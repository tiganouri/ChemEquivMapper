from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple
from collections import deque

import pandas as pd
import pronto

from GEPC.context import PipelineContext
from GEPC.resources import ChebiResourceConfig, resolve_chebi_obo
from GEPC.resources import ReactomeTXTConfig, KEGGTSVConfig  # KEGG cfg unused but ctx loader expects it

from GEPC.steps.step2_chebi.supplementary import parse_chebi_list, sort_chebi_ids
from .schema import Step5Schema

logger = logging.getLogger(__name__)


def _first_nonempty_cell(df: pd.DataFrame, idx, cols: Sequence[str]) -> str:
    for c in cols:
        if c in df.columns:
            v = df.at[idx, c]
            s = "" if v is None else str(v).strip()
            if s and s.lower() != "nan":
                return s
    return ""


def _get_term(onto: pronto.Ontology, chebi_id: str) -> Optional[pronto.Term]:
    try:
        return onto[chebi_id]
    except Exception:
        return None


@dataclass
class Step5IsAGeneralisation:
    """
    Step5: is_a-based generalisation for *previously unmapped* rows.

    - No RO neighbors
    - No KEGG
    - Reactome mapping only (via shared Reactome index)
    - Finds ALL closest mapped ancestors (minimum BFS distance)
    - Produces dataframe columns + a structured report (for PDF)
    """
    repo_root: Path

    # resources
    chebi_cfg: ChebiResourceConfig = ChebiResourceConfig(source="bundled")
    
    _ctx: Optional[PipelineContext] = field(default=None, init=False)

    # behaviour
    max_is_a_depth: int = 20
    excluded_ancestors: Set[str] = field(
        default_factory=lambda: {
            "CHEBI:23367",  # molecular entity
            "CHEBI:24431",  # chemical entity
            "CHEBI:61120",  # nucleobase-containing molecular entity (example)
        }
    )

    # how to decide unmapped rows
    prior_reactome_cols: Sequence[str] = (
        "Reactome_Pathways_Step1",
        "Reactome_Pathways_Step2",
        "Reactome_Pathways_Step3",
        "Reactome_Pathways_Step4",
    )

    # which chebi columns to use as starting point (best available)
    base_chebi_cols: Sequence[str] = ("CHEBI_ID_Step2",)

    # internals
    _onto: Optional[pronto.Ontology] = field(default=None, init=False)
    _reactome_index: Any = field(default=None, init=False)

    def set_shared_resources(
        self,
        *,
        onto: Optional[pronto.Ontology] = None,
        ctx: Optional[PipelineContext] = None,
    ) -> None:
        """Inject already-loaded shared resources (e.g., from Step2/3/4) to avoid reloading."""
        if onto is not None:
            self._onto = onto
        if ctx is not None:
            self._ctx = ctx
            # If ctx already has a reactome index, reuse it directly
            if getattr(ctx, "reactome_index", None) is not None:
                self._reactome_index = ctx.reactome_index

    def _ensure_ontology(self) -> None:
        if self._onto is not None:
            return
        cache_dir = Path(self.repo_root) / ".cache" / "step5_isa_generalisation"
        cache_dir.mkdir(parents=True, exist_ok=True)
        obo_path = resolve_chebi_obo(self.chebi_cfg, cache_dir=cache_dir)
        if obo_path is None:
            raise RuntimeError("Step5: chebi.obo could not be resolved (bundled/local/download).")
        self._onto = pronto.Ontology(str(obo_path))

    def _ensure_reactome_index(self) -> None:
        if self._reactome_index is not None:
            return

        # Prefer an existing ctx.reactome_index if available
        if self._ctx is not None and getattr(self._ctx, "reactome_index", None) is not None:
            self._reactome_index = self._ctx.reactome_index
            return

        # Otherwise create ctx (Reactome only; KEGG cfg required by loader but not used here)
        if self._ctx is None:
            self._ctx = PipelineContext(repo_root=Path(self.repo_root)).load_resources_default(
                reactome_cfg=ReactomeTXTConfig(source="download"),
                kegg_cfg=KEGGTSVConfig(source="bundled"),
                species="Homo sapiens",
            )

        self._reactome_index = self._ctx.reactome_index


    def _chebi_pathways(self, chebi_id: str) -> Set[str]:
        """
        Uses ReactomeIndex. Adjust this ONE line if your ReactomeIndex method name differs.
        """
        self._ensure_reactome_index()
        idx = self._reactome_index

        # Common method patterns:
        if hasattr(idx, "get_pathways"):
            return set(idx.get_pathways(chebi_id))  # type: ignore
        if hasattr(idx, "pathways_for_chebi"):
            return set(idx.pathways_for_chebi(chebi_id))  # type: ignore

        raise AttributeError("ReactomeIndex has no get_pathways/pathways_for_chebi method.")

    def _iter_is_a_parents(self, term: pronto.Term) -> List[pronto.Term]:
        parents: Set[pronto.Term] = set()

        if hasattr(term, "parents"):
            try:
                parents |= set(term.parents)
            except Exception:
                pass

        if not parents and hasattr(term, "superclasses"):
            try:
                for t in term.superclasses(distance=1, with_self=False):
                    parents.add(t)
            except TypeError:
                pass
            except Exception:
                pass

        def _tid(t: pronto.Term) -> str:
            try:
                return t.id or ""
            except Exception:
                return ""

        return sorted(parents, key=_tid)

    def _find_closest_mapped_ancestors(
        self,
        start_chebi: str,
    ) -> List[Tuple[str, List[str], Set[str]]]:
        """
        BFS up is_a parents and return ALL mapped ancestors found at MIN depth.
        Each result is (mapped_chebi, path_ids, pathways).
        Excluded ancestors prune branches.
        """
        self._ensure_ontology()
        assert self._onto is not None

        start_term = _get_term(self._onto, start_chebi)
        if start_term is None:
            return []

        visited: Set[pronto.Term] = {start_term}
        q: deque[Tuple[pronto.Term, int, List[str]]] = deque([(start_term, 0, [start_chebi])])

        found: Dict[str, Tuple[List[str], Set[str]]] = {}
        min_depth_found: Optional[int] = None

        while q:
            node, depth, path_ids = q.popleft()

            if min_depth_found is not None and depth > min_depth_found:
                break

            cid = getattr(node, "id", None)
            if isinstance(cid, str) and cid.startswith("CHEBI:"):
                if cid in self.excluded_ancestors:
                    continue

                pathways = self._chebi_pathways(cid)
                if pathways:
                    if min_depth_found is None:
                        min_depth_found = depth
                    if depth == min_depth_found and cid not in found:
                        found[cid] = (path_ids, pathways)
                    continue

            if min_depth_found is None and depth < self.max_is_a_depth:
                for parent in self._iter_is_a_parents(node):
                    if parent not in visited:
                        visited.add(parent)
                        pid = getattr(parent, "id", None)
                        new_path = path_ids + [pid] if isinstance(pid, str) and pid.startswith("CHEBI:") else path_ids
                        q.append((parent, depth + 1, new_path))

        out: List[Tuple[str, List[str], Set[str]]] = []
        for mapped_chebi, (path, pws) in found.items():
            out.append((mapped_chebi, path, pws))
        out.sort(key=lambda x: x[0])
        return out

    def _row_is_mapped(self, df: pd.DataFrame, idx) -> bool:
        for c in self.prior_reactome_cols:
            if c in df.columns:
                s = "" if df.at[idx, c] is None else str(df.at[idx, c]).strip()
                if s and s.lower() != "nan":
                    return True
        return False

    def run(self, df: pd.DataFrame, generate_report: bool = True, *, ctx: Optional[PipelineContext] = None) -> Tuple[pd.DataFrame, List[Dict[str, Any]]]:
        if ctx is not None:
            self._ctx = ctx

        df_out = df.copy()
        schema = Step5Schema()
        df_out = schema.ensure_columns(df_out)

        self._ensure_ontology()
        self._ensure_reactome_index()

        report_entries: List[Dict[str, Any]] = []

        for idx in df_out.index:
            # Only previously unmapped rows
            if self._row_is_mapped(df_out, idx):
                continue

            base_cell = _first_nonempty_cell(df_out, idx, self.base_chebi_cols)
            base_ids = parse_chebi_list(base_cell)
            if not base_ids:
                continue

            mapped_ancestors: Set[str] = set()
            mapped_pathways: Set[str] = set()

            for base_chebi in sorted(base_ids):
                results = self._find_closest_mapped_ancestors(base_chebi)

                for mapped_chebi, path_ids, pathways in results:
                    mapped_ancestors.add(mapped_chebi)
                    mapped_pathways |= pathways

                    if generate_report:
                        # add labels with term names
                        path_labels: List[str] = []
                        for cid in path_ids:
                            label = cid
                            term = _get_term(self._onto, cid) if self._onto is not None else None
                            if term is not None and getattr(term, "name", None):
                                label = f"{cid} ({term.name})"
                            path_labels.append(label)

                        mapped_name = None
                        term_m = _get_term(self._onto, mapped_chebi) if self._onto is not None else None
                        if term_m is not None and getattr(term_m, "name", None):
                            mapped_name = term_m.name

                        report_entries.append(
                            {
                                "row_index": int(idx),
                                "base_chebi": base_chebi,
                                "path_ids": path_ids,
                                "path_labels": path_labels,
                                "mapped_chebi": mapped_chebi,
                                "mapped_chebi_name": mapped_name,
                                "pathways": sorted(pathways),
                            }
                        )

            if mapped_ancestors and mapped_pathways:
                df_out.at[idx, schema.chebi_generalisation_col] = ",".join(sort_chebi_ids(mapped_ancestors))
                df_out.at[idx, schema.reactome_generalisation_col] = ",".join(sorted(mapped_pathways))
                df_out.at[idx, schema.n_pathways_col] = len(mapped_pathways)

        df_out = schema.normalize(df_out)
        logger.info("Step5IsAGeneralisation completed (Reactome-only, unmapped rows).")
        return df_out, report_entries
