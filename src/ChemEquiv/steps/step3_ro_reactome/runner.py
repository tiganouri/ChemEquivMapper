from __future__ import annotations

import logging
from dataclasses import dataclass, field
from collections import deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pronto
import pandas as pd

from GEPC.context import PipelineContext
from GEPC.resources import ChebiResourceConfig, resolve_chebi_obo
from GEPC.resources import ReactomeTXTConfig, KEGGTSVConfig

from GEPC.steps.step2_chebi.kegg_map import get_chebi_to_kegg_map
from GEPC.steps.step2_chebi.supplementary import (
    parse_chebi_list,
    parse_kegg_list,
    sort_chebi_ids,
    sort_kegg_ids,
)

logger = logging.getLogger(__name__)

RO_EQUIVALENCE_RELS: Sequence[str] = [
    "RO:0018033",  # is conjugate base of
    "RO:0018034",  # is conjugate acid of
    "RO:0018036",  # is tautomer of
    "RO:0018039",  # is enantiomer of
]


def _get_term(onto: pronto.Ontology, chebi_id: str) -> Optional[pronto.Term]:
    try:
        return onto[chebi_id]
    except Exception:
        return None


def _get_relationship(onto: pronto.Ontology, rel_id: str):
    if hasattr(onto, "get_relationship"):
        try:
            return onto.get_relationship(rel_id)
        except Exception:
            pass
    try:
        return onto.get(rel_id)
    except Exception:
        return None


def _term_neighbors_for_rel(term: pronto.Term, rel) -> Set[pronto.Term]:
    targets: Set[pronto.Term] = set()

    # outgoing
    for attr in ("relations", "relationships"):
        if hasattr(term, attr):
            m = getattr(term, attr)
            try:
                targets |= set(m.get(rel, set()))
            except Exception:
                pass

    # incoming
    for attr in ("relation_of", "related_by"):
        if hasattr(term, attr):
            m = getattr(term, attr)
            try:
                targets |= set(m.get(rel, set()))
            except Exception:
                pass

    return targets


@dataclass
class Step3RO:
    """
    Step 3: RO-based ChEBI expansion + KEGG enrichment + pathway mapping via shared mapper.

    Required input columns (from Step2 output):
      - CHEBI_ID_Step2 (comma-separated)
      - KEGG_ID_Step2  (comma-separated, optional)

    Output columns:
      - CHEBI_RO_Equivalents
      - CHEBI_ID_Step3
      - KEGG_ID_Step3
      - ChEBI_ID_Step3_reactome
      - Reactome_Pathways_Step3
      - KEGG_Compounds_Step3
      - KEGG_Pathways_Step3
    """

    repo_root: Path

    # --- ChEBI ontology controls (needed for RO traversal + CHEBI->KEGG) ---
    chebi_cfg: ChebiResourceConfig = ChebiResourceConfig(source="bundled")

    # --- RO traversal ---
    equivalence_rels: Sequence[str] = field(default_factory=lambda: list(RO_EQUIVALENCE_RELS))
    max_depth: int = 20

    # --- internal caches (lazy) ---
    _onto: Optional[pronto.Ontology] = field(default=None, init=False)
    _rel_objs: Optional[List] = field(default=None, init=False)
    _chebi_to_kegg: Optional[Dict[str, Set[str]]] = field(default=None, init=False)
    _equiv_cache: Dict[str, Set[str]] = field(default_factory=dict, init=False)

    def set_shared_resources(
        self,
        *,
        onto: Optional[pronto.Ontology] = None,
        chebi_to_kegg: Optional[Dict[str, Set[str]]] = None,
    ) -> None:
        """
        Inject already-loaded shared resources (e.g., from Step2) to avoid reloading.
        """
        if onto is not None:
            self._onto = onto
        if chebi_to_kegg is not None:
            self._chebi_to_kegg = chebi_to_kegg

    def _cache_dir(self) -> Path:
        return Path(self.repo_root) / ".cache" / "step3_ro"

    def _ensure_ontology(self, *, ctx: Optional[PipelineContext] = None) -> None:
        # If we already have an ontology, do not reload it
        if self._onto is None:
            cache_dir = self._cache_dir()
            cache_dir.mkdir(parents=True, exist_ok=True)

            obo_path = resolve_chebi_obo(self.chebi_cfg, cache_dir=cache_dir)
            if obo_path is None:
                raise RuntimeError("Step3: chebi.obo could not be resolved from bundled/local/download.")

            self._onto = pronto.Ontology(str(obo_path))

        # If we already have CHEBI->KEGG, do not rebuild it
        if self._chebi_to_kegg is not None:
            return

        # If ctx provides chebi_to_kegg, reuse it (no OBO needed)
        if ctx is not None and getattr(ctx, "chebi_to_kegg", None):
            self._chebi_to_kegg = {c: set(v) for c, v in ctx.chebi_to_kegg.items()}
            return

        # Otherwise, build CHEBI->KEGG from ontology (sqlite cached)
        cache_dir = self._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)

        obo_path = resolve_chebi_obo(self.chebi_cfg, cache_dir=cache_dir)
        if obo_path is None:
            raise RuntimeError("Step3: chebi.obo could not be resolved from bundled/local/download.")

        map_db = cache_dir / "chebi_to_kegg.sqlite"
        self._chebi_to_kegg = get_chebi_to_kegg_map(obo_path=obo_path, db_path=map_db, ont=self._onto)

    def _ensure_rel_objs(self, *, ctx: Optional[PipelineContext] = None) -> None:
        if self._rel_objs is not None:
            return

        self._ensure_ontology(ctx=ctx)
        assert self._onto is not None

        rel_objs: List = []
        for rid in self.equivalence_rels:
            rel = _get_relationship(self._onto, rid)
            if rel is not None:
                rel_objs.append(rel)
            else:
                logger.warning(f"Step3: Could not resolve RO relationship '{rid}' from ontology.")

        if not rel_objs:
            raise RuntimeError("Step3: No RO equivalence relationships resolved.")
        self._rel_objs = rel_objs

    def _equivalence_closure_single(self, chebi_id: str, *, ctx: Optional[PipelineContext] = None) -> Set[str]:
        chebi_id = str(chebi_id).strip().upper()
        if not chebi_id.startswith("CHEBI:"):
            return set()

        if chebi_id in self._equiv_cache:
            return set(self._equiv_cache[chebi_id])

        self._ensure_rel_objs(ctx=ctx)
        assert self._onto is not None
        assert self._rel_objs is not None

        start = _get_term(self._onto, chebi_id)
        if start is None:
            self._equiv_cache[chebi_id] = set()
            return set()

        visited: Set[pronto.Term] = {start}
        q: deque[Tuple[pronto.Term, int]] = deque([(start, 0)])

        while q:
            node, d = q.popleft()
            if d >= self.max_depth:
                continue

            for rel in self._rel_objs:
                for nb in _term_neighbors_for_rel(node, rel):
                    if nb not in visited:
                        visited.add(nb)
                        q.append((nb, d + 1))

        out: Set[str] = set()
        for t in visited:
            tid = getattr(t, "id", None)
            if isinstance(tid, str) and tid.startswith("CHEBI:") and tid != chebi_id:
                out.add(tid)

        self._equiv_cache[chebi_id] = set(out)
        return out

    def _expand_ro_equivalents(self, chebi_ids: Set[str], *, ctx: Optional[PipelineContext] = None) -> Set[str]:
        all_eq: Set[str] = set()
        for cid in chebi_ids:
            all_eq |= self._equivalence_closure_single(cid, ctx=ctx)
        return all_eq - chebi_ids

    def _map_chebi_to_kegg(self, chebi_ids: Set[str], *, ctx: Optional[PipelineContext] = None) -> Set[str]:
        self._ensure_ontology(ctx=ctx)
        assert self._chebi_to_kegg is not None
        out: Set[str] = set()
        for cid in chebi_ids:
            out |= self._chebi_to_kegg.get(cid, set())
        return out

    def run(self, 
        df: pd.DataFrame, 
        chebi_col_step2: str="CHEBI_ID_Step2", 
        kegg_col_step2: str="KEGG_ID_Step2", 
        *, 
        ctx: Optional[PipelineContext] = None
    ) -> pd.DataFrame:
        if chebi_col_step2 not in df.columns:
            raise KeyError(f"Step3 requires '{chebi_col_step2}' in the input DataFrame.")

        df_out = df.copy()

        # Prepare output columns
        df_out["CHEBI_RO_Equivalents"] = ""
        df_out["CHEBI_ID_Step3"] = ""
        df_out["KEGG_ID_Step3"] = ""

        # Ensure ctx first (so we can reuse ctx.chebi_to_kegg)
        if ctx is None:
            ctx = PipelineContext(repo_root=Path(self.repo_root)).load_resources_default(
                reactome_cfg=ReactomeTXTConfig(source="download"),
                kegg_cfg=KEGGTSVConfig(source="bundled"),
                species="Homo sapiens",
            )

        # Ensure ontology + maps using ctx (no reload if provided)
        self._ensure_rel_objs(ctx=ctx)
        self._ensure_ontology(ctx=ctx)

        # Provide mapping to ctx (optional but nice for later steps)
        if getattr(ctx, "chebi_to_kegg", None) is None and self._chebi_to_kegg is not None:
            ctx.chebi_to_kegg = {c: sorted(list(v)) for c, v in self._chebi_to_kegg.items()}

        mapper = ctx.build_mapper()

        for idx in df_out.index:
            step2_chebis = parse_chebi_list(df_out.at[idx, chebi_col_step2])
            step2_kegg = (
                parse_kegg_list(df_out.at[idx, kegg_col_step2])
                if kegg_col_step2 in df_out.columns
                else set()
            )

            ro_eq = self._expand_ro_equivalents(step2_chebis, ctx=ctx) if step2_chebis else set()
            step3_chebis = step2_chebis | ro_eq

            df_out.at[idx, "CHEBI_RO_Equivalents"] = ",".join(sort_chebi_ids(ro_eq)) if ro_eq else ""
            df_out.at[idx, "CHEBI_ID_Step3"] = ",".join(sort_chebi_ids(step3_chebis)) if step3_chebis else ""

            derived_kegg = self._map_chebi_to_kegg(step3_chebis, ctx=ctx) if step3_chebis else set()
            kegg_step3 = step2_kegg | derived_kegg
            df_out.at[idx, "KEGG_ID_Step3"] = ",".join(sort_kegg_ids(kegg_step3)) if kegg_step3 else ""

            # Map to Reactome/KEGG pathways via shared mapper (Step3 columns)
            mapper.apply_to_df(
                df_out,
                idx,
                step=3,
                chebi_ids=df_out.at[idx, "CHEBI_ID_Step3"],
                kegg_compound_ids=df_out.at[idx, "KEGG_ID_Step3"],
                dedup_pathways=False,   # keep duplicates (your requirement)
                dedup_compounds=True,
            )

        logger.info("Step3RO completed.")
        return df_out
