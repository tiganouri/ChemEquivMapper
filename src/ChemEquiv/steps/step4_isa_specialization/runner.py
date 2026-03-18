from __future__ import annotations

import logging
from dataclasses import dataclass, field
from collections import deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd
import pronto

from ChemEquiv.context import PipelineContext
from ChemEquiv.resources import (
    ChebiResourceConfig,
    resolve_chebi_obo,
    ReactomeTXTConfig,
    KEGGTSVConfig,
)

from ChemEquiv.steps.step2_chebi.kegg_map import get_chebi_to_kegg_map
from ChemEquiv.steps.step2_chebi.supplementary import (
    parse_chebi_list,
    parse_kegg_list,
    sort_chebi_ids,
    sort_kegg_ids,
)

from .isa_traversal import get_is_a_children_for_ids
from .schema import Step4Schema

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

    for attr in ("relations", "relationships"):
        if hasattr(term, attr):
            m = getattr(term, attr)
            try:
                targets |= set(m.get(rel, set()))
            except Exception:
                pass

    for attr in ("relation_of", "related_by"):
        if hasattr(term, attr):
            m = getattr(term, attr)
            try:
                targets |= set(m.get(rel, set()))
            except Exception:
                pass

    return targets


@dataclass
class Step4IsASpecialization:
    """
    Step4: is_a(children) + RO(children) + mapping via shared mapper.

    Input:
      - CHEBI_ID_Step3
      - KEGG_ID_Step3 (optional)

    Output:
      - CHEBI_is_a_Equivalents
      - CHEBI_ID_Step4
      - KEGG_ID_Step4
      - (mapper outputs with step=4)
    """
    repo_root: Path

    # ontology resource strategy
    chebi_cfg: ChebiResourceConfig = ChebiResourceConfig(source="bundled")

    # traversal depths
    max_is_a_depth: int = 1
    max_equiv_depth: int = 20
    equivalence_rels: Sequence[str] = field(default_factory=lambda: list(RO_EQUIVALENCE_RELS))

    # internal caches
    _onto: Optional[pronto.Ontology] = field(default=None, init=False)
    _rel_objs: Optional[List] = field(default=None, init=False)
    _chebi_to_kegg: Optional[Dict[str, Set[str]]] = field(default=None, init=False)
    _equiv_cache: Dict[str, Set[str]] = field(default_factory=dict, init=False)

    # --- mass filtering ---
    mass_col: str = "Exact mass"
    mass_tolerance: float = 20.0

    _mass_cache: Dict[str, Optional[float]] = field(default_factory=dict, init=False)

    def _parse_float(self, x) -> Optional[float]:
        if x is None:
            return None
        try:
            if isinstance(x, str):
                x = x.strip()
                if not x:
                    return None
            v = float(x)
            if v != v:  # NaN check
                return None
            return v
        except Exception:
            return None

    def _get_chebi_mass(self, chebi_id: str, *, ctx: Optional[PipelineContext] = None) -> Optional[float]:
        """
        Read a ChEBI mass from ontology annotations, cached per CHEBI ID.
    
        Accepts multiple annotation property names because pronto / ontology exports
        may expose mass as e.g.:
          - mass
          - monoisotopic_mass
          - exact_mass
          - ...or other mass-like properties
        """
        chebi_id = str(chebi_id).strip().upper()
        if not chebi_id.startswith("CHEBI:"):
            return None
    
        if chebi_id in self._mass_cache:
            return self._mass_cache[chebi_id]
    
        self._ensure_ontology(ctx=ctx)
        assert self._onto is not None
    
        term = _get_term(self._onto, chebi_id)
        if term is None:
            self._mass_cache[chebi_id] = None
            return None
    
        mass_val: Optional[float] = None
    
        try:
            anns = getattr(term, "annotations", []) or []
    
            # Prefer monoisotopic mass first, then generic mass
            preferred_tokens = [
                "monoisotopic_mass",
                "monoisotopic mass",
                "exact_mass",
                "exact mass",
                "chemrof:mass",
                "mass",
            ]
    
            candidates = []
    
            for ann in anns:
                prop = str(getattr(ann, "property", "")).strip().lower()
                lit = getattr(ann, "literal", None)
    
                if lit is None:
                    continue
    
                lit_val = self._parse_float(lit)
                if lit_val is None:
                    continue
    
                candidates.append((prop, lit_val))
    
            # ranked search
            for token in preferred_tokens:
                for prop, val in candidates:
                    if token in prop:
                        mass_val = val
                        break
                if mass_val is not None:
                    break
    
        except Exception:
            mass_val = None
    
        self._mass_cache[chebi_id] = mass_val
        return mass_val

    def _get_row_reference_mass(
        self,
        df: pd.DataFrame,
        idx,
        step3_chebis: Set[str],
        *,
        ctx: Optional[PipelineContext] = None,
    ) -> Optional[float]:
        # 1) prefer Exact mass column if present/parseable
        if self.mass_col in df.columns:
            m = self._parse_float(df.at[idx, self.mass_col])
            if m is not None:
                return m

        # 2) fallback: use first available ontology mass from step3 chebis
        for cid in sort_chebi_ids(step3_chebis):
            m = self._get_chebi_mass(cid, ctx=ctx)
            if m is not None:
                return m

        return None

    def _within_mass_tolerance(
        self,
        ref_mass: Optional[float],
        candidate_chebi: str,
        *,
        ctx: Optional[PipelineContext] = None,
    ) -> bool:
        # If we don't have a reference mass, we can't apply the criterion -> keep old behavior
        if ref_mass is None:
            return True

        cand_mass = self._get_chebi_mass(candidate_chebi, ctx=ctx)

        # STRICT: if reference mass exists, candidate MUST have mass to be considered
        if cand_mass is None:
            return False

        return abs(cand_mass - float(ref_mass)) <= float(self.mass_tolerance)


    def set_shared_resources(
        self,
        *,
        onto: Optional[pronto.Ontology] = None,
        chebi_to_kegg: Optional[Dict[str, Set[str]]] = None,
    ) -> None:
        """Inject already-loaded shared resources (e.g., from Step2/Step3) to avoid reloading."""
        if onto is not None:
            self._onto = onto
        if chebi_to_kegg is not None:
            self._chebi_to_kegg = chebi_to_kegg

    def _cache_dir(self) -> Path:
        return Path(self.repo_root) / ".cache" / "step4_isa"

    def _ensure_ontology_and_maps(self, *, ctx: Optional[PipelineContext] = None) -> None:
        # Ontology: load only if missing
        if self._onto is None:
            cache_dir = self._cache_dir()
            cache_dir.mkdir(parents=True, exist_ok=True)

            obo_path = resolve_chebi_obo(self.chebi_cfg, cache_dir=cache_dir)
            if obo_path is None:
                raise RuntimeError("Step4: chebi.obo could not be resolved (bundled/local/download).")

            self._onto = pronto.Ontology(str(obo_path))

        # CHEBI->KEGG: if already present, stop
        if self._chebi_to_kegg is not None:
            return

        # Prefer ctx.chebi_to_kegg if available
        if ctx is not None and getattr(ctx, "chebi_to_kegg", None):
            self._chebi_to_kegg = {c: set(v) for c, v in ctx.chebi_to_kegg.items()}
            return

        # Otherwise build/load from sqlite cache (and IMPORTANT: reuse ontology if your Step2 change exists)
        cache_dir = self._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)

        map_db = cache_dir / "chebi_to_kegg.sqlite"

        self._chebi_to_kegg = get_chebi_to_kegg_map(obo_path=obo_path, db_path=map_db, ont=self._onto)

    def _ensure_rel_objs(self, *, ctx: Optional[PipelineContext] = None) -> None:
        if self._rel_objs is not None:
            return

        self._ensure_ontology_and_maps(ctx=ctx)
        assert self._onto is not None

        rel_objs: List = []
        for rid in self.equivalence_rels:
            rel = _get_relationship(self._onto, rid)
            if rel is not None:
                rel_objs.append(rel)
            else:
                logger.warning(f"Step4: Could not resolve RO relationship '{rid}' from ontology.")

        if not rel_objs:
            raise RuntimeError("Step4: No RO equivalence relationships resolved.")
        self._rel_objs = rel_objs

    # RO equivalence closure
    def _equivalence_closure_single(self, chebi_id: str, *, ref_mass: Optional[float] = None, ctx: Optional[PipelineContext] = None) -> Set[str]:
        chebi_id = str(chebi_id).strip().upper()
        if not chebi_id.startswith("CHEBI:"):
            return set()

        if ref_mass is None and chebi_id in self._equiv_cache:
            return set(self._equiv_cache[chebi_id])

        self._ensure_rel_objs(ctx=ctx)
        assert self._onto is not None
        assert self._rel_objs is not None

        start = _get_term(self._onto, chebi_id)

        if start is None:
            if ref_mass is None:
                self._equiv_cache[chebi_id] = set()
            return set()


        visited: Set[pronto.Term] = {start}
        q: deque[Tuple[pronto.Term, int]] = deque([(start, 0)])

        while q:
            node, d = q.popleft()
            if d >= self.max_equiv_depth:
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
                if self._within_mass_tolerance(ref_mass, tid, ctx=ctx):
                    out.add(tid)

        # only cache when no mass-filtering is active
        if ref_mass is None:
            self._equiv_cache[chebi_id] = set(out)

        return out


    def _expand_ro_equivalents(self, chebi_ids: Set[str], *, ref_mass: Optional[float] = None, ctx: Optional[PipelineContext] = None) -> Set[str]:
        all_eq: Set[str] = set()
        for cid in chebi_ids:
            all_eq |= self._equivalence_closure_single(cid, ref_mass=ref_mass, ctx=ctx)
        return all_eq

    def _map_chebi_to_kegg(self, chebi_ids: Set[str], *, ctx: Optional[PipelineContext] = None) -> Set[str]:
        self._ensure_ontology_and_maps(ctx=ctx)
        assert self._chebi_to_kegg is not None
        out: Set[str] = set()
        for cid in chebi_ids:
            out |= self._chebi_to_kegg.get(cid, set())
        return out


    def run(self, df: pd.DataFrame, chebi_col_step3: str="CHEBI_ID_Step3", kegg_col_step3: str="KEGG_ID_Step3", *, ctx: Optional[PipelineContext] = None) -> pd.DataFrame:
        if chebi_col_step3 not in df.columns:
            raise KeyError(f"Step4 requires '{chebi_col_step3}' in input dataframe.")

        schema = Step4Schema(chebi_step3_col=chebi_col_step3, kegg_step3_col=kegg_col_step3)
        df_out = schema.ensure_columns(df)

        # shared mapper (Reactome + KEGG)
        if ctx is None:
            ctx = PipelineContext(repo_root=Path(self.repo_root)).load_resources_default(
                reactome_cfg=ReactomeTXTConfig(source="download"),
                kegg_cfg=KEGGTSVConfig(source="bundled"),
                species="Homo sapiens",
            )

        # preload once
        self._ensure_ontology_and_maps(ctx=ctx)
        self._ensure_rel_objs(ctx=ctx)
        assert self._onto is not None
        assert self._chebi_to_kegg is not None
                
        # only set ctx.chebi_to_kegg if not already there
        if not getattr(ctx, "chebi_to_kegg", None):
            ctx.chebi_to_kegg = {c: sorted(list(ks)) for c, ks in self._chebi_to_kegg.items()}
        
        mapper = ctx.build_mapper()

        for idx in df_out.index:
            step3_chebis = parse_chebi_list(df_out.at[idx, chebi_col_step3])
            step3_keggs = (
                parse_kegg_list(df_out.at[idx, kegg_col_step3])
                if kegg_col_step3 in df_out.columns
                else set()
            )

            if not step3_chebis:
                continue

            ref_mass = self._get_row_reference_mass(df_out, idx, step3_chebis, ctx=ctx)

            # 1) is_a children
            isa_children = get_is_a_children_for_ids(
                chebi_ids=step3_chebis,
                onto=self._onto,
                max_depth=self.max_is_a_depth,
            )

            # 1b) MASS FILTER on is_a children (strict when ref_mass exists)
            if ref_mass is not None and isa_children:
                isa_children = {
                    cid for cid in isa_children
                    if self._within_mass_tolerance(ref_mass, cid, ctx=ctx)
                }

            # 2) RO expansion on children
            isa_children_equiv = (
                self._expand_ro_equivalents(isa_children, ref_mass=ref_mass, ctx=ctx)
                if isa_children else set()
            )

            # 3) new specialized chebis (exclude already present)
            isa_specific = (isa_children | isa_children_equiv) - step3_chebis

            df_out.at[idx, schema.chebi_isa_equiv_col] = (
                ",".join(sort_chebi_ids(isa_specific)) if isa_specific else ""
            )

            # 4) union into step4 chebis
            step4_chebis = step3_chebis | isa_specific
            df_out.at[idx, schema.chebi_step4_col] = (
                ",".join(sort_chebi_ids(step4_chebis)) if step4_chebis else ""
            )

            # 5) kegg union
            derived_kegg = self._map_chebi_to_kegg(step4_chebis, ctx=ctx) if step4_chebis else set()
            step4_keggs = step3_keggs | derived_kegg
            df_out.at[idx, schema.kegg_step4_col] = (
                ",".join(sort_kegg_ids(step4_keggs)) if step4_keggs else ""
            )

            # 6) mapper outputs for step=4
            mapper.apply_to_df(
                df_out,
                idx,
                step=4,
                chebi_ids=df_out.at[idx, schema.chebi_step4_col],
                kegg_compound_ids=df_out.at[idx, schema.kegg_step4_col],
                dedup_pathways=False,
                dedup_compounds=True,
            )

        # normalize CHEBI/KEGG union columns; mapper pathway columns not touched
        df_out = schema.normalize_joined_ids(df_out)
        logger.info("Step4IsASpecialization completed.")
        return df_out
