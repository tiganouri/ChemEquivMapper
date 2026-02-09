from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any
import pronto

import pandas as pd

from .schema import Step2Schema
from .name_index import NameIndexConfig, ensure_name_index_db, load_indices_from_db
from .kegg_map import get_chebi_to_kegg_map
from .supplementary import (
    norm,
    norm_remove_special,
    parse_chebi_list,
    parse_kegg_list,
    sort_chebi_ids,
    sort_kegg_ids,
)

from ChemEquiv.context import PipelineContext
from ChemEquiv.resources import (
    ReactomeTXTConfig,
    KEGGTSVConfig,
    ChebiResourceConfig,
    resolve_chebi_obo,
)


def _first_nonempty(*vals: object) -> str:
    for v in vals:
        s = "" if v is None else str(v).strip()
        if s and s.lower() != "nan" and s != "-":
            return s
    return ""


def _join(ids: List[str]) -> str:
    return ",".join([x for x in ids if x])


def _lookup_from_indices(
    name: str,
    exact: Dict[str, List[str]],
    fuzzy: Dict[str, List[str]],
) -> Tuple[List[str], List[str], List[str]]:
    """
    Returns (exact_hits, fuzzy_hits, all_hits) as lists.
    Uses your existing normalization functions:
      - exact key: norm(name)
      - fuzzy key: norm_remove_special(name)
    """
    if not name:
        return [], [], []
    k1 = norm(name)
    k2 = norm_remove_special(name)

    exact_hits = exact.get(k1, []) if k1 else []
    fuzzy_hits = fuzzy.get(k2, []) if k2 else []

    all_hits = list(exact_hits) + [x for x in fuzzy_hits if x not in exact_hits]
    return exact_hits, fuzzy_hits, all_hits

def _norm_chebi_token(tok: str) -> str:
    s = str(tok).strip().upper()
    if not s or s == "-" or s == "NAN":
        return ""
    if s.isdigit():
        return f"CHEBI:{s}"
    if s.startswith("CHEBI:"):
        return s
    if s.startswith("CHEBI_"):
        return s.replace("CHEBI_", "CHEBI:")
    # if something odd comes in, keep it empty to avoid junk
    return ""

def _parse_chebi_any(x: object) -> List[str]:
    """
    Parses comma-separated values and normalizes tokens to CHEBI: format.
    Works for Step1 (digits) and Step2 (already CHEBI:xxxx).
    """
    if x is None:
        return []
    s = str(x).strip()
    if not s or s.lower() == "nan" or s == "-":
        return []
    toks = [t.strip() for t in s.split(",") if t.strip()]
    out = []
    for t in toks:
        nt = _norm_chebi_token(t)
        if nt:
            out.append(nt)
    return out

@dataclass
class Step2ChebiLookup:
    repo_root: Path

    # where Step1 put names (RefMet output)
    input_name_col: str = "Input name"
    standardized_name_col: str = "Standardized name"

    # ChEBI ontology source strategy
    chebi_cfg: ChebiResourceConfig = ChebiResourceConfig(source="bundled")

    # where to store sqlite indices/caches
    cache_dir: Optional[Path] = None

    # ---- cached heavy resources (loaded once per Step2ChebiLookup instance) ----
    _obo_path: Optional[Path] = field(default=None, init=False, repr=False)
    _onto: Optional[pronto.Ontology] = field(default=None, init=False, repr=False)
    _exact_idx: Optional[Dict[str, List[str]]] = field(default=None, init=False, repr=False)
    _fuzzy_idx: Optional[Dict[str, List[str]]] = field(default=None, init=False, repr=False)
    _chebi_to_kegg: Optional[Dict[str, Set[str]]] = field(default=None, init=False, repr=False)

    def _cache_dir(self) -> Path:
        return Path(self.cache_dir) if self.cache_dir is not None else (Path(self.repo_root) / ".cache" / "step2_chebi")
        
    def _ensure_loaded(self) -> None:
        """
        Load ChEBI ontology + indices + CHEBI->KEGG mapping exactly once per Step2ChebiLookup instance.
        Subsequent runs reuse them (avoids RAM blowups).
        """
        # If already loaded, do nothing
        if (
            self._onto is not None
            and self._exact_idx is not None
            and self._fuzzy_idx is not None
            and self._chebi_to_kegg is not None
        ):
            return

        cache_dir = self._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)

        # 1) Resolve chebi.obo.gz once
        if self._obo_path is None:
            obo_path = resolve_chebi_obo(self.chebi_cfg, cache_dir=cache_dir)
            if obo_path is None:
                raise FileNotFoundError("Could not resolve chebi.obo (bundled/local/download failed).")
            self._obo_path = Path(obo_path)

        # 2) Load ontology once
        if self._onto is None:
            self._onto = pronto.Ontology(str(self._obo_path))

        # 3) Ensure/load indices once
        if self._exact_idx is None or self._fuzzy_idx is None:
            name_db = cache_dir / "chebi_name_index.sqlite"
            nicfg = NameIndexConfig(db_path=name_db)

            ensure_name_index_db(nicfg, obo_path=self._obo_path, force_rebuild=False, onto=self._onto)
            self._exact_idx, self._fuzzy_idx = load_indices_from_db(nicfg)

        # 4) Build/load CHEBI -> KEGG mapping once
        if self._chebi_to_kegg is None:
            kegg_db = cache_dir / "chebi_kegg.sqlite"
            self._chebi_to_kegg = get_chebi_to_kegg_map(
                obo_path=self._obo_path,
                db_path=kegg_db,
                ont=self._onto,
            )

    def run(
        self, 
        meta_step1: pd.DataFrame, 
        *, 
        ctx: Optional[PipelineContext] = None
    ) -> Tuple[pd.DataFrame, Dict[str, Set[str]], pronto.Ontology]:
        """
        Input: Step1 meta dataframe (should have Input name + Standardized name)
        Output:
          - meta_step2 dataframe with Step2 CHEBI/KEGG + mapping columns
          - chebi_to_kegg dict (CHEBI -> set(Cxxxxx))
          - onto (pronto.Ontology) cached + reused across datasets
        """
        out = meta_step1.copy()
        schema = Step2Schema()
        out = schema.ensure_columns(out)
        
        # --- Normalize Step1 IDs so they match Step2 formatting ---
        if "ChEBI_ID" in out.columns:
            out["ChEBI_ID"] = out["ChEBI_ID"].apply(lambda v: _join(sort_chebi_ids(_parse_chebi_any(v))))
        if "KEGG_ID" in out.columns:
            # keep KEGG as simple cleaned comma list
            out["KEGG_ID"] = out["KEGG_ID"].apply(lambda v: _join(sort_kegg_ids(parse_kegg_list(v))))

        # Load heavy resources once per Step2ChebiLookup instance
        self._ensure_loaded()
        onto = self._onto
        exact_idx = self._exact_idx
        fuzzy_idx = self._fuzzy_idx
        chebi_to_kegg = self._chebi_to_kegg
    
        if onto is None or exact_idx is None or fuzzy_idx is None or chebi_to_kegg is None:
            raise RuntimeError("Step2 resources not loaded correctly (ontology/indices/map missing).")

        # ---- 1) Name -> CHEBI lookup using standardized then input name ----
        for idx in out.index:
            std_name = _first_nonempty(out.at[idx, self.standardized_name_col])
            in_name = _first_nonempty(out.at[idx, self.input_name_col])
    
            # lookup both (if present)
            ex_std, fu_std, all_std = _lookup_from_indices(std_name, exact_idx, fuzzy_idx)
            ex_in,  fu_in,  all_in  = _lookup_from_indices(in_name,  exact_idx, fuzzy_idx)
    
            # prefer standardized hits first
            exact_hits = list(ex_std) + [x for x in ex_in if x not in ex_std]
            fuzzy_hits = list(fu_std) + [x for x in fu_in if x not in fu_std]
    
            all_hits = list(all_std)
            for x in all_in:
                if x not in all_hits:
                    all_hits.append(x)
    
            out.at[idx, schema.chebi_exact_col] = _join(exact_hits)
            out.at[idx, schema.chebi_fuzzy_col] = _join(fuzzy_hits)
            out.at[idx, schema.chebi_all_col]   = _join(all_hits)
    
            # Step2 hits (from name index)
            chebi_step2 = set(sort_chebi_ids(_parse_chebi_any(_join(all_hits))))

            # Step1 hits (RefMet) - already normalized above
            chebi_step1 = set(_parse_chebi_any(out.at[idx, "ChEBI_ID"])) if "ChEBI_ID" in out.columns else set()

            # Union: keep everything
            chebi_union = sort_chebi_ids(list(chebi_step1 | chebi_step2))
            out.at[idx, schema.chebi_final_col] = _join(chebi_union)

        
        # ---- 2) Fill KEGG_ID_Step2 using CHEBI_ID_Step2 ----
        for idx in out.index:
            chebis = parse_chebi_list(out.at[idx, schema.chebi_final_col])
    
            kegg_hits: List[str] = []
            for c in sort_chebi_ids(list(chebis)):
                kegg_hits.extend(sorted(list(chebi_to_kegg.get(c, set()))))
    
            kegg_step2 = set(sort_kegg_ids(parse_kegg_list(",".join(kegg_hits))))
            kegg_step1 = set(parse_kegg_list(out.at[idx, "KEGG_ID"])) if "KEGG_ID" in out.columns else set()

            kegg_union = sort_kegg_ids(list(kegg_step1 | kegg_step2))
            out.at[idx, schema.kegg_final_col] = _join(kegg_union)

            # keep aligned for now
            out.at[idx, schema.kegg_exact_col] = out.at[idx, schema.kegg_final_col]
            out.at[idx, schema.kegg_fuzzy_col] = out.at[idx, schema.kegg_final_col]
            out.at[idx, schema.kegg_all_col]   = out.at[idx, schema.kegg_final_col]

    
            # keep aligned for now
            out.at[idx, schema.kegg_exact_col] = out.at[idx, schema.kegg_final_col]
            out.at[idx, schema.kegg_fuzzy_col] = out.at[idx, schema.kegg_final_col]
            out.at[idx, schema.kegg_all_col]   = out.at[idx, schema.kegg_final_col]
    
        out = schema.normalize_joined_ids(out)
    
        # ---- 3) Mapper: map Step2 CHEBI/KEGG -> Reactome/KEGG pathways ----
        if ctx is None:
            ctx = PipelineContext(repo_root=Path(self.repo_root)).load_resources_default(
                reactome_cfg=ReactomeTXTConfig(source="download"),
                kegg_cfg=KEGGTSVConfig(source="bundled"),
                species="Homo sapiens",
            )

        # 6) Mapper: map Step2 CHEBI/KEGG -> Reactome/KEGG pathways
        if ctx is None:
            ctx = PipelineContext(repo_root=Path(self.repo_root)).load_resources_default(
                reactome_cfg=ReactomeTXTConfig(source="download"),
                kegg_cfg=KEGGTSVConfig(source="bundled"),
                species="Homo sapiens",
            )
            
        # Provide CHEBI->KEGG expansion to mapper (as list[str])
        # (Keep it lightweight; if memory is tight you can avoid copying by caching this too.)
        ctx.chebi_to_kegg = {c: sorted(list(ks)) for c, ks in chebi_to_kegg.items()}
        mapper = ctx.build_mapper()
    
        for idx in out.index:
            mapper.apply_to_df(
                out,
                idx,
                step=2,
                chebi_ids=out.at[idx, schema.chebi_final_col],
                kegg_compound_ids=out.at[idx, schema.kegg_final_col],
                dedup_pathways=False,
                dedup_compounds=True,
            )
    
        return out, chebi_to_kegg, onto
