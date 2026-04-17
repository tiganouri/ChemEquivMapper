from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal, List

import re
import shutil
import zipfile

import pandas as pd

from .zip_manager import RefMetZipConfig, resolve_refmet_zip
from .client import RefMetClient
from .cache import SQLiteRefMetCache
from .schema import RefMetSchema

from ...context import PipelineContext


def _normalize_name(x: str) -> str:
    x = "" if x is None else str(x)
    x = x.strip()
    x = re.sub(r"\s+", " ", x)
    return x


def _split_ids_cell(cell: object) -> str:
    """
    Normalize ID cell to comma-separated string.
    The mapper accepts comma-separated strings.
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return ""
    return str(cell).strip()


def _normalize_chebi_cell(cell: object) -> str:
    """
    Normalize ChEBI cell to comma-separated CHEBI:<id> tokens.
    Accepts '89605', 'CHEBI:89605', '89605, 1234', etc.
    """
    s = _split_ids_cell(cell)
    if not s:
        return ""

    parts = [p.strip() for p in s.split(",") if p.strip()]
    out: List[str] = []
    for p in parts:
        if p.isdigit():
            out.append(f"CHEBI:{p}")
        elif p.upper().startswith("CHEBI:"):
            out.append("CHEBI:" + p.split(":", 1)[1].strip())
        else:
            out.append(p)  # leave unknown tokens untouched
    return ",".join(out)


@dataclass
class Step1RefMet:
    """
    Step 1:
    - Input: preprocess metabolite_stats (expects a 'metabolite_id' column)
    - Output: NEW RefMet annotation table (one row per unique input name).
    - Optional: map to Reactome/KEGG pathways using a PRELOADED PipelineContext (standard practice).

    Notes:
    - This step does NOT merge annotations back into metabolite_stats.
    - Any metabolite collapsing/duplicate logic should be handled in preprocess.
    """
    repo_root: Path

    # zip control (optional legacy support)
    refmet_source: Literal["auto", "local", "download"] = "auto"
    refmet_zip_path: Optional[Path] = None
    download_cache_zip: Optional[Path] = None

    # POST settings
    chunk_size: int = 500
    timeout_s: int = 60

    # cache settings
    cache_db: Optional[Path] = None  # if None => repo_root/.cache/refmet_cache.sqlite

    # -------------------------
    # RefMet zip resolution (optional)
    # -------------------------
    def resolve_zip(self) -> Optional[Path]:
        """
        Decide which RefMet zip to use (download with fallback to local).
        Returns zip path or None.
        """
        default_local = self.repo_root / "vendor" / "refmet" / "RefMet_python.zip"
        if self.download_cache_zip is None:
            self.download_cache_zip = self.repo_root / ".cache" / "refmet" / "RefMet_python.zip"

        cfg = RefMetZipConfig(
            source=self.refmet_source,
            default_local_zip=default_local,
            zip_path=self.refmet_zip_path,
            download_cache_zip=self.download_cache_zip,
            timeout_s=30,
        )
        return resolve_refmet_zip(cfg)

    def unzip_refmet(self, zip_path: Path, extract_dir: Path) -> Optional[Path]:
        """
        Unzip safely. Returns extract_dir or None. Never raises.
        Always clears extract_dir first to avoid mixing versions.
        """
        try:
            if extract_dir.exists():
                shutil.rmtree(extract_dir)
            extract_dir.mkdir(parents=True, exist_ok=True)

            with zipfile.ZipFile(zip_path, "r") as z:
                z.extractall(extract_dir)
            return extract_dir
        except Exception:
            return None

    def _cache_path(self) -> Path:
        return self.cache_db if self.cache_db is not None else (self.repo_root / ".cache" / "refmet_cache.sqlite")

    # -------------------------
    # RefMet POST + cache
    # -------------------------
    def run_refmet_post(self, names: List[str]) -> pd.DataFrame:
        """
        RefMet POST with SQLite cache:
          1) load cached hits
          2) POST only missing names (chunked)
          3) normalize schema
          4) store new results to cache
          5) return cached + new

        Never crashes: returns cached-only or empty df if failures.
        """
        names = [_normalize_name(n) for n in names if isinstance(n, str) and n.strip()]
        if not names:
            return pd.DataFrame()

        cache = SQLiteRefMetCache(self._cache_path())
        df_cached, missing = cache.get_many(names)

        if not missing:
            return RefMetSchema().normalize_df(df_cached)

        client = RefMetClient(timeout_s=self.timeout_s)

        all_new = []
        try:
            for i in range(0, len(missing), self.chunk_size):
                chunk = missing[i: i + self.chunk_size]
                df_new = client.fetch(chunk)
                all_new.append(df_new)

            df_new_all = pd.concat(all_new, ignore_index=True) if all_new else pd.DataFrame()
            df_new_all = RefMetSchema().normalize_df(df_new_all)

            cache.put_many(df_new_all, input_col="Input name")

            df_cached = RefMetSchema().normalize_df(df_cached)

            if df_cached.empty:
                return df_new_all
            if df_new_all.empty:
                return df_cached
            return pd.concat([df_cached, df_new_all], ignore_index=True)

        except Exception:
            # if POST fails, still return cached results
            return RefMetSchema().normalize_df(df_cached)

    # -------------------------
    # Pathway mapping (Step1) -- STANDARD: use preloaded PipelineContext
    # -------------------------
    def _map_pathways_step1(self, meta: pd.DataFrame, *, ctx: PipelineContext) -> pd.DataFrame:
        """
        Standard mapping (Step1):
        Uses a PRELOADED PipelineContext (resources resolved once; indices reused).
        Maps:
          - ChEBI -> Reactome pathways
          - KEGG compound -> KEGG pathways
        Writes Step1 columns.
        """
        if meta is None or meta.empty:
            return meta

        if ctx is None or ctx.reactome_index is None or ctx.kegg_index is None:
            raise ValueError(
                "PipelineContext is not loaded. Call ctx.load_resources_default(...) before mapping."
            )

        mapper = ctx.build_mapper()
        out = meta.copy()

        if "ChEBI_ID" not in out.columns:
            out["ChEBI_ID"] = ""
        if "KEGG_ID" not in out.columns:
            out["KEGG_ID"] = ""

        for idx in out.index:
            chebi_cell = _normalize_chebi_cell(out.at[idx, "ChEBI_ID"])
            kegg_cell = _split_ids_cell(out.at[idx, "KEGG_ID"])

            mapper.apply_to_df(
                out,
                idx,
                step=1,
                chebi_ids=chebi_cell,
                kegg_compound_ids=kegg_cell,
                dedup_pathways=False,
                dedup_compounds=True,
            )

        return out

    # -------------------------
    # Public API: metabolite_stats -> Step1 annotation table
    # -------------------------
    def run(
        self,
        metabolite_stats: pd.DataFrame,
        *,
        ctx: Optional[PipelineContext] = None,
        unzip: bool = False,
        output_csv: Optional[Path] = None,
        map_pathways: bool = True,
    ) -> pd.DataFrame:
        """
        Step1:
        - Input: preprocess metabolite_stats (expects a 'metabolite_id' column)
        - Output: NEW RefMet annotation table (one row per unique input name).
        - Standard mapping: if map_pathways=True, ctx MUST be provided and preloaded.
        """
        if metabolite_stats is None or metabolite_stats.empty:
            return metabolite_stats

        df_in = metabolite_stats.copy()

        # 1) Extract names from metabolite_id (preferred), else fall back safely
        if "metabolite_id" in df_in.columns and df_in["metabolite_id"].notna().any():
            names = df_in["metabolite_id"].astype(str).tolist()
        else:
            # fallback: index (only if not default RangeIndex), else first column
            if not isinstance(df_in.index, pd.RangeIndex):
                names = df_in.index.astype(str).tolist()
            else:
                names = df_in.iloc[:, 0].astype(str).tolist()

        names = [_normalize_name(x) for x in names]
        names = [n for n in names if n]

        if not names:
            return pd.DataFrame()

        # 2) Optional unzip (legacy; NOT required for POST)
        if unzip:
            z = self.resolve_zip()
            if z is not None:
                extract_dir = self.repo_root / ".cache" / "refmet" / "RefMet_python_unzipped"
                _ = self.unzip_refmet(z, extract_dir)

        # 3) Deduplicate WITHOUT sorting (preserve first-seen order)
        names_unique = list(dict.fromkeys(names))

        # 4) RefMet POST + cache (schema normalized)
        df_ann = self.run_refmet_post(names_unique)
        df_ann = RefMetSchema().normalize_df(df_ann)  # safe redundancy
        
        if "ChEBI_ID" in df_ann.columns:
            df_ann["ChEBI_ID"] = df_ann["ChEBI_ID"].apply(_normalize_chebi_cell)
    
        # 5) Standard mapping using ctx (required if map_pathways=True)
        if map_pathways and (df_ann is not None) and (not df_ann.empty):
            if ctx is None:
                raise ValueError(
                    "map_pathways=True requires a preloaded PipelineContext. "
                    "Create ctx = PipelineContext(...).load_resources_default(...) and pass it here."
                )
            df_ann = self._map_pathways_step1(df_ann, ctx=ctx)

        # 6) Optional output CSV
        if output_csv is not None:
            output_csv = Path(output_csv)
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            df_ann.to_csv(output_csv, index=False)

        return df_ann
