from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict
import pandas as pd


# Canonical output columns we want downstream (stable API)
REFMET_CANONICAL_COLUMNS: List[str] = [
    "Input name",
    "RefMet_ID",
    "RefMet name",
    "Formula",
    "Exact mass",
    "Super class",
    "Main class",
    "Sub class",
    "PubChem_CID",
    "ChEBI_ID",
    "HMDB_ID",
    "LM_ID",
    "KEGG_ID",
    "INCHI_KEY",
    "SMILES",
]

# If RefMet uses slightly different headers, map them to canonical.
REFMET_COLUMN_ALIASES: Dict[str, str] = {
    "Exact_mass": "Exact mass",
    "RefMet_name": "RefMet name",
    "Super_class": "Super class",
    "Main_class": "Main class",
    "Sub_class": "Sub class",
    "PubChem CID": "PubChem_CID",
    "ChEBI ID": "ChEBI_ID",
    "HMDB ID": "HMDB_ID",
    "InChIKey": "INCHI_KEY",
}


@dataclass(frozen=True)
class RefMetSchema:
    canonical_columns: List[str] = None

    def __post_init__(self):
        if self.canonical_columns is None:
            object.__setattr__(self, "canonical_columns", REFMET_CANONICAL_COLUMNS)

    def normalize_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        - Rename alias columns to canonical
        - Ensure all canonical columns exist (add missing as NA)
        - Order columns: canonical first, then extra columns (if any)
        """
        if df is None or df.empty:
            return pd.DataFrame(columns=self.canonical_columns)

        out = df.copy()

        # rename aliases
        out.rename(columns={c: REFMET_COLUMN_ALIASES.get(c, c) for c in out.columns}, inplace=True)

        # ensure canonical columns exist
        for c in self.canonical_columns:
            if c not in out.columns:
                out[c] = pd.NA

        # reorder: canonical first then extras
        extras = [c for c in out.columns if c not in self.canonical_columns]
        out = out[self.canonical_columns + extras]

        return out
