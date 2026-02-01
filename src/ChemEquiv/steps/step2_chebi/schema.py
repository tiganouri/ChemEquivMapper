from __future__ import annotations

from dataclasses import dataclass
from typing import List

import pandas as pd

from .supplementary import (
    parse_chebi_list,
    parse_kegg_list,
    sort_chebi_ids,
    sort_kegg_ids,
)


@dataclass(frozen=True)
class Step2Schema:
    """
    Step2 output schema helper:
    - ensures output columns exist
    - normalizes/cleans CHEBI/KEGG columns
    """

    chebi_exact_col: str = "CHEBI_Exact_From_Names"
    chebi_fuzzy_col: str = "CHEBI_Fuzzy_From_Names"
    chebi_all_col: str = "CHEBI_From_Names_All"
    chebi_final_col: str = "CHEBI_ID_Step2"

    kegg_exact_col: str = "KEGG_Exact_From_Names"
    kegg_fuzzy_col: str = "KEGG_Fuzzy_From_Names"
    kegg_all_col: str = "KEGG_From_Names_All"
    kegg_final_col: str = "KEGG_ID_Step2"

    def ensure_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        df2 = df.copy()
        for c in [
            self.chebi_exact_col,
            self.chebi_fuzzy_col,
            self.chebi_all_col,
            self.chebi_final_col,
            self.kegg_exact_col,
            self.kegg_fuzzy_col,
            self.kegg_all_col,
            self.kegg_final_col,
        ]:
            if c not in df2.columns:
                df2[c] = ""
        return df2

    def normalize_joined_ids(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalizes the Step2 final ID columns into canonical, deduped, sorted comma-strings.
        """
        df2 = df.copy()

        if self.chebi_final_col in df2.columns:
            df2[self.chebi_final_col] = df2[self.chebi_final_col].apply(
                lambda v: ",".join(sort_chebi_ids(parse_chebi_list(v))) if str(v).strip() else ""
            )

        if self.kegg_final_col in df2.columns:
            df2[self.kegg_final_col] = df2[self.kegg_final_col].apply(
                lambda v: ",".join(sort_kegg_ids(parse_kegg_list(v))) if str(v).strip() else ""
            )

        return df2

    def output_cols(self) -> List[str]:
        return [
            self.chebi_exact_col,
            self.chebi_fuzzy_col,
            self.chebi_all_col,
            self.chebi_final_col,
            self.kegg_exact_col,
            self.kegg_fuzzy_col,
            self.kegg_all_col,
            self.kegg_final_col,
        ]
