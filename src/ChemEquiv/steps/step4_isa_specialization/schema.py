from __future__ import annotations

from dataclasses import dataclass
import pandas as pd

from GEPC.steps.step2_chebi.supplementary import (
    parse_chebi_list,
    parse_kegg_list,
    sort_chebi_ids,
    sort_kegg_ids,
)


@dataclass
class Step4Schema:
    """
    Step4 schema only owns Step4 expansion outputs.
    Mapping outputs (Reactome/KEGG pathways) are written by PathwayMapper.

    Expected input cols:
      - chebi_step3_col (default: CHEBI_ID_Step3)
      - kegg_step3_col  (default: KEGG_ID_Step3)
    """
    chebi_step3_col: str = "CHEBI_ID_Step3"
    kegg_step3_col: str = "KEGG_ID_Step3"

    # step4 outputs
    chebi_isa_equiv_col: str = "CHEBI_is_a_Equivalents"
    chebi_step4_col: str = "CHEBI_ID_Step4"
    kegg_step4_col: str = "KEGG_ID_Step4"

    def ensure_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for col in [self.chebi_isa_equiv_col, self.chebi_step4_col, self.kegg_step4_col]:
            if col not in out.columns:
                out[col] = ""
        return out

    def normalize_joined_ids(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize CHEBI/KEGG CSV cells for the Step4 union columns:
          - parse
          - sort
          - re-join
        Does NOT touch mapper-generated pathway columns.
        """
        out = df.copy()

        # CHEBI_is_a_Equivalents
        out[self.chebi_isa_equiv_col] = out[self.chebi_isa_equiv_col].apply(
            lambda x: ",".join(sort_chebi_ids(parse_chebi_list(x)))
        )

        # CHEBI_ID_Step4
        out[self.chebi_step4_col] = out[self.chebi_step4_col].apply(
            lambda x: ",".join(sort_chebi_ids(parse_chebi_list(x)))
        )

        # KEGG_ID_Step4
        out[self.kegg_step4_col] = out[self.kegg_step4_col].apply(
            lambda x: ",".join(sort_kegg_ids(parse_kegg_list(x)))
        )

        return out
