from __future__ import annotations

from dataclasses import dataclass
import pandas as pd

from ChemEquiv.steps.step2_chebi.supplementary import (
    parse_chebi_list,
    sort_chebi_ids,
)


@dataclass
class Step5Schema:
    # outputs
    chebi_generalisation_col: str = "CHEBI_IsA_Generalisation"
    reactome_generalisation_col: str = "Reactome_Pathways_IsA_Generalisation"
    n_pathways_col: str = "n_pathways_isa_generalisation"

    def ensure_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        if self.chebi_generalisation_col not in out.columns:
            out[self.chebi_generalisation_col] = ""
        if self.reactome_generalisation_col not in out.columns:
            out[self.reactome_generalisation_col] = ""
        if self.n_pathways_col not in out.columns:
            out[self.n_pathways_col] = 0
        return out

    def normalize(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out[self.chebi_generalisation_col] = out[self.chebi_generalisation_col].apply(
            lambda x: ",".join(sort_chebi_ids(parse_chebi_list(x)))
        )
        # Reactome pathways are strings like "R-HSA-..." — keep unique + sorted
        out[self.reactome_generalisation_col] = out[self.reactome_generalisation_col].apply(
            lambda x: ",".join(sorted({p.strip() for p in str(x).split(",") if p.strip()}))
        )
        # n_pathways int
        out[self.n_pathways_col] = out[self.reactome_generalisation_col].apply(
            lambda x: 0 if not str(x).strip() else len([p for p in str(x).split(",") if p.strip()])
        )
        return out
