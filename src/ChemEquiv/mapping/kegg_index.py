from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd


@dataclass
class KEGGIndex:
    """
    In-memory index: KEGG compound ID -> list of KEGG pathway IDs

    Compound IDs:
      - stored normalized as "C00031" (no "KEGG:" prefix)
      - get_pathways() normalizes too so callers can pass "KEGG:C00031" or "C00031"

    Pathway IDs:
      - whatever is in your TSV (typically "map00010", "map00020", etc.)
    """
    compound_to_pathways: Dict[str, List[str]]

    @staticmethod
    def _normalize_compound_id(x: str) -> str:
        x = (x or "").strip()
        if not x:
            return ""
        # Accept KEGG:C00031 or C00031 or cpd:C00031
        if ":" in x:
            x = x.split(":", 1)[1].strip()
        return x

    def get_pathways(self, kegg_compound_id: str) -> List[str]:
        kid = KEGGIndex._normalize_compound_id(kegg_compound_id)
        return self.compound_to_pathways.get(kid, [])

    @staticmethod
    def from_compound_pathway_tsv(tsv_path: Path, assume_header: bool = True) -> "KEGGIndex":
        """
        Build compound -> pathways from a TSV with two columns:
          kegg_compound_id   kegg_pathway_id
        """
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"KEGG TSV not found: {tsv_path}")

        df = pd.read_csv(tsv_path, sep="\t", header=0 if assume_header else None, dtype=str)
        df = df.fillna("")

        # Determine columns
        if assume_header:
            c_comp = df.columns[0]
            c_path = df.columns[1] if len(df.columns) > 1 else df.columns[0]
        else:
            c_comp = df.columns[0]
            c_path = df.columns[1]

        comp_to_paths: Dict[str, List[str]] = {}
        for comp, pw in zip(df[c_comp].astype(str), df[c_path].astype(str)):
            comp = KEGGIndex._normalize_compound_id(comp)
            pw = (pw or "").strip()
            if not comp or not pw:
                continue
            comp_to_paths.setdefault(comp, []).append(pw)

        return KEGGIndex(comp_to_paths)
