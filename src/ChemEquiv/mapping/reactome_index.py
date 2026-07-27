from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd


@dataclass
class ReactomeIndex:
    """
    In-memory index: CHEBI:<id> -> list of Reactome pathway IDs.

    Notes:
    - Preserves duplicates from the source file.
    - Supports Reactome ChEBI2Reactome_All_Levels.txt (headerless TSV).
    """
    chebi_to_pathways: Dict[str, List[str]]

    def get_pathways(self, chebi_id: str) -> List[str]:
        """Return pathways for a CHEBI id. Accepts either '10033' or 'CHEBI:10033'."""
        chebi_id = str(chebi_id).strip()
        if chebi_id.isdigit():
            chebi_id = f"CHEBI:{chebi_id}"
        elif chebi_id.upper().startswith("CHEBI:"):
            chebi_id = "CHEBI:" + chebi_id.split(":", 1)[1]
        return self.chebi_to_pathways.get(chebi_id, [])

    @staticmethod
    def _normalize_chebi(x: str) -> str:
        x = str(x).strip()
        if not x:
            return ""
        if x.isdigit():
            return f"CHEBI:{x}"
        if x.upper().startswith("CHEBI:"):
            return "CHEBI:" + x.split(":", 1)[1]
        return x

    @staticmethod
    def from_reactome_tsv(
        tsv_path: Path,
        species: str = "Homo sapiens",
        assume_header: bool = False,  # Reactome TXT is usually headerless
    ) -> "ReactomeIndex":
        """
        Build CHEBI -> pathways index from a Reactome mapping TXT/TSV.

        Reactome ChEBI2Reactome_All_Levels.txt is typically a headerless TSV with columns:
          0: ChEBI numeric id
          1: Reactome pathway stable id (e.g., R-HSA-..., R-BTA-...)
          2: URL
          3: Pathway name
          4: Evidence
          5: Species (e.g., Homo sapiens)

        If assume_header=True, will attempt to find named columns and fall back to heuristics.
        """
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Reactome mapping not found: {tsv_path}")

        df = pd.read_csv(tsv_path, sep="\t", header=0 if assume_header else None, dtype=str).fillna("")

        # Fast path: headerless, fixed-position Reactome TXT
        if not assume_header:
            if df.shape[1] < 2:
                raise ValueError(f"Reactome mapping has too few columns: {tsv_path} (cols={df.shape[1]})")

            col_chebi, col_pathway = 0, 1
            col_species = 5 if df.shape[1] > 5 else None

            if species and col_species is not None:
                df = df[df[col_species].astype(str).str.strip() == species]

            chebi_to_paths: Dict[str, List[str]] = {}
            for chebi, pw in zip(df[col_chebi].astype(str), df[col_pathway].astype(str)):
                chebi = ReactomeIndex._normalize_chebi(chebi)
                pw = str(pw).strip()
                if not chebi or not pw:
                    continue
                chebi_to_paths.setdefault(chebi, []).append(pw)

            return ReactomeIndex(chebi_to_paths)

        # Header path: try named columns, then heuristics
        col_chebi = None
        col_pathway = None
        col_species = None

        lower = {str(c).lower(): c for c in df.columns}

        # Named columns first
        for cand in ["chebi", "chebi_id", "chebi identifier", "chebiid"]:
            if cand in lower:
                col_chebi = lower[cand]
                break

        for cand in ["pathway", "pathway_id", "reactome", "reactome_id", "pathway identifier"]:
            if cand in lower:
                col_pathway = lower[cand]
                break

        for cand in ["species", "organism"]:
            if cand in lower:
                col_species = lower[cand]
                break

        # Heuristics fallback
        if col_chebi is None:
            for c in df.columns:
                s = df[c].astype(str)
                # accept either CHEBI: prefix or pure numeric IDs
                if (s.str.contains("CHEBI:", regex=False).mean() > 0.05) or (s.str.match(r"^\d+$").mean() > 0.2):
                    col_chebi = c
                    break

        if col_pathway is None:
            for c in df.columns:
                s = df[c].astype(str)
                if (s.str.startswith("R-", na=False).mean() > 0.2):
                    col_pathway = c
                    break

        if col_species is None and species:
            for c in df.columns:
                s = df[c].astype(str)
                if (s == species).mean() > 0.2:
                    col_species = c
                    break

        if col_chebi is None or col_pathway is None:
            raise ValueError(
                f"Could not infer Reactome CHEBI/pathway columns from: {tsv_path}. "
                f"Found chebi={col_chebi}, pathway={col_pathway}, species={col_species}"
            )

        if species and col_species is not None:
            df = df[df[col_species].astype(str).str.strip() == species]

        chebi_to_paths: Dict[str, List[str]] = {}
        for chebi, pw in zip(df[col_chebi].astype(str), df[col_pathway].astype(str)):
            chebi = ReactomeIndex._normalize_chebi(chebi)
            pw = str(pw).strip()
            if not chebi or not pw:
                continue
            chebi_to_paths.setdefault(chebi, []).append(pw)

        return ReactomeIndex(chebi_to_paths)
