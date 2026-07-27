from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Union, Any

from .reactome_index import ReactomeIndex
from .kegg_index import KEGGIndex


def _split_csv_cell(cell: str) -> List[str]:
    if cell is None:
        return []
    s = str(cell).strip()
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def _normalize_chebi(x: str) -> str:
    x = (x or "").strip()
    if not x:
        return ""
    # Accept "89605" or "CHEBI:89605"
    if x.isdigit():
        return f"CHEBI:{x}"
    if x.upper().startswith("CHEBI:"):
        return "CHEBI:" + x.split(":", 1)[1].strip()
    return x


def _normalize_kegg_compound(x: str) -> str:
    x = (x or "").strip()
    if not x:
        return ""
    # Accept KEGG:C00031 or C00031
    if x.startswith("KEGG:"):
        x = x.split(":", 1)[1].strip()
    return x


@dataclass
class MappingResult:
    reactome_chebis: List[str]
    reactome_pathways: List[str]
    kegg_compounds: List[str]
    kegg_pathways: List[str]


class PathwayMapper:
    """
    Reusable mapper for multiple steps.

    - Reactome mapping uses CHEBI -> pathways.
    - KEGG mapping uses KEGG compound -> pathways, with optional CHEBI -> KEGG compound expansion.

    Key defaults:
      - dedup_pathways=False (preserve duplicates; matches your requirement)
      - dedup_compounds=True (avoid repeated lookups; faster)
    """

    def __init__(
        self,
        reactome_index: Optional[ReactomeIndex] = None,
        kegg_index: Optional[KEGGIndex] = None,
        chebi_to_kegg: Optional[Dict[str, List[str]]] = None,
    ):
        self.reactome_index = reactome_index
        self.kegg_index = kegg_index
        self.chebi_to_kegg = chebi_to_kegg or {}

    def map(
        self,
        chebi_ids: Optional[Union[Iterable[str], str]] = None,
        kegg_compound_ids: Optional[Union[Iterable[str], str]] = None,
        *,
        dedup_pathways: bool = False,
        dedup_compounds: bool = True,
    ) -> MappingResult:
        # Allow CSV strings directly (common in step output cells)
        if isinstance(chebi_ids, str):
            chebi_ids = _split_csv_cell(chebi_ids)
        if isinstance(kegg_compound_ids, str):
            kegg_compound_ids = _split_csv_cell(kegg_compound_ids)

        # Normalize inputs
        chebis: List[str] = []
        if chebi_ids:
            for c in chebi_ids:
                nc = _normalize_chebi(str(c))
                if nc:
                    chebis.append(nc)

        keggs: List[str] = []
        if kegg_compound_ids:
            for k in kegg_compound_ids:
                nk = _normalize_kegg_compound(str(k))
                if nk:
                    keggs.append(nk)

        # Expand KEGG from CHEBI if available
        if self.chebi_to_kegg and chebis:
            for c in chebis:
                for k in self.chebi_to_kegg.get(c, []):
                    nk = _normalize_kegg_compound(k)
                    if nk:
                        keggs.append(nk)

        # Dedup compounds (keeps stable order)
        if dedup_compounds:
            seen = set()
            chebis = [x for x in chebis if not (x in seen or seen.add(x))]
            seen = set()
            keggs = [x for x in keggs if not (x in seen or seen.add(x))]

        # Reactome mapping
        reactome_chebis: List[str] = []
        reactome_pathways: List[str] = []
        if self.reactome_index and chebis:
            for c in chebis:
                pws = self.reactome_index.get_pathways(c)
                if pws:
                    reactome_chebis.append(c)
                    # keep duplicates across chebis if same pathway appears multiple times
                    reactome_pathways.extend(pws)

        # KEGG mapping
        kegg_compounds: List[str] = []
        kegg_pathways: List[str] = []
        if self.kegg_index and keggs:
            for k in keggs:
                pws = self.kegg_index.get_pathways(k)
                if pws:
                    kegg_compounds.append(k)
                    kegg_pathways.extend(pws)

        # Optional dedup of pathways (OFF by default)
        if dedup_pathways:
            seen = set()
            reactome_pathways = [x for x in reactome_pathways if not (x in seen or seen.add(x))]
            seen = set()
            kegg_pathways = [x for x in kegg_pathways if not (x in seen or seen.add(x))]

        return MappingResult(
            reactome_chebis=reactome_chebis,
            reactome_pathways=reactome_pathways,
            kegg_compounds=kegg_compounds,
            kegg_pathways=kegg_pathways,
        )

    def apply_to_df(
        self,
        df: Any,
        idx: Any,
        *,
        step: int,
        chebi_ids: Optional[Union[Iterable[str], str]] = None,
        kegg_compound_ids: Optional[Union[Iterable[str], str]] = None,
        dedup_pathways: bool = False,
        dedup_compounds: bool = True,
    ) -> MappingResult:
        """
        Convenience method for steps: perform mapping and write step-specific columns.

        Columns written:
          - ChEBI_ID_Step{step}_reactome
          - Reactome_Pathways_Step{step}
          - KEGG_Compounds_Step{step}
          - KEGG_Pathways_Step{step}
        """
        res = self.map(
            chebi_ids=chebi_ids,
            kegg_compound_ids=kegg_compound_ids,
            dedup_pathways=dedup_pathways,
            dedup_compounds=dedup_compounds,
        )

        df.at[idx, f"ChEBI_ID_Step{step}_reactome"] = ",".join(res.reactome_chebis) if res.reactome_chebis else ""
        df.at[idx, f"Reactome_Pathways_Step{step}"] = ",".join(res.reactome_pathways) if res.reactome_pathways else ""

        df.at[idx, f"KEGG_Compounds_Step{step}"] = ",".join(res.kegg_compounds) if res.kegg_compounds else ""
        df.at[idx, f"KEGG_Pathways_Step{step}"] = ",".join(res.kegg_pathways) if res.kegg_pathways else ""

        return res
