from __future__ import annotations

import re
from typing import List, Set


# -----------------------------
# Normalization helpers
# -----------------------------

def norm(s: object) -> Optional[str]:
    """
    Case-insensitive normalizer for names, tolerating None.

    - Strips leading/trailing whitespace
    - Normalizes internal whitespace to a single space
    - Returns casefolded string
    """
    if not isinstance(s, str):
        return None
    normalized = " ".join(s.strip().split())
    return normalized.casefold() if normalized else None


def norm_remove_special(s: object) -> Optional[str]:
    """
    Trim and collapse whitespace
    Additional normalizer that replaces parentheses, brackets, and asterisks with spaces.
    Collapse whitespace again
    Lowercase

    Used for fuzzy matching when exact match fails.

    This allows matching:
        'PE(16:0/18:1)'  -> 'pe 16:0/18:1'
        'PE 16:0/18:1'   -> 'pe 16:0/18:1'
    """
    if not isinstance(s, str):
        return None

    text = " ".join(s.strip().split())
    text = re.sub(r"[(){}\[\]*]", " ", text)
    text = " ".join(text.split())
    return text.casefold() if text else None

# -----------------------------
# ID parsing + sorting
# -----------------------------
# NOTE:
# Step1/RefMet may return ChEBI as raw numbers (e.g., "89605") or "-"
# Step2 must accept both "CHEBI:89605" and "89605" and normalize to "CHEBI:89605".

_CHEBI_CANON_RE = re.compile(r"(CHEBI:\d+)", re.IGNORECASE)
_CHEBI_NUM_RE = re.compile(r"\b(\d{3,})\b")  # ChEBI numeric IDs are typically >= 3 digits

_KEGG_CPD_RE = re.compile(r"\b(C\d{5})\b", re.IGNORECASE)


def parse_chebi_list(val: object) -> Set[str]:
    """
    Parse ChEBI IDs from a value that may contain:
      - "CHEBI:89605"
      - "89605"
      - "89605, 73704" or "89605;73704" or "89605|73704"
      - "-" / empty / NaN

    Returns canonical uppercase IDs: {"CHEBI:89605", ...}
    """
    if val is None:
        return set()

    s = str(val).strip()
    if not s or s == "-" or s.lower() == "nan":
        return set()

    # unify separators
    s2 = s.replace(";", ",").replace("|", ",")
    out: Set[str] = set()

    # 1) capture canonical CHEBI:xxxxx if present
    for h in _CHEBI_CANON_RE.findall(s2):
        if h:
            out.add(h.upper())

    # 2) if no canonical IDs found, accept raw numbers as ChEBI IDs
    if not out:
        for n in _CHEBI_NUM_RE.findall(s2):
            if n:
                out.add(f"CHEBI:{n}")

    return out


def parse_kegg_list(val: object) -> Set[str]:
    """
    Parse KEGG compound IDs (C00031 etc) from a value that may contain:
      - "C00031"
      - "KEGG:C00031"
      - "C00031, C00022" etc
      - "-" / empty / NaN
    Returns uppercase IDs: {"C00031", ...}
    """
    if val is None:
        return set()

    s = str(val).strip()
    if not s or s == "-" or s.lower() == "nan":
        return set()

    hits = _KEGG_CPD_RE.findall(s.upper())
    return {h.upper() for h in hits if h}


def sort_chebi_ids(ids: Set[str] | List[str]) -> List[str]:
    """Sort CHEBI:xxxxx by numeric part."""
    def key_fn(x: str) -> int:
        m = re.search(r"CHEBI:(\d+)", str(x).upper())
        return int(m.group(1)) if m else 10**18

    return sorted({str(i).upper() for i in ids if i}, key=key_fn)


def sort_kegg_ids(ids: Set[str] | List[str]) -> List[str]:
    """Sort Cxxxxx by numeric part."""
    def key_fn(x: str) -> int:
        m = re.search(r"C(\d{5})", str(x).upper())
        return int(m.group(1)) if m else 10**18

    return sorted({str(i).upper() for i in ids if i}, key=key_fn)
