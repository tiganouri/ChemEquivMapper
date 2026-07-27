"""
- fdr_bh: BH-FDR correction
- fisher_exact_right_tail: enrichment Fisher test
- metabolite_sig_map: build metabolite_id -> 0/1 significance mapping from preprocessing stats
"""

from __future__ import annotations

from typing import Iterable, List, Dict, Tuple

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm


def fdr_bh(pvals: Iterable[float], alpha: float = 0.05) -> List[float]:
    """Benjamini-Hochberg FDR correction (statsmodels). Returns q-values aligned to input order."""
    p = list(pvals)
    if len(p) == 0:
        return []
    return list(sm.stats.multipletests(p, alpha=alpha, method="fdr_bh")[1])


def fisher_exact_right_tail(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """
    One-sided Fisher exact test (right tail / enrichment).

    Contingency table:
        in_path   not_in_path
    sel    a          b
    bg     c          d
    """
    if any(x < 0 for x in (a, b, c, d)):
        return 1.0, 1.0
    table = np.array([[a, b], [c, d]], dtype=int)
    odds, p = stats.fisher_exact(table, alternative="greater")
    return float(odds), float(p)


def metabolite_sig_map(
    metabolite_stats: pd.DataFrame,
    *,
    metabolite_id_col: str = "metabolite_id",
    significant_col: str = "significant",
) -> Dict[str, int]:
    """
    Create dict metabolite_id -> 0/1.
    If duplicates exist for the same metabolite_id, uses OR rule (max).
    """
    if metabolite_id_col not in metabolite_stats.columns:
        raise KeyError(f"metabolite_stats missing '{metabolite_id_col}'")
    if significant_col not in metabolite_stats.columns:
        raise KeyError(f"metabolite_stats missing '{significant_col}'")

    tmp = metabolite_stats[[metabolite_id_col, significant_col]].copy()
    tmp[metabolite_id_col] = tmp[metabolite_id_col].astype(str)
    
    tmp[significant_col] = tmp[significant_col].astype(int)

    s = tmp.groupby(metabolite_id_col)[significant_col].max()
    return {k: int(v) for k, v in s.to_dict().items()}
