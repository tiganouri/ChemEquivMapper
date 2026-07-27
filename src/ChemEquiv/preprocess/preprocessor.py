from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.stats as sp_stats
import statsmodels.api as sm

from .input_table import InputTableConfig, InputTablePreprocessor


# ----------------------------
# Stats helpers (BH-FDR + Fisher)
# ----------------------------

def _fdr_bh(pvals: List[float]) -> List[float]:
    p = np.array([float(x) if np.isfinite(x) else 1.0 for x in pvals], dtype=float)
    p = np.clip(p, 0.0, 1.0)
    return sm.stats.multipletests(p, alpha=0.05, method="fdr_bh")[1].tolist()


def _combine_pvalues_fisher(pvals: List[float]) -> float:
    p = np.array([float(x) for x in pvals if np.isfinite(x)], dtype=float)
    if p.size == 0:
        return 1.0
    p = np.clip(p, 0.0, 1.0)
    return float(sp_stats.combine_pvalues(p, method="fisher")[1])


def _median_nonzero(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    x = x[x > 0]
    if x.size == 0:
        return np.nan
    return float(np.median(x))


def _feature_pvalue(
    x: np.ndarray,
    mask0: np.ndarray,
    mask1: np.ndarray,
    *,
    test: Literal["mwu", "ttest"] = "ttest",
    ttest_equal_var: bool = False,
    min_n: int = 2,
) -> float:
    a = x[mask0]
    b = x[mask1]
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < min_n or len(b) < min_n:
        return 1.0

    if test == "mwu":
        return float(sp_stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)

    return float(sp_stats.ttest_ind(a, b, equal_var=ttest_equal_var, nan_policy="omit").pvalue)


# ----------------------------
# Config + Result
# ----------------------------

@dataclass(frozen=True)
class PreprocessConfig:
    """
    Whole-preprocess configuration.

    ratio_thresh:
      Scale similarity threshold: max(median_nonzero)/min(median_nonzero) <= ratio_thresh => similar

    test:
      "ttest" for 2-sample t-test (recommended default),"mwu" for Mann–Whitney  

    alpha:
      significance threshold applied after BH-FDR

    feature_col:
      feature name column in the input table; if None uses the first column
    """
    feature_col: Optional[str] = None

    ratio_thresh: float = 10.0

    test: Literal["mwu", "ttest"] = "ttest"
    ttest_equal_var: bool = False
    alpha: float = 0.05
    min_n_per_group: int = 2


@dataclass
class PreprocessResult:
    """
    Outputs consumed by Steps + Analysis.

    abundance_feature:
      features x samples (index=feature_id, columns=sample_id)

    feature_meta:
      indexed by feature_id with:
        - feature_id (unique)
        - feature_name (canonicalised; duplicates allowed)
        - row_number (optional)

    metabolite_stats:
      indexed by metabolite_id (= canonical feature_name grouping)
      with p_value, p_adj, significant, n_features, method_used

    scale_report:
      metabolite_id, n_features, median_min, median_max, median_ratio, log10_median_ratio, scales_similar

    metabolite_feature_map:
      metabolite_id -> list of feature_id that contributed
    """
    abundance_feature: pd.DataFrame
    feature_meta: pd.DataFrame
    metabolite_stats: pd.DataFrame
    scale_report: pd.DataFrame
    metabolite_feature_map: pd.DataFrame


# ----------------------------
# Main Preprocessor
# ----------------------------

class Preprocessor:
    """
    One entrypoint that prepares a consistent foundation for:
      - Steps (feature-level)
      - ORA (metabolite-level significance)

    Assumptions:
      - Input feature table is features as rows, samples as columns
      - sample_meta index contains sample IDs matching the sample columns
      - group_col provides 2-group labels (required for current hybrid p-value strategy)
    """

    def __init__(self, cfg: PreprocessConfig = PreprocessConfig()):
        self.cfg = cfg

    def _align_samples(
        self,
        abundance_feature: pd.DataFrame,
        sample_meta: pd.DataFrame,
        *,
        group_col: str,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, List[str]]:
        # sample IDs are abundance_feature columns
        sample_cols = [str(c) for c in abundance_feature.columns]
        if sample_meta.index.name is None:
            # ok; but we require index to be sample ids
            pass

        meta_index = [str(i) for i in sample_meta.index.tolist()]

        common = [c for c in sample_cols if c in set(meta_index)]
        if len(common) == 0:
            raise ValueError(
                "No overlapping sample IDs between feature table columns and sample_meta.index. "
                "Make sure sample_meta.index are the sample IDs."
            )

        dropped_table = [c for c in sample_cols if c not in set(common)]
        dropped_meta = [i for i in meta_index if i not in set(common)]

        if dropped_table:
            warnings.warn(f"Dropping {len(dropped_table)} sample columns not present in metadata. Example: {dropped_table[:5]}")
        if dropped_meta:
            warnings.warn(f"Ignoring {len(dropped_meta)} metadata rows not present in feature table. Example: {dropped_meta[:5]}")

        abundance_feature = abundance_feature[common].copy()
        sample_meta = sample_meta.loc[common].copy()

        if group_col not in sample_meta.columns:
            raise ValueError(f"group_col='{group_col}' not found in sample_meta columns: {sample_meta.columns.tolist()}")

        if sample_meta[group_col].isna().any():
            bad = sample_meta.index[sample_meta[group_col].isna()].tolist()
            raise ValueError(f"Missing group labels for {len(bad)} samples. Example: {bad[:10]}")

        groups = list(pd.unique(sample_meta[group_col].astype(str)))
        if len(groups) != 2:
            raise ValueError(f"Hybrid preprocess currently requires exactly 2 groups; got {len(groups)}: {groups}")

        g0, g1 = groups[0], groups[1]
        mask0 = (sample_meta[group_col].astype(str) == g0).to_numpy()
        mask1 = (sample_meta[group_col].astype(str) == g1).to_numpy()

        return abundance_feature, sample_meta, mask0, mask1, [g0, g1]

    def run(
        self,
        feature_table_path: str | pd.PathLike,
        sample_meta: pd.DataFrame,
        *,
        group_col: str,
        metabolite_feature_map: Optional[pd.DataFrame] = None,
    ) -> PreprocessResult:
        
        # 1) Read + canonicalize + ensure unique feature_id
        reader = InputTablePreprocessor(InputTableConfig(feature_col=self.cfg.feature_col))
        abundance_feature, feature_meta = reader.read(feature_table_path)

        # 2) Align samples with metadata
        abundance_feature, sample_meta2, mask0, mask1, group_names = self._align_samples(
            abundance_feature, sample_meta, group_col=group_col
        )

        # 3) Build metabolite groups by canonical feature_name
        # metabolite_id := canonicalized name (feature_meta["feature_name"])
        # build mapping: metabolite_id -> list(feature_id)
        metab_to_features: Dict[str, List[str]] = {}

        if metabolite_feature_map is not None:
            # expected columns: metabolite_id, feature_ids_used (list) OR metabolite_id, feature_id rows
            m = metabolite_feature_map.copy()
        
            if "metabolite_id" not in m.columns:
                raise ValueError("metabolite_feature_map must contain a 'metabolite_id' column.")
        
            if "feature_ids_used" in m.columns:
                for _, r in m.iterrows():
                    metab = str(r["metabolite_id"])
                    fids = r["feature_ids_used"]
                    if isinstance(fids, str):
                        # if stored as stringified list, try safe parse fallback
                        # (optional: you can skip parsing if you guarantee lists)
                        fids = [fids]
                    metab_to_features[metab] = [str(x) for x in fids]
            elif "feature_id" in m.columns:
                for metab, sub in m.groupby("metabolite_id"):
                    metab_to_features[str(metab)] = [str(x) for x in sub["feature_id"].tolist()]
            else:
                raise ValueError("metabolite_feature_map must contain 'feature_ids_used' OR 'feature_id' column.")
        else:
            # fallback to group by feature_name
            for fid, row in feature_meta.iterrows():
                metab = str(row["feature_name"])
                metab_to_features.setdefault(metab, []).append(str(row["feature_id"]))

        # 4) Compute scale_report + hybrid p-values
        # Use X = samples x features for stats
        X = abundance_feature.T  # samples x features

        scale_rows = []
        stat_rows = []
        map_rows = []

        for metab, fids in metab_to_features.items():
            # gather feature vectors
            feats = [f for f in fids if f in X.columns]
            if len(feats) == 0:
                continue

            sub = X[feats].to_numpy(dtype=float)  # n_samples x n_feats

            # drop metabolite if all NaN or all zeros across all features
            if np.all(~np.isfinite(sub)):
                continue
            if np.nanmax(sub) <= 0:
                continue

            # medians for scale check
            meds = np.array([_median_nonzero(sub[:, j]) for j in range(sub.shape[1])], dtype=float)
            meds = meds[np.isfinite(meds)]

            if meds.size == 0:
                median_min = median_max = median_ratio = log10_ratio = np.nan
                scales_similar = False
            else:
                median_min = float(np.min(meds))
                median_max = float(np.max(meds))
                if median_min <= 0:
                    median_ratio = float("inf")
                else:
                    median_ratio = float(median_max / median_min)
                log10_ratio = float(np.log10(median_ratio)) if np.isfinite(median_ratio) and median_ratio > 0 else np.nan
                scales_similar = bool(np.isfinite(median_ratio) and median_ratio <= self.cfg.ratio_thresh)

            scale_rows.append(
                dict(
                    metabolite_id=metab,
                    n_features=int(len(feats)),
                    median_min=median_min,
                    median_max=median_max,
                    median_ratio=median_ratio,
                    log10_median_ratio=log10_ratio,
                    scales_similar=scales_similar,
                )
            )

            map_rows.append({"metabolite_id": metab, "feature_ids_used": feats})

            # Hybrid p-value logic
            if len(feats) == 1:
                p = _feature_pvalue(
                    sub[:, 0],
                    mask0,
                    mask1,
                    test=self.cfg.test,
                    ttest_equal_var=self.cfg.ttest_equal_var,
                    min_n=self.cfg.min_n_per_group,
                )
                method = "single"
            else:
                if scales_similar:
                    agg = np.nanmean(sub, axis=1)
                    p = _feature_pvalue(
                        agg,
                        mask0,
                        mask1,
                        test=self.cfg.test,
                        ttest_equal_var=self.cfg.ttest_equal_var,
                        min_n=self.cfg.min_n_per_group,
                    )
                    method = "aggregate"
                else:
                    pvals = [
                        _feature_pvalue(
                            sub[:, j],
                            mask0,
                            mask1,
                            test=self.cfg.test,
                            ttest_equal_var=self.cfg.ttest_equal_var,
                            min_n=self.cfg.min_n_per_group,
                        )
                        for j in range(sub.shape[1])
                    ]
                    p = _combine_pvalues_fisher(pvals)
                    method = "fisher"

            stat_rows.append(
                dict(
                    metabolite_id=metab,
                    p_value=float(p) if np.isfinite(p) else 1.0,
                    method_used=method,
                    n_features=int(len(feats)),
                    n_per_group=f"{group_names[0]}:{int(mask0.sum())},{group_names[1]}:{int(mask1.sum())}",
                )
            )

        scale_report = pd.DataFrame(scale_rows).sort_values(
            ["scales_similar", "median_ratio"], ascending=[False, True]
        ).reset_index(drop=True)

        metabolite_stats = pd.DataFrame(stat_rows)
        if metabolite_stats.empty:
            metabolite_stats = pd.DataFrame(
                columns=["metabolite_id", "p_value", "p_adj", "significant", "method_used", "n_features", "n_per_group"]
            )
        else:
            metabolite_stats["p_adj"] = _fdr_bh(metabolite_stats["p_value"].tolist())
            metabolite_stats["significant"] = metabolite_stats["p_adj"] <= self.cfg.alpha
            metabolite_stats = metabolite_stats.sort_values(["p_adj", "p_value"], ascending=[True, True]).reset_index(drop=True)

        metabolite_feature_map = pd.DataFrame(map_rows)

        return PreprocessResult(
            abundance_feature=abundance_feature,
            feature_meta=feature_meta,
            metabolite_stats=metabolite_stats,
            scale_report=scale_report,
            metabolite_feature_map=metabolite_feature_map,
        )
