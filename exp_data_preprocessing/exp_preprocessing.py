from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

# Matplotlib is only needed if you enable plotting
import matplotlib.pyplot as plt


@dataclass
class DatasetPreprocessConfig:
    # filtering / imputation
    max_missing_frac: float = 0.50          # 50% rule (groupwise; keep if ANY group passes)
    min_present_frac: float = 0.50          # (currently not used; kept for backwards compatibility)
    drop_blanks: bool = True
    knn_k: int = 5
    knn_weights: str = "distance"

    na_group_label: str = "NA"              # assign missing group labels to this

    fallback_divisor: float = 1000.0        # fill NaNs with (row_min / fallback_divisor)
    fallback_when_knn_fails: bool = True    # run fallback if KNN throws
    fallback_when_nans_remain: bool = True  # run fallback if KNN completes but NaNs remain

    # log transform (disabled by default; set to "log2" or "log10" to enable)
    log_transform: Optional[str] = None     # {"log2", "log10", None}
    log_pseudocount: float = 1.0            # added before log to avoid log(0)

    # PCA plot controls (plots are optional and can also be controlled per-call)
    pca_n_components: int = 2
    plot_dpi: int = 300


class DatasetPreprocessor:
    """
    Feature-table preprocessing:
      - read tables
      - normalize metadata (sample_id, group)
      - drop blanks (optional)
      - align/subset samples by groups_of_interest
      - groupwise missingness filtering (keep if ANY group passes max_missing_frac)
      - numeric coercion
      - KNN imputation (+ fallback)
      - optional log transform (log2/log10)
      - optional PCA plots (raw, transformed, or side-by-side)
    """

    def __init__(self, cfg: DatasetPreprocessConfig = DatasetPreprocessConfig()):
        self.cfg = cfg

    # ----------------------------
    # IO
    # ----------------------------
    def read_table(self, path: Path, nrows: Optional[int] = None) -> pd.DataFrame:
        suf = path.suffix.lower()
        seps = ["\t"] if suf in [".tsv", ".txt"] else [",", "\t"]
        encodings = ["utf-8", "utf-8-sig", "utf-16", "latin1", "cp1252"]

        last_err: Optional[Exception] = None
        for sep in seps:
            for enc in encodings:
                try:
                    return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False, encoding=enc)
                except Exception as e:
                    last_err = e
        if last_err is None:
            raise RuntimeError("Failed to read table for unknown reasons.")
        raise last_err

    # ----------------------------
    # metadata normalization
    # ----------------------------
    def normalize_metadata(self, meta: pd.DataFrame) -> pd.DataFrame:
        meta = meta.copy()
        cols_lower = {c.lower(): c for c in meta.columns}

        sample_candidates = [
            "sample_id", "sample", "sample name", "samplename",
            "sample_name", "Sample Name", "Sample"
        ]
        sample_col = None
        for k in sample_candidates:
            if k.lower() in cols_lower:
                sample_col = cols_lower[k.lower()]
                break
        if sample_col is None:
            sample_col = meta.columns[0]

        group_col = cols_lower.get("group", None)
        if group_col is None and meta.shape[1] > 1:
            group_col = meta.columns[1]

        rename_map = {sample_col: "sample_id"}
        if group_col is not None:
            rename_map[group_col] = "group"

        meta = meta.rename(columns=rename_map)
        meta["sample_id"] = meta["sample_id"].astype(str).str.strip()

        if "group" in meta.columns:
            # normalize group strings, fill missing with NA label
            meta["group"] = (
                meta["group"]
                .where(~meta["group"].isna(), other=self.cfg.na_group_label)
                .astype(str)
                .str.strip()
            )
            # treat empty/whitespace as NA group too
            meta.loc[meta["group"].str.strip().eq(""), "group"] = self.cfg.na_group_label
        else:
            # if no group col exists, create NA group for all
            meta["group"] = self.cfg.na_group_label

        return meta

    # ----------------------------
    # data cleaning
    # ----------------------------
    def drop_blank_columns(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
        blank_cols = df.columns[df.columns.astype(str).str.contains("blank", case=False, na=False)]
        return df.drop(columns=blank_cols), list(blank_cols)

    def align_and_subset_by_groups(
        self,
        df: pd.DataFrame,
        meta_raw: pd.DataFrame,
        groups_keep: Sequence[str],
        feature_col: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
        df = df.copy()
        meta = self.normalize_metadata(meta_raw)

        if "group" not in meta.columns:
            raise ValueError("Metadata has no 'group' column after normalization.")

        if feature_col is None:
            feature_col = df.columns[0]

        groups_keep = set([str(g).strip() for g in groups_keep])
        meta_sub = meta[meta["group"].isin(groups_keep)].copy()

        # standardize df sample col names
        sample_cols = [c for c in df.columns if c != feature_col]
        df = df.rename(columns={c: str(c).strip() for c in sample_cols})

        wanted = set(meta_sub["sample_id"].tolist())
        overlap = [c for c in df.columns if c != feature_col and c in wanted]

        if len(overlap) == 0:
            raise ValueError(
                "No overlapping sample IDs between df columns and metadata sample_id.\n"
                "Check if df columns are sample IDs."
            )

        df_sub = df[[feature_col] + overlap].copy()

        # align meta to df order
        meta_sub = meta_sub[meta_sub["sample_id"].isin(overlap)].copy()
        meta_sub["sample_id"] = pd.Categorical(meta_sub["sample_id"], categories=overlap, ordered=True)
        meta_sub = meta_sub.sort_values("sample_id").reset_index(drop=True)

        report = {
            "feature_col": feature_col,
            "groups_keep": sorted(list(groups_keep)),
            "n_df_samples_before": len(sample_cols),
            "n_overlap": len(overlap),
            "group_counts": meta_sub["group"].value_counts().to_dict(),
            "df_head_samples": overlap[:10],
        }

        meta_sub = meta_sub.set_index("sample_id")
        return df_sub, meta_sub, report

    def to_numeric_matrix(self, df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
        sample_cols = [c for c in df.columns if c != feature_col]
        X = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        out = pd.concat([df[[feature_col]].astype(str), X], axis=1)
        return out

    # ----------------------------
    # log transform (optional)
    # ----------------------------
    def apply_log_transform(
        self,
        df: pd.DataFrame,
        feature_col: str,
        transform: Optional[str],
        *,
        pseudocount: Optional[float] = None,
    ) -> pd.DataFrame:
        """
        Apply log2 or log10 transform to numeric values.
        Uses pseudocount to avoid log(0).
        """
        if transform is None:
            return df.copy()

        if transform not in {"log2", "log10"}:
            raise ValueError("log_transform must be one of {None, 'log2', 'log10'}")

        if pseudocount is None:
            pseudocount = float(self.cfg.log_pseudocount)

        df = df.copy()
        sample_cols = [c for c in df.columns if c != feature_col]

        X = df[sample_cols].astype(float) + float(pseudocount)

        if transform == "log2":
            X = np.log2(X)
        else:  # log10
            X = np.log10(X)

        df[sample_cols] = X
        return df

    # ----------------------------
    # groupwise missing filter
    # ----------------------------
    def drop_rows_by_missing_fraction_groupwise(
        self,
        df: pd.DataFrame,
        feature_col: str,
        meta_indexed: pd.DataFrame,  # index=sample_id, includes 'group'
        *,
        max_missing_frac: Optional[float] = None,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Keep a feature if ANY group has missing_frac <= max_missing_frac (default cfg.max_missing_frac).

        Missing defined as:
          - NaN after numeric coercion
          - OR empty string/whitespace in raw cells
        """
        if max_missing_frac is None:
            max_missing_frac = self.cfg.max_missing_frac

        if "group" not in meta_indexed.columns:
            raise ValueError("meta_indexed must include a 'group' column (already normalized).")

        df = df.copy()
        sample_cols = [c for c in df.columns if c != feature_col]

        # restrict meta to sample_cols
        meta = meta_indexed.copy()
        meta = meta.loc[meta.index.intersection(sample_cols)].copy()

        # any sample not present in meta? assign NA group
        missing_in_meta = [s for s in sample_cols if s not in meta.index]
        if missing_in_meta:
            add = pd.DataFrame({"group": self.cfg.na_group_label}, index=missing_in_meta)
            meta = pd.concat([meta, add], axis=0)

        groups = meta["group"].astype(str).fillna(self.cfg.na_group_label)
        groups = groups.where(~groups.str.strip().eq(""), other=self.cfg.na_group_label)

        X_raw = df[sample_cols]
        empty_mask = X_raw.astype(str).apply(lambda s: s.str.strip().eq("")).fillna(False)
        X_num = X_raw.apply(pd.to_numeric, errors="coerce")
        nan_mask = X_num.isna()
        missing_mask = nan_mask | empty_mask

        group_levels = list(pd.unique(groups))
        per_group_missing: Dict[str, pd.Series] = {}

        for g in group_levels:
            cols_g = groups.index[groups == g].tolist()
            if not cols_g:
                continue
            miss_frac_g = missing_mask[cols_g].mean(axis=1)
            per_group_missing[str(g)] = miss_frac_g

        keep = pd.Series(False, index=df.index)
        for _, miss_frac_g in per_group_missing.items():
            keep = keep | (miss_frac_g <= float(max_missing_frac))

        out = df.loc[keep].copy()

        report = {
            "rule": "keep if ANY group has missing_frac <= max_missing_frac",
            "max_missing_frac": float(max_missing_frac),
            "rows_before": int(df.shape[0]),
            "rows_after": int(out.shape[0]),
            "rows_dropped": int((~keep).sum()),
            "groups_seen": sorted(list(per_group_missing.keys())),
            "group_sizes": {str(g): int((groups == g).sum()) for g in pd.unique(groups)},
            "per_group_mean_missing_frac": {g: float(s.mean()) for g, s in per_group_missing.items()},
            "per_group_median_missing_frac": {g: float(s.median()) for g, s in per_group_missing.items()},
        }
        return out, report

    def _fallback_impute_min_over_divisor(self, df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
        """
        Fill NaNs per feature(row) with (min_non_missing / divisor).
        If a row has no non-missing values, fill with 0.
        """
        df = df.copy()
        sample_cols = [c for c in df.columns if c != feature_col]

        X = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        row_min = X.min(axis=1, skipna=True)
        fill_vals = (row_min / float(self.cfg.fallback_divisor)).fillna(0.0)

        X_filled = X.copy()
        nan_mask = X_filled.isna()
        X_filled = X_filled.where(~nan_mask, other=np.expand_dims(fill_vals.values, axis=1))

        df[sample_cols] = X_filled
        return df

    def knn_impute_feature_table(
        self,
        df: pd.DataFrame,
        feature_col: str,
        *,
        k: Optional[int] = None,
        weights: Optional[str] = None,
    ) -> pd.DataFrame:
        if k is None:
            k = self.cfg.knn_k
        if weights is None:
            weights = self.cfg.knn_weights

        df = df.copy()
        sample_cols = [c for c in df.columns if c != feature_col]

        X = df[sample_cols].apply(pd.to_numeric, errors="coerce")  # features x samples
        X_t = X.T                                                  # samples x features

        try:
            imputer = KNNImputer(n_neighbors=int(k), weights=str(weights))
            X_imp_t = pd.DataFrame(imputer.fit_transform(X_t), index=X_t.index, columns=X_t.columns)
            df[sample_cols] = X_imp_t.T

            if self.cfg.fallback_when_nans_remain and pd.isna(df[sample_cols]).to_numpy().any():
                df = self._fallback_impute_min_over_divisor(df, feature_col)

            return df

        except Exception:
            if not self.cfg.fallback_when_knn_fails:
                raise
            df[sample_cols] = X
            df = self._fallback_impute_min_over_divisor(df, feature_col)
            return df

    # ----------------------------
    # PCA plotting (optional)
    # ----------------------------
    def _plot_pca_on_axis(
        self,
        ax: plt.Axes,
        df: pd.DataFrame,
        meta_indexed: pd.DataFrame,
        feature_col: str,
        *,
        title: str,
    ) -> None:
        X = df.drop(columns=[feature_col]).T  # samples x features
        pca = PCA(n_components=int(self.cfg.pca_n_components))
        scores = pca.fit_transform(X)

        groups = meta_indexed.loc[X.index, "group"]
        for g in groups.unique():
            idx = (groups == g).to_numpy()
            ax.scatter(scores[idx, 0], scores[idx, 1], label=str(g), alpha=0.8)

        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.legend()

    def save_pca_plot(
        self,
        *,
        df_raw: pd.DataFrame,
        df_log: Optional[pd.DataFrame],
        meta: pd.DataFrame,
        feature_col: str,
        ds_name: str,
        outpath: Path,
        plot_mode: str,
        log_transform: Optional[str],
    ) -> None:
        """
        plot_mode:
          - "none"    : do nothing
          - "raw"     : single PCA on raw
          - "log"     : single PCA on transformed (requires log_transform != None)
          - "compare" : side-by-side raw vs transformed (requires log_transform != None)
        """
        mode = str(plot_mode).lower().strip()
        if mode == "none":
            return

        if mode in {"log", "compare"} and log_transform is None:
            raise ValueError(f"plot_mode='{plot_mode}' requires log_transform to be set (log2/log10).")

        outpath.parent.mkdir(parents=True, exist_ok=True)

        if mode == "raw":
            fig, ax = plt.subplots(1, 1, figsize=(6, 5))
            self._plot_pca_on_axis(ax, df_raw, meta, feature_col, title="Raw")
            fig.suptitle(f"{ds_name} PCA")
            fig.tight_layout()
            fig.savefig(outpath, dpi=int(self.cfg.plot_dpi))
            plt.close(fig)
            return

        if mode == "log":
            assert df_log is not None
            fig, ax = plt.subplots(1, 1, figsize=(6, 5))
            self._plot_pca_on_axis(ax, df_log, meta, feature_col, title=f"{log_transform}")
            fig.suptitle(f"{ds_name} PCA")
            fig.tight_layout()
            fig.savefig(outpath, dpi=int(self.cfg.plot_dpi))
            plt.close(fig)
            return

        if mode == "compare":
            assert df_log is not None
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
            self._plot_pca_on_axis(axes[0], df_raw, meta, feature_col, title="Raw")
            self._plot_pca_on_axis(axes[1], df_log, meta, feature_col, title=f"{log_transform}")
            fig.suptitle(f"{ds_name} PCA (raw vs {log_transform})")
            fig.tight_layout()
            fig.savefig(outpath, dpi=int(self.cfg.plot_dpi))
            plt.close(fig)
            return

        raise ValueError("plot_mode must be one of {'none','raw','log','compare'}")

    # ----------------------------
    # main builder
    # ----------------------------
    def build_preprocessed_datasets(
        self,
        ds_paths: Dict[str, Path],
        group_paths: Dict[str, Path],
        ds_to_meta: Dict[str, str],
        groups_of_interest: Dict[str, Sequence[str]],
        *,
        # Per-call overrides (so you can enable/disable without mutating cfg)
        log_transform: Optional[str] = None,          # None => use cfg.log_transform
        plot_mode: str = "none",                      # "none"|"raw"|"log"|"compare"
        plot_outdir: Optional[Path] = None,           # required if plot_mode != "none"
        plot_filename_suffix: str = "",               # optional
    ) -> Dict[str, Dict[str, Any]]:
        datasets: Dict[str, Dict[str, Any]] = {}

        # resolve transform for this call
        effective_log = self.cfg.log_transform if log_transform is None else log_transform
        if effective_log not in {None, "log2", "log10"}:
            raise ValueError("log_transform must be one of {None, 'log2', 'log10'}")

        for ds_name, ds_path in ds_paths.items():
            meta_key = ds_to_meta[ds_name]
            meta_path = group_paths[meta_key]
            groups_keep = groups_of_interest[meta_key]

            raw_df = self.read_table(ds_path)
            raw_meta = self.read_table(meta_path)

            feature_col = raw_df.columns[0]

            blanks: List[str] = []
            df1 = raw_df
            if self.cfg.drop_blanks:
                df1, blanks = self.drop_blank_columns(df1)

            # align + subset first (needed for groupwise filter)
            df2, meta2, rep_align = self.align_and_subset_by_groups(df1, raw_meta, groups_keep, feature_col=feature_col)

            # groupwise drop sparse rows
            df3, rep_sparse = self.drop_rows_by_missing_fraction_groupwise(df2, feature_col, meta2)

            # numeric + impute
            df3 = self.to_numeric_matrix(df3, feature_col)
            df3 = self.knn_impute_feature_table(df3, feature_col)

            df_raw = df3.copy()
            df_log = self.apply_log_transform(df_raw, feature_col, transform=effective_log)

            # optional plotting
            if str(plot_mode).lower().strip() != "none":
                if plot_outdir is None:
                    raise ValueError("plot_outdir must be provided when plot_mode != 'none'.")

                suffix = f"_{plot_filename_suffix}" if plot_filename_suffix else ""
                outpath = Path(plot_outdir) / f"PCA_{ds_name}{suffix}.png"

                # if no transform configured, df_log is still a copy of df_raw
                # but we enforce semantic correctness: requesting log/compare requires effective_log != None
                self.save_pca_plot(
                    df_raw=df_raw,
                    df_log=(df_log if effective_log is not None else None),
                    meta=meta2,
                    feature_col=feature_col,
                    ds_name=ds_name,
                    outpath=outpath,
                    plot_mode=plot_mode,
                    log_transform=effective_log,
                )

            datasets[ds_name] = {
                "paths": {"data": ds_path, "meta": meta_path, "meta_key": meta_key},
                "groups_keep": list(groups_keep),
                "feature_col": feature_col,
                "raw_df": raw_df,
                "raw_meta": raw_meta,
                "df_raw": df_raw,
                "df_log": (df_log if effective_log is not None else None),
                "meta": meta2,
                "log_transform": effective_log,
                "report_align": rep_align,
                "report_blanks": {"dropped_blank_cols": blanks, "n_dropped": len(blanks)},
                "report_sparse_rows": rep_sparse,
            }

        return datasets
