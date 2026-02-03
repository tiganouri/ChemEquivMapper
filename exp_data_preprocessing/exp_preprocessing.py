from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Any

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


@dataclass
class DatasetPreprocessConfig:
    max_missing_frac: float = 0.50          # 50% rule
    min_present_frac: float = 0.50          # keep if ANY group has >= this present (kept for compatibility)
    drop_blanks: bool = True
    knn_k: int = 5
    knn_weights: str = "distance"

    na_group_label: str = "NA"              # assign missing group labels to this

    fallback_divisor: float = 1000.0        # fill NaNs with (row_min / fallback_divisor)
    fallback_when_knn_fails: bool = True    # run fallback if KNN throws
    fallback_when_nans_remain: bool = True  # run fallback if KNN completes but NaNs remain


@dataclass
class LogTransformConfig:
    """
    Only controls log transform and optional PCA compare plot.

    log_mode:
      - "none" | "log2" | "log10"

    If pca_compare_out_dir is not None AND log_mode != "none",
    a 2-panel PCA figure is saved:
      left = PCA on RAW
      right = PCA on LOG
    """
    log_mode: str = "none"                  # "none" | "log2" | "log10"
    log_offset: float = 1e-9                # avoid log(0)

    # PCA compare output (optional)
    pca_compare_out_dir: Optional[Path] = None
    pca_n_components: int = 2
    pca_dpi: int = 220
    pca_prefix: str = "PCA_COMPARE"

    pca_loadings_out_dir: Optional[Path] = None
    pca_loadings_prefix: str = "PCA_LOADINGS"
    
class DatasetPreprocessor:
    def __init__(
        self,
        cfg: DatasetPreprocessConfig = DatasetPreprocessConfig(),
        log_cfg: LogTransformConfig = LogTransformConfig(),
    ):
        self.cfg = cfg
        self.log_cfg = log_cfg

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
            meta["group"] = (
                meta["group"]
                .astype(str)
                .where(~meta["group"].isna(), other=self.cfg.na_group_label)
            )
            meta["group"] = meta["group"].astype(str).str.strip()
            meta.loc[meta["group"].str.strip().eq(""), "group"] = self.cfg.na_group_label
        else:
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
        Keep a feature if at least one group has missing_frac <= max_missing_frac (default 0.50).

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

        # group labels
        groups = meta["group"].astype(str).fillna(self.cfg.na_group_label)
        groups = groups.where(~groups.str.strip().eq(""), other=self.cfg.na_group_label)

        X_raw = df[sample_cols]

        empty_mask = X_raw.astype(str).apply(lambda s: s.str.strip().eq("")).fillna(False)
        X_num = X_raw.apply(pd.to_numeric, errors="coerce")
        nan_mask = X_num.isna()
        missing_mask = nan_mask | empty_mask  # features x samples boolean

        # compute per-group missing fractions
        group_levels = list(pd.unique(groups))
        per_group_missing: Dict[str, pd.Series] = {}

        for g in group_levels:
            cols_g = groups.index[groups == g].tolist()
            if not cols_g:
                continue
            miss_frac_g = missing_mask[cols_g].mean(axis=1)
            per_group_missing[str(g)] = miss_frac_g

        # keep if ANY group passes the threshold
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

    # ----------------------------
    # imputation
    # ----------------------------
    def _fallback_impute_min_over_divisor(self, df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
        """
        Fill NaNs per feature(row) with (min_non_missing / divisor).
        If a row has no non-missing values (should be rare after filtering), fill with 0.
        """
        df = df.copy()
        sample_cols = [c for c in df.columns if c != feature_col]

        X = df[sample_cols].apply(pd.to_numeric, errors="coerce")  # features x samples

        row_min = X.min(axis=1, skipna=True)
        fill_vals = row_min / float(self.cfg.fallback_divisor)

        # if row_min is NaN (all missing), fallback to 0
        fill_vals = fill_vals.fillna(0.0)

        X_filled = X.copy()
        nan_mask = X_filled.isna()
        # broadcast per-row fill values
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

        X = df[sample_cols].apply(pd.to_numeric, errors="coerce")   # features x samples
        X_t = X.T                                                   # samples x features

        try:
            imputer = KNNImputer(n_neighbors=int(k), weights=str(weights))
            X_imp_t = pd.DataFrame(imputer.fit_transform(X_t), index=X_t.index, columns=X_t.columns)
            df[sample_cols] = X_imp_t.T

            if self.cfg.fallback_when_nans_remain:
                if pd.isna(df[sample_cols]).to_numpy().any():
                    df = self._fallback_impute_min_over_divisor(df, feature_col)

            return df

        except Exception:
            if not self.cfg.fallback_when_knn_fails:
                raise
            # fallback if KNN fails
            df[sample_cols] = X  # ensure numeric-coerced
            df = self._fallback_impute_min_over_divisor(df, feature_col)
            return df

    # ----------------------------
    # LOG transform ONLY
    # ----------------------------
    def apply_log_transform(self, df: pd.DataFrame, feature_col: str) -> Optional[pd.DataFrame]:
        """
        Returns:
          - None if log_mode == "none"
          - otherwise returns a new dataframe (same shape) with log applied to sample columns.
        """
        mode = str(self.log_cfg.log_mode).lower().strip()
        if mode == "none":
            return None

        df2 = df.copy()
        sample_cols = [c for c in df2.columns if c != feature_col]
        X = df2[sample_cols].apply(pd.to_numeric, errors="coerce")

        X = X + float(self.log_cfg.log_offset)

        if mode == "log2":
            X = np.log2(X)
        elif mode == "log10":
            X = np.log10(X)
        else:
            raise ValueError(f"Unknown log_mode: {mode}. Use: none | log2 | log10")

        df2[sample_cols] = X
        return df2

    # ----------------------------
    # PCA compare plot: RAW vs LOG
    # ----------------------------
    def _pca_scores(self, df: pd.DataFrame, feature_col: str, meta: pd.DataFrame) -> Tuple[np.ndarray, pd.Series, List[str]]:
        """
        PCA on samples (columns). No centering/scaling beyond PCA internals.
        Returns:
          scores: (n_samples, 2)
          groups: sample->group
          samples: sample order used
        """
        sample_cols = [c for c in df.columns if c != feature_col]
        X = df[sample_cols].apply(pd.to_numeric, errors="coerce").T  # samples x features

        # align metadata to available samples
        meta2 = meta.loc[meta.index.intersection(X.index)].copy()
        X = X.loc[meta2.index]

        pca = PCA(n_components=int(self.log_cfg.pca_n_components))
        scores = pca.fit_transform(X.to_numpy())
        groups = meta2["group"].astype(str)
        samples = list(meta2.index.astype(str))
        return scores, groups, samples

    def save_pca_compare_plot(
        self,
        *,
        df_raw: pd.DataFrame,
        df_log: pd.DataFrame,
        feature_col: str,
        meta: pd.DataFrame,
        title: str,
        out_path: Path,
    ) -> None:
        """
        Creates a 2-panel PCA figure:
          left = raw PCA
          right = log PCA
        """
        scores_raw, groups_raw, _ = self._pca_scores(df_raw, feature_col, meta)
        scores_log, groups_log, _ = self._pca_scores(df_log, feature_col, meta)

        # Make sure panels use same group order
        uniq_groups = sorted(set(groups_raw.unique()).union(set(groups_log.unique())))

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        fig.suptitle(title)

        # left: RAW
        ax = axes[0]
        ax.set_title("PCA (RAW)")
        for g in uniq_groups:
            m = (groups_raw == g).to_numpy()
            if m.sum() == 0:
                continue
            ax.scatter(scores_raw[m, 0], scores_raw[m, 1], label=g)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.legend(loc="best", fontsize=8)

        # right: LOG
        ax = axes[1]
        ax.set_title(f"PCA ({self.log_cfg.log_mode.upper()})")
        for g in uniq_groups:
            m = (groups_log == g).to_numpy()
            if m.sum() == 0:
                continue
            ax.scatter(scores_log[m, 0], scores_log[m, 1], label=g)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.legend(loc="best", fontsize=8)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=int(self.log_cfg.pca_dpi), bbox_inches="tight")
        plt.close(fig)

    def save_pca_loadings_csv(
        self,
        *,
        df: pd.DataFrame,
        feature_col: str,
        meta: pd.DataFrame,
        out_path: Path,
        n_components: Optional[int] = None,
    ) -> None:
        """
        Saves PCA loadings (per metabolite/feature) to CSV.
    
        Output columns:
          feature (metabolite name), PC1..PCk loadings, plus explained variance (%).
        """
        if n_components is None:
            n_components = int(self.log_cfg.pca_n_components)
    
        sample_cols = [c for c in df.columns if c != feature_col]
    
        # X: samples x features (same orientation as your PCA plotting)
        X = df[sample_cols].apply(pd.to_numeric, errors="coerce").T
    
        # align to metadata sample order (consistent with your PCA plots)
        meta2 = meta.loc[meta.index.intersection(X.index)].copy()
        X = X.loc[meta2.index]
    
        pca = PCA(n_components=int(n_components))
        pca.fit(X.to_numpy())
    
        # feature names (metabolites) in the SAME order as columns of X
        feature_names = df[feature_col].astype(str).tolist()
    
        # loadings: (n_features, n_components)
        loadings = pca.components_.T
    
        cols = [f"PC{i+1}_loading" for i in range(loadings.shape[1])]
        out_df = pd.DataFrame(loadings, columns=cols)
        out_df.insert(0, feature_col, feature_names)
    
        # optional: add explained variance info (handy for later visualization labels)
        evr = pca.explained_variance_ratio_
        for i, v in enumerate(evr, start=1):
            out_df[f"PC{i}_explained_var_ratio"] = float(v)
    
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_path, index=False)
    def make_feature_ids_and_map(
        self,
        df: pd.DataFrame,
        feature_col: str,
        *,
        id_sep: str = "__",
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Turn non-unique feature names into unique feature_id and create:
          - feature_meta: feature_id -> feature_name, row_number
          - metabolite_feature_map: metabolite_id (=feature_name) -> list(feature_id)
    
        Returns:
          df2: same matrix but first column contains feature_id (unique)
          feature_meta: indexed by feature_id
          metabolite_feature_map: metabolite_id -> feature_ids_used
        """
        df2 = df.copy()
        names = df2[feature_col].astype(str).str.strip()
    
        # stable per-row occurrence index
        occ = names.groupby(names).cumcount() + 1
        feature_ids = names + id_sep + occ.astype(str)
    
        df2[feature_col] = feature_ids
    
        feature_meta = pd.DataFrame(
            {
                "feature_id": feature_ids.values,
                "feature_name": names.values,
                "row_number": np.arange(1, len(df2) + 1, dtype=int),
            }
        ).set_index("feature_id")
    
        metabolite_feature_map = (
            feature_meta.reset_index()
            .groupby("feature_name")["feature_id"]
            .apply(list)
            .reset_index()
            .rename(columns={"feature_name": "metabolite_id", "feature_id": "feature_ids_used"})
        )
    
        return df2, feature_meta, metabolite_feature_map

    # ----------------------------
    # main builder
    # ----------------------------
    def build_preprocessed_datasets(
        self,
        ds_paths: Dict[str, Path],
        group_paths: Dict[str, Path],
        ds_to_meta: Dict[str, str],
        groups_of_interest: Dict[str, Sequence[str]],
    ) -> Dict[str, Dict[str, Any]]:
        datasets: Dict[str, Dict[str, Any]] = {}

        for ds_name, ds_path in ds_paths.items():
            meta_key = ds_to_meta[ds_name]
            meta_path = group_paths[meta_key]
            groups_keep = groups_of_interest[meta_key]

            raw_df = self.read_table(ds_path)
            raw_meta = self.read_table(meta_path)

            feature_col = raw_df.columns[0]

            # optional: drop blanks
            blanks: List[str] = []
            df1 = raw_df
            if self.cfg.drop_blanks:
                df1, blanks = self.drop_blank_columns(df1)

            # align + subset first (needed for groupwise filter)
            df2, meta2, rep_align = self.align_and_subset_by_groups(df1, raw_meta, groups_keep, feature_col=feature_col)

            # groupwise drop sparse rows
            df3, rep_sparse = self.drop_rows_by_missing_fraction_groupwise(df2, feature_col, meta2)

            # numeric
            df3 = self.to_numeric_matrix(df3, feature_col)

            # KNN impute remaining missing values
            df3 = self.knn_impute_feature_table(df3, feature_col)

            # make feature_id unique + build maps (for PCA + downstream ORA)
            df3, feature_meta, metabolite_feature_map = self.make_feature_ids_and_map(df3, feature_col)

            # ---- NEW: keep raw + (optional) log ----
            df_raw_clean = df3
            df_log = self.apply_log_transform(df_raw_clean, feature_col)  # None if log_mode == "none"
            
            out: Dict[str, Any] = {
                "paths": {"data": ds_path, "meta": meta_path, "meta_key": meta_key},
                "groups_keep": list(groups_keep),
                "feature_col": feature_col,
                "raw_df": raw_df,
                "raw_meta": raw_meta,
                "df_raw": df_raw_clean,
                "df_log": df_log,          # either DataFrame or None
                "meta": meta2,
                "report_align": rep_align,
                "report_blanks": {"dropped_blank_cols": blanks, "n_dropped": len(blanks)},
                "report_sparse_rows": rep_sparse,
                "log_cfg": {"log_mode": self.log_cfg.log_mode, "log_offset": self.log_cfg.log_offset},
                "feature_meta": feature_meta,
                "metabolite_feature_map" : metabolite_feature_map,
            }

            # ---- NEW: PCA compare plot (only if log exists and output dir is set) ----
            if self.log_cfg.pca_compare_out_dir is not None and df_log is not None:
                out_name = f"{self.log_cfg.pca_prefix}_{ds_name}_groups={'-'.join(sorted(set(groups_keep)))}_log={self.log_cfg.log_mode}.png"
                out_path = Path(self.log_cfg.pca_compare_out_dir) / out_name

                self.save_pca_compare_plot(
                    df_raw=df_raw_clean,
                    df_log=df_log,
                    feature_col=feature_col,
                    meta=meta2,
                    title=f"{ds_name} | PCA compare (RAW vs {self.log_cfg.log_mode})",
                    out_path=out_path,
                )
                out["pca_compare_plot"] = str(out_path)

            # ---- NEW: PCA loadings export (optional) ----
            if self.log_cfg.pca_loadings_out_dir is not None:
                df_for_pca = df_log if df_log is not None else df_raw_clean
                which = self.log_cfg.log_mode if df_log is not None else "raw"
            
                out_name = (
                    f"{self.log_cfg.pca_loadings_prefix}_{ds_name}"
                    f"_groups={'-'.join(sorted(set(groups_keep)))}"
                    f"_matrix={which}_nPC={self.log_cfg.pca_n_components}.csv"
                )
                out_path = Path(self.log_cfg.pca_loadings_out_dir) / out_name
            
                self.save_pca_loadings_csv(
                    df=df_for_pca,
                    feature_col=feature_col,
                    meta=meta2,
                    out_path=out_path,
                    n_components=self.log_cfg.pca_n_components,
                )
                out["pca_loadings_csv"] = str(out_path)

            datasets[ds_name] = out

        return datasets

