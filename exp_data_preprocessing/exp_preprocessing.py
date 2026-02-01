from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Any

import pandas as pd
from sklearn.impute import KNNImputer


@dataclass
class DatasetPreprocessConfig:
    max_missing_frac: float = 0.50          # 50% rule
    min_present_frac: float = 0.50          # keep if ANY group has >= this present
    drop_blanks: bool = True
    knn_k: int = 5
    knn_weights: str = "distance"
    
    na_group_label: str = "NA"              # assign missing group labels to this
    
    fallback_divisor: float = 1000.0        # fill NaNs with (row_min / fallback_divisor)
    fallback_when_knn_fails: bool = True    # run fallback if KNN throws
    fallback_when_nans_remain: bool = True  # run fallback if KNN completes but NaNs remain


class DatasetPreprocessor:
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
                .astype(str)
                .where(~meta["group"].isna(), other=self.cfg.na_group_label)
            )
            meta["group"] = meta["group"].astype(str).str.strip()
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
        Keep a feature if at least one group has present fraction >= min_present_frac.
        Equivalent: keep if exists group with missing_frac <= max_missing_frac (default 0.50).

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
        per_group_present: Dict[str, pd.Series] = {}

        for g in group_levels:
            cols_g = groups.index[groups == g].tolist()
            if not cols_g:
                continue
            miss_frac_g = missing_mask[cols_g].mean(axis=1)
            pres_frac_g = 1.0 - miss_frac_g
            per_group_missing[str(g)] = miss_frac_g
            per_group_present[str(g)] = pres_frac_g

        # keep if ANY group passes the threshold
        # i.e. exists g with missing_frac <= max_missing_frac
        keep = pd.Series(False, index=df.index)
        for g, miss_frac_g in per_group_missing.items():
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

            # groupwise drop sparse rows (>50% missing within all groups; keep if any group passes)
            df3, rep_sparse = self.drop_rows_by_missing_fraction_groupwise(df2, feature_col, meta2)

            # numeric
            df3 = self.to_numeric_matrix(df3, feature_col)

            # KNN impute remaining missing values
            df3 = self.knn_impute_feature_table(df3, feature_col)

            datasets[ds_name] = {
                "paths": {"data": ds_path, "meta": meta_path, "meta_key": meta_key},
                "groups_keep": list(groups_keep),
                "feature_col": feature_col,
                "raw_df": raw_df,
                "raw_meta": raw_meta,
                "df": df3,
                "meta": meta2,
                "report_align": rep_align,
                "report_blanks": {"dropped_blank_cols": blanks, "n_dropped": len(blanks)},
                "report_sparse_rows": rep_sparse,
            }

        return datasets
