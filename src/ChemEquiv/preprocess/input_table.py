from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
import re

import pandas as pd


def _normalize_name(x: object) -> str:
    x = "" if x is None else str(x)
    x = x.strip()
    x = re.sub(r"\s+", " ", x)
    return x


@dataclass(frozen=True)
class InputTableConfig:
    """
    Input reading + canonicalisation options.

    feature_col:
      - If None, uses first column as feature names.
      - Else, uses the named column.
    """
    feature_col: Optional[str] = None

    keep_row_number: bool = True
    duplicate_report: bool = True

    # Ensures stable unique IDs even when names repeat
    make_unique_ids: bool = True
    id_suffix_sep: str = "__"  # "name__2", "name__3", ...


class InputTablePreprocessor:
    """
    Reads an input feature table (features as rows, samples as columns) and guarantees:
      - a stable UNIQUE feature_id index
      - a canonicalised feature_name column (can repeat)
      - a numeric-only abundance table

    Returns:
      abundance_feature: DataFrame (features x samples), index=feature_id, columns=samples
      feature_meta: DataFrame indexed by feature_id with feature_id, feature_name, row_number
    """

    def __init__(self, cfg: InputTableConfig = InputTableConfig()):
        self.cfg = cfg

    def read(self, path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
        path = Path(path)

        if path.suffix.lower() == ".tsv":
            df = pd.read_csv(path, sep="\t")
        elif path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
        elif path.suffix.lower() in [".xlsx", ".xls"]:
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path, sep=None, engine="python")

        if df.shape[1] < 2:
            raise ValueError("Input must have >=2 columns (feature name + >=1 sample column).")

        feature_col = self.cfg.feature_col or df.columns[0]

        df = df.reset_index(drop=True)
        df["_row_number_"] = df.index + 1
        df["_feature_name_norm_"] = df[feature_col].map(_normalize_name).astype(str)

        # Duplicate reporting (kept)
        dup_mask = df["_feature_name_norm_"].duplicated(keep=False)
        if self.cfg.duplicate_report and dup_mask.any():
            dup_df = df.loc[dup_mask, ["_row_number_", "_feature_name_norm_"]]
            print("\n[Preprocess] Duplicate canonical names detected (kept):")
            for name, grp in dup_df.groupby("_feature_name_norm_"):
                rows = grp["_row_number_"].tolist()
                print(f"  - '{name}' at rows: {rows}")
            print("[Preprocess] End duplicate report.\n")

        feature_name = df["_feature_name_norm_"]

        # Build unique feature_id
        if self.cfg.make_unique_ids:
            rank = feature_name.groupby(feature_name).cumcount() + 1
            feature_id = feature_name.where(rank == 1, feature_name + self.cfg.id_suffix_sep + rank.astype(str))
        else:
            feature_id = feature_name

        # Numeric-only abundance (drop feature column + helper cols)
        abundance = df.drop(columns=[feature_col, "_row_number_", "_feature_name_norm_"]).copy()
        abundance = abundance.apply(pd.to_numeric, errors="coerce")
        abundance.index = feature_id

        meta = pd.DataFrame(index=feature_id)
        meta["feature_id"] = feature_id.values
        meta["feature_name"] = feature_name.values
        if self.cfg.keep_row_number:
            meta["row_number"] = df["_row_number_"].values

        return abundance, meta
