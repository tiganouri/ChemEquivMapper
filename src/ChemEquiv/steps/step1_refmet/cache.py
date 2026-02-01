from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
import pandas as pd

REFMET_CACHE_TABLE = "refmet_cache_v1"


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


@dataclass
class SQLiteRefMetCache:
    """
    Tiny SQLite cache:
    key = normalized input name
    value = JSON payload of the RefMet output row for that input name
    """
    db_path: Path

    def _connect(self) -> sqlite3.Connection:
        _ensure_parent_dir(self.db_path)
        con = sqlite3.connect(str(self.db_path))
        con.execute(f"""
            CREATE TABLE IF NOT EXISTS {REFMET_CACHE_TABLE} (
                key TEXT PRIMARY KEY,
                payload_json TEXT NOT NULL
            )
        """)
        return con

    @staticmethod
    def normalize_key(name: str) -> str:
        return " ".join(str(name).strip().lower().split())

    def get_many(self, names: List[str]) -> Tuple[pd.DataFrame, List[str]]:
        """
        Returns:
          df_cached: annotations for cached names (includes "Input name" = original name)
          missing: names not in cache
        """
        keys = [self.normalize_key(n) for n in names]
        missing: List[str] = []
        rows = []

        con = self._connect()
        try:
            for orig_name, k in zip(names, keys):
                cur = con.execute(
                    f"SELECT payload_json FROM {REFMET_CACHE_TABLE} WHERE key=?",
                    (k,)
                )
                hit = cur.fetchone()
                if hit is None:
                    missing.append(orig_name)
                else:
                    payload = json.loads(hit[0])
                    payload["Input name"] = orig_name
                    rows.append(payload)
        finally:
            con.close()

        df_cached = pd.DataFrame(rows) if rows else pd.DataFrame()
        return df_cached, missing

    def put_many(self, df: pd.DataFrame, input_col: str = "Input name") -> None:
        """
        Store RefMet output rows keyed by normalized "Input name".
        """
        if df.empty or input_col not in df.columns:
            return

        con = self._connect()
        try:
            for _, row in df.iterrows():
                name = row.get(input_col)
                if not isinstance(name, str) or not name.strip():
                    continue
                k = self.normalize_key(name)
                payload_json = json.dumps(row.to_dict(), ensure_ascii=False)
                con.execute(
                    f"INSERT OR REPLACE INTO {REFMET_CACHE_TABLE} (key, payload_json) VALUES (?, ?)",
                    (k, payload_json)
                )
            con.commit()
        finally:
            con.close()
