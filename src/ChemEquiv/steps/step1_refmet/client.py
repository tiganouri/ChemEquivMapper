from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, List
from io import StringIO
import time

import pandas as pd
import requests

REFMET_POST_URL = "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minIDS.php"


@dataclass
class RefMetClient:
    timeout_s: int = 60
    max_retries: int = 4
    backoff_s: float = 1.5
    session: Optional[requests.Session] = None

    def _sess(self) -> requests.Session:
        if self.session is None:
            self.session = requests.Session()
        return self.session

    def fetch(self, names: List[str]) -> pd.DataFrame:
        """
        POST a list of names to RefMet and return annotations as DataFrame.
        Returns empty DataFrame if the request succeeds but has no data.
        Raises only after retries; caller can handle gracefully.
        """
        clean = [n for n in names if isinstance(n, str) and n.strip()]
        payload = {"metabolite_name": "\n".join(clean)}

        last_err = None
        for attempt in range(self.max_retries):
            try:
                r = self._sess().post(REFMET_POST_URL, data=payload, timeout=self.timeout_s)
                r.raise_for_status()
                text = r.text.strip()
                if not text:
                    return pd.DataFrame()
                return pd.read_csv(StringIO(text), sep="\t")
            except Exception as e:
                last_err = e
                time.sleep(self.backoff_s * (attempt + 1))

        raise RuntimeError(f"RefMet POST failed after {self.max_retries} retries: {last_err}")
