from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal
import shutil
import zipfile

import requests

REFMET_ZIP_URL = "https://www.metabolomicsworkbench.org/databases/refmet/RefMet_python.zip"


@dataclass
class RefMetZipConfig:
    # "auto": use local if present, else try download, else None
    # "local": only local, no download
    # "download": try download, if fails fall back to local
    source: Literal["auto", "local", "download"] = "auto"

    # default local zip shipped/placed in repo
    default_local_zip: Optional[Path] = None

    # optional user-provided zip path
    zip_path: Optional[Path] = None

    # where to download/cache the zip
    download_cache_zip: Optional[Path] = None

    url: str = REFMET_ZIP_URL
    timeout_s: int = 30


def _valid_zip(path: Path) -> bool:
    try:
        return path.exists() and path.is_file() and zipfile.is_zipfile(path)
    except Exception:
        return False


def _download_zip(url: str, dest: Path, timeout_s: int) -> bool:
    """
    Safe download: returns True if OK, False if anything fails.
    Does NOT raise.
    """
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True, timeout=timeout_s) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                shutil.copyfileobj(r.raw, f)
        tmp.replace(dest)
        return _valid_zip(dest)
    except Exception:
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass
        return False


def resolve_refmet_zip(cfg: RefMetZipConfig) -> Optional[Path]:
    """
    Resolve which RefMet zip to use.
    Guaranteed not to crash on download errors.
    """
    # 1) User-provided zip
    if cfg.zip_path and _valid_zip(cfg.zip_path):
        return cfg.zip_path

    # 2) Local default zip
    local = cfg.default_local_zip
    local_ok = bool(local and _valid_zip(local))

    if cfg.source == "local":
        return local if local_ok else None

    # 3) Download / cached download
    if cfg.source in ("auto", "download") and cfg.download_cache_zip:
        if _valid_zip(cfg.download_cache_zip):
            return cfg.download_cache_zip

        ok = _download_zip(cfg.url, cfg.download_cache_zip, cfg.timeout_s)
        if ok:
            return cfg.download_cache_zip

    # 4) Fallback to local if download failed
    if local_ok:
        return local

    return None
