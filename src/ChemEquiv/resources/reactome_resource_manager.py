from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import requests
from importlib import resources
from importlib.resources import as_file

logger = logging.getLogger(__name__)

ReactomeSource = Literal["download", "bundled", "local"]


@dataclass(frozen=True)
class ReactomeTXTConfig:
    """
    Reactome mapping resolution strategy.

    Default:
      - source="download" -> try url -> if fails, fallback bundled
      - source="local" -> local_txt_path if exists else fallback bundled
      - source="bundled" -> bundled only
    """
    source: ReactomeSource = "download"

    # It downloads reactome's current version, which is the latest as this package is written 
    download_url: str = "https://reactome.org/download/current/ChEBI2Reactome.txt" #https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt" to access all the levels 

    # local option
    local_txt_path: Optional[Path] = None

    # cache location (downloaded file stored here)
    cache_txt_path: Optional[Path] = None

    # bundled default path inside package
    package_name: str = "ChemEquiv"
    bundled_rel_path: str = "data/reactome/ChEBI2Reactome.txt"    #data/reactome/ChEBI2Reactome_All_Levels.txt to access all the levels 

    timeout_s: int = 60


def _bundled_path(cfg: ReactomeTXTConfig) -> Optional[Path]:
    try:
        traversable = resources.files(cfg.package_name).joinpath(cfg.bundled_rel_path)
        with as_file(traversable) as p:
            return Path(p)
    except Exception as e:
        logger.warning(f"Could not resolve bundled Reactome TXT: {e}")
        return None


def resolve_reactome_txt(cfg: ReactomeTXTConfig, cache_dir: Path, download_url: Optional[str] = None,) -> Optional[Path]:
    """
    Resolve Reactome mapping file path with:
      - download then fallback bundled
      - local then fallback bundled
      - bundled only
    Never raises; returns None if nothing available.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_txt = cfg.cache_txt_path or (cache_dir / "ChEBI2Reactome_All_Levels.txt")

    if cfg.source == "local":
        if cfg.local_txt_path and Path(cfg.local_txt_path).exists():
            return Path(cfg.local_txt_path)

        b = _bundled_path(cfg)
        if b and b.exists():
            return b

        logger.warning("Reactome local requested but local file missing and bundled missing.")
        return None

    if cfg.source == "bundled":
        b = _bundled_path(cfg)
        if b and b.exists():
            return b
        logger.warning("Reactome bundled requested but bundled file not found.")
        return None

    if cfg.source == "download":
        try:
            url = download_url or cfg.download_url
            logger.info(f"Downloading Reactome mapping from {url}")
            r = requests.get(url, timeout=cfg.timeout_s)
            r.raise_for_status()
            cache_txt.write_bytes(r.content)
            logger.info(f"Reactome mapping downloaded -> {cache_txt}")
            return cache_txt
        except Exception as e:
            logger.warning(f"Reactome download failed; falling back to bundled. Reason: {e}")

        b = _bundled_path(cfg)
        if b and b.exists():
            return b

        logger.warning("Reactome download failed and bundled file not found.")
        return None

    return None
