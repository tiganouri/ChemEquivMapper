from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import requests
from importlib import resources
from importlib.resources import as_file

logger = logging.getLogger(__name__)

ChebiSource = Literal["download", "bundled", "local"]


@dataclass(frozen=True)
class ChebiResourceConfig:
    """
    ChEBI OBO resolution strategy.

    Default:
      - source="bundled" -> use packaged src/ChemEquiv/data/chebi/chebi.obo.gz
      - source="download" -> download_url -> cache -> fallback bundled
      - source="local" -> local_obo_gz_path -> fallback bundled
    """
    source: ChebiSource = "bundled"

    # Optional download support (nice for tutorial freshness)
    download_url: str = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"

    # local option
    local_obo_gz_path: Optional[Path] = None

    # cache path override
    cache_obo_gz_path: Optional[Path] = None

    # bundled path inside package
    package_name: str = "ChemEquiv"
    bundled_rel_path: str = "data/chebi/chebi.obo.gz"

    timeout_s: int = 120


def _bundled_path(cfg: ChebiResourceConfig) -> Optional[Path]:
    try:
        traversable = resources.files(cfg.package_name).joinpath(cfg.bundled_rel_path)
        with as_file(traversable) as p:
            return Path(p)
    except Exception as e:
        logger.warning(f"Could not resolve bundled ChEBI OBO: {e}")
        return None


def resolve_chebi_obo(cfg: ChebiResourceConfig, cache_dir: Path, download_url: Optional[str] = None) -> Optional[Path]:
    """
    Resolve a usable chebi.obo.gz path.
    Never raises; returns None if nothing is available.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_path = cfg.cache_obo_gz_path or (cache_dir / "chebi.obo.gz")

    # local
    if cfg.source == "local":
        if cfg.local_obo_gz_path and Path(cfg.local_obo_gz_path).exists():
            return Path(cfg.local_obo_gz_path)

        b = _bundled_path(cfg)
        if b and b.exists():
            return b

        logger.warning("ChEBI local requested but local file missing and bundled missing.")
        return None

    # bundled
    if cfg.source == "bundled":
        b = _bundled_path(cfg)
        if b and b.exists():
            return b
        logger.warning("ChEBI bundled requested but bundled file not found.")
        return None

    # download (with fallback bundled)
    if cfg.source == "download":
        try:
            url = download_url or cfg.download_url
            logger.info(f"Downloading ChEBI OBO from {url}")
            r = requests.get(url, timeout=cfg.timeout_s)
            r.raise_for_status()
            cache_path.write_bytes(r.content)
            logger.info(f"ChEBI OBO downloaded -> {cache_path}")
            return cache_path
        except Exception as e:
            logger.warning(f"ChEBI download failed; falling back to bundled. Reason: {e}")

        b = _bundled_path(cfg)
        if b and b.exists():
            return b

        logger.warning("ChEBI download failed and bundled file not found.")
        return None

    return None
