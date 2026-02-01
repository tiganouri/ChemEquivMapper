from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, List, Tuple

import requests
from importlib import resources
from importlib.resources import as_file

logger = logging.getLogger(__name__)

KEGGSource = Literal["download", "bundled", "local"]


@dataclass(frozen=True)
class KEGGTSVConfig:
    """
    KEGG mapping resolution strategy for compound->pathway TSV.

    - source="bundled": use packaged TSV only
    - source="local": use local_tsv_path if provided and exists, else fallback bundled
    - source="download": download from KEGG REST (map pathways) into cache, else fallback bundled

    TSV format:
      kegg_compound_id <tab> kegg_pathway_id
    """
    source: KEGGSource = "bundled"

    # local option
    local_tsv_path: Optional[Path] = None

    # cache location (downloaded TSV stored here)
    cache_tsv_path: Optional[Path] = None

    # bundled default path inside package
    package_name: str = "GEPC"
    bundled_rel_path: str = "data/kegg/kegg_compound_to_pathway.tsv"

    # download tuning (KEGG REST)
    timeout_s: int = 60
    max_retries: int = 4
    sleep_s: float = 0.2
    batch_size: int = 10


def _bundled_path(cfg: KEGGTSVConfig) -> Optional[Path]:
    try:
        traversable = resources.files(cfg.package_name).joinpath(cfg.bundled_rel_path)
        with as_file(traversable) as p:
            return Path(p)
    except Exception as e:
        logger.warning(f"Could not resolve bundled KEGG TSV: {e}")
        return None


def _get_text(url: str, timeout_s: int, max_retries: int) -> str:
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, timeout=timeout_s)
            r.raise_for_status()
            return r.text
        except Exception as e:
            last_err = e
            time.sleep(1.0 * attempt)
    raise RuntimeError(f"GET failed: {url}\nLast error: {last_err}")


def _list_map_pathways(cfg: KEGGTSVConfig) -> List[str]:
    """
    KEGG map pathways: map00010, map01100, ...
    """
    url = "https://rest.kegg.jp/list/pathway"
    text = _get_text(url, cfg.timeout_s, cfg.max_retries)

    ids: List[str] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        left = line.split("\t", 1)[0].strip()  # path:map00010
        if ":" in left:
            left = left.split(":", 1)[1].strip()
        if left.startswith("map") and len(left) == 8:
            ids.append(left)
    return ids


def _fetch_compound_links_for_maps(cfg: KEGGTSVConfig, map_ids: List[str]) -> List[Tuple[str, str]]:
    """
    Fetch (compound_id, map_pathway_id) pairs using:
      https://rest.kegg.jp/link/cpd/map00010
    Supports batching dbentries with '+'.
    """
    pairs: List[Tuple[str, str]] = []

    total = len(map_ids)
    batches = (total + cfg.batch_size - 1) // cfg.batch_size

    for bi in range(batches):
        batch = map_ids[bi * cfg.batch_size : (bi + 1) * cfg.batch_size]
        dbentries = "+".join(batch)
        url = f"https://rest.kegg.jp/link/cpd/{dbentries}"

        text = _get_text(url, cfg.timeout_s, cfg.max_retries)
        time.sleep(cfg.sleep_s)

        for line in text.splitlines():
            if not line.strip():
                continue
            a, b = line.split("\t")
            pw = a.split(":", 1)[1].strip() if ":" in a else a.strip()      # map00010
            comp = b.split(":", 1)[1].strip() if ":" in b else b.strip()   # C00031
            if pw and comp:
                pairs.append((comp, pw))

        if (bi + 1) in (1, batches) or (bi + 1) % 25 == 0:
            logger.info(f"KEGG download progress: batch {bi+1}/{batches}, pairs={len(pairs):,}")

    return pairs


def _download_kegg_tsv(cfg: KEGGTSVConfig, out_path: Path) -> Optional[Path]:
    """
    Build the TSV from KEGG REST and write to out_path.
    Returns out_path on success, None on failure.
    """
    try:
        logger.info("Listing KEGG map pathways...")
        map_ids = _list_map_pathways(cfg)
        logger.info(f"Found {len(map_ids):,} map pathways. Downloading compound links...")

        pairs = _fetch_compound_links_for_maps(cfg, map_ids)
        if not pairs:
            raise RuntimeError("KEGG REST returned zero compound-pathway pairs.")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8", newline="\n") as f:
            f.write("kegg_compound_id\tkegg_pathway_id\n")
            for comp, pw in pairs:
                f.write(f"{comp}\t{pw}\n")

        logger.info(f"KEGG TSV downloaded/built -> {out_path}")
        return out_path
    except Exception as e:
        logger.warning(f"KEGG download failed: {e}")
        return None


def resolve_kegg_tsv(cfg: KEGGTSVConfig, cache_dir: Path) -> Optional[Path]:
    """
    Resolve KEGG compound->pathway TSV path with:
      - bundled / local / download then fallback bundled
    Never raises; returns None if nothing available.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_tsv = cfg.cache_tsv_path or (cache_dir / "kegg_compound_to_pathway.tsv")

    # 1) local mode
    if cfg.source == "local":
        if cfg.local_tsv_path and Path(cfg.local_tsv_path).exists():
            return Path(cfg.local_tsv_path)

        b = _bundled_path(cfg)
        if b and b.exists():
            return b

        logger.warning("KEGG local requested but local TSV missing and bundled missing.")
        return None

    # 2) bundled mode
    if cfg.source == "bundled":
        b = _bundled_path(cfg)
        if b and b.exists():
            return b
        logger.warning("KEGG bundled requested but bundled TSV not found.")
        return None

    # 3) download mode
    if cfg.source == "download":
        p = _download_kegg_tsv(cfg, cache_tsv)
        if p and p.exists():
            return p

        b = _bundled_path(cfg)
        if b and b.exists():
            logger.info("Using bundled KEGG TSV after download failure.")
            return b

        logger.warning("KEGG download failed and bundled TSV not found.")
        return None

    return None
