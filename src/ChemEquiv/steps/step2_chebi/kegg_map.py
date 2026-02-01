from __future__ import annotations

from pathlib import Path
from typing import Dict, Set, Optional, Iterable
import sqlite3
import logging
import re

import pronto

logger = logging.getLogger(__name__)

_KEGG_TOKEN = re.compile(r"\b(C\d{5})\b")


def _init_db(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS chebi_kegg (chebi TEXT NOT NULL, kegg TEXT NOT NULL, PRIMARY KEY (chebi, kegg))"
    )
    cur.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT)")
    con.commit()
    con.close()


def _db_count_pairs(db_path: Path) -> int:
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("SELECT COUNT(*) FROM chebi_kegg")
    n = int(cur.fetchone()[0])
    con.close()
    return n


def _load_db(db_path: Path) -> Dict[str, Set[str]]:
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("SELECT chebi, kegg FROM chebi_kegg")
    out: Dict[str, Set[str]] = {}
    for chebi, kegg in cur.fetchall():
        out.setdefault(chebi, set()).add(kegg)
    con.close()
    return out


def _store_pairs(db_path: Path, pairs: Dict[str, Set[str]]) -> None:
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = [(c, k) for c, ks in pairs.items() for k in ks]
    cur.executemany("INSERT OR IGNORE INTO chebi_kegg (chebi, kegg) VALUES (?, ?)", rows)
    con.commit()
    con.close()


def _extract_kegg_from_xrefs(xrefs: Iterable[object]) -> Set[str]:
    """
    pronto term.xrefs yields objects; we convert to string and regex out Cxxxxx.
    """
    out: Set[str] = set()
    for x in xrefs:
        xs = str(x)
        # keep only KEGG-related xrefs to reduce false matches
        if "KEGG" not in xs.upper():
            continue
        for hit in _KEGG_TOKEN.findall(xs.upper()):
            out.add(hit)
    return out


def build_chebi_to_kegg_map(obo_path: Path) -> Dict[str, Set[str]]:
    """
    Parse chebi.obo via pronto and extract KEGG xrefs.
    Output: "CHEBI:xxxxx" -> {"C00031", ...}
    """
    ont = pronto.Ontology(str(obo_path))
    out: Dict[str, Set[str]] = {}

    for term in ont.terms():
        tid = term.id
        if not tid or not tid.startswith("CHEBI:"):
            continue
        kegg_ids = _extract_kegg_from_xrefs(term.xrefs)
        if kegg_ids:
            out[tid] = kegg_ids

    return out

def build_chebi_to_kegg_map_from_ontology(ont: pronto.Ontology) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    for term in ont.terms():
        tid = term.id
        if not tid or not tid.startswith("CHEBI:"):
            continue
        kegg_ids = _extract_kegg_from_xrefs(term.xrefs)
        if kegg_ids:
            out[tid] = kegg_ids
    return out

def get_chebi_to_kegg_map(obo_path: Optional[Path], db_path: Path, *, ont: Optional[pronto.Ontology] = None) -> Dict[str, Set[str]]:
    """
    Load CHEBI->KEGG mapping from sqlite cache, building it once from obo if needed.
    - If db has pairs, loads and returns.
    - If db empty and obo_path exists, builds then caches.
    - If no obo_path, returns {}.
    """
    _init_db(db_path)

    if _db_count_pairs(db_path) > 0:
        return _load_db(db_path)

    if obo_path is None and ont is None:
        logger.warning("No chebi.obo available; CHEBI->KEGG map will be empty.")
        return {}

    if ont is None:
        logger.info(f"Building CHEBI->KEGG map from ontology: {obo_path}")
        built = build_chebi_to_kegg_map(obo_path)  # existing
    else:
        built = build_chebi_to_kegg_map_from_ontology(ont)
        
    _store_pairs(db_path, built)
    return built
