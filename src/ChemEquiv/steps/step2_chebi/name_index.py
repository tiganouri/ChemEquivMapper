from __future__ import annotations

import logging
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pronto

from .supplementary import norm, norm_remove_special

logger = logging.getLogger(__name__)

NameIndex = Dict[str, List[str]]  # key -> list of CHEBI:xxxx


@dataclass(frozen=True)
class NameIndexConfig:
    """
    SQLite-backed name/synonym indices for ChEBI ontology terms.

    DB tables:
      - exact(key TEXT, chebi TEXT)
      - fuzzy(key TEXT, chebi TEXT)
      - meta(k TEXT PRIMARY KEY, v TEXT)
    """
    db_path: Path


def _init_db(cfg: NameIndexConfig) -> None:
    cfg.db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS exact (key TEXT NOT NULL, chebi TEXT NOT NULL)")
    cur.execute("CREATE TABLE IF NOT EXISTS fuzzy (key TEXT NOT NULL, chebi TEXT NOT NULL)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_exact_key ON exact(key)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_fuzzy_key ON fuzzy(key)")
    cur.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT)")
    con.commit()
    con.close()


def _count_rows(cfg: NameIndexConfig) -> Tuple[int, int]:
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    cur.execute("SELECT COUNT(*) FROM exact")
    n_exact = int(cur.fetchone()[0])
    cur.execute("SELECT COUNT(*) FROM fuzzy")
    n_fuzzy = int(cur.fetchone()[0])
    con.close()
    return n_exact, n_fuzzy


def _clear_indices(cfg: NameIndexConfig) -> None:
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    cur.execute("DELETE FROM exact")
    cur.execute("DELETE FROM fuzzy")
    con.commit()
    con.close()


def _store_pairs(cfg: NameIndexConfig, table: str, pairs: List[Tuple[str, str]], chunk: int = 50000) -> None:
    """
    pairs: [(key, chebi), ...]
    """
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    q = f"INSERT INTO {table} (key, chebi) VALUES (?, ?)"
    for i in range(0, len(pairs), chunk):
        cur.executemany(q, pairs[i : i + chunk])
        con.commit()
    con.close()


def _set_meta(cfg: NameIndexConfig, k: str, v: str) -> None:
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    cur.execute("INSERT OR REPLACE INTO meta (k, v) VALUES (?, ?)", (k, v))
    con.commit()
    con.close()


def _get_meta(cfg: NameIndexConfig, k: str) -> Optional[str]:
    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()
    cur.execute("SELECT v FROM meta WHERE k=?", (k,))
    row = cur.fetchone()
    con.close()
    return row[0] if row else None


def _term_all_names(term: pronto.Term) -> Set[str]:
    """
    Collect primary name + synonyms text from a pronto term.

    pronto term API differs slightly across versions, so we use safe probing.
    """
    out: Set[str] = set()

    # primary label
    name = getattr(term, "name", None)
    if isinstance(name, str) and name.strip():
        out.add(name.strip())

    # synonyms (often set of Synonym objects)
    syns = getattr(term, "synonyms", None)
    if syns:
        try:
            for s in syns:
                # Synonym may have .description or stringify
                desc = getattr(s, "description", None)
                if isinstance(desc, str) and desc.strip():
                    out.add(desc.strip())
                else:
                    ss = str(s).strip()
                    # Sometimes str(s) includes quotes/metadata; still useful
                    if ss:
                        out.add(ss)
        except Exception:
            pass

    return out


def build_name_indices_from_obo(
    obo_path: Path,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Build pairs lists for sqlite:
      exact_pairs = [(norm(name), chebi_id), ...]
      fuzzy_pairs = [(norm_remove_special(name), chebi_id), ...]

    Returns (exact_pairs, fuzzy_pairs).
    """
    onto = pronto.Ontology(str(obo_path))

    exact_pairs: List[Tuple[str, str]] = []
    fuzzy_pairs: List[Tuple[str, str]] = []

    n_terms = 0
    for term in onto.terms():
        tid = getattr(term, "id", None)
        if not isinstance(tid, str) or not tid.startswith("CHEBI:"):
            continue

        n_terms += 1
        names = _term_all_names(term)

        for nm in names:
            k_exact = norm(nm)
            if k_exact:
                exact_pairs.append((k_exact, tid))

            k_fuzzy = norm_remove_special(nm)
            if k_fuzzy:
                fuzzy_pairs.append((k_fuzzy, tid))

    logger.info(f"Built name indices from ontology: terms_used={n_terms}, exact_pairs={len(exact_pairs)}, fuzzy_pairs={len(fuzzy_pairs)}")
    return exact_pairs, fuzzy_pairs

def build_name_indices_from_ontology(onto: pronto.Ontology) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    exact_pairs: List[Tuple[str, str]] = []
    fuzzy_pairs: List[Tuple[str, str]] = []

    n_terms = 0
    for term in onto.terms():
        tid = getattr(term, "id", None)
        if not isinstance(tid, str) or not tid.startswith("CHEBI:"):
            continue
        n_terms += 1
        names = _term_all_names(term)
        for nm in names:
            k_exact = norm(nm)
            if k_exact:
                exact_pairs.append((k_exact, tid))
            k_fuzzy = norm_remove_special(nm)
            if k_fuzzy:
                fuzzy_pairs.append((k_fuzzy, tid))

    logger.info(f"Built name indices from ontology: terms_used={n_terms}, exact_pairs={len(exact_pairs)}, fuzzy_pairs={len(fuzzy_pairs)}")
    return exact_pairs, fuzzy_pairs

def ensure_name_index_db(cfg: NameIndexConfig, 
    obo_path: Path, 
    force_rebuild: bool = False, 
    *, 
    onto: Optional[pronto.Ontology] = None
    ) -> None:

    """
    Ensures sqlite DB exists and has data. Builds if missing/empty or force_rebuild.
    """
    _init_db(cfg)

    n_exact, n_fuzzy = _count_rows(cfg)
    if not force_rebuild and n_exact > 0 and n_fuzzy > 0:
        return

    logger.info(f"Building ChEBI name/synonym indices into sqlite: {cfg.db_path}")
    _clear_indices(cfg)

    if onto is None:
        exact_pairs, fuzzy_pairs = build_name_indices_from_obo(obo_path)  # existing path
    else:
        exact_pairs, fuzzy_pairs = build_name_indices_from_ontology(onto)
    # store
    _store_pairs(cfg, "exact", exact_pairs)
    _store_pairs(cfg, "fuzzy", fuzzy_pairs)

    # minimal provenance
    _set_meta(cfg, "obo_path", str(obo_path))
    _set_meta(cfg, "obo_mtime", str(obo_path.stat().st_mtime))

    logger.info("Name index DB ready.")


def load_indices_from_db(cfg: NameIndexConfig) -> Tuple[NameIndex, NameIndex]:
    """
    Load full indices into memory dicts (fast lookups).
    For very large DBs, this still fits usually (ChEBI size is manageable),
    but later we can switch to on-demand sqlite queries if needed.
    """
    _init_db(cfg)

    exact: NameIndex = {}
    fuzzy: NameIndex = {}

    con = sqlite3.connect(cfg.db_path)
    cur = con.cursor()

    cur.execute("SELECT key, chebi FROM exact")
    for k, cid in cur.fetchall():
        exact.setdefault(k, []).append(cid)

    cur.execute("SELECT key, chebi FROM fuzzy")
    for k, cid in cur.fetchall():
        fuzzy.setdefault(k, []).append(cid)

    con.close()
    return exact, fuzzy
