from __future__ import annotations

from collections import deque
from typing import Optional, Set, Tuple

import pronto


def _safe_get_term(onto: pronto.Ontology, chebi_id: str) -> Optional[pronto.Term]:
    try:
        return onto[chebi_id]
    except Exception:
        return None


def get_is_a_children_for_ids(
    *,
    chebi_ids: Set[str],
    onto: pronto.Ontology,
    max_depth: int,
) -> Set[str]:
    """
    is_a traversal (children / more specific terms).

    For each CHEBI:XXXXX in chebi_ids:
      - traverse subclasses() up to max_depth
      - return all discovered CHEBI children IDs
      - exclude starting IDs

    Notes:
      - defensive for pronto variations
      - returns CHEBI:XXXXX only
    """
    if not chebi_ids or max_depth <= 0:
        return set()

    out: Set[str] = set()

    for cid in chebi_ids:
        cid = str(cid).strip().upper()
        if not cid.startswith("CHEBI:"):
            continue

        start = _safe_get_term(onto, cid)
        if start is None:
            continue

        visited_terms = {start}
        q: deque[Tuple[pronto.Term, int]] = deque([(start, 0)])

        while q:
            node, d = q.popleft()
            if d >= max_depth:
                continue

            # subclasses() is the cleanest (children)
            if hasattr(node, "subclasses"):
                try:
                    for child in node.subclasses():
                        if child in visited_terms:
                            continue
                        visited_terms.add(child)
                        q.append((child, d + 1))
                        tid = getattr(child, "id", None)
                        if isinstance(tid, str) and tid.startswith("CHEBI:") and tid not in chebi_ids:
                            out.add(tid)
                except Exception:
                    # if subclasses() fails in some pronto versions, skip
                    pass

    out -= chebi_ids
    return out
