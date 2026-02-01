"""
analysis.ora (metabolite-centric ORA)

Input:
- meta_df: the stepwise output table (Step1..Step5 columns with IDs)
- metabolite_stats: preprocessing stats table with metabolite_id + significant

Core rule:
- all identifiers (ChEBI/KEGG) derived from a metabolite inherit that metabolite’s significant flag
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple, Literal

import pandas as pd

from ChemEquiv.context import PipelineContext
from .stats import fisher_exact_right_tail, fdr_bh, metabolite_sig_map


@dataclass(frozen=True)
class StepColumnSpec:
    step: int
    chebi_col: str
    kegg_col: str


# Keep Step5 spec: if the ID columns don't exist in a given meta_df,
# build_identifier_table will return empty tables and ORA will just return empty.
DEFAULT_STEP_SPECS: Sequence[StepColumnSpec] = [
    StepColumnSpec(1, "ChEBI_ID", "KEGG_ID"),
    StepColumnSpec(2, "CHEBI_ID_Step2", "KEGG_ID_Step2"),
    StepColumnSpec(3, "CHEBI_ID_Step3", "KEGG_ID_Step3"),
    StepColumnSpec(4, "CHEBI_ID_Step4", "KEGG_ID_Step4"),
]


def _split_ids(x: object) -> List[str]:
    if x is None:
        return []
    s = str(x).strip()
    if not s or s.lower() == "nan" or s == "-":
        return []
    return [t.strip() for t in s.split(",") if t.strip()]


def _norm_chebi(x: str) -> str:
    s = str(x).strip().upper()
    if not s:
        return ""
    if s.isdigit():
        return f"CHEBI:{s}"
    if s.startswith("CHEBI:"):
        return s
    s = s.replace("CHEBI_", "CHEBI:")
    if s.startswith("CHEBI:"):
        return s
    return ""


def _norm_kegg(x: str) -> str:
    return str(x).strip()


class StepwiseORA:
    def __init__(
        self,
        *,
        ctx: PipelineContext,
        step_specs: Sequence[StepColumnSpec] = DEFAULT_STEP_SPECS,
    ) -> None:
        self.ctx = ctx
        self.step_specs = list(step_specs)

        if getattr(ctx, "reactome_index", None) is None:
            raise ValueError("ctx.reactome_index is None; load resources before ORA.")
        if getattr(ctx, "kegg_index", None) is None:
            raise ValueError("ctx.kegg_index is None; load resources before ORA.")

    def build_identifier_table(
        self,
        meta_df: pd.DataFrame,
        metabolite_stats: pd.DataFrame,
        *,
        meta_join_col: str = "Input name",
        stats_join_col: str = "metabolite_id",
        sig_col: str = "significant",
        id_col: str,
        kind: Literal["chebi", "kegg"],
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Returns:
          ids_df: unique identifiers with [identifier, significant]
          shared_df: identifiers used by >1 metabolite (diagnostic)

        Rules:
          - each identifier inherits metabolite significance
          - duplicates collapsed with OR rule (max)
        """
        if meta_join_col not in meta_df.columns:
            raise KeyError(f"meta_df missing '{meta_join_col}'")

        if id_col not in meta_df.columns:
            empty_ids = pd.DataFrame(columns=["identifier", "significant"])
            empty_shared = pd.DataFrame(columns=["identifier", "n_metabolites", "metabolites"])
            return empty_ids, empty_shared

        sig_map = metabolite_sig_map(
            metabolite_stats,
            metabolite_id_col=stats_join_col,
            significant_col=sig_col,
        )

        rows: List[Tuple[str, int, str]] = []
        for _, r in meta_df.iterrows():
            met = str(r[meta_join_col]).strip() if r[meta_join_col] is not None else ""
            if not met or met.lower() == "nan":
                continue

            s = int(sig_map.get(met, 0))  # default 0 if metabolite missing in stats

            for tok in _split_ids(r[id_col]):
                _id = _norm_chebi(tok) if kind == "chebi" else _norm_kegg(tok)
                if _id:
                    rows.append((_id, s, met))

        if not rows:
            empty_ids = pd.DataFrame(columns=["identifier", "significant"])
            empty_shared = pd.DataFrame(columns=["identifier", "n_metabolites", "metabolites"])
            return empty_ids, empty_shared

        tmp = pd.DataFrame(rows, columns=["identifier", "significant", "metabolite_id"])

        # diagnostic: shared IDs across metabolites
        counts = tmp.groupby("identifier")["metabolite_id"].nunique().reset_index(name="n_metabolites")
        shared = counts[counts["n_metabolites"] > 1].copy()
        if not shared.empty:
            ex = (
                tmp.groupby("identifier")["metabolite_id"]
                .apply(lambda x: list(pd.unique(x))[:10])
                .reset_index(name="metabolites")
            )
            shared = shared.merge(ex, on="identifier", how="left")
        else:
            shared = pd.DataFrame(columns=["identifier", "n_metabolites", "metabolites"])

        # OR-collapse at identifier level
        ids_df = tmp.groupby("identifier", as_index=False)["significant"].max()
        ids_df["significant"] = ids_df["significant"].astype(int)

        return ids_df, shared

    # ---------- pathway mapping helpers ----------

    def _reactome_pw2ids(self, chebi_ids: Set[str]) -> Dict[str, Set[str]]:
        """
        Uses ReactomeIndex.get_pathways(chebi_id) if available; otherwise uses chebi_to_pathways dict.
        Your ReactomeIndex has: _normalize_chebi, chebi_to_pathways, get_pathways.
        """
        idx = self.ctx.reactome_index
        pw2ids: Dict[str, Set[str]] = {}

        normalizer = getattr(idx, "_normalize_chebi", None)
        get_pathways = getattr(idx, "get_pathways", None)
        chebi_to_pathways = getattr(idx, "chebi_to_pathways", None)

        for cid in chebi_ids:
            cid_s = str(cid).strip()
            if not cid_s:
                continue

            cid_norm = normalizer(cid_s) if callable(normalizer) else cid_s.upper()

            if callable(get_pathways):
                pws = get_pathways(cid_norm)
            elif isinstance(chebi_to_pathways, dict):
                pws = chebi_to_pathways.get(cid_norm, [])
            else:
                raise AttributeError("ReactomeIndex has neither get_pathways() nor chebi_to_pathways dict.")

            if not pws:
                continue

            for pw in pws:
                pw2ids.setdefault(str(pw), set()).add(cid_norm)

        return pw2ids

    def _kegg_pw2ids(self, kegg_ids: Set[str]) -> Dict[str, Set[str]]:
        """
        Tries KEGGIndex.get_pathways(compound_id) if available, else tries common dict mappings.
        """
        idx = self.ctx.kegg_index
        pw2ids: Dict[str, Set[str]] = {}

        get_pathways = getattr(idx, "get_pathways", None)

        dict_map = None
        for nm in ("compound_to_pathways", "kegg_to_pathways", "kegg_compound_to_pathways"):
            m = getattr(idx, nm, None)
            if isinstance(m, dict):
                dict_map = m
                break

        for kid in kegg_ids:
            kid_s = _norm_kegg(kid)
            if not kid_s:
                continue

            if callable(get_pathways):
                pws = get_pathways(kid_s)
            elif dict_map is not None:
                pws = dict_map.get(kid_s, [])
            else:
                raise AttributeError(
                    "KEGGIndex has neither get_pathways() nor a recognized *to_pathways dict "
                    "(compound_to_pathways / kegg_to_pathways / kegg_compound_to_pathways)."
                )

            if not pws:
                continue

            for pw in pws:
                pw2ids.setdefault(str(pw), set()).add(kid_s)

        return pw2ids

    # ---------- ORA table ----------

    def _ora_table(
        self,
        *,
        foreground: Set[str],
        background: Set[str],
        pw2ids_all: Dict[str, Set[str]],
        min_pathway_size: int = 2,
        fdr_alpha: float = 0.05,
    ) -> pd.DataFrame:
        fg = set(foreground)
        bg = set(background)

        if len(bg) == 0:
            return pd.DataFrame(columns=["pathway_id", "pathway_size", "overlap", "oddsratio", "p_value", "q_value"])

        rows: List[Tuple[str, int, int, float, float]] = []
        pvals: List[float] = []

        for pw, ids_in_pw in pw2ids_all.items():
            ids_in_pw = set(ids_in_pw) & bg
            if len(ids_in_pw) < min_pathway_size:
                continue

            a = len(ids_in_pw & fg)
            if a == 0:
                continue

            b = len(fg) - a
            c = len(ids_in_pw) - a
            d = len(bg) - (a + b + c)

            odds, p = fisher_exact_right_tail(a, b, c, d)
            rows.append((str(pw), len(ids_in_pw), a, odds, p))
            pvals.append(p)

        if not rows:
            return pd.DataFrame(columns=["pathway_id", "pathway_size", "overlap", "oddsratio", "p_value", "q_value"])

        qvals = fdr_bh(pvals, alpha=fdr_alpha)

        out = pd.DataFrame(rows, columns=["pathway_id", "pathway_size", "overlap", "oddsratio", "p_value"])
        out["q_value"] = qvals
        out = out.sort_values(["q_value", "p_value", "overlap"], ascending=[True, True, False]).reset_index(drop=True)
        return out

    # ---------- Run per step + all steps ----------

    def run_step(
        self,
        meta_df: pd.DataFrame,
        metabolite_stats: pd.DataFrame,
        *,
        step: int,
        meta_join_col: str,
        stats_join_col: str = "metabolite_id",
        sig_col: str = "significant",
        min_pathway_size: int = 2,
    ) -> Dict[str, pd.DataFrame]:
        spec = next((s for s in self.step_specs if s.step == step), None)
        if spec is None:
            raise ValueError(f"No StepColumnSpec for step={step}")

        chebi_ids_df, shared_chebi = self.build_identifier_table(
            meta_df, metabolite_stats,
            meta_join_col=meta_join_col,
            stats_join_col=stats_join_col,
            sig_col=sig_col,
            id_col=spec.chebi_col,
            kind="chebi",
        )
        kegg_ids_df, shared_kegg = self.build_identifier_table(
            meta_df, metabolite_stats,
            meta_join_col=meta_join_col,
            stats_join_col=stats_join_col,
            sig_col=sig_col,
            id_col=spec.kegg_col,
            kind="kegg",
        )

        chebi_bg = set(chebi_ids_df["identifier"].astype(str))
        chebi_fg = set(chebi_ids_df.loc[chebi_ids_df["significant"] == 1, "identifier"].astype(str))

        kegg_bg = set(kegg_ids_df["identifier"].astype(str))
        kegg_fg = set(kegg_ids_df.loc[kegg_ids_df["significant"] == 1, "identifier"].astype(str))

        reactome_pw2ids = self._reactome_pw2ids(chebi_bg)
        kegg_pw2ids = self._kegg_pw2ids(kegg_bg)

        reactome_ora = self._ora_table(
            foreground=chebi_fg,
            background=chebi_bg,
            pw2ids_all=reactome_pw2ids,
            min_pathway_size=min_pathway_size,
        )
        kegg_ora = self._ora_table(
            foreground=kegg_fg,
            background=kegg_bg,
            pw2ids_all=kegg_pw2ids,
            min_pathway_size=min_pathway_size,
        )

        return {
            "chebi_ids": chebi_ids_df,
            "kegg_ids": kegg_ids_df,
            "shared_chebi_ids": shared_chebi,
            "shared_kegg_ids": shared_kegg,
            "reactome_ora_chebi": reactome_ora,
            "kegg_ora": kegg_ora,
        }

    def run_all_steps(
        self,
        meta_df: pd.DataFrame,
        metabolite_stats: pd.DataFrame,
        *,
        steps: Optional[Sequence[int]] = None,
        meta_join_col: str,
        stats_join_col: str = "metabolite_id",
        sig_col: str = "significant",
        min_pathway_size: int = 2,
    ) -> Dict[int, Dict[str, pd.DataFrame]]:
        if steps is None:
            steps = [s.step for s in self.step_specs]

        out: Dict[int, Dict[str, pd.DataFrame]] = {}
        for st in steps:
            out[st] = self.run_step(
                meta_df, metabolite_stats,
                step=st,
                meta_join_col=meta_join_col,
                stats_join_col=stats_join_col,
                sig_col=sig_col,
                min_pathway_size=min_pathway_size,
            )
        return out
