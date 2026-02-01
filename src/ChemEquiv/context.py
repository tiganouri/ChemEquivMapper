from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, List

from .mapping.mapper import PathwayMapper
from .mapping.reactome_index import ReactomeIndex
from .mapping.kegg_index import KEGGIndex

from .resources.reactome_resource_manager import ReactomeTXTConfig, resolve_reactome_txt
from .resources.kegg_resource_manager import KEGGTSVConfig, resolve_kegg_tsv


@dataclass
class PipelineContext:
    """
    Shared context that holds heavyweight resources (indices) so they are loaded once
    and reused across steps.
    """
    repo_root: Path
    reactome_index: Optional[ReactomeIndex] = None
    kegg_index: Optional[KEGGIndex] = None
    chebi_to_kegg: Optional[Dict[str, List[str]]] = None

    def build_mapper(self) -> PathwayMapper:
        return PathwayMapper(
            reactome_index=self.reactome_index,
            kegg_index=self.kegg_index,
            chebi_to_kegg=self.chebi_to_kegg or {},
        )

    # -------------------------
    # Reactome
    # -------------------------
    def load_reactome_index(self, reactome_txt_path: Path, *, species: str = "Homo sapiens") -> "PipelineContext":
        """
        Load Reactome CHEBI->pathways index from a tab-separated Reactome TXT.
        IMPORTANT: Reactome TXT (ChEBI2Reactome_All_Levels.txt) is headerless -> assume_header=False.
        """
        reactome_txt_path = Path(reactome_txt_path)
        self.reactome_index = ReactomeIndex.from_reactome_tsv(
            reactome_txt_path,
            species=species,
            assume_header=False,
        )
        return self

    def load_reactome_index_auto(
        self,
        *,
        cfg: ReactomeTXTConfig = ReactomeTXTConfig(),
        cache_dir: Optional[Path] = None,
        species: str = "Homo sapiens",
    ) -> "PipelineContext":
        """
        Resolve Reactome mapping (download/local/bundled) then load index.
        """
        cache_dir = cache_dir or (Path(self.repo_root) / ".cache" / "reactome")
        reactome_path = resolve_reactome_txt(cfg, cache_dir=cache_dir)
        if reactome_path is None:
            raise FileNotFoundError(
                "Could not resolve Reactome mapping file (download/local/bundled all failed)."
            )
        return self.load_reactome_index(reactome_path, species=species)

    # -------------------------
    # KEGG
    # -------------------------
    def load_kegg_index(self, kegg_tsv_path: Path) -> "PipelineContext":
        """
        Load KEGG compound->pathways index from a TSV with columns:
          kegg_compound_id, kegg_pathway_id
        """
        kegg_tsv_path = Path(kegg_tsv_path)
        self.kegg_index = KEGGIndex.from_compound_pathway_tsv(kegg_tsv_path, assume_header=True)
        return self

    def load_kegg_index_auto(
        self,
        *,
        cfg: KEGGTSVConfig = KEGGTSVConfig(),  # default is bundled
        cache_dir: Optional[Path] = None,
    ) -> "PipelineContext":
        """
        Resolve KEGG TSV (bundled/local/download) then load index.
        """
        cache_dir = cache_dir or (Path(self.repo_root) / ".cache" / "kegg")
        kegg_path = resolve_kegg_tsv(cfg, cache_dir=cache_dir)
        if kegg_path is None:
            raise FileNotFoundError(
                "Could not resolve KEGG TSV (download/local/bundled all failed)."
            )
        return self.load_kegg_index(kegg_path)

    # -------------------------
    # One-call convenience
    # -------------------------
    def load_resources_default(
        self,
        *,
        reactome_cfg: ReactomeTXTConfig = ReactomeTXTConfig(),
        kegg_cfg: KEGGTSVConfig = KEGGTSVConfig(),  # default bundled
        species: str = "Homo sapiens",
    ) -> "PipelineContext":
        """
        Load both Reactome and KEGG indices with default configs.
        """
        self.load_reactome_index_auto(cfg=reactome_cfg, species=species)
        self.load_kegg_index_auto(cfg=kegg_cfg)
        return self
