"""
Mapping layer: load indices and map identifiers to pathways (Reactome/KEGG).
"""

from .reactome_index import ReactomeIndex
from .kegg_index import KEGGIndex
from .mapper import PathwayMapper, MappingResult

__all__ = [
    "ReactomeIndex",
    "KEGGIndex",
    "PathwayMapper",
    "MappingResult",
]
