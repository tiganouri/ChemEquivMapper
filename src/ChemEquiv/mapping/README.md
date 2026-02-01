# Mapping layer (Reactome / KEGG)



This folder provides the mapping layer used by the pipeline to translate identifiers (ChEBI IDs and/or KEGG compound IDs) into pathway IDs for Reactome and KEGG.



It consists of:

- in-memory indices (`ReactomeIndex`, `KEGGIndex`) built from mapping tables

- a reusable mapper (`PathwayMapper`) that performs normalization, optional expansion, and writes step-specific outputs back into a DataFrame



The package exports:

- `ReactomeIndex`

- `KEGGIndex`

- `PathwayMapper`, `MappingResult`


## Role in the pipeline



Typical upstream steps produce one or more identifier columns per feature/metabolite (e.g., ChEBI IDs from name resolution, KEGG compound IDs from cross-references). This mapping layer converts those IDs into:

- **Reactome pathway stable IDs** (e.g., `R-HSA-...`)

- **KEGG pathway IDs** (e.g., `map00010`)



These pathway IDs can then be used for enrichment/ORA, reporting, or network analysis.



## Public API



This module exposes the following objects (exported in `__init__.py`): 

- `ReactomeIndex`

- `KEGGIndex`

- `PathwayMapper`

- `MappingResult`



## ReactomeIndex



### What it is

`ReactomeIndex` is an in-memory dictionary:

- key: `CHEBI:<id>`

- value: list of Reactome pathway IDs



Duplicates from the source file are preserved. 



### ID normalization

`get_pathways()` accepts either:

- `"10033"` or `"CHEBI:10033"`

and normalizes to `"CHEBI:10033"` internally.



### Building the index

Use `ReactomeIndex.from_reactome_tsv(...)` to build from a Reactome mapping TXT/TSV.

It supports the common Reactome file `ChEBI2Reactome_All_Levels.txt` which is typically headerless TSV.



Filtering:

- By default, it filters rows by `species="Homo sapiens"` when the species column is present.



Example:



```python

from pathlib import Path

from mapping import ReactomeIndex



reactome_index = ReactomeIndex.from_reactome_tsv(

    tsv_path=Path("ChEBI2Reactome_All_Levels.txt"),

    species="Homo sapiens",

    assume_header=False,
)

pathways = reactome_index.get_pathways("CHEBI:10033")



