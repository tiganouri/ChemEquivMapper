# Step 2 — ChEBI lookup from names + CHEBI→KEGG expansion + pathway mapping

This folder implements **Step 2** of the pipeline. It takes the output of Step 1 (RefMet annotations),
builds fast lookup indices from the ChEBI ontology, and uses them to infer additional ChEBI IDs from
metabolite names (exact + fuzzy matching). It then expands those ChEBI IDs to KEGG compound IDs using
ChEBI xrefs, and finally maps CHEBI/KEGG IDs to Reactome/KEGG pathways (Step2 mapping columns).

Main entrypoint:

- `Step2ChebiLookup`

Exports:

- `Step2ChebiLookup`

## Role in the pipeline

Step 1 uses RefMet to annotate metabolites and may provide ChEBI/KEGG IDs.
Step 2 improves coverage by additionally resolving ChEBI IDs directly from:
- **Standardized name** (preferred)
- **Input name** (fallback)

It produces:
- Step2 ChEBI inference columns
- Step2 KEGG compounds derived from inferred ChEBI
- Step2 pathway columns for Reactome/KEGG via the shared mapper

## Public API

The package exports `Step2ChebiLookup`.

## Inputs

### Input DataFrame (Step1 output)

`Step2ChebiLookup.run(meta_step1, ...)` expects a DataFrame that includes name columns from Step1:

- `Standardized name` (preferred lookup key)
- `Input name` (fallback lookup key) 

The DataFrame index is treated as the row identifier and is preserved.

## Outputs

`run(...)` returns a tuple:

1) `meta_step2`: the input DataFrame plus Step2 ID and pathway columns
2) `chebi_to_kegg`: a dict mapping `"CHEBI:xxxxx" -> { "C00031", ... }` (built or loaded from cache) 


### Step2 CHEBI/KEGG columns

The Step2 schema defines canonical output columns and guarantees they exist:
- `CHEBI_Exact_From_Names`
- `CHEBI_Fuzzy_From_Names`
- `CHEBI_From_Names_All`
- `CHEBI_ID_Step2`

- `KEGG_Exact_From_Names`
- `KEGG_Fuzzy_From_Names`
- `KEGG_From_Names_All`
- `KEGG_ID_Step2` 

Final Step2 ID columns (`CHEBI_ID_Step2`, `KEGG_ID_Step2`) are normalized into canonical,
deduped, sorted comma-separated strings.

### Step2 pathway columns

After Step2 fills `CHEBI_ID_Step2` and `KEGG_ID_Step2`, it uses the shared mapper to write

step-specific pathway outputs:
- `ChEBI_ID_Step2_reactome`
- `Reactome_Pathways_Step2`
- `KEGG_Compounds_Step2`
- `KEGG_Pathways_Step2` 


## What Step 2 does

### 1) Resolve ChEBI ontology file (chebi.obo.gz)


Step2 uses the shared resource manager to resolve `chebi.obo.gz` using a strategy:
- bundled / local / download (default bundled)

### 2) Build or load SQLite-backed name indices


Step2 builds two indices from ChEBI ontology names + synonyms:

- `exact`: key = `norm(name)` (lowercase, trimmed, collapse whitespace)

- `fuzzy`: key = `norm_remove_special(name)` (remove non-alnum, collapse whitespace)


Both map from normalized key -> list of `CHEBI:xxxxx`. 


The indices are stored in SQLite for reuse and loaded into memory for fast lookups.


### 3) Name → ChEBI lookup logic

For each row, Step2 performs ChEBI name lookup using **both** name fields when available:

- `Standardized name` (if present / non-empty)
- `Input name` (if present / non-empty)

The lookup results are merged:
- exact hits: `Standardized name` hits first, then any additional `Input name` hits
- fuzzy hits: `Standardized name` hits first, then any additional `Input name` hits
- combined hits: union of exact+fuzzy across both names (deduplicated, stable order)


It returns:

- exact hits (from `exact` index)
- fuzzy hits (from `fuzzy` index)
- combined hits (exact first, then fuzzy not already in exact)


These are written into the CHEBI Step2 columns.


### 4) Build/load CHEBI → KEGG mapping (SQLite cached)


Step2 extracts KEGG compound IDs from ChEBI ontology term xrefs:

- keeps only xrefs that look KEGG-related
- regex extracts `Cd{5}` tokens 


The resulting mapping is cached in SQLite so it is built only once. 


Step2 then populates:

- `KEGG_ID_Step2` by expanding the inferred `CHEBI_ID_Step2` through this map 


### 5) Map Step2 IDs to Reactome/KEGG pathways



Step2 loads Reactome/KEGG mapping resources using `PipelineContext` and applies the shared mapper.

By default:

- Reactome mapping uses `source="download"` and filters to `"Homo sapiens"`
- KEGG mapping uses `source="bundled"` 

## Usage


```python

from pathlib import Path
from step2 import Step2ChebiLookup


# meta_step1 is the output DataFrame from Step1 (must include name columns)
# meta_step1 = ...

step2 = Step2ChebiLookup(repo_root=Path("."))


meta_step2, chebi_to_kegg = step2.run(meta_step1)



