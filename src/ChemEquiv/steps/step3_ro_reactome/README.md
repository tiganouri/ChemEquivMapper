# Step 3 — RO-based ChEBI expansion + KEGG enrichment + pathway mapping


This folder implements **Step 3** of the pipeline: it expands the ChEBI identifiers obtained in Step 2 by
traversing **RO equivalence relationships** in the ChEBI ontology (e.g., tautomer/enantiomer/conjugate acid/base),
then enriches KEGG compound coverage from the expanded ChEBI set, and finally maps both ChEBI and KEGG IDs
to Reactome/KEGG pathway IDs using the shared mapper.


Main entrypoint:

- `Step3RO`

## Role in the pipeline

Step 2 produces:
- `CHEBI_ID_Step2` (required)
- `KEGG_ID_Step2` (optional)



Step 3 improves pathway mapping coverage by:

1) Expanding ChEBI IDs via RO equivalence relations (ontology graph traversal)
2) Deriving additional KEGG compound IDs from the expanded ChEBI set (using ChEBI xrefs)
3) Mapping Step3 IDs to Reactome/KEGG pathway IDs (Step3 mapping columns)

## Public API

This package exports `Step3RO`. 

## Inputs

`Step3RO.run(df, ...)` expects a pandas DataFrame that contains:

Required:

- `CHEBI_ID_Step2` (comma-separated ChEBI IDs)

Optional:

- `KEGG_ID_Step2` (comma-separated KEGG compound IDs) 

You can override the column names:

- `chebi_col_step2` (default `"CHEBI_ID_Step2"`)
- `kegg_col_step2` (default `"KEGG_ID_Step2"`) 

## Outputs

Step 3 returns a copy of the input DataFrame with additional columns:


New Step3 ID columns:

- `CHEBI_RO_Equivalents` (RO-expanded equivalents only, excluding original Step2 IDs)
- `CHEBI_ID_Step3` (Step2 ChEBI IDs ∪ RO equivalents)
- `KEGG_ID_Step3` (Step2 KEGG IDs ∪ KEGG compounds derived from Step3 ChEBI) 

New Step3 pathway columns (written by the shared mapper with `step=3`):

- `ChEBI_ID_Step3_reactome`
- `Reactome_Pathways_Step3`
- `KEGG_Compounds_Step3`
- `KEGG_Pathways_Step3` 

## What Step 3 does

### 1) Resolve ChEBI ontology and build CHEBI→KEGG map

Step 3 resolves `chebi.obo.gz` via the shared resource resolver (bundled/local/download strategies),
then loads it into `pronto.Ontology`.

It also builds a CHEBI→KEGG map using the Step2 logic and caches it in:
- `<repo_root>/.cache/step3_ro/chebi_to_kegg.sqlite` 

### 2) RO equivalence traversal (graph expansion)

For each `CHEBI:<id>` in `CHEBI_ID_Step2`, Step 3 computes an **equivalence closure** by traversing a set of RO

relationships in the ontology graph, up to a maximum BFS depth (`max_depth`, default 20). 

Default RO equivalence relationships used:

- `RO:0018033` — is conjugate base of
- `RO:0018034` — is conjugate acid of
- `RO:0018036` — is tautomer of
- `RO:0018039` — is enantiomer of 

The output equivalents:

- include only `CHEBI:` IDs
- exclude the starting ID itself
- are cached per input ChEBI ID to speed up repeated runs 

### 3) KEGG enrichment from expanded ChEBI

After RO expansion:

- `CHEBI_ID_Step3 = CHEBI_ID_Step2 ∪ CHEBI_RO_Equivalents`
- derive KEGG compounds from ChEBI xrefs for `CHEBI_ID_Step3`
- `KEGG_ID_Step3 = KEGG_ID_Step2 ∪ derived_kegg`

### 4) Pathway mapping (Step3)

Step 3 loads default mapping resources via `PipelineContext` if a context is not provided:
- Reactome mapping: `source="download"`, species `"Homo sapiens"`
- KEGG mapping: `source="bundled"` 

Then it uses the shared mapper to write Step3 pathway columns with:

- `dedup_pathways=False` (duplicates preserved)
- `dedup_compounds=True` 

## Usage

```python

from pathlib import Path
from step3 import Step3RO

step3 = Step3RO(repo_root=Path("."))

df_step3 = step3.run(
	df_step2,  # must contain CHEBI_ID_Step2
)



