# Step 4 — is_a specialization (children) + RO expansion + pathway mapping


This folder implements **Step 4** of the pipeline: starting from Step 3 IDs, it expands ChEBI identifiers by
traversing **is_a children** (more specific terms) in the ChEBI ontology, optionally applies RO equivalence
expansion to those children, then unions the results into `CHEBI_ID_Step4` and derives `KEGG_ID_Step4`
from ChEBI xrefs. Finally, it maps Step4 IDs to Reactome/KEGG pathways using the shared mapper.

Main entrypoint:

- `Step4IsASpecialization`

Exports:
- `Step4IsASpecialization`

## Role in the pipeline

Step 3 produces:

- `CHEBI_ID_Step3` (required)
- `KEGG_ID_Step3` (optional)

Step 4 increases downstream pathway mapping coverage by adding **more specific ChEBI terms** (children)
and their RO equivalents. This is useful when you want “specialized” terms that may map to additional
pathways.


Step 4 outputs:

- `CHEBI_is_a_Equivalents`: only *new* specialized ChEBIs (children and their RO equivalents), excluding Step 3 IDs
- `CHEBI_ID_Step4`: Step 3 ChEBIs ∪ specialized ChEBIs
- `KEGG_ID_Step4`: Step 3 KEGG ∪ KEGG compounds derived from Step4 ChEBIs
- Mapper pathway columns for `step=4` (Reactome + KEGG)

## Public API

This package exports:
- `Step4IsASpecialization`

## Inputs

`Step4IsASpecialization.run(df, ...)` expects a pandas DataFrame with:

Required:
- `CHEBI_ID_Step3` (comma-separated ChEBI IDs)

Optional:
- `KEGG_ID_Step3` (comma-separated KEGG compound IDs)

You can override column names:
- `chebi_col_step3` (default `"CHEBI_ID_Step3"`)
- `kegg_col_step3` (default `"KEGG_ID_Step3"`) 

## Outputs

Step 4 returns a copy of the input DataFrame with Step4 expansion columns added (if missing):

- `CHEBI_is_a_Equivalents`
- `CHEBI_ID_Step4`
- `KEGG_ID_Step4` 

Plus pathway mapping columns written by the shared mapper (`step=4`):

- `ChEBI_ID_Step4_reactome`
- `Reactome_Pathways_Step4`
- `KEGG_Compounds_Step4`
- `KEGG_Pathways_Step4` 

After mapping, the Step4 CHEBI/KEGG union columns are normalized (parse → sort → re-join). 

## What Step 4 does

### 1) Resolve ChEBI ontology + CHEBI→KEGG map (cached)

Step 4 resolves `chebi.obo.gz` using the shared resolver (bundled/local/download). If it cannot be resolved,
Step 4 raises a runtime error. 

It loads the ontology with `pronto.Ontology` and builds a CHEBI→KEGG map (cached in SQLite):

- cache directory: `<repo_root>/.cache/step4_isa/`

- map DB: `chebi_to_kegg.sqlite` 

### 2) Resolve RO equivalence relationships (cached)

Step 4 resolves RO relationship objects from the ontology and stores them for reuse.
If none can be resolved, Step 4 raises an error.

Default RO relationships used:

- `RO:0018033` — is conjugate base of
- `RO:0018034` — is conjugate acid of
- `RO:0018036` — is tautomer of
- `RO:0018039` — is enantiomer of :contentReference[oaicite:14]{index=14}

### 3) is_a traversal (children / more specific terms)


For each Step3 ChEBI ID, Step 4 traverses **subclasses()** (children) up to `max_is_a_depth`
(default 1) and returns discovered `CHEBI:` children IDs excluding the starting IDs.

### 4) RO expansion on is_a children

If any is_a children are found, Step 4 computes RO equivalence closure **starting from the children**,
not from the original Step3 IDs.

It then defines the “specialized” set:
- `isa_specific = (isa_children ∪ RO_equiv(isa_children)) - step3_chebis` 

This set is written into:
- `CHEBI_is_a_Equivalents` 

### 5) Union into Step4 ID columns

Step 4 unions Step3 IDs with the specialized set:
- `CHEBI_ID_Step4 = step3_chebis ∪ isa_specific` 

Then it derives KEGG compounds from the Step4 ChEBI IDs:
- `KEGG_ID_Step4 = step3_keggs ∪ derived_kegg(step4_chebis)` 

### 6) Pathway mapping for Step4

Step 4 uses a shared `PipelineContext` (or loads one by default) and applies the mapper with `step=4`. 

Default resource strategy when `ctx` is not provided:
- Reactome mapping: `source="download"`, species `"Homo sapiens"`
- KEGG mapping: `source="bundled"` 

Mapping call defaults:

- `dedup_pathways=False` (duplicates preserved)
- `dedup_compounds=True` 

### 7) Normalize Step4 union columns

Finally, `Step4Schema.normalize_joined_ids()` normalizes the Step4 union columns by:

- parsing the CSV cell
- sorting IDs
- joining back with commas

It does **not** touch mapper-generated pathway columns. 

## Usage


```python

from pathlib import Path
from step4 import Step4IsASpecialization

step4 = Step4IsASpecialization(repo_root=Path("."))


df_step4 = step4.run(
	df_step3,  # must contain CHEBI_ID_Step3
)



