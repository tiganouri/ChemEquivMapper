# ChemEquiv.analysis

This subpackage implements **stepwise metabolite-centric Over-Representation Analysis (ORA)** for pathway enrichment using **Reactome (via ChEBI)** and **KEGG** identifiers.

It is designed to operate downstream of a preprocessing/statistical pipeline where **metabolite-level significance** has already been computed.

**Core rule**

> All identifiers (ChEBI / KEGG) derived from a metabolite inherit that metabolite’s significance flag.

This allows enrichment analysis to be performed at multiple annotation steps while preserving metabolite-level statistical decisions.

---

## Package structure

ChemEquiv/analysis/
├── init.py
├── stats.py
└── ora.py

---

## Modules

### `stats.py`

Statistical utilities used by ORA:

- `fdr_bh`  
  Benjamini–Hochberg false discovery rate (FDR) correction.

- `fisher_exact_right_tail`  
  One-sided Fisher exact test for enrichment.

- `metabolite_sig_map`  
  Builds a `metabolite_id → {0,1}` significance mapping using an OR rule if duplicates exist.

---

### `ora.py`

Implements the stepwise ORA workflow:

- `StepColumnSpec`  
  Defines which ChEBI and KEGG columns correspond to each annotation step.

- `DEFAULT_STEP_SPECS`  
  Default mapping for Steps 1–4.

- `StepwiseORA`  
  Main class that:
  - builds identifier tables per step
  - propagates metabolite significance to identifiers
  - performs ORA against Reactome and KEGG pathways

---

## Conceptual workflow

1. Inputs:
   - a metadata table (`meta_df`) with metabolites and stepwise annotation columns
   - a statistics table (`metabolite_stats`) with metabolite-level significance

2. For each annotation step:
   - explode ChEBI / KEGG identifiers
   - normalize identifiers
   - propagate metabolite significance to identifiers
   - collapse duplicate identifiers using an OR rule

3. Build pathway → identifier mappings using:
   - `PipelineContext.reactome_index`
   - `PipelineContext.kegg_index`

4. Perform ORA using:
   - right-tailed Fisher exact test
   - BH-FDR correction

---

## Requirements

- pandas
- numpy
- scipy
- statsmodels

A valid `PipelineContext` **must** be provided and must contain:
- `reactome_index`
- `kegg_index`

If either index is missing, ORA will raise an error.

---

## Input data requirements

### `meta_df`

A pandas DataFrame containing:
- a metabolite identifier column (default: `Input name`)
- step-specific identifier columns

Example columns:

| Column name         | Description              |
|---------------------|--------------------------|
| `ChEBI_ID`          | Step 1 ChEBI identifiers |
| `KEGG_ID`           | Step 1 KEGG identifiers  |
| `CHEBI_ID_Step2`    | Step 2 ChEBI identifiers |
| `KEGG_ID_Step2`     | Step 2 KEGG identifiers  |

- Multiple identifiers may be comma-separated.
- Missing step columns are handled gracefully and produce empty results.

---

### `metabolite_stats`

A pandas DataFrame containing:
- `metabolite_id`
- `significant` (0/1)

If a metabolite appears multiple times, significance is propagated using an OR rule.

---

## Basic usage

```python
from ChemEquiv.analysis import StepwiseORA
from ChemEquiv.context import PipelineContext

ora = StepwiseORA(ctx=pipeline_context)

results = ora.run_step(
    meta_df=meta_df,
    metabolite_stats=metabolite_stats,
    step=1,
    meta_join_col="Input name"
)
