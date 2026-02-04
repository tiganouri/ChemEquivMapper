This folder contains the preprocessing step of the pipeline. It is responsible for reading a feature table, canonicalizing feature names, aligning samples with metadata, and producing consistent feature-level and metabolite-level outputs for downstream analysis steps.



## Role in the pipeline

The preprocess step converts a raw feature table + sample metadata into:

- a clean **feature × sample** abundance matrix with stable unique feature identifiers
- a **feature metadata** table (canonical names, original row numbers)
- **metabolite-level statistics** by grouping duplicate/canonical feature names and applying a hybrid significance strategy suitable for downstream ORA / reporting


## Public API



This module exposes the following objects:

- `InputTableConfig`, `InputTablePreprocessor`
- `PreprocessConfig`, `PreprocessResult`, `Preprocessor`

## Inputs

### Feature table



A file path to a table where:

- **rows** represent features
- **columns** represent samples
- one column contains feature names (by default the first column, or a named column via configuration) 



Supported formats:

- `.csv`, `.tsv`, `.xlsx` / `.xls` (and a best-effort delimiter sniff for other text tables). 

### Sample metadata (`sample_meta`)

A pandas DataFrame where:

- `**sample_meta**.**index**` contains **sample IDs**
- a column `group_col` contains **group labels** for each sample

Important constraints:

- Sample IDs must overlap between feature-table columns and `sample_meta.index`. If there is no overlap, preprocessing raises an error.
- The current “hybrid” strategy requires ****exactly 2 groups****. 


## Outputs

The main entrypoint returns a `PreprocessResult` with:



### 1) `abundance_feature`

A DataFrame of shape **(n_features × n_samples)**:

- index = unique `feature_id`
- columns = sample IDs



### 2) `feature_meta`

A DataFrame indexed by `feature_id` with:
- `feature_id` (unique)
- `feature_name` (canonicalized; duplicates allowed)
- `row_number` (optional)



### 3) `metabolite_stats`

A metabolite-level stats table (metabolite_id is the canonical feature name grouping) including:

- `p_value`
- `p_adj` (BH-FDR)
- `significant` (based on configured alpha)
- `method_used` (single / aggregate / fisher)		
- `n_features`
- `n_per_group` 

### 4) `scale_report`

For each metabolite group:

- number of features
- min/max non-zero median
- median ratio and log10 ratio
- a boolean flag indicating whether feature scales are considered “similar” (ratio <= threshold) 


### 5) `metabolite_feature_map`

A mapping table:
- `metabolite_id`
- `feature_ids_used` (list of feature_ids contributing to that metabolite group) 

## What the preprocess step does

### 1) Canonicalize feature names and create stable unique `feature_id`

Feature names are normalized by trimming and collapsing whitespace, then used to build:

- `feature_name` (canonical; can repeat)
- `feature_id` (unique; duplicates receive a suffix like `__2`, `__3`, …) 

Optional duplicate reporting is printed to stdout when duplicates are detected (duplicates are kept).

The abundance table is coerced to numeric values only (`errors="coerce"`).

### 2) Align samples with metadata

The preprocessor:

- keeps only samples present in both the feature table columns and `sample_meta.index`
- warns about dropped sample columns or metadata rows that do not match

It validates that:

- `group_col` exists in `sample_meta`
- no missing group labels are present



### 3) Group features by metabolite_id (canonical name)

Metabolite groups are formed by using the canonical `feature_name` and collecting all corresponding `feature_id`s.


### 4) Hybrid p-value strategy (per metabolite group)



For each metabolite group, the method depends on how many features are present and whether the feature scales are similar.



Scale similarity is determined using non-zero medians across features:

- `median_ratio = max(median_nonzero) / min(median_nonzero)`

- “similar” if `median_ratio <= ratio_thresh` 



Then:

- **1 feature** → compute p-value for that feature (`single`)

- **>1 feature + similar scales** → aggregate feature abundances by mean across features, then test (`aggregate`)

- **>1 feature + different scales** → compute per-feature p-values and combine using Fisher’s method (`fisher`)



Statistical tests supported:

- 2-sample t-test (`ttest`, default)

- Mann–Whitney U (`mwu`)



Multiple testing correction:

- BH-FDR is applied to metabolite p-values to produce `p_adj`, and `significant` is set based on `alpha`. 



## Usage



### Standalone usage



```python

import pandas as pd

from preprocess import Preprocessor, PreprocessConfig



# sample_meta must have **sample IDs in the index** and a group column

# sample_meta = pd.read_csv("sample_meta.csv", index_col=0)



cfg = PreprocessConfig(

&nbsp;   feature_col=None,      # use first column as feature name column

&nbsp;   test="ttest",          # or "mwu"

&nbsp;   ratio_thresh=10.0,

&nbsp;   alpha=0.05,

&nbsp;   min_n_per_group=2,

)



pp = Preprocessor(cfg)

res = pp.run(
	feature_table_path="feature_table.xlsx",

	sample_meta=sample_meta,

	group_col="Group",

)



abundance_feature = res.abundance_feature

feature_meta = res.feature_meta

metabolite_stats = res.metabolite_stats

scale_report = res.scale_report

metabolite_feature_map = res.metabolite_feature_map

