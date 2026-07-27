# Preprocessing

Feature-table preprocessing utilities for metabolomics-style datasets.

This module provides a **robust, reproducible preprocessing pipeline** for
feature tables with rows = features (e.g. metabolites) and columns = samples.

It is designed for **multi-group experimental designs** and explicitly
separates *feature filtering* from *imputation*.

---

## What the preprocessing does

The preprocessing pipeline applies the following steps **in order**:

1. **Flexible file reading**
   - Automatically detects separators and encodings for `.txt`, `.tsv`, and `.csv` files.

2. **Metadata normalization**
   - Standardizes sample identifiers.
   - Ensures a `group` column exists (missing values are assigned to a configurable label).

3. **Group-aware sample alignment**
   - Retains only samples present in both the feature table and metadata.
   - Subsets the dataset to user-defined groups of interest.

4. **Groupwise missing-value filtering (feature-level)**
   - A feature is **retained if at least one group** has  
     ≥ 50% **non-missing values** (i.e. missing fraction ≤ 0.50).
   - Features failing this criterion are **removed entirely**.
   - **No imputation is performed before this filtering step.**

5. **Numeric coercion**
   - All sample values are coerced to numeric.
   - Non-numeric values are converted to `NaN`.

6. **Missing-value imputation**
   - **Primary method:** KNN-based imputation (sample-wise, across features).
   - **Fallback rule:**  
     If KNN imputation fails *or leaves missing values*, remaining `NaN`s are
     imputed **per feature** using:

     ```
     min(non-missing values of the feature) / 1000
     ```

   - The fallback is applied **only after** the groupwise missing-value filter
     and only to features that passed the filtering step.

---

## Design principles

- **Feature-centric processing**  
  Rows represent features; columns represent samples.

- **Separation of concerns**
  - Feature *validity* is decided before imputation.
  - Imputation is used only to stabilize retained features.

- **Group-aware robustness**
  - Features specific to a single biological group are preserved.
  - Sparse or uninformative features are removed early.

- **Numerical stability**
  - The fallback imputation preserves feature-specific scale and avoids zeros,
    which is critical for downstream log-transforms and multivariate analysis.

---

## Main Components

### `DatasetPreprocessConfig`

Configuration container controlling preprocessing behavior.

```python
from ChemEquiv.exp_data_preprocessing import DatasetPreprocessConfig
