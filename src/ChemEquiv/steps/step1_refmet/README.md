# Step 1 — RefMet annotation + pathway mapping

This folder implements **Step 1** of the pipeline. It takes metabolite-level metadata (one row per canonical metabolite name), annotates metabolite names using **RefMet** (Metabolomics Workbench), caches results locally in SQLite to avoid repeated network calls, and maps RefMet-derived identifiers to **Reactome** and **KEGG** pathways.

The package exports `Step1RefMet`.

---

## Role in the pipeline

Upstream preprocessing produces a metabolite-level table (e.g. `metabolite_stats`) indexed by a **unique, canonical metabolite name**.  
Step 1 enriches that table with:

1. **RefMet standardisation / annotation**
   - RefMet ID
   - Formula
   - ChEBI ID
   - KEGG ID
   - etc.  
   via a POST request with retry logic and SQLite caching.

2. **Pathway mapping**
   - ChEBI → Reactome pathways
   - KEGG compound → KEGG pathways  

   These are written as **Step1-specific columns** in the output DataFrame.

> **Important contract**  
> Step 1 does **not** read feature tables and does **not** handle duplicates.  
> These responsibilities are handled upstream during preprocessing.

---

## Public API

- **`Step1RefMet`** — main entry point

Supporting modules:

- `RefMetClient`  
  POST client for RefMet with retries and backoff
- `SQLiteRefMetCache`  
  SQLite-backed cache keyed by normalized input name
- `RefMetSchema`  
  Normalizes RefMet output into stable canonical columns
- `RefMetZipConfig`, `resolve_refmet_zip`  
  Optional utilities for zip resolution / downloading

---

## Inputs and outputs

### Input: metabolite-level metadata table

`Step1RefMet.run(...)` expects a pandas `DataFrame` that:

- has **one row per metabolite**
- is indexed by a **unique canonical metabolite name**

---

### Output: enriched metadata table

The returned DataFrame:

- preserves the same index
- adds:
  - RefMet annotation columns (canonical schema)
  - Step1 mapping columns for Reactome and KEGG pathways

---

## What Step 1 does

### 1) Normalize metabolite names

Metabolite names are normalized by:

- trimming leading/trailing whitespace
- collapsing internal whitespace

before being sent to RefMet.

---

### 2) RefMet POST annotation (with caching)

Annotations are fetched by POSTing newline-separated metabolite names to the RefMet endpoint.

**Caching behavior:**

- SQLite cache:
  - **key**: normalized input name (lowercased, whitespace-collapsed)
  - **value**: JSON representation of the RefMet row
- Only **missing** names are sent to the network
- Requests are chunked by `chunk_size`
- Results are written back to the cache
- Output is normalized to stable columns via `RefMetSchema`

---

### 3) Merge annotations back into the input table

Step 1 ensures that all canonical RefMet columns exist:

- missing columns are created and filled with `NA`
- values are filled per metabolite index name

---

### 4) Map IDs to pathways (Step1 columns)

Using RefMet-derived identifier columns (`ChEBI_ID`, `KEGG_ID`), Step 1 performs:

- ChEBI → Reactome pathway mapping
- KEGG compound → KEGG pathway mapping

and writes results into Step1-specific output columns using the shared pipeline mapper.

**Conservative mapping rules:**

- ChEBI ID cells are normalized into comma-separated `CHEBI:<id>` tokens
- KEGG ID cells are kept as conservative strings  
  (the mapper supports comma-separated tokens)

---

## Usage

### Minimal example

```python
from pathlib import Path
import pandas as pd

from step1_refmet import Step1RefMet

# metabolite_stats: index must be unique canonical metabolite name
# metabolite_stats = pd.read_csv("metabolite_stats.csv", index_col=0)

step1 = Step1RefMet(repo_root=Path("."))

out = step1.run(
    metabolite_stats,
    unzip=True,                  # optional: resolves/unzips RefMet zip for compatibility
    reactome_download_url=None,  # optional override
)
