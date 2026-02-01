# Step 5 — is_a Generalisation for previously unmapped rows (Reactome-only)

## Purpose

Step 5 is a **recovery step** that runs **only on rows that remain unmapped to Reactome after Steps 1–4**.  
For those rows, it attempts to recover Reactome pathways by traversing **is_a parents upward in the ChEBI ontology** and identifying the **closest Reactome-mapped ancestors**.

---

## Key features

- Reactome-only (no KEGG work in this step)
- No RO-neighbor expansion
- Finds **ALL closest mapped ancestors** (minimum BFS distance)
- Produces:
  - new dataframe columns with generalisation results
  - an optional structured report that can be exported to PDF

---

## When Step 5 runs (row eligibility)

A row is considered **already mapped** if **any** of the following columns is non-empty:

- `Reactome_Pathways_Step1`
- `Reactome_Pathways_Step2`
- `Reactome_Pathways_Step3`
- `Reactome_Pathways_Step4`

Step 5 **only processes rows where all of the above are empty**.

---

## Which ChEBI IDs Step 5 starts from

For each unmapped row, Step 5 uses **only**:

- `CHEBI_ID_Step2`

This cell is:

- parsed into a list of `CHEBI:<id>` tokens
- each token is used as an independent starting point

If `CHEBI_ID_Step2` is missing or empty, the row is **skipped** by Step 5.

---

## Core method: closest mapped ancestors via is_a BFS

For each base ChEBI ID, Step 5 performs a **breadth-first search (BFS)** upward through `is_a` parents and:

- searches up to `max_is_a_depth` (default: `20`)
- stops once the **minimum depth** with mapped ancestors is found
- returns **all mapped ancestors at that same minimum depth**
- prunes overly-generic ancestors via `excluded_ancestors`  
  (default exclusions include:
  - `CHEBI:23367` — *molecular entity*
  - `CHEBI:24431` — *chemical entity*)

An ancestor is considered **mapped** if it has **at least one Reactome pathway** in the shared Reactome index  
(ChEBI → Reactome pathways).

---

## Outputs (new columns added)

Step 5 adds or fills the following columns:

- **`CHEBI_IsA_Generalisation`**  
  Comma-separated list of the closest mapped ancestor ChEBI IDs (sorted).

- **`Reactome_Pathways_IsA_Generalisation`**  
  Comma-separated list of Reactome pathway IDs obtained from those ancestors  
  (unique and sorted).

- **`n_pathways_isa_generalisation`**  
  Integer count of pathways assigned by Step 5.

These values are written **only if** Step 5 finds:
- at least one mapped ancestor **and**
- at least one Reactome pathway.

---

## Report and PDF export (optional)

If `generate_report=True`, Step 5 produces a list of `report_entries` that can be exported to PDF.

Each report entry includes:

- `row_index`
- `base_chebi`
- `path_ids` / `path_labels` describing the `is_a` traversal route
- `mapped_chebi`
- `mapped_chebi_name` (if available)
- `pathways` (list of Reactome pathways)

---

## Usage

### Minimal example

```python
from pathlib import Path
from step5 import Step5IsAGeneralisation, save_generalisation_report_pdf

step5 = Step5IsAGeneralisation(repo_root=Path("."))

df_step5, report_entries = step5.run(
    df_step4,
    generate_report=True
)

# Optional: save PDF report
save_generalisation_report_pdf(
    report_entries,
    pdf_path=Path("step5_is_a_generalisation_report.pdf"),
    dataset_label="MyDataset"
)
