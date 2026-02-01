# Resource managers (Reactome / KEGG / ChEBI)

This folder provides small **resource resolver** utilities that locate (or fetch) reference files required by downstream steps (e.g., metabolite-to-pathway mapping). Each resolver returns a **local file path** to a usable resource, based on a configurable strategy: use a bundled file, use a local file, or download-and-cache (with fallback to bundled).

The module exposes:

- `ReactomeTXTConfig`, `resolve_reactome_txt`
- `KEGGTSVConfig`, `resolve_kegg_tsv`
- `ChebiResourceConfig`, `resolve_chebi_obo`

---

## Role in the pipeline

Many pipelines need curated/reference mapping resources:

- **Reactome**: ChEBI → Reactome pathways (TXT mapping)
- **KEGG**: compound → pathway links (TSV mapping)
- **ChEBI**: ontology file (OBO gz) used for ID resolution / normalization

These utilities provide a consistent way to obtain these files so that the pipeline can be:

- reproducible (bundled mode),
- up-to-date when needed (download mode),
- or controlled by the user (local mode).

---

## Public API

The public API is exported from this package’s `__init__.py`:

- `ReactomeTXTConfig`, `resolve_reactome_txt`
- `KEGGTSVConfig`, `resolve_kegg_tsv`
- `ChebiResourceConfig`, `resolve_chebi_obo`

---

## Common concepts

### Resolution strategy

Each resource uses a `source` field (a `Literal`):

- `"bundled"`: use a file shipped inside the Python package (recommended for reproducibility)
- `"local"`: use a user-provided path, with fallback to bundled
- `"download"`: download into a cache directory, with fallback to bundled (where supported)

All resolvers:

- create the cache directory if needed,
- **never raise** on failure (they return `None` if nothing usable is available).

---

### Bundled resources inside the package

Bundled files are resolved using `importlib.resources` and `as_file(...)` to produce a concrete filesystem `Path`.

By default, the package name and relative paths are:

- package: `GEPC`
- ChEBI: `data/chebi/chebi.obo.gz`
- KEGG: `data/kegg/kegg_compound_to_pathway.tsv`
- Reactome: `data/reactome/ChEBI2Reactome_All_Levels.txt`

If you rename the package, update `package_name` in the configs.

---

## Reactome resource manager

### What it provides

- `ReactomeTXTConfig`
- `resolve_reactome_txt(cfg, cache_dir, download_url=None) -> Optional[Path]`

---

### Default behavior

Reactome defaults to `source="download"` and attempts to download the Reactome mapping file, caching it locally, and falling back to the bundled file if download fails.

Default download URL (as currently implemented):

- `https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt`

---

### Example usage

```python
from pathlib import Path
from resources import ReactomeTXTConfig, resolve_reactome_txt

cfg = ReactomeTXTConfig(source="download")
path = resolve_reactome_txt(cfg, cache_dir=Path(".cache/resources"))

if path is None:
    raise RuntimeError("Reactome mapping not available")
