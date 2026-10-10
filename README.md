---
title: Retina scRNA-seq Pattern Viewer
emoji: 🧬
colorFrom: indigo
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
---

# Retina scRNA-seq Pattern Viewer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23282649.svg)](https://doi.org/10.5281/zenodo.23282649)

An interactive [Plotly Dash](https://dash.plotly.com/) application for exploring
single-cell RNA-seq data with spatial/topographic structure from the developing
retina. It accompanies a manuscript on topographic gene expression in the early
retina (the bioRxiv preprint is linked at the top of the app). Two views are
available:

- **Embedding / spatial map** — color cells on a UMAP or on a **custom embedding
  built from two metadata score axes** (e.g. dorsoventral vs. nasotemporal), by a
  metadata field or by the expression of a searchable gene. Optionally switch to a
  binned ("metacell") view, a per-bin **"% positive (detected)"** map, a **two-gene**
  side-by-side comparison, or the **whole-mount ("flower") reprojection** that warps
  the flat DV/NT plane into a flat-mounted retina (with adjustable projection
  parameters), or a **3D sphere projection**. A selectable **HAA pointer**
  marks the chick high-acuity landmark.
- **Expression by group** — the distribution of a gene's expression (or **any
  continuous metadata variable**) across the categories of an `.obs` column, as a
  violin / box / strip, a gene-set dot plot, or a pseudobulk + per-cell figure.

- **ROI differential expression (exploratory)** — draw polygon ROIs on the raw
  topographic map, or load manuscript figure-region presets (chick Fig 6H,
  human Fig S24), then run pseudobulk DE and inspect a volcano and downloadable table.
- **Mask-normalised smoothing (opt-in)** — the default zero-fill Gaussian smoothing
  matches the manuscript pipeline: invalid bins contribute zero, reducing signal near
  sparse edges. The opt-in mode divides `smooth(values × mask)` by `smooth(mask)`;
  an isolated 100% detected bin stays 100%. Both modes keep masked bins masked and
  use the selected sigma; the HAA pointer uses the same smoothed field.

Views (including all of the above controls) can be shared via an encoded URL, and
each dataset's `.h5ad` is available for download.

## Datasets

The datasets served are declared in [`datasets_config.yml`](datasets_config.yml).
The current configuration ships three cross-species retinal progenitor cell
(RPC) datasets (chick, human, mouse), each with inferred 2D topographic
coordinates (`DV.Score`, `NT.Score`); a full chick retina object spanning all
cell classes (UMAP-navigable, no topographic scores); and the prior mouse
alignment retained as a labeled **legacy** entry:

| ID | Title | Cells | Source / attribution |
|----|-------|-------|----------------------|
| `chick_rpc` | Chick Retinal Progenitor Cells | 29,025 | This study (Cepko Lab), GEO **GSE322831** (public), + 1 reprocessed public library, GEO **GSE142244** (Ghinia Tegla et al. 2020) |
| `chick_full` | Chick Retina — Full (all cell classes) | 85,135 | This study (Cepko Lab), GEO **GSE322831** (public), + 1 reprocessed public library, GEO **GSE142244** (Ghinia Tegla et al. 2020); UMAP-navigable, no DV/NT scores (those apply to `chick_rpc`) |
| `human_rpc` | Human Retinal Progenitor Cells | 21,793 | Reprocessed from GEO **GSE138002** (Lu et al. 2020), **GSE234963** (Dorgau et al. 2024), **GSE246169** (Wohlschlegel et al. 2023) |
| `mouse_rpc` | Mouse Retinal Progenitor Cells | 26,505 | Reprocessed (Cell Ranger 9.0.1 / GRCm39) from GEO **GSE118614** (Clark et al.), **GSE139904** (Balasubramanian et al., control cells only), **GSE149040** (Wu et al.), **GSE122466** (Lo Giudice et al.) |
| `mouse_rpc_legacy` | Mouse Retinal Progenitor Cells (legacy) | 25,202 | **Superseded** — original four-library alignment (older reference) from GEO **GSE139904**, **GSE118614**, which pooled wild-type + Fgfr1/2-mutant cells. Retained only for reproducibility; not for new analysis. |

Human and mouse datasets are reprocessed entirely from publicly available GEO
series; chick data is public at GEO **GSE322831**, with one additional
reprocessed public library. The mouse
dataset was re-aligned (Cell Ranger 9.0.1 / GRCm39-2024-A) to de-contaminate a
mixed-genotype source library; the prior alignment is kept as `mouse_rpc_legacy`
for reproducibility and should not be used for new analysis. Please cite the
original accessions above when reusing these data.

> The `.h5ad` data files themselves are **not** committed to this repository
> (they are git-ignored and can be large). See [Data](#data-not-bundled) below
> for how they are provided at runtime.

To add a dataset, append an entry to `datasets_config.yml` with a `title`,
`description`, `file_path` (relative to the data directory, e.g.
`data/my_dataset.h5ad`), `last_updated`, and optional `metadata`. These four
keys (`title`, `description`, `file_path`, `last_updated`) are required; all other
per-dataset keys are optional. Use a `sha256` checksum to pin new data. A leading
`data/` is stripped before joining under `DATA_DIR`; nested paths are preserved
(e.g. `data/sub/x.h5ad` → `$DATA_DIR/sub/x.h5ad`). Basenames must be unique. The app loads
embeddings (`adata.obsm`), metadata columns (`adata.obs`) and genes
(`adata.var_names`) dynamically from each file -- no code changes are required
to add another dataset (any species; a species switch first tries a case-insensitive gene-name match).

## Configuration

Behavior is controlled by environment variables (all optional):

| Variable | Default | Purpose |
|----------|---------|---------|
| `DASH_DEBUG` | `false` | Set to `true` to enable Dash debug mode + hot reload. Keep `false` in production. |
| `DATA_DIR` | `data` | Directory holding the `.h5ad` files. Dataset paths in `datasets_config.yml` are resolved relative to this directory. |
| `PORT` | `7860` | Port for local and Docker entrypoints (Spaces expects `7860`). |
| `HOST` | `0.0.0.0` | Bind address for local and Docker entrypoints. |
| `HF_DATA_REPO` | _(unset)_ | If set (e.g. `username/retina-scrnaseq-data`), provision configured data lazily and pre-warm in the background. No-op when unset. |
| `HF_DATA_REVISION` | Configured immutable revision | Override the top-level `hf_data_revision` when fetching data. Configured SHA256 checks still apply. |
| `DEG_N_CPUS` | `4` (capped by available CPUs) | Threads used for pseudobulk DE fitting. One DE fit runs process-wide. |
| `HF_DATA_REPO_TYPE` | `dataset` | Type of the `HF_DATA_REPO` (`dataset`, `model`, or `space`). |
| `HF_TOKEN` | _(unset)_ | Token for private data and for `scripts/deploy_space.py`. Set it as a Space *secret*, never in code. |

Logging goes to stderr by default. `python app.py -debug` enables debug-level
logging and `app_debug.log` (rotated at 5 MB, three backups); this is separate from
`DASH_DEBUG`, which enables Dash's development server features.

Per-dataset display knobs (`smooth_sigma`, `min_cells_per_bin`, `color_floor`,
`gene_modules`, `default_gene`, `annotation_colors`, …) live in
[`datasets_config.yml`](datasets_config.yml); other plot defaults are set in
`utils/plotting.py`.

## Data (not bundled)

The `.h5ad` files are intentionally **not** part of the image or repo. There are
two supported ways to provide them at runtime:

1. **Local / mounted:** point `DATA_DIR` at a directory containing the files
   referenced in `datasets_config.yml`. With the default config that is:

   ```
   $DATA_DIR/20250604_chick_RPC.h5ad
   $DATA_DIR/20260620_chick_full.h5ad
   $DATA_DIR/20250604_human_RPC.h5ad
   $DATA_DIR/20260528_mouse_RPC_cr9_e13e16.h5ad
   $DATA_DIR/20250604_mouse_RPC.h5ad          # legacy (superseded)
   ```

2. **Hugging Face repo (used on Spaces):** set `HF_DATA_REPO` to a Hugging Face
   dataset repo holding those files. At startup the app downloads any missing
   files into `DATA_DIR` (see [`utils/data_provision.py`](utils/data_provision.py)).
   For a private data repo, also set `HF_TOKEN` as a Space secret. The config pins
   revision `2c21052aad5203a64b5798474ed2045fca8e3733` and each shipped file's
   SHA256. Downloads and existing local files are verified by streaming SHA256;
   mismatches fail and the bad file is deleted. Verification stamps keyed by
   path, size and mtime are cached beside each file for fast restarts.

## Run locally

```bash
# 1. Create an environment and install dependencies
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. Put your .h5ad files where datasets_config.yml expects them
#    (default: ./data/), or set DATA_DIR to point elsewhere.
mkdir -p data
# cp /path/to/*.h5ad data/

# 3. Run (defaults to http://0.0.0.0:7860)
python app.py

# Optional: enable debug + a different port
DASH_DEBUG=true PORT=8050 python app.py
```

Then open the printed URL in a browser.

## Deploy to Hugging Face Spaces

This app is set up as a **Docker** Space (see the YAML frontmatter at the top of
this README, which Hugging Face reads, and the [`Dockerfile`](Dockerfile)).

1. Create a new Space at <https://huggingface.co/new-space>, choosing **Docker**
   as the SDK. Start it **private** for testing if you like.
2. From a **clean, committed checkout**, run:

   ```bash
   HF_TOKEN=... python scripts/deploy_space.py owner/space-name
   ```

   The script refuses a dirty git tree, exports tracked files with `git archive`,
   stamps the source SHA in `viewer_revision.txt`, and uploads the checkout with
   `huggingface_hub.upload_folder`. The footer displays `viewer <version> · <short sha>`.
   The Space is deployed by this upload, **not by git push**. Hugging Face builds
   the Dockerfile and serves the app on port `7860`.
3. Provide the data via a separate Hugging Face **dataset** repo: upload the
   `.h5ad` files there, then set the Space variables `HF_DATA_REPO`
   (e.g. `username/retina-scrnaseq-data`) and, while that repo is private,
   `HF_TOKEN` (as a *secret*). The `Dockerfile` does **not** copy any data into
   the image.
4. Leave `DASH_DEBUG` unset (or `false`) for the public deployment.

## Project layout

```
app.py                 # Dash app, layout, core callbacks, /download route, entrypoint
datasets_config.yml    # Declares the datasets served (titles, file paths, metadata)
requirements.txt       # Python dependencies
Dockerfile             # Hugging Face Docker Space build
callbacks/             # Dash callbacks (dataset, selection, URL, plotting, status, ROI/DEG)
components/            # Reusable UI pieces (sidebar, status bar)
layouts/               # Page layout (sidebar + main panel)
utils/                 # Data loading, HF provisioning, plotting, state, DE, transforms
docs/                  # Viewer specification and hardening plans
assets/                # Browser ROI drawing code
scripts/               # Clean-checkout Space deployment
tests/                 # Synthetic pytest regressions (no .h5ad files)
```


## Testing

```bash
pip install pytest && pytest tests/
```

Tests use synthetic data and mocked readers/downloads; no `.h5ad` files are needed.

## Citation

Please cite the [preprint](https://doi.org/10.64898/2026.01.04.697548) by Heer N. V.
Joisher, ChangHee Lee, Chaitra Prabhakara, Isabella van der Weide, Yichen Si,
Nicholas Lonfat and Constance Cepko, and the original GEO accessions above.
Machine-readable citation metadata is in [CITATION.cff](CITATION.cff).
The [public analysis repository](https://github.com/chlee-tabin/retina-spatial-scrna-analysis)
contains `scripts/figures/sfig24_area_deg.R` and
`scripts/figures/supptables1_2_area_deg.R`; see manuscript Methods for the
expression-by-group figure. Viewer DE uses `~ library + condition` when estimable
and the configured replicate unit. Its `min_frac=0.5` detection filter is a
viewer-only addition that changes adjusted p-values versus manuscript tables;
the interactive results are exploratory.

## License (MIT)

Viewer code is licensed under [MIT](LICENSE). Cite the sources above when reusing data.
