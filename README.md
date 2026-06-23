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
  parameters).
- **Expression by group** — the distribution of a gene's expression (or **any
  continuous metadata variable**) across the categories of an `.obs` column, as a
  violin / box / strip, a gene-set dot plot, or a pseudobulk + per-cell figure.

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
| `chick_rpc` | Chick Retinal Progenitor Cells | 29,025 | This study (Cepko Lab) + 1 reprocessed public library, GEO **GSE142244** (Emerson et al.) |
| `chick_full` | Chick Retina — Full (all cell classes) | 85,135 | This study (Cepko Lab) + 1 reprocessed public library, GEO **GSE142244** (Emerson et al.); UMAP-navigable, no DV/NT scores (those apply to `chick_rpc`) |
| `human_rpc` | Human Retinal Progenitor Cells | 21,793 | Reprocessed from GEO **GSE138002** (Sridhar et al.), **GSE234963**, **GSE246169** |
| `mouse_rpc` | Mouse Retinal Progenitor Cells | 26,505 | Reprocessed (Cell Ranger 9.0.1 / GRCm39) from GEO **GSE118614** (Clark et al.), **GSE139904** (Balasubramanian et al., control cells only), **GSE149040** (Wu et al.), **GSE122466** (Lo Giudice et al.) |
| `mouse_rpc_legacy` | Mouse Retinal Progenitor Cells (legacy) | 25,202 | **Superseded** — original four-library alignment (older reference) from GEO **GSE139904**, **GSE118614**, which pooled wild-type + Fgfr1/2-mutant cells. Retained only for reproducibility; not for new analysis. |

Human and mouse datasets are reprocessed entirely from publicly available GEO
series; chick data is in-house except one reprocessed public library. The mouse
dataset was re-aligned (Cell Ranger 9.0.1 / GRCm39-2024-A) to de-contaminate a
mixed-genotype source library; the prior alignment is kept as `mouse_rpc_legacy`
for reproducibility and should not be used for new analysis. Please cite the
original accessions above when reusing these data.

> The `.h5ad` data files themselves are **not** committed to this repository
> (they are git-ignored and can be large). See [Data](#data-not-bundled) below
> for how they are provided at runtime.

To add a dataset, append an entry to `datasets_config.yml` with a `title`,
`description`, `file_path` (relative to the data directory, e.g.
`data/my_dataset.h5ad`), `last_updated`, and optional `metadata`. The app loads
embeddings (`adata.obsm`), metadata columns (`adata.obs`) and genes
(`adata.var_names`) dynamically from each file -- no code changes are required
to add another dataset (any species; gene-name casing is taken verbatim from
`var_names`).

## Configuration

Behavior is controlled by environment variables (all optional):

| Variable | Default | Purpose |
|----------|---------|---------|
| `DASH_DEBUG` | `false` | Set to `true` to enable Dash debug mode + hot reload. Keep `false` in production. |
| `DATA_DIR` | `data` | Directory holding the `.h5ad` files. Dataset paths in `datasets_config.yml` are resolved relative to this directory. |
| `PORT` | `7860` | Port the server listens on (Hugging Face Spaces expects `7860`). |
| `HOST` | `0.0.0.0` | Bind address (`0.0.0.0` so the container is reachable). |
| `HF_DATA_REPO` | _(unset)_ | If set (e.g. `username/retina-scrnaseq-data`), download the configured datasets from this Hugging Face repo into `DATA_DIR` at startup. No-op when unset. |
| `HF_DATA_REPO_TYPE` | `dataset` | Type of the `HF_DATA_REPO` (`dataset`, `model`, or `space`). |
| `HF_TOKEN` | _(unset)_ | Access token, only needed while `HF_DATA_REPO` is **private**. Set it as a Space *secret*, never in code. |

Plot defaults (colors, marker size, bin size, etc.) live in
[`config.yml`](config.yml).

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
   For a private data repo, also set `HF_TOKEN` as a Space secret.

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
2. Push this repository to the Space's git remote. Hugging Face builds the image
   from the `Dockerfile` and serves the app on port `7860` (matching `app_port`
   in the frontmatter).
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
config.yml             # Plot/layout defaults
requirements.txt       # Python dependencies
Dockerfile             # Hugging Face Docker Space build
callbacks/             # Dash callbacks (dataset, selection, URL, plotting, status)
components/            # Reusable UI pieces (sidebar, status bar)
layouts/               # Page layout (sidebar + main panel)
utils/                 # Data loading, HF data provisioning, plotting, state, config
```
