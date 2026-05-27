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
single-cell RNA-seq data with spatial/topographic structure from developing
retina. It accompanies a manuscript on topographic gene expression in the early
retina. Users can select a dataset, choose an embedding (UMAP or a **custom
embedding built from two metadata score axes**, e.g. dorsoventral vs.
nasotemporal), color cells by a metadata field or by the expression of a
searchable gene, and optionally switch to a binned ("metacell") density view.
Views can be shared via an encoded URL, and each dataset's `.h5ad` is available
for download.

## Datasets

The datasets served are declared in [`datasets_config.yml`](datasets_config.yml).
The current configuration ships three early-chick retina datasets:

| ID | Title | Description |
|----|-------|-------------|
| `full_dataset` | Chick Full Retinal Dataset | Full early chick retinal dataset (Cepko lab + reprocessed Emmerson lab data) |
| `fabp7_dataset` | Subset: Chick FABP7+ RPCs | FABP7+ retinal progenitor cells with topographic layout |
| `otx2neg_dataset` | Subset: Chick OTX2- Proliferative RPCs | Proliferative RPCs lacking OTX2 expression |

> The `.h5ad` data files themselves are **not** committed to this repository
> (they are git-ignored and can be large). See
> [Data](#data-not-bundled) below for how to provide them at runtime.

To add a dataset, append an entry to `datasets_config.yml` with a `title`,
`description`, `file_path` (relative to the data directory, e.g.
`data/my_dataset.h5ad`), `last_updated`, and optional `metadata`. The app loads
embeddings (`adata.obsm`), metadata columns (`adata.obs`) and genes
(`adata.var_names`) dynamically from each file -- no code changes are required
to add another chick dataset. (Adding **human/mouse** datasets is discussed in
the PR description / "Still needed before public".)

## Configuration

Behavior is controlled by environment variables (all optional):

| Variable | Default | Purpose |
|----------|---------|---------|
| `DASH_DEBUG` | `false` | Set to `true` to enable Dash debug mode + hot reload. Keep `false` in production. |
| `DATA_DIR` | `data` | Directory holding the `.h5ad` files. Dataset paths in `datasets_config.yml` are resolved relative to this directory, so the data can be mounted/uploaded separately instead of baked into the image. |
| `PORT` | `7860` | Port the server listens on (Hugging Face Spaces expects `7860`). |
| `HOST` | `0.0.0.0` | Bind address (`0.0.0.0` so the container is reachable). |

Plot defaults (colors, marker size, bin size, etc.) live in
[`config.yml`](config.yml).

## Data (not bundled)

The `.h5ad` files are intentionally **not** part of the image or repo. Provide
them at runtime by pointing `DATA_DIR` at a directory that contains the files
referenced in `datasets_config.yml`. For example, with the default config the
directory should contain:

```
$DATA_DIR/20240815_full.h5ad
$DATA_DIR/20240815_fabp7.h5ad
$DATA_DIR/20240815_otx2negRPC.h5ad
```

On Hugging Face Spaces, upload the files into the Space (e.g. a `data/`
directory in the Space repo, or a mounted dataset) and set `DATA_DIR`
accordingly.

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

1. Create a new Space at <https://huggingface.co/new-space>, choosing
   **Docker** as the SDK.
2. Push this repository to the Space's git remote (or connect the GitHub repo).
   Hugging Face builds the image from the `Dockerfile` and serves the app on
   port `7860` (matching `app_port` in the frontmatter).
3. Provide the data: upload your `.h5ad` files into the Space and set the
   `DATA_DIR` variable (Space *Settings -> Variables and secrets*) if they are
   not in the default `data/` directory. The `Dockerfile` does **not** copy any
   data into the image.
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
utils/                 # Data loading, plotting, state encoding, monitoring, config
```
