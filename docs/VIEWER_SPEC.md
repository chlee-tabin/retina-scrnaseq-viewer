# Retina scRNA-seq Viewer — Specification

A living spec for the Dash/Plotly single-cell viewer. It defines the view model,
the controls, the dataset-config schema, the shareable-state schema (v2), and the
definition of the "Expression by group" figure. **This document is the contract used
for independent (ultracode) review** — each claim below should hold in the code.

## Goals

From a browser, a user must be able to reproduce the core manuscript analyses:

1. **NPY⁺ (positive-cell) expression across stages** — the 2-panel pseudobulk +
   positive-cell violin figure (the "Figure" style of the group view).
2. **Topographic gene-vs-gene comparison** — two genes side-by-side on the same
   embedding, each on its **own** colour scale.
3. **Working violin / distribution representations** across cell-type groups.

## View model

Two top-level views, switched by **View** (`plot-type`):

- **Embedding / spatial map** (`map`): a 2D embedding (UMAP / PCA / Harmony) or the
  custom DV/NT topographic axes, coloured by gene expression or an obs column.
  Supports binned + Gaussian-smoothed spatial maps, a "% positive (detected)" bin
  colour mode, and two-gene side-by-side.
- **Expression by group** (`group`): distribution of a gene (or gene set) across the
  categories of an obs column. Styles: **Figure** (default), Violin, Box, Strip,
  Dot plot.

## Controls

| Control (id) | View | Purpose |
|---|---|---|
| `plot-type` | both | Map vs Expression-by-group |
| `embedding-select` | map | Embedding (or `custom_embedding` = DV/NT) |
| `custom-x-select` / `custom-y-select` | map | Axes for the custom embedding |
| `enable-binning` | map (custom) | Binned "metacell" spatial view |
| `bin-stat` | map (custom) | Bin colour: `mean` expression vs `frac_pos` (% positive) |
| `enable-smoothing` | map (custom) | Gaussian smoothing of the binned map |
| `bin-number-slider` | map (custom) | Number of bins (2–120) |
| `percentile-slider` | map (custom) | Colour-scale percentile cutoff |
| `color-select` | map | Gene expression or an obs column |
| `gene-select` | map | Gene (when colouring by expression) |
| `compare-genes` / `gene-select-2` | map | Two genes side-by-side |
| `compare-shared-scale` | map | Shared (absolute) vs independent colour scales |
| `viz-mode` | map | Plot order (random/asc/desc) |
| `group-gene-select` | group | Gene for figure/violin/box/strip |
| `group-by-select` | group | Categorical column (x-axis) |
| `group-split-select` | group | Optional second grouping (violin/box/strip) |
| `group-style` | group | Figure / Violin / Box / Strip / Dot plot |
| `group-positive-only` | group (figure) | Violin uses only cells with expr > 0 |
| `group-replicate-select` | group (figure) | Pseudobulk replicate unit (Panel A) |
| `gene-module-select` | group (dotplot) | Gene set for the dot plot |
| `group-value-source` / `group-meta-select` | group | Plot value: gene expression vs a continuous `.obs` variable |
| `custom-projection` | map (custom) | Raw axes vs whole-mount (DV/NT flower) |
| `wm-rho-nt` / `wm-rho-dv` / `wm-gap` / `wm-stretch` / `wm-cuts` | map (flower) | Flower extent, relief width, petal stretch, number of cuts |
| `wm-symmetric` / `wm-pole` / `wm-dewarp` / `wm-pow` / `wm-gap-mode` / `wm-gap-frac` | map (flower) | Advanced flower geometry (defaults reproduce Fig R2.5) |

Style-specific group controls are revealed by `group-style`: the figure controls
show only for `figure`, the module selector only for `dotplot`.

## Dataset config schema (`datasets_config.yml`)

Per dataset (all optional unless noted): `title`, `description`, `file_path`
(required), `default_gene`, `default_embedding`, `smooth_sigma`,
`min_cells_per_bin`, `color_floor`, `gene_modules`, `deg_results_path` (reserved,
currently unused), plus:

- **`annotation_column`** — the default categorical column for grouping/colour.
- **`annotation_order`** — biological (maturation) order; applied to the scatter
  legend, the violin/figure x-axis, and the dot-plot x-axis.
- **`annotation_colors`** — `{category: hex}`; the **same** colour for a cell type in
  every view. Unmapped categories fall back to a distinct qualitative palette
  (`plotly.express.colors.qualitative.Dark24`, replacing the pale `Set3`).
- **`replicate_columns`** — columns combined to form the pseudobulk replicate unit
  (e.g. `[library, genotype]`). Shown as ` × `; the internal key joins them with a
  non-printing `US` (`\x1f`) separator so a value containing `|` cannot mis-split.

## Categorical colouring

When colouring/grouping on `annotation_column`, a cell type renders the **same**
colour in the UMAP scatter, the violin, the figure, and the dot-plot x-axis. The map
`color_discrete_map` colours configured values; any unmapped value still draws from
the Dark24 fallback (nothing renders colourless).

## Two-gene comparison (#1)

Two panels over the same embedding, one gene each.

- **Independent colour scales by default**: each panel's colorbar maximum is its own
  gene's max, so a sparse gene (e.g. CYP26C1, well under 1 % positive — 0.5 % in
  chick_full, 0.9 % in chick_rpc) is **not** flattened by a strong one (e.g. FGF8).
  `compare-shared-scale` switches to one absolute scale.
- Expressing cells are drawn **last** (sorted ascending by expression) so a sparse
  gene's positive cells are not hidden under the many non-expressing cells.

## Expression-by-group "Figure" (#4) — NPY 2-panel definition

Plotly reproduction of `scripts/viewer_full_chick/16b_npy_figure.R` (analysis repo,
PR #19). For a gene `g`, a grouping column, and a replicate unit:

- **Panel A — pseudobulk per (group × replicate)**: one dot per replicate;
  `y = log1p( Σ raw / Σ depth × T )`; dot **size increases with n cells**; grey
  crossbar = mean across replicates; replicates with **< 10 cells dropped**.
  - `T` is the dataset's counts-per-X normalisation target, **detected at load** from
    the per-cell `expm1(.X)` sum (1e4 for CP10K, 1e6 for CPM). Raw counts are
    reconstructed as `raw_i = expm1(X_i) · nCount_i / T` (verified integer-exact on
    chick_full / chick_rpc), so Panel A is correct for any `log1p(CP*)` normalisation
    and stays in the same units as Panel B. A near-integrality check additionally
    guards the depth basis; Panel A is **omitted** (per-cell violin only) when `.X` is
    not a clean `log1p(CP)` (e.g. z-scored, scran-pooled) — no raw-count layer needed.
- **Panel B — per-cell positive-cell violins**: `y =` the per-cell `.X` log-norm value
  (`log1p(CP·T)`, same units as Panel A) for cells with `X > 0`; one violin per group
  (`scalemode="width"`); faint jitter overlaid except for groups with > 2000 positive
  cells; **`n=` = number of positive cells** annotated above each violin. (With the
  positive-cells-only gate off, Panel B shows all cells and `n=` is the per-group cell
  count.)
- x-axis = `annotation_column` in `annotation_order`; colours from
  `annotation_colors`.

If the dataset lacks `nCount_RNA` or a replicate column, Panel A is omitted and only
the per-cell violin panel is drawn.

## Dot plot

Gene set (rows) × groups (columns). Dot **size** = fraction of cells with expr > 0;
dot **colour** = mean expression over all cells in the group (the standard dot-plot
convention). The
x-axis follows `annotation_order`. A default module is preselected so the dot plot
renders immediately when chosen.

## Share state (v2)

Base64-encoded JSON carried in `?state=`. Keys (only those relevant to the current
view are written; absent keys restore to defaults):

```
v=2, dataset, embedding, color, gene, mode, plot_type,
compare_genes, gene2, compare_shared_scale,                       # two-gene
custom_x, custom_y, bins, percentile, enable_binning,             # map binning
  enable_smoothing, bin_stat,
custom_projection, wm_rho_nt, wm_rho_dv, wm_gap, wm_stretch,      # whole-mount flower / sphere (3D)
  wm_cuts, wm_symmetric, wm_dewarp, wm_pow, wm_gap_mode,
  wm_gap_frac, wm_pole,
group_by, group_split, group_gene, group_style, gene_module,      # group view
  group_positive_only, group_replicate, group_value_source, group_meta
```

**Restoration ownership** (each control is written by exactly one restore-capable
callback, and each reads the same URL state):

- `callbacks/url_callbacks.initialize_from_url` → dataset + global controls
  (mode, bins, percentile, enable_binning, enable_smoothing, bin_stat, plot_type,
  compare_genes, gene2, compare_shared_scale). Fires on initial load
  (`prevent_initial_call='initial_duplicate'`) so a pasted link works.
- `app.update_data` → embedding, colour, gene (triggered by the restored dataset).
- `main_callbacks.update_custom_embedding_controls` → custom_x / custom_y.
- `main_callbacks.populate_group_controls` → all group controls.

Every URL-reading owner (`update_data`, `update_custom_embedding_controls`,
`populate_group_controls`) applies the shared state **only when the URL's `dataset`
matches the loaded `dataset_id`**, so a later manual dataset switch is not
re-overridden by the stale `?state=` (which is never cleared from the URL).

The earlier duplicate `initialize_from_url` in `app.py` was removed so no two
callbacks write the same control value.
