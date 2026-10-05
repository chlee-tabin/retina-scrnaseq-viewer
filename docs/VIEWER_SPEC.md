# Retina scRNA-seq Viewer — Specification

A living spec for the Dash/Plotly single-cell viewer. It defines the view model,
the controls, the dataset-config schema, the shareable-state schema (v4), and the
definition of the "Expression by group" figure. **This document is the contract used
for independent review** — each claim below should hold in the code.

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
  colour mode, two-gene side-by-side, flower and 3D sphere projections, and an
  optional chick HAA pointer / uncaptured-nasal overlay. The raw custom map supports
  exploratory polygon ROI pseudobulk differential expression and figure-region presets.
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
| `enable-smoothing` | map (custom) | Gaussian smoothing of the binned map (on/off gate) |
| `smoothing-mode` | map (custom) | Manuscript zero-fill (default) vs opt-in mask-normalised smoothing |
| `smooth-sigma-slider` | map (custom) | Smoothing strength σ (0–4, bin units); per-dataset config = default |
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
| `custom-projection` | map (custom) | Raw axes vs whole-mount flower vs 3D sphere |
| `wm-rho-nt` / `wm-rho-dv` / `wm-gap` / `wm-stretch` / `wm-cuts` | map (flower) | Flower extent, relief width, petal stretch, number of cuts |
| `wm-symmetric` / `wm-pole` / `wm-dewarp` / `wm-pow` / `wm-gap-mode` / `wm-gap-frac` | map (flower) | Advanced flower geometry (defaults come from the current dataset) |
| `haa-mode` | map (custom, chick) | Off / footprint / expression / domain / peak HAA pointer |
| `nasal-gap-layers` / `nasal-gap-frac` / `nasal-gap-dorsal` / `nasal-gap-ventral` | flower/sphere (chick) | Render-only uncaptured-nasal cap reach / depth / arc |
| `sphere-camera-store` | sphere | Persist and share the 3D camera |
| `wm-reset` | flower/sphere | Reset to the current dataset's projection defaults |
| ROI A/B draw / clear, `roi-preset-select` | map (raw custom axes) | Polygon ROIs or manuscript bin-region presets |
| `deg-min-cells` / `run-deg-btn` | ROI DE | Pseudobulk cell floor and explicit run |
| `volcano-lfc-thresh` / `volcano-padj-thresh` / `deg-download-btn` | ROI DE | Display thresholds and result CSV |

Style-specific group controls are revealed by `group-style`: the figure controls
show only for `figure`, the module selector only for `dotplot`.

## Dataset config schema (`datasets_config.yml`)

Required per-dataset keys: `title`, `description`, `file_path`, `last_updated`
(the same contract as README). Other keys are optional, including `metadata`,
`default_gene`, `default_embedding`, `smooth_sigma`, `min_cells_per_bin`,
`color_floor`, `gene_modules`, plus:

- **`sha256`** — streamed checksum; configured downloads and local files must match.
  Every shipped dataset is pinned. Verification is cached by path, size and mtime.
- **`wholemount_params`** — overrides for `utils/wholemount.py::LOCKED_PARAMS`.
  Chick keeps the tuned directional stretch; human/mouse use an empty stretch list.
  Dataset changes and Reset projection apply this dataset's defaults.
- **`figure_regions`** — named `{name, n, nt: [lo, hi], dv: [lo, hi]}` presets.
  Select directly on zero-based bin indices with half-open tests
  `nt[0] <= NT_bin < nt[1]+1` and `dv[0] <= DV_bin < dv[1]+1`.
  Axis maxima have index `n` and are excluded. Rectangles display the outline;
  preset selection uses the bin test, including after a rounded share restore.
- **`haa_marker`** — marker for the optional HAA pointer (chick `CYP26C1`).
- **`nasal_gap`** — render-only uncaptured-nasal cap (`frac`, `layers`, optional
  `dorsal_deg`, `ventral_deg`); absent means no overlay.

Top-level **`hf_data_revision`** pins the Hugging Face data commit; the
`HF_DATA_REVISION` environment variable overrides it. A leading `data/` in
`file_path` is removed before joining under `DATA_DIR`; nested paths are preserved
across loading, provisioning and downloads. Dataset basenames must be unique.

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

Plotly reproduction of the manuscript Methods (expression-by-group figure). For a gene `g`, a grouping column, and a replicate unit:

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

## Spatial smoothing

Zero-fill Gaussian smoothing remains the default to match the manuscript pipeline:
masked bins contribute zero to the convolution and are masked again afterward.
The opt-in `mask_normalised` mode uses `smooth(values × mask) / smooth(mask)` and
re-masks invalid bins; it preserves an isolated valid 100% bin at 100% rather than
attenuating it near empty bins. Sigma and the enable-smoothing gate apply to both
modes, in raw, flower, sphere and comparison views. The HAA marker uses the same mode.

## ROI differential expression (exploratory)

Polygon selection uses raw custom-axis coordinates; manuscript figure presets use
the half-open bin gates above. With one ROI, rest excludes all cells with non-finite
plot coordinates. Two ROIs form a disjoint A/B contrast, with overlap going to the
smaller ROI. The volcano names the actual A/B or region/rest contrast and counts
excluded `padj = NA` genes.

Replicate unit = configured `replicate_columns` (chick: library × genotype;
human/mouse: library); cells missing any replicate label are excluded and reported.
Counts must be verified raw counts or an integer-like log1p(CP) reconstruction.
The design is `~ library + condition` when estimable, otherwise explicitly reported
`~ condition`; fit failures never trigger a silent refit. The `min_frac=0.5` gene
filter is a viewer-only addition that changes padj versus the manuscript's tables.
The [public analysis repository](https://github.com/chlee-tabin/retina-spatial-scrna-analysis)
contains `scripts/figures/sfig24_area_deg.R` and
`scripts/figures/supptables1_2_area_deg.R` (chick Fig 6H, human Fig S24).

## Share state (v4)

Base64-encoded JSON carried in `?state=`. The single owner is
`utils/state.py::STATE_SCHEMA` (component, property, default, coercion, version
introduced, restoration owner). New links capture all controls; omitted keys in
older links reset to schema defaults, clearing a recipient's comparison gene and
shared scale. Dataset-dependent defaults are resolved from server configuration.

```
v=4, dataset, data_version, embedding, color, gene, mode, plot_type,
  smooth_sigma, smoothing_mode,
compare_genes, gene2, compare_shared_scale,                       # two-gene
custom_x, custom_y, bins, percentile, enable_binning,             # map binning
  enable_smoothing, bin_stat,
custom_projection, wm_rho_nt, wm_rho_dv, wm_gap, wm_stretch,      # whole-mount flower / sphere (3D)
  wm_cuts, wm_symmetric, wm_dewarp, wm_pow, wm_gap_mode,
  wm_gap_frac, wm_pole,
ng_layers, ng_frac, ng_dorsal, ng_ventral,                        # uncaptured-nasal cap (flower + sphere)
haa_mode,                                                         # HAA pointer mode (off/footprint/expression/domain/peak)
sphere_cam,                                                       # 3D-sphere viewpoint (camera; both panels)
group_by, group_split, group_gene, group_style, gene_module,      # group view
  group_positive_only, group_replicate, group_value_source, group_meta,
roi, volcano_lfc, volcano_padj, deg_min_cells,                    # ROI DE
```

Added in v0.221 (additive, backward compatible — pre-v0.221 links omit these keys
and restore the R2.5 cap defaults + the default face-on camera): `ng_layers`,
`ng_frac`, `ng_dorsal`, `ng_ventral` (the uncaptured-nasal-cap extent / wedge,
historically written for non-raw whole-mount views); `haa_mode` (the HAA-pointer
definition); and `sphere_cam` (one camera covers both compare panels). v4 captures
these controls in every link, including a null camera for the face-on default.

Added in v0.31 — share schema bumped **v2 → v3** (additive, backward compatible):
`smooth_sigma`, the Gaussian smoothing strength, now a live slider rather than a fixed
per-dataset config value (the config value is just the slider's default on dataset
load). Captured for every new link (top-level, not only the spatial view) so a recipient
who later switches to the map gets the sharer's strength. A **pre-slider link is
identified by its schema version** (`v < 3`): for those, `decode_state` migrates the
per-dataset default that was in effect then (`LEGACY_SMOOTH_SIGMA` in
`utils/smoothing.py` — `human_rpc` **0.5** (the old near-no-op default, since corrected
to 2.0), `chick_rpc`/`mouse_rpc`/`mouse_rpc_legacy` 2.0, `chick_full` 1.5), so the
shared view is reproduced exactly. Identifying old links by version (not by a missing
key) stays correct even if a future change makes the `smooth_sigma` write conditional.
An out-of-range σ from a hand-edited URL is clamped to the slider's [0, 4] domain.

Added in v0.32 — schema **v4** adds `data_version` and `smoothing_mode`.
Dataset ids are stable: replacing data requires a new id or the version notice.
When retaining an id, use a new file basename so older links trigger the notice;
do not replace data in place under the same basename.
`data_version` is the configured file basename, resolved server-side. A mismatch
still restores but displays: "This link was created with an earlier version of this
dataset (<old>); it now shows <new>." Links without the field restore silently.
Missing `smoothing_mode` uses historical zero-fill; pre-v3 links retain the legacy
sigma migration above. Missing nasal-cap keys retain their historical defaults,
and a missing camera uses the face-on default. These are defaults, never the
recipient's previously selected settings.

**Restoration ownership** (option-populating callbacks retain their owners; all
value Outputs and restore values come from the schema):

- `callbacks/url_callbacks.initialize_from_url` → dataset + global controls
  (mode, bins, percentile, enable_binning, enable_smoothing, bin_stat, plot_type,
  compare_genes, gene2, compare_shared_scale, smoothing_mode, ROI/volcano controls, the uncaptured-nasal-cap controls
  ng_layers/ng_frac/ng_dorsal/ng_ventral, and the sphere camera `sphere_cam`).
  Fires on initial load (`prevent_initial_call='initial_duplicate'`) so a pasted
  link works.
- `app.update_data` → embedding, colour, gene, and `smooth-sigma-slider` (triggered
  by the restored dataset; the slider's sigma is resolved URL-value > legacy default
  > config default by `resolve_smooth_sigma`).
- `main_callbacks.update_custom_embedding_controls` → custom_x / custom_y.
- `main_callbacks.populate_group_controls` → group options and values.
- `main_callbacks.reset_wholemount_params` → dataset projection defaults (or the
  matching URL values on dataset load); Reset always restores dataset defaults.

Every URL-reading owner (`update_data`, `update_custom_embedding_controls`,
`populate_group_controls`) applies the shared state **only when the URL's `dataset`
matches the loaded `dataset_id`**, so a later manual dataset switch is not
re-overridden by the stale `?state=` (which is never cleared from the URL).

The earlier duplicate `initialize_from_url` in `app.py` was removed so no two
callbacks write the same control value.
