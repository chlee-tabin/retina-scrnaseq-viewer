from dash import Input, Output, State, callback
from utils.data_loading import load_adata, load_dataset_config, gene_search_options
from utils.plotting import (
    create_scatter_plot,
    create_binned_plot,
    create_group_expression_plot,
    create_group_expression_figure,
    create_dotplot,
    create_dual_gene_figure,
    create_wholemount_binned_figure,
    create_sphere_figure,
    _message_figure,
)
import pandas as pd
import logging
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import state_for_dataset
from utils import wholemount
import scipy.sparse
import numpy as np
import traceback

logger = logging.getLogger(__name__)


def _gene_vector(adata, gene):
    """Return a dense 1D expression vector for `gene` (handles sparse X)."""
    sub = adata[:, gene].X
    if scipy.sparse.issparse(sub):
        return sub.toarray().flatten()
    return np.asarray(sub).flatten()


def _on(value):
    """True iff a dcc.Checklist value contains the 'enabled' flag. Tolerates a
    non-list (e.g. a hand-edited ?state=) without raising on the membership test."""
    return isinstance(value, (list, tuple)) and 'enabled' in value


# Internal separator for composite replicate keys: a rare control char that cannot
# appear in obs string values, so a value that itself contains '|' never mis-splits.
_REP_SEP = '\x1f'


def _replicate_key(cols):
    """Build the dropdown value for a multi-column replicate unit (collision-safe)."""
    return _REP_SEP.join(cols)


def _replicate_series(adata, replicate_value, replicate_columns):
    """Resolve the per-cell pseudobulk replicate label for the group figure.

    `replicate_value` is the dropdown choice: a single obs column, a composite key
    built by `_replicate_key` (several columns), or None -> fall back to the dataset's
    configured `replicate_columns` composite, then a single configured column. A
    requested composite must resolve to ALL of its columns: if any is missing from obs
    the replicate is treated as unavailable (return None -> Panel A omitted) rather
    than silently narrowing to a coarser unit. Returns a numpy array of string labels,
    or None if no usable replicate exists.
    """
    obs = adata.obs

    def join(cols):
        return obs[cols].astype(str).agg(_REP_SEP.join, axis=1).to_numpy()

    if replicate_value and replicate_value in obs.columns:
        return obs[replicate_value].astype(str).to_numpy()
    if replicate_value and _REP_SEP in replicate_value:
        want = replicate_value.split(_REP_SEP)
        present = [c for c in want if c in obs.columns]
        return join(present) if present and len(present) == len(want) else None

    configured = list(replicate_columns or [])
    present = [c for c in configured if c in obs.columns]
    if len(configured) >= 2:
        return join(present) if len(present) == len(configured) else None
    if len(configured) == 1:
        return obs[present[0]].astype(str).to_numpy() if present else None
    return None


def _resolve_group_field(adata, field):
    """(display_name, per-cell string array) for a group-by / split-by selection.

    A composite key built by `_replicate_key` (e.g. library x genotype) resolves to its
    columns joined for display ('libA | donor0'), so a per-library demux donor label
    ('genotype') is only ever grouped WITHIN its library -- never on its own, where the
    repeated donor labels across libraries would be meaningless. A plain obs column
    resolves to itself. Returns (field, None) when the field is missing/unresolvable.
    """
    if not field:
        return field, None
    if _REP_SEP in field:
        cols = field.split(_REP_SEP)
        if any(c not in adata.obs.columns for c in cols):
            return field, None
        disp = ' x '.join(cols)
        ser = adata.obs[cols].astype(str).agg(' | '.join, axis=1).to_numpy()
        return disp, ser
    if field in adata.obs.columns:
        return field, adata.obs[field].astype(str).to_numpy()
    return field, None


def _annotation_style(data_store, column):
    """(color_map, category_order) for a categorical `column`: the configured per-type
    colours + order when `column` is the dataset's annotation column, else (None,
    None) so the default palette/order applies."""
    if column and column == data_store.get('annotation_column'):
        return (data_store.get('annotation_colors') or {},
                data_store.get('annotation_order') or [])
    return None, None


def _is_categorical_series(series):
    """Heuristic used across the viewer: category/object dtype or small-int."""
    return (
        series.dtype.name in ['category', 'object'] or
        (pd.api.types.is_integer_dtype(series) and len(series.unique()) <= 50)
    )


def _wholemount_params(rho_nt, rho_dv, gap_gain, stretch_scale, n_cuts,
                       symmetric=None, dewarp=None, pow_p=None,
                       gap_mode=None, gap_frac=None, pole=None):
    """Merge the Advanced-projection control values over the shipped reviewer preset
    (wholemount.LOCKED_PARAMS) into a params dict for flower_transform. None values keep
    the preset, so the defaults reproduce Fig R2.5 exactly."""
    P = dict(wholemount.LOCKED_PARAMS)
    if rho_nt is not None:
        P['rho_max_nt_deg'] = float(rho_nt)
    if rho_dv is not None:
        P['rho_max_dv_deg'] = float(rho_dv)
    if gap_gain is not None:
        P['gap_gain'] = float(gap_gain)
    # Scale the locked directional stretch bumps (temporal / ventro-temporal / dorso-nasal)
    # by the slider; 0 -> no anatomical stretch, 1 -> the shipped amounts.
    if stretch_scale is not None:
        s = float(stretch_scale)
        P['stretch_bumps'] = tuple((ang, amp * s, wid)
                                   for (ang, amp, wid) in wholemount.LOCKED_PARAMS['stretch_bumps'])
    # Evenly spaced relief cuts; 0 -> a solid disk (no slits).
    if n_cuts is not None:
        k = int(n_cuts)
        P['cut_angles_deg'] = tuple(i * 360.0 / k for i in range(k)) if k > 0 else ()
    # Knobs not pinned by LOCKED_PARAMS (they inherit DEFAULT_PARAMS): exposed so the user
    # can explore them. Each stays at the R2.5 value unless its control overrides it.
    if symmetric is not None:
        P['symmetric'] = bool(symmetric)   # False -> p99 per-axis scaling (true asymmetry)
    if dewarp is not None:
        P['dewarp'] = dewarp               # 'arcsin' (R2.5) | 'pow' | 'none'
    if pow_p is not None:
        P['pow_p'] = float(pow_p)          # exponent for dewarp='pow'
    if gap_mode is not None:
        P['gap_mode'] = gap_mode           # 'deficit' (R2.5) | 'linear'
    if gap_frac is not None:
        P['gap_frac'] = float(gap_frac)    # rip width for gap_mode='linear'
    if pole is not None:
        P['pole'] = pole                   # 'origin' (R2.5) | 'median'
    return P


@callback(
    Output('main-plot', 'figure', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('custom-x-select', 'value'),
     Input('custom-y-select', 'value'),
     Input('color-select', 'value'),
     Input('gene-select', 'value'),
     Input('viz-mode', 'value'),
     Input('bin-number-slider', 'value'),
     Input('percentile-slider', 'value'),
     Input('enable-binning', 'value'),
     Input('enable-smoothing', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search'),
     # New inputs (appended): plot-type switch, expression-by-group controls,
     # and the second-gene comparison controls.
     Input('plot-type', 'value'),
     Input('group-gene-select', 'value'),
     Input('group-by-select', 'value'),
     Input('group-split-select', 'value'),
     Input('group-style', 'value'),
     Input('gene-module-select', 'value'),
     Input('compare-genes', 'value'),
     Input('gene-select-2', 'value'),
     Input('bin-stat', 'value'),
     Input('group-positive-only', 'value'),
     Input('group-replicate-select', 'value'),
     Input('compare-shared-scale', 'value'),
     Input('custom-projection', 'value'),
     Input('wm-rho-nt', 'value'),
     Input('wm-rho-dv', 'value'),
     Input('wm-gap', 'value'),
     Input('wm-stretch', 'value'),
     Input('wm-cuts', 'value'),
     Input('group-value-source', 'value'),
     Input('group-meta-select', 'value'),
     Input('wm-symmetric', 'value'),
     Input('wm-dewarp', 'value'),
     Input('wm-pow', 'value'),
     Input('wm-gap-mode', 'value'),
     Input('wm-gap-frac', 'value'),
     Input('wm-pole', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, custom_x, custom_y, color_by, gene, viz_mode,
                bin_number, percentile, enable_binning, enable_smoothing, selection_data,
                url_search, plot_type, group_gene, group_by, group_split, group_style,
                gene_module, compare_genes, gene2, bin_stat,
                group_positive_only, group_replicate, compare_shared_scale,
                custom_projection, wm_rho_nt, wm_rho_dv, wm_gap, wm_stretch, wm_cuts,
                group_value_source, group_meta,
                wm_symmetric, wm_dewarp, wm_pow, wm_gap_mode, wm_gap_frac, wm_pole):
    logger.debug("update_plot called with parameters:")

    # Initialize treat_as_categorical as False by default
    treat_as_categorical = False

    if not data_store:
        return {}

    try:
        config = load_dataset_config()
        adata = load_adata(data_store['filename'])

        # ---- Expression-by-group view (violin/box/strip/dotplot) ----
        if plot_type == 'group':
            if not group_by:
                return _message_figure("Select a categorical column to group by.")

            value_is_gene = (group_value_source != 'meta')

            # Resolve the grouping field. A composite (library x genotype) demux key maps
            # to a per-cell 'libA | donor0' label, so a per-library demux donor ('genotype')
            # is only ever grouped WITHIN its library, never on its own.
            g_disp, g_series = _resolve_group_field(adata, group_by)
            if g_series is None:
                return _message_figure("The selected grouping column is not available in this dataset.")
            # Per-type colour/order only applies to the dataset's plain annotation column.
            g_cmap, g_corder = _annotation_style(data_store, group_by)

            # Dot plot is a gene-set view; it does not apply to a continuous metadata value.
            if group_style == 'dotplot':
                if not value_is_gene:
                    return _message_figure(
                        "Dot plot summarizes gene sets. Switch 'Plot value' to Gene "
                        "expression, or choose violin / box / strip / figure."
                    )
                gene_modules = data_store.get('gene_modules', {}) or {}
                module_genes = gene_modules.get(gene_module) if gene_module else None
                if module_genes:
                    genes = [g for g in module_genes if g in adata.var_names]
                elif group_gene:
                    genes = [group_gene] if group_gene in adata.var_names else []
                else:
                    return _message_figure("Select a gene or a gene module for the dot plot.")
                if not genes:
                    return _message_figure(
                        "None of the selected genes are present in this dataset."
                    )
                expr_data = {g: _gene_vector(adata, g) for g in genes}
                expr_df = pd.DataFrame(expr_data)
                expr_df[g_disp] = g_series
                return create_dotplot(expr_df, genes, g_disp, category_order=g_corder)

            # Resolve the per-cell value: a gene's log-norm expression, or a continuous
            # .obs variable (QC metric, DV/NT score, ...).
            if value_is_gene:
                if not group_gene:
                    return _message_figure("Select a gene to plot expression by group.")
                if group_gene not in adata.var_names:
                    return _message_figure(
                        f"Gene '{group_gene}' not found in this dataset. "
                        "Try a different gene or check the name/casing."
                    )
                value_name = group_gene
                value_vec = _gene_vector(adata, group_gene)
            else:
                if not group_meta or group_meta not in adata.obs.columns:
                    return _message_figure("Select a continuous metadata variable to plot.")
                value_name = group_meta
                value_vec = pd.to_numeric(adata.obs[group_meta], errors='coerce').to_numpy()

            # Polished 2-panel figure. For a gene, Panel A is the depth-normalised log1p(CP)
            # pseudobulk and Panel B the positive-cell violin; for a continuous metadata
            # variable, Panel A is the per-replicate MEAN of the variable and Panel B the
            # all-cell violin (create_group_expression_figure branches on value_is_gene).
            if group_style == 'figure':
                fig_df = pd.DataFrame({'expr': value_vec})
                fig_df[g_disp] = g_series
                if 'nCount_RNA' in adata.obs.columns:
                    fig_df['ncount'] = adata.obs['nCount_RNA'].to_numpy()
                rep_series = _replicate_series(
                    adata, group_replicate, data_store.get('replicate_columns', []))
                if rep_series is not None:
                    fig_df['replicate'] = rep_series
                positive_only = _on(group_positive_only)
                return create_group_expression_figure(
                    fig_df, value_name, g_disp,
                    color_map=g_cmap, category_order=g_corder,
                    positive_only=positive_only,
                    norm_target=data_store.get('x_norm_target'),
                    value_is_gene=value_is_gene,
                )

            # Violin / box / strip (and the metadata 'figure' fallback).
            style = group_style if group_style in ('violin', 'box', 'strip') else 'violin'
            df = pd.DataFrame({'expr': value_vec})
            df[g_disp] = g_series
            s_disp, s_series = _resolve_group_field(adata, group_split)
            if s_series is not None and s_disp != g_disp:
                df[s_disp] = s_series
                split_col = s_disp
            else:
                split_col = None
            return create_group_expression_plot(
                df, value_name, g_disp, split_by=split_col, style=style,
                color_map=g_cmap, category_order=g_corder, value_is_gene=value_is_gene,
            )

        # ---- Embedding / spatial map view ----
        if not embedding:
            return {}

        # Get coordinates based on embedding type. `embedding_label` is the string the
        # plotting helpers use for axis labels + styling; it diverges from `embedding`
        # only for the whole-mount projection (a 'wholemount' styling token), while
        # `embedding` stays 'custom_embedding' so the binning logic below is unchanged.
        embedding_label = embedding
        proj_params = None   # whole-mount transform params (None unless the flower view)
        sphere_xyz = None    # (x, y, z) on the unit sphere when the 3D projection is active
        if embedding == 'custom_embedding':
            if custom_projection in ('flower', 'sphere'):
                # Whole-mount reprojection of the DV/NT topographic scores: the flat "flower"
                # or its native 3D sphere. Both need the score columns; datasets without them
                # (e.g. the full-retina object) get a clear message.
                if not wholemount.has_scores(adata.obs.columns):
                    return _message_figure(
                        "Whole-mount / sphere projection needs DV.Score and NT.Score, which "
                        "this dataset doesn't have. It's available for the RPC datasets "
                        "(chick / human / mouse retinal progenitor cells)."
                    )
                # Advanced-projection controls merged over the reviewer preset (defaults = R2.5).
                proj_params = _wholemount_params(
                    wm_rho_nt, wm_rho_dv, wm_gap, wm_stretch, wm_cuts,
                    symmetric=_on(wm_symmetric), dewarp=wm_dewarp, pow_p=wm_pow,
                    gap_mode=wm_gap_mode, gap_frac=wm_gap_frac,
                    pole=('median' if wm_pole == 'median' else 'origin'),
                )
                # Fit the pole + per-axis scale ONCE from the cells and reuse it for every
                # warp (scatter cells AND the binned grid), so the views share one basis
                # and the layout stays stable under cell subsetting.
                dv_cells = adata.obs[wholemount.DV_COL].to_numpy()
                nt_cells = adata.obs[wholemount.NT_COL].to_numpy()
                proj_params = {**proj_params,
                               **wholemount.compute_scale_fit(dv_cells, nt_cells, **proj_params)}
                if custom_projection == 'sphere':
                    # Native 3D geometry: reuse the flower's basis but place cells on the
                    # unit sphere. Built below (a 3D figure) once the colour series resolves;
                    # x/y are placeholders for the shared DataFrame / compare guards.
                    sphere_xyz = wholemount.sphere_coords(dv_cells, nt_cells, params=proj_params)
                    x, y = sphere_xyz[0], sphere_xyz[1]
                else:
                    # Per-cell warp (used by the scatter + dual-gene views; the binned view
                    # re-bins in score space below for a stray-free, faithful map).
                    x, y = wholemount.wholemount_coords(dv_cells, nt_cells, params=proj_params)
                x_label = y_label = ''
                embedding_label = 'wholemount'
            else:
                if not custom_x or not custom_y:
                    return {}
                x = adata.obs[custom_x]
                y = adata.obs[custom_y]
                x_label = custom_x
                y_label = custom_y
        else:
            coordinates = adata.obsm[embedding]
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            x_label = f"{embedding}_1"
            y_label = f"{embedding}_2"

        # Binning is a SPATIAL-map concept: only bin the custom DV/NT embedding
        # (when its binning toggle is on). Every other embedding (UMAP/PCA/...)
        # always renders as a per-cell scatter -- this avoids the pixelated
        # heatmap and the useless per-category subplot grid on UMAPs.
        binning_on = (
            embedding == 'custom_embedding'
            and _on(enable_binning)
        )

        # ---- Two genes side-by-side (map view, gene-expression color) ----
        # The 3D sphere renders a single variable; skip the (2D) dual-gene comparison there.
        compare_on = _on(compare_genes)
        if color_by == 'gene_expression' and compare_on and gene and gene2 and sphere_xyz is None:
            missing = [g for g in (gene, gene2) if g not in adata.var_names]
            if missing:
                return _message_figure(
                    f"Gene(s) not found in this dataset: {', '.join(missing)}. "
                    "Check the name/casing."
                )
            smooth_on = _on(enable_smoothing)
            smooth_sigma = float(data_store.get('smooth_sigma', 0) or 0) if smooth_on else 0
            return create_dual_gene_figure(
                np.asarray(x), np.asarray(y),
                [_gene_vector(adata, gene), _gene_vector(adata, gene2)],
                [gene, gene2],
                embedding_label,
                binned=binning_on,
                bin_size=bin_number,
                percentile=percentile,
                smooth_sigma=smooth_sigma,
                min_cells=data_store.get('min_cells_per_bin', 1),
                color_floor=data_store.get('color_floor', 0.05),
                bin_stat=(bin_stat or 'mean'),
                shared_scale=_on(compare_shared_scale),
            )

        # Handle gene expression
        if color_by == 'gene_expression':
            logger.debug("Color by gene_expression selected")
            if not gene:
                logger.debug("No gene selected")
                color_series = None
            elif gene not in adata.var_names:
                # Guard against genes absent from this dataset (e.g. a shared
                # URL, a stale selection, or wrong gene-name casing across
                # species). Return a friendly message instead of a 500.
                logger.info(f"Gene '{gene}' not found in dataset '{data_store['filename']}'")
                return _message_figure(
                    f"Gene '{gene}' not found in this dataset. "
                    "Try a different gene or check the name/casing."
                )
            else:
                color_series = _gene_vector(adata, gene)
                treat_as_categorical = False
        else:
            logger.debug(f"Color by '{color_by}' selected")
            color_series = adata.obs[color_by] if color_by else None
            if color_series is not None:
                logger.debug(f"Color series dtype: {color_series.dtype}")
                treat_as_categorical = _is_categorical_series(color_series)
            logger.debug(f"treat_as_categorical: {treat_as_categorical}")

        # ---- 3D spherical whole-mount: the topographic map on its native geometry ----
        # Renders a complete 3D figure (bypasses the 2D binning/scatter dispatch below).
        # Colour resolution above is shared; categorical annotation keeps its configured map.
        if sphere_xyz is not None:
            s_cmap, s_corder = _annotation_style(data_store, color_by)
            return create_sphere_figure(
                sphere_xyz[0], sphere_xyz[1], sphere_xyz[2],
                color_series, color_by, gene=gene,
                treat_as_categorical=treat_as_categorical,
                color_map=s_cmap, category_order=s_corder,
            )

        # Create DataFrame for plotting
        df = pd.DataFrame({
            'x': x,
            'y': y,
            'color': color_series
        })

        # Create the plot. Only the custom spatial embedding is ever binned;
        # all other embeddings always use the per-cell scatter (F1).
        if binning_on:
            smooth_on = _on(enable_smoothing)
            smooth_sigma = float(data_store.get('smooth_sigma', 0) or 0) if smooth_on else 0
            min_cells = data_store.get('min_cells_per_bin', 1)
            if custom_projection == 'flower' and not treat_as_categorical and color_series is not None:
                # Faithful whole-mount: bin in DV/NT SCORE space, warp, drop strays --
                # instead of histogramming the warped per-cell coordinates.
                base = gene if color_by == 'gene_expression' else str(color_by)
                cl = f"{base} (% detected)" if bin_stat == 'frac_pos' else base
                fig = create_wholemount_binned_figure(
                    dv_cells,
                    nt_cells,
                    np.asarray(color_series, dtype=float),
                    bin_size=bin_number,
                    percentile=percentile,
                    smooth_sigma=smooth_sigma,
                    min_cells=min_cells,
                    color_floor=data_store.get('color_floor', 0.05),
                    bin_stat=(bin_stat or 'mean'),
                    color_label=cl,
                    params=proj_params,
                )
            else:
                fig = create_binned_plot(
                    df, embedding_label, color_by,
                    bin_size=bin_number,
                    percentile=percentile,
                    treat_as_categorical=treat_as_categorical,
                    smooth_sigma=smooth_sigma,
                    min_cells=min_cells,
                    color_floor=data_store.get('color_floor', 0.05),
                    bin_stat=(bin_stat or 'mean'),
                )
        else:
            s_cmap, s_corder = _annotation_style(data_store, color_by)
            fig = create_scatter_plot(
                df, embedding_label, color_by,
                treat_as_categorical=treat_as_categorical,
                color_map=s_cmap, category_order=s_corder,
                plot_order=(viz_mode or 'random'),
            )

        return fig

    except Exception as e:
        logger.error(f"Error updating plot: {str(e)}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        return _message_figure("Internal error rendering this view -- see server logs.")

@callback(
    [Output('custom-x-select', 'options'),
     Output('custom-y-select', 'options'),
     Output('custom-x-select', 'value'),
     Output('custom-y-select', 'value'),
     Output('custom-embedding-container', 'style')],
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_custom_embedding_controls(data_store, embedding, url_search):
    if not data_store or not data_store.get('metadata_cols'):
        return [], [], None, None, {'display': 'none'}
    
    # Get numeric columns for custom embedding
    numeric_cols = [
        col for col, type_ in data_store.get('column_types', {}).items()
        if type_ == 'numeric'
    ]
    
    options = [{'label': col, 'value': col} for col in numeric_cols]
    
    # Set container visibility based on embedding type
    container_style = {'display': 'block'} if embedding == 'custom_embedding' else {'display': 'none'}
    
    # Check if we should initialize from URL state
    if url_search:
        # state_for_dataset returns the shared state only when it belongs to the loaded
        # dataset (a stale ?state= from another dataset must not re-apply on a switch).
        state = state_for_dataset(url_search, data_store.get('dataset_id'))
        if state and state.get('embedding') == 'custom_embedding':
            custom_x = state.get('custom_x')
            custom_y = state.get('custom_y')
            # Only set values if they exist in the options
            valid_values = [opt['value'] for opt in options]
            if custom_x in valid_values and custom_y in valid_values:
                return options, options, custom_x, custom_y, container_style
    
    # Default return if not initializing from URL: pre-select the topographic
    # axes (NT.Score on x, DV.Score on y) when those score columns exist.
    default_x = 'NT.Score' if 'NT.Score' in numeric_cols else None
    default_y = 'DV.Score' if 'DV.Score' in numeric_cols else None
    return options, options, default_x, default_y, container_style

@callback(
    Output('viz-mode', 'style'),
    [Input('embedding-select', 'value'),
     Input('enable-binning', 'value')]
)
def toggle_plot_order(embedding, enable_binning):
    if embedding == 'custom_embedding' and _on(enable_binning):
        return {'display': 'none'}
    return {'display': 'block'}

@callback(
    Output('binning-controls', 'style'),
    [Input('enable-binning', 'value')]
)
def toggle_binning_controls(enable_binning):
    if _on(enable_binning):
        return {'display': 'block'}
    return {'display': 'none'}


@callback(
    Output('raw-axes-controls', 'style'),
    Input('custom-projection', 'value')
)
def toggle_raw_axes(custom_projection):
    # The whole-mount projections (flat flower and 3D sphere) source their coordinates from
    # DV.Score / NT.Score, so the manual X/Y axis pickers are irrelevant there -- hide them.
    if custom_projection in ('flower', 'sphere'):
        return {'display': 'none'}
    return {'display': 'block'}


@callback(
    Output('wholemount-advanced', 'style'),
    Input('custom-projection', 'value')
)
def toggle_wholemount_advanced(custom_projection):
    # The Advanced-projection controls apply to both whole-mount views. On the sphere the
    # extent / de-warp / pole / symmetric knobs reshape the cap; the relief gap / cuts /
    # stretch are flat-layout-only and have no effect there (the sphere has no rips).
    if custom_projection in ('flower', 'sphere'):
        return {'display': 'block'}
    return {'display': 'none'}


# ---- F2: toggle map vs. expression-by-group control panels ----
@callback(
    [Output('map-controls', 'style'),
     Output('group-controls', 'style')],
    Input('plot-type', 'value')
)
def toggle_plot_type(plot_type):
    if plot_type == 'group':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'block'}, {'display': 'none'}


# ---- Within the group view, reveal only the controls the chosen style uses ----
@callback(
    [Output('figure-style-controls', 'style'),
     Output('module-controls', 'style')],
    Input('group-style', 'value')
)
def toggle_group_style_controls(group_style):
    fig_style = {'display': 'block'} if group_style == 'figure' else {'display': 'none'}
    mod_style = {'display': 'block'} if group_style == 'dotplot' else {'display': 'none'}
    return fig_style, mod_style


# ---- Group view: show the gene picker or the continuous-metadata picker ----
@callback(
    [Output('group-gene-block', 'style'),
     Output('group-meta-block', 'style')],
    Input('group-value-source', 'value')
)
def toggle_group_value_source(group_value_source):
    if group_value_source == 'meta':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'block'}, {'display': 'none'}


# ---- Group view: only the distribution styles apply to a continuous metadata value ----
@callback(
    [Output('group-style', 'options'),
     Output('group-style', 'value', allow_duplicate=True)],
    Input('group-value-source', 'value'),
    State('group-style', 'value'),
    prevent_initial_call=True
)
def restrict_group_styles(group_value_source, current_style):
    """The dot plot summarizes a GENE set, so it has no meaning for a single continuous
    metadata variable -- disable it under the metadata value source. The Figure DOES
    generalize (pseudobulk becomes the per-replicate mean of the variable, the violin is
    over all cells), so it stays available."""
    meta = (group_value_source == 'meta')
    options = [
        {'label': 'Figure: pseudobulk + per-cell violin', 'value': 'figure'},
        {'label': 'Violin (all cells)', 'value': 'violin'},
        {'label': 'Box', 'value': 'box'},
        {'label': 'Strip', 'value': 'strip'},
        {'label': 'Dot plot (gene set)', 'value': 'dotplot', 'disabled': meta},
    ]
    if meta and current_style == 'dotplot':
        return options, 'figure'
    return options, (current_style or 'figure')


# ---- F2: populate the group-by / split-by / gene-module dropdowns ----
@callback(
    [Output('group-by-select', 'options'),
     Output('group-by-select', 'value'),
     Output('group-split-select', 'options'),
     Output('group-split-select', 'value'),
     Output('gene-module-select', 'options'),
     Output('gene-module-select', 'value'),
     Output('group-replicate-select', 'options'),
     Output('group-replicate-select', 'value'),
     Output('group-gene-select', 'options', allow_duplicate=True),
     Output('group-gene-select', 'value'),
     Output('group-style', 'value', allow_duplicate=True),
     Output('group-positive-only', 'value', allow_duplicate=True),
     Output('group-meta-select', 'options'),
     Output('group-meta-select', 'value')],
    [Input('data-store', 'data')],
    [State('url', 'search'),
     State('group-value-source', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def populate_group_controls(data_store, url_search, group_value_source):
    empty_split = [{'label': '(none)', 'value': ''}]
    if not data_store:
        return [], None, empty_split, '', [], None, [], None, [], None, 'figure', ['enabled'], [], None

    # A shared URL may carry expression-by-group state; apply it only when the loaded
    # dataset matches the one in the URL, so a later manual dataset switch is not
    # silently re-overridden by stale shared state.
    state = (state_for_dataset(url_search, data_store.get('dataset_id')) if url_search else None) or {}

    column_types = data_store.get('column_types', {}) or {}
    categorical_cols = [c for c, t in column_types.items() if t == 'categorical']
    numeric_cols = [c for c, t in column_types.items() if t == 'numeric']

    # Demultiplex composite (e.g. library x genotype): a per-library demux donor label
    # ('genotype') is meaningless grouped on its own (donor labels repeat across
    # libraries), so offer the (library x genotype) COMPOSITE as a group key and drop the
    # trailing demux sub-label(s) from the standalone group/split choices. The leading
    # column (library) stays groupable on its own.
    replicate_columns = [c for c in (data_store.get('replicate_columns') or [])
                         if c in categorical_cols]
    demux_subcols = set(replicate_columns[1:]) if len(replicate_columns) >= 2 else set()
    composite_key = _replicate_key(replicate_columns) if len(replicate_columns) >= 2 else None
    composite_label = (' x '.join(replicate_columns) + ' (demux donor)') if composite_key else None
    standalone_cats = [c for c in categorical_cols if c not in demux_subcols]

    group_options = ([{'label': composite_label, 'value': composite_key}] if composite_key else []) \
        + [{'label': c, 'value': c} for c in standalone_cats]
    valid_group = ({composite_key} if composite_key else set()) | set(standalone_cats)

    # Default grouping: shared-URL value, else the configured annotation column, else
    # 'annotation_refined', else the first standalone categorical column (NOT the demux
    # composite -- cell type is the natural default).
    ann_col = data_store.get('annotation_column')
    if state.get('group_by') in valid_group:
        group_value = state['group_by']
    elif ann_col in standalone_cats:
        group_value = ann_col
    elif 'annotation_refined' in standalone_cats:
        group_value = 'annotation_refined'
    elif standalone_cats:
        group_value = standalone_cats[0]
    else:
        group_value = composite_key

    # Split-by includes a "(none)" option (value '') so a second grouping is optional.
    split_options = empty_split + group_options
    valid_split = {''} | valid_group
    split_value = state['group_split'] if state.get('group_split') in valid_split else ''

    gene_modules = data_store.get('gene_modules', {}) or {}
    module_options = [
        {'label': f"{name} ({len(genes)} genes)", 'value': name}
        for name, genes in gene_modules.items()
    ]
    # Default to the first module so the dot-plot style renders immediately when
    # chosen (the earlier "not working" was an empty module with no default).
    if state.get('gene_module') in gene_modules:
        module_value = state['gene_module']
    else:
        module_value = next(iter(gene_modules), None)

    # Pseudobulk replicate options: the configured composite (e.g. library x genotype)
    # first, then each standalone categorical column -- excluding the demux sub-label
    # ('genotype'), which on its own pools cells across libraries (the same mistake the
    # group-by guard prevents).
    replicate_options = []
    default_replicate = None
    if composite_key:
        replicate_options.append(
            {'label': ' x '.join(replicate_columns) + ' (demux replicate)', 'value': composite_key})
        default_replicate = composite_key
    replicate_options += [{'label': c, 'value': c} for c in standalone_cats]
    valid_replicate = {opt['value'] for opt in replicate_options}
    # Default replicate: shared-URL value, then the configured composite, then a known
    # sample/batch column, else None (Panel A omitted) rather than an arbitrary first
    # categorical column (which could be a barcode/cluster axis).
    SAMPLE_HINTS = ('library', 'orig.ident', 'sample', 'batch', 'donor', 'genotype')
    if state.get('group_replicate') in valid_replicate:
        replicate_value = state['group_replicate']
    elif default_replicate is not None:
        replicate_value = default_replicate
    else:
        replicate_value = next((c for c in categorical_cols if c.lower() in SAMPLE_HINTS), None)

    # Group-view gene: shared-URL value, else the dataset's default gene. (Options are
    # also kept in sync by update_group_gene_select.)
    genes = data_store.get('genes', []) or []
    gene_options = [{'label': g, 'value': g} for g in sorted(genes)]
    if state.get('group_gene') in genes:
        gene_value = state['group_gene']
    else:
        default_gene = data_store.get('default_gene')
        gene_value = default_gene if default_gene in genes else None

    # Style + positive-cell gate: shared-URL values, else the figure default. Only the
    # gene-set dot plot is invalid under a metadata value source (the Figure generalizes).
    if group_value_source == 'meta':
        valid_styles = {'figure', 'violin', 'box', 'strip'}
    else:
        valid_styles = {'figure', 'violin', 'box', 'strip', 'dotplot'}
    group_style = state['group_style'] if state.get('group_style') in valid_styles else 'figure'
    positive_only = state.get('group_positive_only')
    if not isinstance(positive_only, list):
        positive_only = ['enabled']

    # Continuous-metadata options for the "Metadata" value source: every numeric .obs
    # column. Default to a familiar QC / topographic metric when one is present.
    meta_options = [{'label': c, 'value': c} for c in numeric_cols]
    META_HINTS = ('nCount_RNA', 'nFeature_RNA', 'DV.Score', 'NT.Score')
    if state.get('group_meta') in numeric_cols:
        meta_value = state['group_meta']
    else:
        meta_value = next((c for c in META_HINTS if c in numeric_cols),
                          numeric_cols[0] if numeric_cols else None)

    return (group_options, group_value, split_options, split_value,
            module_options, module_value, replicate_options, replicate_value,
            gene_options, gene_value, group_style, positive_only,
            meta_options, meta_value)


# ---- F2: searchable gene dropdown for the group view (mirrors gene-select) ----
@callback(
    Output('group-gene-select', 'options', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('group-gene-select', 'search_value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_group_gene_select(data_store, search_value):
    if not data_store or 'genes' not in data_store:
        return []
    return gene_search_options(data_store['genes'], search_value)


# ---- F4: show/hide the second-gene dropdown ----
@callback(
    Output('gene-select-2-container', 'style'),
    Input('compare-genes', 'value')
)
def toggle_compare_genes(compare_genes):
    if _on(compare_genes):
        return {'display': 'block'}
    return {'display': 'none'}


# ---- F4: searchable second-gene dropdown (mirrors gene-select) ----
@callback(
    Output('gene-select-2', 'options'),
    [Input('data-store', 'data'),
     Input('gene-select-2', 'search_value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_gene_select_2(data_store, search_value):
    if not data_store or 'genes' not in data_store:
        return []
    return gene_search_options(data_store['genes'], search_value)