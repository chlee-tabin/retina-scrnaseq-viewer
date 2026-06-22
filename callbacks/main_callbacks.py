from dash import Input, Output, State, callback
import json
import base64
from urllib.parse import parse_qs, urlencode
from utils.data_loading import load_adata, load_dataset_config
from utils.plotting import (
    create_scatter_plot,
    create_metacell_plot,
    create_binned_plot,
    create_group_expression_plot,
    create_group_expression_figure,
    create_dotplot,
    create_dual_gene_figure,
)
import pandas as pd
import logging
from utils.processing import create_metacells, calculate_selection_stats
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import parse_url_state
import scipy.sparse
import numpy as np
import traceback
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def _message_figure(text):
    """Return an empty plot displaying a centered, user-friendly message."""
    fig = go.Figure()
    fig.add_annotation(
        text=text,
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=16, color="#666")
    )
    fig.update_layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white"
    )
    return fig


def _gene_vector(adata, gene):
    """Return a dense 1D expression vector for `gene` (handles sparse X)."""
    sub = adata[:, gene].X
    if scipy.sparse.issparse(sub):
        return sub.toarray().flatten()
    return np.asarray(sub).flatten()


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
     Input('compare-shared-scale', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, custom_x, custom_y, color_by, gene, viz_mode,
                bin_number, percentile, enable_binning, enable_smoothing, selection_data,
                url_search, plot_type, group_gene, group_by, group_split, group_style,
                gene_module, compare_genes, gene2, bin_stat,
                group_positive_only, group_replicate, compare_shared_scale):
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

            # Consistent per-type colour + biological order, applied when grouping on
            # the dataset's configured annotation column (else default palette/order).
            g_cmap, g_corder = _annotation_style(data_store, group_by)

            # Dot plot of a gene set (module) x groups. A module selection is
            # self-sufficient (no single gene required); without a module, fall
            # back to the chosen gene as a 1-gene dot plot.
            if group_style == 'dotplot':
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
                expr_df[group_by] = adata.obs[group_by].to_numpy()
                return create_dotplot(expr_df, genes, group_by, category_order=g_corder)

            # Gene-based styles (figure / violin / box / strip) need a single gene.
            if not group_gene:
                return _message_figure("Select a gene to plot expression by group.")
            if group_gene not in adata.var_names:
                return _message_figure(
                    f"Gene '{group_gene}' not found in this dataset. "
                    "Try a different gene or check the name/casing."
                )

            # Polished 2-panel figure: pseudobulk per (group x replicate) + positive-
            # cell log1p violin. Panel A reconstructs raw counts from .X + nCount_RNA.
            if group_style == 'figure':
                fig_df = pd.DataFrame({'expr': _gene_vector(adata, group_gene)})
                fig_df[group_by] = adata.obs[group_by].to_numpy()
                if 'nCount_RNA' in adata.obs.columns:
                    fig_df['ncount'] = adata.obs['nCount_RNA'].to_numpy()
                rep_series = _replicate_series(
                    adata, group_replicate, data_store.get('replicate_columns', []))
                if rep_series is not None:
                    fig_df['replicate'] = rep_series
                positive_only = (isinstance(group_positive_only, (list, tuple))
                                 and 'enabled' in group_positive_only)
                return create_group_expression_figure(
                    fig_df, group_gene, group_by,
                    color_map=g_cmap, category_order=g_corder,
                    positive_only=positive_only,
                    norm_target=data_store.get('x_norm_target', 1e4),
                )

            # Simple violin / box / strip of one gene across groups.
            df = pd.DataFrame({'expr': _gene_vector(adata, group_gene)})
            df[group_by] = adata.obs[group_by].to_numpy()
            split_col = group_split if (group_split and group_split in adata.obs.columns) else None
            if split_col:
                df[split_col] = adata.obs[split_col].to_numpy()
            return create_group_expression_plot(
                df, group_gene, group_by, split_by=split_col, style=group_style,
                color_map=g_cmap, category_order=g_corder,
            )

        # ---- Embedding / spatial map view ----
        if not embedding:
            return {}

        # Get coordinates based on embedding type
        if embedding == 'custom_embedding':
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
            and bool(enable_binning) and 'enabled' in enable_binning
        )

        # ---- Two genes side-by-side (map view, gene-expression color) ----
        compare_on = isinstance(compare_genes, (list, tuple)) and 'enabled' in compare_genes
        if color_by == 'gene_expression' and compare_on and gene and gene2:
            missing = [g for g in (gene, gene2) if g not in adata.var_names]
            if missing:
                return _message_figure(
                    f"Gene(s) not found in this dataset: {', '.join(missing)}. "
                    "Check the name/casing."
                )
            smooth_on = bool(enable_smoothing) and 'enabled' in enable_smoothing
            smooth_sigma = float(data_store.get('smooth_sigma', 0) or 0) if smooth_on else 0
            return create_dual_gene_figure(
                np.asarray(x), np.asarray(y),
                [_gene_vector(adata, gene), _gene_vector(adata, gene2)],
                [gene, gene2],
                embedding,
                binned=binning_on,
                bin_size=bin_number,
                percentile=percentile,
                smooth_sigma=smooth_sigma,
                min_cells=data_store.get('min_cells_per_bin', 1),
                color_floor=data_store.get('color_floor', 0.05),
                bin_stat=(bin_stat or 'mean'),
                shared_scale=(isinstance(compare_shared_scale, (list, tuple))
                              and 'enabled' in compare_shared_scale),
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

        # Create DataFrame for plotting
        df = pd.DataFrame({
            'x': x,
            'y': y,
            'color': color_series
        })

        # Create the plot. Only the custom spatial embedding is ever binned;
        # all other embeddings always use the per-cell scatter (F1).
        if binning_on:
            smooth_on = bool(enable_smoothing) and 'enabled' in enable_smoothing
            smooth_sigma = float(data_store.get('smooth_sigma', 0) or 0) if smooth_on else 0
            min_cells = data_store.get('min_cells_per_bin', 1)
            fig = create_binned_plot(
                df, embedding, color_by,
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
                df, embedding, color_by,
                treat_as_categorical=treat_as_categorical,
                color_map=s_cmap, category_order=s_corder,
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
        state = parse_url_state(url_search)
        # Only honour shared custom-axis state when it belongs to the loaded dataset
        # (a stale ?state= from another dataset must not re-apply on a switch).
        if (state and state.get('embedding') == 'custom_embedding'
                and state.get('dataset') == data_store.get('dataset_id')):
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
    if embedding == 'custom_embedding' and enable_binning and 'enabled' in enable_binning:
        return {'display': 'none'}
    return {'display': 'block'}

@callback(
    Output('binning-controls', 'style'),
    [Input('enable-binning', 'value')]
)
def toggle_binning_controls(enable_binning):
    if enable_binning and 'enabled' in enable_binning:
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
     Output('group-positive-only', 'value', allow_duplicate=True)],
    Input('data-store', 'data'),
    State('url', 'search'),
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def populate_group_controls(data_store, url_search):
    empty_split = [{'label': '(none)', 'value': ''}]
    if not data_store:
        return [], None, empty_split, '', [], None, [], None, [], None, 'figure', ['enabled']

    # A shared URL may carry expression-by-group state; apply it only when the loaded
    # dataset matches the one in the URL, so a later manual dataset switch is not
    # silently re-overridden by stale shared state.
    state = parse_url_state(url_search) if url_search else None
    if not state or state.get('dataset') != data_store.get('dataset_id'):
        state = {}

    column_types = data_store.get('column_types', {}) or {}
    categorical_cols = [c for c, t in column_types.items() if t == 'categorical']
    group_options = [{'label': c, 'value': c} for c in categorical_cols]

    # Default grouping: shared-URL value, else the configured annotation column, else
    # 'annotation_refined', else the first categorical column.
    ann_col = data_store.get('annotation_column')
    if state.get('group_by') in categorical_cols:
        group_value = state['group_by']
    elif ann_col in categorical_cols:
        group_value = ann_col
    elif 'annotation_refined' in categorical_cols:
        group_value = 'annotation_refined'
    elif categorical_cols:
        group_value = categorical_cols[0]
    else:
        group_value = None

    # Split-by includes a "(none)" option (value '') so a second grouping is optional.
    split_options = empty_split + group_options
    valid_split = {''} | set(categorical_cols)
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

    # Pseudobulk replicate options: the configured composite (e.g. library x
    # genotype) first, then each categorical column on its own.
    replicate_columns = [c for c in (data_store.get('replicate_columns') or [])
                         if c in categorical_cols]
    replicate_options = []
    default_replicate = None
    if len(replicate_columns) >= 2:
        composite = _replicate_key(replicate_columns)
        replicate_options.append(
            {'label': ' × '.join(replicate_columns) + ' (demux replicate)', 'value': composite})
        default_replicate = composite
    replicate_options += [{'label': c, 'value': c} for c in categorical_cols]
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

    # Style + positive-cell gate: shared-URL values, else the figure defaults.
    valid_styles = {'figure', 'violin', 'box', 'strip', 'dotplot'}
    group_style = state['group_style'] if state.get('group_style') in valid_styles else 'figure'
    positive_only = state.get('group_positive_only')
    if not isinstance(positive_only, list):
        positive_only = ['enabled']

    return (group_options, group_value, split_options, split_value,
            module_options, module_value, replicate_options, replicate_value,
            gene_options, gene_value, group_style, positive_only)


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
    return _gene_search_options(data_store['genes'], search_value)


# ---- F4: show/hide the second-gene dropdown ----
@callback(
    Output('gene-select-2-container', 'style'),
    Input('compare-genes', 'value')
)
def toggle_compare_genes(compare_genes):
    if compare_genes and 'enabled' in compare_genes:
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
    return _gene_search_options(data_store['genes'], search_value)


def _gene_search_options(genes, search_value):
    """Fuzzy starts-with/contains gene search shared by the extra gene dropdowns.

    Mirrors the logic in callbacks.dataset_callbacks.update_gene_select: when a
    search string is present, starts-with matches sort ahead of contains
    matches; otherwise all genes are returned alphabetically.
    """
    if search_value:
        sv = search_value.lower()
        starts_with = []
        contains = []
        for gene in genes:
            gl = gene.lower()
            if gl.startswith(sv):
                starts_with.append(gene)
            elif sv in gl:
                contains.append(gene)
        matching = sorted(starts_with) + sorted(contains)
        return [{'label': g, 'value': g} for g in matching]
    return [{'label': g, 'value': g} for g in sorted(genes)]