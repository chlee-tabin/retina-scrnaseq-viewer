from dash import Input, Output, State, callback, ctx, no_update
import json
import base64
from urllib.parse import parse_qs, urlencode
import logging
import dash
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import create_share_url, parse_url_state, encode_state

logger = logging.getLogger(__name__)

# Restore controls from a shared URL. This callback owns only the controls no other
# callback restores: it sets the dataset (which triggers update_data, the url-aware
# owner of embedding/colour/gene), plus the global view controls. The custom axes are
# restored by update_custom_embedding_controls and every expression-by-group control
# by populate_group_controls -- each reads the same URL state. prevent_initial_call=
# 'initial_duplicate' lets this fire on the very first render so a pasted link works.
@callback(
    [Output('dataset-select', 'value', allow_duplicate=True),
     Output('viz-mode', 'value', allow_duplicate=True),
     Output('bin-number-slider', 'value', allow_duplicate=True),
     Output('percentile-slider', 'value', allow_duplicate=True),
     Output('enable-binning', 'value', allow_duplicate=True),
     Output('enable-smoothing', 'value', allow_duplicate=True),
     Output('bin-stat', 'value', allow_duplicate=True),
     Output('plot-type', 'value', allow_duplicate=True),
     Output('compare-genes', 'value', allow_duplicate=True),
     Output('gene-select-2', 'value', allow_duplicate=True),
     Output('compare-shared-scale', 'value', allow_duplicate=True)],
    [Input('url', 'search'),
     Input('url', 'pathname')],
    prevent_initial_call='initial_duplicate'
)
@handle_callback_error
@log_callback_info
def initialize_from_url(search, pathname):
    n_out = 11
    triggered_id = ctx.triggered_id
    # Pure pathname navigation (no new shared state): leave controls untouched.
    if triggered_id == 'url.pathname' or not search:
        return (no_update,) * n_out
    state = parse_url_state(search)
    if not state:
        return (no_update,) * n_out

    return (
        state.get('dataset', no_update),
        state.get('mode', 'random'),
        state.get('bins', no_update),
        state.get('percentile', no_update),
        state.get('enable_binning', no_update),
        state.get('enable_smoothing', no_update),
        state.get('bin_stat', no_update),
        state.get('plot_type', no_update),
        state.get('compare_genes', no_update),
        state.get('gene2', no_update),
        state.get('compare_shared_scale', no_update),
    )

@callback(
    [Output('share-url', 'style', allow_duplicate=True),
     Output('share-url', 'value', allow_duplicate=True)],
    Input('share-button', 'n_clicks'),
    [State('dataset-select', 'value'),
     State('embedding-select', 'value'),
     State('color-select', 'value'),
     State('gene-select', 'value'),
     State('viz-mode', 'value'),
     State('custom-x-select', 'value'),
     State('custom-y-select', 'value'),
     State('bin-number-slider', 'value'),
     State('percentile-slider', 'value'),
     State('enable-binning', 'value'),
     State('enable-smoothing', 'value'),
     State('bin-stat', 'value'),
     State('compare-genes', 'value'),
     State('gene-select-2', 'value'),
     State('compare-shared-scale', 'value'),
     State('plot-type', 'value'),
     State('group-gene-select', 'value'),
     State('group-by-select', 'value'),
     State('group-split-select', 'value'),
     State('group-style', 'value'),
     State('gene-module-select', 'value'),
     State('group-positive-only', 'value'),
     State('group-replicate-select', 'value'),
     State('url', 'href')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def share_url(n_clicks, dataset, embedding, color_by, gene, viz_mode,
              custom_x, custom_y, bin_number, percentile, enable_binning,
              enable_smoothing, bin_stat, compare_genes, gene2, compare_shared_scale,
              plot_type, group_gene, group_by, group_split, group_style,
              gene_module, group_positive_only, group_replicate,
              current_url):
    if n_clicks is None:
        return {'display': 'none'}, ''

    # Share-state schema v2 (see docs/VIEWER_SPEC.md). Only the keys relevant to the
    # current view are written, so URLs stay lean; absent keys restore to defaults.
    state_dict = {
        'v': 2,
        'dataset': dataset,
        'embedding': embedding,
        'color': color_by,
        'mode': viz_mode,
        'plot_type': plot_type,
    }

    if color_by == 'gene_expression' and gene:
        state_dict['gene'] = gene

    # Two-gene side-by-side comparison.
    if compare_genes:
        state_dict['compare_genes'] = compare_genes
        if gene2:
            state_dict['gene2'] = gene2
        if compare_shared_scale:
            state_dict['compare_shared_scale'] = compare_shared_scale

    # Spatial-map / binning controls (only meaningful on the custom DV/NT embedding).
    if embedding == 'custom_embedding':
        state_dict['custom_x'] = custom_x
        state_dict['custom_y'] = custom_y
        state_dict['bins'] = bin_number
        state_dict['percentile'] = percentile
        state_dict['enable_binning'] = enable_binning
        state_dict['enable_smoothing'] = enable_smoothing
        state_dict['bin_stat'] = bin_stat

    # Expression-by-group controls.
    if plot_type == 'group':
        state_dict['group_by'] = group_by
        state_dict['group_split'] = group_split
        state_dict['group_gene'] = group_gene
        state_dict['group_style'] = group_style
        state_dict['gene_module'] = gene_module
        state_dict['group_positive_only'] = group_positive_only
        state_dict['group_replicate'] = group_replicate

    base_url = current_url.split('?')[0]
    url_value = create_share_url(base_url, state_dict)

    return {'display': 'block', 'width': '100%', 'marginTop': '10px'}, url_value