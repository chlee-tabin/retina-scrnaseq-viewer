from dash import Input, Output, State, callback, ctx, no_update
import json
import base64
from urllib.parse import parse_qs, urlencode
import logging
import dash
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import create_share_url, parse_url_state, encode_state

logger = logging.getLogger(__name__)

@callback(
    [Output('dataset-select', 'value', allow_duplicate=True),
     Output('embedding-select', 'value', allow_duplicate=True),
     Output('color-select', 'value', allow_duplicate=True),
     Output('viz-mode', 'value', allow_duplicate=True),
     Output('gene-select', 'value', allow_duplicate=True),
     Output('custom-x-select', 'value', allow_duplicate=True),
     Output('custom-y-select', 'value', allow_duplicate=True),
     Output('bin-number-slider', 'value', allow_duplicate=True),
     Output('percentile-slider', 'value', allow_duplicate=True),
     Output('enable-binning', 'value', allow_duplicate=True),
     Output('custom-embedding-container', 'style', allow_duplicate=True),
     Output('gene-select-container', 'style', allow_duplicate=True),
     Output('gene-select', 'options', allow_duplicate=True)],
    [Input('url', 'search'),
     Input('url', 'pathname')],
    [State('dataset-select', 'options'),
     State('dataset-select', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def initialize_from_url(search, pathname, dataset_options, current_dataset):
    triggered_id = ctx.triggered_id if ctx.triggered_id else 'url.search'
    
    if triggered_id == 'url.pathname':
        # Pure pathname navigation (no shared state): leave controls untouched.
        return (no_update,) * 13
    
    if not search:
        # No shared-URL state: leave every control at its layout/dataset default.
        return (no_update,) * 13
    
    state = parse_url_state(search)
    if not state:
        # No shared-URL state: leave every control at its layout/dataset default.
        return (no_update,) * 13
    
    # Get custom embedding values and binning parameters
    custom_x = state.get('custom_x') if state.get('embedding') == 'custom_embedding' else None
    custom_y = state.get('custom_y') if state.get('embedding') == 'custom_embedding' else None
    bin_number = state.get('bins', 50)
    percentile = state.get('percentile', 0.95)
    enable_binning = state.get('enable_binning', [])
    
    # Show custom embedding container if custom embedding is selected
    container_style = {'display': 'block'} if state.get('embedding') == 'custom_embedding' else {'display': 'none'}
    
    # Determine gene select container visibility
    gene_container_style = {'display': 'block'} if state.get('color') == 'gene_expression' else {'display': 'none'}
    
    return (
        state.get('dataset'),
        state.get('embedding'),
        state.get('color'),
        state.get('mode', 'random'),
        state.get('gene'),
        custom_x,
        custom_y,
        bin_number,
        percentile,
        enable_binning,
        container_style,
        gene_container_style,
        []  # Empty options, will be populated by the gene_select callback
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
     State('main-plot', 'relayoutData'),
     State('url', 'href')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def share_url(n_clicks, dataset, embedding, color_by, gene, viz_mode, 
              custom_x, custom_y, bin_number, percentile, enable_binning,
              relay_data, current_url):
    if n_clicks is None:
        return {'display': 'none'}, ''
    
    state_dict = {
        'dataset': dataset,
        'embedding': embedding,
        'color': color_by,
        'mode': viz_mode
    }
    
    if color_by == 'gene_expression' and gene:
        state_dict['gene'] = gene
    
    if embedding == 'custom_embedding':
        state_dict['custom_x'] = custom_x
        state_dict['custom_y'] = custom_y
        state_dict['bins'] = bin_number
        state_dict['percentile'] = percentile
        state_dict['enable_binning'] = enable_binning
    
    if relay_data:
        view_state = {}
        if 'xaxis.range[0]' in relay_data and 'xaxis.range[1]' in relay_data:
            view_state['xrange'] = [relay_data['xaxis.range[0]'], relay_data['xaxis.range[1]']]
        if 'yaxis.range[0]' in relay_data and 'yaxis.range[1]' in relay_data:
            view_state['yrange'] = [relay_data['yaxis.range[0]'], relay_data['yaxis.range[1]']]
        if view_state:
            state_dict['view'] = view_state
    
    base_url = current_url.split('?')[0]
    share_url = create_share_url(base_url, state_dict)
    
    return {'display': 'block', 'width': '100%', 'marginTop': '10px'}, share_url 