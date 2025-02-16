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
     Output('gene-select', 'value', allow_duplicate=True)],
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
        return no_update, no_update, no_update, no_update, no_update
    
    if not search:
        # Return defaults instead of no_update for viz_mode
        return no_update, no_update, no_update, 'random', no_update
    
    state = parse_url_state(search)
    if not state:
        return no_update, no_update, no_update, 'random', no_update
    
    return (
        state.get('dataset'),
        state.get('embedding'),
        state.get('color'),
        state.get('mode', 'random'),  # Default to 'random'
        state.get('gene')
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
     State('main-plot', 'relayoutData'),
     State('url', 'href')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def share_url(n_clicks, dataset, embedding, color_by, gene, viz_mode, relay_data, current_url):
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
    
    if relay_data:
        view_state = {}
        if 'xaxis.range[0]' in relay_data and 'xaxis.range[1]' in relay_data:
            view_state['xrange'] = [relay_data['xaxis.range[0]'], relay_data['xaxis.range[1]']]
        if 'yaxis.range[0]' in relay_data and 'yaxis.range[1]' in relay_data:
            view_state['yrange'] = [relay_data['yaxis.range[0]'], relay_data['yaxis.range[1]']]
        state_dict['view'] = view_state
    
    share_url = create_share_url(current_url.split('?')[0], state_dict)
    return {'display': 'block'}, share_url 