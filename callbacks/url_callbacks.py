from dash import Input, Output, State, callback, ctx, no_update
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import (create_share_url, parse_url_state, share_states, schema_keys,
                         restore_outputs, restore_values, capture_state)
from utils.data_loading import get_dataset_config


# Option-populating owners restore their controls from the same STATE_SCHEMA.
@callback(
    restore_outputs('url'),
    [Input('url', 'search'), Input('url', 'pathname')],
    prevent_initial_call='initial_duplicate',
)
@handle_callback_error
@log_callback_info
def initialize_from_url(search, pathname):
    if not search:
        return (no_update,) * len(schema_keys('url'))
    triggered = getattr(ctx, 'triggered_prop_ids', {})
    if 'url.pathname' in triggered and 'url.search' not in triggered:
        return (no_update,) * len(schema_keys('url'))
    state = parse_url_state(search)
    dataset = get_dataset_config(state.get('dataset')) if state else None
    if dataset is None:
        return (no_update,) * len(schema_keys('url'))
    return restore_values(state, 'url', dataset)


@callback(
    [Output('share-url', 'style', allow_duplicate=True),
     Output('share-url', 'value', allow_duplicate=True)],
    Input('share-button', 'n_clicks'),
    share_states() + [State('url', 'href')],
    prevent_initial_call=True,
)
@handle_callback_error
@log_callback_info
def share_url(n_clicks, *values):
    if not n_clicks:
        return {'display': 'none'}, ''
    controls, current_url = values[:-1], values[-1]
    dataset_id = dict(zip(schema_keys(), controls)).get('dataset')
    dataset = get_dataset_config(dataset_id)
    if dataset is None:
        return no_update, no_update
    # Share-state schema v4: capture every control, including smoothing mode.
    state = capture_state(controls, dataset)
    url = create_share_url(current_url.split('?')[0], state)
    return {'display': 'block', 'width': '100%', 'marginTop': '10px'}, url
