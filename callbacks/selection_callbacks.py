from dash import Input, Output, State, callback
from utils.error_handling import handle_callback_error, log_callback_info

@callback(
    [Output('selection-store', 'data', allow_duplicate=True),
     Output('selection-info', 'children', allow_duplicate=True)],
    Input('main-plot', 'selectedData'),
    State('data-store', 'data'),
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_selection(selected_data, data_store):
    if not selected_data or not data_store:
        return None, "No cells selected"
    
    points = selected_data.get('points', [])
    n_selected = len(points)
    
    selection_data = {
        'indices': [p['pointIndex'] for p in points] if points else []
    }
    
    return selection_data, f"Selected {n_selected} cells" 