from dash import Input, Output, State, callback
import dash
from utils.server_monitoring import get_memory_usage, update_session
import time

@callback(
    Output('server-status', 'children'),
    Input('interval-component', 'n_intervals'),
    State('session-id', 'data'),
    prevent_initial_call=False
)
def update_status(n, session_id):
    active_users = update_session(session_id)
    memory_usage = get_memory_usage()
    return f"Active users: {active_users} | Memory usage: {memory_usage}" 