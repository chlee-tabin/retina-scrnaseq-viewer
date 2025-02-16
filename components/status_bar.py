from dash import html
import dash_bootstrap_components as dbc

def create_status_bar():
    return html.Div(
        dbc.Row([
            dbc.Col(
                html.Small(id='server-status', className='text-muted'),
                width='auto'
            )
        ], justify='end'),
        style={'padding': '2px 10px', 'backgroundColor': '#f8f9fa'}
    ) 