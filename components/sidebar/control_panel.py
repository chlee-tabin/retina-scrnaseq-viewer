from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.config import load_config

def create_control_panel():
    config = load_config()
    
    return html.Div([
        dbc.Card([
            dbc.CardBody([
                html.H5("Visualization Controls", className="card-title"),
                
                html.Label("Embedding:"),
                dcc.Dropdown(
                    id='embedding-select',
                    placeholder="Select embedding"
                ),
                
                html.Label("Color by:", className="mt-3"),
                dcc.Dropdown(
                    id='color-select',
                    placeholder="Select feature"
                ),
                
                html.Label("Visualization Mode:", className="mt-3"),
                dcc.RadioItems(
                    id='viz-mode',
                    options=[
                        {'label': 'Single cells', 'value': 'cells'},
                        {'label': 'Metacells', 'value': 'metacells'}
                    ],
                    value='cells'
                )
            ])
        ])
    ], className="mt-3") 