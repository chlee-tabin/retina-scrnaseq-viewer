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
                
                # Gene selection input
                html.Div([
                    html.Label("Gene:", className="mt-3"),
                    dcc.Dropdown(
                        id='gene-select',
                        placeholder="Type gene name..."
                    )
                ], id='gene-select-container', style={'display': 'none'}),
                
                html.Label("Plot Order:", className="mt-3"),
                dcc.RadioItems(
                    id='viz-mode',
                    options=[
                        {'label': 'Random', 'value': 'random'},
                        {'label': 'Ascending', 'value': 'ordered_asc'},
                        {'label': 'Descending', 'value': 'ordered_desc'}
                    ],
                    value='random'
                )
            ])
        ])
    ], className="mt-3") 