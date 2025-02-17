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
                
                # New container for custom embedding controls with better styling
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Custom Embedding Controls", className="card-subtitle mb-2 text-muted"),
                        html.Label("X-axis:", className="mt-2"),
                        dcc.Dropdown(
                            id='custom-x-select',
                            placeholder="Select X-axis metric"
                        ),
                        html.Label("Y-axis:", className="mt-2"),
                        dcc.Dropdown(
                            id='custom-y-select',
                            placeholder="Select Y-axis metric"
                        )
                    ])
                ], id='custom-embedding-container', 
                   className="mt-3",
                   style={'display': 'none', 'backgroundColor': '#f8f9fa'}),
                
                html.Label("Color by:", className="mt-3"),
                dcc.Dropdown(
                    id='color-select',
                    placeholder="Select feature"
                ),
                
                # Gene selection input with better styling
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Gene Selection", className="card-subtitle mb-2 text-muted"),
                        dcc.Dropdown(
                            id='gene-select',
                            placeholder="Type gene name..."
                        )
                    ])
                ], id='gene-select-container', 
                   className="mt-3",
                   style={'display': 'none', 'backgroundColor': '#f8f9fa'}),
                
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