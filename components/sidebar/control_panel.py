from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.config import load_config

def create_control_panel():
    config = load_config()
    plot_settings = config.get('plot_settings', {})
    
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
                ),
                
                # Add additional controls based on config
                html.Div([
                    html.Label("Plot Settings:", className="mt-3"),
                    dcc.Slider(
                        id='opacity-slider',
                        min=0,
                        max=1,
                        step=0.1,
                        value=plot_settings.get('opacity', 0.7),
                        marks={i/10: str(i/10) for i in range(0, 11, 2)}
                    )
                ] if plot_settings.get('show_opacity_control', True) else [])
            ])
        ])
    ], className="mt-3") 