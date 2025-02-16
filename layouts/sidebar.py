from dash import html, dcc
import dash_bootstrap_components as dbc
from components.sidebar.dataset_info import create_dataset_info
from components.sidebar.control_panel import create_control_panel

def create_sidebar():
    return html.Div([
        html.H1("Single-cell Data Viewer", className="mb-4"),
        dbc.Card([
            dbc.CardBody([
                # Dataset selection
                html.Label("Select Dataset:"),
                dcc.Dropdown(
                    id='dataset-select',
                    placeholder="Choose a dataset"
                ),
                
                # Dataset info panel
                create_dataset_info(),
                
                # Share button and URL
                dbc.Button("Share View", id="share-button", color="primary", className="mt-3"),
                dbc.Input(id="share-url", type="text", style={'display': 'none'}, className="mt-2"),
                
                # Control panel
                create_control_panel()
            ])
        ])
    ], className="p-3") 