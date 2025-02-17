from dash import html, dcc
import dash_bootstrap_components as dbc
from components.sidebar.dataset_info import create_dataset_info
from components.sidebar.control_panel import create_control_panel

def create_sidebar():
    return html.Div([
        # Share View button always visible at top
        dbc.Button("Share View", id="share-button", color="primary", className="mb-3"),
        dbc.Input(id="share-url", type="text", style={'display': 'none'}, className="mt-2"),
        
        # Accordion for collapsible sections
        dbc.Accordion([
            dbc.AccordionItem([
                create_dataset_info()
            ], title="Dataset", item_id="dataset-section"),
            
            dbc.AccordionItem([
                create_control_panel()
            ], title="Visualization Control", item_id="viz-section"),
        ], start_collapsed=False, always_open=True, active_item=["dataset-section", "viz-section"])  # Set both sections to be active by default
    ], className="p-3") 