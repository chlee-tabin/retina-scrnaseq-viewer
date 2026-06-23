from dash import html, dcc
import dash_bootstrap_components as dbc
from components.sidebar.dataset_info import create_dataset_info
from components.sidebar.control_panel import create_control_panel

# Permanent link to the manuscript preprint, shown at the top of the sidebar.
PREPRINT_URL = "https://www.biorxiv.org/content/10.64898/2026.01.04.697548v1"

def create_sidebar():
    return html.Div([
        # Permanent header: link to the manuscript preprint (bioRxiv).
        html.Div(
            html.A(
                "📄 Manuscript preprint (bioRxiv)",
                href=PREPRINT_URL, target="_blank", rel="noopener noreferrer",
                className="fw-semibold",
            ),
            className="mb-2",
        ),
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