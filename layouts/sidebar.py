from dash import html, dcc
import dash_bootstrap_components as dbc
from components.sidebar.dataset_info import create_dataset_info
from components.sidebar.control_panel import create_control_panel
from utils.version import VIEWER_VERSION, viewer_revision

# Viewer version (bump on each released change) + manuscript preprint link, shown at
# the top of the sidebar.
PREPRINT_URL = "https://doi.org/10.64898/2026.01.04.697548"
SOURCE_URL = "https://github.com/chlee-tabin/retina-scrnaseq-viewer"
ANALYSIS_URL = "https://github.com/chlee-tabin/retina-spatial-scrna-analysis"


def _link(text, href):
    return html.A(text, href=href, target="_blank", rel="noopener noreferrer")

def create_sidebar():
    return html.Div([
        # Permanent header: viewer name + version, with a link to the manuscript preprint.
        html.Div([
            html.Div(f"Retinal scRNA-seq topographic viewer v{VIEWER_VERSION}",
                     className="fw-semibold"),
            html.A(
                "📄 Manuscript preprint (bioRxiv)",
                href=PREPRINT_URL, target="_blank", rel="noopener noreferrer",
                className="small",
            ),
        ], className="mb-2"),
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
        ], start_collapsed=False, always_open=True, active_item=["dataset-section", "viz-section"]),
        html.Footer([
            f"viewer {VIEWER_VERSION} · {viewer_revision()} · ",
            _link("Source code (MIT license)", SOURCE_URL), " · ",
            _link("Analysis code", ANALYSIS_URL),
        ], className='small text-muted mt-3'),
    ], className="p-3")
