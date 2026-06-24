from dash import html, dcc
import dash_bootstrap_components as dbc
from components.deg_panel import create_deg_panel

def create_main_panel():
    return html.Div([
        # ROI-DEG state: polygon vertices ({A:[[x,y]...], B:[...]}), the active draw mode
        # ('A'/'B'/None), the computed results, and the CSV download target.
        dcc.Store(id='roi-vertices-store'),
        dcc.Store(id='roi-draw-mode-store'),
        dcc.Store(id='deg-results-store'),
        dcc.Download(id='deg-download'),

        dcc.Loading(
            id="loading-upload",
            type="circle",
            children=[
                html.Div(id="loading-output"),
                dcc.Graph(
                    id='main-plot',
                    style={'height': '85vh'}  # Make the plot larger
                ),
                html.Div(id='selection-info', className='mt-3')
            ]
        ),

        # ROI differential-expression section (drawn on the map above; hidden in group view).
        create_deg_panel(),
    ], className="p-3")
