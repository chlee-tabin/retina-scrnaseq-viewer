from dash import html, dcc
import dash_bootstrap_components as dbc
from components.deg_panel import create_deg_panel

def create_main_panel():
    return html.Div([
        # ROI-DEG state: captured ROIs, computed results, and the CSV download target.
        dcc.Store(id='deg-roi-store'),
        dcc.Store(id='deg-results-store'),
        dcc.Download(id='deg-download'),

        # Map / group view (the main-plot). Hidden when the ROI-DEG view is active.
        html.Div([
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
            )
        ], id='map-view-container'),

        # ROI differential-expression view (two topographic panels + volcano).
        create_deg_panel(),
    ], className="p-3")
