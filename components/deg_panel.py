"""Main-panel section for the ROI differential-expression view (shown when the
sidebar View is set to "ROI differential expression"). Two topographic panels --
lasso ROI A in the left, ROI B in the right -- then a resolved-selection recap and
the volcano below. Hidden by default; toggled by callbacks.deg_callbacks.toggle_deg_view.
"""
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

# Expose lasso + box select on the modebar (and make lasso the default drag mode in the
# figure layout) so a region can be drawn without first hunting through the toolbar.
_ROI_GRAPH_CONFIG = {'modeBarButtonsToAdd': ['select2d', 'lasso2d'], 'displaylogo': False}


def create_deg_panel():
    return html.Div([
        dbc.Alert([
            html.B("ROI differential expression (pseudobulk). "),
            "Lasso a region in each panel. ", html.B("One ROI"), " = region vs. the rest; ",
            html.B("two ROIs"), " = A vs. B. An overlapping cell is assigned to the ",
            "smaller ROI. The test is a negative-binomial pseudobulk DE (pydeseq2) over the "
            "demux replicate unit (library x genotype) — the interactive form of the "
            "manuscript area-DEG volcano. Available where the dataset has DV/NT scores and "
            "a configured replicate unit.",
        ], color="light", className="border small mb-2"),

        dbc.Row([
            dbc.Col([
                html.Div("ROI A", className="fw-semibold text-center"),
                dcc.Graph(id='roi-plot-a', style={'height': '42vh'}, config=_ROI_GRAPH_CONFIG),
                html.Div(id='roi-a-info', className='small text-muted text-center'),
            ], width=6),
            dbc.Col([
                html.Div("ROI B  (leave empty → contrast vs. rest)",
                         className="fw-semibold text-center"),
                dcc.Graph(id='roi-plot-b', style={'height': '42vh'}, config=_ROI_GRAPH_CONFIG),
                html.Div(id='roi-b-info', className='small text-muted text-center'),
            ], width=6),
        ]),

        html.Div([
            dbc.Button("Run differential expression", id='run-deg-btn',
                       color='primary', className='me-2'),
            dbc.Button("Clear ROIs", id='clear-roi-btn', color='secondary',
                       outline=True, className='me-2'),
            dbc.Button("Download CSV", id='deg-download-btn', color='success', outline=True),
        ], className='my-2'),

        html.Div(id='deg-status', className='small my-2'),
        dcc.Loading(dcc.Graph(id='deg-recap-plot', style={'height': '34vh'}), type='circle'),
        dcc.Loading(dcc.Graph(id='volcano-plot', style={'height': '55vh'}), type='circle'),
        dash_table.DataTable(
            id='deg-table', page_size=15, sort_action='native', filter_action='native',
            style_table={'overflowX': 'auto'},
            style_cell={'fontSize': 12, 'fontFamily': 'monospace', 'textAlign': 'left'},
            style_header={'fontWeight': 'bold'},
        ),
    ], id='deg-view-container', style={'display': 'none'}, className='p-3')
