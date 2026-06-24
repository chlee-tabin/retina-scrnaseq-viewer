"""ROI differential-expression section, shown below the spatial map (not a separate
view). The user draws one or two polygon ROIs directly on the map above (click to add
vertices, double-click to close; handled by assets/roi_draw.js), then runs a pseudobulk
NB-GLM DE over the demux replicate unit. Hidden in the "Expression by group" view.
"""
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc


def create_deg_panel():
    return html.Div([
        html.Hr(),
        html.H5("Differential expression (ROI)", className="mb-1"),
        dbc.Alert([
            "Click ", html.B("Draw ROI A"), ", then click on the map above to place polygon "
            "vertices and ", html.B("double-click to close"), ". One ROI = region vs. the "
            "rest; add ", html.B("ROI B"), " for A vs. B (an overlapping cell goes to the "
            "smaller ROI). Selection is point-in-polygon on the DV/NT scores, so it works on "
            "the per-cell or binned/smoothed map. The test is a negative-binomial pseudobulk "
            "DE (pydeseq2) over the demux replicate unit (library × genotype) — the "
            "interactive form of the manuscript area-DEG volcano.",
        ], color="light", className="border small mb-2"),

        dbc.Row([
            dbc.Col(dbc.ButtonGroup([
                dbc.Button("Draw ROI A", id='roi-draw-a-btn', color='primary', outline=True),
                dbc.Button("Draw ROI B", id='roi-draw-b-btn', color='warning', outline=True),
                dbc.Button("Clear", id='roi-clear-btn', color='secondary', outline=True),
            ]), width="auto"),
            dbc.Col(html.Span("Min cells/pseudobulk:", className="small me-1"),
                    width="auto", className="d-flex align-items-center"),
            dbc.Col(dcc.Input(id='deg-min-cells', type='number', value=50, min=5, step=5,
                              debounce=True, style={'width': '80px'}),
                    width="auto", className="d-flex align-items-center"),
            dbc.Col(html.Span(id='roi-draw-status', className='small text-muted'),
                    className="d-flex align-items-center"),
        ], className="g-2 mb-2 align-items-center"),

        html.Div([
            dbc.Button("Run differential expression", id='run-deg-btn',
                       color='primary', className='me-2'),
            dbc.Button("Download CSV", id='deg-download-btn', color='success', outline=True),
        ], className='my-2'),

        # Animated indeterminate activity bar (pydeseq2's fit is one long opaque step);
        # shown/hidden + the Run button disabled by clientside callbacks in deg_callbacks.
        html.Div(
            dbc.Progress(value=100, striped=True, animated=True, color="info",
                         label="Computing differential expression… (large ROIs ~1 min)",
                         style={'height': '22px'}),
            id='deg-progress-wrap', style={'display': 'none'}, className='my-2'),

        html.Div(id='deg-status', className='small my-2'),
        dcc.Graph(id='deg-recap-plot', style={'height': '34vh'}),
        dcc.Graph(id='volcano-plot', style={'height': '55vh'}),
        dash_table.DataTable(
            id='deg-table', page_size=15, sort_action='native', filter_action='native',
            style_table={'overflowX': 'auto'},
            style_cell={'fontSize': 12, 'fontFamily': 'monospace', 'textAlign': 'left'},
            style_header={'fontWeight': 'bold'}),

        # Hidden sinks for the clientside drawing callbacks (setMode / syncVerts / redraw).
        html.Div(id='roi-mode-dummy', style={'display': 'none'}),
        html.Div(id='roi-verts-dummy', style={'display': 'none'}),
        html.Div(id='roi-verts-dummy2', style={'display': 'none'}),
    ], id='deg-section', className='mt-2')
