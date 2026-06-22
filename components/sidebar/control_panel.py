from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.config import load_config

def create_control_panel():
    config = load_config()

    return html.Div([
        dbc.Card([
            dbc.CardBody([
                html.H5("Visualization Controls", className="card-title"),

                # Top-level view switch: the spatial/embedding map vs. the
                # "expression by group" (violin/box/strip/dotplot) view.
                html.Label("View:"),
                dcc.RadioItems(
                    id='plot-type',
                    options=[
                        {'label': 'Embedding / spatial map', 'value': 'map'},
                        {'label': 'Expression by group', 'value': 'group'},
                    ],
                    value='map',
                    className="mb-3",
                    labelStyle={'display': 'block'},
                ),

                # ---- MAP CONTROLS (embedding / spatial map view) ----
                html.Div([
                    html.Label("Embedding:"),
                    dcc.Dropdown(
                        id='embedding-select',
                        value='custom_embedding',  # default to the topographic DV/NT view
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
                            ),
                            html.Label("Enable Binning:", className="mt-2"),
                            dcc.Checklist(
                                id='enable-binning',
                                options=[{'label': 'Show binned view', 'value': 'enabled'}],
                                value=['enabled'],  # default: binned ("metacell") spatial view
                                className="mb-2"
                            ),
                            html.Div([
                                html.Label("Bin color:", className="mt-2"),
                                dcc.RadioItems(
                                    id='bin-stat',
                                    options=[
                                        {'label': 'Mean expression', 'value': 'mean'},
                                        {'label': '% positive (detected)', 'value': 'frac_pos'},
                                    ],
                                    value='mean',  # 'frac_pos' -> topographic map of fraction of cells expressing
                                    className="mb-2",
                                ),
                                dcc.Checklist(
                                    id='enable-smoothing',
                                    options=[{'label': 'Gaussian smoothing', 'value': 'enabled'}],
                                    value=['enabled'],  # default: smoothed, matching the paper's spatial images
                                    className="mb-2"
                                ),
                                html.Label("Number of bins:", className="mt-2"),
                                dcc.Slider(
                                    id='bin-number-slider',
                                    min=2,
                                    max=120,
                                    step=1,
                                    value=50,
                                    marks={i: str(i) for i in [2, 20, 40, 60, 80, 100, 120]},
                                ),
                                html.Label("Percentile cutoff:", className="mt-2"),
                                dcc.Slider(
                                    id='percentile-slider',
                                    min=0,
                                    max=1,
                                    step=0.01,
                                    value=0.95,
                                    marks={i/10: str(i/10) for i in range(0, 11)},
                                )
                            ], id='binning-controls', style={'display': 'none'}),
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
                                placeholder="Type gene name...",
                                searchable=True,
                                optionHeight=35,
                                persistence=False,  # Ensure fresh search results each time
                                clearable=True
                            ),

                            # Optional second gene for side-by-side comparison.
                            dcc.Checklist(
                                id='compare-genes',
                                options=[{'label': 'Compare a second gene', 'value': 'enabled'}],
                                value=[],  # default: single gene
                                className="mt-2 mb-1",
                            ),
                            html.Div([
                                dcc.Dropdown(
                                    id='gene-select-2',
                                    placeholder="Type second gene name...",
                                    searchable=True,
                                    optionHeight=35,
                                    persistence=False,
                                    clearable=True,
                                ),
                                # Independent per-gene colour scales by default so a
                                # weak gene (CYP26C1) is not washed out by a strong one
                                # (FGF8); tick to compare on one absolute scale.
                                dcc.Checklist(
                                    id='compare-shared-scale',
                                    options=[{'label': ' Shared colour scale (absolute compare)',
                                              'value': 'enabled'}],
                                    value=[],
                                    className="mt-2",
                                ),
                            ], id='gene-select-2-container', style={'display': 'none'}),
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
                    ),
                ], id='map-controls'),

                # ---- GROUP CONTROLS (expression-by-group view) ----
                html.Div([
                    dbc.Card([
                        dbc.CardBody([
                            html.H6("Gene", className="card-subtitle mb-2 text-muted"),
                            dcc.Dropdown(
                                id='group-gene-select',
                                placeholder="Type gene name...",
                                searchable=True,
                                optionHeight=35,
                                persistence=False,
                                clearable=True,
                            ),
                        ])
                    ], className="mt-1", style={'backgroundColor': '#f8f9fa'}),

                    html.Label("Group by (categorical):", className="mt-3"),
                    dcc.Dropdown(
                        id='group-by-select',
                        placeholder="Select grouping column",
                    ),

                    html.Label("Split by (optional):", className="mt-3"),
                    dcc.Dropdown(
                        id='group-split-select',
                        placeholder="(none)",
                    ),

                    html.Label("Style:", className="mt-3"),
                    dcc.RadioItems(
                        id='group-style',
                        options=[
                            {'label': 'Figure: pseudobulk + positive-cell violin', 'value': 'figure'},
                            {'label': 'Violin (all cells)', 'value': 'violin'},
                            {'label': 'Box', 'value': 'box'},
                            {'label': 'Strip', 'value': 'strip'},
                            {'label': 'Dot plot (gene set)', 'value': 'dotplot'},
                        ],
                        value='figure',  # default: the polished 2-panel NPY-style figure
                        labelStyle={'display': 'block'},
                    ),

                    # Figure-style options (shown only for the 'figure' style):
                    # Panel A pseudobulk replicate unit + Panel B positive-cell gate.
                    html.Div([
                        dcc.Checklist(
                            id='group-positive-only',
                            options=[{'label': ' Violin: positive cells only (expr > 0)',
                                      'value': 'enabled'}],
                            value=['enabled'],
                            className="mt-2 mb-1",
                        ),
                        html.Label("Pseudobulk replicate (Panel A):", className="mt-2"),
                        dcc.Dropdown(
                            id='group-replicate-select',
                            placeholder="Replicate unit for pseudobulk dots",
                        ),
                    ], id='figure-style-controls'),

                    # Gene-module selector (shown only for the 'dotplot' style).
                    html.Div([
                        html.Label("Gene module (dot plot):", className="mt-3"),
                        dcc.Dropdown(
                            id='gene-module-select',
                            placeholder="Select a gene module",
                        ),
                    ], id='module-controls', style={'display': 'none'}),
                ], id='group-controls', style={'display': 'none'}),
            ])
        ])
    ], className="mt-3")
