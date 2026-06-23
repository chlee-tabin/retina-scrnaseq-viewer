from dash import html, dcc
import dash_bootstrap_components as dbc

def create_control_panel():
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
                        {'label': 'ROI differential expression (volcano)', 'value': 'deg'},
                    ],
                    value='map',
                    className="mb-3",
                    labelStyle={'display': 'block'},
                ),

                # ---- ROI-DEG CONTROLS (shown only for the volcano view) ----
                html.Div([
                    html.Label("Min cells per pseudobulk:"),
                    dcc.Input(id='deg-min-cells', type='number', value=50, min=5, step=5,
                              debounce=True, style={'width': '100%'}),
                    html.Div(
                        "Pseudobulks below this cell count are dropped (manuscript floor "
                        "= 50). Lower it for small ROIs, at the cost of noisier estimates.",
                        className="small text-muted mt-1"),
                ], id='deg-controls', style={'display': 'none'}, className="mb-3"),

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

                            # Projection of the topographic scores: the raw X-vs-Y plane,
                            # or the whole-mount ("flower") reprojection that warps the flat
                            # DV/NT plane into a flat-mounted retina (reviewer-response
                            # transform; needs DV.Score + NT.Score).
                            html.Label("Projection:", className="mt-2"),
                            dcc.RadioItems(
                                id='custom-projection',
                                options=[
                                    {'label': 'Raw axes (X vs Y)', 'value': 'raw'},
                                    {'label': 'Whole-mount (DV/NT flower)', 'value': 'flower'},
                                ],
                                value='raw',
                                className="mb-2",
                                labelStyle={'display': 'block'},
                            ),

                            # Raw-axis pickers, hidden under the whole-mount projection
                            # (which derives its coordinates from DV.Score / NT.Score).
                            html.Div([
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
                            ], id='raw-axes-controls'),

                            # Advanced whole-mount projection parameters (shown only under
                            # the flower projection). Defaults reproduce the shipped
                            # reviewer figure; every knob feeds wholemount.flower_transform.
                            html.Div([
                                html.Hr(className="my-2"),
                                html.Label("Advanced projection", className="fw-semibold small"),
                                html.Label("Nasotemporal extent (deg):", className="mt-1 small"),
                                dcc.Slider(id='wm-rho-nt', min=40, max=110, step=1, value=82,
                                           marks={40: '40', 82: '82', 110: '110'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label("Dorsoventral extent (deg):", className="mt-2 small"),
                                dcc.Slider(id='wm-rho-dv', min=40, max=110, step=1, value=64,
                                           marks={40: '40', 64: '64', 110: '110'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label("Relief gap (deficit mode):", className="mt-2 small"),
                                dcc.Slider(id='wm-gap', min=0, max=1.5, step=0.05, value=1.0,
                                           marks={0: '0', 1: '1', 1.5: '1.5'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label("Petal stretch:", className="mt-2 small"),
                                dcc.Slider(id='wm-stretch', min=0, max=1.5, step=0.05, value=1.0,
                                           marks={0: '0', 1: '1', 1.5: '1.5'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label("Number of cuts (petals):", className="mt-2 small"),
                                dcc.Slider(id='wm-cuts', min=0, max=8, step=1, value=4,
                                           marks={0: '0', 2: '2', 4: '4', 6: '6', 8: '8'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),

                                # Scaling / geometry knobs (inherit DEFAULT_PARAMS; every
                                # default below reproduces Fig R2.5). Each is always shown so
                                # the layout effect is explorable; labels note when one only
                                # applies in a particular mode.
                                html.Hr(className="my-2"),
                                dcc.Checklist(
                                    id='wm-symmetric',
                                    options=[{'label': ' Symmetric per-side scaling (equalize opposing petals)',
                                              'value': 'enabled'}],
                                    value=['enabled'],   # R2.5; untick -> p99 keeps true NT/DV asymmetry
                                    className="small mb-1",
                                ),
                                html.Label("Pole (score origin):", className="mt-1 small"),
                                dcc.RadioItems(
                                    id='wm-pole',
                                    options=[
                                        {'label': 'Origin (0,0) — HAA pole', 'value': 'origin'},
                                        {'label': 'Median of cells', 'value': 'median'},
                                    ],
                                    value='origin',
                                    className="small mb-1",
                                    labelStyle={'display': 'block'},
                                ),
                                html.Label("Radial de-warp:", className="mt-1 small"),
                                dcc.RadioItems(
                                    id='wm-dewarp',
                                    options=[
                                        {'label': 'arcsin (spherical-cap)', 'value': 'arcsin'},
                                        {'label': 'power', 'value': 'pow'},
                                        {'label': 'none (linear)', 'value': 'none'},
                                    ],
                                    value='arcsin',
                                    className="small mb-1",
                                    labelStyle={'display': 'block'},
                                ),
                                html.Label("Power exponent (de-warp = power):", className="mt-1 small"),
                                dcc.Slider(id='wm-pow', min=1.0, max=3.0, step=0.1, value=1.6,
                                           marks={1: '1', 1.6: '1.6', 2: '2', 3: '3'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label("Relief mode:", className="mt-2 small"),
                                dcc.RadioItems(
                                    id='wm-gap-mode',
                                    options=[
                                        {'label': 'deficit (curvature-true rifts)', 'value': 'deficit'},
                                        {'label': 'linear (V-notch)', 'value': 'linear'},
                                    ],
                                    value='deficit',
                                    className="small mb-1",
                                    labelStyle={'display': 'block'},
                                ),
                                html.Label("Linear rip width (relief = linear):", className="mt-1 small"),
                                dcc.Slider(id='wm-gap-frac', min=0.0, max=1.0, step=0.05, value=0.5,
                                           marks={0: '0', 0.5: '0.5', 1: '1'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                            ], id='wholemount-advanced', style={'display': 'none'},
                               className="mb-2"),

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
                    # Plot value: a gene's expression, or any continuous .obs variable
                    # (QC metric, DV/NT score, ...).
                    html.Label("Plot value:"),
                    dcc.RadioItems(
                        id='group-value-source',
                        options=[
                            {'label': 'Gene expression', 'value': 'gene'},
                            {'label': 'Metadata (continuous)', 'value': 'meta'},
                        ],
                        value='gene',
                        className="mb-2",
                        labelStyle={'display': 'block'},
                    ),

                    # Gene picker (shown for the 'gene' value source).
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
                    ], id='group-gene-block'),

                    # Continuous-metadata picker (shown for the 'meta' value source).
                    html.Div([
                        html.Label("Metadata variable (continuous):", className="mt-1"),
                        dcc.Dropdown(
                            id='group-meta-select',
                            placeholder="Select a continuous .obs variable",
                        ),
                    ], id='group-meta-block', style={'display': 'none'}),

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
