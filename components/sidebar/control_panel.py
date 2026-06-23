from dash import html, dcc
import dash_bootstrap_components as dbc


def _help(help_id, header, text):
    """A clickable '?' badge that opens a Popover with a longer explanation on CLICK
    (dbc trigger='legacy': click to open, click away to dismiss). `help_id` must be unique.
    Click-to-open (not hover) avoids the help-cursor glyph and lets the text be read."""
    return html.Span([
        dbc.Badge("?", id=help_id, color="info", pill=True,
                  className="ms-1", style={'cursor': 'pointer'}),
        dbc.Popover(
            [dbc.PopoverHeader(header), dbc.PopoverBody(text)],
            target=help_id, trigger="legacy", placement="right",
            style={'maxWidth': '340px'}),
    ])


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
                                    {'label': 'Whole-mount sphere (3D)', 'value': 'sphere'},
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
                                # Geometry knobs (inherit DEFAULT_PARAMS; every default below
                                # reproduces Fig R2.5). The controls in THIS group shape both the
                                # flat flower AND the 3D sphere (they set the spherical-cap angle,
                                # pole and scaling); the "flat-layout relief" group further down is
                                # flower-only and is hidden under the sphere projection.
                                html.Div("Cap extent, pole & scaling — affect flower and sphere.",
                                         className="text-muted fst-italic", style={'fontSize': '11px'}),
                                html.Label([
                                    "Nasotemporal extent (deg):",
                                    _help('help-rho-nt', "Nasotemporal extent",
                                          "Max colatitude (angle from the HAA pole) the most extreme "
                                          "nasal/temporal cells reach on the cap. Default 82° reproduces "
                                          "Fig R2.5. 90° = the equator (a hemisphere); a real retina lines "
                                          "~72–76% of the eye = a cap of ~116–121°, so the slider runs to "
                                          "150°. NOTE: with the arcsin de-warp the effect saturates near "
                                          "90° and folds past it (110° → a SMALLER cap than 90°) — to use "
                                          "the >90° range switch de-warp to power or none."),
                                ], className="mt-1 small"),
                                dcc.Slider(id='wm-rho-nt', min=40, max=150, step=1, value=82,
                                           marks={40: '40', 90: '90', 120: '120', 150: '150'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                html.Label([
                                    "Dorsoventral extent (deg):",
                                    _help('help-rho-dv', "Dorsoventral extent",
                                          "Same as nasotemporal extent but for the dorsal/ventral axis. "
                                          "Default 64° — smaller than the 82° N-T, which is why the map is "
                                          "elongated nasotemporally. 90° = equator; ~116–121° ≈ a real "
                                          "retina's coverage. Same arcsin saturation/fold caveat near and "
                                          "past 90°; use power/none de-warp to push past the equator."),
                                ], className="mt-2 small"),
                                dcc.Slider(id='wm-rho-dv', min=40, max=150, step=1, value=64,
                                           marks={40: '40', 90: '90', 120: '120', 150: '150'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),
                                dcc.Checklist(
                                    id='wm-symmetric',
                                    options=[{'label': ' Symmetric per-side scaling (equalize opposing petals)',
                                              'value': 'enabled'}],
                                    value=['enabled'],   # R2.5; untick -> p99 keeps true NT/DV asymmetry
                                    className="small mb-1 mt-2",
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
                                html.Label([
                                    "Radial de-warp:",
                                    _help('help-dewarp', "Radial de-warp",
                                          "How a cell's score-radius r (0 at the pole, 1 at the rim) becomes "
                                          "colatitude ρ. arcsin: ρ=arcsin(r·sin(extent)) — treats r as the "
                                          "flattened (orthographic) image of a sphere and inverts it; the "
                                          "spherical-cap model, but its response to 'extent' saturates near "
                                          "90° and cannot exceed it. power: ρ=extent·r^p (see Power exponent). "
                                          "none: ρ=extent·r (linear). Use power or none to reach >90°."),
                                ], className="mt-1 small"),
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
                                html.Label([
                                    "Power exponent (de-warp = power):",
                                    _help('help-pow', "Power exponent",
                                          "Only used when de-warp = power. ρ = extent · r^p. p=1 is linear; "
                                          "p>1 pushes cells toward the pole (compresses the centre, expands "
                                          "the rim); p<1 the opposite. Default 1.6. Unlike arcsin, 'extent' "
                                          "scales ρ directly here, so raising it past 90° genuinely enlarges "
                                          "the cap toward the realistic retinal extent (~116–121°)."),
                                ], className="mt-1 small"),
                                dcc.Slider(id='wm-pow', min=1.0, max=3.0, step=0.1, value=1.6,
                                           marks={1: '1', 1.6: '1.6', 2: '2', 3: '3'},
                                           tooltip={'placement': 'bottom', 'always_visible': False}),

                                # Flat-layout relief: the orange-peel cuts / wedge gaps / petal
                                # stretch that turn the flat disk into a dissected whole-mount.
                                # These are purely a FLATTENING device -- a sphere has no curvature
                                # deficit to relieve -- so the sphere projection hides this group
                                # (toggle_wholemount_relief). Defaults reproduce Fig R2.5.
                                html.Div([
                                    html.Hr(className="my-2"),
                                    html.Label("Flat-layout relief (flower only)",
                                               className="fw-semibold small"),
                                    html.Label("Number of cuts (petals):", className="mt-1 small"),
                                    dcc.Slider(id='wm-cuts', min=0, max=8, step=1, value=4,
                                               marks={0: '0', 2: '2', 4: '4', 6: '6', 8: '8'},
                                               tooltip={'placement': 'bottom', 'always_visible': False}),
                                    html.Label("Relief gap (deficit mode):", className="mt-2 small"),
                                    dcc.Slider(id='wm-gap', min=0, max=1.5, step=0.05, value=1.0,
                                               marks={0: '0', 1: '1', 1.5: '1.5'},
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
                                    html.Label("Petal stretch:", className="mt-2 small"),
                                    dcc.Slider(id='wm-stretch', min=0, max=1.5, step=0.05, value=1.0,
                                               marks={0: '0', 1: '1', 1.5: '1.5'},
                                               tooltip={'placement': 'bottom', 'always_visible': False}),
                                ], id='wholemount-relief-controls'),
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
