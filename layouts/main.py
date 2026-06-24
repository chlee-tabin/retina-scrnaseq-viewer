from dash import html, dcc
import dash_bootstrap_components as dbc

def create_main_panel():
    return html.Div([
        dcc.Loading(
            id="loading-upload",
            type="circle",
            children=[
                html.Div(id="loading-output"),
                # Holds the user's last 3D-sphere camera so it persists across gene/colour
                # re-renders (dcc.Graph's uirevision doesn't reliably hold the 3D camera).
                dcc.Store(id='sphere-camera-store'),
                dcc.Graph(
                    id='main-plot',
                    style={'height': '85vh'}  # Make the plot larger
                ),
                html.Div(id='selection-info', className='mt-3')
            ]
        )
    ], className="p-3") 