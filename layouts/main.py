from dash import html, dcc
import dash_bootstrap_components as dbc

def create_main_panel():
    return html.Div([
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
    ], className="p-3") 