from dash import html
import dash_bootstrap_components as dbc
from utils.error_handling import handle_callback_error

def create_dataset_info():
    return html.Div([
        dbc.Card([
            dbc.CardBody([
                html.Div(id='dataset-info'),
                # Add statistics panel
                html.Div(id='dataset-stats', className="mt-3")
            ])
        ], className="mt-3")
    ])

@handle_callback_error
def format_dataset_info(dataset, n_cells=None, stats=None):
    """Helper function to format dataset information"""
    info_elements = [
        html.H5(dataset['title'], className='card-title'),
        html.P(dataset['description'], className='card-text'),
        html.P([
            html.Strong("Last Updated: "),
            dataset['last_updated']
        ], className='card-text'),
    ]
    
    if n_cells:
        info_elements.append(html.P([
            html.Strong("Number of Cells: "),
            f"{n_cells:,}"
        ], className='card-text'))
    
    if stats:
        info_elements.extend([
            html.Hr(),
            html.H6("Dataset Statistics:", className="mt-3"),
            html.P([
                html.Strong("Mean Genes per Cell: "),
                f"{stats.get('mean_genes', 'N/A'):,.1f}"
            ], className='card-text'),
            html.P([
                html.Strong("Mean Counts per Cell: "),
                f"{stats.get('mean_counts', 'N/A'):,.1f}"
            ], className='card-text')
        ])
    
    return info_elements 