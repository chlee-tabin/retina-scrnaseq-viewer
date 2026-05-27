from dash import html, dcc
import dash_bootstrap_components as dbc
import logging
from utils.error_handling import handle_callback_error
from utils.data_loading import load_dataset_config, validate_datasets

logger = logging.getLogger(__name__)

def create_dataset_info():
    # Populate the dataset dropdown at layout-build time and pre-select the
    # first dataset, so the app opens on a dataset instead of a blank state.
    try:
        datasets = validate_datasets(load_dataset_config())
        dataset_options = [{'label': d['title'], 'value': i} for i, d in datasets.items()]
        default_dataset = next(iter(datasets), None)
    except Exception as e:
        logger.error(f"Could not pre-load dataset options: {e}")
        dataset_options, default_dataset = [], None

    return html.Div([
        # Dataset selector dropdown
        html.Label("Dataset:"),
        dcc.Dropdown(
            id='dataset-select',
            options=dataset_options,
            value=default_dataset,
            placeholder="Select dataset",
            className="mb-3"
        ),
        
        # Dataset info card
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