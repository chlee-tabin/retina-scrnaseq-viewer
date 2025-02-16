from dash import Input, Output, State, callback
import dash_bootstrap_components as dbc
from dash import html
from utils.data_loading import load_dataset_config, validate_datasets, load_adata
from utils.error_handling import handle_callback_error, log_callback_info
from utils.config import load_config
import logging

logger = logging.getLogger(__name__)

@callback(
    Output('dataset-select', 'options', allow_duplicate=True),
    Input('url', 'pathname'),
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def initialize_dataset_dropdown(pathname):
    config = load_dataset_config()
    datasets = validate_datasets(config)
    
    return [
        {'label': data['title'], 'value': dataset_id}
        for dataset_id, data in datasets.items()
    ]

@callback(
    [Output('data-store', 'data', allow_duplicate=True),
     Output('embedding-select', 'options', allow_duplicate=True),
     Output('color-select', 'options', allow_duplicate=True),
     Output('loading-output', 'children', allow_duplicate=True)],
    Input('dataset-select', 'value'),
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_data(dataset_id):
    if not dataset_id:
        return None, [], [], ""
    
    config = load_dataset_config()
    dataset = config['datasets'][dataset_id]
    adata = load_adata(dataset['file_path'])
    
    data_store = {
        'filename': dataset['file_path'],
        'n_cells': adata.n_obs,
        'embeddings': list(adata.obsm.keys()),
        'metadata_cols': list(adata.obs.columns)
    }
    
    embedding_options = [{'label': emb, 'value': emb} for emb in data_store['embeddings']]
    color_options = [{'label': col, 'value': col} for col in data_store['metadata_cols']]
    
    return data_store, embedding_options, color_options, ""

@callback(
    Output('dataset-info', 'children', allow_duplicate=True),
    [Input('dataset-select', 'value'),
     Input('data-store', 'data')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_dataset_info(dataset_id, data_store):
    if not dataset_id:
        return ""
    
    config = load_dataset_config()
    dataset = config['datasets'][dataset_id]
    
    return dbc.Card([
        dbc.CardBody([
            html.H5(dataset['title'], className='card-title'),
            html.P(dataset['description'], className='card-text'),
            html.P([
                html.Strong("Last Updated: "),
                dataset['last_updated']
            ], className='card-text'),
            html.P([
                html.Strong("Number of Cells: "),
                f"{data_store['n_cells']:,}" if data_store else "Loading..."
            ], className='card-text')
        ])
    ]) 