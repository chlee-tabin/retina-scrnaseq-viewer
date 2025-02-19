from dash import Input, Output, State, callback
import dash_bootstrap_components as dbc
from dash import html
from utils.data_loading import load_dataset_config, validate_datasets, load_adata
from utils.error_handling import handle_callback_error, log_callback_info
from utils.config import load_config
import logging
import pandas as pd

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
     Output('color-select', 'value', allow_duplicate=True),
     Output('gene-select', 'value', allow_duplicate=True),
     Output('loading-output', 'children', allow_duplicate=True)],
    [Input('dataset-select', 'value'),
     State('color-select', 'value'),
     State('gene-select', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_data(dataset_id, current_color, current_gene):
    if not dataset_id:
        return None, [], [], None, None, ""
    
    try:
        config = load_dataset_config()
        dataset = config['datasets'][dataset_id]
        adata = load_adata(dataset['file_path'])
        
        # Determine column types
        column_types = {}
        for col in adata.obs.columns:
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                column_types[col] = 'numeric'
            else:
                column_types[col] = 'categorical'
        
        # Create embedding list including custom embedding option
        available_embeddings = list(adata.obsm.keys())
        
        data_store = {
            'filename': dataset['file_path'],
            'n_cells': adata.n_obs,
            'embeddings': available_embeddings + ['custom_embedding'],  # Add custom_embedding to the list
            'metadata_cols': list(adata.obs.columns),
            'genes': list(adata.var_names),
            'column_types': column_types
        }
        
        # Create embedding options with custom embedding as first option
        embedding_options = [
            {'label': 'Custom embedding by meta scores', 'value': 'custom_embedding'}
        ] + [
            {'label': emb, 'value': emb} for emb in available_embeddings
        ]
        
        color_options = [
            {'label': 'Gene Expression', 'value': 'gene_expression'}
        ]
        color_options.extend([
            {'label': col, 'value': col} 
            for col in data_store['metadata_cols']
        ])
        
        valid_color_values = ['gene_expression'] + data_store['metadata_cols']
        color_value = (url_color if url_color in valid_color_values 
                      else current_color if current_color in valid_color_values 
                      else None)
        
        gene_value = (url_gene if url_gene in data_store['genes']
                     else current_gene if current_gene in data_store['genes']
                     else None)
        
        return data_store, embedding_options, color_options, color_value, gene_value, ""
        
    except Exception as e:
        error_message = f"Error loading data: {str(e)}"
        logger.error(error_message)
        return None, [], [], None, None, error_message

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
    
    # Create the download link path
    download_path = f"/download/{dataset['file_path']}"
    
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
            ], className='card-text'),
            # Add download link
            html.P([
                html.Strong("Download: "),
                html.A(
                    "Download .h5ad file",
                    href=download_path,
                    download=dataset['file_path'].split('/')[-1],
                    className="btn btn-outline-primary btn-sm"
                )
            ], className='card-text')
        ])
    ])

@callback(
    [Output('gene-select-container', 'style', allow_duplicate=True),
     Output('gene-select', 'options', allow_duplicate=True)],
    [Input('color-select', 'value'),
     Input('data-store', 'data'),
     Input('gene-select', 'search_value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_gene_select(color_value, data_store, search_value):
    logger.debug("update_gene_select called with:")
    logger.debug(f"color_value: {color_value}")
    logger.debug(f"search_value: {search_value}")

    if color_value != 'gene_expression' or not data_store or 'genes' not in data_store:
        logger.debug("Hiding gene-select container")
        return {'display': 'none'}, []

    genes = data_store['genes']
    
    if search_value:
        search_value = search_value.lower()
        # First, find all matching genes
        matching_genes = []
        starts_with = []
        contains = []
        
        for gene in genes:
            gene_lower = gene.lower()
            if gene_lower.startswith(search_value):
                starts_with.append(gene)
            elif search_value in gene_lower:
                contains.append(gene)
        
        # Sort and combine matches
        matching_genes = sorted(starts_with) + sorted(contains)
        
        # Create options only for matching genes
        gene_options = [{'label': gene, 'value': gene} for gene in matching_genes]
    else:
        # When no search, show all genes in alphabetical order
        gene_options = [{'label': gene, 'value': gene} for gene in sorted(genes)]
    
    logger.debug(f"Generated {len(gene_options)} gene options")
    if search_value:
        logger.debug(f"First few matches for '{search_value}': {[opt['label'] for opt in gene_options[:5]]}")
    
    return {'display': 'block'}, gene_options 