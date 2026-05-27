from dash import Input, Output, State, callback
import dash_bootstrap_components as dbc
from dash import html
from utils.data_loading import load_dataset_config, validate_datasets, load_adata
from utils.error_handling import handle_callback_error, log_callback_info
from utils.config import load_config
import logging
import pandas as pd

logger = logging.getLogger(__name__)

# NOTE: dataset-select options, data-store loading, and dataset-info rendering
# (including the download link) are handled by the authoritative callbacks in
# app.py (initialize_dataset_dropdown, update_data, update_dataset_info). The
# earlier copies here were duplicates that referenced helpers which were never
# implemented (load_data_store / get_embedding_options / get_color_options) and
# built a broken download path; they have been removed. The gene-search
# callback below is unique and remains the sole owner of its outputs.

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