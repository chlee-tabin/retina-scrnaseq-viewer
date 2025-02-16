import scanpy as sc
import logging
import yaml
from pathlib import Path
from functools import lru_cache

logger = logging.getLogger(__name__)

@lru_cache(maxsize=2)
def load_adata(filename):
    logger.info(f"Starting to load {filename}")
    adata = sc.read_h5ad(filename)
    logger.info(f"Successfully loaded {filename} with {adata.n_obs} cells")
    return adata

def load_dataset_config(config_path='datasets_config.yml'):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def validate_datasets(config):
    """Validate that all dataset files exist and are accessible"""
    valid_datasets = {}
    for id, dataset in config['datasets'].items():
        file_path = Path(dataset['file_path'])
        if file_path.exists() and file_path.suffix == '.h5ad':
            valid_datasets[id] = dataset
            logger.info(f"Validated dataset: {dataset['title']}")
        else:
            logger.warning(f"Dataset file not found or invalid: {file_path}")
    return valid_datasets 