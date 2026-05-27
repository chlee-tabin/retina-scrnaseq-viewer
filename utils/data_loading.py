import os
import scanpy as sc
import logging
import yaml
from pathlib import Path
from functools import lru_cache

logger = logging.getLogger(__name__)


def resolve_data_path(file_path):
    """Resolve a dataset path against the configurable DATA_DIR.

    Dataset paths in datasets_config.yml are written relative to a ``data/``
    directory (e.g. ``data/20240815_full.h5ad``). DATA_DIR (default ``data``)
    lets the data live elsewhere -- mounted or uploaded separately rather than
    baked into the image -- without editing the config. A leading ``data/`` (or
    ``data\\``) segment is stripped before rejoining under DATA_DIR; absolute
    paths are honored as-is.
    """
    data_dir = os.getenv("DATA_DIR", "data")
    p = str(file_path)
    if os.path.isabs(p):
        return p
    parts = p.replace("\\", "/").split("/", 1)
    if parts[0] == "data" and len(parts) > 1:
        p = parts[1]
    return os.path.join(data_dir, p)


@lru_cache(maxsize=2)
def load_adata(filename):
    resolved = resolve_data_path(filename)
    logger.info(f"Starting to load {resolved}")
    adata = sc.read_h5ad(resolved)
    logger.info(f"Successfully loaded {resolved} with {adata.n_obs} cells")
    return adata

def load_dataset_config(config_path='datasets_config.yml'):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def validate_datasets(config):
    """Validate that all dataset files exist and are accessible"""
    valid_datasets = {}
    for id, dataset in config['datasets'].items():
        file_path = Path(resolve_data_path(dataset['file_path']))
        if file_path.exists() and file_path.suffix == '.h5ad':
            valid_datasets[id] = dataset
            logger.info(f"Validated dataset: {dataset['title']}")
        else:
            logger.warning(f"Dataset file not found or invalid: {file_path}")
    return valid_datasets 