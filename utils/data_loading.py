import os
import anndata as ad
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
    # Provision the file on first access (no-op if already on disk or if
    # HF_DATA_REPO is unset). Lets the app start without blocking on downloads;
    # the dataset is fetched the moment it is actually selected.
    try:
        from utils.data_provision import ensure_one_dataset
        ensure_one_dataset(filename)
    except Exception as e:  # noqa: BLE001 - fall through to the read, which will error clearly
        logger.error(f"On-demand provisioning failed for {filename}: {e}")
    resolved = resolve_data_path(filename)
    logger.info(f"Starting to load {resolved}")
    adata = ad.read_h5ad(resolved)
    logger.info(f"Successfully loaded {resolved} with {adata.n_obs} cells")
    return adata

def load_dataset_config(config_path='datasets_config.yml'):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def validate_datasets(config):
    """List the datasets the app should offer.

    Datasets are provisioned lazily (downloaded on first selection in
    ``load_adata`` -> ``ensure_one_dataset``), so a dataset is listed as long as
    its configured path is an ``.h5ad`` -- we deliberately do NOT require the file
    to already be on disk. This lets the dropdown populate immediately at startup
    instead of being silently emptied while the background download runs (and lets
    the server bind its port without waiting on any download).
    """
    valid_datasets = {}
    for id, dataset in config['datasets'].items():
        file_path = Path(resolve_data_path(dataset['file_path']))
        if file_path.suffix == '.h5ad':
            valid_datasets[id] = dataset
        else:
            logger.warning(f"Dataset '{id}' has a non-.h5ad path, skipping: {file_path}")
    return valid_datasets


def choose_default_embedding(available_embeddings, obs_columns,
                             config_default=None, url_embedding=None):
    """Pick the embedding a dataset should open on.

    Precedence: an embedding named in the shared URL > the dataset's configured
    ``default_embedding`` > the topographic custom view when DV.Score & NT.Score
    both exist (the spatial RPC datasets) > the first UMAP in obsm (so a dataset
    WITHOUT spatial scores, e.g. the full-retina object, opens on its UMAP rather
    than a blank custom-embedding plot) > ``'custom_embedding'`` as a last resort.
    """
    options = list(available_embeddings) + ['custom_embedding']
    if url_embedding in options:
        return url_embedding
    if config_default in options:
        return config_default
    if 'DV.Score' in obs_columns and 'NT.Score' in obs_columns:
        return 'custom_embedding'
    umaps = [e for e in available_embeddings if 'umap' in e.lower()]
    if umaps:
        return umaps[0]
    return available_embeddings[0] if available_embeddings else 'custom_embedding'


def gene_search_options(genes, search_value):
    """Fuzzy starts-with/contains gene search shared by the gene dropdowns.

    With a search string present, starts-with matches sort ahead of contains matches;
    otherwise all genes are returned alphabetically. Returns Dash dropdown option dicts.
    Lives here (a leaf util both callback modules already import) so neither callback
    module has to import a private helper from the other.
    """
    if search_value:
        sv = search_value.lower()
        starts_with = []
        contains = []
        for gene in genes:
            gl = gene.lower()
            if gl.startswith(sv):
                starts_with.append(gene)
            elif sv in gl:
                contains.append(gene)
        matching = sorted(starts_with) + sorted(contains)
        return [{'label': g, 'value': g} for g in matching]
    return [{'label': g, 'value': g} for g in sorted(genes)]