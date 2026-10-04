import os
import anndata as ad
import logging
import yaml
import threading
import numpy as np
import scipy.sparse as sp
from utils.validation import obs_column_types
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
    p = str(file_path).replace('\\', '/')
    if os.path.isabs(p):
        return p
    parts = p.split("/", 1)
    if parts[0] == "data" and len(parts) > 1:
        p = parts[1]
    return os.path.join(data_dir, p)


def load_dataset_config(config_path='datasets_config.yml'):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    basenames = {}
    for dataset_id, dataset in config['datasets'].items():
        name = os.path.basename(dataset['file_path'].replace('\\', '/'))
        if name in basenames:
            raise ValueError(f"Duplicate dataset basename '{name}': {basenames[name]} and {dataset_id}")
        basenames[name] = dataset_id
    return config


def _dataset_cache_size():
    return max(1, sum(Path(ds['file_path']).suffix == '.h5ad'
                      for ds in load_dataset_config()['datasets'].values()))


@lru_cache(maxsize=_dataset_cache_size())
def _cached_load_adata(filename):
    # Provision the file on first access (no-op if already on disk or if
    # HF_DATA_REPO is unset). Lets the app start without blocking on downloads;
    # the dataset is fetched the moment it is actually selected.
    from utils.data_provision import ensure_one_dataset
    resolved = ensure_one_dataset(filename)
    logger.info(f"Starting to load {resolved}")
    adata = ad.read_h5ad(resolved)
    logger.info(f"Successfully loaded {resolved} with {adata.n_obs} cells")
    return adata

_load_locks = {}
_load_locks_guard = threading.Lock()


def configured_file_paths():
    return {ds['file_path'] for ds in load_dataset_config()['datasets'].values()}


def store_dataset_id(data_store):
    """Dataset id from the client-held data-store; a forged non-dict store is unknown."""
    return data_store.get('dataset_id') if isinstance(data_store, dict) else None


def get_dataset_config(dataset_id):
    """Resolve identity/settings only from the server's configuration."""
    if not isinstance(dataset_id, str):
        return None
    return load_dataset_config()['datasets'].get(dataset_id)


def load_adata(filename):
    if not isinstance(filename, str) or filename not in configured_file_paths():
        raise ValueError('Unknown dataset file.')
    with _load_locks_guard:
        lock = _load_locks.setdefault(filename, threading.Lock())
    # Lookup the LRU after acquiring the file lock: a concurrent cold miss rechecks
    # the populated cache rather than reading again. Other files load independently.
    with lock:
        return _cached_load_adata(filename)


# Preserve the public cache introspection/invalidation API and per-file load lock.
load_adata.cache_info = _cached_load_adata.cache_info
load_adata.cache_parameters = _cached_load_adata.cache_parameters


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


def gene_search_options(genes, search_value, selected_values=None):
    """Fuzzy starts-with/contains gene search shared by the gene dropdowns.

    With a search string present, starts-with matches sort ahead of contains matches;
    no search results are offered until a character is typed; at most 200 options are
    returned.  A current selection is retained even before a search begins, because Dash
    clears a Dropdown value that is absent from its options.
    Lives here (a leaf util both callback modules already import) so neither callback
    module has to import a private helper from the other.
    """
    if isinstance(selected_values, str):
        selected_values = [selected_values]
    elif not isinstance(selected_values, (list, tuple, set)):
        selected_values = []
    selected = []
    for value in selected_values:
        gene = match_gene(value, genes)
        if gene and gene not in selected:
            selected.append(gene)

    if isinstance(search_value, str) and search_value.strip():
        sv = search_value.strip().lower()
        starts_with = []
        contains = []
        for gene in genes:
            gl = gene.lower()
            if gl.startswith(sv):
                starts_with.append(gene)
            elif sv in gl:
                contains.append(gene)
        matching = sorted(starts_with) + sorted(contains)
    else:
        matching = []
    # Put selected genes first so a selected gene survives the 200-option cap even
    # when it does not match the current search text.
    choices = selected + [gene for gene in matching if gene not in selected]
    return [{'label': gene, 'value': gene} for gene in choices[:200]]


def match_gene(gene, genes):
    """Exact name first, then species-correct casing, before choosing a default."""
    if not isinstance(gene, str):
        return None
    if gene in genes:
        return gene
    return next((g for g in genes if g.casefold() == gene.casefold()), None)

def detect_cp_target(adata, n_sample=256):
    """Detect the counts-per-X normalisation target of a log1p-normalised ``.X``.

    For CP-normalised data every cell's ``expm1(.X)`` sums to the same target (1e4 for
    CP10K, 1e6 for CPM). Sample cells, sum ``expm1`` across genes, and return the
    median when the sums are tightly clustered; return None otherwise (raw counts,
    z-scored, scran-pooled, or otherwise not a clean log1p(CP)) so the pseudobulk panel
    is omitted rather than reconstructed from an unknown normalisation.
    """
    try:
        n = adata.n_obs
        if n == 0:
            return None
        idx = np.unique(np.linspace(0, n - 1, min(n_sample, n)).astype(int))
        X = adata.X[idx]
        if sp.issparse(X):
            sums = np.asarray(np.expm1(X.toarray()).sum(axis=1)).ravel()
        else:
            sums = np.expm1(np.asarray(X)).sum(axis=1).ravel()
        sums = sums[np.isfinite(sums) & (sums > 0)]
        if sums.size < max(5, 0.5 * len(idx)):
            return None
        med = float(np.median(sums))
        if med <= 0:
            return None
        # A genuine single-target CP normalisation makes every cell's sum EQUAL to the
        # target (only ~1e-7 float round-off). A robust 0.5-99.5 percentile band stays
        # ~0 for such data but blows past 2% for any minority block on a different scale
        # (e.g. a CPM arm is 100x off and, even at ~1% of cells, dominates the
        # reconstructed pseudobulk) -- so omit Panel A. Non-finite sums are filtered
        # above; a stray aberrant cell in the 0.5% tails is tolerated.
        lo, hi = np.percentile(sums, [0.5, 99.5])
        if (float(hi) - float(lo)) / med > 0.02:
            return None
        return med
    except Exception as e:  # noqa: BLE001 - best-effort; omit Panel A on any failure
        logger.warning(f"CP-target detection failed ({e}); pseudobulk panel omitted")
        return None

@lru_cache(maxsize=8)
def dataset_column_types(filename):
    """Per-dataset obs column typing, computed once (it scans every obs column)."""
    return obs_column_types(load_adata(filename).obs)


@lru_cache(maxsize=_dataset_cache_size())
def dataset_norm_target(filename):
    """Server-derived normalisation, never a value from dcc.Store."""
    return detect_cp_target(load_adata(filename))


def _clear_dataset_caches():
    """Invalidate the AnnData cache together with everything derived from it."""
    _cached_load_adata.cache_clear()
    dataset_column_types.cache_clear()
    dataset_norm_target.cache_clear()


load_adata.cache_clear = _clear_dataset_caches


def load_dataset_state(data_store):
    """Return (adata, authoritative settings/columns); unknown ids do no file access."""
    dataset_id = store_dataset_id(data_store)
    dataset = get_dataset_config(dataset_id)
    if dataset is None:
        return None, {}
    adata = load_adata(dataset['file_path'])
    settings = dict(dataset)
    settings.update(dataset_id=dataset_id, filename=dataset['file_path'],
                    column_types=dataset_column_types(dataset['file_path']),
                    metadata_cols=list(adata.obs.columns), genes=list(adata.var_names),
                    x_norm_target=dataset_norm_target(dataset['file_path']))
    return adata, settings
