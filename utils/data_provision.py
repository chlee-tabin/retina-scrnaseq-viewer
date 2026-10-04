"""Provisioning of dataset files from a Hugging Face repo.

The ``.h5ad`` datasets are intentionally not committed to this repo or baked
into the Docker image (see ``.dockerignore`` / ``Dockerfile``). On Hugging Face
Spaces we keep them in a separate Hugging Face repo and download the files
referenced by ``datasets_config.yml`` into ``DATA_DIR``.

Provisioning is **lazy + non-blocking**: files are downloaded on first access
(``ensure_one_dataset`` from ``load_adata``) and, optionally, pre-warmed in a
background thread (``ensure_datasets``). Crucially this is NOT done synchronously
at import time, so the web server can bind its port and pass the platform health
check immediately instead of blocking ~30-40 s on the full download (which was
preventing the Space from being promoted after a restart).

This is a **no-op unless ``HF_DATA_REPO`` is set**, so local development (where
the files already live in ``DATA_DIR``) is unaffected.

Environment variables
----------------------
``HF_DATA_REPO``       e.g. ``"username/retina-scrnaseq-data"`` (required to activate)
``HF_DATA_REPO_TYPE``  ``"dataset"`` (default), ``"model"`` or ``"space"``
``DATA_DIR``           target directory (default ``"data"``)
``HF_TOKEN``           access token; only needed while the HF repo is private
"""
import os
import logging
import threading

logger = logging.getLogger(__name__)

# One lock per file basename so concurrent callers (the background pre-warm thread
# and an on-demand load) never download the same file at the same time.
_file_locks = {}
_locks_guard = threading.Lock()


def _lock_for(name):
    with _locks_guard:
        if name not in _file_locks:
            _file_locks[name] = threading.Lock()
        return _file_locks[name]


def ensure_one_dataset(file_path):
    """Ensure a single dataset file exists in ``DATA_DIR``, downloading it from
    ``HF_DATA_REPO`` on demand if missing. Returns the local path.

    Safe to call repeatedly and concurrently: a no-op when the file is already
    present, and serialized per-file so two callers can't race the same download.
    If ``HF_DATA_REPO`` is unset the local path is returned as-is (local dev).
    """
    from utils.data_loading import configured_file_paths
    if not isinstance(file_path, str) or file_path not in configured_file_paths():
        raise ValueError('Unknown dataset file.')
    data_dir = os.getenv("DATA_DIR", "data")
    name = os.path.basename(str(file_path))
    dest = os.path.join(data_dir, name)
    if not name:
        return dest
    if os.path.exists(dest):
        return dest

    repo = os.getenv("HF_DATA_REPO")
    if not repo:
        # Local development: nothing to download; caller handles a missing file.
        return dest

    with _lock_for(name):
        if os.path.exists(dest):  # re-check now that we hold the lock
            return dest
        try:
            from huggingface_hub import hf_hub_download
        except ImportError:
            logger.error(
                "HF_DATA_REPO is set but huggingface_hub is not installed; "
                "cannot provision %s.", name
            )
            return dest
        repo_type = os.getenv("HF_DATA_REPO_TYPE", "dataset")
        token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")
        os.makedirs(data_dir, exist_ok=True)
        logger.info("Downloading '%s' from %s repo '%s' ...", name, repo_type, repo)
        try:
            hf_hub_download(
                repo_id=repo,
                filename=name,
                repo_type=repo_type,
                local_dir=data_dir,
                token=token,
            )
            logger.info("Downloaded '%s' -> %s", name, dest)
        except Exception as e:  # noqa: BLE001 - log and continue; missing data filtered later
            logger.error("Failed to download '%s' from '%s': %s", name, repo, e)
    return dest


def ensure_datasets(config):
    """Best-effort pre-warm of every configured dataset file that is missing.

    Intended to be run in a BACKGROUND thread (see app.py) so it never blocks the
    web server from binding its port. Each file is also fetched on demand by
    ``load_adata`` -> ``ensure_one_dataset`` if a user selects it before the
    pre-warm reaches it, so nothing breaks if this is still in progress.
    """
    repo = os.getenv("HF_DATA_REPO")
    if not repo:
        logger.info(
            "HF_DATA_REPO not set; skipping Hugging Face data provisioning "
            "(expecting data already present in DATA_DIR)."
        )
        return

    datasets = (config or {}).get("datasets", {}) or {}
    seen = set()
    for ds in datasets.values():
        name = os.path.basename(str(ds.get("file_path", "")))
        if name and name not in seen:
            seen.add(name)
            ensure_one_dataset(ds['file_path'])
