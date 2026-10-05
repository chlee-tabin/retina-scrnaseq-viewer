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

Downloads require ``HF_DATA_REPO``; existing local files are also verified against
their configured SHA256, with verification stamps cached for fast restarts.

Environment variables
----------------------
``HF_DATA_REPO``       e.g. ``"username/retina-scrnaseq-data"`` (required to activate)
``HF_DATA_REPO_TYPE``  ``"dataset"`` (default), ``"model"`` or ``"space"``
``HF_DATA_REVISION``   overrides the configured immutable data revision
``DATA_DIR``           target directory (default ``"data"``)
``HF_TOKEN``           access token; only needed while the HF repo is private
"""
import os
import logging
import threading
import hashlib
import json
from pathlib import Path

logger = logging.getLogger(__name__)

# One lock per file basename so concurrent callers (the background pre-warm thread
# and an on-demand load) never download the same file at the same time.
_file_locks = {}
_locks_guard = threading.Lock()
_verified_files = {}


def verify_sha256(path, expected):
    """Stream once per (path, size, mtime); persist the stamp for fast restarts."""
    path = Path(path)
    stamp = path.stat()
    key = (str(path.resolve()), stamp.st_size, stamp.st_mtime_ns)
    sidecar = Path(str(path) + '.sha256.json')
    saved = _verified_files.get(key)
    if saved is None:
        try:
            record = json.loads(sidecar.read_text())
            if isinstance(record, dict) and record.get('stamp') == list(key):
                saved = record.get('sha256')
        except (OSError, ValueError):
            pass
    if saved != expected:
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        saved = digest.hexdigest()
        if saved != expected:
            path.unlink()
            sidecar.unlink(missing_ok=True)
            logger.error('SHA256 mismatch for %s: expected %s, got %s', path.name, expected, saved)
            raise ValueError(f"SHA256 mismatch for dataset {path.name}; bad file deleted.")
        try:
            sidecar.write_text(json.dumps({'stamp': list(key), 'sha256': saved}))
        except OSError:
            logger.debug('Could not persist checksum stamp for %s', path.name)
    _verified_files[key] = saved
    return True


def _lock_for(name):
    with _locks_guard:
        if name not in _file_locks:
            _file_locks[name] = threading.Lock()
        return _file_locks[name]


def ensure_one_dataset(file_path):
    """Ensure a single dataset file exists in ``DATA_DIR``, downloading it from
    ``HF_DATA_REPO`` on demand if missing. Returns the local path.

    Safe to call repeatedly and concurrently: verification is cached when the file
    is unchanged, and serialized per-file so callers can't race a download.
    If ``HF_DATA_REPO`` is unset a missing local path is returned to the reader.
    """
    from utils.data_loading import load_dataset_config, resolve_data_path
    config = load_dataset_config()
    dataset = next((ds for ds in config['datasets'].values()
                    if ds['file_path'] == file_path), None) if isinstance(file_path, str) else None
    if dataset is None:
        raise ValueError('Unknown dataset file.')
    data_dir = os.getenv("DATA_DIR", "data")
    name = os.path.basename(str(file_path))
    dest = resolve_data_path(file_path)
    repo = os.getenv("HF_DATA_REPO")
    with _lock_for(dest):
        if not repo:
            # Local development: files in DATA_DIR are the developer's own (possibly
            # synthetic) copies -- never checksum-verify or delete them.
            return dest
        if os.path.exists(dest):
            if not dataset.get('sha256'):
                return dest
            try:
                verify_sha256(dest, dataset['sha256'])
                return dest
            except ValueError:
                # A provisioned copy that no longer matches the pinned data was deleted
                # by verify_sha256; fetch the pinned revision again below.
                logger.warning("Re-downloading '%s' after checksum mismatch", name)
        from huggingface_hub import hf_hub_download
        repo_type = os.getenv("HF_DATA_REPO_TYPE", "dataset")
        token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")
        os.makedirs(data_dir, exist_ok=True)
        relative = Path(dest).resolve().relative_to(Path(data_dir).resolve()).as_posix()
        logger.info("Downloading '%s' from %s repo '%s' ...", name, repo_type, repo)
        hf_hub_download(
            repo_id=repo,
            filename=relative,
            repo_type=repo_type,
            revision=os.getenv('HF_DATA_REVISION') or config.get('hf_data_revision'),
            local_dir=data_dir,
            token=token,
        )
        if dataset.get('sha256'):
            verify_sha256(dest, dataset['sha256'])
        logger.info("Downloaded and verified '%s'", name)
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
            try:
                ensure_one_dataset(ds['file_path'])
            except Exception:
                logger.exception("Failed to provision '%s'", name)
