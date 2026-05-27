"""Optional startup provisioning of dataset files from a Hugging Face repo.

The ``.h5ad`` datasets are intentionally not committed to this repo or baked
into the Docker image (see ``.dockerignore`` / ``Dockerfile``). On Hugging Face
Spaces we instead keep them in a separate Hugging Face repo and download the
files referenced by ``datasets_config.yml`` into ``DATA_DIR`` at startup.

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

logger = logging.getLogger(__name__)


def ensure_datasets(config):
    """Download any configured dataset files that are missing from ``DATA_DIR``.

    Parameters
    ----------
    config : dict
        Parsed ``datasets_config.yml`` (expects a ``datasets`` mapping whose
        entries have a ``file_path``).
    """
    repo = os.getenv("HF_DATA_REPO")
    if not repo:
        logger.info(
            "HF_DATA_REPO not set; skipping Hugging Face data provisioning "
            "(expecting data already present in DATA_DIR)."
        )
        return

    data_dir = os.getenv("DATA_DIR", "data")
    repo_type = os.getenv("HF_DATA_REPO_TYPE", "dataset")
    token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")

    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        logger.error(
            "HF_DATA_REPO is set but huggingface_hub is not installed; "
            "cannot provision data. Add huggingface_hub to requirements.txt."
        )
        return

    os.makedirs(data_dir, exist_ok=True)

    # Collect the unique basenames referenced by the config.
    datasets = (config or {}).get("datasets", {}) or {}
    wanted = []
    for ds in datasets.values():
        name = os.path.basename(str(ds.get("file_path", "")))
        if name and name not in wanted:
            wanted.append(name)

    for name in wanted:
        dest = os.path.join(data_dir, name)
        if os.path.exists(dest):
            logger.info(f"Dataset already present, skipping download: {dest}")
            continue
        logger.info(f"Downloading '{name}' from {repo_type} repo '{repo}' ...")
        try:
            hf_hub_download(
                repo_id=repo,
                filename=name,
                repo_type=repo_type,
                local_dir=data_dir,
                token=token,
            )
            logger.info(f"Downloaded '{name}' -> {dest}")
        except Exception as e:  # noqa: BLE001 - log and continue; missing data is filtered later
            logger.error(f"Failed to download '{name}' from '{repo}': {e}")
