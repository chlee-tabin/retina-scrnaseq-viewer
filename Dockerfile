# Hugging Face Docker Space image for the Retina scRNA-seq Pattern Viewer.
# Build: docker build -t retina-viewer .
# Run:   docker run -p 7860:7860 -v /path/to/h5ads:/app/data retina-viewer
#
# NOTE: this image deliberately does NOT bundle the .h5ad datasets. Provide them
# at runtime via a volume/upload and (optionally) the DATA_DIR env var. See the
# "Data (not bundled)" section of README.md.

FROM python:3.11-slim

# System libs commonly needed by the scientific stack (scanpy/anndata/h5py/scipy).
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Hugging Face Spaces run containers as a non-root user (UID 1000).
RUN useradd -m -u 1000 appuser

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=7860 \
    HOST=0.0.0.0 \
    DATA_DIR=data \
    DASH_DEBUG=false

WORKDIR /app

# Install Python dependencies first to leverage Docker layer caching.
COPY requirements.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt

# Copy application code (data is intentionally excluded; see note above).
COPY . .

# Ensure a default data directory exists and is writable; mount real data here.
RUN mkdir -p /app/data && chown -R appuser:appuser /app

USER appuser

EXPOSE 7860

# app.py reads HOST/PORT/DASH_DEBUG from the environment.
CMD ["python", "app.py"]
