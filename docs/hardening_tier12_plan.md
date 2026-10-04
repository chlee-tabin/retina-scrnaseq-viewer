# Public-release hardening — Tier 1 (abuse/security) + Tier 2 (wrong results)

Spec for one PR. Source: a whole-repo review (Codex + multi-agent) before the repo is made
public. The app is a Plotly Dash viewer served publicly (Hugging Face Space, gunicorn
1 worker × 4 threads, 16 GB). Every Dash callback input can be forged by a client, so
**slider/Input min/max are not a security boundary** — validate on the server.

Keep the diff minimal and in the style of the surrounding code. Reuse existing helpers
(`utils/smoothing.py::_coerce_sigma` is the model for clamping). No new dependencies.

## Tier 1 — abuse / security

### T1.1 Server-side bounds for every numeric/structured client input
Add small coerce helpers (one module, e.g. extend `utils/state.py` or a new
`utils/validation.py`) and apply them at BOTH entry points: share-link restore
(`callbacks/url_callbacks.py::initialize_from_url`) and the compute callbacks that consume
the values (`callbacks/main_callbacks.py::update_plot` and any other callback that passes
these to numpy/scipy; `callbacks/deg_callbacks.py`). Non-numeric / non-finite → fallback.
Ranges = the UI's own ranges:
- `bins` / bin_number: int in [2, 120] (UI default 50)
- `percentile`: float in [0, 1] (default 0.95)
- smoothing sigma: reuse `_coerce_sigma` ([0, 4])
- `wm_rho_nt`, `wm_rho_dv`: [40, 150]; `wm_cuts`: int [0, 8]; `wm_gap`: [0, 1.5];
  `wm_gap_frac`: [0, 1]; `wm_stretch`: [0, 1.5]; check `components/sidebar/control_panel.py`
  for every other slider/Input restored from the URL (`wm_*`, `ng_*`, volcano thresholds,
  sphere camera) and clamp each to its declared min/max.
- `deg_min_cells`: int, clamp to [5, 100000]; `0`/missing → 50 (do not silently turn a
  valid user value into 50; only invalid ones).
- ROI polygons (URL `roi` and the ROI store): dict with only keys A/B; each a list of
  [x, y] finite-number pairs; at most 500 vertices per polygon; anything else → treated
  as no ROI (with no exception).
- Categorical colour / group-by column: refuse (fall back to default) columns with more
  than 200 categories; apply the same cardinality cap where `app.py` and
  `callbacks/main_callbacks.py` decide which obs columns are offered as categorical, so
  ID-like columns (`barcode`, `cell_id`) are never offered and cannot be forced via URL.

### T1.2 Dataset file identity is server-side only
Today `data-store['filename']` (client-held `dcc.Store`) is passed to `load_adata`, to HF
provisioning, and to the download link. Change so the server only ever loads/provisions/
serves files declared in `datasets_config.yml`:
- Callbacks resolve the file path (and `replicate_columns`, `x_norm_target`-relevant
  settings, `figure_regions`, etc.) from the dataset id via `load_dataset_config()` on the
  server; never from a client-supplied filename/path. Unknown dataset id → no-op/error
  message, no file access.
- `utils/data_provision.py`: only download configured `file_path`s.
- `/download/<path>` route in `app.py`: only serve basenames of configured datasets
  (keep the existing DATA_DIR confinement as defence in depth).
- Error messages shown to users must not echo server filesystem paths or tracebacks
  (log them server-side instead; see `utils/error_handling.py`).

### T1.3 Concurrency limits
- `utils/data_loading.py::load_adata`: concurrent cold misses for the same file must load
  once (per-file `threading.Lock`, double-checked around the cache). Keep `lru_cache`
  semantics otherwise.
- DE runs (`callbacks/deg_callbacks.py` → `utils/deg.py`): at most ONE DE fit at a time
  process-wide (`threading.Semaphore(1)`, non-blocking acquire). If busy, return a clear
  user message ("Another differential-expression run is in progress; please retry in a
  minute.") instead of queueing.

## Tier 2 — wrong results shown to users

### T2.4 Heatmap hover text is transposed
`utils/plotting.py` (~L295-310 categorical facets, ~L405-420 continuous): `z=H.T` but
`customdata=hover_text` is built `[i][j]` over `H` (x-major). Build hover data in the same
orientation as `z` (i.e. transpose it). Check every `go.Heatmap` in `utils/plotting.py`
and `utils/wholemount.py` for the same mismatch.

### T2.5 Percentile clipping must affect colour only
Currently each gene's grid values are clipped before plotting and before choosing the
shared maximum (`utils/plotting.py` ~L516 and ~L1405), so hover values are capped and
equal true means can render differently on a "shared" scale. Change to: keep the data
(z and hover) unclipped; compute the per-gene clip value (the percentile) and apply it as
the colour-scale upper limit (`zmax`/`cmax`). For the shared scale, use one common
`zmax` = max of the per-gene clip values, applied to all panels. The "Max:" badge (if
any) must still report the colour-scale maximum it reports today.

### T2.6 Signed continuous metadata
When colouring by a continuous metadata column that has negative values (e.g. DV.Score,
NT.Score), do not force `zmin=0` / `cmin=0` (`utils/plotting.py` ~L392, sphere bins in
`utils/wholemount.py`); use the data range (or a symmetric range around 0). Gene
expression (non-negative) keeps `zmin=0`.

### T2.7 Differential expression correctness (`utils/deg.py`, `callbacks/deg_callbacks.py`)
- `get_raw_counts`: verify the count source. For `.raw`, sample up to ~10k nonzero values
  and require non-negative and integer-like (|v - round(v)| < 1e-2, the same tolerance as
  the pseudobulk panel check in `utils/plotting.py` ~L1556-1565). For the expm1×nCount
  reconstruction, apply the same check to the reconstructed sample. Remove the
  "last resort: assume .X already holds counts" branch. If no verified count source →
  raise `DEGError("raw counts unavailable for this dataset; differential expression is
  disabled")`.
- Replicate labels: cells whose value is missing (NaN/None/empty) in ANY replicate
  column must be excluded from pseudobulking (never stringified into a 'nan' replicate),
  and the number excluded reported in the status message. Same fix in the group-figure
  pseudobulk path (`callbacks/main_callbacks.py` ~L116).
- Design fallback: keep the manuscript design `~ library + condition` (matches Fig 6H
  method). Remove the broad `except Exception` downgrade to `~condition`. Instead, BEFORE
  fitting, check identifiability: if every library contributes pseudobulks to only one
  side (library fully confounded with condition), fit `~ condition` and say so explicitly
  in the status ("library term not estimable: each library is on one side only").
  Otherwise fit `~ library + condition`; if pydeseq2 raises, surface it as a DE error
  (logged server-side), do not silently refit.
- Leave `min_frac` behaviour unchanged (documented later).

### T2.8 Citations (README.md, datasets_config.yml descriptions/data_source, any UI text)
Verified against GEO/PubMed:
- GSE138002 → "Lu et al. 2020" (not "Sridhar et al.")
- GSE142244 → "Ghinia Tegla et al. 2020" (not "Emerson et al.")
- GSE234963 → add "Dorgau et al. 2024"
- GSE246169 → add "Wohlschlegel et al. 2023"
Replace every occurrence; keep accession numbers as they are.

### T2.9 UI wording
- Visible text must not cite reviewer-response figure numbers: replace "Fig R2.5"
  in labels/tooltips/buttons (e.g. "↺ Reset projection to Fig R2.5 default" →
  "↺ Reset projection to default"; "Default 82° reproduces Fig R2.5" → "Default 82°").
  Code comments may keep provenance notes.
- DE help text (`components/deg_panel.py` ~L4, ~L20) and docstring
  (`callbacks/deg_callbacks.py` ~L6): do not claim "library × genotype" for all datasets;
  say the replicate unit is the dataset's configured replicate columns (chick:
  library × genotype; human/mouse: library), ideally rendered from the dataset's
  `replicate_columns`.
- Replace the unprofessional comment in `utils/plotting.py` ~L41 ("awful").

## Verification (must pass before you finish)
1. Write `tests/test_hardening.py`, runnable as `python tests/test_hardening.py` (plain
   asserts + `if __name__ == "__main__"`, no pytest), covering with synthetic inputs:
   coerce helpers (out-of-range, NaN, strings, huge ints, malformed ROI, >500 vertices);
   heatmap hover orientation (asymmetric synthetic grid: hovered value at (x,y) equals z
   at (x,y)); shared-scale clipping (two genes [1,1,1,10] vs [10,10,10,10] at
   percentile 0.5 → equal values get equal z, hover unclipped); signed metadata zmin < 0;
   `get_raw_counts` refuses log-normalised `.raw` and a no-counts object, accepts integer
   `.raw` and a valid log1p(CP)+nCount reconstruction; NaN replicate excluded; confounded
   design → `~condition` with explicit message; load_adata concurrent misses → one load
   (mock the reader); download route refuses unconfigured names.
2. `python tests/test_hardening.py`, `python test_nasal_overlay.py`,
   `python test_smooth_sigma_resolve.py` all pass.
3. `python -c "import app"` succeeds (imports, callback registration).

## Constraints
- Work only inside this worktree. No git commit/push, no branch changes. Leave changes
  uncommitted.
- Do not start the server. Do not read or modify any .h5ad or anything outside the worktree.
- Python: use `/Users/chlee/repos/_retina_viewer_venv/bin/python` for all commands.
- Final message: per spec item, what changed (file:line) and anything you could not do
  or judged differently, with reasons.
