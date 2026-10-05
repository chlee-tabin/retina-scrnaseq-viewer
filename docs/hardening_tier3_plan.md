# Public-release hardening — Tier 3 (robustness, reproducibility, docs)

Second PR, on top of the merged Tier 1+2 hardening (`docs/hardening_tier12_plan.md` —
reuse its validation helpers; do not duplicate them). Keep diffs minimal and in the style
of the surrounding code. No new runtime dependencies (pytest may be added as a dev/test
dependency only).

## A. Share links and data identity
A1. **Data version in share state.** Add a `data_version` field to every share link: the
    configured `file_path` basename of the dataset being viewed (from
    `datasets_config.yml`, resolved server-side). On restore, if the link's
    `data_version` differs from the dataset id's current file, still restore but show a
    visible notice: "This link was created with an earlier version of this dataset
    (<old>); it now shows <new>." Links without `data_version` (older links) restore
    silently as today. Policy (write it in docs/VIEWER_SPEC.md): dataset ids are stable;
    replacing a dataset's data requires a new id or this notice.
A2. **Single owner for the share-state schema.** Create one `STATE_SCHEMA` table in
    `utils/state.py` (key → component id, property, default, coerce function, schema
    version introduced). Generate the share-link `State`s, the restore `Output`s, and the
    restore values from it (replace the hardcoded `n_out = 34` positional tuple in
    `callbacks/url_callbacks.py`). Migrations live in `decode_state` (e.g. v2→v3 fills
    `smooth_sigma` via the existing `LEGACY_SMOOTH_SIGMA`). Fix the "schema v2" comment
    next to `'v': 3`. Keep every existing link working (round-trip test).
A3. **Omitted keys restore defaults.** When restoring a link, keys absent from the link
    reset their control to the schema default (not `no_update`), so a single-gene link
    does not keep the recipient's comparison gene / shared scale. Exception: keys
    introduced after the link's schema version keep the documented legacy behaviour.

## B. Data provisioning and paths
B1. **Pin the data.** Add to `datasets_config.yml` a top-level
    `hf_data_revision: "2c21052aad5203a64b5798474ed2045fca8e3733"` and per-dataset `sha256`:
    - 20250604_chick_RPC.h5ad 7e5f363aeb58f0c9093d73cb6f82532c36c51a130bd9a4d1edd25274ea0308ea
    - 20250604_human_RPC.h5ad b9d0253f12f4a8d494c46f1d521c678de519dcdcef829a7a4a9612fb534dd6fd
    - 20250604_mouse_RPC.h5ad bc04f446cc5c5a0ddc4a0067b1547734427fc9f53bc74a3ee680b651c435a2a5
    - 20260528_mouse_RPC_cr9_e13e16.h5ad 3fe08118d68a2b2bc0d2b0d5ded6ec93cd01b7e69e407058dd498a22aa23a2ef
    - 20260620_chick_full.h5ad 62bd0d860c8e44c4e4b4947dc19c205ba158bbe0d7574ed042e6d4575c52679d
    `utils/data_provision.py`: pass `revision=` to `hf_hub_download` (env
    `HF_DATA_REVISION` overrides the config); after download verify sha256 (streamed) and
    fail loudly on mismatch (delete the bad file). Files already present locally are
    verified once and the result cached by (path, size, mtime) so restarts stay fast.
B2. **One path rule.** Loading, provisioning and the download route must resolve the
    same path for a config `file_path` (including nested `data/sub/x.h5ad`); duplicate
    basenames across datasets are rejected at config load with a clear error.
B3. Pin every entry in `requirements.txt` exactly (currently `scipy>=`,
    `huggingface_hub>=` are open) to the versions in the deployed environment — read them
    from `/Users/chlee/repos/_retina_viewer_venv` (`pip freeze`) only for packages
    already listed.

## C. Server / ops
C1. Per-browser session id (generate in a clientside callback or `uuid4` via a layout
    function) so the active-user count is real; poll interval 5 s → 30 s.
C2. Dockerfile: honour `PORT`/`HOST` (shell-form CMD, default 7860/0.0.0.0).
C3. Logging: stdout/stderr only by default; file logging (`RotatingFileHandler`, 5 MB ×
    3) only with the `-debug` flag. Do not log the whole `data_store`.
C4. `handle_callback_error` (`utils/error_handling.py`): after logging, raise
    `PreventUpdate` (or return `no_update` per output) instead of a bare `None`; never
    send tracebacks to the browser.
C5. Gene search: no options until ≥1 character typed; cap results at 200.
C6. Add `scripts/deploy_space.py`: uploads a clean checkout to the HF Space with
    `huggingface_hub.upload_folder` (token from env `HF_TOKEN`; refuses if the git tree
    is dirty), and writes the git SHA into an env/file the app shows in the footer
    ("viewer <version> · <short sha>"). Update the README deploy section to match (the
    Space is NOT deployed by git push).
C7. CI: `.github/workflows/tests.yml` running `pytest tests/` on push/PR (Python version
    from the Dockerfile). Tests must not need the .h5ad data.

## D. Analysis behaviour (decided by the author)
D1. **Smoothing:** keep the current zero-fill smoothing as the DEFAULT (it matches the
    manuscript pipeline). Add an opt-in "mask-normalised smoothing" toggle
    (normalised convolution: smooth(values·mask) / smooth(mask), masked bins stay
    masked); include it in the share state; document the difference in README + spec.
D2. **Whole-mount projection:** add a per-dataset `wholemount_params:` config key
    (overrides for `utils/wholemount.py::LOCKED_PARAMS`). chick datasets keep the current
    tuned values; human/mouse default to NO stretch bumps (stretch list empty / neutral).
    The "Reset projection" button resets to the current dataset's defaults.
D3. Figure-region presets (`figure_regions`): select cells with the half-open interval
    test on the binned DV/NT indices directly (as documented in the config), not via the
    polygon point-in-path test.
D4. Volcano: x-axis label from the actual contrast ("log2FC (ROI A / ROI B)",
    "(ROI A / rest)", "(ROI B / rest)"); subtitle notes "N genes with padj = NA not shown"
    when N > 0. In "vs rest" mode, "rest" excludes cells with non-finite plot coordinates.
D5. Group-figure replicate fallback (`callbacks/main_callbacks.py` ~L1007): remove
    `'genotype'` from `SAMPLE_HINTS` and search only standalone categorical columns.
D6. One categorical-column helper used by both `app.py` and
    `callbacks/main_callbacks.py` (non-numeric incl. pandas StringDtype, or small-int).
D7. Species switch: if the current gene is absent in the new dataset, try a
    case-insensitive match (VAX1→Vax1) before falling back; when a shared link's gene is
    missing, show "gene X not found in this dataset; showing Y".
D8. One `DEFAULT_MIN_CELLS_PER_BIN` constant (app.py ~L309 falls back to 1,
    main_callbacks ~L459/658 to 5). Config: give `chick_full` an explicit
    `smooth_sigma: 1.5`; add a one-line comment why human uses `min_cells_per_bin: 2`
    (lower cell density per bin); bump the config header `version`/`last_updated`.
D9. `utils/data_loading.py`: size the in-memory cache to the number of configured
    datasets that can be selected (or memory-aware), keeping the per-file load lock.

## E. Documentation and citation (outside readers)
E1. README: feature list adds ROI differential expression (polygon ROI → pseudobulk
    volcano, figure-region presets, exploratory), 3D sphere projection, HAA pointer,
    mask-normalised smoothing; layout block lists `docs/`, `assets/`, `scripts/`,
    `tests/`, deg callbacks. Env-var table adds `DEG_N_CPUS`, `HF_DATA_REVISION`; note
    `PORT`/`HOST` apply to Docker too (after C2); document `-debug`.
E2. README chick rows cite GEO **GSE322831** (public). Add "Citation" and "License (MIT)"
    sections: preprint https://doi.org/10.64898/2026.01.04.697548 (Joisher, Lee,
    Prabhakara, van der Weide, Si, Lonfat, Cepko), pointer to CITATION.cff.
    `layouts/sidebar.py` preprint URL → version-less https://doi.org/10.64898/2026.01.04.697548.
E3. CITATION.cff: add `version` (app version string), `date-released`, and a
    `preferred-citation` (type: article; the preprint above with all 7 authors in
    manuscript order: Heer N. V. Joisher, ChangHee Lee, Chaitra Prabhakara, Isabella van
    der Weide, Yichen Si, Nicholas Lonfat, Constance Cepko; doi 10.64898/2026.01.04.697548).
E4. docs/VIEWER_SPEC.md: drop "ultracode"; add config keys `figure_regions` (half-open
    bin semantics), `haa_marker`, `nasal_gap`, `wholemount_params`, `sha256`,
    `hf_data_revision`; share-state keys `roi`, `volcano_lfc`, `volcano_padj`,
    `deg_min_cells`, `data_version`, smoothing mode; controls incl. ROI/DEG/haa/nasal-gap/
    sphere; a short DE section (design `~ library + condition`, replicate unit =
    configured `replicate_columns`, `min_frac` gene filter is a viewer-only addition that
    changes padj vs the manuscript's tables); README and spec agree on required config
    keys; drop `deg_results_path` everywhere (dead).
E5. Replace references to internal-only paths/PRs ("analysis repo, PR #19",
    `scripts/viewer_full_chick/16b_npy_figure.R`, `scripts/figures/fig6h_sfig23_area_deg.R`,
    `scripts/analysis/area_significant_deg.R`) with the public analysis repo
    https://github.com/chlee-tabin/retina-spatial-scrna-analysis and its script names
    `scripts/figures/sfig24_area_deg.R` / `scripts/figures/supptables1_2_area_deg.R`;
    where no public script exists, cite "manuscript Methods" instead.
E6. `app.py` browser tab title → "Retina scRNA-seq Pattern Viewer".
E7. Config/help text: "14 individual donor eyes / libraries" → "14 libraries"; any
    remaining "Fig 6H"/"Fig S23" preset names → v15 numbering: chick preset "Fig 6H",
    human preset "Fig S24" (the human area DEG is Figure S24 in the revised manuscript).

## F. Tests
F1. Move `test_nasal_overlay.py`, `test_smooth_sigma_resolve.py`, and
    `tests/test_hardening.py` under `tests/` as pytest-compatible tests (keep the
    `__main__` runners); remove personal venv paths from docstrings; README "Testing"
    section: `pip install pytest && pytest tests/`.
F2. New tests: share-state round trip for every schema key + v2 link migration +
    omitted-key reset + data_version notice; sha256 verify (tiny temp file); path rule
    incl. nested + duplicate basename rejection; mask-normalised smoothing (isolated valid
    100% bin stays 100%); wholemount per-dataset params; figure-region interval selection
    edges; volcano label per mode; species-switch case-insensitive gene match.

## Verification
- `/Users/chlee/repos/_retina_viewer_venv/bin/python -m pytest tests/` passes (install
  pytest into that venv if missing — allowed).
- `python -c "import app"` succeeds.

## Constraints
Work only inside this worktree; no git commit/push/branch changes; do not start the
server; do not read/modify .h5ad files or anything outside the worktree (except
installing pytest into the named venv). Final message: per item, what changed
(file:line) and anything skipped or done differently, with reasons.
