"""Pseudobulk differential expression for user-drawn ROIs.

The interactive generalization of the manuscript area-DEG volcano
(scripts/figures/fig6h_sfig23_area_deg.R in the analysis repo): instead of a fixed
DV/NT grid quadrant, the user lassos one or two regions of interest on the
topographic map, and we run the SAME class of test the figure used.

Pipeline, mirroring the figure:
  1. Resolve two (possibly overlapping) ROIs into disjoint per-cell side labels.
     One ROI  -> region vs. the rest (the figure's `area1` contrast).
     Two ROIs -> A vs. B; an overlapping cell goes to the SMALLER ROI and the
     larger ROI cedes it, so no cell feeds both pseudobulk groups (which would
     bias the contrast).
  2. Sum RAW counts into per-(side x replicate) pseudobulks, where the replicate
     unit is the demux embryo (library x genotype for chick) -- the same unit the
     figure pseudobulked on. Drop pseudobulks with < min_cells cells.
  3. Keep genes detected in >= min_frac of pseudobulks.
  4. Test side A vs B with a negative-binomial GLM (pydeseq2 -- the Python analog
     of the figure's glmGamPoi) under design ~library + condition when the library
     batch term is estimable, else ~condition. Positive log2FC = enriched in A.

Numbers will not be bit-identical to the R figure (a different NB engine), but the
methodology matches: pseudobulk, library-adjusted, region contrast.
"""
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Matches _REP_SEP in callbacks.main_callbacks: the library x genotype replicate key
# is joined on this control char (cannot appear in obs string values).
SEP = '\x1f'

# pydeseq2 parallelises gene fitting; keep the worst-case (large ROI, ~15k genes) under
# the gunicorn request timeout. Conservative cap so it plays nice on a small HF Space.
_N_CPUS = max(1, min(int(os.getenv('DEG_N_CPUS', '4')), (os.cpu_count() or 1)))


class DEGError(Exception):
    """Base for expected, user-facing DEG refusals (caught in the callback)."""


class NoReplicate(DEGError):
    pass


class InsufficientReplicates(DEGError):
    def __init__(self, n_a, n_b, min_cells):
        self.n_a, self.n_b, self.min_cells = n_a, n_b, min_cells
        super().__init__(f"{n_a} A / {n_b} B pseudobulks >= {min_cells} cells")


def resolve_rois(idx_a, idx_b, n_obs):
    """Disjoint per-cell side labels for two (possibly overlapping) ROIs.

    Returns (labels, info). `labels` is an object array of 'A' / 'B' / '' (unused).
    Smaller ROI keeps contested cells; the larger cedes the overlap. With exactly
    one ROI populated, the other side is the rest = every remaining cell.
    """
    a = np.zeros(n_obs, dtype=bool)
    b = np.zeros(n_obs, dtype=bool)
    if idx_a:
        a[np.asarray(idx_a, dtype=int)] = True
    if idx_b:
        b[np.asarray(idx_b, dtype=int)] = True
    na, nb = int(a.sum()), int(b.sum())
    labels = np.full(n_obs, '', dtype=object)

    if na and nb:
        overlap = a & b
        # Smaller ROI keeps the contested cells; the larger cedes them.
        if na <= nb:
            b = b & ~overlap
        else:
            a = a & ~overlap
        labels[a] = 'A'
        labels[b] = 'B'
        info = {'mode': 'A_vs_B', 'n_overlap': int(overlap.sum())}
    elif na or nb:
        # One ROI vs the rest. Keep the drawn region as the 'A' side so the DE direction is
        # stable (positive log2FC = enriched in the drawn region), but record which button the
        # user actually used ('solo') so the recap/status/volcano name + colour it correctly --
        # otherwise a B-only draw is reported as "ROI A".
        roi = a if na else b
        labels[roi] = 'A'
        labels[~roi] = 'B'   # B == the rest
        info = {'mode': 'A_vs_rest', 'n_overlap': 0, 'solo': 'A' if na else 'B'}
    else:
        info = {'mode': None, 'n_overlap': 0}

    info['n_A'] = int((labels == 'A').sum())
    info['n_B'] = int((labels == 'B').sum())
    return labels, info


def polygon_to_indices(verts, xs, ys):
    """Indices of cells whose (xs[i], ys[i]) fall inside the closed polygon `verts`
    ([[x, y], ...] in the SAME coordinate space as xs/ys, e.g. NT.Score × DV.Score).

    This is what makes the ROI work on ANY rendering of that space (per-cell scatter,
    binned heatmap, smoothed map): selection is a point-in-polygon test against the
    cells' own coordinates, independent of how the view is drawn. A polygon with < 3
    vertices selects nothing.
    """
    if not verts or len(verts) < 3:
        return []
    from matplotlib.path import Path
    path = Path(np.asarray(verts, dtype=float))
    pts = np.column_stack([np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)])
    return np.where(path.contains_points(pts))[0].tolist()


def replicate_labels(obs, replicate_columns):
    """Per-cell replicate label from the configured columns (e.g. library x genotype),
    joined on SEP. None if no configured column is present (-> pseudobulk impossible)."""
    cols = [c for c in (replicate_columns or []) if c in obs.columns]
    if not cols:
        return None
    return obs[cols].astype(str).agg(SEP.join, axis=1).to_numpy()


def get_raw_counts(adata, norm_target=None):
    """(counts cells x genes, gene_names) of RAW counts.

    Prefer adata.raw (true integer counts, e.g. chick). Otherwise reconstruct from a
    log1p(CP) .X via expm1 * nCount / norm_target -- the SAME reconstruction the group
    figure's pseudobulk panel uses -- keeping the matrix sparse. Counts are rounded to
    integers at the pseudobulk level (aggregate_pseudobulk), as pydeseq2 requires.
    """
    if adata.raw is not None:
        return adata.raw.X, np.asarray(adata.raw.var_names)
    X = adata.X
    if norm_target and 'nCount_RNA' in adata.obs.columns:
        nc = adata.obs['nCount_RNA'].to_numpy().astype(float)
        scale = nc / float(norm_target)
        if sp.issparse(X):
            r = sp.csr_matrix(X, copy=True).astype(float)
            r.data = np.expm1(r.data)               # log1p(CP) -> CP (sparsity preserved)
            r = sp.diags(scale) @ r                 # row-scale by depth/target
        else:
            r = np.expm1(np.asarray(X, dtype=float)) * scale[:, None]
        return r, np.asarray(adata.var_names)
    # Last resort: assume .X already holds counts.
    return X, np.asarray(adata.var_names)


def aggregate_pseudobulk(counts, genes, replicate, side, min_cells=50):
    """Sum raw counts into per-(side, replicate) pseudobulks.

    counts: cells x genes (sparse or ndarray), raw counts. replicate: per-cell label.
    side: per-cell 'A'/'B'/''. Returns (pb_df samples x genes, meta_df). meta columns:
    side, replicate, library, n_cells -- where library is the leading SEP-component of
    the replicate key (the batch covariate; == replicate for a single-column key).
    Pseudobulks with < min_cells cells are dropped.
    """
    counts = counts.tocsr() if sp.issparse(counts) else np.asarray(counts)
    keep = np.isin(side, ('A', 'B'))
    rows = np.where(keep)[0]
    if rows.size == 0:
        return pd.DataFrame(), pd.DataFrame()

    grp_key = np.array([f"{side[i]}{SEP}{replicate[i]}" for i in rows])
    uniq, inv = np.unique(grp_key, return_inverse=True)
    inv = np.asarray(inv).ravel()
    # group-indicator (n_groups x n_kept) @ counts_kept -> n_groups x genes
    G = sp.csr_matrix((np.ones(rows.size), (inv, np.arange(rows.size))),
                      shape=(uniq.size, rows.size))
    sub = counts[rows]
    pb = G @ sub
    pb = np.asarray(pb.todense()) if sp.issparse(pb) else np.asarray(pb)
    n_cells = np.asarray(G.sum(axis=1)).ravel().astype(int)

    parts = [u.split(SEP, 2) for u in uniq]   # [side, lib, rest...] ; lib = leading rep col
    sides = np.array([p[0] for p in parts])
    reps = np.array([SEP.join(p[1:]) for p in parts])
    libs = np.array([p[1] if len(p) > 1 else p[0] for p in parts])
    display_idx = np.array([u.replace(SEP, '|') for u in uniq])

    meta = pd.DataFrame({'side': sides, 'replicate': reps, 'library': libs,
                         'n_cells': n_cells}, index=display_idx)
    ok = (meta['n_cells'] >= min_cells).to_numpy()
    meta = meta[ok]
    pb = np.rint(pb[ok]).astype(np.int64)
    pb_df = pd.DataFrame(pb, index=meta.index, columns=np.asarray(genes))
    return pb_df, meta


def filter_genes(pb_df, min_frac=0.5):
    """Keep genes detected (>0) in >= min_frac of pseudobulk samples."""
    if pb_df.shape[0] == 0:
        return pb_df
    frac = (pb_df > 0).mean(axis=0)
    return pb_df.loc[:, frac >= min_frac]


def choose_design(meta):
    """Pick '~library + condition' when the library batch term is estimable and not
    confounded with side; else '~condition'. Returns (design_str, reason)."""
    libs = meta['library']
    n = len(meta)
    n_lib = libs.nunique()
    # A library must span BOTH sides for `library` to be a batch term rather than a
    # proxy for side. And we need residual df: params = intercept + (n_lib-1) + 1(cond).
    spans_both = any(set(meta.loc[libs == L, 'side']) >= {'A', 'B'} for L in libs.unique())
    if n_lib >= 2 and spans_both and n > n_lib + 1:
        return "~library + condition", f"library-adjusted ({n_lib} libraries)"
    if n_lib < 2:
        reason = "single library"
    elif not spans_both:
        reason = "library confounded with ROI"
    else:
        reason = "too few replicates to adjust for library"
    return "~condition", reason


def run_deg(pb_df, meta, min_frac=0.5, min_cells=50):
    """Pseudobulk NB-GLM DE (pydeseq2): side A vs B. Returns (results_df, info).

    results_df: columns gene, baseMean, log2FoldChange, pvalue, padj (positive
    log2FoldChange = enriched in ROI A). Raises InsufficientReplicates when a side has
    < 2 pseudobulks (DE is not meaningful) -- the cross-species honesty guardrail.
    """
    pb_df = filter_genes(pb_df, min_frac)
    per_side = meta['side'].value_counts()
    n_a, n_b = int(per_side.get('A', 0)), int(per_side.get('B', 0))
    if n_a < 2 or n_b < 2:
        raise InsufficientReplicates(n_a, n_b, min_cells)
    if pb_df.shape[1] == 0:
        raise DEGError("No genes passed the detection filter.")

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    md = meta[['side', 'library']].rename(columns={'side': 'condition'}).copy()
    design, reason = choose_design(meta)
    try:
        dds = DeseqDataSet(counts=pb_df, metadata=md, design=design,
                           n_cpus=_N_CPUS, quiet=True)
        dds.deseq2()
    except Exception:
        # Full-rank / estimability backstop: drop the library covariate and retry.
        design, reason = "~condition", reason + " (library term not estimable)"
        dds = DeseqDataSet(counts=pb_df, metadata=md, design="~condition",
                           n_cpus=_N_CPUS, quiet=True)
        dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "A", "B"], n_cpus=_N_CPUS, quiet=True)
    st.summary()
    res = (st.results_df[["baseMean", "log2FoldChange", "pvalue", "padj"]]
           .reset_index().rename(columns={'index': 'gene'}))
    info = {'design': design, 'design_reason': reason,
            'n_pb_A': n_a, 'n_pb_B': n_b, 'n_genes_tested': int(res.shape[0])}
    return res, info


def deg_from_labels(adata, replicate_columns, labels, min_cells=50, min_frac=0.5,
                    norm_target=None):
    """End-to-end DE from resolved per-cell side `labels`. Returns (results_df, info)."""
    rep = replicate_labels(adata.obs, replicate_columns)
    if rep is None:
        raise NoReplicate()
    counts, genes = get_raw_counts(adata, norm_target)
    pb_df, meta = aggregate_pseudobulk(counts, genes, rep, labels, min_cells=min_cells)
    n_pb_a = int((meta.get('side') == 'A').sum()) if len(meta) else 0
    n_pb_b = int((meta.get('side') == 'B').sum()) if len(meta) else 0
    if n_pb_a < 2 or n_pb_b < 2:
        raise InsufficientReplicates(n_pb_a, n_pb_b, min_cells)
    return run_deg(pb_df, meta, min_frac=min_frac, min_cells=min_cells)


# --------------------------------------------------------------------------- #
# Self-check: `python -m utils.deg` (or `python utils/deg.py`). Asserts the
# overlap rule and that a spiked gene is recovered through the full pipeline,
# via both the .raw path (chick) and the log1p(CP) reconstruction path.
def _demo():
    import anndata as ad

    # 1) Overlap resolution: smaller ROI keeps the contested cells.
    labels, info = resolve_rois(idx_a=[0, 1, 2], idx_b=[2, 3, 4, 5, 6], n_obs=10)
    assert labels[2] == 'A', "overlap must go to the smaller ROI (A)"
    assert info['mode'] == 'A_vs_B' and info['n_overlap'] == 1
    assert list(labels[[3, 4, 5, 6]]) == ['B'] * 4 and info['n_A'] == 3 and info['n_B'] == 4
    # one ROI -> vs rest
    lab1, info1 = resolve_rois([0, 1], [], 6)
    assert info1['mode'] == 'A_vs_rest' and (lab1 == 'B').sum() == 4

    # 2) Full pipeline on a synthetic dataset with a clean spike-in.
    rng = np.random.default_rng(0)
    n_cells, n_genes = 600, 40
    counts = rng.poisson(20, size=(n_cells, n_genes)).astype(np.int64)
    # 3 libraries x 2 genotypes spanning the whole field (so library spans both sides)
    obs = pd.DataFrame({
        'library': rng.choice(['L1', 'L2', 'L3'], n_cells),
        'genotype': rng.choice(['d0', 'd1'], n_cells),
        'NT.Score': rng.uniform(-1, 1, n_cells),
        'DV.Score': rng.uniform(-1, 1, n_cells),
    })
    obs['nCount_RNA'] = counts.sum(1)
    region_a = obs['NT.Score'].to_numpy() > 0.3        # ROI A = nasal stripe
    counts[region_a, 0] += 150                         # gene g0 strongly up in A
    var = pd.DataFrame(index=[f"g{j}" for j in range(n_genes)])

    # polygon-in-score-space selection reproduces the boolean region (the ROI mechanism)
    poly = [[0.3, -1.1], [1.1, -1.1], [1.1, 1.1], [0.3, 1.1]]
    poly_idx = polygon_to_indices(poly, obs['NT.Score'], obs['DV.Score'])
    assert set(poly_idx) == set(np.where(region_a)[0]), "polygon selection must match NT>0.3"

    adata = ad.AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    adata.raw = adata                                  # exercise the .raw path
    lab, _ = resolve_rois(np.where(region_a)[0].tolist(),
                          np.where(~region_a)[0].tolist(), n_cells)
    res, info = deg_from_labels(adata, ['library', 'genotype'], lab, min_cells=20)
    g0 = res.set_index('gene').loc['g0']
    assert g0['log2FoldChange'] > 1 and g0['padj'] < 0.05, f"spike not recovered: {g0.to_dict()}"
    assert 'library' in info['design'], f"expected library-adjusted design, got {info['design']}"
    print(f"[.raw path] design={info['design']!r}  g0 lfc={g0['log2FoldChange']:.2f} "
          f"padj={g0['padj']:.1e}  genes={info['n_genes_tested']}")

    # 3) Reconstruction path: log1p(CP10K) .X, no .raw.
    cp = counts / counts.sum(1, keepdims=True) * 1e4
    adata2 = ad.AnnData(X=np.log1p(cp).astype(np.float32), obs=obs.copy(), var=var)
    res2, info2 = deg_from_labels(adata2, ['library', 'genotype'], lab,
                                  min_cells=20, norm_target=1e4)
    g0b = res2.set_index('gene').loc['g0']
    assert g0b['log2FoldChange'] > 1 and g0b['padj'] < 0.05, f"recon spike: {g0b.to_dict()}"
    print(f"[recon path] g0 lfc={g0b['log2FoldChange']:.2f} padj={g0b['padj']:.1e}")
    print("deg.py self-check passed.")


if __name__ == '__main__':
    _demo()
