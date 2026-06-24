"""Whole-mount ("flower" / orange-peel) reprojection of the (DV.Score, NT.Score) map.

A standalone, pure-numpy port of `flower_reproj/flower_reproject.py:flower_transform`
from the analysis repo (retina-spatial-scrna-analysis), shipped there as a reviewer
response: the flat DV/NT score plane is warped to resemble a dissected, flat-mounted
retinal whole-mount, so peripheral relationships are not artefactually compressed.

Geometric model (verbatim from the source):
  - centre at the HAA pole (score origin); scale each axis robustly; go to polar
    (r, theta) with NT -> horizontal (nasal +x, right) and DV -> vertical (dorsal +y, up).
  - treat the score-radius as orthographic foreshortening of a spherical cap and invert
    it (`dewarp='arcsin'`: rho = arcsin(r * sin(rho_max))) so the periphery expands.
  - lay out azimuthal-equidistant: display radius R = rho (colatitude).
  - cut K radial slits ("gores"); within each gore open a wedge gap carrying the
    curvature deficit, so the petals separate toward the rim (orange-peel relief).
  - anisotropy: rho_max larger toward NT, so the nasal/temporal petals run longer.

This module provides the per-cell coordinate transform (flower_transform, a pure function
of the two score arrays + numpy) plus compute_scale_fit -- a once-per-render fit of the
pole and per-axis scale FROM THE CELLS, reused across the scatter and binned warps so the
two views share one basis. The matplotlib binning/rendering in the source script is not
ported -- the Dash viewer reuses its own binning/smoothing pipeline on the warped (X, Y).

Orientation (FISH flat-mount convention): nasal = right, dorsal = top, temporal = left,
ventral = bottom; HAA at centre.
"""
import numpy as np
import pandas as pd

# Coordinate-relevant defaults from flower_reproject.DEFAULT_PARAMS (the render-only
# missing_nasal_* keys are dropped; see the module docstring).
DEFAULT_PARAMS = dict(
    pole=(0.0, 0.0),               # HAA pole = score origin (balanced D/V & N/T module scores)
    symmetric=True,                # equalize N/T and D/V extents via per-direction scaling
    scale="p99",                   # used only when symmetric=False
    rho_max_dv_deg=64.0,           # cap half-angle toward dorsal/ventral
    rho_max_nt_deg=82.0,           # cap half-angle toward nasal/temporal (stretched wider)
    dewarp="arcsin",               # "arcsin" | "pow" | "none"
    pow_p=1.6,
    cut_angles_deg=(0, 90, 180, 270),  # slits at cardinals -> 4 lobes at diagonals (DN,DT,VT,VN)
    gap_gain=0.85,                 # deficit-mode rip width
    gap_mode="linear",             # "linear" (visible V-notch from rip-start) | "deficit" (curvature-true)
    gap_frac=0.5,                  # linear mode: max gap as fraction of the half-gore width (at rim)
    rho_join_deg=18.0,             # keep lobes joined within this colatitude (solid HAA center) = rip-start
    stretch_bumps=(),              # selective directional radial stretch: [(angle_deg, amplitude, width_deg), ...]
    rot_deg=0.0,
)

# Parameters of the shipped reviewer figure (flower_reproj LOCKED_PARAMS,
# Fig_R2.5_flower_reprojection.png). The `missing_nasal` hatched cap is a render-only
# overlay (not a coordinate of any cell), so it is intentionally absent here.
# Stretch-bump angles: temporal=180 deg, ventro-temporal=225 deg, dorso-nasal=45 deg.
LOCKED_PARAMS = dict(
    rho_max_nt_deg=82.0,
    rho_max_dv_deg=64.0,
    gap_gain=1.0,
    gap_mode="deficit",
    rho_join_deg=14.0,
    dewarp="arcsin",
    cut_angles_deg=(0, 90, 180, 270),
    stretch_bumps=((180.0, 0.60, 50.0), (225.0, 0.85, 50.0), (45.0, 0.40, 50.0)),
)


def _resolve_pole(dv, nt, pole):
    if pole == "median":
        return float(np.nanmedian(dv)), float(np.nanmedian(nt))
    if pole == "origin":
        return 0.0, 0.0
    return float(pole[0]), float(pole[1])


def _resolve_scale(x, y, scale):
    if scale == "p99":
        return (float(np.nanpercentile(np.abs(x), 99)) or 1.0,
                float(np.nanpercentile(np.abs(y), 99)) or 1.0)
    return float(scale[0]), float(scale[1])


def _sym_factors(v, q=97.0):
    """Per-side scale factors (sp, sn): the q-th percentile of the positive and the
    abs-negative values. Exposed so the factors can be computed ONCE from the cells and
    reused when warping the bin grid -- keeping the scatter and binned views on one basis
    (see compute_scale_fit)."""
    v = np.asarray(v, float)
    pos = v[v > 0]; neg = -v[v < 0]
    sp = (float(np.nanpercentile(pos, q)) if pos.size else 1.0) or 1.0
    sn = (float(np.nanpercentile(neg, q)) if neg.size else 1.0) or 1.0
    return sp, sn


def _apply_sym(v, sp, sn):
    """Apply per-side factors: positives / sp, abs-negatives / sn -> ~[-1, 1]."""
    v = np.asarray(v, float)
    return np.where(v >= 0, v / sp, v / sn)


def _sym_scale(v, q=97.0):
    """Scale + and - sides independently to their q-th percentile -> ~[-1,1] symmetric.
    Equalizes opposing petals (e.g. nasal vs temporal) despite asymmetric score ranges."""
    return _apply_sym(v, *_sym_factors(v, q))


def flower_transform(dv, nt, **p):
    """Map (DV.Score, NT.Score) -> (X, Y). Returns (X, Y, diag) with per-point gore/rho."""
    P = {**DEFAULT_PARAMS, **p}
    dv = np.asarray(dv, float); nt = np.asarray(nt, float)
    dv0, nt0 = _resolve_pole(dv, nt, P["pole"])
    x = nt - nt0          # NT -> horizontal (nasal +x, right)
    y = dv - dv0          # DV -> vertical   (dorsal +y, up)
    if P.get("symmetric", True):
        # Reuse cell-derived per-side factors when injected (compute_scale_fit), so the
        # bin-grid warp shares the cells' basis; otherwise fit to this array.
        sf = P.get("sym_factors")
        if sf is not None:
            (spx, snx), (spy, sny) = sf
            xs, ys = _apply_sym(x, spx, snx), _apply_sym(y, spy, sny)
        else:
            xs, ys = _sym_scale(x), _sym_scale(y)
        sx = sy = None
    else:
        sx, sy = _resolve_scale(x, y, P["scale"])
        xs, ys = x / sx, y / sy
    th = np.arctan2(ys, xs) + np.deg2rad(P["rot_deg"])
    r = np.clip(np.hypot(xs, ys), 0.0, 1.0)

    w_nt = np.cos(th) ** 2
    rho_max = np.deg2rad(P["rho_max_nt_deg"]) * w_nt + np.deg2rad(P["rho_max_dv_deg"]) * (1.0 - w_nt)
    if P["dewarp"] == "arcsin":
        rho = np.arcsin(np.clip(r * np.sin(rho_max), -1.0, 1.0))
    elif P["dewarp"] == "pow":
        rho = rho_max * (r ** P["pow_p"])
    else:
        rho = rho_max * r
    R = rho.copy()

    cuts = np.deg2rad(np.sort(np.asarray(P["cut_angles_deg"], float) % 360.0))
    K = len(cuts)
    phi = th % (2 * np.pi)
    gore = np.zeros_like(r, dtype=int)
    if K > 0:
        rj = np.deg2rad(P["rho_join_deg"])
        if P.get("gap_mode", "linear") == "deficit":
            with np.errstate(invalid="ignore", divide="ignore"):
                deficit = 2 * np.pi * (1.0 - np.where(rho > 1e-9, np.sin(rho) / rho, 1.0))
            deficit = np.clip(deficit * P["gap_gain"], 0.0, 2 * np.pi * 0.9)
            delta = deficit / (2 * K)
        else:  # "linear": V-notch that opens from the rip-start (rho_join) out to the rim
            half_gore = np.pi / K
            frac = np.clip((rho - rj) / np.maximum(rho_max - rj, 1e-6), 0.0, 1.0)
            delta = P.get("gap_frac", 0.5) * half_gore * frac
        thm = th % (2 * np.pi)
        thm = np.where(thm < cuts[0], thm + 2 * np.pi, thm)
        edges = np.concatenate([cuts, [cuts[0] + 2 * np.pi]])
        gore = np.clip(np.searchsorted(edges, thm, side="right") - 1, 0, K - 1)
        c0 = edges[gore]; c1 = edges[gore + 1]
        width = c1 - c0
        t = (thm - c0) / np.where(width > 0, width, 1.0)
        span = np.clip(width - 2 * delta, 1e-3, None)
        phi = c0 + delta + t * span
    # selective directional stretch (anatomical direction = th): elongate chosen sectors
    if P.get("stretch_bumps"):
        st = np.ones_like(R)
        for ang, amp, wid in P["stretch_bumps"]:
            d = np.angle(np.exp(1j * (th - np.deg2rad(ang))))   # signed angular diff in (-pi, pi]
            st = st + amp * np.exp(-0.5 * (d / np.deg2rad(wid)) ** 2)
        R = R * st
    X = R * np.cos(phi); Y = R * np.sin(phi)
    # `theta` is the raw azimuth (pre-gore, pre-stretch); `phi` is after the relief gaps.
    # The sphere view reuses (rho, theta) -- its native, un-ripped geometry (see sphere_coords).
    return X, Y, dict(r=r, rho=rho, theta=th, phi=phi, gore=gore, pole=(dv0, nt0), scale=(sx, sy))


def compute_scale_fit(dv, nt, **p):
    """Pole + per-axis scale derived ONCE from the cells, to be merged into the params and
    reused for every flower_transform call (the cells AND the bin-grid vertices/centres), so
    the scatter and binned whole-mount share one basis and the layout is stable under cell
    subsetting. Returns a dict to splat into the flower_transform kwargs.

    Without this, each call re-fits its scale to whatever array it receives -- the cells for
    the scatter, the regular grid for the binned map -- so the two views' petals diverge.
    """
    P = {**DEFAULT_PARAMS, **p}
    dv = np.asarray(dv, float); nt = np.asarray(nt, float)
    dv0, nt0 = _resolve_pole(dv, nt, P["pole"])
    x = nt - nt0
    y = dv - dv0
    fit = {"pole": (dv0, nt0)}
    if P.get("symmetric", True):
        fit["sym_factors"] = (_sym_factors(x), _sym_factors(y))
    else:
        fit["scale"] = _resolve_scale(x, y, P["scale"])
    return fit


# Default obs column names carrying the topographic scores (consistent across the
# chick / human / mouse RPC datasets in datasets_config.yml).
DV_COL = "DV.Score"
NT_COL = "NT.Score"


def wholemount_coords(dv, nt, params=None):
    """Convenience wrapper: (DV, NT) score arrays -> (X, Y) whole-mount coordinates.

    Uses the shipped LOCKED_PARAMS by default (the reviewer-figure look). Non-finite
    score rows propagate to NaN coordinates, which the viewer's scatter/binning then
    drop -- so a dataset where only RPCs carry scores simply shows those cells.
    """
    P = dict(LOCKED_PARAMS if params is None else params)
    X, Y, _ = flower_transform(dv, nt, **P)
    return np.asarray(X, dtype=float), np.asarray(Y, dtype=float)


def sphere_coords(dv, nt, params=None):
    """(DV, NT) score arrays -> (x, y, z) on the UNIT SPHERE: the topographic map on its
    native (near-spherical) geometry, of which the flat flower is the azimuthal-equidistant
    projection (display radius = colatitude). The colatitude `rho` (the arcsin-dewarped
    spherical-cap angle) and azimuth `theta` come straight from flower_transform, so the
    sphere shares the flower's pole + per-axis scale basis when the same `params` (with a
    compute_scale_fit) are passed.

    The flat layout's relief cuts / wedge gap / petal stretch are FLATTENING devices for a
    curved surface (they relieve the curvature deficit); on the sphere there is no deficit,
    so they are intentionally ignored here -- the surface is continuous and un-ripped.

    Orientation matches the flower: HAA pole at the top (+z, rho=0); nasal -> +x, dorsal ->
    +y. Cells only reach rho <= rho_max (~64-82 deg), so they cover a polar CAP; the
    uncaptured periphery and the nasal gap are simply empty sphere (no cell there).
    """
    P = dict(LOCKED_PARAMS if params is None else params)
    _, _, diag = flower_transform(dv, nt, **P)
    rho = np.asarray(diag['rho'], dtype=float)
    th = np.asarray(diag['theta'], dtype=float)
    sin_rho = np.sin(rho)
    return sin_rho * np.cos(th), sin_rho * np.sin(th), np.cos(rho)


def marker_footprint_center(dv, nt, present, *, bin_size=50, min_cells=5):
    """Geometric centre of a sparse marker's footprint in (DV, NT) score space: the
    equal-weight centroid of the score-grid BINS that contain >=1 'present' cell (and clear the
    cell floor). Returns (dv, nt) or None.

    Bin-PRESENCE, not per-cell density: each occupied bin counts once, so the centre is the
    middle of where the marker is detected (its "bull's-eye" / ring centre) and is NOT dragged
    toward whichever side happens to have more detected cells. This is the robust landmark for a
    marker seen in <1% of cells (chick CYP26C1 -> HAA): a per-cell median/mean, or a per-bin
    expression-weighted centre, chases the densest / brightest specks (which sparsity scatters
    to the periphery), whereas the footprint centroid stays at the spatial centre. `min_cells`
    matches the displayed binned map so the landmark sits on bins the user actually sees.
    """
    dv = np.asarray(dv, dtype=float); nt = np.asarray(nt, dtype=float)
    present = np.asarray(present, dtype=bool)
    fin = np.isfinite(dv) & np.isfinite(nt)
    if not fin.any():
        return None
    nb = int(bin_size)
    cnt, ne, de = np.histogram2d(nt[fin], dv[fin], bins=nb)          # x=NT, y=DV (as the map)
    pcnt, _, _ = np.histogram2d(nt[fin & present], dv[fin & present], bins=[ne, de])
    foot = (pcnt > 0) & (cnt >= max(1, int(min_cells)))
    if not foot.any():
        return None
    ntc = (ne[:-1] + ne[1:]) / 2; dvc = (de[:-1] + de[1:]) / 2
    NTc, DVc = np.meshgrid(ntc, dvc, indexing='ij')
    return float(DVc[foot].mean()), float(NTc[foot].mean())


# HAA-pointer modes (label shown in the control + the caption description). The marker is
# sparse and ring-ish, so "the centre" legitimately depends on weighting -- expose the choice.
HAA_MODE_DESC = {
    'footprint':  'detection-footprint centre (centroid of marker+ bins)',
    'expression': 'expression-weighted centre (binned-mean centre of mass)',
    'domain':     'expression-domain centre (bins ≥ 50% of peak; FISH-calibrated gate)',
    'peak':       'expression peak (top bins, ≥ 95th-pctile of bin means)',
}


def haa_center(dv, nt, expr, mode='domain', *, bin_size=50, min_cells=5, frac=0.5,
               peak_pctile=95.0, smooth_sigma=0.0):
    """(dv, nt) centre of a marker's HAA in score space, by one of several definitions; None if
    undefined. The marker (chick: CYP26C1) is detected in <1% of cells and ring-ish, so the
    estimators legitimately differ -- the viewer exposes `mode` so users can explore:

      'footprint'  equal-weight centroid of marker+ bins  -> geometric centre of detection
      'expression' binned-mean centre of mass             -> centre of mass of the signal
      'domain'     centroid of bins >= frac*peak, peak = `peak_pctile` of per-bin means
                   -> the FISH-calibrated spatial-domain gate of retina-spatial-scrna-analysis
                      (item2_fishcalib_perbin_domains: frac=0.5, peak=95th pctile) -- robust
      'peak'       centroid of the top bins (>= peak_pctile of per-bin means) -> the hot spot

    All modes bin on the (NT,DV) score grid (x=NT, y=DV) and gate on bins with >= min_cells, so
    the landmark sits on bins the displayed map actually shows. `expr` is the marker's per-cell
    value (log-norm); footprint uses only presence (expr>0). `smooth_sigma` Gaussian-smooths the
    per-bin mean the SAME way the displayed binned map does -- pass the dataset's smooth_sigma so
    the expression/domain/peak centre lands on the bright core the user actually SEES (computing
    on the raw, unsmoothed means leaves the diamond a few bins off the smoothed blob)."""
    dv = np.asarray(dv, dtype=float); nt = np.asarray(nt, dtype=float); expr = np.asarray(expr, dtype=float)
    if mode == 'footprint':
        return marker_footprint_center(dv, nt, expr > 0, bin_size=bin_size, min_cells=min_cells)
    fin = np.isfinite(dv) & np.isfinite(nt) & np.isfinite(expr)
    if not fin.any():
        return None
    nb = int(bin_size)
    esum, ne, de = np.histogram2d(nt[fin], dv[fin], bins=nb, weights=np.clip(expr[fin], 0, None))
    cnt, _, _ = np.histogram2d(nt[fin], dv[fin], bins=[ne, de])
    valid = cnt >= max(1, int(min_cells))
    with np.errstate(invalid='ignore', divide='ignore'):
        bm = np.where(valid, esum / np.maximum(cnt, 1.0), np.nan)
    if smooth_sigma and smooth_sigma > 0:
        # Mirror _binned_mean: fill invalid with 0, smooth, then re-mask -> the displayed field.
        from scipy.ndimage import gaussian_filter
        bm = gaussian_filter(np.where(valid, np.nan_to_num(bm), 0.0), sigma=float(smooth_sigma))
        bm = np.where(valid, bm, np.nan)
    finite_pos = bm[np.isfinite(bm) & (bm > 0)]
    if finite_pos.size == 0:
        return None
    ntc = (ne[:-1] + ne[1:]) / 2; dvc = (de[:-1] + de[1:]) / 2
    NTc, DVc = np.meshgrid(ntc, dvc, indexing='ij')
    if mode == 'expression':
        w = np.where(np.isfinite(bm), np.clip(bm, 0.0, None), 0.0)
    else:  # 'domain' or 'peak': threshold on the robust (percentile) peak of the bin means
        peak = float(np.percentile(finite_pos, float(peak_pctile)))
        thr = frac * peak if mode == 'domain' else peak
        w = (np.isfinite(bm) & (bm >= thr)).astype(float)
    if w.sum() <= 0:
        return None
    return float((DVc * w).sum() / w.sum()), float((NTc * w).sum() / w.sum())


def has_scores(obs_columns):
    """True iff both topographic score columns are present (so the projection applies).

    `obs_columns` is typically a pandas Index (``adata.obs.columns``); build the set from
    it directly -- ``Index or []`` raises ("truth value of an Index is ambiguous").
    """
    cols = set(obs_columns) if obs_columns is not None else set()
    return DV_COL in cols and NT_COL in cols


if __name__ == "__main__":
    # ponytail: one runnable check for the sphere geometry (no test framework in this repo).
    # Run with `python -m utils.wholemount`.
    dv = np.array([0.0, 1.0, -1.0, 0.5, -0.7, 0.0])
    nt = np.array([0.0, 0.8, 0.3, -0.6, -0.2, 1.0])
    x, y, z = sphere_coords(dv, nt)
    assert np.allclose(x**2 + y**2 + z**2, 1.0), "points must lie on the unit sphere"
    assert np.allclose([x[0], y[0], z[0]], [0, 0, 1]), "HAA pole (score origin) -> +z top"
    xn, yn, _ = sphere_coords(np.array([0.0]), np.array([1.0]))   # pure nasal (+NT)
    assert xn[0] > 0 and abs(yn[0]) < 1e-9, "nasal -> +x (matches flower orientation)"
    xd, yd, _ = sphere_coords(np.array([1.0]), np.array([0.0]))   # pure dorsal (+DV)
    assert yd[0] > 0 and abs(xd[0]) < 1e-9, "dorsal -> +y (matches flower orientation)"
    print("wholemount.sphere_coords self-check OK")
