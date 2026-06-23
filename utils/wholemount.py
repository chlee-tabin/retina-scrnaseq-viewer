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

This module contains ONLY the per-cell coordinate transform (a pure function of the two
score arrays + numpy); the matplotlib binning/rendering in the source script is not
ported -- the Dash viewer reuses its own binning/smoothing pipeline on the warped (X, Y).

Orientation (FISH flat-mount convention): nasal = right, dorsal = top, temporal = left,
ventral = bottom; HAA at centre.
"""
import numpy as np
import pandas as pd

# Defaults copied verbatim from flower_reproject.DEFAULT_PARAMS.
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
    return X, Y, dict(r=r, rho=rho, phi=phi, gore=gore, pole=(dv0, nt0), scale=(sx, sy))


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


def has_scores(obs_columns):
    """True iff both topographic score columns are present (so the projection applies).

    `obs_columns` is typically a pandas Index (``adata.obs.columns``); build the set from
    it directly -- ``Index or []`` raises ("truth value of an Index is ambiguous").
    """
    cols = set(obs_columns) if obs_columns is not None else set()
    return DV_COL in cols and NT_COL in cols
