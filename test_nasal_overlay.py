"""Geometry self-check for the exploratory uncaptured-nasal sphere overlay
(plotting._missing_nasal_overlay). Synthetic cells only -- no data file needed.

Run: ~/repos/_retina_viewer_venv/bin/python test_nasal_overlay.py
"""
import numpy as np

from utils import wholemount as wm
from utils.plotting import _missing_nasal_overlay, _missing_nasal_overlay_flat


def _disc_cells(n=80, rmax=0.9, nt_shift=0.0):
    """A dense square grid of (DV, NT) scores clipped to a disc -> a populated polar cap.
    nt_shift moves the disc along NT (nasal +): a negative shift = temporal-only, no nasal edge."""
    g = np.linspace(-1.0, 1.0, n)
    dv, nt = np.meshgrid(g, g)
    dv, nt = dv.ravel(), nt.ravel()
    keep = (dv ** 2 + (nt - nt_shift) ** 2) <= rmax ** 2
    return dv[keep], nt[keep] + nt_shift


def _rho_az(x, y, z):
    return np.arccos(np.clip(z, -1, 1)), np.arctan2(y, x)   # colatitude, azimuth


def test_off_switch():
    dv, nt = _disc_cells()
    assert _missing_nasal_overlay(dv, nt, None, frac=0.0) == []          # frac 0 -> nothing
    assert _missing_nasal_overlay(dv, nt, None, frac=0.13, layers=0) == []   # 0 layers -> off
    assert _missing_nasal_overlay(np.array([]), np.array([]), None, frac=0.13) == []


def test_wedge_shrinks_with_reach():
    # Narrowing the dorsal+ventral reach must keep strictly fewer band tiles (shrink to a
    # sub-arc of the nasal rim) -- the "shrink to certain locations" control.
    dv, nt = _disc_cells()
    wide = _missing_nasal_overlay(dv, nt, None, frac=0.13, layers=2, dorsal_deg=45, ventral_deg=55)
    narrow = _missing_nasal_overlay(dv, nt, None, frac=0.13, layers=2, dorsal_deg=15, ventral_deg=15)
    assert wide and narrow, "both wedges should still draw something"
    assert len(np.asarray(narrow[0].x)) < len(np.asarray(wide[0].x)), "narrow wedge should have fewer verts"


def test_band_geometry():
    dv, nt = _disc_cells()
    traces = _missing_nasal_overlay(dv, nt, None, frac=0.13, layers=2, bin_size=40)
    assert len(traces) == 2, "expect a Mesh3d band + a line outline"
    mesh = traces[0]
    vx, vy, vz = (np.asarray(a, float) for a in (mesh.x, mesh.y, mesh.z))

    # 1. every band vertex lies on the unit sphere.
    norm = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    assert np.allclose(norm, 1.0, atol=1e-9), f"verts off the sphere: norm range {norm.min()}..{norm.max()}"

    # 2. the band is on the NASAL side (+x) -- centroid x > 0.
    assert vx.mean() > 0, f"band not nasal: mean x = {vx.mean()}"

    # 3. the band sits BEYOND the cap's nasal rim: its max colatitude exceeds the
    #    nasal-wedge cells' max colatitude (extrusion to higher rho actually happened).
    cx, cy, cz = wm.sphere_coords(dv, nt)
    rho_c, az_c = _rho_az(cx, cy, cz)
    wedge = (np.cos(az_c) > 0) & (np.sin(az_c) <= np.cos(az_c)) & (np.sin(az_c) >= -1.4 * np.cos(az_c))
    rho_b, _ = _rho_az(vx, vy, vz)
    assert rho_b.max() > rho_c[wedge].max(), (
        f"band max rho {np.degrees(rho_b.max()):.1f}deg not beyond nasal rim "
        f"{np.degrees(rho_c[wedge].max()):.1f}deg")

    # 4. the band's inner edge is ANCHORED at the rim (not floating): its min rho is within
    #    a few degrees of some nasal-wedge cell's colatitude.
    assert rho_b.min() <= rho_c[wedge].max() + np.radians(5), "band detached from the cap rim"
    print(f"OK: {len(rho_b)} band verts, rho {np.degrees(rho_b.min()):.0f}..{np.degrees(rho_b.max()):.0f}deg, "
          f"nasal rim ~{np.degrees(rho_c[wedge].max()):.0f}deg")


def test_flat_overlay():
    # The flat-flower twin: one hatched fill trace, tiles nasal (+x in display) and extruded
    # beyond the captured cap's nasal extent; layers=0 / temporal-only -> nothing.
    dv, nt = _disc_cells()
    traces = _missing_nasal_overlay_flat(dv, nt, None, frac=0.13, layers=2, bin_size=40)
    assert len(traces) == 1 and traces[0].fill == 'toself', "expect one hatched fill trace"
    tx = np.array([v for v in traces[0].x if v is not None], float)
    assert tx.size and tx.mean() > 0, "flat tiles should be nasal (+x)"
    cx, _ = wm.wholemount_coords(dv, nt)
    assert np.nanmax(tx) > np.nanmax(cx[np.isfinite(cx)]), "tiles should extrude past the cap rim"
    assert _missing_nasal_overlay_flat(dv, nt, None, frac=0.13, layers=0) == []
    assert _missing_nasal_overlay_flat(*_disc_cells(nt_shift=-0.8), None, frac=0.13, layers=2, bin_size=40) == []


def test_temporal_only_has_no_band():
    # Cells pushed temporally (NT negative) -> no populated bin in the nasal wedge -> no band.
    dv, nt = _disc_cells(nt_shift=-0.8)
    assert _missing_nasal_overlay(dv, nt, None, frac=0.13, layers=2, bin_size=40) == []


if __name__ == "__main__":
    test_off_switch()
    test_wedge_shrinks_with_reach()
    test_band_geometry()
    test_flat_overlay()
    test_temporal_only_has_no_band()
    print("all geometry checks passed")
