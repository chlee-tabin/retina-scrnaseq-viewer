// Click-to-click polygon ROI drawing for the differential-expression view.
//
// Plotly has no native click-to-vertex polygon tool (lasso/draw are all drag), so we
// attach our own click / double-click listeners: each click converts its pixel position
// to data coordinates (NT.Score / DV.Score) and adds a polygon vertex; a double-click
// closes the current polygon. Vertices are mirrored into the roi-vertices-store, where
// the server reads them for a point-in-polygon cell selection -- so the ROI works on
// ANY rendering of that coordinate space (per-cell scatter, binned heatmap, smoothed).
window.dash_clientside = window.dash_clientside || {};

(function () {
  var S = { mode: null, verts: { A: [], B: [] }, gd: null, pending: [],
            onClick: null, onDbl: null, boundGd: null };

  var COL = { A: { line: 'rgba(44,127,184,0.95)', fill: 'rgba(44,127,184,0.15)' },
              B: { line: 'rgba(217,95,14,0.95)',  fill: 'rgba(217,95,14,0.15)' } };

  // ponytail: 280ms click/double-click disambiguation window — a touch over the typical
  // inter-click gap so a closing double-click cancels its own two pending vertices;
  // raise it if fast closes still leave a stray vertex.
  var COMMIT_MS = 280;

  function gdEl() {
    var el = document.getElementById('main-plot');
    return el ? (el.querySelector('.js-plotly-plot') || el) : null;
  }

  // The x/y axis pair whose pixel extent contains (px, py). Generalises to subplots
  // (e.g. the two-gene comparison); null when the click is outside every plot area
  // (margins, modebar, title) so such clicks never drop a stray vertex.
  function axesAt(gd, px, py) {
    var fl = gd._fullLayout, xa = null, ya = null;
    for (var k in fl) {
      var ax = fl[k];
      if (!ax || !ax._id || ax._offset == null || ax._length == null) continue;
      if (ax._id.charAt(0) === 'x' && px >= ax._offset && px <= ax._offset + ax._length) xa = ax;
      if (ax._id.charAt(0) === 'y' && py >= ax._offset && py <= ax._offset + ax._length) ya = ax;
    }
    return (xa && ya) ? [xa, ya] : null;
  }

  function pixelToData(gd, evt) {
    var bb = gd.getBoundingClientRect();
    var px = evt.clientX - bb.left, py = evt.clientY - bb.top;
    var pair = axesAt(gd, px, py);
    if (!pair) return null;
    return [pair[0].p2d(px - pair[0]._offset), pair[1].p2d(py - pair[1]._offset)];
  }

  function pathStr(v) {
    var s = 'M ' + v[0][0] + ',' + v[0][1];
    for (var i = 1; i < v.length; i++) s += ' L ' + v[i][0] + ',' + v[i][1];
    if (v.length > 2) s += ' Z';
    return s;
  }

  // react-plotly imports plotly as a module, so window.Plotly may be undefined until the
  // first Graph mounts. Resolve it lazily and no-op until it's ready.
  function plotly() { return window.Plotly || null; }

  // The x/y subplot pairs in the figure: [{xref:'x',yref:'y'}, {xref:'x2',yref:'y2'}, ...].
  // "Compare a second gene" shows the same DV/NT space twice, so we draw the overlay on every
  // pair -- the ROI then appears on whichever panel was clicked (and its twin) rather than
  // always being relaid out onto the left panel.
  function subplotRefs(gd) {
    var fl = gd._fullLayout, xs = [];
    for (var k in fl) {
      var ax = fl[k];
      if (ax && ax._id && ax._id.charAt(0) === 'x' && xs.indexOf(ax._id) < 0) xs.push(ax._id);
    }
    if (!xs.length) xs = ['x'];
    return xs.map(function (x) { return { xref: x, yref: 'y' + x.slice(1) }; });
  }

  function draw() {
    var P = plotly();
    var gd = gdEl();
    // Re-fetch the live gd and require it to be fully initialised: relayout on a div that
    // is still (re)rendering throws "_guiEditing" (e.g. the redraw-on-figure hook racing a
    // main-plot re-render). Bail quietly; the next sync redraws.
    if (!gd || !gd._fullLayout || !P) return;
    S.gd = gd;
    var refs = subplotRefs(gd);
    var shapes = [];
    ['A', 'B'].forEach(function (k) {
      var v = S.verts[k];
      if (!v || !v.length) return;
      refs.forEach(function (r) {
        shapes.push({ type: 'path', path: pathStr(v), xref: r.xref, yref: r.yref, layer: 'above',
          line: { color: COL[k].line, width: 2 },
          fillcolor: v.length > 2 ? COL[k].fill : 'rgba(0,0,0,0)' });
      });
    });
    try { P.relayout(gd, { shapes: shapes }); } catch (e) { /* gd mid-render; a later sync redraws */ }
  }

  function sync() {
    window.dash_clientside.set_props('roi-vertices-store',
      { data: JSON.parse(JSON.stringify(S.verts)) });
  }

  function onClick(evt) {
    if (!S.mode) return;
    var pt = pixelToData(S.gd, evt);
    if (!pt) return;
    var m = S.mode, entry = {};
    // Delay each vertex commit so the two clicks of a closing double-click can be
    // cancelled. Distinct vertices land at different spots (no dblclick) and commit.
    entry.timer = setTimeout(function () {
      S.verts[m].push(pt); draw(); sync();
      S.pending = S.pending.filter(function (e) { return e !== entry; });
    }, COMMIT_MS);
    S.pending.push(entry);
  }

  function onDbl(evt) {
    if (!S.mode) return;
    evt.preventDefault();
    evt.stopPropagation();
    // A double-click fires two click events at the same spot. Cancel just the second one and
    // let the first commit (via its timer, in click order) as the closing vertex -- so
    // following the displayed "double-click the final vertex" workflow keeps that corner
    // instead of dropping it. Earlier pending single clicks still commit normally.
    var dup = S.pending.pop();
    if (dup) clearTimeout(dup.timer);
    window.dash_clientside.set_props('roi-draw-mode-store', { data: null });   // close
  }

  function detach() {
    if (S.boundGd && S.onClick) {
      S.boundGd.removeEventListener('click', S.onClick, true);
      S.boundGd.removeEventListener('dblclick', S.onDbl, true);
    }
    S.onClick = S.onDbl = null;
    S.boundGd = null;
  }

  // Attach the click / double-click listeners to a graph div, detaching from any previous one.
  function bindTo(gd) {
    detach();
    S.onClick = onClick; S.onDbl = onDbl;
    S.boundGd = gd;
    gd.addEventListener('click', S.onClick, true);
    gd.addEventListener('dblclick', S.onDbl, true);
  }

  window.dash_clientside.roi = {
    // Attach / detach the drawing listeners as the draw mode changes.
    setMode: function (mode) {
      detach();
      S.gd = gdEl();
      S.mode = mode || null;
      var P = plotly();
      if (!S.gd) return '';
      if (S.mode) {
        bindTo(S.gd);
        if (P) P.relayout(S.gd, { dragmode: false });   // clicks place vertices, no box/lasso/zoom
      } else if (P) {
        P.relayout(S.gd, { dragmode: 'zoom' });          // restore normal map interaction on exit
      }
      return '';
    },
    // React to external store changes (e.g. the Clear button empties the vertices).
    syncVerts: function (v) {
      S.verts = (v && v.A && v.B) ? v : { A: [], B: [] };
      var gd = gdEl();
      if (gd) S.gd = gd;
      // A map re-render (view / control change) while drawing replaces the Plotly div and
      // drops our listeners; re-bind to the live div so clicks keep placing vertices.
      if (S.mode && gd && gd !== S.boundGd) {
        bindTo(gd);
        var P = plotly(); if (P) { try { P.relayout(gd, { dragmode: false }); } catch (e) {} }
      }
      draw();
      return '';
    }
  };
})();
