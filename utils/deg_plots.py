"""Volcano + ROI-recap figures for the interactive ROI differential-expression view.

The volcano mirrors the manuscript area-DEG styling (x = log2FC, y = -log10 adj p,
flag |lfc| > 1 & adj p < 0.05, label the top genes); the recap shows the resolved
disjoint A/B assignment on the topographic map, so it is transparent which cells were
contrasted (the "illustrate what got selected, above the volcano" requirement).
"""
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from utils.plotting import create_scatter_plot, _message_figure


def create_volcano_figure(res, lfc_thresh=1.0, padj_thresh=0.05, top_n=12, subtitle=None,
                          foreground_label='A', mode='A_vs_B'):
    """Volcano from a results frame with gene / log2FoldChange / padj columns.
    Positive log2FC = enriched in the drawn ('A'-side) region -- named by `foreground_label`
    ('A' normally, 'B' for a B-only draw)."""
    res = res.copy()
    for c in ('padj', 'log2FoldChange'):
        res[c] = pd.to_numeric(res[c], errors='coerce')
    n_na = int(res['padj'].isna().sum())
    contrast = f"ROI {foreground_label} / rest" if mode == 'A_vs_rest' else 'ROI A / ROI B'
    axis_label = f"log2FC ({contrast})"
    na_note = f"{n_na:,} genes with padj = NA not shown" if n_na else ''
    r = res.dropna(subset=['padj', 'log2FoldChange']).copy()
    if r.empty:
        fig = _message_figure("No genes passed filtering for this contrast.")
        fig.update_layout(xaxis_title=axis_label, title=na_note)
        return fig

    # -log10(padj), with padj==0 clipped to the smallest positive value for a finite y.
    pos = r['padj'][r['padj'] > 0]
    floor = float(pos.min()) if not pos.empty else 1e-300
    r['nlp'] = -np.log10(r['padj'].clip(lower=floor))
    sig = (r['log2FoldChange'].abs() > lfc_thresh) & (r['padj'] < padj_thresh)

    fig = go.Figure()
    for mask, color, name in ((~sig, 'lightgrey', 'ns'), (sig, '#d6604d', 'sig')):
        sub = r[mask]
        fig.add_trace(go.Scattergl(
            x=sub['log2FoldChange'], y=sub['nlp'], mode='markers',
            marker=dict(color=color, size=6 if name == 'sig' else 5),
            name=name, text=sub['gene'],
            hovertemplate='%{text}<br>log2FC=%{x:.2f}<br>-log10 padj=%{y:.2f}<extra></extra>'))

    fig.add_vline(x=lfc_thresh, line=dict(color='grey', dash='dash', width=1))
    fig.add_vline(x=-lfc_thresh, line=dict(color='grey', dash='dash', width=1))
    # Only draw the significance line for a positive threshold -- padj_thresh == 0 (which
    # flags nothing) would make -log10 infinite and break the render.
    if padj_thresh and padj_thresh > 0:
        fig.add_hline(y=-np.log10(padj_thresh), line=dict(color='grey', dash='dash', width=1))

    # Label the most significant genes (by padj) among the flagged set.
    for _, row in r[sig].nsmallest(top_n, 'padj').iterrows():
        fig.add_annotation(x=row['log2FoldChange'], y=row['nlp'], text=row['gene'],
                           showarrow=False, font=dict(size=10), yshift=9)

    if subtitle is None:
        subtitle = (f"|log2FC| &gt; {lfc_thresh:g} · adj p &lt; {padj_thresh:g} "
                    f"· {int(sig.sum()):,} significant")
    if na_note:
        subtitle = f"{subtitle} · {na_note}" if subtitle else na_note
    title = f"Volcano — positive log2FC = enriched in ROI {foreground_label}"
    if subtitle:
        title += f"<br><sub>{subtitle}</sub>"
    fig.update_layout(
        title=title, xaxis_title=axis_label,
        yaxis_title="-log10 adjusted p", template='plotly_white',
        height=520, showlegend=False, margin=dict(t=70))
    return fig


def create_recap_figure(x, y, labels, mode, x_label='NT.Score', y_label='DV.Score',
                        foreground_label='A'):
    """Scatter coloured by the resolved ROI assignment (A / B-or-rest / unused) -- the
    disjoint selection that was actually contrasted, in the embedding's axes. `foreground_label`
    names the drawn 'A'-side region: 'A' normally, 'B' for a B-only draw so it shows + colours
    as ROI B (matching the orange polygon on the map) rather than being mislabelled ROI A."""
    fg = f'ROI {foreground_label}'
    other = 'rest' if mode == 'A_vs_rest' else 'ROI B'
    disp = np.where(labels == 'A', fg, np.where(labels == 'B', other, '(unused)'))
    df = pd.DataFrame({'x': np.asarray(x), 'y': np.asarray(y), 'color': disp})
    cmap = {'ROI A': '#2c7fb8', 'ROI B': '#d95f0e', 'rest': '#bdbdbd', '(unused)': '#e8e8e8'}
    order = ['(unused)', 'rest', 'ROI B', 'ROI A']   # ROIs drawn last (on top)
    fig = create_scatter_plot(
        df, f'{x_label} vs {y_label}', 'ROI assignment',
        treat_as_categorical=True, color_map=cmap, category_order=order,
        plot_order='none')
    fig.update_layout(title="Selected cells (resolved A / B assignment)", margin=dict(t=50))
    return fig
