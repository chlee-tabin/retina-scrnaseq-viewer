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


def create_volcano_figure(res, lfc_thresh=1.0, padj_thresh=0.05, top_n=12, subtitle=None):
    """Volcano from a results frame with gene / log2FoldChange / padj columns.
    Positive log2FC = enriched in ROI A."""
    res = res.copy()
    for c in ('padj', 'log2FoldChange'):
        res[c] = pd.to_numeric(res[c], errors='coerce')
    r = res.dropna(subset=['padj', 'log2FoldChange']).copy()
    if r.empty:
        return _message_figure("No genes passed filtering for this contrast.")

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
    fig.add_hline(y=-np.log10(padj_thresh), line=dict(color='grey', dash='dash', width=1))

    # Label the most significant genes (by padj) among the flagged set.
    for _, row in r[sig].nsmallest(top_n, 'padj').iterrows():
        fig.add_annotation(x=row['log2FoldChange'], y=row['nlp'], text=row['gene'],
                           showarrow=False, font=dict(size=10), yshift=9)

    title = "Volcano — positive log2FC = enriched in ROI A"
    if subtitle:
        title += f"<br><sub>{subtitle}</sub>"
    fig.update_layout(
        title=title, xaxis_title="log2 fold change (A / B)",
        yaxis_title="-log10 adjusted p", template='plotly_white',
        height=520, showlegend=False, margin=dict(t=70))
    return fig


def create_recap_figure(x, y, labels, mode, x_label='NT.Score', y_label='DV.Score'):
    """Scatter coloured by the resolved ROI assignment (A / B-or-rest / unused) -- the
    disjoint selection that was actually contrasted, in the embedding's axes."""
    other = 'rest' if mode == 'A_vs_rest' else 'ROI B'
    disp = np.where(labels == 'A', 'ROI A', np.where(labels == 'B', other, '(unused)'))
    df = pd.DataFrame({'x': np.asarray(x), 'y': np.asarray(y), 'color': disp})
    cmap = {'ROI A': '#2c7fb8', 'ROI B': '#d95f0e', 'rest': '#bdbdbd', '(unused)': '#e8e8e8'}
    order = ['(unused)', 'rest', 'ROI B', 'ROI A']   # ROIs drawn last (on top)
    fig = create_scatter_plot(
        df, f'{x_label} vs {y_label}', 'ROI assignment',
        treat_as_categorical=True, color_map=cmap, category_order=order,
        plot_order='none')
    fig.update_layout(title="Selected cells (resolved A / B assignment)", margin=dict(t=50))
    return fig
