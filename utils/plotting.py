import plotly.express as px
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def _message_figure(text):
    """Return an empty figure with a centered, user-friendly message.

    Mirrors callbacks.main_callbacks._message_figure but lives here so the
    plotting helpers have a self-contained fallback (no callback import).
    """
    fig = go.Figure()
    fig.add_annotation(
        text=text,
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=16, color="#666"),
    )
    fig.update_layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white",
    )
    return fig


# A distinct, colour-blind-aware qualitative palette for categorical colouring when
# a dataset does not configure an explicit `annotation_colors` map. Replaces plotly's
# pale, fast-recycling Set3 (the source of the "awful" full-chick categorical scheme).
CATEGORICAL_FALLBACK = px.colors.qualitative.Dark24


def _resolve_order(values, category_order=None, by_series=None):
    """Stable category order: a configured `category_order` first (present categories
    only), then the rest alphabetically. With no configured order, order by descending
    median of `by_series` when given, else alphabetically.
    """
    present = list(pd.unique(pd.Series(values).astype(str)))
    if category_order:
        head = [str(c) for c in category_order if str(c) in present]
        tail = sorted(c for c in present if c not in set(head))
        return head + tail
    if by_series is not None:
        med = (pd.DataFrame({'g': pd.Series(values).astype(str),
                             'v': np.asarray(by_series, dtype=float)})
               .groupby('g', observed=True)['v'].median()
               .sort_values(ascending=False))
        return med.index.tolist()
    return sorted(present)


def _resolve_color_map(order, color_map=None):
    """Resolve a full {category -> colour} map over `order`: keep the configured
    `color_map` colours and assign each unmapped category a deterministic fallback
    (advancing the palette by the number already assigned) so the scatter, violin and
    figure all colour an unmapped category identically -- even under a partial
    `annotation_colors`. Every view is passed this FULL map (no color_discrete_sequence),
    so plotly never applies its own fallback.
    """
    full = dict(color_map or {})
    for cat in (order or []):
        if cat not in full:
            full[cat] = CATEGORICAL_FALLBACK[len(full) % len(CATEGORICAL_FALLBACK)]
    return full


def _categorical_px_kwargs(col_name, color_map=None, order=None):
    """plotly-express kwargs for consistent categorical colour + order. When the render
    `order` is known we resolve a FULL colour map (configured colours + a deterministic
    fallback per unmapped category) and pass only that -- so every view (scatter, violin,
    figure) colours a given category identically even with a partial `annotation_colors`.
    """
    kw = {}
    if order:
        kw['category_orders'] = {col_name: order}
        kw['color_discrete_map'] = _resolve_color_map(order, color_map)
    else:
        kw['color_discrete_sequence'] = CATEGORICAL_FALLBACK
        if color_map:
            kw['color_discrete_map'] = dict(color_map)
    return kw


def _plot_order_perm(values, plot_order, treat_as_categorical):
    """Row permutation implementing the sidebar 'Plot Order' control for a scatter.

    Plotly draws points in data order, so the LAST point is on top. Hence:
      'ordered_asc'  -> sort ascending by value: the HIGHEST values are drawn last (on
                        top). The useful default for a sparse gene (e.g. NPY) -- the few
                        positive cells sit above the many zero cells instead of being
                        buried under them.
      'ordered_desc' -> sort descending: the lowest values are drawn on top.
      'random'       -> deterministic shuffle (FIXED seed) so points don't re-jump on
                        every unrelated callback, while still breaking the
                        acquisition-order occlusion of the raw cell order.
    Any other value (or <=1 point) -> identity (original order).

    For categorical colour the sort key is the string label; for continuous colour it is
    the numeric value (non-numeric/NaN coerced, sorted to the end of the ascending order).
    """
    n = len(values)
    if n <= 1 or plot_order not in ('ordered_asc', 'ordered_desc', 'random'):
        return np.arange(n)
    if plot_order == 'random':
        return np.random.RandomState(0).permutation(n)
    if treat_as_categorical:
        keys = pd.Series(values).astype(str).to_numpy()
    else:
        keys = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy()
    perm = np.argsort(keys, kind='stable')
    if plot_order == 'ordered_desc':
        perm = perm[::-1]
    return perm


def create_scatter_plot(df, embedding, color_by, treat_as_categorical=False, selection_data=None,
                        color_map=None, category_order=None, plot_order='random'):
    """
    Create a scatter plot with proper styling based on embedding type
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing 'x', 'y', and 'color' columns
    embedding : str
        Name of the embedding (affects axis labels and styling)
    color_by : str
        Name of the column used for coloring
    treat_as_categorical : bool
        Whether to treat the color column as categorical
    selection_data : dict
        Dictionary containing selection indices
    """
    # For custom embedding, use the actual column names
    if ' vs ' in embedding:
        x_label, y_label = embedding.split(' vs ')
    else:
        x_label = f'{embedding}_1'
        y_label = f'{embedding}_2'
    
    # Plot Order: a row permutation drawn on top last (computed on the ORIGINAL colour
    # values, before any categorical string-cast below).
    perm = _plot_order_perm(df['color'].to_numpy(), plot_order, treat_as_categorical)

    px_kwargs = dict(
        labels={'x': x_label, 'y': y_label, 'color': color_by},
        title=f'Single-cell visualization - {embedding}',
        hover_data=None,
    )
    if treat_as_categorical:
        df = df.copy()
        df['color'] = df['color'].astype(str)
        order = _resolve_order(df['color'], category_order=category_order)
        px_kwargs.update(_categorical_px_kwargs('color', color_map=color_map, order=order))
    else:
        px_kwargs['color_continuous_scale'] = 'viridis'

    # Apply the Plot Order (no-op when perm is the identity). iloc keeps the colour/x/y
    # rows aligned; px.scatter then draws in this new order.
    df = df.iloc[perm]

    fig = px.scatter(df, x='x', y='y', color='color', **px_kwargs)
    
    # Apply different styling based on embedding type
    if 'umap' in embedding.lower():
        fig.update_layout(
            plot_bgcolor='white',
            xaxis=dict(
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                scaleanchor="y",
                scaleratio=1,
            ),
            yaxis=dict(
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                scaleanchor="x",
                scaleratio=1,
            )
        )
    else:
        fig.update_layout(
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
            )
        )
    
    # Update selection styling if provided. Selection indices reference the ORIGINAL
    # cell order (pointIndex into the un-reordered data), so map them through `perm` to
    # the post-Plot-Order trace positions, keeping the highlight on the right cells.
    if selection_data and selection_data.get('indices'):
        inv = np.empty(len(perm), dtype=np.int64)
        inv[perm] = np.arange(len(perm))
        n = len(inv)
        sel = [int(inv[i]) for i in selection_data['indices']
               if isinstance(i, (int, np.integer)) and 0 <= i < n]
        fig.update_traces(
            selectedpoints=sel,
            selected=dict(marker=dict(opacity=1)),
            unselected=dict(marker=dict(opacity=0.1))
        )

    return fig

def create_metacell_plot(H, xedges, yedges, embedding):
    """
    Create a density-based visualization of cells
    """
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    
    if ' vs ' in embedding:
        x_label, y_label = embedding.split(' vs ')
    else:
        x_label = f'{embedding}_1'
        y_label = f'{embedding}_2'
    
    fig = px.imshow(
        H.T,
        x=x_centers,
        y=y_centers,
        labels={'x': x_label, 'y': y_label},
        title=f'Density visualization - {embedding}',
        aspect='equal'
    )
    
    # Apply UMAP-specific styling
    if 'umap' in embedding.lower():
        fig.update_layout(
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False)
        )
    
    return fig

def create_binned_plot(df, embedding, color_by, bin_size=50, percentile=0.95, treat_as_categorical=False, smooth_sigma=0, min_cells=1, color_floor=0.05, bin_stat='mean'):
    """
    Create a binned visualization of cells
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing 'x', 'y', and 'color' columns
    embedding : str
        Name of the embedding
    color_by : str
        Name of the column used for coloring
    bin_size : int
        Number of bins for both x and y axes
    percentile : float
        Percentile cutoff for color scale
    treat_as_categorical : bool
        Whether to treat the color column as categorical
    """
    if ' vs ' in embedding:
        x_label, y_label = embedding.split(' vs ')
    else:
        x_label = f'{embedding}_1'
        y_label = f'{embedding}_2'

    if treat_as_categorical:
        # Get unique categories and their counts
        categories = df['color'].unique()
        n_cols = min(4, len(categories))  # Allow up to 4 columns instead of 3
        n_rows = (len(categories) + n_cols - 1) // n_cols
        
        # Create subplot grid with tighter spacing
        fig = make_subplots(
            rows=n_rows, 
            cols=n_cols,
            subplot_titles=[str(cat) for cat in categories],
            horizontal_spacing=0.05,  # Reduced from 0.1
            vertical_spacing=0.05     # Reduced from 0.1
        )
        
        # Create a histogram for each category
        for idx, category in enumerate(categories):
            row = idx // n_cols
            col = idx % n_cols
            
            mask = df['color'] == category
            category_df = df[mask]
            
            H, xedges, yedges = np.histogram2d(
                category_df['x'], 
                category_df['y'], 
                bins=bin_size
            )
            
            # Keep original counts for hover text
            H_display = H.copy()
            
            # Create mask for empty bins and set them to NaN for grey color
            mask = H_display > 0
            H_display[~mask] = np.nan
            
            x_centers = (xedges[:-1] + xedges[1:]) / 2
            y_centers = (yedges[:-1] + yedges[1:]) / 2
            
            # Create hover data with bin indices and cell counts
            hover_text = [[
                f'Category: {category}<br>'
                f'Bin: ({i},{j})<br>'
                f'x: {x_centers[i]:.2f}<br>'
                f'y: {y_centers[j]:.2f}<br>'
                f'Cells: {int(H[i,j])}'
                for j in range(H.shape[1])
            ] for i in range(H.shape[0])]
            
            # Add heatmap for this category
            fig.add_trace(
                go.Heatmap(
                    z=H_display.T,
                    x=x_centers,
                    y=y_centers,
                    colorscale='viridis',
                    customdata=hover_text,
                    hovertemplate='%{customdata}',
                    hoverongaps=False,
                    showscale=False  # Hide individual colorbars
                ),
                row=row + 1,
                col=col + 1
            )
        
        # For categorical values, use fixed dimensions that scale reasonably with the number of categories
        base_height = 300  # Base height per row
        base_width = 350   # Base width per column
        
        # Calculate total dimensions with some reasonable limits
        max_rows = 4  # Maximum number of rows before scrolling
        n_rows = min(len(categories), max_rows)
        n_cols = (len(categories) + max_rows - 1) // max_rows  # Ceiling division
        
        total_height = base_height * n_rows
        total_width = base_width * n_cols
        
        # Create the figure with appropriate dimensions
        fig.update_layout(
            title=f'Binned visualization by {color_by} - {embedding}',
            showlegend=False,
            height=total_height,
            width=total_width,
            margin=dict(t=60, l=60, r=20, b=20),
            grid=dict(rows=n_rows, columns=n_cols)
        )
        
        # Update axes for all subplots with equal aspect ratio
        for i in range(1, n_rows + 1):
            for j in range(1, n_cols + 1):
                if (i-1) * n_cols + j <= len(categories):  # Only update existing subplots
                    fig.update_xaxes(
                        showgrid=False, 
                        zeroline=False, 
                        title_text=x_label,
                        constrain='domain',
                        row=i, 
                        col=j
                    )
                    fig.update_yaxes(
                        showgrid=False, 
                        zeroline=False, 
                        title_text=y_label, 
                        autorange=True,
                        scaleanchor=f"x{(i-1)*n_cols + j}",  # Link to corresponding x-axis
                        scaleratio=1,
                        row=i, 
                        col=j
                    )
    
    else:
        # For continuous values (including gene expression). The per-bin mean,
        # smoothing, and color-scale math are factored into _binned_mean so the
        # dual-gene side-by-side figure uses the identical computation.
        H, H_counts, xedges, yedges, color_vmax = _binned_mean(
            df['x'], df['y'], df['color'],
            bin_size=bin_size,
            percentile=percentile,
            smooth_sigma=smooth_sigma,
            min_cells=min_cells,
            color_floor=color_floor,
            stat=bin_stat,
        )
        is_pct = (bin_stat == 'frac_pos')
        color_label = f"{color_by} (% detected)" if is_pct else color_by

        x_centers = (xedges[:-1] + xedges[1:]) / 2
        y_centers = (yedges[:-1] + yedges[1:]) / 2

        # Create bin indices for hover text
        bin_indices_x = np.arange(len(x_centers))
        bin_indices_y = np.arange(len(y_centers))
        bin_indices_x, bin_indices_y = np.meshgrid(bin_indices_x, bin_indices_y)

        fig = px.imshow(
            H.T,  # Transpose to match the expected orientation
            x=x_centers,
            y=y_centers,
            labels={'x': x_label, 'y': y_label, 'color': color_label},
            title=f'Binned visualization - {embedding}',
            color_continuous_scale='viridis',
            zmin=0,
            zmax=color_vmax,
            aspect='equal'
        )

        # '% detected' is stored as a fraction in [0, 1]; render the colorbar ticks as
        # percentages (e.g. 0.2 -> "20%") so the scale matches the "(% detected)" label
        # and the hover read-out instead of showing bare fractions.
        if is_pct:
            fig.update_coloraxes(colorbar=dict(tickformat='.0%'))

        # Create hover data with bin indices and cell counts
        hover_text = [[
            f'Bin: ({i},{j})<br>'
            f'x: {x_centers[i]:.2f}<br>'
            f'y: {y_centers[j]:.2f}<br>'
            + (f'% detected: {H[i,j]*100:.1f}%<br>' if is_pct else f'Value: {H[i,j]:.2f}<br>')
            + f'Cells: {int(H_counts[i,j])}'
            for j in range(H.shape[1])
        ] for i in range(H.shape[0])]
        
        # Update traces with hover template and fix y-axis orientation
        fig.update_traces(
            customdata=hover_text,
            hovertemplate='%{customdata}',
            hoverongaps=False  # Disable hover for empty bins
        )
        
        # Fix y-axis orientation (minus at bottom, plus at top)
        fig.update_yaxes(autorange=True)  # This will ensure proper orientation

        # For continuous values (including gene expression), use a fixed pixel height
        fig.update_layout(
            title=f'Binned visualization - {embedding}',
            showlegend=False,
            height=800,  # Fixed pixel height instead of viewport units
            margin=dict(t=60, l=60, r=20, b=20)
        )

    return fig


def _binned_mean(x, y, vals, bin_size=50, percentile=0.95, smooth_sigma=0, min_cells=1,
                 color_floor=0.05, stat='mean'):
    """Compute a per-bin statistic of a continuous value over a 2D (x, y) grid.

    This is the shared spatial-binning math used by both create_binned_plot's
    continuous branch and create_dual_gene_figure, so the smoothing, masking and
    color-scale clipping stay identical across the single- and dual-gene views.

    Parameters
    ----------
    stat : {'mean', 'frac_pos'}
        'mean' -> per-bin mean of `vals` (default). 'frac_pos' -> per-bin
        fraction of cells with vals > 0 (the topographic "% positive / detected"
        map), in [0, 1].

    Returns
    -------
    H_mean : 2D float array
        Per-bin statistic (NaN for bins below the cell floor), clipped to vmax.
    H_counts : 2D float array
        Cells per bin.
    xedges, yedges : 1D arrays
        Bin edges from np.histogram2d.
    vmax : float
        Color-scale maximum.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    vals = np.asarray(vals, dtype=float)

    # Drop non-finite rows before histogramming: np.histogram2d raises on NaN/inf
    # coordinates ("autodetected range ... is not finite"), and a NaN value would
    # poison its bin's weighted sum. Custom-axis columns (QC metrics, subset-defined
    # scores) legitimately contain NaN, so this path is reachable via the axis picker.
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(vals)
    x, y, vals = x[finite], y[finite], vals[finite]
    if x.size == 0:
        nb = int(bin_size)
        edges = np.linspace(0.0, 1.0, nb + 1)
        empty = np.full((nb, nb), np.nan, dtype=float)
        return empty, np.zeros_like(empty), edges, edges, float(color_floor)

    # Weight by the value (mean) or by a detected-indicator (fraction positive).
    weights = (vals > 0).astype(float) if stat == 'frac_pos' else vals
    H, xedges, yedges = np.histogram2d(x, y, bins=bin_size, weights=weights)
    H_counts, _, _ = np.histogram2d(x, y, bins=bin_size)

    # Bins with at least `min_cells` cells are valid; the rest are greyed.
    valid = H_counts >= max(1, int(min_cells))

    # Per-bin statistic (mean, or fraction of cells positive); NaN elsewhere.
    H_mean = np.full(H.shape, np.nan, dtype=float)
    np.divide(H, H_counts, out=H_mean, where=valid)

    # Optional Gaussian smoothing (mirrors the analysis pipeline's spatial
    # images): fill invalid bins with 0, smooth, then re-mask them to NaN so
    # empty regions stay grey instead of bleeding into the smoother.
    if smooth_sigma and smooth_sigma > 0:
        filled = np.where(valid, H_mean, 0.0)
        H_mean = gaussian_filter(filled, sigma=float(smooth_sigma))
        H_mean[~valid] = np.nan

    # Color scale from the NONZERO valid bins (matches the analysis pipeline,
    # which clips over expressing pixels rather than all bins). For 'frac_pos'
    # the floor is small (so low-but-real detection still renders) and vmax is
    # capped at 1.0 (it is a fraction).
    floor = 0.01 if stat == 'frac_pos' else float(color_floor)
    pos = H_mean[valid & np.isfinite(H_mean) & (H_mean > 0)]
    if pos.size and percentile < 1.0:
        vmax = float(np.nanpercentile(pos, percentile * 100))
    elif pos.size:
        vmax = float(np.nanmax(pos))
    else:
        vmax = 0.0
    vmax = max(vmax, floor)
    if stat == 'frac_pos':
        vmax = min(vmax, 1.0)
    H_mean = np.where(H_mean > vmax, vmax, H_mean)

    return H_mean, H_counts, xedges, yedges, vmax


def create_group_expression_plot(df, gene, group_by, split_by=None, style='violin',
                                 color_map=None, category_order=None):
    """Expression of one gene across the categories of a .obs column.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain an 'expr' column (log-normalized expression), the
        `group_by` column, and optionally the `split_by` column.
    gene : str
        Gene name (for titles/axis labels).
    group_by : str
        Categorical column to group on (x-axis).
    split_by : str or None
        Optional second categorical column used for color/legend grouping.
    style : {'violin', 'box', 'strip'}
        Plot style. (The 'dotplot' style is handled by create_dotplot.)
    """
    if df is None or len(df) == 0 or group_by not in df.columns:
        return _message_figure("No data to plot for this selection.")

    # Drop cells with no group label or no expression value (avoids a stray
    # 'nan' category from astype(str), matching create_dotplot's behavior).
    plot_df = df[df[group_by].notna() & df['expr'].notna()].copy()
    if len(plot_df) == 0:
        return _message_figure("No data to plot for this selection.")

    # Stringify the categorical columns so plotly keeps the explicit order and
    # treats them as discrete.
    plot_df[group_by] = plot_df[group_by].astype(str)

    # Group order: a configured cell-type order if given, else descending median
    # expression so the most-expressing categories read left-to-right.
    order = _resolve_order(plot_df[group_by], category_order=category_order,
                           by_series=plot_df['expr'])

    # Colour by the split column if one is chosen; otherwise colour each group with
    # its fixed palette colour (consistent with the UMAP / dot-plot views).
    if split_by and split_by in plot_df.columns:
        plot_df[split_by] = plot_df[split_by].astype(str)
        color_arg = split_by
        color_kwargs = {'color_discrete_sequence': CATEGORICAL_FALLBACK,
                        'category_orders': {group_by: order}}
    else:
        color_arg = group_by
        color_kwargs = _categorical_px_kwargs(group_by, color_map=color_map, order=order)

    common = dict(
        x=group_by,
        y='expr',
        color=color_arg,
        title=f"{gene} expression by {group_by}",
        labels={'expr': f"{gene} (log-normalized)", group_by: group_by},
        **color_kwargs,
    )

    if style == 'box':
        # On ALL cells a zero-inflated gene's box collapses to a flat line at 0
        # (median = Q1 = Q3 = 0): only the whisker "range" shows, no IQR box. Box the
        # DETECTED cells (expr > 0) so the quartile box is meaningful; the title flags it.
        box_df = plot_df[plot_df['expr'] > 0]
        if len(box_df) == 0:
            return _message_figure(f"No {gene}-positive cells to box-plot for this selection.")
        box_common = {**common, 'title': f"{gene} expression by {group_by} (detected cells, expr > 0)"}
        fig = px.box(box_df, points=False, **box_common)
    elif style == 'strip':
        fig = px.strip(plot_df, **common)
    else:  # 'violin' (default)
        fig = px.violin(plot_df, box=True, points=False, **common)

    fig.update_layout(
        plot_bgcolor='white',
        height=600,
        margin=dict(t=60, l=60, r=20, b=120),
        xaxis=dict(tickangle=-40),
    )
    return fig


def create_dotplot(expr_df, genes, group_by, category_order=None):
    """Scanpy-style dot plot: gene set (rows) x categories (columns).

    Dot size encodes the fraction of cells with detectable expression
    (expr > 0); dot color encodes the mean expression across ALL cells in the
    group (the scanpy convention -- mean over all cells, not only expressing
    ones).

    Parameters
    ----------
    expr_df : pandas.DataFrame
        One column per gene in `genes` (log-normalized expression) plus the
        `group_by` column.
    genes : list of str
        Genes to show (y-axis).
    group_by : str
        Categorical column (x-axis).
    """
    genes = [g for g in (genes or []) if g in expr_df.columns]
    if not genes or group_by not in expr_df.columns:
        return _message_figure("No genes available for the dot plot.")

    records = []
    grouped = expr_df.groupby(group_by, observed=True)
    for grp, sub in grouped:
        for g in genes:
            col = sub[g].to_numpy()
            n = col.size
            frac = float((col > 0).sum()) / n if n else 0.0
            mean_expr = float(col.mean()) if n else 0.0
            records.append({
                group_by: str(grp),
                'gene': g,
                'fraction': frac,
                'mean_expr': mean_expr,
            })

    dot_df = pd.DataFrame.from_records(records)
    # Keep the requested gene order top-to-bottom (reverse so the first gene is
    # at the top of the y-axis).
    gene_order = list(reversed(genes))
    # Order the category (x) axis by the configured cell-type order when available.
    x_order = _resolve_order(dot_df[group_by], category_order=category_order)

    fig = px.scatter(
        dot_df,
        x=group_by,
        y='gene',
        size='fraction',
        color='mean_expr',
        color_continuous_scale='viridis',
        size_max=18,
        category_orders={'gene': gene_order, group_by: x_order},
        labels={
            'fraction': 'Fraction detected',
            'mean_expr': 'Mean expression',
            'gene': 'Gene',
        },
        title=f"Expression dot plot by {group_by}",
    )
    fig.update_layout(
        plot_bgcolor='white',
        height=max(400, 60 * len(genes) + 200),
        margin=dict(t=60, l=120, r=20, b=120),
        xaxis=dict(tickangle=-40),
    )
    return fig


def _dual_colorbar_title(name, shared_scale, is_pct):
    """Colorbar title for a two-gene panel: the per-gene name by default; a single
    neutral title under a shared scale; '(% detected)' when the bin stat is frac_pos."""
    if shared_scale:
        return '% detected' if is_pct else 'expression'
    return f"{name} (% detected)" if is_pct else str(name)


def create_dual_gene_figure(x, y, vals_list, names, embedding, binned=False,
                            bin_size=50, percentile=0.95, smooth_sigma=0,
                            min_cells=1, color_floor=0.05, shared_scale=False,
                            bin_stat='mean'):
    """Two genes side-by-side over the same embedding for visual comparison.

    Each panel shows one gene's expression. When `binned` (custom spatial
    embedding + binning enabled) each panel is a per-bin-mean heatmap computed
    with the shared _binned_mean helper; otherwise each panel is a per-cell
    Scattergl colored by expression. By default each gene gets its OWN colour scale
    (so a sparse gene like CYP26C1 is not flattened by a strong one like FGF8); pass
    shared_scale=True to put both panels on one absolute scale.

    Parameters
    ----------
    x, y : array-like
        Embedding coordinates (same for both panels).
    vals_list : list of array-like
        Per-cell expression for each of the two genes.
    names : list of str
        Gene names (used as subplot titles).
    embedding : str
        Embedding name (axis labels / styling).
    binned : bool
        Heatmap (True) vs per-cell scatter (False).
    """
    if ' vs ' in embedding:
        x_label, y_label = embedding.split(' vs ')
    else:
        x_label = f'{embedding}_1'
        y_label = f'{embedding}_2'

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    fig = make_subplots(rows=1, cols=2, subplot_titles=[str(n) for n in names],
                        horizontal_spacing=0.14)

    # By default each gene gets its OWN colour scale + colorbar, so a weakly
    # expressed gene (e.g. CYP26C1) is not flattened by a strongly expressed one
    # (e.g. FGF8) under a shared maximum -- the reason CYP26C1 "did not show up".
    # shared_scale=True restores a single scale for absolute side-by-side comparison.
    cbx = [0.43, 1.005]  # colorbar x-positions for the left / right panels

    if binned:
        grids = []
        vmaxes = []
        for vals in vals_list:
            H_mean, H_counts, xedges, yedges, vmax = _binned_mean(
                x, y, vals,
                bin_size=bin_size, percentile=percentile,
                smooth_sigma=smooth_sigma, min_cells=min_cells,
                color_floor=color_floor, stat=bin_stat,
            )
            grids.append((H_mean, xedges, yedges))
            vmaxes.append(vmax)
        shared_vmax = max(vmaxes) if vmaxes else float(color_floor)
        is_pct = (bin_stat == 'frac_pos')

        for idx, (H_mean, xedges, yedges) in enumerate(grids):
            x_centers = (xedges[:-1] + xedges[1:]) / 2
            y_centers = (yedges[:-1] + yedges[1:]) / 2
            zmax = shared_vmax if shared_scale else vmaxes[idx]
            fig.add_trace(
                go.Heatmap(
                    z=H_mean.T,
                    x=x_centers,
                    y=y_centers,
                    colorscale='viridis',
                    zmin=0,
                    zmax=zmax,
                    hoverongaps=False,
                    showscale=(idx == 1) if shared_scale else True,
                    colorbar=dict(title=_dual_colorbar_title(names[idx], shared_scale, is_pct),
                                  x=cbx[idx], len=1.0, thickness=12,
                                  # fraction [0,1] -> percentage ticks for the % detected map
                                  tickformat='.0%' if is_pct else None),
                ),
                row=1, col=idx + 1,
            )
    else:
        # Per-cell scatter. Each panel uses its own maximum unless shared_scale.
        per_max = []
        for vals in vals_list:
            arr = np.asarray(vals, dtype=float)
            per_max.append(float(np.nanmax(arr)) if arr.size and np.isfinite(arr).any() else 0.0)
        shared_vmax = max(per_max + [float(color_floor)])

        for idx, vals in enumerate(vals_list):
            cmax = shared_vmax if shared_scale else max(per_max[idx], float(color_floor))
            v = np.asarray(vals, dtype=float)
            # Draw expressing cells LAST (on top) so a sparse gene's few positive
            # cells are not hidden under the many non-expressing cells.
            order = np.argsort(v, kind='stable')
            cb_title = _dual_colorbar_title(names[idx], shared_scale, is_pct=False)
            fig.add_trace(
                go.Scattergl(
                    x=x[order], y=y[order],
                    mode='markers',
                    marker=dict(
                        size=3,
                        color=v[order],
                        colorscale='viridis',
                        cmin=0,
                        cmax=cmax,
                        showscale=(idx == 1) if shared_scale else True,
                        colorbar=dict(title=cb_title, x=cbx[idx], len=1.0, thickness=12),
                    ),
                    showlegend=False,
                    hoverinfo='skip',
                ),
                row=1, col=idx + 1,
            )

    is_umap = 'umap' in embedding.lower()
    for col in (1, 2):
        fig.update_xaxes(
            title_text=x_label,
            showgrid=False, zeroline=False,
            showticklabels=not is_umap,
            row=1, col=col,
        )
        fig.update_yaxes(
            title_text=y_label,
            showgrid=False, zeroline=False,
            showticklabels=not is_umap,
            scaleanchor=('x' if col == 1 else 'x2'),
            scaleratio=1,
            row=1, col=col,
        )

    fig.update_layout(
        plot_bgcolor='white',
        height=600,
        margin=dict(t=60, l=60, r=20, b=60),
    )
    return fig


def create_group_expression_figure(
    df, gene, group_by, *, color_map=None, category_order=None,
    positive_only=True, min_cells_pseudobulk=10, jitter_skip_threshold=2000,
    norm_target=1e4,
):
    """NPY-style two-panel expression figure for one gene across cell groups.

    Plotly reproduction of scripts/viewer_full_chick/16b_npy_figure.R.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per cell. Required column 'expr' = per-cell log1p(CP_norm_target)
        (the h5ad .X). Optional 'ncount' (= nCount_RNA, the per-cell library size)
        and 'replicate' enable Panel A; without them only the per-cell violin panel
        is drawn. Plus the `group_by` categorical column.
    gene, group_by : str
        Gene name (titles/axis); categorical column for the x-axis (e.g. cell type).
    color_map, category_order : dict / list or None
        Fixed {category -> hex} and explicit left-to-right order (e.g. maturation
        order), so each group keeps its dataset-wide colour.
    positive_only : bool
        Panel B shows only cells with expr > 0 (gene-positive). Panel A pseudobulk
        always uses every cell in the group (it is depth-normalised).
    min_cells_pseudobulk : int
        Drop a (group x replicate) pseudobulk dot with fewer than this many cells.
    jitter_skip_threshold : int
        Suppress the per-cell jitter overlay for groups with more positive cells
        than this (keeps the violin readable for the big progenitor classes).
    norm_target : float or None
        The counts-per-X normalisation target of .X (1e4 for CP10K, 1e6 for CPM),
        detected per dataset at load. Panel A reconstructs raw counts as
        raw_i = expm1(expr_i) * ncount_i / norm_target and plots
        log1p( sum(raw)/sum(depth) * norm_target ) -- the SAME log1p(CP) units as
        Panel B. If None (e.g. .X is z-scored / not a clean log1p(CP)), Panel A is
        omitted rather than computed from an unknown normalisation.
    """
    if df is None or len(df) == 0 or group_by not in df.columns:
        return _message_figure("No data to plot for this selection.")
    d = df[df[group_by].notna() & df['expr'].notna()].copy()
    if len(d) == 0:
        return _message_figure("No data to plot for this selection.")
    d[group_by] = d[group_by].astype(str)

    order = _resolve_order(d[group_by], category_order=category_order, by_series=d['expr'])
    pos_index = {g: i for i, g in enumerate(order)}
    n_groups = len(order)

    full_cmap = _resolve_color_map(order, color_map)

    def gcolor(g):
        return full_cmap[g]  # _resolve_color_map covers every category in `order`

    # Panel A (pseudobulk) needs a replicate, a per-cell depth, AND a known CP
    # normalisation target so raw counts can be reconstructed. norm_target carries the
    # detected value (None when .X is not a clean log1p(CP*), e.g. z-scored).
    has_cols = ('ncount' in d.columns and 'replicate' in d.columns
                and d['ncount'].notna().any() and d['replicate'].notna().any())
    draw_panel_a = has_cols and bool(norm_target)
    # Reason Panel A is omitted DESPITE the data carrying a replicate + depth (so the
    # user can tell "not recognised" from "no replicate configured"); surfaced below.
    pb_omit_reason = (".X is not a recognised log1p(CP) normalisation"
                      if (has_cols and not norm_target) else None)
    if draw_panel_a:
        # Verify raw_i = expm1(X)*nCount/norm_target lands on integers for this gene;
        # if .X was normalised on a depth other than nCount_RNA, omit Panel A rather
        # than draw unverified pseudobulk.
        nc = pd.to_numeric(d['ncount'], errors='coerce').to_numpy()
        raw = np.expm1(d['expr'].to_numpy()) * nc / float(norm_target)
        raw_pos = raw[np.isfinite(raw) & (d['expr'].to_numpy() > 0)]
        if raw_pos.size and float(np.max(np.abs(raw_pos - np.round(raw_pos)))) > 1e-2:
            draw_panel_a = False
            pb_omit_reason = ".X does not appear normalised on nCount_RNA"

    rows = 2 if draw_panel_a else 1
    titles = []
    if draw_panel_a:
        titles.append(f"Pseudobulk {gene} per population x replicate "
                      f"(n = replicates passing the >={int(min_cells_pseudobulk)}-cell gate; dot size ~ n cells)")
    pos_label = f"{gene}+ cells only" if positive_only else "all cells"
    titles.append(f"Per-cell {gene} ({pos_label}), by population")

    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True, vertical_spacing=0.10,
                        subplot_titles=titles)
    row_b = 2 if draw_panel_a else 1

    # ---- Panel A: pseudobulk per (group x replicate); dot size grows with n cells ----
    if draw_panel_a:
        cells = d[d['ncount'].notna()].copy()
        cells['raw'] = np.expm1(cells['expr'].to_numpy()) * cells['ncount'].to_numpy() / float(norm_target)
        cells = cells[np.isfinite(cells['raw'])]  # drop any non-finite reconstruction
        grouped = cells.groupby([group_by, 'replicate'], observed=True)
        pb_all = grouped.agg(raw_sum=('raw', 'sum'), depth=('ncount', 'sum'),
                             n=('expr', 'size')).reset_index()
        pb_all = pb_all[pb_all['depth'] > 0].copy()              # usable-depth replicate units
        pb_df = pb_all[pb_all['n'] >= int(min_cells_pseudobulk)].copy()
        # Replicate-dot counts per population: K = usable units, k = units passing the
        # >=min_cells gate. k < K is exactly why populations show different numbers of
        # pseudobulk dots, so both are surfaced as a per-population label below.
        n_total_g = pb_all.groupby(group_by, observed=True).size().to_dict()
        n_shown_g = pb_df.groupby(group_by, observed=True).size().to_dict()
        pb_df['pb'] = np.log1p(pb_df['raw_sum'] / pb_df['depth'] * float(norm_target))
        n_ref = max(1, int(pb_df['n'].max())) if len(pb_df) else 1
        # diameter grows with sqrt(n); +6px floor keeps small replicates visible
        # (so size increases with n but is not strictly area-proportional).
        pb_df['px'] = 6.0 + 18.0 * np.sqrt(pb_df['n'] / n_ref)
        for g in order:
            sub = pb_df[pb_df[group_by] == g]
            if not len(sub):
                continue
            i = pos_index[g]
            jx = i + np.random.RandomState(i + 1).uniform(-0.16, 0.16, len(sub))
            fig.add_trace(go.Scatter(
                x=jx, y=sub['pb'], mode='markers',
                marker=dict(size=sub['px'], color=gcolor(g), line=dict(width=0), opacity=0.8),
                showlegend=False,
                customdata=np.stack([sub['replicate'].astype(str), sub['n']], axis=-1),
                hovertemplate=(f"{g}<br>replicate=%{{customdata[0]}}<br>"
                               f"pseudobulk {gene}=%{{y:.3f}}<br>n cells=%{{customdata[1]}}<extra></extra>"),
            ), row=1, col=1)
            m = float(sub['pb'].mean())
            fig.add_trace(go.Scatter(
                x=[i], y=[m], mode='markers',
                marker=dict(symbol='line-ew', size=24, color='rgba(60,60,60,0.9)',
                            line=dict(width=2.5, color='rgba(60,60,60,0.9)')),
                showlegend=False,
                hovertemplate=f"{g}<br>mean across replicates={m:.3f}<extra></extra>",
            ), row=1, col=1)
        # Per-population replicate-dot count (mirrors Panel B's n=), so the varying number
        # of pseudobulks -- driven by the >=min_cells gate -- is explicit. "n=k" when all
        # usable units pass; "n=k/K" when K-k were gated out (k can be 0).
        if len(pb_df):
            a_top = float(pb_df['pb'].max()); a_bot = float(pb_df['pb'].min())
            a_lab = a_top + 0.12 * max(a_top - a_bot, 0.5)
            for g in order:
                k = int(n_shown_g.get(g, 0)); K = int(n_total_g.get(g, 0))
                label = f"n={k}" if k == K else f"n={k}/{K}"
                fig.add_trace(go.Scatter(
                    x=[pos_index[g]], y=[a_lab], mode='text', text=[label],
                    textfont=dict(size=9, color='grey'), showlegend=False, hoverinfo='skip',
                ), row=1, col=1)
        cp_label = "log1p(CP10K)" if abs(float(norm_target) - 1e4) < 1.0 else f"log1p(CP/{float(norm_target):g})"
        fig.update_yaxes(title_text=f"Pseudobulk {gene}<br>{cp_label}", row=1, col=1)

    # ---- Panel B: per-cell violins of (positive) cells ----
    pos_cells = d[d['expr'] > 0].copy() if positive_only else d.copy()
    if len(pos_cells) == 0 and not draw_panel_a:
        return _message_figure(f"No {gene}-positive cells to plot for this selection.")
    n_pos = pos_cells.groupby(group_by, observed=True)['expr'].size().to_dict() if len(pos_cells) else {}
    y_top = float(pos_cells['expr'].max()) if len(pos_cells) else 1.0
    for g in order:
        i = pos_index[g]
        sub = pos_cells[pos_cells[group_by] == g]
        if len(sub):
            fig.add_trace(go.Violin(
                x=np.full(len(sub), i), y=sub['expr'].to_numpy(),
                width=0.85, scalemode='width', points=False,
                line=dict(color='rgba(50,50,50,0.6)', width=1),
                fillcolor=gcolor(g), opacity=0.55, showlegend=False,
                hovertemplate=f"{g}<br>{gene}=%{{y:.3f}}<extra></extra>",
            ), row=row_b, col=1)
            if len(sub) <= int(jitter_skip_threshold):
                jx = i + np.random.RandomState(1000 + i).uniform(-0.14, 0.14, len(sub))
                fig.add_trace(go.Scatter(
                    x=jx, y=sub['expr'], mode='markers',
                    marker=dict(size=2.5, color='rgba(20,20,20,0.22)'),
                    showlegend=False, hoverinfo='skip',
                ), row=row_b, col=1)
        # Always annotate the positive-cell count (n=0 for groups with none) so a
        # zero-positive group reads as an explicit n=0 rather than a silent gap.
        fig.add_trace(go.Scatter(
            x=[i], y=[y_top * 1.06], mode='text',
            text=[f"n={int(n_pos.get(g, 0))}"],
            textfont=dict(size=10, color='grey'),
            showlegend=False, hoverinfo='skip',
        ), row=row_b, col=1)
    fig.update_yaxes(title_text=f"{gene} (log-norm)<br>{pos_label}", row=row_b, col=1)
    if len(pos_cells) == 0:
        # draw_panel_a is True here (else we returned above): Panel A is informative
        # but no cell is positive -- label the empty violin panel instead of blank.
        xr = 'x domain' if row_b == 1 else f'x{row_b} domain'
        yr = 'y domain' if row_b == 1 else f'y{row_b} domain'
        fig.add_annotation(text=f"No {gene}-positive cells in any population",
                           xref=xr, yref=yr, x=0.5, y=0.5, showarrow=False,
                           font=dict(size=14, color='#888'))

    if not draw_panel_a and pb_omit_reason:
        # The default "figure" style promises a pseudobulk panel; when it is omitted
        # for a recognised reason (the data had a replicate but .X is not a clean CP
        # normalisation), say so rather than silently showing a single panel.
        fig.add_annotation(text=f"Pseudobulk panel omitted: {pb_omit_reason}",
                           xref='paper', yref='paper', x=0.5, y=1.0, yanchor='bottom',
                           showarrow=False, font=dict(size=11, color='#a06000'))

    for r in range(1, rows + 1):
        fig.update_xaxes(
            tickmode='array', tickvals=list(range(n_groups)), ticktext=order,
            tickangle=-40, range=[-0.6, n_groups - 0.4], showgrid=False, row=r, col=1,
        )
    fig.update_layout(
        plot_bgcolor='white',
        violinmode='overlay',
        height=780 if draw_panel_a else 480,
        margin=dict(t=60, l=80, r=30, b=140),
    )
    return fig