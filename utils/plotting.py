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

def create_scatter_plot(df, embedding, color_by, treat_as_categorical=False, selection_data=None):
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
    
    fig = px.scatter(
        df, x='x', y='y', color='color',
        labels={'x': x_label, 'y': y_label},
        title=f'Single-cell visualization - {embedding}',
        color_discrete_sequence=px.colors.qualitative.Set3 if treat_as_categorical else None,
        color_continuous_scale='viridis' if not treat_as_categorical else None,
        hover_data=None
    )
    
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
    
    # Update selection styling if provided
    if selection_data and selection_data.get('indices'):
        fig.update_traces(
            selectedpoints=selection_data['indices'],
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

def create_binned_plot(df, embedding, color_by, bin_size=50, percentile=0.95, treat_as_categorical=False, smooth_sigma=0, min_cells=1, color_floor=0.05):
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
        )

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
            labels={'x': x_label, 'y': y_label, 'color': color_by},
            title=f'Binned visualization - {embedding}',
            color_continuous_scale='viridis',
            zmin=0,
            zmax=color_vmax,
            aspect='equal'
        )
        
        # Create hover data with bin indices and cell counts
        hover_text = [[
            f'Bin: ({i},{j})<br>'
            f'x: {x_centers[i]:.2f}<br>'
            f'y: {y_centers[j]:.2f}<br>'
            f'Value: {H[i,j]:.2f}<br>'
            f'Cells: {int(H_counts[i,j])}'
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


def _binned_mean(x, y, vals, bin_size=50, percentile=0.95, smooth_sigma=0, min_cells=1, color_floor=0.05):
    """Compute the per-bin mean of a continuous value over a 2D (x, y) grid.

    This is the shared spatial-binning math used by both create_binned_plot's
    continuous branch and create_dual_gene_figure, so the smoothing, masking and
    color-scale clipping stay identical across the single- and dual-gene views.

    Returns
    -------
    H_mean : 2D float array
        Per-bin mean (NaN for bins below the cell floor), clipped to vmax.
    H_counts : 2D float array
        Cells per bin.
    xedges, yedges : 1D arrays
        Bin edges from np.histogram2d.
    vmax : float
        Color-scale maximum (>= color_floor).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    vals = np.asarray(vals, dtype=float)

    H, xedges, yedges = np.histogram2d(x, y, bins=bin_size, weights=vals)
    H_counts, _, _ = np.histogram2d(x, y, bins=bin_size)

    # Bins with at least `min_cells` cells are valid; the rest are greyed.
    valid = H_counts >= max(1, int(min_cells))

    # Mean expression per valid bin (NaN elsewhere).
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
    # which clips over expressing pixels rather than all bins). The floor keeps
    # near-empty genes from being auto-stretched by the smoother into spurious
    # "signal".
    pos = H_mean[valid & np.isfinite(H_mean) & (H_mean > 0)]
    if pos.size and percentile < 1.0:
        vmax = float(np.nanpercentile(pos, percentile * 100))
    elif pos.size:
        vmax = float(np.nanmax(pos))
    else:
        vmax = 0.0
    vmax = max(vmax, float(color_floor))
    H_mean = np.where(H_mean > vmax, vmax, H_mean)

    return H_mean, H_counts, xedges, yedges, vmax


def create_group_expression_plot(df, gene, group_by, split_by=None, style='violin'):
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

    # Order groups by descending median expression so the most-expressing
    # categories read left-to-right.
    order = (
        df.groupby(group_by, observed=True)['expr']
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )
    order = [str(c) for c in order]

    # Stringify the categorical columns so plotly keeps the explicit order and
    # treats them as discrete.
    plot_df = df.copy()
    plot_df[group_by] = plot_df[group_by].astype(str)
    color_arg = None
    if split_by and split_by in plot_df.columns:
        plot_df[split_by] = plot_df[split_by].astype(str)
        color_arg = split_by

    common = dict(
        x=group_by,
        y='expr',
        color=color_arg,
        category_orders={group_by: order},
        title=f"{gene} expression by {group_by}",
        labels={'expr': f"{gene} (log-normalized)", group_by: group_by},
    )

    if style == 'box':
        fig = px.box(plot_df, points=False, **common)
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


def create_dotplot(expr_df, genes, group_by):
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

    fig = px.scatter(
        dot_df,
        x=group_by,
        y='gene',
        size='fraction',
        color='mean_expr',
        color_continuous_scale='viridis',
        size_max=18,
        category_orders={'gene': gene_order},
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


def create_dual_gene_figure(x, y, vals_list, names, embedding, binned=False,
                            bin_size=50, percentile=0.95, smooth_sigma=0,
                            min_cells=1, color_floor=0.05):
    """Two genes side-by-side over the same embedding for visual comparison.

    Each panel shows one gene's expression. When `binned` (custom spatial
    embedding + binning enabled) each panel is a per-bin-mean heatmap computed
    with the shared _binned_mean helper; otherwise each panel is a per-cell
    Scattergl colored by expression. The two panels share one color scale (the
    combined max across both genes) so the comparison is fair.

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
                        horizontal_spacing=0.08)

    if binned:
        # Compute both bin grids first so the shared color scale spans both.
        grids = []
        vmaxes = []
        for vals in vals_list:
            H_mean, H_counts, xedges, yedges, vmax = _binned_mean(
                x, y, vals,
                bin_size=bin_size, percentile=percentile,
                smooth_sigma=smooth_sigma, min_cells=min_cells,
                color_floor=color_floor,
            )
            grids.append((H_mean, xedges, yedges))
            vmaxes.append(vmax)
        shared_vmax = max(vmaxes) if vmaxes else float(color_floor)

        for idx, (H_mean, xedges, yedges) in enumerate(grids):
            x_centers = (xedges[:-1] + xedges[1:]) / 2
            y_centers = (yedges[:-1] + yedges[1:]) / 2
            fig.add_trace(
                go.Heatmap(
                    z=H_mean.T,
                    x=x_centers,
                    y=y_centers,
                    colorscale='viridis',
                    zmin=0,
                    zmax=shared_vmax,
                    hoverongaps=False,
                    showscale=(idx == 1),  # one shared colorbar on the right
                    colorbar=dict(title='expr') if idx == 1 else None,
                ),
                row=1, col=idx + 1,
            )
    else:
        # Per-cell scatter; share the color scale across panels.
        finite_max = 0.0
        for vals in vals_list:
            arr = np.asarray(vals, dtype=float)
            if arr.size and np.isfinite(arr).any():
                finite_max = max(finite_max, float(np.nanmax(arr)))
        shared_vmax = max(finite_max, float(color_floor))

        for idx, vals in enumerate(vals_list):
            fig.add_trace(
                go.Scattergl(
                    x=x, y=y,
                    mode='markers',
                    marker=dict(
                        size=3,
                        color=np.asarray(vals, dtype=float),
                        colorscale='viridis',
                        cmin=0,
                        cmax=shared_vmax,
                        showscale=(idx == 1),
                        colorbar=dict(title='expr') if idx == 1 else None,
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