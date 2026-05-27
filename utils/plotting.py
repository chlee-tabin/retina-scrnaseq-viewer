import plotly.express as px
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
from plotly.subplots import make_subplots
import plotly.graph_objects as go

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
        # For continuous values (including gene expression)
        H, xedges, yedges = np.histogram2d(
            df['x'], df['y'], 
            bins=bin_size,
            weights=df['color'] if 'color' in df else None
        )
        H_counts, _, _ = np.histogram2d(df['x'], df['y'], bins=bin_size)
        
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
        # which clips over expressing pixels rather than all bins). The floor
        # keeps genes with little or no expression from being auto-stretched by
        # the smoother into spurious "signal": a near-zero map (e.g. a gene in a
        # handful of cells) renders near-blank, consistent with the raw view,
        # instead of a bright blob.
        pos = H_mean[valid & np.isfinite(H_mean) & (H_mean > 0)]
        if pos.size and percentile < 1.0:
            vmax = float(np.nanpercentile(pos, percentile * 100))
        elif pos.size:
            vmax = float(np.nanmax(pos))
        else:
            vmax = 0.0
        vmax = max(vmax, float(color_floor))
        H_mean = np.where(H_mean > vmax, vmax, H_mean)

        H = H_mean
        color_vmax = vmax

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