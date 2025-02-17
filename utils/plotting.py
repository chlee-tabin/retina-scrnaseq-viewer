import plotly.express as px
import pandas as pd
import numpy as np

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