import plotly.express as px
import pandas as pd
import numpy as np

def create_scatter_plot(df, embedding, color_by, treat_as_categorical=False):
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
    
    # Add custom styling for UMAP plots
    if 'umap' in embedding.lower():
        fig.update_layout(
            plot_bgcolor='white',
            xaxis=dict(
                showgrid=False,
                showticklabels=False,
                scaleanchor="y",
                scaleratio=1,
            ),
            yaxis=dict(
                showgrid=False,
                showticklabels=False,
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
    
    return fig

def create_metacell_plot(H, xedges, yedges, embedding):
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    
    return px.imshow(
        H.T,
        x=x_centers,
        y=y_centers,
        labels={'x': f'{embedding}_1', 'y': f'{embedding}_2'},
        title=f'Metacell visualization - {embedding}',
        aspect='equal'
    ) 