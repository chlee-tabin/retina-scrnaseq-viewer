import numpy as np
from scipy.stats import binned_statistic_2d
from scipy.ndimage import gaussian_filter

def create_metacells(x, y, values=None, bins=50, smoothing=1.0):
    """Create metacells from single-cell data"""
    if values is None:
        # If no values provided, count cells per bin
        H, xedges, yedges = np.histogram2d(x, y, bins=bins)
    else:
        # If values provided, calculate mean per bin
        H, xedges, yedges, _ = binned_statistic_2d(
            x, y, values, 
            statistic='mean', 
            bins=bins
        )
    
    # Apply Gaussian smoothing if requested
    if smoothing > 0:
        H = gaussian_filter(H, sigma=smoothing)
    
    return H, xedges, yedges

def calculate_selection_stats(adata, indices):
    """Calculate statistics for selected cells"""
    if not indices:
        return {}
    
    selected_cells = adata[indices]
    stats = {
        'n_cells': len(indices),
        'mean_counts': np.mean(selected_cells.X.sum(axis=1)),
        'mean_genes': np.mean((selected_cells.X > 0).sum(axis=1))
    }
    
    # Add categorical variable statistics
    for col in selected_cells.obs.select_dtypes(include=['category']).columns:
        value_counts = selected_cells.obs[col].value_counts()
        stats[f'{col}_distribution'] = value_counts.to_dict()
    
    return stats 