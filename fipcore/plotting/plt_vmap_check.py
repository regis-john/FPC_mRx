"""
================================================================================
Plotting functions for MMS VDF data.

Functions:
- plt_vmap_check: Validate the binning of velocity vectors into a velocity map.

Author: Regis John
Created: 2026-03-26
================================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plt_vmap_check(vv, vmap, edges_norm, vth_mean, comp_x = 2, comp_y = 0, 
                    sweep_idx=0, xlabel='v_x', ylabel='v_y', coord_str='FAC',
                    xlim=1.5, ylim=1.5):
    """
    Validate the binning of velocity vectors into a velocity map.

    Plots a scatter plot of velocity vectors with x and y components
    colored by their vmap bin indices. Overlaid on the scatter plot
    is a grid of hollow rectangles representing the bin boundaries.

    Parameters:
    - vv: np.array of shape (3, N, T); velocity vectors
    - vmap: np.array of shape (3, N, T); vmap bin indices
    - edges_norm: np.array of shape (n_vbins); normalized bin edges
    - vth_mean: float; mean thermal velocity for normalization
    - comp_x, comp_y: int; components of velocity vectors to plot
    - sweep_idx: int; index of sweep to plot
    - xlabel, ylabel: str; labels for x and y axes
    - coord_str: str; string describing coordinate system
    - xlim, ylim: float; limits of x and y axes
    """
    # Extract velocities and normalize by vth_mean
    v_y_norm = vv[comp_y, :, sweep_idx] / vth_mean
    v_x_norm = vv[comp_x, :, sweep_idx] / vth_mean
    
    # Extract vmap indices
    m_y = vmap[comp_y, :, sweep_idx]
    m_x  = vmap[comp_x, :, sweep_idx]
    
    # Filter for valid binned points only
    valid = (m_y >= 0) & (m_x >= 0)
    vp_norm_val = v_y_norm[valid]
    vz_norm_val = v_x_norm[valid]
    m_p_val = m_y[valid]
    m_z_val = m_x[valid]
    
    # Normalize Grid Edges
    dv_norm = edges_norm[1] - edges_norm[0]
    
    # Plotting
    fig, ax = plt.subplots(figsize=(8, 8))
    color_map = plt.get_cmap('tab20', 20)    # Get 20 discrete tab20 colors
    
    # Identify unique occupied bins to avoid redundant drawing
    unique_bins = np.unique(np.column_stack((m_z_val, m_p_val)), axis=0)
    
    for b_z, b_p in unique_bins:
        # Use modulo to cycle through 20 colors based on bin indices
        color_idx = (b_z + b_p) % 20
        edge_color = color_map(color_idx)
        
        # Draw a hollow rectangle representing the bin boundary
        # edges_norm[idx] is the left/bottom edge
        rect = patches.Rectangle(
            (edges_norm[b_z], edges_norm[b_p]), dv_norm, dv_norm,
            linewidth=1.0, edgecolor=edge_color, fill=False, alpha=0.8, zorder=2
        )
        ax.add_patch(rect)
        
    # Scatter the points using same color logic as the boxes
    point_colors = [(b_z + b_p) % 20 for b_z, b_p in zip(m_z_val, m_p_val)]
    ax.scatter(vz_norm_val, vp_norm_val, c=point_colors, cmap='tab20', s=8, zorder=3, edgecolors='none')
    
    # Formatting
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f'Velocity Map Binning Validation (Sweep {sweep_idx}) in {coord_str}')
    
    # Limits and Aspect
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()