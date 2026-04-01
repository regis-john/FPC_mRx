"""
================================================================================
12 Panel Plotting functions for MMS VDF and FPC data.

Functions:
- imagecont_vdf_fpc: Generate 12 panel subplots of VDF and FPC components.
- imagecont_vdf_fpc_panel: Helper function for imagecont_vdf_fpc.

Author: Regis John
Created: 2026-03-31
================================================================================
"""

import h5py
from pyspedas import time_string
import matplotlib.pyplot as plt
import os

from src.utils.helper_utils import (mms_name_make, slice3d_to_2d, 
                                    global_vmin_vmax_vdf, global_vmin_vmax_fpc)
from src.plotting.imagecont_vdf import imagecont_vdf
from src.plotting.imagecont_fpc import imagecont_fpc


def imagecont_vdf_fpc(trange, species, bin_width_frac, time_index=0,
                          coord_type='both', suptitle_sfx=None, dpi_val=100, 
                          axis_fntsz=16, cbar_ht=0.050, save_fig=False, 
                          outdir_sfx='', show_tind=True, **kwargs):
    """
    Master function to generate 12 panel subplots of VDF and FPC components for 
    a given time index with the VDFs on the top row and the FPCs on the bottom 
    three rows.

    Parameters:
    - trange (list of str): start time, end time] in the format:
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss'].
    - species (str): species name ('e' or 'i').
    - bin_width_frac (float): bin width as a fraction of thermal velocity.
    - time_index (int): time index for the plot.
    - coord_type (str): coordinate system to use ('fac', 'lmn', or 'both').
    - suptitle_sfx (str or None): suffix for the suptitle.
    - dpi_val (int): DPI value for the plot.
    - axis_fntsz (int): font size for the axis labels.
    - cbar_ht (float): height of the colorbar.
    - save_fig (bool): if True, save the figure.
    - outdir_sfx (str): suffix for the output directory.
    - show_tind (bool): if True, show the time index in the plot at the top right corner.

    Returns:
    - fig (matplotlib.figure.Figure): The figure object.
    - ax (matplotlib.axes.Axes): The axes object.
    """
    if suptitle_sfx is None:
        suptitle_sfx = f"with {bin_width_frac} binwidth"
    
    data_dir = "data"
    hfile_pref = f'mms_fpc_{species}_{bin_width_frac:.2f}'
    hfile = os.path.join(data_dir, mms_name_make(hfile_pref, trange[0], trange[1]))
    # hfile = mms_name_make(hfile_pref, trange[0], trange[1])
    
    with h5py.File(hfile, "r") as f:
        time = f["meta"]["time"][:]
        species = f["meta"]["species"].asstr()[()]
        probe = f["meta"]["probe"].asstr()[()]
        title_time = time_string(time[time_index])[:-3]
        species_name = {'e':'Electron', 'i':'Ion'}

        if coord_type in ("fac", "both"):
            suptitle = (f"{species_name[species]} VDF and FPC components (FAC) " 
                        f"{suptitle_sfx}\n MMS{probe} — {title_time}")
            vdf_fac = f["fac"]["vdf_vol"][:, :, :, time_index]
            c_fac = f["fac"]["c_fac_vol"][:, :, :, :, time_index]
            xc_fac = f["fac"]["binc_facn"][:]
            yc_fac = f["fac"]["binc_facn"][:]
            zc_fac = f["fac"]["binc_facn"][:]
            je_fac = f["fac"]["JE_tot"][:, time_index]
            vdf_xy, vdf_yz, vdf_xz = slice3d_to_2d(vdf_fac)
            cx_xy, cx_yz, cx_xz = slice3d_to_2d(c_fac[0])
            cy_xy, cy_yz, cy_xz = slice3d_to_2d(c_fac[1])
            cz_xy, cz_yz, cz_xz = slice3d_to_2d(c_fac[2])
            vdf_2d_fac = (vdf_xy, vdf_xz, vdf_yz.T) # The order is xy, xz and zy
            cx_2d_fac = (cx_xy, cx_xz, cx_yz.T)
            cy_2d_fac = (cy_xy, cy_xz, cy_yz.T)
            cz_2d_fac = (cz_xy, cz_xz, cz_yz.T)
            ax_pairs_fac = ((xc_fac, yc_fac), (xc_fac, zc_fac), (zc_fac, yc_fac))
            axlabels_fac = (r"$v_{\perp 1}/v_{te}$", r"$v_{\perp 2}/v_{te}$", 
                                r"$v_{\parallel}/v_{te}$")
            axlabel_pairs_fac = ((axlabels_fac[0], axlabels_fac[1]), # xy, xz, zy
                                (axlabels_fac[0], axlabels_fac[2]), 
                                (axlabels_fac[2], axlabels_fac[1]))
            row_labels_fac = (rf"$\mathrm{{log}}(f_{{{species}}})$", 
                    r"$C_{\mathrm{E}\perp 1}$", r"$C_{\mathrm{E}\perp 2}$", 
                    r"$C_{\mathrm{E}\parallel}$")
            diag_labels_fac = (rf"$J_{{\perp 1}}E_{{\perp 1}} = {je_fac[0]:.2e}$", 
                               rf"$J_{{\perp 2}}E_{{\perp 2}} = {je_fac[1]:.2e}$", 
                               rf"$J_{{\parallel}}E_{{\parallel}} = {je_fac[2]:.2e}$")
            fig = imagecont_vdf_fpc_panel(vdf_2d_fac, cx_2d_fac, cy_2d_fac, cz_2d_fac, 
                        ax_pairs_fac, axlabel_pairs_fac, row_labels_fac, diag_labels_fac, 
                        dpi_val, axis_fntsz, cbar_ht)
            fig.suptitle(suptitle, fontsize = 24)
            
            if show_tind:
                fig.text(0.98, 0.98, f"t: {time_index}", ha="right", va="top", fontsize=18)

            if save_fig:
                plot_name = f"mms_fpc_{species}_fac_{bin_width_frac:.2f}_{time_index}.png"
                sfx_str = f'_{outdir_sfx}' if outdir_sfx else ''
                output_dir = f'plots_fac{sfx_str}'
                os.makedirs(output_dir, exist_ok=True)
                fig.savefig(os.path.join(output_dir, plot_name), dpi=dpi_val, bbox_inches='tight')
                plt.close(fig)
            else:
                plt.show()


        if coord_type in ("lmn", "both"):
            suptitle_sfx
            suptitle = (f"{species_name[species]} VDF and FPC components (LMN) " 
                            f"{suptitle_sfx}\n MMS{probe} — {title_time}")
            vdf_lmn = f["lmn"]["vdf_vol"][:, :, :, time_index]
            c_lmn = f["lmn"]["c_lmn_vol"][:, :, :, :, time_index]
            yc_lmn = f["lmn"]["binc_lmnn"][:]
            xc_lmn = f["lmn"]["binc_lmnn"][:]
            zc_lmn = f["lmn"]["binc_lmnn"][:]
            je_lmn = f["lmn"]["JE_tot"][:, time_index]
            vdf_xy, vdf_yz, vdf_xz = slice3d_to_2d(vdf_lmn)
            cx_xy, cx_yz, cx_xz  = slice3d_to_2d(c_lmn[0])
            cy_xy, cy_yz, cy_xz  = slice3d_to_2d(c_lmn[1])
            cz_xy, cz_yz, cz_xz  = slice3d_to_2d(c_lmn[2])
            vdf_2d_lmn = (vdf_xy, vdf_xz, vdf_yz.T) # The order is xy, xz and zy
            cx_2d_lmn = (cx_xy, cx_xz, cx_yz.T)
            cy_2d_lmn = (cy_xy, cy_xz, cy_yz.T)
            cz_2d_lmn = (cz_xy, cz_xz, cz_yz.T)
            ax_pairs_lmn = ((xc_lmn, yc_lmn), (xc_lmn, zc_lmn), (zc_lmn, yc_lmn))
            axlabels_lmn = (r"$v_L$", r"$v_M$", r"$v_N$")
            axlabel_pairs_lmn = ((axlabels_lmn[0], axlabels_lmn[1]), 
                                (axlabels_lmn[0], axlabels_lmn[2]), 
                                (axlabels_lmn[2], axlabels_lmn[1]))
            row_labels_lmn = (rf"$\mathrm{{log}}(f_{{{species}}})$", 
                    r"$C_{\mathrm{E}_L}$", r"$C_{\mathrm{E}_M}$", 
                    r"$C_{\mathrm{E}_N}$")
            diag_labels_lmn = (rf"$J_{{L}}E_{{L}} = {je_lmn[0]:.2e}$", 
                               rf"$J_{{M}}E_{{M}} = {je_lmn[1]:.2e}$", 
                               rf"$J_{{N}}E_{{N}} = {je_lmn[2]:.2e}$")
            fig = imagecont_vdf_fpc_panel(vdf_2d_lmn, cx_2d_lmn, cy_2d_lmn, cz_2d_lmn, 
                        ax_pairs_lmn, axlabel_pairs_lmn, row_labels_lmn, diag_labels_lmn, 
                        dpi_val, axis_fntsz, cbar_ht)
            fig.suptitle(suptitle, fontsize = 24)

            if show_tind:
                fig.text(0.98, 0.98, f"t: {time_index}", ha="right", va="top", fontsize=18)

            if save_fig:
                plot_name = f"mms_fpc_{species}_lmn_{bin_width_frac:.2f}_{time_index}.png"
                sfx_str = f'_{outdir_sfx}' if outdir_sfx else ''
                output_dir = f'plots_lmn{sfx_str}'
                os.makedirs(output_dir, exist_ok=True)
                fig.savefig(os.path.join(output_dir, plot_name), dpi=dpi_val, bbox_inches='tight')
                plt.close(fig)
            else:
                plt.show()


def imagecont_vdf_fpc_panel(vdf_2d, cx_2d, cy_2d, cz_2d, ax_pairs, axlabels_pairs, 
                            row_labels, diag_labels, dpi_val=100, axis_fntsz=16, 
                            cbar_ht=0.047):
    """
    Helper function to the above imagecont_vdf_fpc_panel function which plots 
    a 4x3 panel of VDF and FPC components.

    Parameters:
    - vdf_2d (tuple of length 3 of 2D arrays): The VDF data.
    - cx_2d (tuple of length 3 of 2D arrays): The Cx FPC component data.
    - cy_2d (tuple of length 3 of 2D arrays): The Cy FPC component data.
    - cz_2d (tuple of length 3 of 2D arrays): The Cz FPC component data.
    - ax_pairs (tuple of tuples of length 3 of 2-tuples of 1D arrays): The pairs of axes for the plots.
    - axlabels_pairs (tuple of tuples of length 3 of 2-tuples of strings): The labels for the axes.
    - row_labels (tuple of strings): The labels for the rows.
    - diag_labels (tuple of strings): The labels for the diagonal plots.
    - dpi_val (int): The DPI for the figure.
    - axis_fntsz (int): The fontsize for the axis labels.
    - cbar_ht (float): The height of the colorbars.

    Returns:
    - fig (matplotlib.figure.Figure): The figure object.
    """
    fig, axes = plt.subplots(4, 3, figsize=(18, 18), dpi=dpi_val, constrained_layout=True)

    # --- Row 1: VDF top row ---
    for i in range(3):
        vmin, vmax = global_vmin_vmax_vdf(vdf_2d[i])
        ax1, ax2 = ax_pairs[i]
        _, ylabel = axlabels_pairs[i]
        im1, _ = imagecont_vdf(ax1, ax2, vdf_2d[i], xlabel=None, ylabel=ylabel, 
                           ax=axes[0,i], vmin=vmin, vmax=vmax, axis_fntsz=axis_fntsz)
        height, width = len(ax2), len(ax1)
        cbar = fig.colorbar(im1, ax=axes[0,i], fraction=cbar_ht*height/width, pad=0.08)
        cbar.set_label(r"[m$^{-3}$]")
        axes[0,0].text(-0.4, 0.85, row_labels[0], transform=axes[0,0].transAxes, 
                va='center', ha='left', rotation='horizontal', fontsize=20)


    # --- Row 2: Cx ---
    for j in range(3):
        vmin, vmax = global_vmin_vmax_fpc(cx_2d[j])
        ax1, ax2 = ax_pairs[j]
        _, ylabel = axlabels_pairs[j]
        im2, _ = imagecont_fpc(ax1, ax2, cx_2d[j], xlabel=None, ylabel=ylabel, 
                           ax=axes[1,j], vmin=vmin, vmax=vmax, axis_fntsz=axis_fntsz)
        height, width = len(ax2), len(ax1)
        cbar = fig.colorbar(im2, ax=axes[1,j], fraction=cbar_ht*height/width, pad=0.005)
        cbar.set_label(r"[W m$^{-3}$]")
        axes[1,0].text(-0.4, 0.85, row_labels[1], transform=axes[1,0].transAxes, 
                va='center', ha='left', rotation='horizontal', fontsize=20)
        axes[1,2].text(0.05, 0.95, diag_labels[0], transform=axes[1,2].transAxes, 
                       ha="left", va="top", fontsize=14, color="black")
        

    # --- Row 3: Cy ---
    for k in range(3):
        vmin, vmax = global_vmin_vmax_fpc(cy_2d[k])
        ax1, ax2 = ax_pairs[k]
        _, ylabel = axlabels_pairs[k]
        im3, _ = imagecont_fpc(ax1, ax2, cy_2d[k], xlabel=None, ylabel=ylabel, 
                           ax=axes[2,k], vmin=vmin, vmax=vmax, axis_fntsz=axis_fntsz)
        height, width = len(ax2), len(ax1)
        cbar = fig.colorbar(im3, ax=axes[2,k], fraction=cbar_ht*height/width, pad=0.005)
        cbar.set_label(r"[W m$^{-3}$]")
        axes[2,0].text(-0.4, 0.85, row_labels[2], transform=axes[2,0].transAxes, 
                va='center', ha='left', rotation='horizontal', fontsize=20)
        axes[2,1].text(0.05, 0.95, diag_labels[1], transform=axes[2,1].transAxes, 
                       ha="left", va="top", fontsize=14, color="black")

    # --- Row 4: Cz ---
    for l in range(3):
        vmin, vmax = global_vmin_vmax_fpc(cz_2d[l])
        ax1, ax2 = ax_pairs[l]
        xlabel, ylabel = axlabels_pairs[l]
        im4, _ = imagecont_fpc(ax1, ax2, cz_2d[l], xlabel=xlabel, ylabel=ylabel, 
                           ax=axes[3,l], vmin=vmin, vmax=vmax, axis_fntsz=axis_fntsz)
        height, width = len(ax2), len(ax1)
        cbar = fig.colorbar(im4, ax=axes[3,l], fraction=cbar_ht*height/width, pad=0.005)
        cbar.set_label(r"[W m$^{-3}$]")
        axes[3,0].text(-0.4, 0.85, row_labels[3], transform=axes[3,0].transAxes, 
                va='center', ha='left', rotation='horizontal', fontsize=20)
        axes[3,0].text(0.05, 0.95, diag_labels[2], transform=axes[3,0].transAxes, 
                       ha="left", va="top", fontsize=14, color="black")
        
    return fig