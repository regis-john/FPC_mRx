"""
================================================================================
Runnable script to save a series of 12 panel subplots of VDF and FPC components 
for different time indices.

Author: Regis John
Created: 2026-03-31
================================================================================
"""

from src.plotting.imagecont_vdf_fpc import imagecont_vdf_fpc

trange = ['2015-12-09/05:03:55', '2015-12-09/05:03:59']
t_list = [58, 59, 60, 61, 62, 63, 64]
for t_ind in t_list:
    imagecont_vdf_fpc(trange, species='e', bin_width_frac=0.25, time_index = t_ind,
            coord_type='both', suptitle_sfx='', dpi_val=100, axis_fntsz=16, 
            cbar_ht=0.050, save_fig=True, outdir_sfx='', show_tind=True)
    print(f'Files saved for time_index = {t_ind}')