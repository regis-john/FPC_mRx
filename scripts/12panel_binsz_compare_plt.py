"""
================================================================================
Runnable script to save a series of 12 panel subplots of VDF and FPC components 
for different bin widths.

Author: Regis John
Created: 2026-03-31
================================================================================
"""

from src.plotting.imagecont_vdf_fpc import imagecont_vdf_fpc

trange = ['2015-12-09/05:03:55', '2015-12-09/05:03:59']
binsz = [0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
for bin_width_frac in binsz:
    imagecont_vdf_fpc(trange, species='e', bin_width_frac = bin_width_frac, 
                    time_index=60, coord_type='both', save_fig=True, show_tind=False)
    print(f'Files saved for bin_width_frac = {bin_width_frac}')