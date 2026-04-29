"""
================================================================================
Runnable script to fold the fpc components and save a series of 12 panel 
subplots of VDF and folded FPC components.

Author: Regis John
Created: 2026-04-24
================================================================================
"""

import numpy as np
import os
import h5py
from fipcore import mms_name_make

def compute_fpc_fold(trange, species, bin_width_frac, time_index=0, 
                     coord_type="fac"):

    # Read in the computed fpc:
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_dir = os.path.join(project_root, "data")
    print(data_dir)
    hfile_pref = f'mms_fpc_{species}_{bin_width_frac:.2f}'
    hfile = os.path.join(data_dir, mms_name_make(hfile_pref, trange[0], trange[1]))
    print(hfile)
    # with h5py.File(hfile, "r") as f:
    #     fpc_dat = f["fac"]["c_fac_vol"][:, :, :, :, time_index]
    #     print(fpc_dat.shape)
    # pass