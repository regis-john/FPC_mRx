"""
================================================================================
Functions for processing MMS VDF data.

Functions:
- process_vdf: Process the raw velocity distribution function (VDF) data.
- bin_vdf_3d: Bin VDF data into 3D velocity space.
- bin_vdf_vol: Bin VDF data into 3D velocity space with volume weighting.
- bin_vdf_avg: Bin VDF data into 3D velocity space with averaging.

Author: Regis John
Created: 2026-03-26
================================================================================
"""

import numpy as np
from pyspedas import get_data
from scipy.interpolate import interp1d


def process_vdf(vdf_raw_var, vdf_err_var, dq_flags_var, vvol):
    """
    Process the raw velocity distribution function (VDF) data.

    Parameters
    ----------
    vdf_raw_var : str
        tplot variable of raw VDF data (time, phi, theta, energy)
    vdf_err_var : str
        tplot variable of error in raw VDF data (time, phi, theta, energy)
    dq_flags_var : str
        tplot variable of data quality flags (time, phi, theta, energy)
    vvol : numpy.ndarray
        Volume element (16384, n_sweeps) for weighting

    Returns
    -------
    vdf_raw : numpy.ndarray
        Processed VDF data without volume element weighting in s^3/m^6 (bins, time)
    vdf_vol : numpy.ndarray
        Processed VDF data with volume weighting in m^-3 (bins, time)

    """
    # Get Data
    vdf_rawi = get_data(vdf_raw_var).y # order (time, phi, theta, energy)
    time_x = get_data(vdf_raw_var).times
    vdf_err = get_data(vdf_err_var).y
    dq_flags = get_data(dq_flags_var).y
    
    # Flatten to (ntime, 16384)
    nbin = 32*16*32
    ntime = len(time_x)
    vdf_rawi = vdf_rawi.reshape(ntime, nbin)
    vdf_err = vdf_err.reshape(ntime, nbin)
    
    # Quality Check & Temporal Interpolation
    good_idx = np.where(dq_flags == 0)[0]
    bad_idx = np.where(dq_flags != 0)[0]
    
    if len(bad_idx) > 0 and len(good_idx) > 1:
        # Interpolate every bin across the time dimension
        itp = interp1d(time_x[good_idx], vdf_rawi[good_idx, :], axis=0, 
                       kind='linear', fill_value="extrapolate")
        vdf_rawi[bad_idx, :] = itp(time_x[bad_idx])
        
        # Also interpolate error for the noise calculation
        itp_err = interp1d(time_x[good_idx], vdf_err[good_idx, :], axis=0, 
                           kind='linear', fill_value="extrapolate")
        vdf_err[bad_idx, :] = itp_err(time_x[bad_idx])

    # Counts = round((vdf/sigma)^2)
    with np.errstate(divide='ignore', invalid='ignore'): # Ignore divide by 0, keep console clean
        counts = np.round(np.square(vdf_rawi / vdf_err))

        # Set 0-counts to 0
        both_zero = (vdf_rawi == 0) & (vdf_err == 0)
        counts[both_zero] = 0

        # Replace 0-count or negative vdf with NaN
        mask = (counts <= 0) | (vdf_rawi <= 0) | np.isnan(vdf_rawi)
        vdf_rawi[mask] = np.nan

    # Volume Weighting (Interleaved)
    # f is (ntime, 16384), vvol is (16384, n_sweeps)
    vdf_weighted = np.zeros_like(vdf_rawi)

    ntime = vdf_rawi.shape[0]
    n_sweeps = vvol.shape[1]   # 1 or 2 depending on interleaving
    # vvol is (16384, 1) for both non-interleaved and instantaneous
    ip = np.arange(ntime) % n_sweeps
    vdf_weighted = vdf_rawi * vvol[:, ip].T

    # Scale to SI units (s^3/m^6)
    vdf_vol = (vdf_weighted * 1e12).T
    vdf_raw = (vdf_rawi * 1e12).T


    # Final Output: (bins, time)
    return vdf_raw, vdf_vol


def bin_vdf_3d(vdf_dat, vmap, n_vx, n_vy, n_vz, interleave_check=True):
    """
    Bin VDF data into a 3D grid.

    Parameters
    ----------
    vdf_dat : numpy.ndarray
        Processed VDF data(bins, time)
    vmap : numpy.ndarray
        Mapping of velocity bins in normalized velocity space (3, n_vx*n_vy*n_vz, n_sweeps)
    n_vx : int
        Number of bins along 1st direction
    n_vy : int
        Number of bins along 2nd direction
    n_vz : int
        Number of bins along 3rd direction
    interleave_check : bool
        Flag to check if the data is interleaved, default=True

    Returns
    -------
    vdf_binned : numpy.ndarray
        Binned VDF data with shape: (n_vx, n_vy, n_vz, ntime)
    bin_npts : numpy.ndarray
        Number of points in each bin (n_vx, n_vy, n_vz, ntime)
    """
    # Extracting time array
    _, ntime = vdf_dat.shape

    # Output arrays
    vdf_binned = np.zeros((n_vx, n_vy, n_vz, ntime), dtype=float)
    bin_npts  = np.zeros((n_vx, n_vy, n_vz, ntime), dtype=int)

    # Loop only over time
    for t in range(ntime):

        # Extract bin indices for this time
        ix = vmap[0, :, t]   # v⊥1
        iy = vmap[1, :, t]   # v⊥2
        iz = vmap[2, :, t]   # v∥

        # Mask out invalid bins
        valid = (ix >= 0) & (iy >= 0) & (iz >= 0)

        ixv = ix[valid]
        iyv = iy[valid]
        izv = iz[valid]

        # Values to accumulate
        vals = vdf_dat[valid, t]

        # Scatter-add into 3D grid
        np.add.at(vdf_binned[..., t], (ixv, iyv, izv), vals)
        np.add.at(bin_npts[..., t],  (ixv, iyv, izv), 1)

    return vdf_binned, bin_npts

def bin_vdf_vol(*args, **kwargs):
    """
    This is a wrapper around the `bin_vdf_3d` function to bin vdf data that is 
    volume-weighted.

    Parameters:
    vdf_dat: numpy.ndarray of shape (npts, ntime); VDF data
    vmap: numpy.ndarray of shape (3, npts); Mapping of velocity bins in 
        normalized velocity space
    n_vx, n_vy, n_vz: int; number of bins in each velocity dimension
    interleave_check: bool, Flag to check if the data is interleaved, default=True

    Returns:
    vdf_vol_binned: numpy.ndarray of shape (n_vx, n_vy, n_vz, ntime)
        binned VDF volume weighted data
    bin_npts: numpy.ndarray of shape (n_vx, n_vy, n_vz, ntime)
        number of points in each bin
    """
    vdf_vol_binned, bin_npts = bin_vdf_3d(*args, **kwargs)
    return vdf_vol_binned, bin_npts


def bin_vdf_avg(*args, **kwargs):
    """
    This is a wrapper around the `bin_vdf_3d` function to bin vdf data that is 
    not volume-weighted and hence weight it by averaging.

    Parameters
    ----------
    vdf_dat: numpy.ndarray of shape (npts, ntime); VDF data
    vmap: numpy.ndarray of shape (3, npts); Mapping of velocity bins in 
        normalized velocity space
    n_vx, n_vy, n_vz: int; number of bins in each velocity dimension
    interleave_check: bool, Flag to check if the data is interleaved, default=True

    Returns
    -------
    vdf_avg_binned: numpy.ndarray of shape (n_vx, n_vy, n_vz, ntime)
        binned VDF averaged data
    bin_npts: numpy.ndarray of shape (n_vx, n_vy, n_vz, ntime)
        number of points in each bin
    """
    vdf_raw_binned, bin_npts = bin_vdf_3d(*args, **kwargs)
    with np.errstate(divide='ignore', invalid='ignore'):
        return vdf_raw_binned / bin_npts, bin_npts