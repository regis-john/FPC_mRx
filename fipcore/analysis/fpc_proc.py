"""
=================================================================================
Functions for processing MMS FPC data.

Functions:
- compute_cprime: Compute the C' correlation term.
- bin_cprime_3d: Bin C' correlation term data into a 3D grid.
- bin_cprime_vol: Bin C' correlation term data into a 3D grid with volume weighting.
- bin_cprime_avg: Bin C' correlation term data into a 3D grid with averaging.
- jEtot_cprime: Compute the total energy density J.E from C' term.
- compute_fpc: Compute the FPC term from the C' correlation term.

Author: Regis John
Created: 2026-03-26
================================================================================
"""
import numpy as np
from scipy.constants import elementary_charge as q_e
from pyspedas import get_data

def compute_cprime(vv, efdat, vdf_dat, species='e'):
    """
    Compute the C' correlation term.

    Parameters:
    vv (numpy.ndarray): Velocity vector data.
    efdat (str): Electric field tplot variable.
    vdf_dat (numpy.ndarray): Velocity distribution function data.
    species (str, optional): Species ('e' or 'ion'), default='e'.

    Returns:
    cprime (numpy.ndarray): C' correlation term.
    """
    CHARGE = {'e': -q_e, 'ion': +q_e}
    q = CHARGE[species]
    
    # Extract electric field data
    times, ef = get_data(efdat)
    # Make E broadcastable: (3, 1, ntime) 
    ef = ef.T[:, None, :]
    ntime = ef.shape[2]

    # Broadcasting vv to all times
    if vv.shape[2] == 2:
        ip = np.arange(ef.shape[2]) % 2
        vv_time = vv[:, :, ip]     # (3, 16384, 133)
    elif vv.shape[2] == ntime:
        vv_time = vv
    elif vv.shape[2] == 1:
        vv_time = np.repeat(vv, ntime, axis=2)

    return q * vv_time * vdf_dat * ef  # (3, 16384, ntime)


def bin_cprime_3d(cprime, vmap, n_vx, n_vy, n_vz):
    """
    Bin C' correlation term data into a 3D grid.

    Parameters:
    cprime (numpy.ndarray): C' correlation term data.
    vmap (numpy.ndarray): Mapping of velocity bins in normalized velocity space.
    n_vx, n_vy, n_vz (int): Number of bins in each velocity dimension.

    Returns:
    Cprime_binned (numpy.ndarray): Binned C' correlation term data.
    bin_npts (numpy.ndarray): Number of points in each bin.
    """
    _, _, ntime = cprime.shape  # (3, 16384, ntime)

    # Output arrays
    Cprime_binned = np.zeros((3, n_vx, n_vy, n_vz, ntime), dtype=float)
    bin_npts      = np.zeros((n_vx, n_vy, n_vz, ntime), dtype=int)

    # Loop over time
    for t in range(ntime):

        # Extract bin indices for this time
        ix = vmap[0, :, t]
        iy = vmap[1, :, t]
        iz = vmap[2, :, t]

        # Mask invalid bins
        valid = (ix >= 0) & (iy >= 0) & (iz >= 0)

        ixv = ix[valid]
        iyv = iy[valid]
        izv = iz[valid]

        # Values to accumulate: shape (3, n_valid)
        vals = cprime[:, valid, t]

        # Scatter-add for each LMN component
        for j in range(3):
            np.add.at(Cprime_binned[j, ..., t], (ixv, iyv, izv), vals[j])
        
        np.add.at(bin_npts[..., t], (ixv, iyv, izv), 1)

    return Cprime_binned, bin_npts


def bin_cprime_vol(*args, **kwargs):
    """
    Wrapper function of `bin_cprime_3d` to bin C' with volume weighting.

    Parameters:
    *args, **kwargs : Passed to `bin_cprime_3d`

    Returns:
    cprime_vol (numpy.ndarray): Binned C' correlation term data with volume weighting.
    bin_npts (numpy.ndarray): Number of points in each bin.
    """
    cprime_vol, bin_npts = bin_cprime_3d(*args, **kwargs)
    return cprime_vol, bin_npts

def bin_cprime_avg(*args, **kwargs):
    """
    Wrapper function of `bin_cprime_3d` to bin C' with averaging.

    Parameters:
    *args, **kwargs : Passed to `bin_cprime_3d`

    Returns:
    cprime_avg (numpy.ndarray): Binned C' correlation term data with averaging.
    bin_npts (numpy.ndarray): Number of points in each bin.
    """
    cprime_raw, bin_npts = bin_cprime_3d(*args, **kwargs)
    with np.errstate(divide='ignore', invalid='ignore'):
        return cprime_raw / bin_npts, bin_npts


def jEtot_cprime(bcpfac_vol):
    """
    Compute the total energy density J.E from C' term.

    Parameters:
    bcprime_vol (numpy.ndarray): Binned C' correlation term data with volume weighting.

    Returns:
    jEtot (float): Total energy density transfer rate in W m^-3.
    """
    # Sum over the 3 velocity-space dimensions
    return np.nansum(bcpfac_vol, axis=(1, 2, 3))


def compute_fpc(cprime, binc_x, binc_y, binc_z):
    """
    Compute the FPC term from the C' correlation term.

    Parameters:
    cprime (numpy.ndarray): C' correlation term.
    binc_x, binc_y, binc_z (numpy.ndarray): Bin centers in x, y, and z.

    Returns:
    fpc (numpy.ndarray): FPC term.
    """
    vx_4d = binc_x[:, None, None, None]
    vy_4d = binc_y[None, :, None, None]
    vz_4d = binc_z[None, None, :, None]

    dvx = np.median(np.diff(binc_x))
    dvy = np.median(np.diff(binc_y))
    dvz = np.median(np.diff(binc_z))
    
    cpx, cpy, cpz = cprime

    dcpx_dvx = np.gradient(cpx, dvx, axis=0)
    dcpy_dvy = np.gradient(cpy, dvy, axis=1)
    dcpz_dvz = np.gradient(cpz, dvz, axis=2)

    cx = -0.5 * vx_4d * dcpx_dvx + 0.5 * cpx
    cy = -0.5 * vy_4d * dcpy_dvy + 0.5 * cpy
    cz = -0.5 * vz_4d * dcpz_dvz + 0.5 * cpz

    return np.stack([cx, cy, cz], axis=0)