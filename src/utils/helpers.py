"""
General utility functions for MMS data processing.

Includes time-series manipulation (sampling), metadata generation, 
and file I/O operations for research data.

Author: Regis John
Created: 2026-03-25
"""

import numpy as np
from pyspedas import tres, avg_data


def downsample_cad(tplot_var, ref_tplot_var, trange, newname=None, **kwargs):
    """
    Downsample a tplot variable to match the cadence of a reference tplot variable.

    Parameters
    ----------
    - tplot_var: str
        Name of the tplot variable to downsample.
    - ref_tplot_var: str
        Name of the reference tplot variable to match the cadence of.
    - trange: int/float (optional, but recommended)
        Time range over which to downsample the tplot variable.
    - newname: str
        Name of the downsampled tplot variable. If not set, the default is
        to add a '_dwnsmple' to the input name.
    - **kwargs: dict
        Additional keyword arguments to be passed to the avg_data() function.

    Returns
    -------
    - outname: str
        Name of the downsampled tplot variable.
    """
    # Get target cadence
    des_cadence = tres(ref_tplot_var)

    # Name of downsampled tplot variable
    outname = newname if newname else f'{tplot_var}_dwnsmple'

    # Downsample tplot variable
    return avg_data(tplot_var, trange=trange, res=des_cadence,
                             newname=outname, **kwargs)


def is_2darr_rep(arr, atol=1e-6):
    """
    Checks if a 2D array is a replicated array of the first row.

    Parameters
    ----------
    - arr: numpy.ndarray
        2D array to check for replication.
    - atol: float (optional)
        Absolute tolerance used in the comparison.
        Default is 1e-6.

    Returns
    -------
    - check: bool
        True if the array is a replicated array of the first row, False otherwise.

    Notes
    -----
    - The function works by replicating the first row of the array
      vertically to match the shape of the original array.
    - It then checks if the replicated array is close to the original
      array using numpy's allclose() function.
    """
    arr0 = arr[0, :]  # Takes the first row (time t=0) with shape (32,) across all energy bins.
    test_arr = np.tile(arr0, (arr.shape[0], 1))  # Replicates arr0 vertically to arr.shape[0] creating (N, 32)
    check = np.allclose(arr, test_arr, atol)  # 1 if arrays are close, 0 if not

