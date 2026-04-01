"""
================================================================================
Functions for Lorentz Transforming the MMS Electric Field to the Reconnection 
Frame in both FAC and LMN frames.

Functions:
- evec_to_recx: Lorentz Transform the Electric Field to the Reconnection Frame
- evec_to_fac: Rotate the Electric Field to the FAC Frame
- evec_to_lmn: Rotate the Electric Field to the LMN Frame

Author: Regis John
Created: 2026-03-30
================================================================================
"""
import numpy as np
from pyspedas import get_data, store_data, set_units, set_coords
from src.utils.coord_utils import rotate_to_fac, rotate_to_lmn

def evec_to_recx(e_field_var, b_field_var, v_bulk_mean, newname=None):
    """
    Lorentz Transfor the electric field to the reconnection frame.
    
    Parameters:
    e_field_var (str): Name of the tplot variable containing the electric field data.
    b_field_var (str): Name of the tplot variable containing the magnetic field data.
    v_bulk_mean (numpy.ndarray): Mean bulk velocity in km/s.
    newname (str, optional): Name of the output tplot variable, defaults to "evec_lor".
    
    Returns:
    out_name (str): Name of the output tplot variable.
    """
    # Get data from tplot variables
    _, e_data = get_data(e_field_var)
    times, b_data = get_data(b_field_var)
    
    # Cross product <v_i> × B(t) for all t
    v_cross_b = np.cross(v_bulk_mean, b_data)
    
    # Convert km/s * nT of vxB to [mV/m] using 1e-3 scaling factor
    ef_lor = e_data + (v_cross_b * 1e-3)
    
    # Define output name
    out_name = newname if newname else "evec_lor"
    
    # Store the result back as a tplot variable
    store_data(out_name, data={'x': times, 'y': ef_lor}, 
               attr_dict={'y': {'units': 'mV/m', 'label': "Electric Field (mV/m) in Reconnection Frame."}})
    set_coords(out_name, 'Recx Frame')
    set_units(out_name, 'mV/m')
    
    return out_name


def evec_to_fac(e_lor_var, fac_mat_var, newname=None):
    """
    Rotate the lorentz transformed electric field to the FAC frame.

    Parameters:
    e_lor_var (str): Name of the tplot variable containing the  electric field data.
    fac_mat_var (str): Name of the tplot variable containing the FAC rotation matrix.
    newname (str, optional): Name of the output tplot variable, defaults to "evec_lor_fac".

    Returns:
    out_name (str): Name of the output tplot variable.
    """
    times, e_fac = rotate_to_fac(e_lor_var, fac_mat_var)
    
    # Store the result back as a tplot variable
    out_name = newname if newname else "evec_lor_fac"
    store_data(out_name, data={'x': times, 'y': e_fac}, 
               attr_dict={'y': {'units': 'mV/m', 'label': "Electric Field (mV/m) FAC."}})
    set_coords(out_name, 'FAC')
    set_units(out_name, 'mV/m')
    
    return out_name


def evec_to_lmn(e_lor_var, lmn_mat_var, newname=None):
    """
    Rotate the lorentz transformed electric field to the LMN frame.

    Parameters:
    e_lor_var (str): Name of the tplot variable containing the electric field data.
    lmn_mat_var (str): Name of the tplot variable containing the LMN rotation matrix.
    newname (str, optional): Name of the output tplot variable, defaults to "evec_lor_lmn".

    Returns:
    out_name (str): Name of the output tplot variable.
    """
    times, e_lmn = rotate_to_lmn(e_lor_var, lmn_mat_var)

    # Store result
    out_name = newname if newname else "evec_lor_lmn"
    store_data(out_name, data={'x': times, 'y': e_lmn}, 
               attr_dict={'y': {'units': 'mV/m', 'label': "Electric Field (mV/m) LMN."}})
    set_coords(out_name, 'LMN')
    set_units(out_name, 'mV/m')

    return out_name