import numpy as np
from pyspedas import get_data, store_data, set_coords, set_units

def scpot_corr(nrgy, scpot_dwnsmple, q=-1):
    # Extract arrays from tplot variables
    times, nrgy_arr = get_data(nrgy)
    _, scpot_arr = get_data(scpot_dwnsmple)

    # Number of energy bins (second dimension of nrgy array)
    n_bins = nrgy_arr.shape[1]

    # Expand scpot array from (N,) to (N, 1) to (N, 32)
    scpot_reshaped = np.repeat(scpot_arr[:, np.newaxis], n_bins, axis=1)

    # Apply correction: E_corrected = E_nominal + q * scpot
    nrgy_corr_arr = nrgy_arr + q * scpot_reshaped

    # Name of output tplot variables 
    nrgy_name = f'{nrgy}_corr'
    scpot_name = f'{scpot_dwnsmple}_reshaped'

    # Repacking back to tplot variable
    store_data(nrgy_name, data={'x': times, 'y': nrgy_corr_arr}, attr_dict={'y': {'units': 'eV', 'label': "Energy (eV) (GSE)"}})
    set_coords(nrgy_name, 'GSE')
    set_units(nrgy_name, 'eV')

    store_data(scpot_name, data={'x': times, 'y': scpot_reshaped}, attr_dict={'y': {'units': 'V', 'label': "Reshaped SC potential (V) (GSE)"}})
    set_coords(scpot_name, 'GSE')
    set_units(scpot_name, 'V')
    
    return [nrgy_name, scpot_name]