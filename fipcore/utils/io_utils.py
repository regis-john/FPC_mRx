"""
================================================================================
Input/Output utilities for MMS data and results.

Functions:
- mms_name_make: Generate a standardized filename based on the given parameters.
- h5sav: Save data to an .h5 file and overwrites existing datasets.

Author: Regis John
Created: 2026-03-26
================================================================================
"""
import h5py

def mms_name_make(prefix, tstart, tend, ext="h5"):
    """
    Generate a standardized filename based on the given parameters.

    Parameters
    ----------
    prefix : str
        Prefix of the filename.
    tstart : str
        Start time of the data segment in the format 'YYYY-MM-DD/hh:mm:ss'.
    tend : str
        End time of the data segment in the format 'YYYY-MM-DD/hh:mm:ss'.
    ext : str, optional
        File extension of the generated filename. Default is 'h5'.

    Returns
    -------
    str
        The standardized filename.

    """
    date = tstart.split('/')[0].replace('-', '')
    s_time = tstart.split('/')[1].replace(':', '')
    e_time = tend.split('/')[1].replace(':', '')
    return f"{prefix}_{date}_{s_time}_{e_time}.{ext}"


def h5sav(hfile, groups_dict):
    """
    Saves data to an .h5 file and overwrites existing datasets.

    Parameters:
    href (str): Path to the HDF5 file to be written.
    groups_dict (dict): Dictionary where the keys are group names and the values 
        are inner dictionary containing the data to be written.

    Returns:
    None
    The function prints out the current file structure after writing is complete.
    """
    with h5py.File(hfile, "a") as f:
        for group_name, data_dict in groups_dict.items():

            grp = f[group_name] if group_name in f else f.create_group(group_name)

            for key, data in data_dict.items():
                if key in grp:
                    del grp[key]      # overwrite cleanly
                grp.create_dataset(key, data=data)

    # --- Verification printout ---
    with h5py.File(hfile, "r") as f: 
        print("\nHDF5 write complete. Current file structure:\n") 
        for group_name in f.keys(): 
            print(f"[{group_name}] -> {list(f[group_name].keys())}")