import h5py
import numpy as np

def load_box_data(h5file, sim, redshift, field):
    """
    Load a single dataset for a given simulation and redshift.
    """
    zkey = f'z_{int(redshift)}'

    with h5py.File(h5file, 'r') as f:

        if sim not in f:
            raise KeyError(f"Simulation '{sim}' not found in {h5file}")

        sim_grp = f[sim]

        if zkey not in sim_grp:
            raise KeyError(f"Redshift group '{zkey}' not found for sim '{sim}'")

        zgrp = sim_grp[zkey]

        if 'empty' in zgrp.attrs and zgrp.attrs['empty']:
            return None

        if field not in zgrp:
            raise KeyError(f"Field '{field}' not found at z={redshift}")

        return zgrp[field][...]


def load_box_data_dict(h5file, sim, redshift, include_attrs=True):
    """
    Load ALL datasets for a given simulation and redshift.

    Parameters
    ----------
    h5file : str
        Path to HDF5 file
    sim : str
        Simulation name
    redshift : int or float
        Redshift to load
    include_attrs : bool, optional
        If True, include group attributes (e.g. bh_mass, hsml)

    Returns
    -------
    dict or None
        Dictionary of datasets (and attributes if requested),
        or None if snapshot is marked empty.
    """

    zkey = f'z_{int(redshift)}'

    data = {}

    with h5py.File(h5file, 'r') as f:

        if sim not in f:
            raise KeyError(f"Simulation '{sim}' not found in {h5file}")

        sim_grp = f[sim]

        if zkey not in sim_grp:
            raise KeyError(f"Redshift group '{zkey}' not found for sim '{sim}'")

        zgrp = sim_grp[zkey]

        if 'empty' in zgrp.attrs and zgrp.attrs['empty']:
            return None

        # Load datasets
        for key in zgrp.keys():
            data[key] = zgrp[key][...]

        # Optionally load attributes
        if include_attrs:
            data['attrs'] = {k: zgrp.attrs[k] for k in zgrp.attrs.keys()}

    return data
