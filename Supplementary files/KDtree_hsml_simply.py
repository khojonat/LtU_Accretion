#!/usr/bin/env python3
"""
Compute Hsml-like radii for BH particles based on their nearest gas neighbors.

For each BH particle, this script finds the radius that encloses N gas neighbors
(e.g. N=64) using a cKDTree query, for multiple redshifts. Results are stored
in an HDF5 file with one group per redshift.
"""

import sys
import h5py
import numpy as np
from scipy.spatial import cKDTree
from LtU_get_property import get_particle_property_LTU
from arepo_package import *

# ==============================================================
# User parameters
# ==============================================================
Sim = 'Bondi_constrained_AGN_fewseeds_stellar'
basePath = f'/project/torrey-group/jkho/LtU_accretion/Constrained/{Sim}/output/'
N_neighbors = 64

# Array of redshifts to process
desired_redshifts = np.arange(6,15,1)

# Output file
outfile = f"output/hsmls/BH_gas_hsml_N{N_neighbors}_{Sim}.hdf5"

# ==============================================================
# Open HDF5 file (append mode)
# ==============================================================
with h5py.File(outfile, "a") as f:

    # Store some global metadata once
    f.attrs.setdefault("Simulation", Sim)
    f.attrs.setdefault("N_neighbors", N_neighbors)

    # ==========================================================
    # Loop over redshifts
    # ==========================================================
    for z in desired_redshifts:

        z_label = f"z_{z:.2f}"

        if z_label in f:
            print(f"Group {z_label} already exists — skipping.")
            continue

        print(f"\nLoading coordinates from {Sim} at redshift {z}...")

        gas_coords, zout = get_particle_property_LTU(
            basePath,
            desired_property='Coordinates',
            p_type=0,
            desired_redshift=z
        )

        BH_coords, zout = get_particle_property_LTU(
            basePath,
            desired_property='Coordinates',
            p_type=5,
            desired_redshift=z
        )

        boxsize = get_box_size(basePath)

        # ======================================================
        # KDTree query
        # ======================================================
        print("Building KDTree for gas particles ...")
        gas_tree = cKDTree(gas_coords, boxsize=boxsize + 1)

        print(f"Finding distances to {N_neighbors} nearest gas neighbors ...")
        distances, _ = gas_tree.query(BH_coords, k=N_neighbors)

        if N_neighbors > 1:
            hsml = distances[:, -1]
        else:
            hsml = distances

        # ======================================================
        # Write to HDF5
        # ======================================================
        grp = f.create_group(z_label)

        grp.create_dataset(
            "hsml",
            data=hsml,
            compression="gzip",
            compression_opts=4
        )

        grp.attrs["redshift"] = float(zout)
        grp.attrs["N_BH"] = hsml.size
        grp.attrs["boxsize"] = boxsize

        print(
            f"Saved {hsml.size} hsml values for z={z:.2f} "
            f"(mean={np.mean(hsml):.3f})"
        )

print(f"\nAll results written to {outfile}")
print("Done.")
