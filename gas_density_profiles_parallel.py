import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as cons
import h5py
import os
from matplotlib.colors import LightSource, LogNorm
from scipy.spatial import KDTree
from multiprocessing import Pool, cpu_count

from torreylabtools.visualization import contour_makepic as makepic
from torreylabtools.util import calc_hsml as calc_hsml
from torreylabtools.util import naming as namingx
from torreylabtools import kdtree_smoothing as kda
from LtU_get_property import get_particle_property_LTU
import arepo_package as arepo_package
from density_profile import fibonacci_sphere, calc_density


# ==============================================================
# Configuration
# ==============================================================
Filepath = '/project/torrey-group/jkho/LtU_accretion'
Box = 'Zooms'  # 'Constrained' or 'Low_mass_seeds'
SimPaths = [
    ['Bondi_zoom_boost', 'FF_zoom', 'modFF_zoom'],
    ['Bondi_zoom_AGN', 'FF_zoom_AGN', 'modFF_zoom_AGN']
]
outputpath = 'output'

desired_redshifts = np.arange(12, 5, -1)
header = arepo_package.load_snapshot_header(
    f'{Filepath}/{Box}/{SimPaths[0][0]}/{outputpath}', desired_redshifts[0]
)
DMMass = header['MassTable'][1]
a = 1 / (1 + desired_redshifts)
h = 0.6774

# ==============================================================
# Load particle data (sequential; mostly I/O)
# ==============================================================
def load_gas_data():
    Bondi_mass_noAGN, Bondi_pos_noAGN = [], []
    FF_mass_noAGN, FF_pos_noAGN = [], []
    modFF_mass_noAGN, modFF_pos_noAGN = [], []

    Bondi_mass_AGN, Bondi_pos_AGN = [], []
    FF_mass_AGN, FF_pos_AGN = [], []
    modFF_mass_AGN, modFF_pos_AGN = [], []

    for ii in range(len(SimPaths)):
        for sim in SimPaths[ii]:
            basePath = f'{Filepath}/{Box}/{sim}/{outputpath}/'

            for iii in range(len(desired_redshifts)):
                GasMass = (
                    get_particle_property_LTU(basePath, 'Masses', p_type=0, desired_redshift=desired_redshifts[iii])[0]
                    * 1e10 / h
                )  # Msun
                GasPos = (
                    get_particle_property_LTU(basePath, 'Coordinates', p_type=0, desired_redshift=desired_redshifts[iii])[0]
                    * a[iii] / h
                )  # kpc

                BHMasses, _ = arepo_package.get_subhalo_property(basePath, 'SubhaloBHMass', desired_redshifts[iii])
                target = np.argmax(BHMasses)
                offset = arepo_package.get_subhalo_property(basePath, 'SubhaloOffsetType', desired_redshifts[iii])[0][target]
                length = arepo_package.get_subhalo_property(basePath, 'SubhaloLenType', desired_redshifts[iii])[0][target]

                GalaxyGasPos = GasPos[offset[0]:offset[0] + length[0]]
                GalaxyGasMass = GasMass[offset[0]:offset[0] + length[0]]

                # Assign based on sim type
                if ii == 0:
                    if sim == SimPaths[ii][0]:
                        Bondi_mass_noAGN.append(GalaxyGasMass)
                        Bondi_pos_noAGN.append(GalaxyGasPos)
                    elif sim == SimPaths[ii][1]:
                        FF_mass_noAGN.append(GalaxyGasMass)
                        FF_pos_noAGN.append(GalaxyGasPos)
                    elif sim == SimPaths[ii][2]:
                        modFF_mass_noAGN.append(GalaxyGasMass)
                        modFF_pos_noAGN.append(GalaxyGasPos)
                elif ii == 1:
                    if sim == SimPaths[ii][0]:
                        Bondi_mass_AGN.append(GalaxyGasMass)
                        Bondi_pos_AGN.append(GalaxyGasPos)
                    elif sim == SimPaths[ii][1]:
                        FF_mass_AGN.append(GalaxyGasMass)
                        FF_pos_AGN.append(GalaxyGasPos)
                    elif sim == SimPaths[ii][2]:
                        modFF_mass_AGN.append(GalaxyGasMass)
                        modFF_pos_AGN.append(GalaxyGasPos)

    return (Bondi_mass_noAGN, Bondi_pos_noAGN, FF_mass_noAGN, FF_pos_noAGN,
            modFF_mass_noAGN, modFF_pos_noAGN,
            Bondi_mass_AGN, Bondi_pos_AGN, FF_mass_AGN, FF_pos_AGN,
            modFF_mass_AGN, modFF_pos_AGN)


# ==============================================================
# Density calculation (parallel)
# ==============================================================
rmin = 0.1     # kpc
rmax = 50      # kpc
sphere_samples = 200
radial_samples = 50
DesNgb = 32

all_r = np.logspace(np.log10(rmin), np.log10(rmax), radial_samples)


def compute_densities_for_redshift(i):
    """Compute all density profiles for a single redshift index."""
    Bondi_dens_noAGN, FF_dens_noAGN, modFF_dens_noAGN = [], [], []
    Bondi_dens_AGN, FF_dens_AGN, modFF_dens_AGN = [], [], []

    Bondi_center_noAGN = np.median(Bondi_pos_noAGN[i], axis=0)
    FF_center_noAGN = np.median(FF_pos_noAGN[i], axis=0)
    modFF_center_noAGN = np.median(modFF_pos_noAGN[i], axis=0)

    Bondi_center_AGN = np.median(Bondi_pos_AGN[i], axis=0)
    FF_center_AGN = np.median(FF_pos_AGN[i], axis=0)
    modFF_center_AGN = np.median(modFF_pos_AGN[i], axis=0)

    for r in all_r:
        points = fibonacci_sphere(sphere_samples, r)

        # --- No AGN cases ---
        Bondi_dens_noAGN.append(calc_density(Bondi_pos_noAGN[i] - Bondi_center_noAGN, Bondi_mass_noAGN[i], points, DesNgb))
        FF_dens_noAGN.append(calc_density(FF_pos_noAGN[i] - FF_center_noAGN, FF_mass_noAGN[i], points, DesNgb))
        modFF_dens_noAGN.append(calc_density(modFF_pos_noAGN[i] - modFF_center_noAGN, modFF_mass_noAGN[i], points, DesNgb))

        # --- AGN cases ---
        Bondi_dens_AGN.append(calc_density(Bondi_pos_AGN[i] - Bondi_center_AGN, Bondi_mass_AGN[i], points, DesNgb))
        FF_dens_AGN.append(calc_density(FF_pos_AGN[i] - FF_center_AGN, FF_mass_AGN[i], points, DesNgb))
        modFF_dens_AGN.append(calc_density(modFF_pos_AGN[i] - modFF_center_AGN, modFF_mass_AGN[i], points, DesNgb))

    return (
        i,
        np.array(Bondi_dens_noAGN),
        np.array(Bondi_dens_AGN),
        np.array(FF_dens_noAGN),
        np.array(FF_dens_AGN),
        np.array(modFF_dens_noAGN),
        np.array(modFF_dens_AGN),
    )


# ==============================================================
# Main execution
# ==============================================================
if __name__ == "__main__":
    (Bondi_mass_noAGN, Bondi_pos_noAGN, FF_mass_noAGN, FF_pos_noAGN,
     modFF_mass_noAGN, modFF_pos_noAGN,
     Bondi_mass_AGN, Bondi_pos_AGN, FF_mass_AGN, FF_pos_AGN,
     modFF_mass_AGN, modFF_pos_AGN) = load_gas_data()

    nproc = min(len(desired_redshifts), cpu_count())
    print(f"Using {nproc} processes to compute densities across redshifts...")

    with Pool(processes=nproc) as pool:
        results = pool.map(compute_densities_for_redshift, range(len(desired_redshifts)))

    # Sort results by redshift index
    results.sort(key=lambda x: x[0])

    Bondi_density_noAGN = [r[1] for r in results]
    Bondi_density_AGN = [r[2] for r in results]
    FF_density_noAGN = [r[3] for r in results]
    FF_density_AGN = [r[4] for r in results]
    modFF_density_noAGN = [r[5] for r in results]
    modFF_density_AGN = [r[6] for r in results]

    # ==============================================================
    # Save results
    # ==============================================================
    with h5py.File('density_profiles.hdf5', 'w') as f:
        f.create_dataset('Bondi_density_noAGN', data=np.array(Bondi_density_noAGN))
        f.create_dataset('Bondi_density_AGN', data=np.array(Bondi_density_AGN))
        f.create_dataset('FF_density_noAGN', data=np.array(FF_density_noAGN))
        f.create_dataset('FF_density_AGN', data=np.array(FF_density_AGN))
        f.create_dataset('modFF_density_noAGN', data=np.array(modFF_density_noAGN))
        f.create_dataset('modFF_density_AGN', data=np.array(modFF_density_AGN))

    print(" Density profiles saved to density_profiles.hdf5")



