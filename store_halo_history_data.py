import numpy as np
import h5py
import arepo_package

# All Bondi model, constrained, zoom, and low mass
basePaths = ['/standard/torrey-group/jkho/LtU_Accretion/Constrained/Bondi_constrained_AGN_fewseeds_stellar/output/',
            '/standard/torrey-group/jkho/LtU_Accretion/Zoom/Bondi_zoom_AGN/output/',
            '/standard/torrey-group/jkho/LtU_Accretion/Small_Uniform/Bondi_lowmass_AGN_fewseeds_z127/output/',
            '/standard/torrey-group/jkho/LtU_Accretion/Rainer/FF_Rainer/output/']

desired_redshifts = np.arange(20, 5, -1)
h = 0.6774

Con_Mass = []
Zoom_Mass = []
Low_Mass = []
R_Mass = []

Con_StellarMass = []
Zoom_StellarMass = []
Low_StellarMass = []
R_StellarMass = []

for i in range(len(desired_redshifts)):

    for ii in range(len(basePaths)):

        basePath = basePaths[ii]

        GroupMasses = (
            arepo_package.get_group_property(
                basePath, 'GroupMass', desired_redshifts[i]
            )[0] * 1e10 / h
        )

        GroupMassType = (
            arepo_package.get_group_property(
                basePath, 'GroupMassType', desired_redshifts[i]
            )[0] * 1e10 / h
        )

        StellarMasses = GroupMassType[:, 4]

        if ii == 1 or ii == 3:
            mask = GroupMassType[:, 2] / GroupMassType[:, 1] < 0.01
            GroupMasses = GroupMasses[mask]
            StellarMasses = StellarMasses[mask]

        target = np.argmax(GroupMasses)

        if ii == 0:
            Con_Mass.append(GroupMasses[target])
            Con_StellarMass.append(StellarMasses[target])

        elif ii == 1:
            Zoom_Mass.append(GroupMasses[target])
            Zoom_StellarMass.append(StellarMasses[target])

        elif ii == 2:
            Low_Mass.append(GroupMasses[target])
            Low_StellarMass.append(StellarMasses[target])

        elif ii == 3:
            R_Mass.append(GroupMasses[target])
            R_StellarMass.append(StellarMasses[target])

# Convert to arrays
Con_Mass = np.array(Con_Mass)
Zoom_Mass = np.array(Zoom_Mass)
Low_Mass = np.array(Low_Mass)
R_Mass = np.array(R_Mass)

Con_StellarMass = np.array(Con_StellarMass)
Zoom_StellarMass = np.array(Zoom_StellarMass)
Low_StellarMass = np.array(Low_StellarMass)
R_StellarMass = np.array(R_StellarMass)

# Save to HDF5
outfile = "output/halo_growth_histories.hdf5"

with h5py.File(outfile, "w") as f:
    f.create_dataset("redshift", data=desired_redshifts)

    f.create_dataset("Con_Mass", data=Con_Mass)
    f.create_dataset("Zoom_Mass", data=Zoom_Mass)
    f.create_dataset("Low_Mass", data=Low_Mass)
    f.create_dataset("R_Mass", data=R_Mass)

    f.create_dataset("Con_StellarMass", data=Con_StellarMass)
    f.create_dataset("Zoom_StellarMass", data=Zoom_StellarMass)
    f.create_dataset("Low_StellarMass", data=Low_StellarMass)
    f.create_dataset("R_StellarMass", data=R_StellarMass)

print(f"Saved halo growth histories to {outfile}")