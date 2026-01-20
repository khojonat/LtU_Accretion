import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LogNorm
import numpy as np
import os
import h5py
from astropy import constants as c
from astropy import units as u
import matplotlib as mpl
from torreylabtools import helpers
from LtU_get_property import get_particle_property_LTU
import arepo_package as arepo_package

# -------------------------------------------------------
# Paths and box definitions
# -------------------------------------------------------

Filepath = '/project/torrey-group/jkho/LtU_accretion'
Boxes = ['Zooms', 'Constrained', 'Low_mass_seeds']
outputpath = 'output'

# Redshifts
redshifts = np.linspace(25, 6, 20) # redshifts = np.linspace(12, 6, 7)
a = 1.0 / (1.0 + redshifts)
h = 0.6774

# Mapping Simpaths index → HDF5 group
group_labels = ["AGN_S", "AGN_NS", "NAGN_S", "NAGN_NS"]

# -------------------------------------------------------
# Loop through boxes
# -------------------------------------------------------

for Box in Boxes:

    # ---------------------------------------------------
    # Define Simpaths for this Box
    # ---------------------------------------------------
    if Box == 'Zooms':
        Simpaths = [
            ['Bondi_zoom_AGN',                 'FF_zoom_AGN',                 'modFF_zoom_AGN'],                 # AGN_S
            ['Bondi_zoom_AGN_0.1stellar',      'FF_zoom_AGN_0.1stellar',      'modFF_zoom_AGN_0.1stellar'],      # AGN_NS
            ['Bondi_zoom_boost_noAGN_stellar', 'FF_zoom_noAGN_stellar',        'modFF_zoom_noAGN_stellar'],       # NAGN_S
            ['Bondi_zoom_noAGN_0.1stellar',    'FF_zoom_noAGN_0.1stellar',     'modFF_zoom_noAGN_0.1stellar']     # NAGN_NS
        ]

    elif Box == 'Low_mass_seeds':
        Simpaths = [
            ['Bondi_lowmass_AGN_fewseeds_z127',            'FF_lowmass_AGN_fewseeds_z127',            'modFF_lowmass_AGN_fewseeds_z127'],
            ['Bondi_lowmass_AGN_fewseeds_lowstellar_z127', 'FF_lowmass_AGN_fewseeds_lowstellar_z127', 'modFF_lowmass_AGN_fewseeds_lowstellar_z127'],
            ['Bondi_lowmass_noAGN_fewseeds_z127',          'FF_lowmass_noAGN_fewseeds_z127',          'modFF_lowmass_noAGN_fewseeds_z127'],
            ['Bondi_lowmass_noAGN_fewseeds_nostellar_z127','FF_lowmass_noAGN_fewseeds_nostellar_z127','modFF_lowmass_noAGN_fewseeds_nostellar_z127']
        ]

    elif Box == 'Constrained':
        Simpaths = [
            ['Bondi_constrained_AGN_fewseeds_stellar',   'FF_constrained_AGN_fewseeds_stellar',   'modFF_constrained_AGN_fewseeds_stellar'],
            ['Bondi_constrained_AGN_fewseeds_0.1stellar','FF_constrained_AGN_fewseeds_0.1stellar','modFF_constrained_AGN_fewseeds_0.1stellar'],
            ['Bondi_constrained_noAGN_fewseeds',         'FF_constrained_noAGN_fewseeds',         'modFF_constrained_noAGN_fewseeds'],
            ['Bondi_constrained_noAGN_0.1stellar',       'FF_constrained_noAGN_0.1stellar',       'modFF_constrained_noAGN_0.1stellar']
        ]

    # -------------------------------------------------------
    # Prepare output HDF5 file
    # -------------------------------------------------------

    outname = f"output/{Box}_25.hdf5"
    print(f"\n### Writing output for {Box} → {outname}",flush = True)

    hf = h5py.File(outname, "w")

    # -------------------------------------------------------
    # Loop through the 4 categories (AGN_S, AGN_NS, ...)
    # -------------------------------------------------------
    for gi, label in enumerate(group_labels):

        print(f"\n### Processing group: {label}",flush = True)
        label_group = hf.create_group(label)

        # Each element: [Bondi_path, FF_path, modFF_path]
        bondi_path, ff_path, modff_path = Simpaths[gi]

        subnames = ["Bondi", "FF", "modFF"]
        simset = [bondi_path, ff_path, modff_path]

        # ---------------------------------------------------
        # Loop over Bondi, FF, modFF
        # ---------------------------------------------------
        for subname, sim in zip(subnames, simset):

            print(f"  Loading: {sim}", flush = True)
            subgroup = label_group.create_group(subname)

            basePath = f"{Filepath}/{Box}/{sim}/{outputpath}/"

            # Lists per redshift
            Masses = []
            Mdots = []
            Accretion_growth = []
            Edd = []
            Progs = []
            # rho = []

            # ---------------------------------------------------
            # Load data across redshifts
            # ---------------------------------------------------
            for i, z in enumerate(redshifts):

                BH_Mass  = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift=z)
                BH_Mdot  = get_particle_property_LTU(basePath,'BH_Mdot',p_type=5, desired_redshift=z)
                BH_QM    = get_particle_property_LTU(basePath,'BH_CumMassGrowth_QM',p_type=5, desired_redshift=z)
                BH_RM    = get_particle_property_LTU(basePath,'BH_CumMassGrowth_RM',p_type=5, desired_redshift=z)
                BH_Edd   = get_particle_property_LTU(basePath,'BH_MdotEddington',p_type=5, desired_redshift=z)
                BH_Progs = get_particle_property_LTU(basePath,'BH_Progs',p_type=5, desired_redshift=z)
                # BH_rho   = get_particle_property_LTU(basePath,'SubfindDensity',p_type=5, desired_redshift=z)
                    
                if BH_Mass is None:
                                        
                    Masses.append(0)
                    Mdots.append(0)
                    Accretion_growth.append(0)
                    Edd.append(0)
                    Progs.append(0)
                    
                    continue

                most_massive_ind = np.argmax(BH_Mass[0])

                Masses.append(BH_Mass[0][most_massive_ind] * 1e10 / h)
                Mdots.append(BH_Mdot[0][most_massive_ind] * 1e10 / 0.978)
                Accretion_growth.append((BH_QM[0][most_massive_ind] + BH_RM[0][most_massive_ind]) * 1e10 / h)
                Edd.append(BH_Edd[0][most_massive_ind] * 1e10 / 0.978)
                Progs.append(BH_Progs[0][most_massive_ind])
                # rho.append(BH_rho[0][most_massive_ind] * (1e10/h) / (a[i]/h)**3)

            # ---------------------------------------------------
            # Store datasets
            # ---------------------------------------------------
            subgroup.create_dataset("Masses", data=np.array(Masses))
            subgroup.create_dataset("Mdots",  data=np.array(Mdots))
            subgroup.create_dataset("Accretion_growth", data=np.array(Accretion_growth))
            subgroup.create_dataset("Edd", data=np.array(Edd))
            subgroup.create_dataset("Progs", data=np.array(Progs))
            # subgroup.create_dataset("rho", data=np.array(rho))

    hf.close()
    print(f"### Finished writing {outname}\n", flush = True)

print("\nAll boxes completed successfully.", flush = True)
