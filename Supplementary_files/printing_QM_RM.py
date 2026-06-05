import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LogNorm
import numpy as np
import sys 
import os
import h5py
from torreylabtools import helpers
import arepo_package
from astropy import constants as c
from astropy import units as u
from LtU_get_property import get_particle_property_LTU
from astropy.cosmology import Planck18 as cosmo

mpl.use('agg')
mpl.rcParams['text.usetex'] = True

os.makedirs("Plots/Specific", exist_ok=True)

accretion_models = ['Bondi','FF','modFF']

Filepath = '/project/torrey-group/jkho/LtU_accretion'  # '/project/torrey-group/jkho/LtU_accretion/Low_mass_seeds' 
Boxes = ['Low_mass_seeds','Constrained','Zooms']
# Simpath = 'Bondi_lowmass_noAGN_z127'# 'FF_zoom' # 'Bondi_lowmass_noAGN_fewseeds_z127' # 'Bondi_lowmass_noAGN_z127' 
outputpath = 'output'  # 'Bondi_seed3.19_nofeedback_fewseeds'

h = 0.6774
GAMMA=5./3
GAMMA_MINUS1=GAMMA-1
redshifts = np.linspace(15,6,10)
a = 1/(1+redshifts)
bondi_boost = 100

for Box in Boxes:

    path = f'{Filepath}/{Box}/'

    Simpaths = [entry.name for entry in os.scandir(path) if entry.is_dir()]

    for Simpath in Simpaths:

        try:

            matched = [m for m in accretion_models if m in Simpath]

            if 'modFF' in matched:
                accretion_model = 'modFF'
            elif len(matched) == 1:
                accretion_model = matched[0]
            else:
                raise ValueError(f"Could not determine accretion model from {Simpath}")
        
            basePath = f'{Filepath}/{Box}/{Simpath}/{outputpath}/'
    
            Masses = []
            BH_Cum_QM = []
            BH_Cum_RM = []
            
            for i in range(len(redshifts)):
                
                BH_Mass = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift = redshifts[i])
                BH_QM = get_particle_property_LTU(basePath,'BH_CumEgyInjection_QM',p_type=5, desired_redshift = redshifts[i])
                BH_RM = get_particle_property_LTU(basePath,'BH_CumEgyInjection_RM',p_type=5, desired_redshift = redshifts[i])
                
                if len(BH_Mass[0]) == 0:
                    raise ValueError(f"No BH particles in snapshot {redshifts[i]} for {Simpath}")

                most_massive_ind = np.argmax(BH_Mass[0])

                BH_Cum_QM.append(BH_QM[0][most_massive_ind] * (1e10/h) * (a[i]/h)**2 / (0.978/h)**2 ) # Msun*kpc^2/Gyr^2
                BH_Cum_RM.append(BH_RM[0][most_massive_ind] * (1e10/h) * (a[i]/h)**2 / (0.978/h)**2 ) # Msun*kpc^2/Gyr^2

            print(f"{Simpath} QM: {BH_Cum_QM}",flush=True)
            print(f"{Simpath} RM: {BH_Cum_RM}",flush=True)


        except Exception as e:
            print(f"Error in {Box}/{Simpath}: {e}")
    
