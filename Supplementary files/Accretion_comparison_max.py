import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LogNorm
import numpy as np
import os
import h5py
# import helpers
import arepo_package.arepo_package as arepo_package
from astropy import constants as c
from astropy import units as u
from LtU_get_property import get_particle_property_LTU
from astropy.cosmology import Planck18 as cosmo

# Sim_folder = '/standard/torrey-group/jkho/LtU_Accretion/Bondi'
# output_folder = 'SFMFGM5_FOF10_LW300_seed5.00_Bondi_DFD_3_nofeedback'
Sim_folder = 'Low_mass_FF_nofeedback_fewseeds'
output_folder=  'FF_seed3.19_nofeedback_fewseeds'
basePath = f'{Sim_folder}/{output_folder}/'
Omega_baryon = 0.0486
h = 0.6774
d = 0.140 # kpc

header = arepo_package.load_snapshot_header(basePath,6)
boxsize = header['BoxSize']
DM_mass = header['MassTable'][1]
baryons_mass = 5 * Omega_baryon * DM_mass * 1e10/h # Solar masses, 5 times because that's what it turns out to be...
print('Baryon mass:',baryons_mass)

GAMMA=5./3
GAMMA_MINUS1=GAMMA-1
redshifts = np.linspace(15,6,10)
a = 1/(1+redshifts)
bondi_boost = 1
BH_seed_mass = 1e3

G = c.G.to(u.km**3/(u.Msun*u.s**2)).value
light_speed = c.c.to(u.km/u.s).value
kpc2km = u.kpc.to(u.km)
sec_per_Gyr = u.Gyr.to(u.s)

cells_per_soft = 3 # Shooting for an average of ~30 cells
gas_soft_max = 0.0625 # kpc
max_cells_in_d = np.floor((d / gas_soft_max)**3) * cells_per_soft # Number of gas cells we can fit in a sphere of size ~d^3
accretion_mass = max_cells_in_d*baryons_mass # Total mass within radius d of BH, in solar massses

print(f'Max # of cells: {max_cells_in_d}')

specific_accretion_ff = sec_per_Gyr*0.001 * np.sqrt(G/BH_seed_mass) * accretion_mass/(d*kpc2km)**(3/2)
specific_accretion_ff_mod = sec_per_Gyr*100*np.sqrt(2) * G/light_speed * accretion_mass/(d*kpc2km)**2 

print(f'Freefall specific accretion rate: {specific_accretion_ff} Msun/Gyr')
print(f'Modified freefall specific accretion rate: {specific_accretion_ff_mod} Msun/Gyr')



# f_env_Bondi = np.array(SoundSpeed_space)**3/(4*np.pi*G**2*np.array(Densities)/c.kpc.to(u.km).value**3) * (1/u.Gyr.to(u.s)) # Units: Msun*Gyr

# # Time-redshift sampling
# t_vals = cosmo.age(redshifts).to(u.Gyr).value 

# specific_Bondi = []

# for i in range(len(redshifts)):

#     rate = BH_seed_mass/f_env_Bondi[i] # Using seed mass at all times
    
#     specific_Bondi.append(rate*bondi_boost)


