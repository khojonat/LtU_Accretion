import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LogNorm
import numpy as np
import sys 
import os
import h5py
# import helpers
import arepo_package.arepo_package as arepo_package
from astropy import constants as c
from astropy import units as u
from LtU_get_property import get_particle_property_LTU
from astropy.cosmology import Planck18 as cosmo

Sim_folder = 'Low_mass_Bondi_nofeedback'
output_folder = 'Bondi_seed3.19_nofeedback'
basePath = f'{Sim_folder}/{output_folder}/'

header = arepo_package.load_snapshot_header(basePath,6)
boxsize = header['BoxSize']
h = 0.6774

GAMMA=5./3
GAMMA_MINUS1=GAMMA-1
redshifts = np.linspace(15,6,10)
a = 1/(1+redshifts)
bondi_boost = 100

Masses = []
Mdots = []
SoundSpeed_space = []
Densities = []
Eddington_Rates = []

for i in range(len(redshifts)):
    
    BH_Mass = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift = redshifts[i])
    BH_Mdot = get_particle_property_LTU(basePath,'BH_Mdot',p_type=5, desired_redshift = redshifts[i])
    BH_rhos = get_particle_property_LTU(basePath,'BH_Density',p_type=5, desired_redshift = redshifts[i])
    BH_U = get_particle_property_LTU(basePath,'BH_U',p_type=5, desired_redshift = redshifts[i])
    BH_Eddington = get_particle_property_LTU(basePath,'BH_MdotEddington',p_type=5, desired_redshift = redshifts[i])
    
    most_massive_ind = np.argmax(BH_Mass[0])
    
    SoundSpeed_space.append(np.sqrt(GAMMA * GAMMA_MINUS1 * BH_U[0][most_massive_ind])) # Units: km/s
    Densities.append(BH_rhos[0][most_massive_ind] * (1e10/h)/(a[i]/h)**3) # Units: solar masses/kpc^3
    Masses.append(BH_Mass[0][most_massive_ind]*1e10/h)
    Mdots.append(BH_Mdot[0][most_massive_ind]*1e10/0.978)
    Eddington_Rates.append(BH_Eddington[0][most_massive_ind]*1e10/0.978)

G = c.G.to(u.km**3/(u.Msun*u.s**2)).value
f_env = np.array(SoundSpeed_space)**3/(4*np.pi*G**2*np.array(Densities)/c.kpc.to(u.km).value**3) * (1/u.Gyr.to(u.s)) # Units: Msun*Gyr

# Time-redshift sampling
t_vals = cosmo.age(redshifts).to(u.Gyr).value 

specific = []

for i in range(len(redshifts)):

    rate = Masses[i]/f_env[i]
    
    specific.append(rate*bondi_boost)


fig,axs = plt.subplots(1,2,figsize=(10,4))

Accretion_rates = [np.min([np.array(specific[i]),Eddington_Rates[i]/np.array(Masses[i])]) for i in range(len(redshifts))]

axs[0].scatter(redshifts,Accretion_rates,label=r'$M_{\rm seed} = 1e5 M_\odot$',color='green')
axs[0].set_yscale('log')
axs[0].set_xlabel('Redshift')
axs[0].set_ylabel(r'Specifc accretion rate [$\rm 1/Gyr$]',size=15)
axs[0].set_title('Environmental')

axs[1].scatter(redshifts,np.array(Mdots)/np.array(Masses),color='green')
axs[1].set_yscale('log')
axs[1].set_xlabel('Redshift')
axs[1].set_ylabel(r'Specific accretion rate [$\rm 1/Gyr$]',size=15)
axs[1].set_title('Explicit')

# for ax in axs:
#     ax.set_ylim(1e3,5e15)

fig.suptitle(f'{Sim_folder}',y=0.95)
fig.legend(fontsize=10,loc=(0.575,0.75))

fig.tight_layout()
fig.savefig(f'{Sim_folder}/Plots/Bondi_specific_accretion_comparison.png')
