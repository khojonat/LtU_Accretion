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
from LtU_get_property import get_particle_property_LTU
from astropy import constants as c
from astropy import units as u
import matplotlib as mpl
import arepo_package.arepo_package as arepo_package
from astropy import constants as c
from astropy import units as u

Simpath = 'Low_mass_FF_nofeedback'
outputpath = 'FF_seed3.19_nofeedback_r140'
basePath = f'{Simpath}/{outputpath}/'

header = arepo_package.load_snapshot_header(basePath,6)
boxsize = header['BoxSize']
h = 0.6774

GAMMA=5./3
GAMMA_MINUS1=GAMMA-1
redshifts = np.linspace(20,2,10)
a = 1/(1+redshifts)

Masses = []
Mdots = []
SoundSpeed_space = []
Densities = []
all_accretion = []
accretion_masses = []
accretion_pos = []

for i in range(len(redshifts)):

    d = 0.140 # kpc

    BH_Mass = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift = redshifts[i])
    BH_Mdot = get_particle_property_LTU(basePath,'BH_Mdot',p_type=5, desired_redshift = redshifts[i])
    BH_rhos = get_particle_property_LTU(basePath,'BH_Density',p_type=5, desired_redshift = redshifts[i])
    BH_U = get_particle_property_LTU(basePath,'BH_U',p_type=5, desired_redshift = redshifts[i])

    Gas_Mass = get_particle_property_LTU(basePath,'Masses',p_type=0, desired_redshift = redshifts[i])
    Gas_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=0, desired_redshift = redshifts[i])
    BH_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=5, desired_redshift = redshifts[i])
    
    most_massive_ind = np.argmax(BH_Mass[0])

    BH_pos = BH_Pos[0][most_massive_ind] * a[i]/h # kpc
    Gas_pos = Gas_Pos[0] * a[i]/h # kpc

    distances = np.linalg.norm(Gas_pos - BH_pos,axis=1)
    Contributing_gas_mask = distances < d # Only gas cells within a distance d of the BH contribute to accretion
    
    # print(f"Minimum distance: {np.min(distances)} kpc")
    
    if np.sum(Gas_Mass[0][Contributing_gas_mask] * 1e10/h) > 0:
        Total_mgas = Gas_Mass[0][Contributing_gas_mask] * 1e10/h # Msun
    else:
        Total_mgas = Gas_Mass[0][np.argmin(distances)] * 1e10/h
    
    SoundSpeed_space.append(np.sqrt(GAMMA * GAMMA_MINUS1 * BH_U[0][most_massive_ind])) # Units: km/s
    Densities.append(BH_rhos[0][most_massive_ind] * (1e10/h)/(a[i]/h)**3) # Units: solar masses/kpc^3
    Masses.append(BH_Mass[0][most_massive_ind]*1e10/h)
    Mdots.append(BH_Mdot[0][most_massive_ind]*1e10/0.978) # Solar masses/Gyr
    all_accretion.append( (np.sum(BH_QM_Cumgrowth[0])+np.sum(BH_QM_Cumgrowth[0])) * 1e10/h) 
    accretion_masses.append(Total_mgas)
    accretion_pos.append(distances[Contributing_gas_mask])

G = c.G.to(u.km**3/(u.Msun*u.s**2)).value
light_speed = c.c.to(u.km/u.s).value
kpc2km = u.kpc.to(u.km)
sec_per_Gyr = u.Gyr.to(u.s)

specific_accretion = sec_per_Gyr*np.array([100*np.sqrt(2) * G/light_speed * np.sum(accretion_masses[i])/(d*kpc2km)**2 for i in range(len(redshifts))])

plt.scatter(redshifts,np.array(Mdots)/np.array(Masses),label = r'Explicit $\dot{M}/M_{\rm BH}$')
plt.scatter(redshifts,specific_accretion,label = 'Environmentally calculated')
plt.yscale('log')
plt.ylabel('Specific Accretion rate',size=15)
plt.xlabel('Redshift',size=15)
plt.legend()
plt.savefig(f'{Simpath}/Plots/Specific_accretion_comparison_ff_r140.png')
