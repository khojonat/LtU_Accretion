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
file_path = '.' # '/standard/torrey-group/jkho/LtU_Accretion/mod_Freefall' # '.' # 
Sim_folder = 'FF_noAGN' # 'L12.5n192' # 'modFF_Rainer' # 'FF_noAGN_nostellar' #  '/standard/torrey-group/jkho/LtU_Accretion/Freefall/Low_mass_FF_nofeedback'
output_folder= 'FF_nomod_seed3.19_nofeedback_fewseeds' # 'output' # 'FF_seed3.19_nofeedback_r140' 'FF_seed3.19_noAGN_nostellar_fewseeds' #
basePath = f'{file_path}/{Sim_folder}/{output_folder}/'

header = arepo_package.load_snapshot_header(basePath,6)
boxsize = header['BoxSize']
h = 0.6774

GAMMA=5./3
GAMMA_MINUS1=GAMMA-1
redshifts = np.linspace(12,5,8)
a = 1/(1+redshifts)
bondi_boost = 1
BH_seed_masses = np.array([1e3,1e4,1e5])

AllBHmasses = []
Masses = []
Mdots = []
Progs = []
QM_growth = []
RM_growth = []
Eddington_Rates = []
SoundSpeed_space = []
Densities = []
all_progs = []
all_accretion = []
accretion_masses = []
accretion_pos = []

for i in range(len(redshifts)):

    # d = 0.14 # 0.375 # kpc
    d = 0.0015 # 0.000375 # Mpc
    
    BH_Mass = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift = redshifts[i])
    BH_Mdot = get_particle_property_LTU(basePath,'BH_Mdot',p_type=5, desired_redshift = redshifts[i])
    BH_Progs = get_particle_property_LTU(basePath,'BH_Progs',p_type=5, desired_redshift = redshifts[i])
    BH_QM_Cumgrowth = get_particle_property_LTU(basePath,'BH_CumMassGrowth_QM',p_type=5, desired_redshift = redshifts[i])
    BH_RM_Cumgrowth = get_particle_property_LTU(basePath,'BH_CumMassGrowth_RM',p_type=5, desired_redshift = redshifts[i])
    BH_Eddington = get_particle_property_LTU(basePath,'BH_MdotEddington',p_type=5, desired_redshift = redshifts[i])
    BH_rhos = get_particle_property_LTU(basePath,'BH_Density',p_type=5, desired_redshift = redshifts[i])
    BH_U = get_particle_property_LTU(basePath,'BH_U',p_type=5, desired_redshift = redshifts[i])

    Gas_Mass = get_particle_property_LTU(basePath,'Masses',p_type=0, desired_redshift = redshifts[i])
    Gas_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=0, desired_redshift = redshifts[i])
    BH_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=5, desired_redshift = redshifts[i])
    
    most_massive_ind = np.argmax(BH_Mass[0])

    BH_pos = BH_Pos[0][most_massive_ind] * a[i]/h # kpc, Mpc for Rainer
    Gas_pos = Gas_Pos[0] * a[i]/h # kpc, Mpc for Rainer

    distances = np.linalg.norm(Gas_pos - BH_pos,axis=1)
    Contributing_gas_mask = distances < d # Only gas cells within a distance d of the BH contribute to accretion
    
    # print(f"Minimum distance: {np.min(distances)} kpc")
    
    if np.sum(Gas_Mass[0][Contributing_gas_mask] * 1e10/h) > 0:
        Total_mgas = Gas_Mass[0][Contributing_gas_mask] * 1e10/h # Msun
    else:
        Total_mgas = Gas_Mass[0][np.argmin(distances)] * 1e10/h
    
    AllBHmasses.append(BH_Mass[0]*1e10/h)
    SoundSpeed_space.append(np.sqrt(GAMMA * GAMMA_MINUS1 * BH_U[0][most_massive_ind])) # Units: km/s
    Densities.append(BH_rhos[0][most_massive_ind] * (1e10/h)/(a[i]/h)**3) # Units: solar masses/kpc^3
    Masses.append(BH_Mass[0][most_massive_ind]*1e10/h)
    Mdots.append(BH_Mdot[0][most_massive_ind]*1e10/0.978) # Solar masses/Gyr
    Progs.append(BH_Progs[0][most_massive_ind])
    QM_growth.append(BH_QM_Cumgrowth[0][most_massive_ind]*1e10/h)
    RM_growth.append(BH_RM_Cumgrowth[0][most_massive_ind]*1e10/h)
    Eddington_Rates.append(BH_Eddington[0][most_massive_ind]*1e10/0.978)
    all_progs.append(np.sum(BH_Progs[0]-1))
    all_accretion.append( (np.sum(BH_QM_Cumgrowth[0])+np.sum(BH_QM_Cumgrowth[0])) * 1e10/h) 
    accretion_masses.append(Total_mgas)
    accretion_pos.append(distances[Contributing_gas_mask])

    print(f'Total number of gas cells being accreted: {len(Gas_Mass[0][Contributing_gas_mask])}')
    print(f'Average mass of gas cells being accreted: {np.mean(Gas_Mass[0][Contributing_gas_mask])*1e10/h} Msun')

# Convert from Mpc to kpc if using Rainer parameter file
if 'Rainer' or 'L12.5n192' in Sim_folder:
    unit_time = (u.Mpc/(u.km/u.s)).to(u.Gyr) # Units: Gyr/internal units
    Mdots = np.array(Mdots) * 0.978/unit_time # Units: Msun/Gyr
    Eddington_Rates = np.array(Eddington_Rates) * 0.978/unit_time
    d *= 1e3 # Convert to kpc
    Densities = np.array(Densities) * 1e-9 

# Establish constants and convert to same units: Gyr, Msun, kpc
G = c.G.to(u.kpc**3/(u.Msun*u.Gyr**2)).value 
light_speed = c.c.to(u.kpc/u.Gyr).value 
kpc2km = u.kpc.to(u.km)
sec_per_Gyr = u.Gyr.to(u.s)
Densities = np.array(Densities)
SoundSpeed_space = np.array(SoundSpeed_space)* sec_per_Gyr/kpc2km # new units: kpc/Gyr

# Accretion parameters
A_ff = 1e-3
A_modff = 1e2
alpha_ff = 0
alpha_modff = -0.5

for BH_seed_mass in BH_seed_masses:

    tff = (d**3/(G*BH_seed_mass))**0.5 # Gyr
    R_s = 2*G*BH_seed_mass/light_speed**2
    
    eta_ff = A_ff*(d/R_s)**alpha_ff
    eta_modff = A_modff*(d/R_s)**alpha_modff
    
    specific_accretion_ff = np.array([eta_ff*np.sum(accretion_masses[i]) for i in range(len(redshifts))])/(tff * BH_seed_mass) 
    specific_accretion_ff_mod = np.array([eta_modff*np.sum(accretion_masses[i]) for i in range(len(redshifts))])/(tff * BH_seed_mass)
    
    f_env_Bondi = (SoundSpeed_space)**3/(4*np.pi*G**2*Densities) # Units: Msun*Gyr
    
    specific_Bondi = []
    
    for i in range(len(redshifts)):
    
        rate = BH_seed_mass/f_env_Bondi[i] # Using seed mass at all times
        
        specific_Bondi.append(rate*bondi_boost)
    
    # print(f'Actual rates/100: {np.array(Mdots)/(100*np.array(Masses))}')
    # print(f'Specific accretion for FF: {specific_accretion_ff}\n')
    
    fig,axs = plt.subplots(1,figsize=(6,5))
    
    axs.scatter(redshifts,np.array(Mdots)/np.array(Masses),label = 'Actual rates')
    # axs.scatter(redshifts,np.array(Mdots)/(100*np.array(Masses)),label = 'Actual rates/100')
    axs.scatter(redshifts,specific_accretion_ff,label='ff')
    axs.scatter(redshifts,specific_accretion_ff_mod,label='mod ff')
    axs.scatter(redshifts,specific_Bondi,label='Bondi')
    axs.plot(redshifts,np.array(Mdots)/np.array(Masses))
    # axs.plot(redshifts,np.array(Mdots)/(100*np.array(Masses)),ls = '--', lw = 4)
    axs.plot(redshifts,specific_accretion_ff)
    axs.plot(redshifts,specific_accretion_ff_mod)
    axs.plot(redshifts,specific_Bondi)
    axs.set_yscale('log')
    axs.set_xlabel('Redshift',size=15)
    axs.set_ylabel(r'Specifc accretion rate [$\rm 1/Gyr$]',size=15)
    axs.set_title(rf'$M_\bullet$ = {BH_seed_mass} $M_\odot$',size=17)
    
    fig.suptitle(f'{Sim_folder}',y=0.95)
    fig.legend(fontsize=10)
    
    fig.tight_layout()
    fig.savefig(f'{Sim_folder}/Plots/{Sim_folder}_M{BH_seed_mass}_comparison.png') # {file_path}
