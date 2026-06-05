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

    if Box == 'Low_mass_seeds':
        d = 0.140 # kpc
    elif Box == 'Constrained':
        d = 0.375
    elif Box == 'Zooms':
        d = 0.140

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
            Mdots = []
            SoundSpeed_space = []
            Densities = []
            Eddington_Rates = []
            accretion_masses = []
            accretion_pos = []
            BH_Cum_QM = []
            BH_Cum_RM = []
            
            for i in range(len(redshifts)):
                
                BH_Mass = get_particle_property_LTU(basePath,'BH_Mass',p_type=5, desired_redshift = redshifts[i])
                BH_Mdot = get_particle_property_LTU(basePath,'BH_Mdot',p_type=5, desired_redshift = redshifts[i])
                BH_rhos = get_particle_property_LTU(basePath,'BH_Density',p_type=5, desired_redshift = redshifts[i])
                BH_U = get_particle_property_LTU(basePath,'BH_U',p_type=5, desired_redshift = redshifts[i])
                BH_QM = get_particle_property_LTU(basePath,'BH_CumEgyInjection_QM',p_type=5, desired_redshift = redshifts[i])
                BH_RM = get_particle_property_LTU(basePath,'BH_CumEgyInjection_RM',p_type=5, desired_redshift = redshifts[i])
                Gas_Mass = get_particle_property_LTU(basePath,'Masses',p_type=0, desired_redshift = redshifts[i])
                Gas_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=0, desired_redshift = redshifts[i])
                BH_Pos = get_particle_property_LTU(basePath,'Coordinates',p_type=5, desired_redshift = redshifts[i])
                BH_Eddington = get_particle_property_LTU(basePath,'BH_MdotEddington',p_type=5, desired_redshift = redshifts[i])

                if len(BH_Mass[0]) == 0:
                    raise ValueError(f"No BH particles in snapshot {redshifts[i]} for {Simpath}")

                most_massive_ind = np.argmax(BH_Mass[0])

                BH_pos = BH_Pos[0][most_massive_ind] * a[i]/h # kpc
                Gas_pos = Gas_Pos[0] * a[i]/h # kpc
            
                distances = np.linalg.norm(Gas_pos - BH_pos,axis=1)
                Contributing_gas_mask = distances < d # Only gas cells within a distance d of the BH contribute to accretion
                
                if np.sum(Gas_Mass[0][Contributing_gas_mask] * 1e10/h) > 0:
                    Total_mgas = Gas_Mass[0][Contributing_gas_mask] * 1e10/h # Msun
                else:
                    Total_mgas = Gas_Mass[0][np.argmin(distances)] * 1e10/h
                
                SoundSpeed_space.append(np.sqrt(GAMMA * GAMMA_MINUS1 * BH_U[0][most_massive_ind])) # Units: km/s
                Densities.append(BH_rhos[0][most_massive_ind] * (1e10/h)/(a[i]/h)**3) # Units: solar masses/kpc^3
                Masses.append(BH_Mass[0][most_massive_ind]*1e10/h)
                Mdots.append(BH_Mdot[0][most_massive_ind]*1e10/0.978) # Solar masses/Gyr
                # all_accretion.append( (np.sum(BH_QM_Cumgrowth[0])+np.sum(BH_QM_Cumgrowth[0])) * 1e10/h) 
                accretion_masses.append(Total_mgas)
                accretion_pos.append(distances[Contributing_gas_mask])
                Eddington_Rates.append(BH_Eddington[0][most_massive_ind]*1e10/0.978)
                BH_Cum_QM.append(BH_QM[0][most_massive_ind] * (1e10/h) * (a[i]/h)**2 / (0.978/h)**2 ) # Msun*kpc^2/Gyr^2
                BH_Cum_RM.append(BH_RM[0][most_massive_ind] * (1e10/h) * (a[i]/h)**2 / (0.978/h)**2 ) # Msun*kpc^2/Gyr^2

            print(f"{Simpath} QM: {BH_Cum_QM}",flush=True)
            print(f"{Simpath} RM: {BH_Cum_RM}",flush=True)
            
            if accretion_model == 'Bondi':
                
                G = c.G.to(u.km**3/(u.Msun*u.s**2)).value 

                specific_accretion = []

                f_env = np.array(SoundSpeed_space)**3/(4*np.pi*G**2*np.array(Densities)/c.kpc.to(u.km).value**3) * (1/u.Gyr.to(u.s))
            
                for i in range(len(redshifts)):
                
                    rate = Masses[i]/f_env[i]
                    
                    specific_accretion.append(rate*bondi_boost)

                Accretion_rates = [np.min([np.array(specific_accretion[i]),Eddington_Rates[i]/np.array(Masses[i])]) for i in range(len(redshifts))]
    
                fig,axs = plt.subplots(4,1,figsize = (8,10),sharex = True)
    
                axs[0].plot(redshifts,np.array(Mdots)/np.array(Masses),label = r'Explicit $\dot{M}/M_{\rm BH}$')
                axs[0].plot(redshifts,Accretion_rates,label = 'Environmentally calculated',alpha=0.5)
                axs[0].set_ylabel('Specific Accretion rate [1/Gyr]',size=10)
                axs[0].legend(fontsize=10)
                
                axs[1].plot(redshifts,Densities)
                axs[1].set_ylabel(r'Gas density $\rm [M_\odot/kpc^3]$',size=10)
                axs[1].set_xlabel('Redshift',size=10)
    
                axs[2].plot(redshifts,SoundSpeed_space,color='green')
                axs[2].set_ylabel(r'Sound speed [$\rm km/s$]',size=10)
                
                axs[3].plot(redshifts,BH_Cum_QM,label = 'QM')
                axs[3].plot(redshifts,BH_Cum_RM,label = 'RM')
                axs[3].set_ylabel(r'Cumulative energy injection [$M_\odot * kpc^2/Gyr^2$]',size=10)
                axs[3].legend(fontsize=10)

                for ax in axs:
                    ax.set_yscale('log')
    
                fig.tight_layout()
                plt.subplots_adjust(hspace=0)
    
                fig.savefig(f'Plots/Specific/{Simpath}_specific.png')

                plt.close(fig)
    
                continue
                    
            elif accretion_model == 'FF':
    
                G = c.G.to(u.km**3/(u.Msun*u.Gyr**2)).value
                light_speed = c.c.to(u.km/u.Gyr).value
                kpc2km = u.kpc.to(u.km)

                specific_accretion = 0.001 * np.sqrt(G/(np.array(Masses)))/(d*kpc2km)**(3/2) * np.array([np.sum(accretion_masses[i]) for i in range(len(redshifts))])
    
            elif accretion_model == 'modFF':
    
                G = c.G.to(u.km**3/(u.Msun*u.s**2)).value
                light_speed = c.c.to(u.km/u.s).value
                kpc2km = u.kpc.to(u.km)
                sec_per_Gyr = u.Gyr.to(u.s)
    
                specific_accretion = sec_per_Gyr * np.array([100*np.sqrt(2) * G/light_speed * np.sum(accretion_masses[i]) / (d*kpc2km)**2 for i in range(len(redshifts))])
    
            fig,axs = plt.subplots(3,1,figsize = (8,8),sharex = True)

            Accretion_rates = [np.min([np.array(specific_accretion[i]),Eddington_Rates[i]/np.array(Masses[i])]) for i in range(len(redshifts))]
    
            axs[0].plot(redshifts,np.array(Mdots)/np.array(Masses),label = r'Explicit $\dot{M}/M_{\rm BH}$')
            axs[0].plot(redshifts,Accretion_rates,label = 'Environmentally calculated')
            axs[0].set_ylabel('Specific Accretion rate [1/Gyr]',size=10)
            axs[0].legend(fontsize=10)
            
            densities = np.array([np.sum(accretion_masses[i]) for i in range(len(redshifts))])/((4/3)*np.pi*(d)**3)
            axs[1].plot(redshifts,densities)
            axs[1].set_ylabel(r'Gas density within d $\rm [M_\odot/kpc^3]$',size=10)
            axs[1].set_xlabel('Redshift',size=10)

            axs[2].plot(redshifts,BH_Cum_QM,label = 'QM')
            axs[2].plot(redshifts,BH_Cum_RM,label = 'RM')
            axs[2].set_ylabel(r'Cumulative energy injection [$M_\odot * kpc^2/Gyr^2$]',size=10)
            axs[2].legend(fontsize=10)

            for ax in axs:
                ax.set_yscale('log')
    
            fig.tight_layout()
            fig.savefig(f'Plots/Specific/{Simpath}_specific.png')

            plt.close(fig)

        except Exception as e:
            print(f"Error in {Box}/{Simpath}: {e}")
    
