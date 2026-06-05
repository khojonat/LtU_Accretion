import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LogNorm
import numpy as np
import sys 
import os
import h5py
# from torreylabtools import helpers
import arepo_package_notebook as arepo_package
from astropy import constants as cons
from astropy import units as u
from LtU_get_property_notebook import get_particle_property_LTU
from astropy.cosmology import Planck18 as cosmo
from sklearn.linear_model import LinearRegression

def post_process_accretion(basePath=None, alpha_desired = -1.5, redshifts= np.linspace(15,6,10), constant_rho=None, Mseed=1e5,Rainer=False):
    ''' Post process accretion rates for the 'ff' model for a given sim and choice of alpha '''

    # Accretion radius
    d = 0.140 # kpc

    # Setting unit conversions
    kpc2km = u.kpc.to(u.km)
    Gyr_per_s = u.s.to(u.Gyr)
    G = cons.G.to(u.km**3/(u.Msun*u.s**2)).value
    light_speed = cons.c.to(u.km/u.s).value
    
    # Calculating coefficient based on choice of exponent
    if Rainer:
        # Need to implement Rainer's method here?
        # ( (2 G M_0) / (d c^2) )^alpha = (9.57e-10)^alpha
        M0 = 1e7*u.Msun.cgs
        factor = (2 * cons.G.cgs * M0/(d * u.kpc.cgs * cons.c.cgs**2))**alpha_desired
        A_pred = 0.001 * factor # Scaling based off of ff model for a 10^7 Msun BH
    else:
        m = (2 - (-3)) / (-0.5 - 0)
        b = -3
        A_pred = 10**(m * alpha_desired + b)

    if constant_rho is not None:

        tff = Gyr_per_s * np.sqrt( (d*kpc2km)**3/(G*Mseed) )
        R_s = 2*G*Mseed/light_speed**2 / kpc2km

        Total_masses = constant_rho * (4/3) * np.pi * d**3 # rho in Msun/kpc^3

        specific_accretion = A_pred * (d/R_s)**alpha_desired * Total_masses / tff
        
        return np.nan, specific_accretion # Redshift is N/A in this case
        
    else:
        
        header = arepo_package.load_snapshot_header(basePath,6)
        boxsize = header['BoxSize']
        h = 0.6774
        
        # Initializing values
        GAMMA=5./3
        GAMMA_MINUS1=GAMMA-1
        a = 1/(1+redshifts)
    
        n = len(redshifts)
        Masses = np.empty(n)
        Mdots = np.empty(n)
        SoundSpeed_space = np.empty(n)
        Densities = np.empty(n)
        accretion_masses = []
        accretion_pos = []
        
        for i in range(len(redshifts)):
        
            print(f"Currently loading redshift {redshifts[i]}")
                
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
            
            SoundSpeed_space[i] = (np.sqrt(GAMMA * GAMMA_MINUS1 * BH_U[0][most_massive_ind])) # Units: km/s
            Densities[i] = (BH_rhos[0][most_massive_ind] * (1e10/h)/(a[i]/h)**3) # Units: solar masses/kpc^3
            Masses[i] = (BH_Mass[0][most_massive_ind]*1e10/h)
            Mdots[i] = (BH_Mdot[0][most_massive_ind]*1e10/0.978) # Solar masses/Gyr
            accretion_masses.append(Total_mgas)
            accretion_pos.append(distances[Contributing_gas_mask])
        
        Total_masses = np.array([m.sum() for m in accretion_masses])
        tff = Gyr_per_s * np.sqrt( (d*kpc2km)**3/(G*Masses) )
        R_s = 2*G*Masses/light_speed**2 / kpc2km
        
        # Calculate specific accretion rate
        specific_accretion = A_pred * (d/R_s)**alpha_desired * Total_masses / tff
        
        return redshifts, specific_accretion



