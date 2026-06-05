import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as cons
import h5py
import os
from matplotlib.colors import LightSource, LogNorm
from scipy.spatial import KDTree
 
from torreylabtools.visualization import contour_makepic as makepic
from torreylabtools.util import calc_hsml as calc_hsml
from torreylabtools.util import naming as namingx
from torreylabtools import kdtree_smoothing as kda
from LtU_get_property import get_particle_property_LTU
import arepo_package as arepo_package
from density_profile import fibonacci_sphere, calc_density

Filepath = '/project/torrey-group/jkho/LtU_accretion' 
Box = 'Constrained' # 'Low_mass_seeds' # 'Zooms' #
SimPaths = [['Bondi_constrained_AGN_fewseeds_stellar','Bondi_constrained_noAGN_fewseeds_boost','Bondi_constrained_AGN_fewseeds_0.1stellar', 'Bondi_constrained_noAGN_0.1stellar',],
            ['FF_constrained_AGN_fewseeds_stellar','FF_constrained_noAGN_fewseeds','FF_constrained_AGN_fewseeds_0.1stellar','FF_constrained_noAGN_0.1stellar'],
            ['modFF_constrained_AGN_fewseeds_stellar','modFF_constrained_noAGN_fewseeds','modFF_constrained_AGN_fewseeds_0.1stellar','modFF_constrained_noAGN_0.1stellar']] 

#Zooms: [['Bondi_zoom_boost','FF_zoom','modFF_zoom'],['Bondi_zoom_AGN','FF_zoom_AGN','modFF_zoom_AGN']] 
outputpath = 'output'
Ptype = 1 # 0

print('Loading DM Mass ...',flush = True)

desired_redshifts=np.arange(20,6,-1)
header = arepo_package.load_snapshot_header(f'{Filepath}/{Box}/{SimPaths[0][0]}/{outputpath}',desired_redshifts[0])
a = 1/(1 + desired_redshifts)
h = 0.6774
DMMass = header['MassTable'][1] * 1e10/h # Msun

# For density profiles:
rmin = 0.1    #the radius at which to start calculating the density
rmax = 50     #the radius at which to stop calculating the density
sphere_samples = 200 #number of measurement samples for each radius
radial_samples = 50  #number of radii to measure density 
DesNgb = 32   #how many particles to use when estimating density
all_r = np.logspace(np.log10(rmin), np.log10(rmax), radial_samples)

print('Loading masses and positions ...',flush = True)

with h5py.File(f'output/Constrained/rho_prof_p{Ptype}_MMH.hdf5', 'w') as f:

    f.create_dataset('Radii',data=all_r)
    
    for ii in range(len(SimPaths)):
    
        for sim in SimPaths[ii]:
            
            basePath = f'{Filepath}/{Box}/{sim}/{outputpath}/'

            # For all redshift density profiles for this sim:
            dens_z = []
    
            for iii in range(len(desired_redshifts)):
    
                PartPos = get_particle_property_LTU(basePath,'Coordinates', p_type=Ptype,
                                                    desired_redshift = desired_redshifts[iii])[0] * a[iii]/h # kpc
                
                if Ptype == 0:
                    PartMass = get_particle_property_LTU(basePath,'Masses', p_type=Ptype, 
                                                    desired_redshift = desired_redshifts[iii])[0] * 1e10/h # Msun
                elif Ptype == 1:
                    PartMass = DMMass * np.ones(len(PartPos))
    
                # Selecting subhalo with most massive BH
                # BHMasses,o=arepo_package.get_subhalo_property(basePath,'SubhaloBHMass',desired_redshifts[iii])
                # target = np.argmax(BHMasses)
    
                # Selecting most massive halo, hopefully has a BH
                HaloMasses,o=arepo_package.get_group_property(basePath,'GroupMass',desired_redshifts[iii])
                target = np.argmax(HaloMasses)
    
                BH_Masses,o = arepo_package.get_group_property(basePath,'GroupBHMass',desired_redshifts[iii])
    
                if BH_Masses[target] == 0:
                    print(f'Warning!! {sim} sim at z = {desired_redshifts[iii]} has no BHs')

                offset = arepo_package.get_group_property(basePath,'GroupOffsetType',desired_redshifts[iii])[0][target]
                length = arepo_package.get_group_property(basePath,'GroupLenType',desired_redshifts[iii])[0][target]
                
                GalaxyGasPos = PartPos[offset[Ptype]:offset[Ptype]+length[Ptype]]
                GalaxyGasMass = PartMass[offset[Ptype]:offset[Ptype]+length[Ptype]]
                # GalaxyGasVel = GasVel[offset[0]:offset[0]+length[0]]
    
                Gal_center = np.median(GalaxyGasPos,axis=0)

                # For all the density values at this redshift
                dens = []
    
                for r in all_r:
        
                    points = fibonacci_sphere(sphere_samples, r)
    
                    density = calc_density(GalaxyGasPos - Gal_center, GalaxyGasMass, points, DesNgb)
    
                    dens.append(density)
    
                dens_z.append(dens)
    
                print(f"{sim} at z = {desired_redshifts[iii]} had {len(GalaxyGasMass)} gas particles in the subhalo with the largest BH", flush = True)

            # Creating new dataset for this sim
            f.create_dataset(f'{sim}_density', data=dens_z)
            

print('Success!!!', flush = True)

