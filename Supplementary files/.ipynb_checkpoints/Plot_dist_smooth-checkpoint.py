import h5py 
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from illustris_python import illustris_python as il
import astropy.units as u
import astropy.constants as cons
from LtU_get_property import get_particle_property_LTU
import arepo_package as arepo_package
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm

torreylabtools_path = '/sfs/gpfs/tardis/home/yja6qa/FIRE_MW_suite/torreylabtools'
sys.path.insert(0, str(torreylabtools_path))
from torreylabtools import kdtree_smoothing as kd


# Constrained

basePath = f'/project/torrey-group/jkho/LtU_accretion/Constrained/Bondi_constrained_AGN_fewseeds_stellar/output/'
z = 6
Ptype = 1


header = arepo_package.load_snapshot_header(basePath,z)
print(header.keys())
BoxSize = header['BoxSize']
DM_Mass = header['MassTable'][1]
h = 0.6774
DM_Pos,z = get_particle_property_LTU(basePath,'Coordinates',p_type=1, desired_redshift = z)

a = 1/(1+z)
DM_Pos *= a/h
DM_Pos_Constrained = DM_Pos - np.max(DM_Pos)/2 # Centering
DM_Masses_Constrained = np.ones(len(DM_Pos)) * DM_Mass * 1e10/h

HaloMasses,o=arepo_package.get_group_property(basePath,'GroupMass',z)
HaloPos,o=arepo_package.get_group_property(basePath,'GroupPos',z)
HaloRvir,o=arepo_package.get_group_property(basePath,'Group_R_Crit200',z) # Using same radius for all halos
target = np.argmax(HaloMasses)

Constrained_rvir = HaloRvir[target] * a/h
Constrained_Pos = HaloPos[target] * a/h - np.max(DM_Pos)/2
Constrained_target_pos = DM_Pos_Constrained[np.linalg.norm(DM_Pos_Constrained - Constrained_Pos,axis = 1) < Constrained_rvir]
Constrained_target_mass = DM_Masses_Constrained[np.linalg.norm(DM_Pos_Constrained - Constrained_Pos,axis = 1) < Constrained_rvir] 


# Zoom

basePath = f'/project/torrey-group/jkho/LtU_accretion/Zooms/Bondi_zoom_AGN/output/'
z = 6

header = arepo_package.load_snapshot_header(basePath,z)
print(header.keys())
BoxSize = header['BoxSize']
DM_Mass = header['MassTable'][1]
h = 0.6774
DM_Pos,z = get_particle_property_LTU(basePath,'Coordinates',p_type=1, desired_redshift = z)

a = 1/(1+z)
DM_Pos *= a/h
DM_Pos_Zooms = DM_Pos - np.max(DM_Pos)/2 # Centering
DM_Masses_Zooms = np.ones(len(DM_Pos)) * DM_Mass * 1e10/h

HaloMasses,o=arepo_package.get_group_property(basePath,'GroupMass',z)
HaloPos,o=arepo_package.get_group_property(basePath,'GroupPos',z)
HaloRvir,o=arepo_package.get_group_property(basePath,'Group_R_Crit200',z) # Using same radius for all halos
target = np.argmax(HaloMasses)

Zooms_rvir = HaloRvir[target] * a/h
Zooms_Pos = HaloPos[target] * a/h - np.max(DM_Pos)/2
Zooms_target_pos = DM_Pos_Zooms[np.linalg.norm(DM_Pos_Zooms - Zooms_Pos,axis = 1) < Constrained_rvir]
Zooms_target_mass = DM_Masses_Zooms[np.linalg.norm(DM_Pos_Zooms - Zooms_Pos,axis = 1) < Constrained_rvir]

# Small Uniform

basePath = f'/project/torrey-group/jkho/LtU_accretion/Low_mass_seeds/Bondi_lowmass_AGN_fewseeds_z127/output/'
z = 6

header = arepo_package.load_snapshot_header(basePath,z)
print(header.keys())
BoxSize = header['BoxSize']
DM_Mass = header['MassTable'][1]
h = 0.6774
DM_Pos,z = get_particle_property_LTU(basePath,'Coordinates',p_type=1, desired_redshift = z)

a = 1/(1+z)
DM_Pos *= a/h
DM_Pos_Small = DM_Pos - np.max(DM_Pos)/2 # Centering
DM_Masses_Small = np.ones(len(DM_Pos)) * DM_Mass * 1e10/h

HaloMasses,o=arepo_package.get_group_property(basePath,'GroupMass',z)
HaloPos,o=arepo_package.get_group_property(basePath,'GroupPos',z)
HaloRvir,o=arepo_package.get_group_property(basePath,'Group_R_Crit200',z) # Using same radius for all halos
target = np.argmax(HaloMasses)

Small_rvir = HaloRvir[target] * a/h
Small_Pos = HaloPos[target] * a/h - np.max(DM_Pos)/2
Small_target_pos = DM_Pos_Small[np.linalg.norm(DM_Pos_Small - Small_Pos,axis = 1) < Constrained_rvir]
Small_target_mass = DM_Masses_Small[np.linalg.norm(DM_Pos_Small - Small_Pos,axis = 1) < Constrained_rvir]



# =======================
# Toggle: 1 = same colorbar scale, 0 = independent scales
same_scale = 1
cmap = 'viridis'
# =======================

fig, axs = plt.subplots(1, 3, figsize=(12, 3.25))

# ---------------------------------------------------
# Compute global vmin/vmax if same scale is requested
# ---------------------------------------------------
if same_scale:
    all_weights = np.concatenate([
        Constrained_target_mass,
        Zooms_target_mass,
        Small_target_mass
    ])
    vmin = 1e4 # np.min(all_weights[all_weights > 0])
    vmax = 1e10 # np.max(all_weights)
else:
    vmin = vmax = None

bins = 250

# ---------
# Panel 1
# ---------

image = kd.simple_makepic(
    Constrained_target_pos, Constrained_target_mass, image_size = np.max(Constrained_target_pos), 
    pixel_size=np.max(Constrained_target_pos)/1000
)

axs[0].imshow(image,cmap=cmap,norm=LogNorm(vmin=vmin, vmax=vmax))

p = Circle(
    (Constrained_Pos[0], Constrained_Pos[1]),
    Constrained_rvir,
    fill=False,
    edgecolor='red',
    linewidth=1.5
)
axs[0].add_artist(p)
axs[0].axis('off')
# fig.colorbar(h1[3], ax=axs[0])

# ---------
# Panel 2
# ---------

image = kd.simple_makepic(
    Zooms_target_pos, Zooms_target_mass, image_size = np.max(Zooms_target_pos), 
    pixel_size=np.max(Zooms_target_pos)/1000
)

axs[1].imshow(image,cmap=cmap,norm=LogNorm(vmin=vmin, vmax=vmax))

p = Circle(
    (Zooms_Pos[0], Zooms_Pos[1]),
    Zooms_rvir,
    fill=False,
    edgecolor='red',
    linewidth=1.5
)
axs[1].add_artist(p)
axs[1].axis('off')
# fig.colorbar(h2[3], ax=axs[1])

# ---------
# Panel 3
# ---------

image = kd.simple_makepic(
    Small_target_pos, Small_target_mass, image_size = np.max(Small_target_pos), 
    pixel_size=np.max(Small_target_pos)/1000
)

axs[2].imshow(image,cmap=cmap,norm=LogNorm(vmin=vmin, vmax=vmax))

p = Circle(
    (Small_Pos[0], Small_Pos[1]),
    Small_rvir,
    fill=False,
    edgecolor='red',
    linewidth=1.5
)
axs[2].add_artist(p)
axs[2].axis('off')

cbar = fig.colorbar(im, ax=axs)
cbar.set_label(r'$\rm M_\odot~/~kpc$', rotation=270, labelpad=15, fontsize=15)

# Titles
axs[0].set_title('Large Halo ICs', size=15)
axs[1].set_title('Medium Halo ICs', size=15)
axs[2].set_title('Small Halo ICs', size=15)

# fig.tight_layout()
fig.savefig('Plots/Halo_distributions_smoothed.png')

