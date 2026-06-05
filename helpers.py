import sys, os
import numpy as np
import matplotlib.pyplot as plt

# sys.path.insert(1,'~/Python/torreylabtoolsPy3')
sys.path.append("/home/yja6qa/LtU_Accretion/torreylabtools")

from util import calc_hsml # Formerly torreylabtoolsPy3
from visualization import contour_makepic as makepic

plt.rcParams['text.usetex']        = True
plt.rcParams['font.family']        = 'serif'
plt.rcParams['font.size']          = 20

fs_og = 24
plt.rcParams['font.size'] = fs_og
plt.rcParams['axes.linewidth']  = 2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.minor.visible'] = 'true'
plt.rcParams['ytick.minor.visible'] = 'true'
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.0
plt.rcParams['ytick.minor.width'] = 1.0
plt.rcParams['xtick.major.size']  = 7.5
plt.rcParams['ytick.major.size']  = 7.5
plt.rcParams['xtick.minor.size']  = 3.5
plt.rcParams['ytick.minor.size']  = 3.5
plt.rcParams['xtick.top']   = True
plt.rcParams['ytick.right'] = True


def check_save_dir(save_dir):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

def Curti_MZR():
    
    ms = 10**np.arange(8.0, 11.5, 0.1)
    
    Z0 = 8.793
    gamma = 0.28
    beta  = 1.2
    M0    = 1.00E+01 **(10.02)
    
    Zs = Z0 - gamma/beta *np.log10( 1 + (ms/M0)**(-beta) )
    ms = np.log10(ms)
    return ms, Zs


def calc_incl(pos0, vel0, m0, ri, ro):
    rpos = np.sqrt(pos0[:,0]**2.000E+00 +
                   pos0[:,1]**2.000E+00 +
                   pos0[:,2]**2.000E+00 )
    rpos = rpos[~np.isnan(rpos)]
    idx  = (rpos > ri) & (rpos < ro)
    pos  = pos0[idx]
    vel  = vel0[idx]
    m    =   m0[idx]
        
    hl = np.cross(pos, vel)
    L  = np.array([np.multiply(m, hl[:,0]),
                   np.multiply(m, hl[:,1]),
                   np.multiply(m, hl[:,2])])
    L  = np.transpose(L)
    L  = np.array([np.sum(L[:,0]),
                   np.sum(L[:,1]),
                   np.sum(L[:,2])])
    Lmag  = np.sqrt(L[0]**2.000E+00 +
                    L[1]**2.000E+00 +
                    L[2]**2.000E+00 )
    Lhat  = L / Lmag
    incl  = np.array([np.arccos(Lhat[2]), np.arctan2(Lhat[1], Lhat[0])])
    # incl *= 1.800E+02 / np.pi
    if   incl[1]  < 0.000E+00:
         incl[1] += 2*np.pi
    elif incl[1]  > 2*np.pi:
         incl[1] -= 2*np.pi
    return incl

def trans(arr0, incl0):
    arr      = np.copy( arr0)
    incl     = np.copy(incl0)
    # deg2rad  = np.pi / 1.800E+02
    # incl    *= deg2rad
    arr[:,0] = -arr0[:,2] * np.sin(incl[0]) + (arr0[:,0] * np.cos(incl[1]) + arr0[:,1] * np.sin(incl[1])) * np.cos(incl[0])
    arr[:,1] = -arr0[:,0] * np.sin(incl[1]) + (arr0[:,1] * np.cos(incl[1])                                                )
    arr[:,2] =  arr0[:,2] * np.cos(incl[0]) + (arr0[:,0] * np.cos(incl[1]) + arr0[:,1] * np.sin(incl[1])) * np.sin(incl[0])
    del incl
    return arr

def center(pos0, centpos, boxsize = None):
    pos       = np.copy(pos0)
    pos[:,0] -= centpos[0]
    pos[:,1] -= centpos[1]
    pos[:,2] -= centpos[2]
    if (boxsize != None):
        pos[:,0][pos[:,0] < (-boxsize / 2.000E+00)] += boxsize
        pos[:,0][pos[:,0] > ( boxsize / 2.000E+00)] -= boxsize
        pos[:,1][pos[:,1] < (-boxsize / 2.000E+00)] += boxsize
        pos[:,1][pos[:,1] > ( boxsize / 2.000E+00)] -= boxsize
        pos[:,2][pos[:,2] < (-boxsize / 2.000E+00)] += boxsize
        pos[:,2][pos[:,2] > ( boxsize / 2.000E+00)] -= boxsize
    return pos

def make_map(pos, mass):
    
    hsml = calc_hsml.get_particle_hsml( pos[:,0], pos[:,1], pos[:,2], DesNgb=32  )

    n_pixels = 720
    massmap,image = makepic.contour_makepic( pos[:,0], pos[:,1], pos[:,2], hsml, mass,
        xlen = rmax,
        pixels = n_pixels, set_aspect_ratio = 1.0,
        set_maxden = 1.0e10, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
        set_dynrng = 1.0e4  )
    
    return massmap

def center_and_box_wrap(pos, mass, center, boxsize, h):
    for ijk in range(3):
        pos[:,ijk] -= center[ijk]
        pos[ pos[:,ijk] >  1.0*boxsize/2.0 , ijk ] -= boxsize
        pos[ pos[:,ijk] < -1.0*boxsize/2.0 , ijk ] += boxsize

    rad = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
    rad = rad * 1000

    return rad, mass * 1.00E+10 / h
