import numpy as np
from tqdm import tqdm
from scipy.spatial import KDTree, cKDTree

def gaussian_kernel(q, h, truncate_at=10):
    W = np.zeros_like(q)
    mask = q <= truncate_at
    W[mask] = np.exp(-q[mask]**2)
    normalization = 1 / ((np.pi) ** 1.5 * h ** 3)
    W *= normalization
    return W
    
def tophat_kernel(q, h):
    W = np.zeros_like(q)
    W[q <= 1] = 1
    norm = 3 / (4 * np.pi * h**3)
    W *= norm
    return W
    

def cubic_spline_kernel(u, h):
    u = np.asarray(u)
    norm = 8 / (np.pi * h**3)
    W = np.zeros_like(u)

    mask1 = u < 0.5
    mask2 = (u >= 0.5) & (u < 1.0)

    W[mask1] = norm * (1 - 6 * u[mask1]**2 + 6 * u[mask1]**3)
    W[mask2] = norm * (2 * (1 - u[mask2])**3)

    return W

def get_neighbors(tree, grid_positions, DesNgb):
    dist_matrix, idx_matrix = tree.query(grid_positions, k=DesNgb)
    return dist_matrix, idx_matrix

def contour_makepic(
    particle_xyz, weight,
    image_size=10, pixel_size=0.1
):
    
    Npix = int(np.ceil(2 * image_size / pixel_size))
    x = np.linspace(-image_size, image_size, Npix)
    y = np.linspace(-image_size, image_size, Npix)
    z = np.linspace(-image_size, image_size, Npix)
    
    Npix = len(x)
    grid_shape = (Npix, Npix, Npix)
    grid_points = np.array(np.meshgrid(x, y, z, indexing='ij')).reshape(3, -1).T
    
    image = get_map(
        particle_xyz[:, 0], particle_xyz[:, 1], particle_xyz[:, 2], weight,
        grid_points[:,0], grid_points[:, 1], grid_points[:, 2], Npix
    )
    
    return image


def get_map(particle_x, particle_y, particle_z, weight, 
            grid_x, grid_y, grid_z, Npix,
            DesNgb=32,Hmax=33.0, kernel='cubic-spline', two_d=False):
    
    N_grid = len(grid_x)
    N_mass = len(particle_z)
    
    # Prepare KDTree of gas particles
    particle_pos = np.vstack((particle_x, particle_y, particle_z)).T
    weight = np.array(weight)
    
    tree = cKDTree(particle_pos)

    # Store density results
    dens = np.zeros(N_grid)
    
    grid_positions = np.vstack((grid_x, grid_y, grid_z)).T

    kernel_map = {
        'gaussian': gaussian_kernel,
        'cubic-spline': cubic_spline_kernel,
        'tophat': tophat_kernel
    }

    if kernel not in kernel_map:
        raise ValueError("Choose either 'gaussian', 'cubic-spline', or 'tophat'")

    kernel_fn = kernel_map[kernel]
    
    dists, idxs = get_neighbors(tree, grid_positions, DesNgb)

    h = 1.04 * dists[:, -1]
    print("h",np.min(h), np.max(h), np.std(h))
    
    for i in tqdm(range(N_grid)):
        grid_pos = grid_positions[i]
        h_guess  = h[i]
        h2       = h_guess ** 2
        hinv     = 1.0 / h_guess

        neighbor_idxs = idxs[i]
        particle_pos_neighbors  = particle_pos[neighbor_idxs]
        weight_neighbors = weight[neighbor_idxs]

        displacements = grid_pos - particle_pos_neighbors
        r2 = np.sum(displacements**2, axis=1)

        mask = r2 < h2
        if not np.any(mask):
            dens[i] = 0.0
            continue

        r2 = r2[mask]
        weight_neighbors = weight_neighbors[mask]

        r = np.sqrt(r2)
        u = r * hinv
        W = kernel_fn(u, h_guess)

        dens[i] = np.sum(weight_neighbors * W)
    
    if two_d: 
        return dens
    else: ## default mode
        dens = dens.reshape(Npix, Npix, Npix)
        y = np.unique(grid_y)  # should be length Npix # Changed from z to y
        dens = np.trapz(dens, x=y, axis=1) # Changed x to y from z and axis to 1 from 2
        # dens.reshape(Npix, Npix)
        
    return dens

def simple_makepic(
    particle_xyz, weight,
    image_size=10, pixel_size=0.1
):
    Npix = int(np.ceil(2 * image_size / pixel_size))
    x = np.linspace(-image_size, image_size, Npix)
    z = np.linspace(-image_size, image_size, Npix)
    
    xx, zz = np.meshgrid(x, z, indexing='ij')
    yy     = np.zeros_like(xx)
    grid_points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    particle_xyz_copy = particle_xyz.copy()
    
    particle_xyz_copy[:, 1] = np.zeros_like(particle_xyz_copy[:, 0]) ## compress everything along y axis
    
    image = simple_get_map(
        particle_xyz_copy[:, 0], particle_xyz_copy[:, 2], weight,
        grid_points[:,0], grid_points[:,2], Npix,
    )
    
    return image
    
def simple_get_map(particle_x, particle_y, weight,
                   grid_x, grid_y, Npix,
                   DesNgb=32, Hmax=33.0, kernel='cubic-spline'):
    particle_z = np.zeros_like(particle_x)
    grid_z     = np.zeros_like(grid_x)
    
    N_grid = len(grid_x)
    N_mass = len(particle_z)
    
    dens = get_map(
        particle_x, particle_y, particle_z, weight,
        grid_x, grid_y, grid_z, Npix,
        DesNgb=DesNgb, Hmax=Hmax, kernel=kernel, two_d=True
    )
    
    dens = dens.reshape(Npix, Npix)
    
    return dens