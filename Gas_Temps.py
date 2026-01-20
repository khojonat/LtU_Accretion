import sys
import matplotlib
matplotlib.use('agg')
import numpy as np
import os
import h5py
from astropy import constants as cons
from astropy import units as u

from LtU_get_property import get_particle_property_LTU
import arepo_package as arepo_package

# ----------------------
# User parameters
# ----------------------
d = 0.140  # Zoom radius

Filepath = '/project/torrey-group/jkho/LtU_accretion'
outputpath = 'output'

Boxes = {
    'Zooms': [
        ['Bondi_zoom_AGN','Bondi_zoom_boost_noAGN_stellar','Bondi_zoom_AGN_0.1stellar','Bondi_zoom_noAGN_0.1stellar'],
        ['FF_zoom_AGN','FF_zoom_noAGN_stellar','FF_zoom_AGN_0.1stellar','FF_zoom_noAGN_0.1stellar'],
        ['modFF_zoom_AGN','modFF_zoom_noAGN_stellar','modFF_zoom_AGN_0.1stellar','modFF_zoom_noAGN_0.1stellar']
    ],
    'Constrained': [
        ['Bondi_constrained_AGN_fewseeds_stellar','Bondi_constrained_noAGN_fewseeds_boost','Bondi_constrained_AGN_fewseeds_0.1stellar','Bondi_constrained_noAGN_0.1stellar'],
        ['FF_constrained_AGN_fewseeds_stellar','FF_constrained_noAGN_fewseeds','FF_constrained_AGN_fewseeds_0.1stellar','FF_constrained_noAGN_0.1stellar'],
        ['modFF_constrained_AGN_fewseeds_stellar','modFF_constrained_noAGN_fewseeds','modFF_constrained_AGN_fewseeds_0.1stellar','modFF_constrained_noAGN_0.1stellar']
    ],
    'Low_mass_seeds': [
        ['Bondi_lowmass_AGN_fewseeds_z127','Bondi_lowmass_noAGN_fewseeds_z127','Bondi_lowmass_AGN_fewseeds_lowstellar_z127','Bondi_lowmass_noAGN_fewseeds_lowstellar_z127'],
        ['FF_lowmass_AGN_fewseeds_z127','FF_lowmass_noAGN_fewseeds_z127','FF_lowmass_AGN_fewseeds_lowstellar_z127','FF_lowmass_noAGN_fewseeds_lowstellar_z127'],
        ['modFF_lowmass_AGN_fewseeds_z127','modFF_lowmass_noAGN_fewseeds_z127','modFF_lowmass_AGN_fewseeds_lowstellar_z127','modFF_lowmass_noAGN_fewseeds_lowstellar_z127']
    ]
}

redshifts = np.arange(6, 17)
a = 1 / (1 + redshifts)
h = 0.6774

# ----------------------
# Main loop
# ----------------------
for box in Boxes:

    h5name = f'output/Temps/{box}_gas_BH_properties.hdf5'
    print(f'Writing {h5name}')

    with h5py.File(h5name, 'w') as fbox:
        fbox.attrs['redshifts'] = redshifts

        for sim_group in Boxes[box]:
            for sim in sim_group:

                print(f'  Processing {sim}')
                basePath = f'{Filepath}/{box}/{sim}/{outputpath}/'

                sim_grp = fbox.create_group(sim)

                for i, z in enumerate(redshifts):

                    zgrp = sim_grp.create_group(f'z_{z}')

                    # ----------------------
                    # Load BH properties
                    # ----------------------
                    BH_Mass = get_particle_property_LTU(
                        basePath, 'BH_Mass', p_type=5, desired_redshift=z
                    )[0] * 1e10 / h

                    BH_Hsml = get_particle_property_LTU(
                        basePath, 'BH_Hsml', p_type=5, desired_redshift=z
                    )[0] * a[i] / h

                    BH_pos = get_particle_property_LTU(
                        basePath, 'Coordinates', p_type=5, desired_redshift=z
                    )[0] * a[i] / h

                    if len(BH_Mass) == 0 or np.all(np.isnan(BH_Mass)):
                        zgrp.attrs['empty'] = True
                        continue

                    most_massive_BH = np.argmax(BH_Mass)
                    hsml = BH_Hsml[most_massive_BH]
                    MMBH_coord = BH_pos[most_massive_BH]

                    # ----------------------
                    # Load gas properties
                    # ----------------------
                    Gas_pos = get_particle_property_LTU(
                        basePath, 'Coordinates', p_type=0, desired_redshift=z
                    )[0] * a[i] / h

                    Gas_Mass = get_particle_property_LTU(
                        basePath, 'Masses', p_type=0, desired_redshift=z
                    )[0] * 1e10 / h

                    Gas_SFR = get_particle_property_LTU(
                        basePath, 'StarFormationRate', p_type=0, desired_redshift=z
                    )[0]

                    U = get_particle_property_LTU(
                        basePath, 'InternalEnergy', p_type=0, desired_redshift=z
                    )[0] * u.km**2 / u.s**2

                    x_e = get_particle_property_LTU(
                        basePath, 'ElectronAbundance', p_type=0, desired_redshift=z
                    )[0]

                    gas_in_hsml = np.linalg.norm(Gas_pos - MMBH_coord, axis=1) < hsml

                    gas_mass_in_hsml = Gas_Mass[gas_in_hsml]
                    gas_pos_in_hsml = Gas_pos[gas_in_hsml]
                    gas_sfr_in_hsml = Gas_SFR[gas_in_hsml]
                    gas_u_in_hsml = U[gas_in_hsml]
                    gas_xe_in_hsml = x_e[gas_in_hsml]

                    # ----------------------
                    # Temperature
                    # ----------------------
                    XH = 0.76
                    gamma = 5 / 3
                    mu = 4 / (1 + 3*XH + 4*XH*gas_xe_in_hsml) * cons.m_p
                    T = (gamma - 1) * gas_u_in_hsml / cons.k_B * mu
                    T = T.to(u.K).value

                    # ----------------------
                    # Write datasets
                    # ----------------------
                    zgrp.create_dataset('gas_masses', data=gas_mass_in_hsml)
                    zgrp.create_dataset('gas_coords', data=gas_pos_in_hsml)
                    zgrp.create_dataset('gas_sfr_in_hsml', data=gas_sfr_in_hsml)
                    zgrp.create_dataset('gas_temps', data=T)
                    zgrp.create_dataset('bh_coords', data=MMBH_coord)

                    zgrp.attrs['hsml'] = hsml
                    zgrp.attrs['bh_mass'] = BH_Mass[most_massive_BH]
