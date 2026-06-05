#!/usr/bin/env python3
import numpy as np
import h5py
import arepo_package
from pathlib import Path
from LtU_get_property import get_particle_property_LTU

# -------------------------
# Configuration
# -------------------------
FILEPATH = Path('/project/torrey-group/jkho/LtU_accretion')
OUTPUTPATH = 'output'
BOXES = ['Low_mass_seeds', 'Zooms', 'Constrained']
REDSHIFTS = np.linspace(20, 6, 15)

HUBBLE = 0.6774
MASS_CONV = 1e10 / HUBBLE
MDOT_CONV = 1e10 / 0.978

OUTFILE = "output/lenient_bh_growth_SMHM10.hdf5"

# -------------------------
# Simulation dictionary
# -------------------------
SIM_DICT = {
    'Low_mass_seeds': [
        ['Bondi_lowmass_AGN_fewseeds_z127',
         'FF_lowmass_AGN_fewseeds_z127',
         'modFF_lowmass_AGN_fewseeds_z127'],

        ['Bondi_lowmass_AGN_stellar_lenient',
         'FF_lowmass_AGN_stellar_lenient',
         'modFF_lowmass_AGN_stellar_lenient']
    ],

    'Zooms': [
        ['Bondi_zoom_AGN','FF_zoom_AGN','modFF_zoom_AGN'],
        ['Bondi_zoom_AGN_stellar_lenient',
         'FF_zoom_AGN_stellar_lenient',
         'modFF_zoom_AGN_stellar_lenient']
    ],

    'Constrained': [
        ['Bondi_constrained_AGN_fewseeds_stellar',
         'FF_constrained_AGN_fewseeds_stellar',
         'modFF_constrained_AGN_fewseeds_stellar'],
        
        ['Bondi_constrained_AGN_stellar_lenient_HMSM10',
         'FF_constrained_AGN_stellar_lenient_HMSM10',
         'modFF_constrained_AGN_stellar_lenient_HMSM10']
    ]
}

# -------------------------
# Helper: load BH properties
# -------------------------
def load_bh(basePath, redshift):
    """
    Load BH mass, accretion rate, and cumulative quasar growth
    for the most massive BH at a given redshift.
    """
    BH_Mass = get_particle_property_LTU(
        basePath,'BH_Mass',p_type=5,desired_redshift=redshift
    )[0]

    BH_Mdot = get_particle_property_LTU(
        basePath,'BH_Mdot',p_type=5,desired_redshift=redshift
    )[0]

    BH_QM = get_particle_property_LTU(
        basePath,'BH_CumMassGrowth_QM',p_type=5,desired_redshift=redshift
    )[0]

    BH_QM_E = get_particle_property_LTU(
        basePath,'BH_CumEgyInjection_QM',p_type=5,desired_redshift=redshift
    )[0]

    ind = np.argmax(BH_Mass)

    return (
        BH_Mass[ind] * MASS_CONV,
        BH_Mdot[ind] * MDOT_CONV,
        BH_QM[ind] * MASS_CONV,
        BH_QM_E[ind]
    )

# -------------------------
# Main extraction
# -------------------------
def main():

    with h5py.File(OUTFILE, "w") as f:

        # Save global metadata
        f.attrs['HubbleParam'] = HUBBLE
        f.create_dataset("Redshifts", data=REDSHIFTS)

        for box in BOXES:

            print(f"\nProcessing box: {box}")

            box_group = f.create_group(box)

            for mode_index, sims in enumerate(SIM_DICT[box]):

                mode_name = "strict" if mode_index == 0 else "lenient"
                mode_group = box_group.create_group(mode_name)

                for model_index, sim in enumerate(sims):

                    if model_index == 0:
                        model_name = "Bondi"
                    elif model_index == 1:
                        model_name = "FreeFall"
                    else:
                        model_name = "modFreeFall"

                    print(f"  {mode_name} | {model_name}")

                    basePath = f'{FILEPATH}/{box}/{sim}/{OUTPUTPATH}/'

                    masses = np.empty(len(REDSHIFTS))
                    mdots  = np.empty(len(REDSHIFTS))
                    qm     = np.empty(len(REDSHIFTS))
                    qme    = np.empty(len(REDSHIFTS))

                    for i, z in enumerate(REDSHIFTS):
                        masses[i], mdots[i], qm[i], qme[i] = load_bh(basePath, z)

                    model_group = mode_group.create_group(model_name)
                    model_group.create_dataset("Mass", data=masses)
                    model_group.create_dataset("Mdot", data=mdots)
                    model_group.create_dataset("QM_growth", data=qm)
                    model_group.create_dataset("QM_energy", data=qme)

    print(f"\nSaved to {OUTFILE}")


# -------------------------
# Run as script
# -------------------------
if __name__ == "__main__":
    main()
