# LtU_Accretion
Repo for code related to my LtU ([Learning the Universe](https://learning-the-universe.org/)) accretion project testing Bondi, ff, and modff accretion models in different environments and with different AGN and stellar feedback strengths ([paper for ff, modff models](https://www.aanda.org/articles/aa/full_html/2025/08/aa54174-25/aa54174-25.html)). The arepo_package file (used to help load in simulation data) depends on having [illustris python](https://github.com/illustristng/illustris_python) installed in the top level directory.

The code used to generate the 7 figures in my paper can be found in LtU_Accretion_Paper_Figs.ipynb, and the corresponding figures in Plots.

Supplementary_files and Supplementary_Plots contain additional scripts, notebooks, and images I produced and used throughout the duration of this project that did not make it in the paper, and are by no means organized: observer be warned!

The output folder contains data pulled from the various simulations ran for this project (which, due to size, cannot be stored directly here). The notebooks make use of this data to quickly produce the various paper plots, as this is much faster than reloading the data from the simulations every time.

Some of the scripts used to pull and store this simulation data exist in the top level directory, with corresponding slurm submit scripts living in the slurm_scripts folder.

Currently, this paper is submitted (or soon to be submitted) to ApJ.
