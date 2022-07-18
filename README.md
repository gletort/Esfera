# Esfera
Simulations of nuclear droplets diffusion + fusion.

3D agent-based simulations implemented in C++ used in our paper [Al Jord et al., bioRxiv] to simulate diffusion of nuclear speckles with possibility of fusion upon contact. Obstacles can be added to mimick nuclear confinement and can be condensed around the nucleolus to reproduce DNA condensation in NSN to SN transition in mouse oocytes.

## Parameters
The template configuration files to reproduce the simulations presented in our paper are available in the [parameters folder](./parameters).  

## Running simulations
`Esfera` is a C++ software that can be runned through the command line. A Makefile is available in this repository to compile it and it can be runned with the corresponding configuration file (an `xml` file) that contains the parameters of the simulation. The output files (`txt` files) containing the position of each sphere at the output times are created in a `output` folder in the folder where the simulation is runned.

## Visualization
To visualize the output of the simulations, we used [Paraview](https://www.paraview.org/). `Paraview` configuration files that allow to display directly the output files are given in this repository in the folder [paraviewfiles](./paraviewfiles).

## References
**Cytoplasmic forces functionally reorganize nuclear condensates in oocytes**
<br>Adel Al Jord, Gaelle Letort, Adrien Eichmuller, Soline Chanet, Jean-René Huynh, Nir S. Gov, Raphael Voituriez, Marie-Emilie Terret, Marie-Helene Verlhac
<br>[bioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.434387v2) doi:10.1101/2021.03.15.434387

If you use this code or an adaptation of it, please cite our paper, or use this DOI:
[![DOI](https://zenodo.org/badge/429844430.svg)](https://zenodo.org/badge/latestdoi/429844430)

## Remarks
This code was developped in the [Terret/Verlhac](https://www.college-de-france.fr/site/en-cirb/Terret-Verlhac.htm "website team") team at the Centre de Recherche Interdisciplinaire en Biologie in the College de France.
`Esfera` is freely available open-source under the GNU GPL v3.0 License (see License file). 
