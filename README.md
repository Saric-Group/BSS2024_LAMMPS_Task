# BSS2024_LAMMPS_Task

This repository contains bash and Python scripts to:
- Download and compile the February 2016 stable release of LAMMPS
- Setup and run Molecular Dynamics simulations of a patchy nanoparticle being wrapped by a model lipid membrane in different scenarios
- Analyse and extract information from the simulation trajectories

[LAMMPS_installation](LAMMPS_installation) contains the [bash script](LAMMPS_installation/build_lammps.sh) to download and install LAMMPS, which will be placed in a folder [lammps](LAMMPS_installation/lammps) (not part of this repository) after installation.

## How to do things

In order to run the simulations on your laptop you will need to first download and compile the right version of LAMMPS. To do so, from [LAMMPS_installation](LAMMPS_installation) type ‘bash build_lammps.sh‘

Once LAMMPS is properly installed you can then set up your ...
