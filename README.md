# BSS2024_LAMMPS_Task

This repository contains bash and Python scripts to:
- Download and compile the February 2016 stable release of LAMMPS
- Setup and run Molecular Dynamics simulations of a patchy nanoparticle being wrapped by a model lipid membrane in different scenarios
- Analyse and extract information from the simulation trajectories

[LAMMPS_installation](LAMMPS_installation) contains the [bash script](LAMMPS_installation/build_lammps.sh) to download and install LAMMPS, which will be placed in a folder [lammps](LAMMPS_installation/lammps) (not part of this repository) after installation.

## How to do things

In order to run the simulations on your laptop you will need to first download and compile the right version of LAMMPS. To do so, from [LAMMPS_installation](LAMMPS_installation) type `bash build_lammps.sh`

Once LAMMPS is properly installed you can then set up your simulations by:

1. Create a directory where you will be running all the replicas (different simulations with different RNG seeds for a given set of parameters) for a given scenario. For example `~/scenario1` by running `mkdir ~/scenario1`
2. Copy the [Inputs](Inputs) folder into that directory. In this example you can just run `cp -r ./Inputs ~/scenario1` from [here](./)
3. Edit the [InputData.txt](Inputs/InputData.txt) according to the scenario you want to explore.
    [InputData.txt](Inputs/InputData.txt) contains the details of the binding sites on the nanoparticle. Besides the header, the file contains one line per binding site specifying the position (coordinates) and binding interaction strength.
   Each line must contain the following columns:
     1) ID of the binding site (increasing from 1)
     2) type of the binding site (increasing from 3 as there are already two types of beads in the system besides the binding sites: membrane and nanoparticle)
     3) Epsilon: binding potential strength (in kT units)
     4) x,y,z coordinates of the binding site (relative to the nanoparticle's center)
5. Go to your set directory (‘~/scenario1‘) in this example. You can do this by running `cd ~/scenario1`
6. Run the [Python setup script](MakeLammpsInput.py) from the set directory giving it a seed for the random number generator from the command line: `python3 ${path_to_this_repository}/MakeLammpsInput.py 12345` for example, to set it up with seed `12345`
7. Go to the seed folder. Run `cd sd12345`
8. Run the LAMMPS executable from the seed folder: `${path_to_this_repository}/LAMMPS_installation/lammps/src/lmp_serial -in in.local`
9. Repeat this procedure for all the different replicas you might want to run in order to get proper statistics.
10. Run the [analysis](analysis) from this repository by typing `python3 analyse_trajectory.py ~/scenario1` for this particular example from the [analysis](analysis) folder. This will output a .txt file with the budding times for each seed in `scenario1`.
