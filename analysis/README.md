Contains analysis scripts (Python) to extract budding times from LAMMPS simulations.

### Output

The script is setup to run on a folder containing several replicas (different random number generator seeds) for the same scenario (number of binding sites, spatial distribution on the nanoparticle and binding strength values for each). 

It will output two .txt files. The first, `budding_frames.txt`, contains, for each seed (row), the frame of the simulation at which full budding (the nanoparticle is completely wrapped by the membrane and detaches from the flat surface) happens. If it does not happen within the trajectory it will output `NaN` for this seed. The second, `wrapping.txt`, contains, for each seed (row), the final wrapping degree. This is defined as the fraction of the nanoparticle's surface that is covered by the membrane.


### Installation and running

Will need some specific Python modules to run smoothly.

In order to set things up please create a virtual environent named BSS24 with conda by running

    conda create --name BSS24 python=3.9

Then activate it

    conda activate BSS24 

And install the required packages

    pip3 install ovito
    pip3 install freud-analysis
    pip3 install pandas
    

