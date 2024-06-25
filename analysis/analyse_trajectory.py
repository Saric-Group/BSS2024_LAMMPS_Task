import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import numpy as np
import pandas as pd
import os
import glob
import sys
from ovito.io import import_file
import freud


def extract(path, rmax = 1.5):
    # Load trajectory
    p = import_file('%s/output.xyz'%(path))
    # Get number of frames
    nf = p.source.num_frames

    # Initialize cluster object
    cl = freud.cluster.Cluster()
    # Initialize cluster properties object
    cl_props = freud.cluster.ClusterProperties()

    # Initialize list to store number of clusters per frame and frame numbers
    ncl = []
    frs = []
    
    # Loop over frames
    for f in range(nf):
        # Extract data from pipeline frame and keep only membrane positions
        d = p.compute(f)
        pos = np.array(d.particles['Position'])
        typ = np.array(d.particles['Particle Type'])
        mem = pos[typ == 1]
        # Define the box in Freud from the dataframe
        box = freud.Box.from_matrix(d.cell[...])
        # Compute clusters and cluster properties using Freud
        cl.compute((box, mem), neighbors={'r_max': rmax})
        cl_props.compute((box, mem), cl.cluster_idx)
        # Extract cluster sizes
        sizes = cl_props.sizes
        # Store number of clusters and frame number
        ncl.append(len(sizes))
        frs.append(f)
    
    # Define pandas dataframe with number of clusters per frame and frame numbers
    data = pd.DataFrame({'frame':frs, 'ncl':ncl})

    # Save dataframe to file and return it
    data.to_csv('%s/clusters.txt'%(path))
    return data

if __name__ == '__main__':

    # Define path to trajectory from command line argument
    path = sys.argv[1]

    # Extract clusters from trajectory
    data = extract(path)

    # Print frame where budding occurs
    print("Budding occurs at frame %d"%(data[data['ncl'] == 2].iloc[0]['frame']))
