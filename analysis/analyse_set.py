import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import numpy as np
import pandas as pd
import os
import glob
import sys
from ovito.io import import_file
import freud

def get_wrapping(p, f):
    d = p.compute(f)
    positions = np.array(d.particles['Position'])
    types = np.array(d.particles['Particle Type'])
    box = freud.Box.from_matrix(d.cell[...])
    np_pos = positions[types == 2]
    mem_pos = positions[types == 1]
    aabb = freud.locality.AABBQuery(box, mem_pos)
    neighbors = aabb.query(np_pos, {'r_max': 5.5})
    # neighbors = [len(n) for n in neighbors]
    neighbors = sum(1 for _ in neighbors)
    # area= np.pi*0.5**2*len(neighbors)
    area= np.pi*0.5**2*neighbors
    spharea = 4*np.pi*4.5**2
    wrapping = area/spharea
    return wrapping

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
        sizes = sizes[sizes > 10]
        # Store number of clusters and frame number
        ncl.append(len(sizes))
        frs.append(f)

    final_wrapping = get_wrapping(p, nf)
    
    # Define pandas dataframe with number of clusters per frame and frame numbers
    data = pd.DataFrame({'frame':frs, 'ncl':ncl})

    # Save dataframe to file and return it
    data.to_csv('%s/clusters.txt'%(path))
    return data, final_wrapping

def analyse_seeds(path, rmax = 1.5):
    # Change directory to path
    r = os.chdir(path)
    # Get list of seeds
    seeds = glob.glob('sd*')
    # Initialize list to store budding frames
    buds = []
    wraps = []
    # Loop over seeds
    for seed in seeds:
        # Extract clusters from seed
        data, final_wrapping = extract(seed, rmax)
        wraps.append(final_wrapping)
        # Check if budding occurs
        if 2 in data['ncl'].values:
            # If budding occurs, store frame number
            buds.append(data[data['ncl'] == 2].iloc[0]['frame'])
        else:
            # If budding does not occur, store NaN
            buds.append(np.nan)
    # Build pandas dataframe with budding frames
    buds = pd.DataFrame({'seed':seeds,'frame':buds})
    # Save dataframe to file and return it
    buds.to_csv('budding_frames.txt')
    # Build pandas dataframe with wrapping
    wraps = pd.DataFrame({'seed':seeds,'wrapping':wraps})
    # Save dataframe to file and return it
    wraps.to_csv('wrapping.txt')
    return buds, wraps

if __name__ == '__main__':

    # Define path to parameter set (set of simulations for the same parameters but different RNG seeds)
    gpath = sys.argv[1]

    # Extract clusters from trajectory
    # data = extract(path)
    data_bud, data_wrap = analyse_seeds(gpath)

    # Print frame where budding occurs
    print()
    print()
    print(data_bud)
    print()
    print()
    print(data_wrap)
    print()
    print()
