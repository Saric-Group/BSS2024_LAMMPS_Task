import numpy as np
import random
import math
import sys
import os
import pandas



def main():
    
    InputPath = os.path.join(str(os.getcwd()))
    
    print("Reading data ... ")
    f = pandas.read_csv(InputPath+"/Inputs/InputData.txt", header = [0],sep = ' ')
    IDs = f.ID
    types = f.type
    epslions = f.Epsilon # this is the interaction strength of patches (could be changed to have different values for different patches)
    x,y,z = f.x,f.y,f.z    
    NumPatches= len(IDs)
    
    
    #f2 = open(InputPath+"/Inputs/InputData2.txt",'w')
    #f2.write('ID type Epsilon x y z\n')
    #for i in range(len(IDs)):
    #    f2.write(str(IDs[i])+' '+str(types[i])+' '+str(epslions[i])+' '+str(x[i])+' '+str(y[i])+' '+str(float(z[i]-6.5))+'\n')
    #f2.close()
    
    
    
    
    seed = int(sys.argv[1])

    r = os.mkdir('%s/sd%d'%(InputPath,seed))
        
    print("Writing input file ... ")
    
    f_in = open(InputPath + '/sd%s/in.local'%(seed),'w') # this is our input files directory
    f_in.write('''
    
# set up our simulation environment
dimension		3
units			lj
atom_style		hybrid sphere dipole
boundary		p p p\n''')
    f_in.write('read_data			"'+InputPath + '/sd%d/data"\n'%(seed)) # here we specify where the data is for the initial configurations of our particles
    f_in.write('''
    #group particles according to their types

group		mem			type 1
group		vehicle		type 2
group		ligand 		type >= 3
group		np      	type >= 2


#give our particle a small kick towards the membrane

velocity	np	set 0 0 -2

#membrane parameters (see Yuan 2011)
variable	rc			equal	2.6
variable	rmin		equal	1.122462
variable	rclus		equal	2.0
variable	mu			equal	3
variable	zeta		equal	4
variable	eps			equal	4.34
variable	sigma		equal	1.00
variable	theta0_11	equal	0
variable 	memsize equal "count(mem)"

# nanoparticle parameters

variable	peps		equal	2.2

# set up additional variables to be computed during the simulation (see lammps wiki for full explainations)

#compute 		1 		all 		pair lj/cut epair
compute 		cls 	all 		cluster/atom ${rclus}
compute 		ct 		all 		temp/sphere
compute_modify 	ct 		extra 		${memsize}

# set up the pair style for the membrane

pair_style	hybrid		membrane ${rc}	lj/cut 5.04

# membrane-nanoparticle interactions, each is set separately so any distribution of ligand strengths is possible

pair_coeff		*	*	lj/cut		0.0				0.0	0.0
pair_coeff		1	2	lj/cut		100				4.0	4.45
''')
    for n in range(NumPatches):
        f_in.write('pair_coeff		1	'+str(types[n])+'	lj/cut		'+str(epslions[n])+'			1	1.8\n')
    
    f_in.write('''#we set the interaction to zero at its cutoff distance, otherwise we will have a discontinuity

pair_modify		pair 	lj/cut	shift yes

#membrane-membrane interactions (these can be changed by messing with the variables at the top)

pair_coeff		1	1	membrane ${eps} ${sigma} ${rmin} ${rc} ${zeta} ${mu} ${theta0_11}

neigh_modify	every 1	delay 1	exclude group np np

# we set up the integration parameters
''')
    
    f_in.write('fix			fLANG		all		langevin 1.0 1.0 1.0 '+str(seed)+' zero yes omega yes\n') ## reading seed!
    f_in.write('''
fix			fNPH		mem		nph/sphere	x 0.0 0.0 10.0	y 0.0 0.0 10.0 couple xy update dipole dilate all
fix_modify	fNPH		temp ct press thermo_press
fix			fRIGID		np		rigid/nve	group 1 np

#output settings, changing peps will change the output file name as well, change this by removing ${peps} from the dump file name
''')

    f_in.write('dump			coords all custom 100 output.xyz id type x y z c_cls')
    f_in.write('''  
dump_modify	coords sort id

thermo_style	custom	step pe ke etotal

# set up our timestep and runtime

timestep       0.01
thermo         100
run            50000''') ### IF you want to change the length of the simulation, do it here
        
    f_in.close()
               
        


    print("Writing data file ... ")
    
    f_dat = open(InputPath + '/sd%d/data'%(seed),'w')

    # Reading the original data file
    f_ogdat = open(InputPath+"/Inputs//Membrane.data")
    lines = f_ogdat.readlines()
    f_ogdat.close()

    # Box lines
    blines = lines[7:10]
    # Membrane positions lines
    mlines = lines[41:2944]

    # Writing the data file
    f_dat.write('LAMMPS data file\n\n')
    f_dat.write(str(2901+NumPatches)+' atoms\n')
    f_dat.write(str(2+NumPatches)+' atom types\n')
    for line in blines:
        f_dat.write(line)
    f_dat.write('\n')
    f_dat.write('Masses\n\n')
    for i in range(NumPatches+2):
        f_dat.write(str(i+1)+' 1\n')
    f_dat.write('\n')
    for line in mlines:
        f_dat.write(line)
    for n in range(NumPatches):
        f_dat.write(str(2902+n)+ ' ' + str(types[n])+' ' +str(x[n])+' ' +str(y[n])+ ' ' +str(z[n]+6.5)+'  1 1 0   0 0 0\n') ## Inputing postitions of patches
    f_dat.close()
        
        
    
  


if __name__ == "__main__":
    main()
    
    
