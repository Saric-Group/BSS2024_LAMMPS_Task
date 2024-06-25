#!/bin/bash

#this script requires wget, python 2.7, virtualenv and pip
#OS-X has trouble with the gpu package of lammps, install XCode8 (not 9) and CUDA

#commented out lines allow for the building of a gpu accelerated version. I will clean this process up at some point...

STARTDIR=$(pwd)
cd "$( dirname "${BASH_SOURCE[0]}" )"
WDIR=$(pwd)

LAMMPSDIR=./lammps

if [ ! -d $LAMMPSDIR ]; then

        if ! [ -x "$(command -v wget)" ]; then
                echo "you do not have wget installed on this computer, please download lammps Feb16 version manually and place it in a directory called 'lammps' at the same level as this file"
                echo "copy this link into your browser: http://lammps.sandia.gov/tars/lammps-16Feb16.tar.gz"
                echo "extract this tar and copy the resulting file to the lammps directory"
        else
                wget -qO- https://download.lammps.org/tars/lammps-16Feb2016.tar.gz | tar xvz 
                mv lammps* $LAMMPSDIR
        fi
fi

#copy src files to lammps folder
if [ ! -f $LAMMPSDIR/src/lmp_mpi ]; then
        cp -rf src/*.h $LAMMPSDIR/src
        cp -rf src/*.cpp $LAMMPSDIR/src
	cd $LAMMPSDIR/src
	make clean-all
	make no-all
	make yes-dipole yes-rigid yes-molecule
	make -j4 mpi
	cd "$WDIR"
fi

cd "${STARTDIR}"
echo "done!"
