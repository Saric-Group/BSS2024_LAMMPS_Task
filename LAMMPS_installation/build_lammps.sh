#!/bin/bash

STARTDIR=$(pwd)
cd "$( dirname "${BASH_SOURCE[0]}" )"
WDIR=$(pwd)

LAMMPSDIR=./lammps

if [ ! -d $LAMMPSDIR ]; then

        if ! [ -x "$(command -v curl)" ]; then
                echo "you do not have curl installed on this computer, please download lammps Feb16 version manually and place it in a directory called 'lammps' at the same level as this file"
                echo "copy this link into your browser: https://download.lammps.org/tars/lammps-16Feb2016.tar.gz"
                echo "extract this tar and copy the resulting file to the lammps directory"
        else
                curl -O https://download.lammps.org/tars/lammps-16Feb2016.tar.gz
		tar -xzvf lammps-16Feb2016.tar.gz -C $LAMMPSDIR
                # mv lammps* $LAMMPSDIR
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
	make serial
	cd "$WDIR"
fi

cd "${STARTDIR}"
echo "done!"
