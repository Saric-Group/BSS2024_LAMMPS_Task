# build a mixed-precision version of lammps for cpu

VERSION=stable_23Jun2022_update3
BUILD_NAME=serial_BSS24
wget https://github.com/lammps/lammps/archive/${VERSION}.tar.gz
#rm -rf lammps-_${VERSION}
tar zxvf ${VERSION}.tar.gz
cd lammps-${VERSION}

cp ../pair_membrane.* ./src/

mkdir build_${BUILD_NAME}
cd build_${BUILD_NAME}

module purge
module load intel/19.1.1.217
module load intel-mpi/intel/2019.7

cmake3 -D CMAKE_INSTALL_PREFIX=$HOME/.local \
-D LAMMPS_MACHINE=${BUILD_NAME} \
-D ENABLE_TESTING=no \
-D BUILD_MPI=no \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_COMPILER=icpc \
-D CMAKE_CXX_FLAGS_RELEASE="-Ofast -xHost -qopenmp -DNDEBUG" \
-D PKG_MOLECULE=yes \
-D PKG_RIGID=yes \
-D PKG_KSPACE=yes -D FFT=MKL -D FFT_SINGLE=yes \
-D PKG_EXTRA-PAIR=yes \
-D PKG_DIPOLE=yes \
../cmake

make -j 10
make install
