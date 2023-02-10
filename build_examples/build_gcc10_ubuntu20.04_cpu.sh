#!/bin/bash
#################################################################
# "mpif90" is assumed to be in your PATH and points to 
# your chosen MPI library/compiler.
#################################################################
#################################################################
# Please set the location of the HDF5 include & library files. 
# Make sure the HDF5 LIBRARY is COMPILED with 
# the SAME COMPILER used here, and is in the run-time environment.
#################################################################

# NOTE: you must first install the hdf5 packages from the package manager.

HDF5_INCLUDE_DIR="/usr/include/hdf5/serial"
HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"

##################################################################
# Please set the HDF5 linker flags to match the installed version.
##################################################################

HDF5_LIB_FLAGS="-lhdf5_serial_fortran -lhdf5_serialhl_fortran -lhdf5_serial -lhdf5_serial_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
###########################################################################

FFLAGS="-O3 -march=native -fopenmp -ftree-parallelize-loops=${OMP_NUM_THREADS} -fallow-argument-mismatch"

###########################################################################
###########################################################################
###########################################################################

HIPFT_HOME=$PWD

pushd ${HIPFT_HOME}/src >> /dev/null
if [ -e Makefile ]; then
  \rm Makefile
fi 
sed \
  -e "s#<FFLAGS>#${FFLAGS}#g" \
  -e "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" \
  -e "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" \
  -e "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" \
  Makefile.template > Makefile
echo "make 1>build.log 2>build.err"
make clean
make 1>build.log 2>build.err

echo "cp ${HIPFT_HOME}/src/hipft ${HIPFT_HOME}/bin/hipft"
\cp ${HIPFT_HOME}/src/hipft ${HIPFT_HOME}/bin/hipft

