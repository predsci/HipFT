#!/bin/bash
#################################################################
# Enter your MPI compiler (typically "mpif90").
#################################################################

FC="mpif90"

#################################################################
# Please set the location of the HDF5 include & library files. 
# Make sure the HDF5 LIBRARY is COMPILED with 
# the SAME COMPILER used here, and is in the run-time environment.
#################################################################

HDF5_INCLUDE_DIR="${PS_EXT_DEPS_HOME}/hdf5/include"
HDF5_LIB_DIR="${PS_EXT_DEPS_HOME}/hdf5/lib"

##################################################################
# Please set the HDF5 linker flags to match the installed version.
##################################################################

HDF5_LIB_FLAGS="-lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
###########################################################################

FFLAGS="-O3 -march=native -acc=gpu -mp=gpu -stdpar=gpu -gpu=cc60,cc70,cc75,cc80,cc86,cc89,cc90,nomanaged -Minfo=accel"

###########################################################################
# Specify src filename:  Use hipft_gcc.f90 for GCC, otherwise use hipft.f90
###########################################################################

SRCFILE="hipft.f90"

###########################################################################
###########################################################################
###########################################################################

HIPFT_HOME=$PWD

echo "Entering src directory..."
pushd ${HIPFT_HOME}/src >> /dev/null
echo "Removing old Makefile..."
if [ -e Makefile ]; then
  \rm Makefile
fi 
echo "Generating Makefile from Makefile.template..."
sed \
  -e "s#<FC>#${FC}#g" \
  -e "s#<FFLAGS>#${FFLAGS}#g" \
  -e "s#<SRCFILE>#${SRCFILE}#g" \
  -e "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" \
  -e "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" \
  -e "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" \
  Makefile.template > Makefile
echo "Compiling code..."
make clean 1>/dev/null 2>/dev/null ; make 1>build.log 2>build.err
echo "Copying hipft executable to: ${HIPFT_HOME}/bin/hipft"
\cp ${HIPFT_HOME}/src/hipft ${HIPFT_HOME}/bin/hipft
echo "Build complete!  Please add the following to your shell:"
echo "export PATH=${HIPFT_HOME}/bin:\$PATH"

