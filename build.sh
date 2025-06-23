#!/bin/bash
#################################################################
#
# Build script for HipFT.
#
# This script requires a configuration file.
#
# See the configure files in the "/conf" 
# folder for examples.
#
# We welcome contributions to the /conf folder for
# various systems and configurations.
#
#################################################################

codename="HipFT"
codeexe="hipft"

if [[ $# -lt 1 ]]; then
  echo "ERROR: Please specify a build configuration file:"
  echo "  Example:   ./build.sh conf/gcc_cpu_ubuntu.conf"
  exit 1
fi

if [[ $# -gt 1 ]]; then
  echo "ERROR: Please specify a build configuration file:"
  echo "  Example:   ./build.sh conf/gcc_cpu_ubuntu.conf"
  exit 1
fi

conf_file=$1

while read field_name field_value; do
  # Skips blank lines and comments.
  [[ -z "${field_name}" || "${field_name:0:1}" = "#" ]] && continue

  field_name="${field_name:0:${#field_name}-1}"
  if [[ "${field_name}" = "FC" ]]; then
    FC="${field_value}"
  elif [[ "${field_name}" = "FRTFLAGS" ]]; then
    FRTFLAGS="${field_value}"
  elif [[ "${field_name}" = "HDF5_INCLUDE_DIR" ]]; then
    HDF5_INCLUDE_DIR="${field_value}"
  elif [[ "${field_name}" = "HDF5_LIB_DIR" ]]; then
    HDF5_LIB_DIR="${field_value}"
  elif [[ "${field_name}" = "HDF5_LIB_FLAGS" ]]; then
    HDF5_LIB_FLAGS="${field_value}"
  elif [[ "${field_name}" = "USING_GCC" ]]; then
    USING_GCC="${field_value}"
#    if [[ "${field_value}" = "1" ]];then
#      if [[ "${OMP_NuM_THREADS" NOT SET ]];then
#	echo "WARNING! The multithreading of HipFT when using GCC needs OMP_NUM_THREADS set at compile time.  It is not currently set. It is recommeneded to set it, otherwise HipFT will be run in serial for each MPI proccess (realization block).  Setting it to 1 locally to allow build."
#	OMP_NUM_THREADS = "1"
#      fi
#    fi
  fi
done < "${conf_file}"

unset field_value
unset field_name

CODE_HOME=$PWD

cX="\033[0m"
cR="\033[1;31m"
cB="\033[1;34m"
cG="\033[32m"
cC="\033[1;96m"
cM="\033[35m"
cY="\033[1;93m"
Bl="\033[1;5;96m"
echo="echo -e"

${echo} "${cG}=== Starting $codename Build ===${cX}"
${echo} "==> Entering src directory..."
pushd ${CODE_HOME}/src > /dev/null
${echo} "==> Removing old Makefile..."
if [ -e Makefile ]; then
  \rm Makefile
fi 
${echo} "==> Generating Makefile from Makefile.template..."
sed \
  -e "s#<FC>#${FC}#g" \
  -e "s#<FRTFLAGS>#${FRTFLAGS}#g" \
  -e "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" \
  -e "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" \
  -e "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" \
  -e "s#<USING_GCC>#${USING_GCC}#g" \
  Makefile.template > Makefile
${echo} "==> Compiling code..."
make clean 1>/dev/null 2>/dev/null ; make 1>build.log 2>build.err
if [ ! -e $codeexe ]; then
  ${echo} "${cR}!!> ERROR!  $codeeexe executable not found.  Build most likely failed."
  ${echo} "            Contents of src/build.err:"
  cat build.err
  ${echo} "${cX}"
  exit 1
fi
$echo "==> Moving $codeexe executable to: ${CODE_HOME}/bin/$codeexe"
make install
${echo} "${cG}==> Build complete!${cX}"
${echo}      "    Please add the following to your shell startup (e.g. .bashrc, .profile, etc.):"
${echo} "${cC}    export PATH=${CODE_HOME}/bin:\$PATH${cX}"

