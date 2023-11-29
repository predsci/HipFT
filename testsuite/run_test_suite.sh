#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  ${echo} "${cR}==> Caught CTRL-C, shutting down!${cX}"
  exit 1
}

function display_help {
echo "       _    _ _       ______ _______ "
echo "      | |  | (_)     |  ____|__   __|"
echo "      | |__| |_ _ __ | |__     | |"
echo "      |  __  | | '_ \|  __|    | |"
echo "      | |  | | | |_) | |       | |"
echo "      |_|  |_|_| .__/|_|       |_|"
echo "               | |"
echo "               |_|"
echo ""
echo "      TEST SUITE v1.0.0"
echo "USAGE:   ./run_test_suite.sh"
echo ""
echo "By default, the above command will run the test suite on all"
echo "default tests using the 'hipft' executable from the '../bin/'"
echo "folder (assuming it has been built)."
echo ""
echo "OPTIONS:"
echo ""
echo "Non-flag options '-opt=' are specified as '-opt=<OPTION>'."
echo ""
echo "-nochecksetup   Don't check the environment."
echo ""
echo "-hipftexe=      Use this to run the testsuite on a specific hipft executable."
echo "                This should be a full path and is useful for development tests."
echo ""
echo "-test=          Comma-seperated list of subset of tests to run."
echo "                Also can be used to run non-standard/experimental tests."
echo ""
echo "-nocleanup      By default, only the initial and final conditions of a run are "
echo "                kept in the run folders. Set this to keep the full run."
echo ""
echo "-clear          Remove all previous test suite runs and results - CAREFUL!"
echo ""
echo "-norun          Does not run the MAS code.  Checks for a previous run and compares if found."
echo ""
echo "-nocompare      Do not compare the runs to their reference runs."
echo ""
echo "-nocolor        Disable color text output."
echo ""
}

#Set number of processors to use for testsuite:
#(All tests use the same number of procs for now).
NUM_PROCS=1
norun=0
nocompare=0
novis=0
nocleanup=0
nochecksetup=0
nocolor=0
setrefdata=0
algnum=3
hipftexe="hipft"

AVAIL_TEST_RUNS_LIST="
diffuse_soccer
diffuse_advect_soccer
diffuse_dipole
"
#advect_t_blob
#advect_p_blob
#diffuse_advect_flows_map
#"

TEST_RUNS_LIST=${AVAIL_TEST_RUNS_LIST}

for i in "$@"
do
case $i in
    -norun)
    norun=1
    ;;
    -nocompare)
    nocompare=1
    ;;
    -nocleanup)
    nocleanup=1
    ;;
    -nochecksetup)
    nochecksetup=1
    ;;
    -nocolor)
    nocolor=1
    ;;
    -setrefdata)
    setrefdata=1
    ;;
    -test=*)
    TEST_RUNS_LIST="${i#*=}"
    TEST_RUNS_LIST="${TEST_RUNS_LIST//','/' '}"
    ;;
    -hipftexe=*)
    hipftexe="${i#*=}"
    ;;
    -h)
    display_help
    exit 0
    ;;
    --help)
    display_help
    exit 0
    ;;    
    *)
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "ERROR!  Unknown option: $i"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    display_help
    exit 1    
    ;;
esac
done

# ****** Get test suite parameters ******
WD=${PWD}
BINDIR=${WD}/../bin
SRCDIR=${WD}/../src
ROOTDIR=${WD}/..
TSLOG=${RESULTSDIR}/testsuite.log

if [ ${nocolor} == 0 ]
then
  cX="\033[0m"
  cR="\033[1;31m"
  cB="\033[1;34m"
  cG="\033[32m"
  cC="\033[1;96m"
  cM="\033[35m"
  cY="\033[1;93m"
  Bl="\033[1;5;96m"
  echo="echo -e"
else
  cX=
  cR=
  cB=
  cG=
  cC=
  cM=
  cY=
  Bl=
  echo="echo"
fi

#rm ${TSLOG} 1>/dev/null 2>/dev/null

#touch ${TSLOG}

${echo} " "
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "${cY}      _    _ _       ______ _______                                   ${cX}"
${echo} "${cY}     | |  | (_)     |  ____|__   __|                                  ${cX}"
${echo} "${cY}     | |__| |_ _ __ | |__     | |                                     ${cX}"
${echo} "${cY}     |  __  | | '_ \|  __|    | |                                     ${cX}"
${echo} "${cY}     | |  | | | |_) | |       | |                                     ${cX}"
${echo} "${cY}     |_|  |_|_| .__/|_|       |_|                                     ${cX}"
${echo} "${cY}              | |                                                     ${cX}"
${echo} "${cY}              |_|                                                     ${cX}"
${echo} "${cY}  ___ ____ ____ ___    ____ _  _ _ ___ ____                           ${cX}"
${echo} "${cY}   |  |___ [__   |     [__  |  | |  |  |___                           ${cX}"
${echo} "${cY}   |  |___ ___]  |     ___] |__| |  |  |___                           ${cX}"
${echo} "${vY}                                                                      ${cX}"
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "Welcome to the HIPFT test suite!"
#
# ****** Test for correct prerequisites and environment ******
#
${echo} "Checking file structure..."

if [ ! -e ${epremexe} ]
  then
  ${echo} "${cR}!!!> ERROR! HipFT binary executable missing!${cX}"
  exit 1
fi

#
# ****** Test for correct prerequisites and environment ******
#
if [ ${nochecksetup} == 0 ]
then
  ${echo} "Checking software requirements..."
  #Check that python is installed:
  PTEST=$(which python3)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cR}==> ERROR! Python3 does not seem to be installed!${cX}"
    ${echo} "${cR}    Please install it.${cX}"
    exit 1
  fi
  ${echo} "${cG}==> Python3 is installed!${cX}"
 #
 # Check for required packages.
 #
  PYTHON_PKG_LIST="argparse
  sys
  numpy
  "

  for pypkg in $PYTHON_PKG_LIST
  do
    python3 -c "import ${pypkg}" 2>/dev/null
    pychk=$?
    if [ $pychk -eq 1 ]; then
      ${echo} "${cR}==> ERROR! Missing required package ${pypkg}.  Please install it and try again.${cX}"
      exit 1
    fi
    ${echo} "==>        package ${cG}${pypkg}${cX} found!"
  done
#
# Check that the HipFT bin directory is in the user's path, if not, add it.
#
  ${echo} "==> Checking PATH...."
  PTEST=$(which hipft)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cY}==> WARNING: HipFT is not in the PATH!${cX}"
    ${echo} "${cY}==> Appending ${BINDIR} to PATH...${cX}"
    export PATH="${BINDIR}:${PATH}"
  fi
  PTEST=$(which hipft)
  if [ -z "${PTEST}" ]; then
    ${echo} "${cR}==> ERROR! HipFT bin PATH problem!${cX}"
    exit 1
  fi
  ${echo} "${cG}==> Using HIPFT binary: $PTEST ${cX}"

  if [ -z "${OMP_NUM_THREADS}" ]; then
    ${echo} "${cY}==> WARNING: OMP_NUM_THREADS is not set.${cX}"
    ${echo} "${cY}==> If running on CPUs, run may not be parallelized.${cX}"
    ${echo} "${cY}==> For GCC, OMP_NUM_THREADS needs to be set at compile time!${cX}"
  fi
    
  ${echo} "${cG}==> Everything seems OK to run HipFT test suite!${cX}"
fi
#
# ****** Check that user is running the tests and really wants to set ref data:
#
if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 1 ]
then
  ${echo} "${cR} ==> ERROR!  You are trying to set reference data without running the tests!${cX}"
  exit 1
fi
if [ ${setrefdata} -eq 1 ]
then
  ${echo} "${cR}╔═╗┌─┐┌┬┐  ╔╗╔┌─┐┬ ┬  ╦═╗┌─┐┌─┐┌─┐┬─┐┌─┐┌┐┌┌─┐┌─┐  ╔╦╗┌─┐┌┬┐┌─┐${cX}"
  ${echo} "${cR}╚═╗├┤  │   ║║║├┤ │││  ╠╦╝├┤ ├┤ ├┤ ├┬┘├┤ ││││  ├┤    ║║├─┤ │ ├─┤${cX}"
  ${echo} "${cR}╚═╝└─┘ ┴   ╝╚╝└─┘└┴┘  ╩╚═└─┘└  └─┘┴└─└─┘┘└┘└─┘└─┘  ═╩╝┴ ┴ ┴ ┴ ┴${cX}"
  read -p "==> Setting new reference data after runs complete...are you SURE? (y/n):" yn
  if [ ${yn} = "y" ]
  then
    read -p "==> Are you really really really SURE?! (y/n):" yn
    if [ ${yn} = "n" ]
    then
      ${echo} "==> ${cR}Aborting!${cX}"
      exit 0
    fi
  else
    ${echo} "==> ${cR}Aborting!${cX}"
    exit 0
  fi
fi

########################################################################
########################################################################
##
## ****** Start loop through test problems ******
##
########################################################################
########################################################################

Ti=0
for TESTNAME in ${TEST_RUNS_LIST}
do
  Ti=$((${Ti}+1))
#
# ****** Make sure test is in the test suite:
#
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  testok=0
  for run_test in ${AVAIL_TEST_RUNS_LIST}; do
    if [[ ${run_test} == ${TESTNAME} ]]; then
      testok=1
      break
    fi
  done
  if [ ${testok} -eq 0 ]
  then
    ${echo} "${cR}==> TEST ${cX}${cM}${TESTNAME}${cX}${cR} is not a valid test in the testsuite!${cX}"
    continue
  fi
  ${echo} "STARTING TEST ${cM}${TESTNAME}${cX}"
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  ${echo} "==> Gathering test information..."
# ****** Get directories:
  RUNDIR=${WD}/${TESTNAME}/run
  REFDIR=${WD}/${TESTNAME}/reference
  INPUTDIR=${WD}/${TESTNAME}/input

  if [ ! -d ${RUNDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Run directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${REFDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Reference directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${INPUTDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Input directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -e ${INPUTDIR}/hipft.in ]
  then
    ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} does not have an input file!${cX}"
    exit 1
  fi

  cd ${RUNDIR}

  if [ ${norun} == 0 ]
  then

    if [ -e ${RUNDIR}/hipft_run_parameters_used.out ]
    then
      ${echo} "==> Clearing run directory..."
      rm -fr ${RUNDIR}/*
    fi

    ${echo} "==> Copying input files..."
    cp ${INPUTDIR}/* ${RUNDIR}/ 2>/dev/null

    cd ${RUNDIR}

    ${echo} "======================================================="
    ${echo} "${cB}==> RUNNING HIPFT${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running hipft with command:"
    ${echo} "==> mpiexec -np $NUM_PROCS ${hipftexe} hipft.in 1>hipft.log 2>hipft.err"
    mpiexec -np $NUM_PROCS ${hipftexe} hipft.in 1>hipft.log 2>hipft.err
  fi

  # Check that a completed run exists in the run folder
  if [ ! -e ${RUNDIR}/hipft_brmap_final.h5 ]
  then
    if [ ${norun} == 1 ]
    then
      ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} did not run correctly!${cX}"
      ${echo} "${cR}!!!> Check the run folder: ${RUNDIR} ${cX}"
      exit 1
    fi
  fi
  
  if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 0 ]
  then
    ${echo} "${cR}=======================================================${cX}"
    ${echo} "${cR}==> SETTING REFERENCE DATA FOR RUN${cX}"
    ${echo} "${cR}=======================================================${cX}"

    ${echo} "${cR}==> Removing old reference data...${cX}"
    rm -fr ${REFDIR}/* 2>/dev/null
    ${echo} "${cR}==> Copying current run data into reference directory...${cX}"
    cp ${RUNDIR}/* ${REFDIR}/
    rm -f ${REFDIR}/*.h5
  fi

  # Get timing data:
  TIME_RUN_TMP=($(grep "Wall clock time" ${RUNDIR}/hipft.log))
  TIME_RUN_TMP=${TIME_RUN_TMP[3]}

  TIME_REF_TMP=($(grep "Wall clock time" ${REFDIR}/hipft.log))
  TIME_REF_TMP=${TIME_REF_TMP[3]}

  SPEEDUP_TMP=`python3 -c "print(${TIME_REF_TMP}/${TIME_RUN_TMP})"`

  TIME_RUN[${Ti}]=$(printf "%8.3f" ${TIME_RUN_TMP})
  TIME_REF[${Ti}]=$(printf "%8.3f" ${TIME_REF_TMP})
  SPEEDUP[${Ti}]=$(printf "%5.2f" ${SPEEDUP_TMP})

  ${echo} "${cG}==> Test completed in:               ${TIME_RUN[${Ti}]} seconds.${cX}"
  ${echo} "${cG}==> Speedup compared to reference run: ${SPEEDUP[${Ti}]} x${cX}"

# Generate images etc.
#  if [ ${novis} == 0 ]
#  then
#    ${echo} "======================================================="
#    ${echo} "${cB}==> GENERATING IMAGES OF RUN${cX}"
#    ${echo} "======================================================="
#  fi


#
# ****** Compare run data.
#
  if [ ${nocompare} == 0 ]
  then
    ${echo} "======================================================="
    ${echo} "${cB}==> COMPARING RUN DATA TO REFERENCE DATA${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running comparison..."

    PASS_FAIL[${Ti}]=$(hipft_compare_run_diags.py ${REFDIR} ${RUNDIR})

    if [ "${PASS_FAIL[${Ti}]}" = "0" ]
    then
      ${echo} "${cG}==> Test seems to have PASSED!${cX}"
    else
      ${echo} "${cR}==> Test seems to have FAILED!${cX}"
      ${echo} "${cR}==> ${PASS_FAIL[${Ti}]} ${cX}"
      PASS_FAIL[${Ti}]=1
    fi
  fi
#
# ****** Cleanup run data.
#
  if [ ${nocleanup} == 0 ]
  then
    if [ ${norun} == 0 ]
    then
      ${echo} "======================================================="
      ${echo} "${cB}==> CLEANING RUN DATA${cX}"
      ${echo} "======================================================="
      ${echo} "==> Removing files from run data..."
      rm -fr ${RUNDIR}/*
    fi
  fi
done
${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"


# Display summary and timing results.
if [ ${nocompare} == 0 ]
then
  ${echo} "${cY}=================================================================${cX}"
  ${echo} "${cY}Summary of test results:${cX}"
  ${echo} "${cY}=================================================================${cX}"

  Ti=0
    ${echo} "$(printf "%-25s  %9s  %8s  %8s  %7s" "Test name" "PASS/FAIL" "Run-time" "Ref-time" "Speedup")"
  ${echo} "${cY}=================================================================${cX}"    
  for TESTNAME in ${TEST_RUNS_LIST}
  do
    Ti=$((${Ti}+1))
    passfailcomp=( ${PASS_FAIL[${Ti}]} )
    if [ "${passfailcomp[0]}" = "1" ]; then
      pf="${cR}FAIL     ${cX}"
    else
      pf="${cG}PASS     ${cX}"
    fi
    ${echo} "$(printf "%-25s  %9s  %8s  %8s  %7s" "${TESTNAME}" "${pf}" "${TIME_RUN[${Ti}]}" "${TIME_REF[${Ti}]}" "${SPEEDUP[${Ti}]}")"
  done
fi

${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"

exit 0



