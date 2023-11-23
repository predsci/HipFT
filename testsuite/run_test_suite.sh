#!/bin/bash

# HIPFT TESTSUITE
# Two main categories of tests:
#   - Those with exact analytic solutions.
#   - Those that wil be compared to previous runs.

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

hipftexe="${PWD}/../src/hipft"
echo $hipftexe

AVAIL_TEST_RUNS_LIST="
diffuse_soccer
"
#advect_tp_blob
#advect_t_blob
#advect_p_blob
#advect_p_diffuse_soccer
#advect_tp_diffuse_assim_atten_flowac_r3_map
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
    -novis)
    novis=1
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
    *)
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%%%% Warning!  Unknown option! %%%%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    ;;
esac
done

# ****** Get test suite parameters ******
WD=${PWD}
RESULTSDIR=${WD}/results
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

trap ctrl_c INT

function ctrl_c() {
  ${echo} "${cR}==> Caught CTRL-C, shutting down!${cX}"
  exit 1
}

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
${echo} "Welcome to the HIPFT test suite!!!"
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
    PTEST=$(which python)
    if [ -z "${PTEST}" ]
    then
      ${echo} "${cR}==> ERROR! Python does not seem to be installed!${cX}"
      ${echo} "${cR}    Please install it.${cX}"
      exit 1
    else
      PBIN=python
    fi
  else
    PBIN=python3
  fi
  ${echo} "${cG}==> Python is installed!${cX}"
 #
 # Check for required packages.
 #
  PYTHON_PKG_LIST="setuptools
  argparse
  numpy
  jdcal
  bottleneck
  h5py
  matplotlib"

  for pypkg in $PYTHON_PKG_LIST
  do
    $PBIN -c "import ${pypkg}" 2>/dev/null
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
  PTEST=$(which hipft_plot_histories_compare.py)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cY}==> ERROR HipFT BIN folder is not in the PATH!${cX}"
    ${echo} "${cY}==> Appending ${BINDIR} to PATH...${cX}"
    export PATH="${BINDIR}:${PATH}"
  fi
  PTEST=$(which hipft_plot_histories_compare.py)
  if [ -z "${PTEST}" ]; then
    ${echo} "${cR}==> ERROR! HipFT bin PATH problem!${cX}"
    exit 1
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

#
# ****** Remove previous comparison results if making new ones.
#
if [ ${nocompare} == 0 ]
then
  rm ${RESULTSDIR}/files/* 2>/dev/null
  rm ${RESULTSDIR}/files/* 2>/dev/null
  rm ${RESULTSDIR}/hipft_testsuite_results.txt 2>/dev/null
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

    if [ -e ${RUNDIR}/hipft.out ]
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
  if [ ! -e ${RUNDIR}/hipft.out ]
  then
    if [ ${norun} == 1 ]
    then
      ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} did not run correctly!${cX}"
      ${echo} "${cR}!!!> Check the run folder: ${RUNDIR} ${cX}"
      exit 1
    fi
  fi

#  run_completed_test=$(grep "RUN COMPLETE" ${RUNDIR}/eprem.log)

#  if [ -z "${run_completed_test}" ]
#  then
#   ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} did not seem to run correctly!${cX}"
#   exit 1
#  fi

  # Generate images etc.
#  if [ ${novis} == 0 ]
#  then
#    ${echo} "======================================================="
#    ${echo} "${cB}==> GENERATING IMAGES OF RUN${cX}"
#    ${echo} "======================================================="
#  fi

  if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 0 ]
  then
    ${echo} "${cR}=======================================================${cX}"
    ${echo} "${cR}==> SETTING REFERENCE DATA FOR RUN${cX}"
    ${echo} "${cR}=======================================================${cX}"

    ${echo} "${cR}==> Removing old reference data...${cX}"
    rm -fr ${REFDIR}/* 2>/dev/null
    ${echo} "${cR}==> Copying current run data into reference directory...${cX}"
    cp ${RUNDIR}/* ${REFDIR}/
  fi
#
# ****** Compare run data.
#
  if [ ${nocompare} == 0 ]
  then
    ${echo} "======================================================="
    ${echo} "${cB}==> COMPARING RUN DATA TO REFERENCE DATA${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running comparison script on all obs files..."
    result_file=${RESULTSDIR}/files/${TESTNAME}_hipft_compare_run.txt
    rm $result_file 2>/dev/null
    touch $result_file
    echo "=======================================================" >> $result_file
    echo "${TESTNAME} RESULTS" >> $result_file
    echo "=======================================================" >> $result_file
    hipft_compare_run.sh ${REFDIR} ${RUNDIR} >> $result_file
    #See if run passed:
    PASS_FAIL[${Ti}]=$(tail -n 2 $result_file | head -n 1)
    passfailcomp=( ${PASS_FAIL[${Ti}]} )

    if [ "${passfailcomp[0]}" = "FAIL" ]
    then
      ${echo} "${cR}==> Test seems to have FAILED!${cX}"
    else
      ${echo} "${cG}==> Test seems to have PASSED!${cX}"
    fi
    ${echo} "==> Adding comparison to summary file..."
    cat $result_file >> ${RESULTSDIR}/hipft_testsuite_results.txt
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
  ${echo} "${cY}==================================================${cX}"
  ${echo} "${cY}Summary of test results:${cX}"
  ${echo} "${cY}==================================================${cX}"
  echo "==================================================" > ${RESULTSDIR}/passfail.txt
  echo "Summary of test results:" >> ${RESULTSDIR}/passfail.txt
  echo "==================================================" >> ${RESULTSDIR}/passfail.txt

  Ti=0
  for TESTNAME in ${TEST_RUNS_LIST}
  do
    Ti=$((${Ti}+1))
    passfailcomp=( ${PASS_FAIL[${Ti}]} )
    if [ "${passfailcomp[0]}" = "FAIL" ]; then
      ${echo} "$(printf "%-25s %-25s" ${TESTNAME} ${cR}"${PASS_FAIL[${Ti}]}"${cX})"
      echo "$(printf "%-25s %-25s" ${TESTNAME} "${PASS_FAIL[${Ti}]}")">> ${RESULTSDIR}/passfail.txt
    else
      ${echo} "$(printf "%-25s %-25s" ${TESTNAME} ${cG}"${PASS_FAIL[${Ti}]}"${cX})"
      echo "$(printf "%-25s %-25s" ${TESTNAME} "${PASS_FAIL[${Ti}]}")" >> ${RESULTSDIR}/passfail.txt
    fi
  done
  echo "=======================================================" >> ${RESULTSDIR}/passfail.txt
  cat ${RESULTSDIR}/hipft_testsuite_results.txt >> ${RESULTSDIR}/passfail.txt
  mv ${RESULTSDIR}/passfail.txt ${RESULTSDIR}/hipft_testsuite_results.txt
fi

${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"

exit 0



