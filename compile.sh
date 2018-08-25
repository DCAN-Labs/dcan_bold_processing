#! /usr/bin/env bash
set -e
# compiles all project matlab code using the matlab compiler tool

function usage() {
   echo "usage: `basename $0` <Matlab Compiler>"
}

if [ $# -eq 1 ]; then
  MCC_FILE="$1"
fi

if [ ! -e $MCC_FILE ]; then
  echo "unable to locate matlab compiler"
  usage
  exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# set paths
CODE="${DIR}/matlab_code"
BIN="${DIR}/bin"
HCP_MATLAB_PATH="${CODE}/hcp_matlab"
FRAMEWISE_DISPLACEMENT_PATH="${CODE}/framewise_displacement"
SCRIPTS_PATH="${CODE}/scripts"

rm -f ${DIR}/bin/*

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o filtered_movement_regressors    "${CODE}"/filtered_movement_regressors.m    -a ${HCP_MATLAB_PATH} -a ${SCRIPTS_PATH}

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o dcan_signal_processsing         "${CODE}"/DCAN_Signal_Processing.m          -a ${HCP_MATLAB_PATH} -a ${SCRIPTS_PATH} -a ${FRAMEWISE_DISPLACEMENT_PATH}

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o analyses_v2                     "${CODE}"/analyses_v2.m                     -a ${HCP_MATLAB_PATH} -a ${SCRIPTS_PATH} -a ${FRAMEWISE_DISPLACEMENT_PATH}

mv filtered_movement_regressors dcan_signal_processsing analyses_v2 run_filtered_movement_regressors.sh run_dcan_signal_processsing.sh run_analyses_v2.sh "$BIN"/

#add MCR_CACHE_ROOT to all run scripts for Exahead processing
sed -i '/exe_dir=`dirname "$0"`/a if [ ! -d $TMPDIR/$USER ]; then\n    mkdir $TMPDIR/$USER\nfi\nexport MCR_CACHE_ROOT=$TMPDIR/$USER' ${BIN}/run_*.sh

