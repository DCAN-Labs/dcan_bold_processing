#! /usr/bin/env bash
# compiles all 

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
FRAMEWISE_DISPLACEMENT="${CODE}/framewise_displacement"
ADDITIONAL_SCRIPTS="${CODE}/scripts"

rm -f ${DIR}/bin/*

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o "${BIN}"/filtered_movement_regressors    "${CODE}"/filtered_movement_regressors.m    -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS}

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o "${BIN}"/dcan_signal_processsing         "${CODE}"/dcan_signal_processing.m          -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS} -a ${FRAMEWISE_DISPLACEMENT_PATH}

${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o "${BIN}"/analyses_v2                     "${CODE}"/analyses_v2.m                     -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS} -a ${FRAMEWISE_DISPLACEMENT_PATH}

#add MCR_CACHE_ROOT to all run scripts for Exahead processing
sed -i '/exe_dir=`dirname "$0"`/a if [ ! -d $TMPDIR/$USER ]; then\n    mkdir $TMPDIR/$USER\nfi\nexport MCR_CACHE_ROOT=$TMPDIR/$USER' "${BIN}"/run_*.sh

