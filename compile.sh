#! /usr/bin/env bash
# compiles all

function usage() {
   echo "usage: `basename $0` <Matlab Compiler>"
}

echo STARTED:
date

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

fmr="rm -f ${DIR}/bin/run_filtered_movement_regressors.sh ; ${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o filtered_movement_regressors ${CODE}/filtered_movement_regressors.m -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS}"
dsp="rm -f ${DIR}/bin/run_dcan_signal_processsing.sh      ; ${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o dcan_signal_processsing      ${CODE}/dcan_signal_processing.m       -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS} -a ${FRAMEWISE_DISPLACEMENT}"
av2="rm -f ${DIR}/bin/run_analyses_v2.sh                  ; ${MCC_FILE} -v -m -R -singleCompThread -R -nodisplay -o analyses_v2                  ${CODE}/analyses_v2.m                  -a ${HCP_MATLAB_PATH} -a ${ADDITIONAL_SCRIPTS} -a ${FRAMEWISE_DISPLACEMENT}"

echo $fmr
eval $fmr
echo $dsp
eval $dsp
echo $av2
eval $av2

mv filtered_movement_regressors dcan_signal_processsing analyses_v2 run_*.sh ${BIN}/

#add MCR_CACHE_ROOT to all run scripts for Exahead processing
sed -i '/  shift 1/a \  RANDHASH=`cat /dev/urandom | tr -cd "a-f0-9" | head -c 8`\n  export MCR_CACHE_ROOT=$TMPDIR/$USER/$RANDHASH\n  mkdir -p $MCR_CACHE_ROOT' ${BIN}/run_*.sh

echo ENDED:
date
