#!/bin/bash 

# paths to dependencies #
path_to_label_files="/home/exacloud/lustre1/fnl_lab/ROI_sets/Surface_schemes/Human/"
path_to_movment_regressor_check="/home/exacloud/lustre1/fnl_lab/code/movmnt_regressor_check/movmnt_regressor_check.py"
FSL_DIR="/home/exacloud/lustre1/fnl_lab/code/bin/fsl/fsl/"
FSLDIR="/home/exacloud/lustre1/fnl_lab/code/bin/fsl/fsl/"
octave="/home/exacloud/lustre1/fnl_lab/code/bin/octave/bin/octave"
wb_command="/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench/bin_rh_linux64/wb_command" 
PATH=$FSLDIR/bin:$PATH
export PATH FSLDIR
. ${FSLDIR}/etc/fslconf/fsl.sh

framewise_disp_path="/home/exacloud/lustre1/fnl_lab/code/framewise_displacement"
HCP_Mat_Path="/home/exacloud/lustre1/fnl_lab/code/HCP_Matlab"
Matlab_Runtime_Env="/mnt/lustre1/fnl_lab/code/external/utilities/Matlab2016bRuntime/v91"
FNL_preproc_dir=`dirname ${0}`
FNL_preproc_version="FNL_preproc_v2"

# To mimic the old FNL pipelines
motion_filename='motion_numbers.txt'
skip_seconds=7
brain_radius_in_mm=50
expected_contiguous_frame_count=5

# frame displacement th to calculate beta coefficients for regression
fd_th=0.2

# Define filter parameters
bp_order=2 #band pass filter order
lp_Hz=0.009 # low pass frequency, Hz
hp_Hz=0.080 # high pass frequency, Hz

# Define constants
vent_lt_L=4 # white matter lower threshold Left
vent_ut_L=4 # white matter upper threshold Left
vent_lt_R=43 # white matter lower threshold Right
vent_ut_R=43 # white matter upper threshold Right

wm_lt_R=2950 # ventricles lower threshold Right
wm_ut_R=3050 # ventricles upper threshold Right
wm_lt_L=3950 # ventricles lower threshold Left
wm_ut_L=4050 # ventricles upper threshold Left
