#!/bin/bash 
###########################################################################
#	
#	Template for environment setup for FNL_preproc
#	Save a copy as setup_env.sh, and update paths to those specific to
#	your system. Note there are currently also hard-coded config variables.
#
###########################################################################

# paths to dependencies #
path_to_label_files="*/ROI_sets/Surface_schemes/Human/"
path_to_movment_regressor_check="*/movmnt_regressor_check/movmnt_regressor_check.py"
FSL_DIR="*/fnl_lab/code/bin/fsl/fsl/"
octave="*/octave3.8/bin/octave"
wb_command="*/workbench/bin_linux64/wb_command" 
framewise_disp_path="*/framewise_displacement"
HCP_Mat_Path="*/HCP_Matlab"
# To mimic the old FNL pipelines
motion_filename='motion_numbers.txt'
skip_seconds=5
brain_radius_in_mm=50
expected_contiguous_frame_count=5
Matlab_Runtime_Env="/mnt/lustre1/fnl_lab/code/external/utilities/Matlab2016bRuntime/v91"
FNL_preproc_dir=`dirname ${0}`
FNL_preproc_version="FNL_preproc_v2"

# frame displacement th to calculate beta coefficients for regression
fd_th=0.3

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
