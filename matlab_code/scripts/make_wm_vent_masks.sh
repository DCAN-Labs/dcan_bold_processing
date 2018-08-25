#!/bin/bash 
set -e

# Requirements for this script
#  installed versions of: FSL  
#  environment: FSLDIR

################################################ SUPPORT FUNCTIONS ##################################################

Usage() {
  echo "`basename $0`: Script to make masks for ventricles and white matter"
  echo " "
  echo "Usage: `basename $0` --workingdir=<working dir>"
  echo "	     --wm_mask_L=<name for the wm Left hemisphere mask>"
  echo "	     --wm_mask_R=<name for the wm Right hemisphere mask>"
  echo "	     --wm_mask=<name for the wm mask>"
  echo "	     --vent_mask_L=<name for the vent Left hemisphere mask>"
  echo "	     --vent_mask_R=<name for the vent Reft hemisphere mask>"
  echo "	     --vent_mask=<name for the vent mask>"
  echo "	     --wm_lt_L=<white matter lower threshold Left>"
  echo "	     --wm_ut_L=<white matter upper threshold Left>"
  echo "	     --wm_lt_R==<white matter lower threshold Right>"
  echo "	     --wm_ut_R=<white matter upper threshold Right>"
  echo "	     --vent_lt_R=<ventricles lower threshold Right>"
  echo "	     --vent_ut_R=<ventricles upper threshold Right>"
  echo "	     --vent_lt_L=<ventricles lower threshold Left>"
  echo "	     --vent_ut_L=<ventricles upper threshold Left>"
  echo "	     --segBrain=<name of the segmented brain's file>"
  echo "	     --segBrainDir=<path to the segmented's brain file>"
  echo "	     --wm_mask_eroded=<name for the eroded wm mask>"
  echo "	     --vent_mask_eroded=<name for the eroded vent mask>"
  echo " "
}

# function for parsing options
getopt1() {
    sopt="$1"
    shift 1
    for fn in $@ ; do
	if [ `echo $fn | grep -- "^${sopt}=" | wc -w` -gt 0 ] ; then
	    echo $fn | sed "s/^${sopt}=//"
	    return 0
	fi
    done
}

defaultopt() {
    echo $1
}

################################################### OUTPUT FILES #####################################################

# Outputs (in $WD): 
#         NB: all these images are in standard space 
#             but at the specified resolution (to match the fMRI - i.e. low-res)
#     ${T1wImageFile}.${FinalfMRIResolution}  
#     ${FreeSurferBrainMaskFile}.${FinalfMRIResolution}
#     ${BiasFieldFile}.${FinalfMRIResolution}  
#     Scout_gdc_MNI_warp     : a warpfield from original (distorted) scout to low-res MNI
#
# Outputs (not in either of the above):
#     ${OutputTransform}  : the warpfield from fMRI to standard (low-res)
#     ${OutputfMRI}       
#     ${JacobianOut}
#     ${ScoutOutput}
#          NB: last three images are all in low-res standard space

################################################## OPTION PARSING #####################################################

# Just give usage if no arguments specified
if [ $# -eq 0 ] ; then Usage; exit 0; fi
# check for correct options

# parse arguments
wm_mask_L=`getopt1 "--wm_mask_L" $@`  # "$1"
wm_mask_R=`getopt1 "--wm_mask_R" $@`  # "$2"
wm_mask=`getopt1 "--wm_mask" $@`  # "$3"
vent_mask_L=`getopt1 "--vent_mask_L" $@`  # "$4"
vent_mask_R=`getopt1 "--vent_mask_R" $@`  # "$5"
vent_mask=`getopt1 "--vent_mask" $@`  # "$6"
wm_lt_L=`getopt1 "--wm_lt_L" $@`  # "$7"
wm_ut_L=`getopt1 "--wm_ut_L" $@`  # "$8"
wm_lt_R=`getopt1 "--wm_lt_R" $@`  # "$9"
wm_ut_R=`getopt1 "--wm_ut_R=" $@`  # "${10}"
vent_lt_R=`getopt1 "--vent_lt_R" $@`  # "${11}"
vent_ut_R=`getopt1 "--vent_ut_R" $@`  # "${12}"
vent_lt_L=`getopt1 "--vent_lt_L" $@`  # "${13}"
vent_ut_L=`getopt1 "--vent_ut_L" $@`  # "${14}"
segBrain=`getopt1 "--segBrain" $@`  # "${15}"
segBrainDir=`getopt1 "--segBrainDir" $@`  # "${16}"
wm_mask_eroded=`getopt1 "--wm_mask_eroded" $@`  # "${17}"
vent_mask_eroded=`getopt1 "--vent_mask_eroded" $@`  # "${18}"

echo " "
echo " START: making masks"

########################################## DO WORK ########################################## 
echo " making wm mask"
fslmaths ${segBrainDir}/${segBrain} -thr ${wm_lt_R} -uthr ${wm_ut_R} ${segBrainDir}/${wm_mask_R}
fslmaths ${segBrainDir}/${segBrain} -thr ${wm_lt_L} -uthr ${wm_ut_L} ${segBrainDir}/${wm_mask_L}
fslmaths ${segBrainDir}/${wm_mask_R} -add ${segBrainDir}/${wm_mask_L} -bin ${segBrainDir}/${wm_mask}
fslmaths ${segBrainDir}/${wm_mask} -kernel gauss 2 -ero ${segBrainDir}/${wm_mask_eroded}


echo " making vent mask"
fslmaths ${segBrainDir}/${segBrain} -thr ${vent_lt_R} -uthr ${vent_ut_R} ${segBrainDir}/${vent_mask_R}
fslmaths ${segBrainDir}/${segBrain} -thr ${vent_lt_L} -uthr ${vent_ut_L} ${segBrainDir}/${vent_mask_L}
fslmaths ${segBrainDir}/${vent_mask_R} -add ${segBrainDir}/${vent_mask_L} -bin ${segBrainDir}/${vent_mask}
fslmaths ${segBrainDir}/${vent_mask} -kernel gauss 2 -ero ${segBrainDir}/${vent_mask_eroded}

echo "clean up"
rm ${segBrainDir}/${wm_mask_R}
rm ${segBrainDir}/${wm_mask_L}
rm ${segBrainDir}/${wm_mask}
rm ${segBrainDir}/${vent_mask_R}
rm ${segBrainDir}/${vent_mask_L}
rm ${segBrainDir}/${vent_mask}

echo " "
echo "END: masks are done"
#echo " END: `date`" >> $WD/log.txt

########################################## QA STUFF ########################################## 


