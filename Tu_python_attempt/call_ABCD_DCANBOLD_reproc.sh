#!/bin/bash

# activate the python virtual environment (the packages are stored here)
source /data/wheelock/data1/people/Cindy/BCP/SystemMaturity/neuroHarmonize/bin/activate # if looping over subject this can be run once only

data_path='/data/Daenerys/ABCD/data/abcd_collection3165/derivatives/abcd-hcp-pipeline/'
results_path='/data/wheelock/data1/datasets/ABCD/derivatives/' 

subject_id='sub-NDARINV0A4P0LWM' # CHANGE THIS
task_id='task-rest01' # CHANGE THIS
json_path="${data_path}/${subject_id}/ses-baselineYear1Arm1/files/MNINonLinear/Results/${task_id}/DCANBOLDProc_v4.0.0/DCANBOLDProc_v4.0.0_mat_config.json"

echo $json_path

SECONDS=0
python3 DCANBOLD_reprocess.py --data-path $data_path --results-path $results_path --subject-id $subject_id --task-id $task_id --json-path $json_path --GSR 0 --savefig 1 # without GSR (24 params, movement-related regressors only)
echo Takes $SECONDS seconds to run. # You can comment this out if you don't want to estimate time

SECONDS=0
python3 DCANBOLD_reprocess.py --data-path $data_path --results-path $results_path --subject-id $subject_id --task-id $task_id --json-path $json_path --GSR 1 --savefig 0 # with GSR (30 params)
echo Takes $SECONDS seconds to run. # You can comment this out if you don't want to estimate time


# N.B.: json path example for EEG-fMRI:
# /data/wheelock/data1/datasets/AdultEEGfMRI/DockerOutput/sub-01/ses-01/files/MNINonLinear/Results/ses-01_task-rest_run-01/DCANBOLDProc_v4.0.0/DCANBOLDProc_v4.0.0_mat_config.json