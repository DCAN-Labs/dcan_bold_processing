#!/bin/bash

# activate the python virtual environment (the packages are stored here)
source /data/wheelock/data1/people/Cindy/BCP/SystemMaturity/neuroHarmonize/bin/activate # if looping over subject this can be run once only

SECONDS=0
data_path='/data/Daenerys/ABCD/data/abcd_collection3165/derivatives/abcd-hcp-pipeline/'
results_path='/data/wheelock/data1/datasets/ABCD/derivatives/' 

# To-do make the script accept command line arguments so that it can loop subjects
subject_id='sub-NDARINV0A4P0LWM' # CHANGE THIS
task_id='task-rest01' # CHANGE THIS
 
python3 DCANBOLD_reprocess.py --data-path $data_path --results-path $results_path --subject-id $subject_id --task-id $task_id --GSR 0 # without GSR (24 params, movement-related regressors only)
python3 DCANBOLD_reprocess.py --data-path $data_path --results-path $results_path --subject-id $subject_id --task-id $task_id --GSR 1 # with GSR (30 params)

echo Takes $SECONDS seconds to run. # You can comment this out if you don't want to estimate time