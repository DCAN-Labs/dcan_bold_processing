# dcan signal processing

\*\*
This is a repository for the dcan labs bold signal processing. It is 
forked from FNL\_preproc and is meant to take its place.
\*\*

This code repository consists of python wrappers and matlab scripts
for signal processing of the bold signal extracted from fMRI data.
This program is designed for the explicit output data of the HCP
fMRI pipeline or its DCAN derivatives. It is not designed with other 
preprocessed data in mind, so use at your own peril.


## installation

Installation requires use of the matlab compiler tool distributed with 
matlab, or acquiring an up-to-date version of the binaries available upon 
request.

git clone git@gitlab.com:Fair\_lab/dcan\_signal\_processing.git

cd dcan\_signal\_processing

./compile.sh \<matlab compiler path\>


## dcan\_bold\_proc.py

main wrapper for signal processing scripts.

```{bash}

usage: dcan_bold_proc.py [-h] [-v] [--setup] --subject SUBJECT --task TASK
                         [--output-folder OUTPUT_FOLDER]
                         [--legacy-tasknames]
                         [--filter-order FILTER_ORDER] [--lower-bpf LOWER_BPF]
                         [--upper-bpf UPPER_BPF] [--fd-threshold FD_THRESHOLD]
                         [--skip-seconds SKIP_SECONDS]
                         [--contiguous-frames CONTIGUOUS_FRAMES]
                         [--brain-radius BRAIN_RADIUS]
                         [--motion-filter-type {notch,lp}]
                         [--motion-filter-order MOTION_FILTER_ORDER]
                         [--band-stop-min BAND_STOP_MIN]
                         [--band-stop-max BAND_STOP_MAX]
                         [--motion-filter-option MOTION_FILTER_OPTION]
                         [--teardown] [--tasklist TASKLIST] [--physio PHYSIO]

Wraps the compiled DCAN Signal Processing Matlab script, version: 4.0.0.
Runs in 3 main modes:  [setup], [task], and [teardown].

[setup]: creates white matter and ventricular masks for regression, must be
         run prior to task.

[task]: computes fd numbers [1][2], runs regressions on a given task/fmri [3]
        and outputs a corrected dtseries, along with motion numbers in an
        hdf5 (.mat) formatted file.

[teardown]: concatenates any resting state runs into a single dtseries, and
            parcellates all final tasks.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print the software name and version
  --setup               prepare white matter and ventricle masks, must be run
                        prior to individual task runs.
  --subject SUBJECT     subject/participant id
  --task TASK           name of fmri data as used in the dcan fmri pipeline.
                        For bids data it is set to "task-NAME"
  --output-folder OUTPUT_FOLDER
                        output folder which contains all files produced by the
                        dcan fmri-pipeline. Used for setting up standard
                        inputs and outputs
  --legacy-tasknames
                        parse input task names as done in dcan_bold_processing <= 4.0.4.
                        use this flag if the input task filenames use the older DCAN HCP 
                        pipeline filename convention in which run index is appended to 
                        task name, e.g. task-myTask01 instead of task-myTask_run-01. 

bold signal filtering:
  bold signal filtering parameters.

  --filter-order FILTER_ORDER
                        number of filter coefficients for butterworth bandpass
                        filter.
  --lower-bpf LOWER_BPF
                        lower cut-off frequency (Hz) for the butterworth
                        bandpass filter.
  --upper-bpf UPPER_BPF
                        upper cut-off frequency (Hz) for the butterworth
                        bandpass filter.

framewise displacement:
  parameters related to computation of framewise displacment (FD)

  --fd-threshold FD_THRESHOLD
                        upper framewise displacement threshold for use in
                        signal regression.
  --skip-seconds SKIP_SECONDS
                        number of seconds to cut off the beginning of fmri
                        time series.
  --contiguous-frames CONTIGUOUS_FRAMES
                        number of contigious frames for power 2014 fd
                        thresholding.
  --brain-radius BRAIN_RADIUS
                        radius of brain for computation of rotational
                        displacements
  --motion-filter-type {notch,lp}
                        type of band-stop filter to use for removing
                        respiratory artifact from motion regressors. Current
                        options are 'notch' for a notch filter or 'lp' for a
                        lowpass filter.
  --motion-filter-order MOTION_FILTER_ORDER
                        number of filter coeffecients for the band-stop
                        filter.
  --band-stop-min BAND_STOP_MIN
                        lower frequency (bpm) for the band-stop motion filter.
  --band-stop-max BAND_STOP_MAX
                        upper frequency (bpm) for the band-stop motion filter.
  --motion-filter-option MOTION_FILTER_OPTION
                        determines direction(s) in which to filter respiratory
                        artifact. Default is all directions.
  --physio PHYSIO       input .tsv file containing physio data to
                        automatically determine motion filter parameters.
                        Columns, start time, and frequency will also need to
                        be specified. NOT IMPLEMENTED.

final concatenation:
  final stage parameters for after setup and tasks. Concatenates, parcellates,
  and saves combined FD numbers.

  --teardown            flag to run final concatenation steps. After tasks
                        have completed, concatenate resting state data and
                        parcellate.
  --tasklist TASKLIST   comma delimited tasks to be concatenated, pass in
                        argument multiple times to add more task lists. Also
                        determines which tasks will be parcellated, so a
                        single task may be input to parcellate it. Required
                        for this stage. May be a list of one.
```

## Overview

The script is run with calls to three "modes":

### --setup

creates ventricular and wm masks from the anatomical segmentations, and 
computes mean time courses in these classes for use in bold regression later 
on.

### --task TASKNAME

for each fmri run, this script is called to perform regressions on 
motion, ventricular and white matter signals, as well as mean signal 
regression. Framewise displacement (FD) is calculated on the motion numbers 
and used for regression, but also saved for available use in FD thresholding.
If motion band-stop parameters are specified, the motion numbers are first 
filtered in each spatial dimension then FD is computed. The resulting time 
series is saved along with a 'grayplot' displaying relevant time series data.

### --teardown

concatenates any resting state data which shares the same bids task name, 
also concatenates any FD numbers and saves out a matlab (hdf5) file with 
various FD threshold masks computed.

## References

[1] Damien A. Fair, Oscar Miranda-Dominguez, et al. Correction of respiratory artifacts in MRI head motion estimates. NeuroImage, Volume 208, 2020, [doi:10.1016/j.neuroimage.2013.08.048](https://doi.org/10.1016/j.neuroimage.2019.116400).

[2] Power J, et al. Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage [Internet]. Elsevier Inc.; 2014 Jan 1 [cited 2014 Jul 9];84:32041. doi:10.1016/j.neuroimage.2013.08.048

[3] Friston KJ, et al. Movement-related effects in fMRI time-series. Magn Reson Med [Internet]. 1996;35(3):34655. doi:10.1002/mrm.1910350312

