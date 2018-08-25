# dcan signal processing

\*\*
This is a repository for the dcan labs bold signal processing. It is 
forked from FNL\_preproc and is meant to take its place.
\*\*

This code repository consists of python wrappers and matlab scripts
for signal processing of the bold signal extracted from fMRI data.
This program is designed for the explicit output data of the HCP
fMRI pipeline or its DCAN derivative. It is not intended for use with
other preprocessed data, so use at your own peril.

## dcan\_signal\_processing.py

main wrapper for signal processing scripts.

```{bash}
usage: dcan_signal_processing.py [-h] --subject SUBJECT --task TASK [OPTIONS]

Wraps the compiled DCAN Signal Processing Matlab script, version: v4.0.0. Runs
in 3 main modes: [setup], [task], and [teardown]. [setup]: creates white
matter and ventricular masks for regression, must be run prior to task.
[task]: runs regressions on a given task/fmri and outputs a corrected
dtseries, along with power 2014 motion numbers in an hdf5 (.mat) format file.
[teardown]: concatenates any resting state runs into a single dtseries.

optional arguments:
  -h, --help            show this help message and exit
  --subject SUBJECT     subject/participant id (default: None)
  --task TASK           name of fmri data as used in the dcan fmri pipeline.
                        For bids data it is set to "task-NAME" (default: None)
  --output-folder OUTPUT_FOLDER
                        output folder which contains all files produced by the
                        dcan fmri-pipeline. Used for setting up standard
                        inputs and outputs (default: None)
  --fd-threshold FD_THRESHOLD
                        upper frame-wise displacement threshold for use in
                        signal regression. (default: 0.3)
  --filter-order FILTER_ORDER
                        number of filter coefficients for bold bandpass
                        filter. (default: 2)
  --lower-bpf LOWER_BPF
                        lower cut-off frequency (Hz) for the butterworth
                        bandpass filter. (default: 0.009)
  --upper-bpf UPPER_BPF
                        upper cut-off frequency (Hz) for the butterworth
                        bandpass filter. (default: 0.08)
  --motion-filter-type MOTION_FILTER_TYPE
                        type of band-stop filter to use for removing
                        respiratory artifact from motion regressors. Current
                        options are 'notch' for a notch filter or 'lp' for a
                        lowpass filter, or None for singleband or slow TR
                        data. (default: None)
  --physio PHYSIO       input .tsv file containing physio data to
                        automatically determine motion filter parameters.
                        Columns, start time, and frequency will also need to
                        be specified. NOT IMPLEMENTED. (default: None)
  --motion-filter-option MOTION_FILTER_OPTION
                        determines direction(s) in which to filter respiratory
                        artifact. (default: 5)
  --motion-filter-order MOTION_FILTER_ORDER
                        number of filter coeffecients for the band-stop
                        filter. (default: 4)
  --band-stop-min BAND_STOP_MIN
                        lower frequency (bpm) for the band-stop motion filter.
                        (default: None)
  --band-stop-max BAND_STOP_MAX
                        upper frequency (bpm) for the band-stop motion filter.
                        (default: None)
  --skip-seconds SKIP_SECONDS
                        number of seconds to cut off the beginning of fmri
                        time series. (default: 5)
  --setup               prepare white matter and ventricle masks, must be run
                        prior to individual task runs. (default: False)
```

## generate\_qc\_images.py

main wrapper for quality control images which are fed into executive
summary.

