#!/usr/bin/env python3

__prog__ = 'DCANBOLDProc'
__version__ = 'Tu1.0'
__doc__ = \
"""Python-based BOLD signal processing adapted from https://github.com/DCAN-Labs/dcan_bold_processing, version: %s.
# original documentation: section "DCANBOLD Processing" on https://cdnis-brain.readthedocs.io/infant-doc/
Step 1: Make white matter and ventricular masks, create mean signal from those masks for each run
Step 2: Run BOLD processing including demean/detrend, nuissance regression and bandpass filtering, also generates the FD file (respiratory filtered) 
Step 3: Concatenate BOLD and FD for all runs

""" % __version__
__references__ = \
"""References
----------
[1] Fair DA, Miranda-Dominguez O, et al. Correction of respiratory artifacts
in MRI head motion estimates. bioRxiv [Internet]. 2018 Jan 1; Available from:
http://biorxiv.org/content/early/2018/06/07/337360.abstract
[2] Power J, et al. Methods to detect, characterize, and remove motion
artifact in resting state fMRI. Neuroimage [Internet]. Elsevier Inc.; 2014
Jan 1 [cited 2014 Jul 9];84:32041. doi: 10.1016/j.neuroimage.2013.08.048
[3] Friston KJ, et al. Movement-related effects in fMRI time-series. Magn
Reson Med [Internet]. 1996;35(3):34655. doi: 10.1002/mrm.1910350312
"""

import argparse
import os
import subprocess
import json

here = os.path.dirname(os.path.realpath(__file__))

def _cli():
    """
    command line interface
    :return:
    """
    parser = generate_parser()
    args = parser.parse_args()

    kwargs = {
        'subject': args.subject if args.subject else os.path.basename(os.path.dirname(args.input_folder)),
        'folderkeyword': args.folderkeyword,
        'input_folder':args.input_folder,
        'output_folder': args.output_folder,
        'fmri_res': args.fmri_res,
        'roi_res': args.roi_res,
        'no_aparc': args.no_aparc,
        'config':args.config
    }

    return interface(**kwargs)

def generate_parser(parser=None):
    """
    generates argument parser for this program.
    :param parser: if set, args are added to this parser.
    :return: ArgumentParser
    """
    if not parser:
        parser = argparse.ArgumentParser(
            prog='dcan_bold_proc_Tu.py',
            description=__doc__,
            epilog=__references__,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-v', '--version', action='version',
        version='%s_v%s' % (__prog__, __version__),
        help='print the software name and version'
    )

    parser.add_argument(
        '--subject', required=False, help='subject/participant id' # Used to name the masks, default to the input folder name
    )
    
    parser.add_argument(
        '--folderkeyword', required=False, help='task folders keyword',default='' # Used to loop through tasks, if left blank this will loop through all subfolders under /MNINonlinear/Results
    )

    parser.add_argument(
        '--output-folder', required = True,
        help='output folder which contains all files produced by the dcan '
             'fmri-pipeline.  Used for setting up standard outputs'
    )

    parser.add_argument(
        '--input-folder', required=False,
        help='input folder which contains all files produced by the dcan '
             'fmri-pipeline.  Used for setting up standard inputs (incl. *nii.gz, *Atlas_dtseries.nii, *Movement_Regressors.txt), will use the output_folder if not specified'
    )

    parser.add_argument(
        '--config',help='.json file that stores the bold processing parameters',default='./config_20240913.json'
    )

    setup =  parser.add_argument_group(
        'wm/vent regressors',  'options for obtaining white matter and '
                               'ventricle signal regressors')
    setup.add_argument(
        '--fmri-res', type=float, default=2.,
        help='isotropic resolution (mm) for final fmri volume. Default is 2.'
    )
    setup.add_argument(
        '--roi-res', type=float, default=2.,
         help='isotropic resolution (mm) for vent/wm roi volumes. Default is '
              '2.'
    )
    setup.add_argument(
        '--no-aparc', action='store_true',
        help='flag to use original freesurfer LR white matter labels instead '
             'of parcellated labels.'
    )
    return parser

def interface(subject, output_folder, input_folder=None,folderkeyword='', fmri_res=2., roi_res=2., no_aparc=False,
              **kwargs):
    """
    Generates white matter and ventricular masks, calculates mean signal
    in ventricles and white matter for each task folder, then run the motion confound processing

    :param subject: subject id
    :param output_folder: base output files folder for fmri pipeline
    :param input_folder: base input files folder for fmri pipeline
    :param folderkeyword: name of the keyword that specifies string contained in all folders we wanted to process (e.g. 'rfMRI')
    :param kwargs: additional parameters.  Can be used to override default
    paths of inputs and outputs.
    :return:
    """
    # name should only reflect release version, not filter usage.
    version_name = '%s_v%s' % (__prog__, __version__)
    if not input_folder:
        print('No input folder provided, using output folder as input folder.')
        input_folder = output_folder
    
    # Load the configuration json file
    with open(output_spec['config'], 'r') as file:
        json_input = json.load(file)   

    for task in os.listdir(os.path.join(input_folder, 'MNINonLinear','Results')):
        task_path = os.path.join(input_folder, 'MNINonLinear','Results',task)    
        print(task_path)
        if os.path.isdir(task_path) and folderkeyword in task:
            print(f"Processing task: {task}")
            
            # standard input and output folder locations. # Change this part if you want your input/output folder structure to differ
            input_spec = {
                'dtseries': os.path.join(input_folder, 'MNINonLinear', 'Results',
                                        task, '%s_Atlas.dtseries.nii' % task),
                'fmri_volume': os.path.join(input_folder, 'MNINonLinear', 'Results',
                                            task, '%s.nii.gz' % task),
                'movement_regressors': os.path.join(input_folder, 'MNINonLinear',
                                                    'Results', task,
                                                    'Movement_Regressors.txt'),
                'segmentation': os.path.join(input_folder, 'MNINonLinear', 'ROIs',
                                            'wmparc.%g.nii.gz' % roi_res)
            }
            input_spec.update(kwargs.get('input_spec', {}))
            output_spec = {
                'config': os.path.join(output_folder, 'MNINonLinear', 'Results',
                                    version_name,
                                    '%s_mat_config.json' % version_name),
                'output_dtseries_basename': '%s_%s_Atlas' % (task, version_name), # this was in input_spec but I think it makes more sense in the output_spec # JCT 20240908
                'output_motion_numbers': os.path.join(output_folder, 'MNINonLinear',
                                                    'Results', task, version_name,
                                                    'motion_numbers.txt'),
                'output_timecourses': os.path.join(output_folder, version_name,
                                            'analyses_v2','timecourses'),
                'result_dir': os.path.join(output_folder, 'MNINonLinear', 'Results',
                                        task, version_name),
                'summary_folder': os.path.join(output_folder,
                                            'summary_%s' % version_name),
                'vent_mask': os.path.join(output_folder, 'MNINonLinear',
                                        'vent_%gmm_%s_mask_eroded.nii.gz' % \
                                            (fmri_res, subject)),
                'vent_mean_signal': os.path.join(output_folder, 'MNINonLinear',
                                                'Results', task, version_name,
                                                '%s_vent_mean.txt' % task),
                'wm_mask': os.path.join(output_folder, 'MNINonLinear',
                                        'wm_%gmm_%s_mask_eroded.nii.gz' % \
                                            (fmri_res, subject)),
                'wm_mean_signal': os.path.join(output_folder, 'MNINonLinear',
                                            'Results', task, version_name,
                                            '%s_wm_mean.txt' % task)
            }
            output_spec.update(kwargs.get('output_spec', {}))
            
            # Check for inputs
            check_files_exist(input_spec)

            if json_input['GSR']: # Don't worry about white matter and ventricular signal if we are not using it for motion confound regression.
                if ((os.path.isfile(output_spec['wm_mean_signal'])) and (os.path.isfile(output_spec['vent_mean_signal']))):
                    print("White matter and ventricular mean signal found.")
                elif not ((os.path.isfile(output_spec['wm_mask'])) and (os.path.isfile(output_spec['vent_mask']))):
                    print("Generating the white matter and ventricular masks.")
                    if no_aparc:
                        label_override = dict(wm_lt_R=2, wm_ut_R=2, wm_lt_L=41,
                                                wm_ut_L=41)
                    else:
                        label_override = {}
                    
                    # create white matter and ventricle masks for regression (will be used for all runs)
                    make_masks(input_spec['segmentation'], output_spec['wm_mask'],
                                output_spec['vent_mask'], fmri_res=fmri_res,
                                roi_res=roi_res, **label_override)
                
                # Create the result_dir if it's not there
                if not os.path.exists(output_spec['result_dir']):
                    os.mkdir(output_spec['result_dir'])

                # Obtain mean ventricular and white matter signals
                if not os.path.exists(output_spec['wm_mask']):
                    mean_roi_signal(input_spec['fmri_volume'], output_spec['wm_mask'],
                                    output_spec['wm_mean_signal'], fmri_res, roi_res)
                if not os.path.exists(output_spec['vent_mask']):
                    mean_roi_signal(input_spec['fmri_volume'], output_spec['vent_mask'],
                                    output_spec['vent_mean_signal'], fmri_res, roi_res)        
            
            # Search for the path of Connectome Workbench
            path_wb_c = subprocess.run(['which', 'wb_command'], stdout=subprocess.PIPE, text=True)
            path_wb_c = path_wb_c.stdout.strip()

            json_input['path_wb_c'] = path_wb_c # we can prob do it with nibabel

            # Run signal processing
            dcan_signal_processing_Tu(json_input) # This function was not actually finished, it was originally meant to be the .py version of the dcan_signal_processing_developing.ipynb 
            # but I realized that they have a newer public version https://github.com/PennLINC/xcp_d which basically does what I am doing here so I'm not planning to finish it.

def check_files_exist(input_spec):
    missing_files = [file for file in input_spec.values() if not os.path.isfile(file)]
    if missing_files:
        print(f"Missing files: {', '.join(missing_files)}")
    # assert not missing_files, f"Missing files: {', '.join(missing_files)}"

def make_masks(segmentation, wm_mask_out, vent_mask_out, **kwargs):

    """
    generates ventricular and white matter masks from a Desikan/FreeSurfer
    segmentation file.  label constraints may be overridden.
    :param segmentation: Desikan/FreeSurfer spec segmentation nifti file.
    Does not need to be a cifti but must have labels according to FS lookup
    table, including cortical parcellations.
    :param wm_mask_out: binary white matter mask.
    :param vent_mask_out: binary ventricular mask.
    :param kwargs: dictionary of label value overrides.  You may override
    default label number bounds for white matter and ventricle masks in the
    segmentation file.
    :return: None
    """

    wd = os.path.dirname(wm_mask_out)
    # set parameter defaults
    defaults = dict(wm_lt_R=2950, wm_ut_R=3050, wm_lt_L=3950, wm_ut_L=4050,
                    vent_lt_R=43, vent_ut_R=43, vent_lt_L=4, vent_ut_L=4,
                    roi_res=2)
    # set temporary filenames
    tempfiles = {
        'wm_mask_L': os.path.join(wd, 'tmp_left_wm.nii.gz'),
        'wm_mask_R': os.path.join(wd, 'tmp_right_wm.nii.gz'),
        'vent_mask_L': os.path.join(wd, 'tmp_left_vent.nii.gz'),
        'vent_mask_R': os.path.join(wd, 'tmp_right_vent.nii.gz'),
        'wm_mask': os.path.join(wd, 'tmp_wm.nii.gz'),
        'vent_mask': os.path.join(wd, 'tmp_vent.nii.gz')
    }
    # inputs and outputs
    iofiles = {
        'segmentation': segmentation,
        'wm_mask_out': wm_mask_out,
        'vent_mask_out': vent_mask_out
    }
    # command pipeline
    cmdlist = [
        'fslmaths {segmentation} -thr {wm_lt_R} -uthr {wm_ut_R} {wm_mask_R}',
        'fslmaths {segmentation} -thr {wm_lt_L} -uthr {wm_ut_L} {wm_mask_L}',
        'fslmaths {wm_mask_R} -add {wm_mask_L} -bin {wm_mask}',
        'fslmaths {wm_mask} -kernel gauss {roi_res:g} -ero {wm_mask_out}',
        'fslmaths {segmentation} -thr {vent_lt_R} -uthr {vent_ut_R} '
        '{vent_mask_R}',
        'fslmaths {segmentation} -thr {vent_lt_L} -uthr {vent_ut_L} '
        '{vent_mask_L}',
        'fslmaths {vent_mask_R} -add {vent_mask_L} -bin {vent_mask}',
        'fslmaths {vent_mask} -kernel gauss {roi_res:g} -ero {vent_mask_out}'
    ]

    # get params
    defaults.update(kwargs)
    kwargs.update(defaults)
    kwargs.update(iofiles)
    kwargs.update(tempfiles)
    # format and run commands
    for cmdfmt in cmdlist:
        cmd = cmdfmt.format(**kwargs)
        subprocess.call(cmd.split())
    # cleanup
    for key in tempfiles.keys():
        os.remove(tempfiles[key])

def mean_roi_signal(fmri, mask, output, fmri_res=2., roi_res=2.):
    """
    :param fmri: path to fmri nifti
    :param mask: path to mask/roi nifti
    :param output: output text file of time series of mean values within the
    mask/roi
    :return: None
    """
    cmd = 'fslmeants -i {fmri} -o {output} -m {mask}'
    if fmri_res != roi_res:
        resamplecmd = 'flirt -interp nearestneighbour -in {mask} -ref ' \
                      '{fmri} -applyxfm -out {mask}'
        resamplecmd = resamplecmd.format(fmri=fmri, output=output, mask=mask,
                                         fmri_res=fmri_res)
        subprocess.call(resamplecmd.split())
    cmd = cmd.format(fmri=fmri, output=output, mask=mask)
    subprocess.call(cmd.split())

def float_or_None(x):
    if x.lower() == 'none':
        return None
    else:
        return float(x)

if __name__ == '__main__':
    _cli()
