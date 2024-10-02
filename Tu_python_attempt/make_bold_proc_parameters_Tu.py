#!/usr/bin/env python3

__prog__ = 'DCANBOLDProc'
__version__ = 'Tu1.0'
__doc__ = \
'''
Instead of making a parameter for each run, we make a parameter file that will be applied to all files so that we can share the parameters with the paper.
Default is to save to the .json file to the current directory unless given the '--output-folder' argument.
'''

import argparse
import json
import os
from datetime import datetime

here = os.path.dirname(os.path.realpath(__file__))
today_date = datetime.today().strftime('%Y%m%d')

def _cli():
    """
    command line interface
    :return:
    """
    parser = generate_parser()
    args = parser.parse_args()

    kwargs = {
        'output_folder': args.output_folder,
        'fd_threshold': args.fd_threshold,
        'contiguous_frames': args.contiguous_frames,
        'filter_order': args.filter_order,
        'lower_bpf': args.lower_bpf,
        'upper_bpf': args.upper_bpf,
        'motion_filter_type': args.motion_filter_type,
        'motion_filter_option': args.motion_filter_option,
        'motion_filter_order': args.motion_filter_order,
        'band_stop_min': args.band_stop_min,
        'band_stop_max': args.band_stop_max,
        'skip_seconds': args.skip_seconds,
        'brain_radius': args.brain_radius,
        'FD_type':args.FD_type,
        'GSR':args.GSR
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
            prog='make_bold_proc_parameters.py',
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-v', '--version', action='version',
        version='%s_v%s' % (__prog__, __version__),
        help='print the software name and version'
    )
    
    parser.add_argument(
    '--output-folder',default=here,
    help='output folder to store this parameter files')

    bold_filter = parser.add_argument_group(
        'bold signal filtering',  'bold signal filtering parameters.')
    bold_filter.add_argument(
        '--filter-order', type=int, default=2,
        help='number of filter coefficients for butterworth bandpass filter.'
    )
    bold_filter.add_argument(
        '--lower-bpf', type=float, default=0.009,
        help='lower cut-off frequency (Hz) for the butterworth bandpass '
             'filter.'
    )
    bold_filter.add_argument(
        '--upper-bpf', type=float, default=0.080,
        help='upper cut-off frequency (Hz) for the butterworth bandpass '
             'filter.'
    )
    fd = parser.add_argument_group(
        'framewise displacement', 'parameters related to computation of '
                                  'framewise displacment (FD)')
    fd.add_argument(
        '--fd-threshold', type=float, default=0.3,
        help='upper framewise displacement threshold for use in signal '
             'regression.'
    )
    fd.add_argument(
        '--skip-seconds', type=int, default=5,
        help='number of seconds to cut off the beginning of fmri time series.'
    )
    fd.add_argument(
        '--contiguous-frames', type=int, default=5,
        help='number of contigious frames for power 2014 fd thresholding.'
    )
    fd.add_argument(
        '--brain-radius', type=int, default=5,
        help='radius of brain (in mm) for computation of rotational displacements'
    )
    fd.add_argument(
        '--FD-type', type=int, default=1,
        help='FD type (1 = L1-norm, 2 = L2-norm)'
    )
    fd.add_argument(
        '--GSR', type=int, default=1,
        help='Global Signal Regression (1 = YES, 0 = NO)'
    )
    fd.add_argument(
        '--motion-filter-type', choices=['notch','lp'], default=None,
        help='type of band-stop filter to use for removing respiratory '
             'artifact from motion regressors. Current options are \'notch\' '
             'for a notch filter or \'lp\' for a lowpass filter.'
    )
    fd.add_argument(
       '--motion-filter-order', type=int, default=4,
       help='number of filter coeffecients for the band-stop filter.'
    )
    fd.add_argument(
        '--band-stop-min', type=float_or_None,
        help='lower frequency (bpm) for the band-stop motion filter.'
    )
    fd.add_argument(
        '--band-stop-max', type=float_or_None,
        help='upper frequency (bpm) for the band-stop motion filter.'
    )
    fd.add_argument(
       '--motion-filter-option', type=int, default=5,
       help='determines direction(s) in which to filter respiratory '
            'artifact. Default is all directions.'
    )
    return parser

def interface(output_folder, fd_threshold=None,
              filter_order=None, lower_bpf=None, upper_bpf=None,
              motion_filter_type=None, motion_filter_option=None,
              motion_filter_order=None, band_stop_min=None,
              band_stop_max=None, skip_seconds=None, brain_radius=None,
              contiguous_frames=None, FD_type=None, GSR=1, **kwargs):
    """
    main function with 3 modes:
        setup, task, and teardown.

    setup:
    generates white matter and ventricular masks.

    task:
    Runs filtered movement regressors, calculates mean signal
    in ventricles and white matter, then calls dcan signal processing matlab
    script.

    teardown:
    concatenates resting state data and creates parcellated time series.
    :param output_folder: base output files folder to store this json file
    :param legacy_tasknames: support legacy tasknames with run index appended to task, e.g. "task-rest01")
    :param task: name of task
    :param fd_threshold: threshold for use in signal regression
    :param filter_order: order of bold signal bandpass filter
    :param lower_bpf: lower limit of bold signal bandpass filter
    :param upper_bpf: upper limit of bold signal bandpass filter
    :param motion_filter_type: type of bandstop filter for filtering motion
    regressors.  Default: 'notch'
    :param motion_filter_option: dimensions along which to filter motion.
    Default: 1 1 1 1 1 1 (all translations and rotations)
    :param motion_filter_order: bandstop filter order
    :param band_stop_min: lower limit of motion bandstop filter
    :param band_stop_max: upper limit of motion bandstop filter
    :param skip_seconds: number of seconds to cut of beginning of task.
    :param brain_radius: radius for estimation of angular motion regressors
    :param contiguous_frames: minimum contigious frames for fd thresholding.
    :param kwargs: additional parameters.  Can be used to override default
    paths of inputs and outputs.
    :return:
    """
    # name should only reflect release version, not filter usage.
    version_name = '%s_v%s' % (__prog__, __version__)

    output_spec = {
        'config': os.path.join(output_folder, 
                               'config_%s.json' % today_date),
    }
    output_spec.update(kwargs.get('output_spec', {}))

    # check integrity of filter parameters:
    if lower_bpf and upper_bpf:
        assert lower_bpf < upper_bpf, \
            'lower bandpass limit exceeds upper limit.'
    if band_stop_min and band_stop_max:
        assert band_stop_min < band_stop_max, \
            'lower bandstop limit exceeds upper limit.'    
    
    # run signal processing on dtseries
    bold_processing_input = {
        'bp_order': filter_order,
        'lp_Hz': lower_bpf,
        'hp_Hz': upper_bpf,
        'fd_th': fd_threshold,
        'skip_seconds': skip_seconds,
        'motion_filter_type':motion_filter_type,
        'motion_filter_order':motion_filter_order,
        'motion_filter_option':motion_filter_option,
        'band_stop_min':band_stop_min,
        'band_stop_max':band_stop_max,
        'brain_radius_mm':brain_radius, # in mm
        'contiguous_frames':contiguous_frames,
        'FD_type':FD_type,
        'GSR':GSR,
        'version_name':version_name,
    }
    # write input json for matlab script
    with open(output_spec['config'], 'w') as fd:
        json.dump(bold_processing_input, fd, sort_keys=True, indent=4)
    print('Write json to be used for all files')

def float_or_None(x):
    if x.lower() == 'none':
        return None
    else:
        return float(x)

if __name__ == '__main__':
    _cli()
