#!/usr/bin/env python3

import argparse
import os
import subprocess


def _cli():
    parser = generate_parser()
    args = parser.parse_args()

    return interface()


def generate_parser(parser=None):
    # this script may be deprecated.
    if not parser:
        parser = argparse.ArgumentParser(
            prog='concatenate_dtseries.py',
            description="""
            Concatenates a set of preprocessed dense time series files, 
            creates an hdf5 (.mat) structure describing power 2014 frame 
            censoring, and parcellates them based on the roi-sets folder 
            contents.  The folder specification for roi-sets outlined in the 
            README.
            """
        )
    parser.add_argument('--subject')
    parser.add_argument('--tasklist', required=True)
    parser.add_argument('--brain-radius', type=int,
                        help='approximate brain radius (mm) for subject(s).  '
                             'Used to estimate frame-wise displacement caused '
                             'by rotations.'
                        )
    parser.add_argument('--min-contiguous-frames', type=int,
                        help='minimum number of adjacent frames for motion '
                             'censor.  Frames which meet the frame-wise '
                             'displacement threshold but for less than this '
                             'number contiguously, will be censored.')
    parser.add_argument('--roi-sets')

    return parser


def interface():
    pass
