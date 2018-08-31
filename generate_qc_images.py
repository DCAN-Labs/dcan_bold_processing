#!/usr/bin/env python3

import argparse
import os
import re
import subprocess


def _cli():
    parser = generate_parser()
    args = parser.parse_args()

    return interface()


def generate_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser(
            prog='generate_qc_images.py',
            description="""
            Generates a set of qc images for use in executive summary.
            """
        )
    parser.add_argument()
    parser.add_argument()
    parser.add_argument()
    parser.add_argument()
    parser.add_argument()

    return parser


def interface():
    pass


def scene_from_template(image_spec, scene_template, output):
    with open(scene_template) as fd:
        scene = fd.read()
    for name, image in image_spec:
        image_name = os.path.basename(image)
        scene = re.sub('%s_PATH' % name, image, scene)
        scene = re.sub('%s_NAME' % name, image_name, scene)
    with open(output, 'w') as fd:
        fd.write(scene)


def slice_to_png(scene, slice_number, output, **kwargs):
    defaults = {
        'dimx': 900,
        'dimy': 800,
        'stdout': subprocess.DEVNULL,
        'stderr': subprocess.DEVNULL
    }
    kwargs.update(defaults)
    kwargs.update({
        'scene': scene,
        'output': output,
        'slice_number': slice_number
    })
    cmdfmt = '{wb_command} -show-scene {scene} {slice_number} {output} ' \
             '{dimx} {dimy}'
    cmd = cmdfmt.format(**kwargs)
    subprocess.call(cmd, stderr=kwargs['stderr'], stdout=kwargs['stdout'])


def dense_png():
    # loop over slice_to_png
    pass


def silhouetted_mosaic():
    # creates X_in_Y files ("slices" tool wrapper)
    pass


if __name__ == '__main__':
    _cli()
