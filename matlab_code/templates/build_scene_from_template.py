#!/usr/bin/env python2.7

import os
import sys
import argparse
import re
import shutil


"""
input arguments: template scene file
                 list of template fields that need to be modified
                 list of paths to corresponding subject images

1) Create dictionary of template fields to subject paths
2) Read template scene file into memory
3) Make all necessary substitutions
4) Write the subject's scene file

"""




def build_scene_from_template(scene_file, paths_dict):
    with open (scene_file, 'r+') as scene:
        data = scene.read()
        for field in paths_dict:
            data = re.sub(field + "_NAME_and_PATH", paths_dict[field], data)
            data = re.sub(field + "_NAME", paths_dict[field].split('/')[-1], data)
        scene.seek(0)
        scene.write(data)
        scene.truncate()
        scene.close()

def create_scene_from_template(scene_file):
    with open (scene_file, 'r') as scene:
        data = scene.read()
        
        

def main(argc=sys.argv):
    arg_parser = argparse.ArgumentParser(description="Create subject specific scene file from template")
    arg_parser.add_argument('-o', '--output-dir', metavar='OUTPUT_DIR', action='store', required=True,
                            help=("Full path to the subjects processed directory"),
                            dest='output_dir')
    args = arg_parser.parse_args()
    output_dir = args.output_dir
    subject_id = output_dir.split('/')[-1]

    template_fields = ['TX_IMG', 'R_PIAL', 'L_PIAL', 'R_WHITE', 'L_WHITE']

    subject_paths = [output_dir + "/MNINonLinear/T1w_restore.nii.gz",
                     output_dir + "/MNINonLinear/fsaverage_LR32k/" + subject_id + ".R.white.32k_fs_LR.surf.gii",
                     output_dir + "/MNINonLinear/fsaverage_LR32k/" + subject_id + ".R.pial.32k_fs_LR.surf.gii",
                     output_dir + "/MNINonLinear/fsaverage_LR32k/" + subject_id + ".L.white.32k_fs_LR.surf.gii",
                     output_dir + "/MNINonLinear/fsaverage_LR32k/" + subject_id + ".L.pial.32k_fs_LR.surf.gii"]

    paths_dict = dict(zip(template_fields, subject_paths))

    scene_file = output_dir + "/t1_169_scene.scene"
    if os.path.exists(scene_file):
        os.remove(scene_file)
    shutil.copy2("/home/exacloud/lustre1/fnl_lab/projects/andersperrone/FNL_preproc_v3/templates/parasagittal_Tx_169_template.scene", scene_file)

    build_scene_from_template(scene_file, paths_dict)

if __name__ == "__main__":
    sys.exit(main())



