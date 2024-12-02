"""
This script takes the output from ABCD-HCP-pipeline, load the minimally preprocessed CIFTI file and apply confound removal
Adapted from the MATLAB script dcan_bold_processing.m (https://github.com/DCAN-Labs/dcan_bold_processing/blob/main/matlab_code/dcan_signal_processing.m)
All the dependencies are included in this file (no other script is required)

Potential issue: the GSR version is slightly different from the saved file in the original folder (ABCD from DCAN) potentially because of the differences in interpolation and padding?

Args:
Use python DCANBOLD_reprocess.py -h to see the arguments

Requirements:
Python 3 (3.8.10)
numpy (1.24.4)
scipy (1.10.1)
matplotlib (3.7.5)
nibabel (5.2.1)

Author:
    Jiaxin Cindy Tu (tu.j@wustl.edu)

Date: 2024-12-01

ABCD File Structure on Daenerys:
/data/Daenerys/ABCD/data/abcd_collection3165/derivatives/abcd-hcp-pipeline/
- <subject_id>
    -ses-baselineYear1Arm1
        -files
            -MNINonLinear
                -Results
                    - <task_id> e.g. task-rest01
                        - DCANBOLDProc_v4.0.0
                            - DCANBOLDProc_v4.0.0_mat_config.json <- we will use this to get all required files                              
"""

# Dependency packages
import os
import json
from pathlib import Path
from re import search
import argparse
import logging
import numpy as np
import scipy
import matplotlib.pyplot as plt
import nibabel as nib


def main():
    '''
    Main function to run the processing
    '''
    parser = argparse.ArgumentParser(description="confound removal for BOLD fMRI data")
    parser.add_argument(
        "--subject-id",
        type=str,
        help="The folder name for subject, e.g. sub-xxxxxx",
        required=True,
    )
    parser.add_argument(
        "--task-id",
        type=str,
        help="The folder name for run e.g. task-rest01",
        required=True,
    )
    parser.add_argument(
        "--data-path",
        type=str,
        help="DCAN ABCD-HCP-pipeline output folder",
        required=True,
    )
    parser.add_argument(
        "--results-path",
        type=str,
        help="results folder to save cleaned BOLD data",
        required=True,
    )
    parser.add_argument(
        "--json-path",
        type=str,
        help="full path to json file for parameters and input file names",
        required=True,
    )
    parser.add_argument(
        "--savefig",
        type=int,
        help="save the diagnostic plots (1) or not (0), default = 1",
        default=1
    )
    parser.add_argument(
        "--GSR",
        type=int,
        help="include Global Signal Regression (CSF, WM, GM, and their derivatives) [1] or not [0]",
        required=True,
    )
    parser.add_argument(
        "--FD-type", type=int, default=1, help="L1 [1] (default) or L2 [2]"
    )
    parser.add_argument(
        "--brain-radius-mm",
        type=int,
        default=50,
        help="ballpark estimate of brain radius in mm, default = 50, N.B. this is for FD and movement regressor calculation", 
        # but it is possible (but unlikely because they hardcoded 50 mm in the previous versions of the processing) that they used a different brain radius for the file_mov_reg when applying bandstop filters but we don't know that number,
        # and even so the effect should be negligible"
    )
    args = parser.parse_args()
    subject_id = args.subject_id
    task_id = args.task_id
    data_path = args.data_path
    results_path = args.results_path

    # Part I: determine the parameters

    # Load the original json file
    with open(args.json_path, "r") as file:
        json_input = json.load(file)

    # Update paths
    json_input["FD_type"] = args.FD_type
    json_input["brain_radius_mm"] = args.brain_radius_mm
    json_input["GSR"] = args.GSR
    if json_input["GSR"]:
        json_input["result_dir"] = os.path.join(
            results_path, subject_id, task_id, "DCANBOLDProc_v4.0.0", "gsr"
        )
    else:
        json_input["result_dir"] = os.path.join(
            results_path, subject_id, task_id, "DCANBOLDProc_v4.0.0", "nogsr"
        )
    json_input["path_cii"] = str.replace(json_input["path_cii"], "/output/", data_path)
    json_input["file_mov_reg"] = str.replace(
        json_input["file_mov_reg"], "/output/", data_path
    )  # This should be the filtered movement parameters (columns 3-6 in degrees) if you use the filtered version
    json_input["file_vent"] = str.replace(
        json_input["file_vent"], "/output/", data_path
    )
    json_input["file_wm"] = str.replace(json_input["file_wm"], "/output/", data_path)
    json_input["config"] = os.path.join(
        json_input["result_dir"], "DCANBOLDProc_v4.0.0_mat_config.json"
    )
    json_input.pop("path_ex_sum")
    json_input.pop("path_wb_c")

    # Save the parameters in the results folder
    directory_path = Path(json_input["result_dir"])
    directory_path.mkdir(parents=True, exist_ok=True)
    with open(json_input["config"], "w") as fd:
        json.dump(json_input, fd, sort_keys=True, indent=4)

    logging.basicConfig(
        filename=os.path.join(json_input["result_dir"], "output.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info("Now processing subject" + subject_id + ", run: " + task_id)

    # Part II: Load and prepare minimally preprocessed data and regressors

    img = nib.load(json_input["path_cii"])
    data = img.get_fdata()
    DVAR_pre_reg = calculate_dvars_from_cifti(data)

    # Global signals (Line 52,58,83 in dcan_bold_processing.m)
    if json_input["GSR"]:
        wm = np.loadtxt(json_input["file_wm"])
        vent = np.loadtxt(json_input["file_vent"])
        WB = np.mean(data, axis=1)  # Global signal in gray matter
        Glob = np.concatenate(
            (wm.reshape(-1, 1), vent.reshape(-1, 1), WB.reshape(-1, 1)), axis=1
        )
        dGlob = np.vstack(([0, 0, 0], np.diff(Glob, axis=0)))
        
        if args.savefig:
            plt.figure(figsize=(8, 4))
            plt.plot(Glob, linewidth=1, label=("white matter", "ventricle", "grayordinate"))
            plt.title("Global signal")
            plt.legend()
            plt.savefig(
                os.path.join(json_input["result_dir"], "globalsignal_trace.png"),
                format="png",
                dpi=300,
            )

    # Get TR (Line 34 in dcan_signal_processing.m)
    TR = nib.load(json_input["path_cii"]).header.get_axis(0).step  # In seconds
    if TR > 20:  # handle niftis that store TR in ms
        TR /= 1000
    logging.info("TR = " + str(TR))

    # Load movement regressors
    assert os.path.isfile(
        json_input["file_mov_reg"]
    ), "movement regressor file not found"
    filtered_mov_regressors = np.loadtxt(json_input["file_mov_reg"])
    FR = make_friston_regressors(
        filtered_mov_regressors, json_input["brain_radius_mm"]
    )  # 24-param Friston Regressors

    # Calculate FD  (Line 73 in dcan_bold_processing.m)
    FD, meanFD = calc_FD(FR[:,:6], FD_type=json_input["FD_type"])
    FD_file_name = os.path.join(
        json_input["result_dir"], str.replace(Path(json_input["file_mov_reg"]).name,'.txt','_FD.txt')
    )
    np.savetxt(FD_file_name, FD)
    logging.info("Mean FD = " + str(meanFD))
    # Now state the frames to keep for the processing (Line 78 in dcan_bold_processing.m)
    keepframe = FD <= json_input["fd_th"]
    skip_frames = np.int8(np.floor(json_input["skip_seconds"] / TR))
    keepframe[:skip_frames] = False
    logging.info("Retained frames = " + str(np.sum(keepframe)))

    if args.savefig:
        plt.figure(figsize=(8,4))
        for j in np.where(keepframe==0)[0]:
            plt.axvline(x=j, color=[0.5,0.5,0.5], alpha=0.5)
        plt.plot(FD,linewidth=1)
        plt.title('FD')
        plt.savefig(os.path.join(json_input['result_dir'],'FD_trace.png'), format='png', dpi = 300)

    # Concatenate regressors
    if json_input["GSR"]:
        R = np.concatenate(
            (Glob, dGlob, FR), axis=1
        )  # with GSR # N.B. this uses 30 parameter, for 36 parameter it would also involve the square of Glob and dGlob
    else:
        R = FR  # without GSR

    # Part III: Regression and bandpass filter to give cleaned BOLD data

    # Demean/detrend regressors and data
    R = R - np.mean(R[keepframe, :], axis=0)
    R = detrend_manual(R, keepframe)
    data_dd = data - np.mean(data[keepframe, :], axis=0)
    data_dd = detrend_manual(data_dd, keepframe)

    if args.savefig:
        # Plot grayplots
        crange=[-200, 200] # I'm using +/-2% as Jonathan Power's papers, but DCAN used +/-6%
        plt.figure(figsize=(8,4))
        im = plt.imshow(data_dd.transpose(), aspect="auto", cmap="gray")
        plt.colorbar(location="right")
        im.set_clim(crange)
        plt.plot(np.where(keepframe==0)[0],np.repeat(0,sum(keepframe==0)),'r|')
        plt.yticks([])
        plt.xlabel('TR')
        plt.savefig(
            os.path.join(json_input["result_dir"], "grayplots_all.png"),
            format="png",
            dpi=300,
        )

        plt.figure(figsize=(8,4))
        X = data_dd.copy()
        X[np.where(keepframe==1)[0]]=np.NaN
        im = plt.imshow(X.transpose(), aspect="auto", cmap="gray")
        plt.colorbar(location="right")
        im.set_clim(crange)
        plt.yticks([])
        plt.xlabel('TR')
        plt.savefig(
            os.path.join(json_input["result_dir"], "grayplots_removed.png"),
            format="png",
            dpi=300,
        )

        plt.figure(figsize=(8,4))
        X = data_dd.copy()
        X[np.where(keepframe==0)[0]]= np.NaN
        im = plt.imshow(X.transpose(), aspect="auto", cmap="gray")
        plt.colorbar(location="right")
        im.set_clim(crange)
        plt.yticks([])
        plt.xlabel('TR')
        plt.savefig(
            os.path.join(json_input["result_dir"], "grayplots_retained.png"),
            format="png",
            dpi=300,
        )

    # Nuisance Regression (Line 121 in dcan_signal_processing.m)
    b, _, _, _ = np.linalg.lstsq(R[keepframe, :], data_dd[keepframe, :], rcond=None)
    data_postreg = data_dd - R @ b
    DVAR_post_reg = calculate_dvars_from_cifti(data_postreg)

    # Linear interpolation before bandpass filter (Line 144 in dcan_signal_processing.m)
    # Interpolate for missing frames
    x = np.where(keepframe)[0]
    x_removed = np.where(keepframe == 0)[0]
    x_outsidebound = (x_removed < x[0]) | (x_removed > x[-1])
    y_removed = np.apply_along_axis(
        lambda col: np.interp(x_removed, x, col[keepframe]), axis=0, arr=data_postreg
    )
    # Replace extrapolated points with the mean of retained low-motion data
    y_mean = np.mean(data_postreg[keepframe, :], axis=0)
    y_removed[x_outsidebound, :] = y_mean

    data_interpolated = data_postreg.copy()
    data_interpolated[keepframe == 0, :] = y_removed

    # Bandpass filter with manual zero padding
    fs = 1 / TR  # Sampling frequency
    fNy = fs / 2  # Nyquist frequency
    b_filt, a_filt = scipy.signal.butter(
        json_input["bp_order"] / 2,
        np.array([json_input["lp_Hz"], json_input["hp_Hz"]]) / fNy,
        "bandpass",
    )

    # Zero-pad the data for filtering by concatenating rows of zeros on either side of the data
    padding = np.zeros_like(
        data_interpolated
    )  # Create a padding array of the same shape as Rr_int
    pad_amt = padding.shape[0]  # Number of rows to pad

    # Concatenate padding rows on top and bottom of Rr_int
    temp = np.vstack((padding, data_interpolated, padding))

    # Apply the filtfilt function (zero-phase filtering)
    data_filtered = scipy.signal.filtfilt(
        b_filt, a_filt, temp, axis=0, padtype=None
    )  
    
    # Apply filtering along the rows (axis=0)
    data_filtered = data_filtered[pad_amt:-pad_amt]
    DVAR_post_filter = calculate_dvars_from_cifti(data_filtered)

    if args.savefig:
        plt.figure(figsize=(8, 4))
        plt.plot(DVAR_pre_reg, linewidth=1, label="DVARS pre regression", color="b")
        plt.plot(DVAR_post_reg, linewidth=1, label="DVARS post regression", color="r")
        plt.plot(DVAR_post_filter, linewidth=1, label="DVARS post filtered", color="g")
        plt.xlabel("TR")
        plt.legend()
        plt.savefig(
            os.path.join(json_input["result_dir"], "DVARS_trace.png"), format="png", dpi=300
        )

    if json_input["GSR"]:
        logging.info("Saving cleaned data (w/GSR)")
    else:
        logging.info("Saving cleaned data (wo/GSR)")

    # Save filtered data
    ax1 = img.header.get_axis(0)
    ax2 = img.header.get_axis(1)
    header = (ax1, ax2)
    output_img = nib.cifti2.cifti2.Cifti2Image(np.single(data_filtered), header)
    output_img.to_filename(
        os.path.join(
            json_input["result_dir"],
            json_input["FNL_preproc_CIFTI_basename"] + ".dtseries.nii",
        )
    )


# All dependency functions
def calculate_dvars_from_cifti(data):
    """
    This function calculates DVARS (Derivative of Variance) based on grayordinates (WM and non-brain excluded).

    Parameters:
    data (ndarray): 2D numpy array with shape (tr, g), where g represents the number of grayordinates and tr is the number of time points.

    Returns:
    dvars (float): The calculated DVARS value.
    """
    # Check size and transpose if needed
    num_timepoints, num_grayordinates = data.shape
    if num_grayordinates < num_timepoints:
        data = data.T
        print("data transposed due to timepoints > grayordinates, double check input")

    # Calculate differences across timepoints
    data_diff = np.diff(data, axis=0)

    # Calculate DVARS as the root mean square of the differences
    dvars = np.hstack((np.nan, np.sqrt(np.mean(data_diff**2, axis=1))))

    return dvars


def calc_FD(R, FD_type=1):
    '''
    This function calculates framewise displacement (Power et al. 2012 Neuroimage)
    The columns 3-6 (angular displacement) is assumed to be already converted to mm before passed in this function
    '''

    dR = np.diff(R, axis=0)  # First-order derivative
    ddR = np.diff(dR, axis=0)  # Second-order derivative
    if FD_type == 1:
        # L1-norm - sum of absolute values of first-order derivatives
        FD = np.sum(np.absolute(dR), axis=1)
        meanFD = np.mean(FD)
        FD = np.hstack(
            (
                np.zeros(
                    1,
                ),
                FD,
            )
        )  # Pad zeros to make it the same length as the original data
    elif FD_type == 2:
        # L2-norm - sum of absolute values of second-order derivatives
        FD = np.sum(np.absolute(ddR), axis=1)
        meanFD = np.mean(FD)
        FD = np.hstack(
            (
                np.zeros(
                    2,
                ),
                FD,
            )
        )  # Pad zeros to make it the same length as the original data
    return FD, meanFD


def make_friston_regressors(R, hd_mm):
    """
    This function takes a matrix `MR` of 6 degrees of freedom (DOF) movement correction
    parameters and calculates the corresponding 24 Friston regressors.

    Parameters:
    -----------
    MR : numpy array of shape (r, c)
        A matrix where r is the number of time points and c are the 6 DOF movement regressors.
        If the number of columns is more than 6, only the first 6 columns are considered.

    hd_mm : float, optional
        The head radius in mm. Default is 50 mm.

    Returns:
    --------
    FR : numpy array of shape (r, 24)
        A matrix containing 24 Friston regressors.
    """
    MR = R[:, :6]
    MR[:, 3:] = MR[:, 3:] * np.pi * hd_mm / 180
    # Calculate the first part of the Friston regressors (MR and MR^2)
    FR = np.hstack([MR, MR**2])

    # Create a dummy array for the temporal derivatives (lagged version of FR)
    dummy = np.zeros_like(FR)
    dummy[1:, :] = FR[:-1, :]  # shift FR by one time step
    dummy[0, :] = 0  # set the first row to 0

    # Concatenate the original FR and the lagged version
    FR = np.hstack([FR, dummy])
    return FR


def calc_power_spectrum(data):
    assert (
        len(data.shape) == 1 or data.shape[1] == 1
    ), "Input has to be a column vector or 1D array"
    N = data.shape[0]
    yf = scipy.fftpack.fft(data)
    power_spectrum = 2.0 / N * np.abs(yf[0 : N // 2])
    return power_spectrum


def detrend_manual(data, keepframe):
    '''
    Remove linear trends
    '''
    detrended_data = data.copy()
    time_points = np.where(keepframe)[0]
    time_points_all = np.array(range(keepframe.shape[0]))

    # Create the design matrix for linear regression (constant + linear term)
    X = np.vstack(
        [time_points, np.ones(len(time_points))]
    ).T  # Shape (len(keepframe), 2)
    Xall = np.vstack([time_points_all, np.ones(len(time_points_all))]).T

    # Perform the linear regression for all columns at once using least squares
    # Y is the data[keepframe, :] with shape (len(keepframe), n)
    Y = data[keepframe, :]

    # Compute the least squares solution to find the slope and intercept for each column
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)  # beta shape is (2, n)

    # Calculate the trend for each column using the coefficients
    trend = Xall @ beta  # Shape (len(keepframe), n)

    # Subtract the trend from the data at the all indices
    detrended_data -= trend

    return detrended_data


def calc_etasquared(a, b):
    """
    Calculate eta squared based on Cohen 2008 Neuroimage.

    Parameters:
    a : np.ndarray
        First input array.
    b : np.ndarray
        Second input array.

    Returns:
    etasquared : np.ndarray
        Array of eta squared values.
    """

    # Ensure inputs are at least 2D
    if a.ndim == 1:
        a = a[:, np.newaxis]
    if b.ndim == 1:
        b = b[:, np.newaxis]

    assert a.shape[0] == b.shape[0], "input size mismatch"

    cols_a = a.shape[1]
    cols_b = b.shape[1]
    etasquared = np.full((cols_b, cols_a), np.nan)

    for ia in range(cols_a):
        for ib in range(cols_b):
            aa = a[:, ia]
            bb = b[:, ib]

            m = (aa + bb) / 2
            Mbar = np.nanmean(m)

            SSwithin = np.nansum((aa - m) ** 2 + (bb - m) ** 2)
            SStotal = np.nansum((aa - Mbar) ** 2 + (bb - Mbar) ** 2)
            etasquared[ib, ia] = 1 - SSwithin / SStotal

    return etasquared

if __name__ == "__main__":
    main()
