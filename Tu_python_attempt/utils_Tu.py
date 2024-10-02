# Wanted to use the make_mask() function to replace fslmaths but the Gaussian filter seems to be problematic

import nibabel as nib
import numpy as np
import scipy.ndimage as ndimage

def threshold_image(img_data, lower_thr, upper_thr):
    return np.where((img_data >= lower_thr) & (img_data <= upper_thr), img_data, 0)

def binarize_image(img_data):
    return np.where(img_data > 0, 1, 0)

def gaussian_smooth(img_data, sigma_voxels):
    return ndimage.gaussian_filter(img_data, sigma=sigma_voxels)

def erode_image(img_data, sigma_voxels):
    # Create a 3D grid for the Gaussian kernel (size chosen based on sigma)
    kernel_size = int(6 * sigma_voxels)  # Typically, 6 * sigma is enough to cover the Gaussian
    grid = np.ogrid[
        -kernel_size // 2: kernel_size // 2 + 1,
        -kernel_size // 2: kernel_size // 2 + 1,
        -kernel_size // 2: kernel_size // 2 + 1
    ]

    # Create the Gaussian kernel in 3D (adjust for 2D if needed)
    gaussian_kernel = np.exp(-(grid[0]**2 + grid[1]**2 + grid[2]**2) / (2.0 * sigma_voxels**2))
    return ndimage.binary_erosion(img_data, structure=gaussian_kernel)
    #return ndimage.binary_erosion(img_data, structure=np.ones((sigma_voxels, sigma_voxels, sigma_voxels)))

def make_masks_scipy(segmentation, wm_mask_out, vent_mask_out, roi_res, **kwargs):
    # Set default values for thresholds
    defaults = dict(
        wm_lt_R=2950, wm_ut_R=3050, wm_lt_L=3950, wm_ut_L=4050,
        vent_lt_R=43, vent_ut_R=43, vent_lt_L=4, vent_ut_L=4,
        roi_res=2
    )
    
    # Update the defaults with any user-provided overrides (from kwargs)
    defaults.update(kwargs)

    # Load the segmentation image
    segmentation_img = nib.load(segmentation)
    segmentation_data = segmentation_img.get_fdata()

    # Get voxel dimensions for smoothing
    voxel_dims = segmentation_img.header.get_zooms()
    # Convert millimeter-based smoothing (roi_res) to voxel space
    sigma_voxels = np.int8(roi_res / voxel_dims[0])  # Convert mm to voxels

    # Right White Matter Mask
    wm_mask_R = threshold_image(segmentation_data, defaults['wm_lt_R'], defaults['wm_ut_R'])

    # Left White Matter Mask
    wm_mask_L = threshold_image(segmentation_data, defaults['wm_lt_L'], defaults['wm_ut_L'])

    # Combine and binarize the white matter masks
    wm_mask = binarize_image(wm_mask_R + wm_mask_L)

    # Smooth and erode the white matter mask
    #wm_mask_smooth = gaussian_smooth(wm_mask, sigma_voxels)
    wm_mask_eroded = erode_image(wm_mask, sigma_voxels=sigma_voxels).astype(np.int8)
    print(np.sum(wm_mask_eroded))
    # Save the final white matter mask
    wm_mask_out_img = nib.Nifti1Image(wm_mask_eroded, segmentation_img.affine)
    nib.save(wm_mask_out_img, wm_mask_out)

    # Right Ventricular Mask
    vent_mask_R = threshold_image(segmentation_data, defaults['vent_lt_R'], defaults['vent_ut_R'])

    # Left Ventricular Mask
    vent_mask_L = threshold_image(segmentation_data, defaults['vent_lt_L'], defaults['vent_ut_L'])

    # Combine and binarize the ventricular masks
    vent_mask = binarize_image(vent_mask_R + vent_mask_L)

    # Smooth and erode the ventricular mask
    #vent_mask_smooth = gaussian_smooth(vent_mask, sigma_voxels)
    vent_mask_eroded = erode_image(vent_mask, sigma_voxels=sigma_voxels).astype(np.int8)
    
    # Save the final ventricular mask
    vent_mask_out_img = nib.Nifti1Image(vent_mask_eroded, segmentation_img.affine)
    nib.save(vent_mask_out_img, vent_mask_out)

