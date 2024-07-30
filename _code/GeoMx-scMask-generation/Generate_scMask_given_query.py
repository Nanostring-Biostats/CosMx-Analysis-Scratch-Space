"""
Copyright ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.
Author: Lidan Wu
Date: 2024-06-11
Version: v1.0.0
"""

import argparse
import sys
import os
import re
import subprocess
import logging
from datetime import date
import numpy as np
import cv2
from cellpose import io
from skimage import measure, morphology, exposure
from scipy.ndimage import label, find_objects


version_str = f"""
pipeline version: \t v1.0.0\n"""


# setup logging for full script
def main_logger_setup(logger_path:str):
    """ Set up logging of full pipeline"""
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[
                            logging.FileHandler(logger_path),
                            logging.StreamHandler(sys.stdout)
                        ])
    logger = logging.getLogger(__name__)
    logger.info(f"Writing log output of mask generation pipeline to {logger_path}")
    logger.info(version_str)
    return logger

    
# rearrange channel order
def rearrange_channels(image_path: str, channel_order: list) -> np.ndarray:
    """
    Rearrange the channels of a multi-channel image according to a user-defined order.

    Parameters:
    - image_path (str): Path to the input multi-channel image.
    - channel_order (list): User-defined order for channels.

    Returns:
    - np.ndarray: Image with channels rearranged.

    Raises:
    - ValueError: If the input image has fewer channels than required by channel_order.
    """
    # Read the multi-frame image from the file path using cellpose.io
    image = io.imread(image_path)
    
    if image is None:
        raise ValueError(f"Failed to read: {image_path}")
    if image.ndim != 3:
        raise ValueError(f"Not a multi-channel single-frame image: {image_path}")
    else:
        if np.argmin(image.shape) !=0:
            raise ValueError(f"Color must come first, smallest dimension found at Postion {np.argmin(image.shape)+1} in: {image_path}")
        
    # assume color at first position
    # Check if the image has the required number of channels
    if image.shape[0] < max(channel_order) + 1:
        raise ValueError(f"Input image has {image.shape[0]} channels, but channel_order requires at least {max(channel_order) + 1} channels.")

    # Reorder the channels according to the specified order
    reordered_image = image[channel_order, ...]

    return reordered_image
    
    
    
    
# get positive-stained object labels
def get_clean_object_labels(image: np.ndarray, method: str = 'triangle', thresh: int = 0, 
                            fill_holes: bool =  False, kernel_size: int = 5,
                            min_size: int = 50) -> np.ndarray:
    """
    Convert a gray-scale image to a segmented label image using auto-threshold and clean up small objects.

    Parameters:
    - image (np.ndarray): Input 2D gray-scale image.
    - method (str): Thresholding method ('otsu' or 'triangle').
    - thresh (int): Minimal value in 8-bit images for positive regions 
    - fill_holes (bool): Flag to fill holes of binary mask before further filtering on size, output round borders if True  
    - kernel_size (int): Kernel size to do morphology closing on binary masks when fill_holes = True
    - min_size (int): Minimum size of objects to keep.


    Returns:
    - np.ndarray: label image for positive-stain region with small objects removed.

    Raises:
    - ValueError: If an invalid thresholding method is specified.
    """
    # Check if the input image is 2D
    if image.ndim != 2:
        raise ValueError("Input image must be a 2D gray-scale image.")

    # auto-threshold with Otsu or Triangle only works for 8-bit image 
    # Convert images to 16-bit and then rescaled to 8-bit
    image = exposure.rescale_intensity(image.astype('uint16'))
    image = (image/65535 * 255).astype('uint8')   

    # Apply auto-threshold, triangle is more permissive 
    if method == 'otsu':
        _, binary_mask = cv2.threshold(image, thresh, 255, cv2.THRESH_OTSU)
    elif method == 'triangle':
        _, binary_mask = cv2.threshold(image, thresh, 255, cv2.THRESH_TRIANGLE)
    else:
        raise ValueError("Invalid method. Choose 'otsu' or 'triangle'.")

    # Convert binary mask to boolean
    binary_mask = binary_mask > 0

    # fill holes using morphology closing (dilation and then erosion)
    # alternative methods: cv2.findContours, cv2.floodFill
    if fill_holes:
        # Define the kernel size, the size depends on the size of holes to be filled
        # kernel size must be a positive odd number, plus one if even
        kernel_size = max([3, int(kernel_size/2)*2+1])
        kernel = np.ones((kernel_size, kernel_size), np.uint8)
        # Perform morphological closing
        binary_mask = cv2.morphologyEx(binary_mask, cv2.MORPH_CLOSE, kernel)

    # Remove small objects
    cleaned_mask = morphology.remove_small_objects(binary_mask, min_size=min_size)
    
    # Create a labeled version of mask to identify connected components
    labeled_mask, _ = label(cleaned_mask)

    return labeled_mask
    
    

# connectivity operation of two labels 
def group_cells_by_intersection(L1: np.ndarray, L2: np.ndarray, area_cutoff: int = 50) -> (np.ndarray, np.ndarray):
    """
    Find the touching and non-touching cells between two label images based on intersection area.

    Parameters:
    - L1 (np.ndarray): First set of cell labels.
    - L2 (np.ndarray): Second set of cell labels.
    - area_cutoff (int): The minimum intersection area to consider cells as touching.

    Returns:
    - (np.ndarray, np.ndarray): Tuple containing two binary masks:
        1. Combined region of cells with intersection area larger than cutoff.
        2. Cells in L1 with intersection area smaller than cutoff.
    """
    # Ensure L1 and L2 have the same shape
    if L1.shape != L2.shape:
        raise ValueError("L1 and L2 must have the same shape.")

    # Initialize the masks
    combined_touching_mask = np.zeros_like(L1, dtype=bool)
    non_touching_mask = np.zeros_like(L1, dtype=bool)

    # Find properties of labeled regions in L1 and L2
    regions_L1 = measure.regionprops(L1)
    regions_L2 = measure.regionprops(L2)

    # Create a set to keep track of touching labels in L1
    touching_labels_L1 = set()

    # Iterate through each region in L1
    for region_L1 in regions_L1:
        label_L1 = region_L1.label
        coords_L1 = region_L1.coords

        # Create a binary mask for the current region in L1
        mask_L1 = np.zeros_like(L1, dtype=bool)
        mask_L1[coords_L1[:, 0], coords_L1[:, 1]] = True

        # Initialize a flag to check if the region touches any region in L2
        is_touching = False

        # Check for intersection area with each region in L2
        for region_L2 in regions_L2:
            coords_L2 = region_L2.coords

            # Create a binary mask for the current region in L2
            mask_L2 = np.zeros_like(L2, dtype=bool)
            mask_L2[coords_L2[:, 0], coords_L2[:, 1]] = True

            # Compute the intersection area
            intersection_area = np.sum(mask_L1 & mask_L2)

            if intersection_area > area_cutoff:
                is_touching = True
                # Add the touching region from L2 to the combined mask
                combined_touching_mask[coords_L2[:, 0], coords_L2[:, 1]] = True

        if is_touching:
            # Mark the region in L1 as touching
            combined_touching_mask[coords_L1[:, 0], coords_L1[:, 1]] = True
            touching_labels_L1.add(label_L1)
        else:
            # Mark the region in L1 as non-touching
            non_touching_mask[coords_L1[:, 0], coords_L1[:, 1]] = True
    
    # make sure no overlapping pixel between 2 masks
    # if overlapping, remove from both masks 
    flagBad = combined_touching_mask & non_touching_mask
    combined_touching_mask = np.where(flagBad, False, combined_touching_mask)
    non_touching_mask = np.where(flagBad, False, non_touching_mask)

    return combined_touching_mask, non_touching_mask


    
def get_arg_parser():
    """ Parses command line arguments for main function"""
    parser = argparse.ArgumentParser(description="Generate 2 sets of binary masks for positive-vs-negative stained single cells from multi-channel images")
    
    # settings for IO
    input_img_args = parser.add_argument_group("Input Image Arguments")
    input_img_args.add_argument("--in_dir", default=".", type=str, help="absolute path to folder containing multi-channel images to evaluate")
    input_img_args.add_argument("--out_dir", default=None, type=str, help="output parent folder for intermediate results and final masks. (default to use input image folder)")
    input_img_args.add_argument("--cyto_chan", default=1,type=int, help="input channel index for cytoplasm or membrane stain of cells (starting from 1)")
    input_img_args.add_argument("--query_chan", default=2,type=int, help="input channel index for the stain in query to define positive vs. negative staining (starting from 1)")
    input_img_args.add_argument("--nuc_chan", default=0,type=int, help="input channel index for nuclear stain of cells (starting from 1); default to 0 to exclude nuclear stain from cell segmentation")
    
    # settings for cell segmentation
    model_args = parser.add_argument_group("Cell Segmentation Model Arguments")
    model_args.add_argument("--cell_seg_model", default="cyto3",type=str, help="cell segmentation model in use, ok to pass a file path of custom model")
    model_args.add_argument("--cell_diameter", default=30.,type=float, help="median cell diameter of input images, in pixel unit")
    model_args.add_argument("--min_cell_area", default=15,type=int, help="minimal area of a valid cell, in squared pixel unit; can turn off with -1")
    
    # settings for defining positive-stained cells 
    pos_label_args = parser.add_argument_group("Positive-stained cell Arguments")
    pos_label_args.add_argument("--thresh_method", default="triangle",type=str, help="auto-threshold method to define positive-stained object in query channel, use either `triangle` or `otsu`")
    pos_label_args.add_argument("--fill_holes", action="store_true", help="fill holes in positive-stained object before area filtering; if True, use 0.1x of input cell diameter as kernel size and output roundish borders")
    pos_label_args.add_argument("--min_positive_area", default=500.,type=float, help="minimal area of a positive-stained object in query channel, in squared pixel unit; recommend to be 0.5x of expected cell area")
    pos_label_args.add_argument("--min_intersect_area", default=50.,type=float, help="minimal intersection area between a positive-stained cell with nearest positive-stained object in query channel, in squared pixel unit; recommend to be 0.1x of expected cell area")
    pos_label_args.add_argument("--clean_export", action="store_true", help="export masks into file only when there are cells selected in either case")
    
    return parser

def main():
    """ Run pipeline from command line"""
    
    # get input argument
    args = get_arg_parser().parse_args()
	
    # set up 
    inDir =  r'%s' % args.in_dir
    if args.out_dir is None:
        outDir = inDir
    else:
        outDir = r'%s' % args.out_dir
    pretrained_model = r'%s' % args.cell_seg_model
    os.makedirs(outDir, exist_ok=True)
    
    log_file = date.today().strftime("%Y%m%d")
    log_file = os.path.join(outDir, f'{log_file}_pipeline_run_log.txt')
    logger = main_logger_setup(logger_path = log_file)

    # cell segmentation command 
    command = 'python -m cellpose --pretrained_model "' + pretrained_model + '" --save_tif --use_gpu --verbose --no_npy --diameter ' + str(args.cell_diameter) + ' --min_size ' + str(args.min_cell_area)
    
    if args.nuc_chan <1:
        # no nuclear for segmentation
        cpDir = os.path.join(outDir, '2chan_res')
        chan_order = [args.cyto_chan - 1, args.query_chan - 1]
        command = command  +' --chan 1' 
    else:
        # include nuclear for segmentation 
        cpDir = os.path.join(outDir, '3chan_res')
        chan_order = [args.cyto_chan - 1, args.query_chan - 1, args.nuc_chan -1]
        command = command  + ' --chan 1 --chan2 3' 

    # rearrange the input images to new channel order of 1-cyto, 2-query, 3-nuc (optional)
    image_names = io.get_image_files(inDir, mask_filter="_masks")
    nimg = len(image_names)	
    if(nimg <1):
        raise ValueError(f"No images found in: {inDir}")
    
    logger.info("\n")
    logger.info(f"Start channel reorder with new order = [{', '.join(map(str, chan_order))}] , output folder named `{os.path.basename(cpDir)}`\n")
    os.makedirs(cpDir, exist_ok=True)
    for im in image_names:
        newFile = os.path.join(cpDir, os.path.splitext(os.path.basename(im))[0]+'.tiff')
        if os.path.isfile(newFile):
           logger.info(f"File exists, skip reorder for: {os.path.basename(newFile)}")
        else:
            newImg = rearrange_channels(image_path = im, channel_order = chan_order) 
            # save as tiff, channel first  
            io.imsave(newFile, newImg)
    newImg_names = io.get_image_files(cpDir, mask_filter="_masks")
    
    # run cell segmentation, overwriting if labels exist
    labelDir = os.path.join(cpDir, 'cellLabels')
    os.makedirs(labelDir, exist_ok=True)
    command = command  +' --dir "' + cpDir + '" --savedir "' + labelDir + '"'
    # check if labels exist, skip if all presents, overwrite all if partial presents
    lab_names = [os.path.join(labelDir, os.path.splitext(os.path.basename(im))[0]+'_cp_masks.tif') for im in newImg_names]
    lab_exist = [os.path.exists(labImg) for labImg in lab_names]
    
    logger.info("\n")
    if all(lab_exist):
        logger.info(f"All images have existing cell labels, skip cell segmentation.\n")
    else:
        logger.info(f"Start cell segmentation, overwriting {sum(lab_exist)} existing labels.\n")
        # copy outputs to separate log file 
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(result.stdout)
        if result.stderr:
            logger.error(result.stderr)
            
    
    # loop thourgh each new images and get 2 sets of masks 
    posDir = os.path.join(outDir, 'postive_ROI_masks')
    negDir = os.path.join(outDir, 'negative_cell_masks')
    os.makedirs(posDir, exist_ok=True)
    os.makedirs(negDir, exist_ok=True)
    
    # kernel size for fill holes in query stain, use 0.1x of input cell diameter  
    kernel_size = round(args.cell_diameter * 0.1)
    # minimal intensity for postiive query stain in 8-bit scale
    min_query = 0
    min_positive_area = round(args.min_positive_area)
    min_intersect_area = round(args.min_intersect_area)

    logger.info("\n")
    logger.info("Start segmentation of positive-stained objects in query channel and do conectivity analysis")
    logger.info(
                ">>>> query-stained objects: thresh_method = %s, fill_holes = %s, min_positive_area = %d"
                % (args.thresh_method, args.fill_holes, min_positive_area))
    logger.info(
                ">>>> positive-stained cells: min_intersect_area = %d between a cell and overlapping positive-stained object"
                % (min_intersect_area))
    if args.clean_export:
        logger.info(">>>> No export on empty masks")
    
    for n in range(len(newImg_names)):
        # get cell label, skip if no segmented cells
        labImg = io.imread(lab_names[n])
        if labImg.max() <1:
            logger.info(f"No cell identified, skip: {os.path.basename(newImg_names[n])}")
            continue
        else:
            # get object labels for query stain 
            img = io.imread(newImg_names[n])
            queryLab = get_clean_object_labels(image = img[1, ...], # query stain at 2nd channel 
                                               method= args.thresh_method, thresh = min_query, 
                                               fill_holes = args.fill_holes, kernel_size = kernel_size, 
                                               min_size= min_positive_area)
            # check connectivity of cell labels with query-stained objects
            Mask1, Mask2 = group_cells_by_intersection(L1 = labImg, L2 = queryLab, area_cutoff = min_intersect_area) 
            # convert bool matrix to 8-bit (0-255) mask for export
            Mask1 = (Mask1 * 255).astype(np.uint8)
            Mask2 = (Mask2 * 255).astype(np.uint8)
            
            if not args.clean_export:
                io.imsave(os.path.join(posDir, os.path.splitext(os.path.basename(newImg_names[n]))[0]+'_pos_1.tif'), Mask1)
                io.imsave(os.path.join(negDir, os.path.splitext(os.path.basename(newImg_names[n]))[0]+'_neg_2.tif'), Mask2)
            else:
                # only export when mask contain positive pixel 
                if Mask1.max() >0:
                    io.imsave(os.path.join(posDir, os.path.splitext(os.path.basename(newImg_names[n]))[0]+'_pos_1.tif'), Mask1)
                if Mask2.max() >0:
                    io.imsave(os.path.join(negDir, os.path.splitext(os.path.basename(newImg_names[n]))[0]+'_neg_2.tif'), Mask2)
    
    logger.info("\n")
    logger.info(f"Process completed, output at: {outDir}\n")
      
	
if __name__ == '__main__':
    sys.exit(main())	
    
    
    
