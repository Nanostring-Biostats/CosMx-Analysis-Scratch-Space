# ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.

#!/usr/bin/env python
# coding: utf-8
"""
Application: make_composite.py
Author: Vikram Kohli, PSS
Updater: Jacob Hanimann
Version: 1.2.4

Description:
This script creates composite images from layered morphology 2D images.
It extracts layers from TIFF files, converts them to 8-bit (and autocontrasts them),
and generates screen composite images.
The output folders include the chosen file format as a suffix so that processing can be
stopped and later resumed; existing (and valid) files are skipped.

User inputs:
    clipping    - Histogram clipping percentage (adjusts image contrast).
    user_format - Output image format (allowed: jpg, png, tif).

Output:
    raw_<format>                   - Extracted images from the TIFF files.
    8bit_<format>                  - 8-bit converted images.
    8bit_autocontrast_<format>     - Autocontrasted 8-bit images.
    composite_<format>             - Composite images.
    composite_autocontrast_<format> - Autocontrasted composite images.
"""

# Loading required libraries
from PIL import Image, ImageSequence, ImageOps, ImageChops
import numpy as np
import glob
import os
import sys
import re
import shutil
from time import perf_counter
import concurrent.futures

# Allowed image formats and constants
IMAGE_FORMATS = ['jpg', 'png', 'tif']
FOLDER_NAMES = {
    "raw": "raw",
    "8bit": "8bit",
    "8bit_autocontrast": "8bit_autocontrast",
    "composite": "composite",
    "composite_autocontrast": "composite_autocontrast"
}
COLORS = ['cyan', 'red','yellow', 'blue', 'magenta'] # Composite color scheme (channel order)
# COLORS = ['green', 'yellow', 'gray', 'red', 'blue'] 
PATTERN = r"([\w]+_[\w]+_[\w]+_[\w]+_[\w]+_[\w]+_[F]+[\d]*)"
REGEX = re.compile(PATTERN, flags=re.IGNORECASE)
COMPRESS_VALUE = 3  # Lossless file compression value

class ArgumentException(Exception):
    """Exception class for handling argument errors."""
    pass

def check_arguments():
    """Prompt the user for clipping percentage and output format."""
    clipping_input = input('Please specify a clipping percentage (integer or float): ').strip()
    try:
        clipping = float(clipping_input)
    except ValueError:
        print('Error: Clipping value must be a number. Exiting.')
        sys.exit(1)
    
    user_format = input('Please specify one of the allowed image format types [jpg, png, tif]: ').lower().strip()
    if user_format not in IMAGE_FORMATS:
        raise ArgumentException('Error: Allowed formats include: jpg, png, tif. Exiting.')
    
    return clipping, user_format

def create_folders(base_dir, user_format):
    """
    Create the necessary output folders with the file format as a suffix.
    For example, if user_format is 'jpg', folders will be 'raw_jpg', '8bit_jpg', etc.
    """
    updated_folder_names = { key: f"{folder}_{user_format}" for key, folder in FOLDER_NAMES.items() }
    paths = {key: os.path.join(base_dir, folder) for key, folder in updated_folder_names.items()}
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
    print("Folders are ready:", ", ".join(paths.values()))
    return paths

def check_tif_integrity(filepath):
    """Verify the integrity of a TIFF file by iterating through its frames."""
    try:
        with Image.open(filepath) as img:
            for _ in ImageSequence.Iterator(img):
                pass
        return True
    except Exception as e:
        print(f"WARNING: Corrupted TIFF detected: {filepath}: {e}")
        return False

def check_image_integrity(filepath):
    """Verify the integrity of an image file by attempting to open it."""
    try:
        with Image.open(filepath) as img:
            img.verify()
        return True
    except Exception as e:
        print(f"WARNING: Corrupted image detected: {filepath}: {e}")
        return False

def file_exists_and_valid(path):
    """Return True if file exists and is not corrupted."""
    return os.path.exists(path) and check_image_integrity(path)

def process_tiff_file(image_file, raw_folder, user_format):
    """
    Extract image layers from a TIFF file if its filename matches the expected pattern.
    For each channel, skip processing if the expected output file already exists and is valid.
    """
    base_name = os.path.splitext(os.path.basename(image_file))[0]
    if not REGEX.match(base_name):
        return 0  # Skip file if filename does not match the pattern
    
    if not check_tif_integrity(image_file):
        print(f"Skipping corrupted TIFF file: {image_file}")
        return 1  # Count as corrupted

    fov_num = base_name.split('_')[-1]
    try:
        with Image.open(image_file) as image:
            for channel, layer in enumerate(ImageSequence.Iterator(image)):
                output_filename = f"{fov_num}_ch{channel}_raw.{user_format}"
                output_path = os.path.join(raw_folder, output_filename)
                if file_exists_and_valid(output_path):
                    continue  # Skip if already processed and valid
                if user_format == 'jpg':
                    layer.point(lambda value: value * (1. / 256)).convert('L').save(output_path, compress_type=COMPRESS_VALUE)
                else:
                    layer.save(output_path, compress_type=COMPRESS_VALUE)
    except Exception as e:
        print(f"Error processing file {image_file}: {e}")
    return 0

def force_8bit(image_obj):
    """Convert an image to 8-bit using numpy for fast vectorized operations."""
    array = np.array(image_obj)
    max_val = array.max()
    if max_val == 0:
        reduced_bit = array * 255.0
    else:
        reduced_bit = (array / max_val) * 255.0
    return Image.fromarray(reduced_bit.astype('uint8'))

def convert_8bit_worker(raw_image_file, user_format, clipping, paths):
    """
    Convert a raw extracted image to 8-bit and generate an autocontrasted version.
    Skip conversion if output files exist and are valid.
    """
    try:
        base_name = os.path.basename(raw_image_file)
        name_8bit = base_name.replace('_raw', '_8bit')
        name_autocontrast = base_name.replace('_raw', '_8bit_autocontrast')
        bit_path = os.path.join(paths['8bit'], name_8bit)
        auto_path = os.path.join(paths['8bit_autocontrast'], name_autocontrast)
        
        if file_exists_and_valid(bit_path) and file_exists_and_valid(auto_path):
            return
        
        if user_format == 'jpg':
            shutil.copy2(raw_image_file, bit_path)
            with Image.open(raw_image_file) as image:
                autocontrasted = ImageOps.autocontrast(image, cutoff=clipping)
                autocontrasted.save(auto_path, compress_type=COMPRESS_VALUE)
        else:
            with Image.open(raw_image_file) as image:
                img8 = force_8bit(image).convert('L')
                img8.save(bit_path, compress_type=COMPRESS_VALUE)
                autocontrasted = ImageOps.autocontrast(img8, cutoff=clipping)
                autocontrasted.save(auto_path, compress_type=COMPRESS_VALUE)
    except Exception as e:
        print(f"Error converting {raw_image_file}: {e}")

def write_composite(image_files, colors):
    """
    Create a screen composite image from a list of image files using colorization.
    Colors are applied in reverse order.
    """
    try:
        out = Image.open(image_files[-1])
        out = ImageOps.colorize(out, black='black', white=colors[-1])
        for image_file, color in zip(image_files[::-1][1:], colors[::-1][1:]):
            img = Image.open(image_file)
            img = ImageOps.colorize(img, black='black', white=color)
            out = ImageChops.screen(out, img)
        return out
    except Exception as e:
        print(f"Error creating composite: {e}")
        return None

def extract_channel(filename):
    """
    Extract the channel number from a filename.
    Assumes filenames like "F123_ch0_8bit.jpg".
    Returns an integer channel number or -1 if not found.
    """
    match = re.search(r'_ch(\d+)_', filename)
    return int(match.group(1)) if match else -1

def composite_worker(fov, user_format, paths, colors, clipping):
    """
    Generate composite images for a given field-of-view (fov) by processing the 8-bit images.
    Skip composite generation if output files already exist and are valid.
    """
    try:
        comp_filename = f"{fov}_composite.{user_format}"
        auto_comp_filename = f"{fov}_composite_autocontrast.{user_format}"
        comp_path = os.path.join(paths['composite'], comp_filename)
        auto_comp_path = os.path.join(paths['composite_autocontrast'], auto_comp_filename)
        
        if file_exists_and_valid(comp_path) and file_exists_and_valid(auto_comp_path):
            return
        
        search_pattern = os.path.join(paths['8bit'], f"{fov}_*.{user_format}")
        # Sort the image files based on the channel number extracted from the filename.
        image_files = sorted(glob.glob(search_pattern), key=extract_channel)
        if not image_files:
            return
        composite_img = write_composite(image_files, colors)
        if composite_img is None:
            return
        composite_img.save(comp_path, compress_type=COMPRESS_VALUE)
        auto_comp = ImageOps.autocontrast(composite_img, cutoff=clipping)
        auto_comp.save(auto_comp_path, compress_type=COMPRESS_VALUE)
    except Exception as e:
        print(f"Error in composite_worker for {fov}: {e}")

def main():
    # Clear the screen (optional)
    os.system('cls' if os.name == 'nt' else 'clear')
    
    base_dir = os.getcwd()
    print('********************************')
    print('   Composite script started     ')
    print('********************************\n')
    
    try:
        clipping, user_format = check_arguments()
    except ArgumentException as e:
        print(e)
        sys.exit(1)
    
    print(f'\nClipping value: {clipping}% and export format: {user_format}')
    print('\nPreparing necessary folders...')
    paths = create_folders(base_dir, user_format)
    
    start_time = perf_counter()
    print(f'\nProcessing TIFF files in {base_dir}...')
    
    tif_files = glob.glob('*.tif') + glob.glob('*.TIF')
    corrupted_total = 0
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_tiff_file, tif, paths['raw'], user_format)
                   for tif in tif_files]
        for future in concurrent.futures.as_completed(futures):
            corrupted_total += future.result()
    if corrupted_total > 0:
        print(f"\nTotal corrupted TIFF files skipped: {corrupted_total}")
    
    print('\nConverting raw files to 8-bit and creating autocontrasted images...')
    raw_files = glob.glob(os.path.join(paths['raw'], f"*_raw.{user_format}"))
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(convert_8bit_worker, raw_file, user_format, clipping, paths)
                   for raw_file in raw_files]
        concurrent.futures.wait(futures)
    
    print('\nCreating composite images...')
    bit_files = glob.glob(os.path.join(paths['8bit'], f"*.{user_format}"))
    fov_set = {os.path.basename(bf).split('_')[0] for bf in bit_files}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(composite_worker, fov, user_format, paths, COLORS, clipping)
                   for fov in fov_set]
        concurrent.futures.wait(futures)
    
    end_time = perf_counter()
    print(f'\nFinished. Total run time: {int(end_time - start_time)} sec')

if __name__ == '__main__':
    main()
