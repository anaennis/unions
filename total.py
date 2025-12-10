#!/usr/bin/env python3
import argparse
import subprocess
import time
import re
import numpy as np
import os
import shutil
import logging
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter
import pandas as pd
import matplotlib.pyplot as plt

# Import the CONFIGuration file
import CONFIG

# Set up logging configuration
logging.basicConfig(filename=CONFIG.LOG_FILENAME, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s: %(message)s')


# Helper functions

def setup_directories(folder_name):
    """Creates the main data folder and the plots directory if they don't exist."""
    try:
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
            logging.info(f"Created data folder: {folder_name}")
        
        # Get data directory for plots
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, CONFIG.DATA_DIR)
        plot_dir = os.path.join(data_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
            logging.info(f"Created plots folder: {plot_dir}")
            
    except Exception as e:
        logging.error(f"Failed to create directories: {e}")
        raise



def in_UNIONS(RA_input, DEC_input, listing_file):
    """
    Checks if a galaxy's coordinates fall within an existing CFIS/UNIONS tile.

    Parameters
    ----------
    ra_input : float
        Right Ascension of the target (degrees).
    dec_input : float
        Declination of the target (degrees).
    listing_file : str
        Filename containing the list of UNIONS tiles.

    Returns
    -------
    tuple
        (check_if_tile: bool, xxx_cen: int, yyy_cen: int)
    """

    # Load tile list

    try:
        tiles = np.loadtxt(listing_file, dtype=str)
    except FileNotFoundError:
        logging.error(f"Tile listing file not found: {listing_file}")
        return False, 0, 0
    

    # Calculate tile centers and distance
    xxx = np.empty(len(tiles), dtype=int)
    yyy = np.empty(len(tiles), dtype=int)
    distance = np.empty(len(tiles), dtype=float)
    DEC_cen = np.empty(len(tiles),dtype=float)
    RA_cen = np.empty(len(tiles),dtype=float)


    for ii in range(len(tiles)):
        xxx[ii], yyy[ii] = np.array(tiles[ii].split('.'), dtype=int)
        
        DEC_cen[ii] = (yyy[ii] / 2) - 90
        RA_cen[ii] = xxx[ii] / (2 * np.cos(DEC_cen[ii] * np.pi / 180))

        distance[ii] = np.sqrt(np.power((DEC_input-DEC_cen[ii]),2)+np.power((RA_input-RA_cen[ii]),2))

    ch = np.where(distance==min(distance))
        
    # Use [0] since np.where returns an array of indices
    xxx_cen = xxx[ch][0]
    yyy_cen = yyy[ch][0]

    # Use [0] since np.where returns an array of indices

    check_if_tile = np.min(distance) < 0.5    
       
    return check_if_tile,xxx_cen,yyy_cen

# -- Main Pipeline Stages -- 

def download_data(ra, dec, xxx_input, yyy_input, folder):
    """
    Downloads CFIS tiles using vcp (VOSpace client).
    
    Parameters
    ----------
    ra : float
        Right Ascension (deg).
    dec : float
        Declination (deg).
    xxx_input : str
        Tile RA index (formatted '000').
    yyy_input : str
        Tile DEC index (formatted '000').
    folder : str
        Local directory path for the downloads.
        
    Raises
    ------
    subprocess.CalledProcessError
        If VOSpace download or local file check fails.
    """
    logging.info(f"Attempting download for tile CFIS.{xxx_input}.{yyy_input}")
    
    file_u = f'CFIS.{xxx_input}.{yyy_input}.u.fits'
    file_r = f'CFIS.{xxx_input}.{yyy_input}.r.fits'
    file_g_whigs = f'calexp-CFIS_{xxx_input}_{yyy_input}.fits' 
    file_g_final = f'CFIS.{xxx_input}.{yyy_input}.g.fits'


    try:
        import subprocess

        # Define file locations in VOSpace
        vospace_u = f'vos:cfis/tiles_DR5/{file_u}'
        vospace_r = f'vos:cfis/tiles_DR5/{file_r}'
        vospace_g = f'vos:cfis/whigs/stack_images_CFIS_scheme/{file_g_whigs}'
        # Download commands
        subprocess.run(['vcp', vospace_u, folder], check=True, capture_output=True)
        subprocess.run(['vcp', vospace_r, folder], check=True, capture_output=True)
        subprocess.run(['vcp', vospace_g, folder], check=True, capture_output=True)
        
        logging.info(f'Downloaded {file_u}, {file_r}, and {file_g_whigs}')

        # Standardize G-band FITS 
        old_g_path = os.path.join(folder, file_g_whigs)
        new_g_path = os.path.join(folder, file_g_final)
        r_path = os.path.join(folder, file_r) # Path to r-band file for header template
        
        with fits.open(old_g_path) as hdg, fits.open(r_path) as hdr_template:
            # Copy data from HDU 1 of the WHIGS file to a standard FITS structure
            hdr_template[0].data = hdg[1].data
            hdr_template[0].header['FILTER'] = 'G'
            hdr_template[0].header['GAIN'] = CONFIG.DEFAULT_GAIN
            hdr_template.writeto(new_g_path, overwrite=True)
        
        os.remove(old_g_path) # Clean up original WHIGS file
        logging.info(f"Standardized G-band file to {file_g_final}")

    except subprocess.CalledProcessError as e:
        logging.error(f'VOSpace command failed: {e.cmd}. Output: {e.stdout.decode()} {e.stderr.decode()}')
        if os.path.exists(folder):
            shutil.rmtree(folder)
            logging.info(f'Deleted incomplete folder: {folder}')
        raise
    
    except Exception as e:
        logging.error(f'An unexpected error occurred during download/standardization: {e}')
        if os.path.exists(folder):
            shutil.rmtree(folder)
            logging.info(f'Deleted incomplete folder: {folder}')
        raise


def process_and_detect(folder, gal_id, ra_sky, dec_sky, distance_mpc):
    """
    Performs galaxy light subtraction and runs source detection (using SExtractor).

    Parameters
    ----------
    folder : str
        Directory containing the FITS files.
    gal_id : int
        PGC ID of the galaxy.
    ra_sky : float
        Galaxy RA (deg).
    dec_sky : float
        Galaxy DEC (deg).
    distance_mpc : float
        Distance to the galaxy (Mpc).
    """
    logging.info(f"Starting detection phase for PGC {gal_id}")
    
    # Change to the folder containing the FITS files
    original_dir = os.getcwd()
    os.chdir(folder)

    pattern = r"CFIS\.(\d+)\.(\d+)\.([urgiz]).fits"
    fits_files = sorted([f for f in os.listdir('.') if f.endswith('.fits') and f.startswith('CFIS')])

    # Check if all processed images already exist
    expected_filters = ['r', 'u', 'g']
    processed_images = {}
    all_processed = True
    
    for fits_file in fits_files:
        match = re.match(pattern, fits_file)
        if not match:
            continue
        xxx, yyy, flt = match.groups()
        final_image_name = f'{flt}_fin_{xxx}.{yyy}.fits'
        processed_images[flt] = final_image_name
        
        if flt in expected_filters and not os.path.exists(final_image_name):
            all_processed = False
    
    if all_processed and all(processed_images.get(f) for f in expected_filters):
        logging.info(f"All processed images already exist. Skipping image processing for PGC {gal_id}.")
    else:
        # 1. Image Processing (Light Subtraction)
        logging.info(f"Processing images for PGC {gal_id}...")
        for fits_file in fits_files:
            match = re.match(pattern, fits_file)
            if not match:
                continue
            xxx, yyy, flt = match.groups()
            
            final_image_name = f'{flt}_fin_{xxx}.{yyy}.fits'
            processed_images[flt] = final_image_name

            if os.path.exists(final_image_name):
                logging.info(f"Processed image {final_image_name} already exists. Skipping filtering.")
                continue
            
            try:
                with fits.open(fits_file) as hdul:
                    data = hdul[0].data
                    
                    # --- Light Subtraction using CONFIG parameters ---
                    filtered_large = median_filter(data, size=CONFIG.FILTER_LARGE_SCALE)
                    residual_large = data - filtered_large
                    
                    filtered_small = median_filter(residual_large, size=CONFIG.FILTER_SMALL_SCALE)
                    final_residual = residual_large - filtered_small
                    
                    # Write the processed FITS file
                    hdul[0].data = final_residual
                    hdul[0].header['PGC'] = gal_id
                    hdul[0].header['LIGHTSUB'] = f'MEDIAN_{CONFIG.FILTER_LARGE_SCALE[0]}_{CONFIG.FILTER_SMALL_SCALE[0]}'
                    hdul[0].header['DISTANCE'] = distance_mpc
                    hdul.writeto(final_image_name, overwrite=True)
                    logging.info(f"Generated filtered image: {final_image_name}")
                    
            except Exception as e:
                logging.error(f"Error processing image {fits_file}: {e}")
                continue

    # 2. Source Detection 
    r_img = processed_images.get('r')
    u_img = processed_images.get('u')
    g_img = processed_images.get('g')
    
    if not (r_img and u_img and g_img):
        logging.error("Missing one or more processed FITS files (r, u, g). Cannot run DOLPHOT.")
        os.chdir(original_dir)
        return

    cat_name = f'full_catalog_{gal_id}_{ra_sky:.4f}_{dec_sky:.4f}.cat'
    
    # DOLPHOT requires image names without .fits extension
    detection_img_base = r_img.replace('.fits', '')
    g_img_base = g_img.replace('.fits', '')
    u_img_base = u_img.replace('.fits', '')
    
    # Use DOLPHOT for detection 
    dolphot_command = [
        'dolphot', cat_name, '-p' + CONFIG.DOLPHOT_PARAM_FILE, 
        'img1_file=' + detection_img_base,
        'img2_file=' + g_img_base, 
        'img3_file=' + u_img_base
    ]

    try:
        subprocess.run(dolphot_command, check=True, capture_output=True)
        logging.info(f'Successfully ran DOLPHOT. Catalog created: {cat_name}')
    except subprocess.CalledProcessError as e:
        logging.error(f"DOLPHOT failed: {e.cmd}. Error:\n{e.stderr.decode()}")
        
    os.chdir(original_dir) # Return to the original directory

def filter_globular_clusters(catalog_path, x_center_pix, y_center_pix, gal_id, distance_mpc):
    """
    Applies photometric and morphological cuts to select globular cluster candidates.
    
    Parameters
    ----------
    catalog_path : str
        Path to the source extracted catalog (CSV/TXT).
    x_center_pix, y_center_pix : float
        Pixel coordinates of the galaxy center.
    gal_id : int
        PGC ID for logging.
    distance_mpc : float
        Galaxy distance for calculating physical size limits.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame of GC candidates.
    """
    logging.info(f"Starting GC candidate filtering for PGC {gal_id}")
    
    # Define columns based on the default SExtractor output (needs to match your .param file)
    columns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000', 
               'MAG_R', 'MAGERR_R', 'MAG_U', 'MAGERR_U', 'MAG_G', 'MAGERR_G', 
               'CLASS_STAR']
    
    try:
        cat = pd.read_csv(catalog_path, delim_whitespace=True, comment='#', names=columns)
        
    except Exception as e:
        logging.error(f"Failed to read catalog at {catalog_path}: {e}")
        return pd.DataFrame()

    # Calculate distance from the galaxy center in pixels
    cat['RAD_PIX'] = np.sqrt(
        np.power(cat['X_IMAGE'] - x_center_pix, 2) + 
        np.power(cat['Y_IMAGE'] - y_center_pix, 2)
    )
    
    # Calculate colors
    cat['u_r'] = cat['MAG_U'] - cat['MAG_R']
    cat['g_r'] = cat['MAG_G'] - cat['MAG_R']

 # --- Core GC Selection Cuts using CONFIG parameters ---
    
    # 1. Morphological Cut (Stellarity index)
    cat_filtered = cat[cat['CLASS_STAR'] > CONFIG.GC_CLASS_STAR_THRESHOLD]
    
    # 2. Photometric Cuts 
    cat_filtered = cat_filtered[
        (cat_filtered['U_R'] > CONFIG.GC_U_R_MIN) & (cat_filtered['U_R'] < CONFIG.GC_U_R_MAX) &
        (cat_filtered['G_R'] > CONFIG.GC_G_R_MIN) & (cat_filtered['G_R'] < CONFIG.GC_G_R_MAX)
    ]
    
    # 3. Quality Cuts
    cat_filtered = cat_filtered[
        (cat_filtered['MAGERR_R'] < CONFIG.GC_MAG_ERR_MAX) & 
        (cat_filtered['MAG_R'] < CONFIG.GC_MAG_R_MAX) & (cat_filtered['MAG_R'] > CONFIG.GC_MAG_R_MIN)
    ]
    
    # 4. Radial Cut (in pixels)
    kpc_to_arcsec = (1 / (distance_mpc * 1e6)) * 206265 # arcsec per 1 pc
    radius_kpc_arcsec = CONFIG.GC_RADIAL_DISTANCE_KPC * 1000 * kpc_to_arcsec
    radius_kpc_pix = radius_kpc_arcsec / CONFIG.PIXEL_SCALE_CFIS
    
    cat_filtered = cat_filtered[cat_filtered['RAD_PIX'] < radius_kpc_pix]
    
    logging.info(f"Total objects detected: {len(cat)}. GC candidates found: {len(cat_filtered)}")
    
    return cat_filtered




    
def generate_plots(folder, gal_id, ra_sky, dec_sky, distance_mpc, final_catalog):
    """
    Generates key plots for the analysis: Color-Magnitude Diagram and radial distribution.
    
    Parameters
    ----------
    folder : str
        Directory containing the data.
    gal_id : int
        PGC ID.
    ra_sky, dec_sky : float
        Galaxy coordinates.
    distance_mpc : float
        Galaxy distance.
    final_catalog : pd.DataFrame
        The filtered DataFrame of GC candidates.
    """
    
    if final_catalog.empty:
        logging.warning(f"No GC candidates to plot for PGC {gal_id}.")
        return

    # Get absolute path to plots directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, CONFIG.DATA_DIR)
    plot_dir = os.path.join(data_dir, 'plots')
    
    # --- Plot 1: Color-Magnitude Diagram (CMD) ---
    fig_cmd, ax_cmd = plt.subplots(figsize=(7, 6))
    ax_cmd.scatter(final_catalog['U_R'], final_catalog['G_R'], s=15, c='darkorange', alpha=0.7)
    
    # Annotate the expected GC region boundaries from CONFIG
    ax_cmd.axvline(x=CONFIG.GC_U_R_MIN, color='blue', linestyle=':', alpha=0.5)
    ax_cmd.axvline(x=CONFIG.GC_U_R_MAX, color='blue', linestyle=':', alpha=0.5)
    ax_cmd.axhline(y=CONFIG.GC_G_R_MIN, color='red', linestyle=':', alpha=0.5)
    ax_cmd.axhline(y=CONFIG.GC_G_R_MAX, color='red', linestyle=':', alpha=0.5)


    ax_cmd.set_xlabel('u - r (mag)', fontsize=12)
    ax_cmd.set_ylabel('g - r (mag)', fontsize=12)
    ax_cmd.set_title(f'PGC {gal_id} CMD (GC Candidates)', fontsize=14)
    ax_cmd.grid(True, linestyle=':', alpha=0.6)
    fig_cmd.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'cmd_pgc{gal_id}.png'), dpi=150)
    plt.close(fig_cmd)
    logging.info(f"Generated CMD plot for PGC {gal_id}.")

    # --- Plot 2: Radial Distribution ---
    fig_rad, ax_rad = plt.subplots(figsize=(7, 6))
    
    # Convert radial distance from pixels to arcsec
    final_catalog['RAD_ARCSEC'] = final_catalog['RAD_PIX'] * CONFIG.PIXEL_SCALE_CFIS
    
    # Bin the radial distances
    ax_rad.hist(final_catalog['RAD_ARCSEC'], bins=20, histtype='step', 
                color='navy', linewidth=2, label=f'N={len(final_catalog)}')
    
    ax_rad.set_xlabel('Projected Radial Distance (arcsec)', fontsize=12)
    ax_rad.set_ylabel('Number of GCs', fontsize=12)
    ax_rad.set_title(f'PGC {gal_id} Radial Distribution', fontsize=14)
    ax_rad.legend()
    ax_rad.set_yscale('log')
    fig_rad.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'radial_pgc{gal_id}.png'), dpi=150)
    plt.close(fig_rad)
    logging.info(f"Generated radial distribution plot for PGC {gal_id}.")

 

def main():
    """
    Main function to parse arguments and execute the pipeline stages.
    """
    parser = argparse.ArgumentParser(
        description="Globular Cluster Candidate Search Pipeline. Requires VOSpace client (vcp) and SExtractor (sex) installed."
    )
    parser.add_argument(
        '--catalogue', 
        type=str, 
        required=True,
        help="Path to the input catalogue file. Must contain columns: RA (deg), Dec (deg), PGC_ID (int), Distance (Mpc)."
    )
    parser.add_argument(
        '--skip-download', 
        action='store_true',
        help="Skip the download step (use if files are already present locally)."
    )
    parser.add_argument(
        '--skip-detect', 
        action='store_true',
        help="Skip the detection and light subtraction step (use if catalogs are already generated)."
    )
    
    args = parser.parse_args()
    
    print(f"--- Starting Globular Cluster Pipeline (Logging to {CONFIG.LOG_FILENAME}) ---")
    
    start_time = time.time()

    # Load the catalogue data
    try:
        data = np.loadtxt(args.catalogue)
    except FileNotFoundError:
        logging.critical(f"Input catalogue file not found: {args.catalogue}")
        print("CRITICAL ERROR: Input catalogue not found.")
        return
    except Exception as e:
        logging.critical(f"Error loading catalogue: {e}")
        print("CRITICAL ERROR: Failed to load catalogue.")
        return

    num_galaxies = len(data)
    logging.info(f"Loaded {num_galaxies} galaxies from catalogue.")

    for i in range(num_galaxies):
        coords_ra, coords_dec, gal_id, dist_mpc = data[i]
        gal_id = int(gal_id)
        
        logging.info(f"--- Processing Galaxy {i+1}/{num_galaxies}: PGC {gal_id} at ({coords_ra:.4f}, {coords_dec:.4f}) ---")
        
        # 1. Check Tile Coverage
        is_in_tile = False
        xxx_input_int, yyy_input_int = 0, 0
        try:
            # Iterate through the configured tile listing files
            for listing_file in CONFIG.TILE_LISTINGS:
                is_in_tile, xxx_input_int, yyy_input_int = in_UNIONS(coords_ra, coords_dec, listing_file)
                if is_in_tile:
                    break # Found a tile, stop searching

            if not is_in_tile:
                logging.info(f"PGC {gal_id} not in any known CFIS tile. Skipping.")
                continue
                
            xxx_input = '{:03d}'.format(xxx_input_int)
            yyy_input = '{:03d}'.format(yyy_input_int)

            if yyy_input_int == 0 and yyy_input != '000':
                yyy_input = '000'
            
            # Get absolute path for data directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(script_dir, CONFIG.DATA_DIR)
            os.makedirs(data_dir, exist_ok=True)
            
            folder = os.path.join(data_dir, f'CFIS.{xxx_input}.{yyy_input}')
            catalog_path = os.path.join(folder, f'full_catalog_{gal_id}_{coords_ra:.4f}_{coords_dec:.4f}.cat')

        except Exception as e:
            logging.error(f"Error during tile check for PGC {gal_id}: {e}")
            continue

        # 2. Download Data
        if not args.skip_download and not os.path.exists(folder):
            try:
                setup_directories(folder)
                download_data(coords_ra, coords_dec, xxx_input, yyy_input, folder)
            except Exception as e:
                logging.error(f"Download failed for PGC {gal_id}. Skipping detection/analysis.")
                continue

        # 3. Detect Sources and Filter Light
        if not args.skip_detect and os.path.exists(folder):
            try:
                process_and_detect(folder, gal_id, coords_ra, coords_dec, dist_mpc) 
            except Exception as e:
                logging.error(f"Detection failed for PGC {gal_id}: {e}. Skipping analysis.")
                continue

        # 4. Analysis and Plotting
        if os.path.exists(catalog_path):
            try:
                # --- WCS CALCULATION: CRITICAL FIX ---
                # To find the true center, we need WCS from one of the FITS files.
                # Assuming 'r' band is available in the folder.
                r_fits_path = os.path.join(folder, f'r_fin_{xxx_input}.{yyy_input}.fits')
                if not os.path.exists(r_fits_path):
                    r_fits_path = os.path.join(folder, f'CFIS.{xxx_input}.{yyy_input}.r.fits')

                if os.path.exists(r_fits_path):
                    with fits.open(r_fits_path) as hdul:
                        w = WCS(hdul[0].header)
                        # Convert galaxy RA/DEC to pixel coordinates
                        x_center_pix, y_center_pix = w.all_world2pix(coords_ra, coords_dec, 0)
                else:
                    logging.warning(f"FITS file not found for WCS, using hardcoded center (5000, 5000).")
                    x_center_pix, y_center_pix = 5000, 5000 
                # --- END WCS CALCULATION ---
                
                final_catalog = filter_globular_clusters(
                    catalog_path, 
                    x_center_pix, 
                    y_center_pix, 
                    gal_id, 
                    dist_mpc
                )
                
                generate_plots(folder, gal_id, coords_ra, coords_dec, dist_mpc, final_catalog)
                
            except Exception as e:
                logging.error(f"Analysis/Plotting failed for PGC {gal_id}: {e}")
                
        else:
             logging.warning(f"Catalog not found for PGC {gal_id}. Cannot run analysis.")


    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"\n--- Pipeline finished. Total time: {elapsed_time:.2f} seconds ---")


if __name__ == "__main__":
    # Create plots directory in data folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, CONFIG.DATA_DIR)
    plots_dir = os.path.join(data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    main()