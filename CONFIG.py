# --- FILE PATHS ---

# Directory for data downloads and outputs (relative to script location or absolute)
DATA_DIR = '../data'

# Tile listing files for coverage check (searched in order)

TILE_LISTINGS = [
    'UNIONS_CFIS_LSB_r_DR5.txt']

# Logging configuration
LOG_FILENAME = 'analysis_pipeline.log'


# --- ASTROMETRIC & PHOTOMETRIC CONSTANTS ---

PIXEL_SCALE_CFIS = 0.187  # Image scale in arcsec/pixel for CFIS data
DEFAULT_MAG_ZEROPOINT = 30.0 # Placeholder value for seeing FWHM in arcsec in case the FITS header does not have one
DEFAULT_SEEING_FWHM = 1.0 # Placeholder, should also be read from header


# --- IMAGE PROCESSING (Light Subtraction/Filtering) ---

# Median filter kernel size for removing large-scale galaxy light
FILTER_LARGE_SCALE = (50, 50) 
# Median filter kernel size for removing small-scale background noise/residual
FILTER_SMALL_SCALE = (8, 8)


# --- DETECTION PARAMETERS (DOLPHOT) ---
# DOLPHOT parameter values (these map to entries in a DOLPHOT parameter file)
Nimg = 3
img_RAPER = 6.0
img_RPSF = 13
img_APRAD = 10
# Background aperture(s) usually provided as two integers (inner, outer)
img_APSKY = (10, 15)
img_RCentroid = 2
PSFPhot = 1
FSat = 0.999
ZERO = 25.0

# Default instrumental gain to write into FITS headers when needed
DEFAULT_GAIN = 1.0

# DOLPHOT parameter file (used when running dolphot)
DOLPHOT_PARAM_FILE = 'dolphot.param'


# --- GLOBULAR CLUSTER (GC) CANDIDATE FILTERING CUTS ---

# Morphological Cut (CLASS_STAR index: 1.0 is point source/star)
GC_CLASS_STAR_THRESHOLD = 0.8

# Photometric Color Cuts (u-r vs g-r space)
GC_U_R_MIN = 0.5
GC_U_R_MAX = 2.5
GC_G_R_MIN = 0.0
GC_G_R_MAX = 1.0

# Quality/Depth Cuts
GC_MAG_ERR_MAX = 0.1   # Maximum magnitude error in R-band
GC_MAG_R_MAX = 25.0    # Brightest magnitude cut
GC_MAG_R_MIN = 19.0    # Faintest magnitude cut (to avoid noise)

# Radial Cut (Physical size from galaxy center)
GC_RADIAL_DISTANCE_KPC = 10.0 # Maximum radial distance from galaxy center in kpc