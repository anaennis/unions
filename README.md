# UNIONS Globular Cluster Detection Pipeline

An automated pipeline for detecting and analyzing globular cluster candidates in CFIS/UNIONS survey data using DOLPHOT photometry.

## Overview

This pipeline processes multi-band CFIS imaging data to identify globular cluster candidates around target galaxies. It performs galaxy light subtraction, runs DOLPHOT for source detection and photometry, and applies selection cuts to isolate GC candidates based on their photometric and morphological properties.

## Features

- **Automated tile identification**: Checks if target galaxies fall within CFIS survey coverage
- **Data download**: Retrieves u, r, and g-band FITS images from VOSpace
- **Galaxy light subtraction**: Uses dual median filtering to remove large-scale galaxy light and background
- **DOLPHOT photometry**: Performs PSF-fitted photometry on filtered images
- **GC selection**: Applies color, magnitude, morphological, and radial cuts to identify GC candidates
- **Visualization**: Generates color-magnitude diagrams and radial distribution plots

## Requirements

- Python 3.7+
- VOSpace client (`vcp` command)
- DOLPHOT
- Python packages: `numpy`, `astropy`, `scipy`, `pandas`, `matplotlib`

## Configuration

Edit `CONFIG.py` to customize:
- Data storage location (`DATA_DIR`)
- Tile listing files (`TILE_LISTINGS`)
- Image filtering parameters (`FILTER_LARGE_SCALE`, `FILTER_SMALL_SCALE`)
- DOLPHOT parameters (aperture sizes, PSF settings, zeropoints)
- GC selection criteria (color cuts, magnitude limits, radial extent)

## Usage

### Basic Run
```bash
./total.py --catalogue galaxy_catalogue.txt
```

### Skip Download (if data already present)
```bash
./total.py --catalogue galaxy_catalogue.txt --skip-download
```

### Skip Detection (if catalogs already generated)
```bash
./total.py --catalogue galaxy_catalogue.txt --skip-detect
```

## Input Catalogue Format

The input catalogue should be a whitespace-separated text file with 4 columns:
```
RA(deg)  DEC(deg)  PGC_ID  Distance(Mpc)
150.12   2.21      12345   10.5
185.72   -3.46     67890   15.2
```

- No header row
- One galaxy per line
- Comments starting with `#` are allowed

## Output Structure

Data and outputs are organized in the `../data/` directory:
```
data/
├── CFIS.XXX.YYY/              # Tile folder (one per CFIS tile)
│   ├── CFIS.XXX.YYY.r.fits    # Downloaded r-band image
│   ├── CFIS.XXX.YYY.u.fits    # Downloaded u-band image
│   ├── CFIS.XXX.YYY.g.fits    # Downloaded g-band image
│   ├── r_fin_XXX.YYY.fits     # Filtered r-band image
│   ├── u_fin_XXX.YYY.fits     # Filtered u-band image
│   ├── g_fin_XXX.YYY.fits     # Filtered g-band image
│   └── full_catalog_*.cat     # DOLPHOT output catalog
└── plots/
    ├── cmd_pgc12345.png       # Color-magnitude diagram
    └── radial_pgc12345.png    # Radial distribution histogram
```

## Pipeline Stages

1. **Tile Coverage Check**: Verifies if galaxy coordinates fall within available CFIS tiles
2. **Data Download**: Retrieves FITS images from VOSpace (if not present)
3. **Image Processing**: Applies median filtering to subtract galaxy light
4. **Source Detection**: Runs DOLPHOT photometry on filtered images
5. **GC Selection**: Filters catalog based on photometric and morphological criteria
6. **Visualization**: Generates diagnostic plots

## Notes

- The pipeline automatically skips image processing if filtered images already exist
- DOLPHOT parameters are synchronized with `CONFIG.py` using `check_dolphot_params.py`
- All file paths are configured to keep data separate from the code directory

## Logging

Pipeline activity is logged to `analysis_pipeline.log` with detailed information about each processing stage, errors, and warnings.

---

For questions or issues, please refer to the inline documentation in `total.py` and `CONFIG.py`.
