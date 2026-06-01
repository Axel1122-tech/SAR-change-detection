# Sentinel-1 SAR Change Detection — Open-Pit Mine Stockpile Monitoring

A reproducible Python workflow for detecting and quantifying surface change at industrial sites using Sentinel-1 GRD SAR imagery. Demonstrated on the **Chuquicamata open-pit copper mine, Chile** (22°18'S, 68°54'W).

---

## Overview

This project implements a full SAR change detection pipeline from raw Sentinel-1 GRD acquisition through to volumetric change estimation and validation against operational ground truth. The methodology is designed for industrial site monitoring where optical imagery is unavailable or insufficient (cloud cover, night acquisitions, surface moisture sensitivity).

**Primary use case:** Concentrate stockpile volume change estimation at open-pit mining operations.

---

## Workflow

```
Raw Sentinel-1 GRD
        │
        ▼
1. Preprocessing (SNAP / snappy)
   - Orbit file application
   - Thermal noise removal
   - Radiometric calibration → σ⁰
   - Speckle filtering (Refined Lee)
   - Terrain correction (Copernicus GLO-30 DEM)
        │
        ▼
2. AOI Definition & Validation
   - Optical cross-reference (Sentinel-2 / Google Earth)
   - Stockpile polygon delineation
   - Documentation of AOI rationale
        │
        ▼
3. Multi-temporal Stack Construction
   - Co-registration of image pairs
   - dB conversion of σ⁰ values
   - Time series extraction within AOI
        │
        ▼
4. Change Detection
   - Log-ratio differencing between consecutive pairs
   - Statistical significance thresholding
   - Change event classification (accumulation / drawdown)
        │
        ▼
5. Volumetric Estimation
   - Backscatter-to-volume conversion
   - Bulk density assumption documentation
   - Error bound propagation
        │
        ▼
6. Validation
   - Cross-reference against public operational data
   - Lead time analysis (SAR detection vs. public reporting)
```

---

## Repository Structure

```
sentinel1-change-detection/
├── notebooks/
│   ├── 01_preprocessing.ipynb          # SNAP preprocessing via snappy
│   ├── 02_aoi_definition.ipynb         # AOI delineation and validation
│   ├── 03_change_detection.ipynb       # Log-ratio change detection pipeline
│   ├── 04_volumetric_estimation.ipynb  # Backscatter to volume conversion
│   └── 05_validation.ipynb             # Ground truth cross-reference
├── src/
│   ├── preprocessing.py                # Preprocessing utilities
│   ├── change_detection.py             # Change detection functions
│   ├── volumetric.py                   # Volumetric conversion utilities
│   └── visualization.py               # Plotting and map output functions
├── data/
│   └── sample/                         # Sample AOI GeoJSON and metadata
├── outputs/                            # Generated figures and CSVs
├── docs/
│   └── methodology.md                  # Full methodology documentation
├── environment.yml                     # Conda environment
├── requirements.txt
└── README.md
```

---

## Study Area

**Chuquicamata Copper Mine, Chile**
- Coordinates: 22°18'S, 68°54'W
- One of the world's largest open-pit copper mines (operated by Codelco)
- Sentinel-1 IW coverage: multiple passes per week (ascending + descending)
- Concentrate stockpile area clearly identifiable in both SAR and optical imagery

The stockpile AOI was defined using Sentinel-2 RGB and SWIR composites to identify the concentrate storage area independently of the SAR signal, then validated against publicly available Google Earth imagery.

---

## Key Technical Decisions

| Decision | Choice | Rationale |
|---|---|---|
| DEM for terrain correction | Copernicus GLO-30 | Higher accuracy than SRTM in steep Andean terrain |
| Speckle filter | Refined Lee | Preserves edges better than standard Lee; important for stockpile boundary definition |
| Change metric | Log-ratio (dB difference) | Symmetric, normally distributed under no-change hypothesis; standard in SAR literature |
| Significance threshold | 2σ of stable reference pixels | Empirically derived from site-specific no-change areas; avoids hardcoded thresholds |
| Volumetric conversion | Backscatter proxy with documented bulk density | Conservative approach with explicit error bounds; absolute values treated as order-of-magnitude |

---

## Limitations

- Sentinel-1 IW GRD spatial resolution (~10x10m) limits sub-stockpile volumetric precision
- Backscatter is sensitive to surface moisture; rainfall events can produce false change signals
- Layover and foreshortening effects in steep terrain adjacent to the pit can affect AOI boundary pixels
- Volumetric estimates carry compounding uncertainty from resolution, bulk density assumptions, and geometric conversion

All limitations are documented explicitly in `docs/methodology.md` and within each notebook.

---

## Requirements

```
Python >= 3.9
rasterio >= 1.3
GDAL >= 3.4
geopandas >= 0.12
numpy >= 1.23
scipy >= 1.9
matplotlib >= 3.6
snappy (ESA SNAP Python interface) — see installation notes
```

See `environment.yml` for full conda environment specification.

---

## Usage

```bash
# Clone the repo
git clone https://github.com/axelfranke/sentinel1-change-detection.git
cd sentinel1-change-detection

# Create environment
conda env create -f environment.yml
conda activate sar-change-detection

# Run notebooks in order
jupyter lab notebooks/
```

Sentinel-1 GRD data can be downloaded free of charge from [Copernicus Open Access Hub](https://scihub.copernicus.eu/) or via the [Alaska Satellite Facility](https://asf.alaska.edu/).

---

## Author

**Axel Franke** — GIS & Remote Sensing Specialist  
Sentinel-1 SAR processing | SNAP | Google Earth Engine | Python geospatial stack  
[LinkedIn](https://linkedin.com/in/axelfranke) | [Upwork](https://upwork.com)
