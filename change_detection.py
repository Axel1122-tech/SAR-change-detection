"""
change_detection.py
-------------------
SAR backscatter change detection utilities for Sentinel-1 GRD time series.

Implements log-ratio differencing with empirical significance thresholding
derived from stable reference pixels. Designed for industrial site monitoring
(open-pit mines, stockpile areas, port terminals).

Author: Axel Franke
"""

import numpy as np
import rasterio
from rasterio.mask import mask
from pathlib import Path
from typing import Union, Optional
import geopandas as gpd
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Conversion utilities
# ---------------------------------------------------------------------------

def linear_to_db(sigma_linear: np.ndarray) -> np.ndarray:
    """Convert linear sigma-naught backscatter to dB scale.

    Parameters
    ----------
    sigma_linear : np.ndarray
        Linear backscatter values (sigma-naught from radiometric calibration).

    Returns
    -------
    np.ndarray
        Backscatter in dB (10 * log10(sigma_linear)).
        NoData pixels (<=0) are set to NaN.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        db = 10.0 * np.log10(np.where(sigma_linear > 0, sigma_linear, np.nan))
    return db


def db_to_linear(sigma_db: np.ndarray) -> np.ndarray:
    """Convert dB backscatter back to linear scale."""
    return 10.0 ** (sigma_db / 10.0)


# ---------------------------------------------------------------------------
# AOI extraction
# ---------------------------------------------------------------------------

def extract_aoi_values(
    image_path: Union[str, Path],
    aoi_geojson: Union[str, Path],
    band: int = 1,
    to_db: bool = True,
) -> np.ndarray:
    """Extract backscatter values within an AOI polygon.

    Parameters
    ----------
    image_path : str or Path
        Path to a GeoTIFF (terrain-corrected, calibrated Sentinel-1 GRD).
    aoi_geojson : str or Path
        Path to a GeoJSON file containing the AOI polygon.
    band : int
        Band index to read (1-indexed). Default 1 (VV polarisation).
    to_db : bool
        If True, convert linear sigma-naught values to dB before returning.

    Returns
    -------
    np.ndarray
        Flattened array of valid (non-NaN, non-masked) backscatter values
        within the AOI.
    """
    aoi = gpd.read_file(aoi_geojson)

    with rasterio.open(image_path) as src:
        # Reproject AOI to image CRS if necessary
        if aoi.crs != src.crs:
            aoi = aoi.to_crs(src.crs)

        geoms = [geom.__geo_interface__ for geom in aoi.geometry]
        masked_data, _ = mask(src, geoms, crop=True, nodata=np.nan)
        values = masked_data[band - 1].astype(np.float32)

    valid = values[np.isfinite(values)]

    if to_db:
        valid = linear_to_db(valid)
        valid = valid[np.isfinite(valid)]

    return valid


def extract_aoi_mean(
    image_path: Union[str, Path],
    aoi_geojson: Union[str, Path],
    band: int = 1,
    to_db: bool = True,
) -> float:
    """Return the mean backscatter value within an AOI.

    Convenience wrapper around extract_aoi_values. Used to build
    the time series DataFrame for a stack of acquisitions.
    """
    values = extract_aoi_values(image_path, aoi_geojson, band=band, to_db=to_db)
    return float(np.nanmean(values))


# ---------------------------------------------------------------------------
# Log-ratio change detection
# ---------------------------------------------------------------------------

def log_ratio(db_t1: float, db_t2: float) -> float:
    """Compute log-ratio change between two dB backscatter values.

    Log-ratio = dB_t2 - dB_t1

    Positive values indicate backscatter increase (accumulation signal).
    Negative values indicate backscatter decrease (drawdown signal).

    Parameters
    ----------
    db_t1, db_t2 : float
        Mean dB backscatter at time 1 and time 2.

    Returns
    -------
    float
        Log-ratio change in dB.
    """
    return db_t2 - db_t1


def build_time_series(
    image_paths: list,
    acquisition_dates: list,
    aoi_geojson: Union[str, Path],
    reference_geojson: Optional[Union[str, Path]] = None,
    band: int = 1,
) -> pd.DataFrame:
    """Build a backscatter time series DataFrame for a stack of images.

    Parameters
    ----------
    image_paths : list of str or Path
        Ordered list of terrain-corrected GeoTIFF paths.
    acquisition_dates : list of str
        Corresponding acquisition dates (ISO format: YYYY-MM-DD).
    aoi_geojson : str or Path
        AOI polygon over the feature of interest (e.g. stockpile).
    reference_geojson : str or Path, optional
        Stable reference area polygon for empirical threshold derivation.
        If None, threshold is set using a fixed 2 dB heuristic.
    band : int
        Band index (1 = VV, 2 = VH for dual-pol products).

    Returns
    -------
    pd.DataFrame
        Columns: date, backscatter_db, log_ratio, pct_change, significant
    """
    records = []
    for path, date in zip(image_paths, acquisition_dates):
        mean_db = extract_aoi_mean(path, aoi_geojson, band=band, to_db=True)
        records.append({"date": pd.to_datetime(date), "backscatter_db": mean_db})

    df = pd.DataFrame(records).sort_values("date").reset_index(drop=True)

    # Compute log-ratio and percentage change between consecutive pairs
    df["log_ratio"] = df["backscatter_db"].diff()
    df["pct_change"] = (
        (db_to_linear(df["backscatter_db"]) - db_to_linear(df["backscatter_db"].shift(1)))
        / db_to_linear(df["backscatter_db"].shift(1))
        * 100
    )

    # Derive significance threshold from reference pixels or use fixed heuristic
    if reference_geojson is not None:
        ref_changes = []
        for path in image_paths[1:]:
            # Use reference area to characterise no-change distribution
            ref_vals = extract_aoi_values(path, reference_geojson, band=band, to_db=True)
            ref_changes.append(np.nanstd(ref_vals))
        threshold = float(np.mean(ref_changes) * 2.0)  # 2-sigma threshold
    else:
        threshold = 2.0  # conservative fixed threshold (dB)

    df["significant"] = df["log_ratio"].abs() > threshold
    df["direction"] = np.where(
        df["log_ratio"] > 0, "accumulation",
        np.where(df["log_ratio"] < 0, "drawdown", "stable")
    )
    df["threshold_used_db"] = threshold

    return df


# ---------------------------------------------------------------------------
# Volumetric estimation
# ---------------------------------------------------------------------------

def backscatter_to_volume(
    log_ratio_db: float,
    aoi_area_m2: float,
    bulk_density_t_m3: float = 1.8,
    sensitivity_db_per_m: float = 0.5,
) -> dict:
    """Convert a backscatter change to an estimated volumetric change.

    This is an approximate conversion. Backscatter is sensitive to surface
    roughness, moisture, and incidence angle as well as volume. The result
    should be treated as an order-of-magnitude estimate with explicit
    uncertainty, not a precision measurement.

    Parameters
    ----------
    log_ratio_db : float
        Log-ratio backscatter change (dB) over the AOI.
    aoi_area_m2 : float
        Area of the AOI in square metres.
    bulk_density_t_m3 : float
        Assumed bulk density of the material in t/m3.
        Copper concentrate: ~1.8 t/m3 (document source in methodology).
    sensitivity_db_per_m : float
        Empirical sensitivity parameter: dB change per metre of height change.
        Derived from literature or calibration data. Default 0.5 dB/m is
        conservative for open stockpile surfaces.

    Returns
    -------
    dict with keys:
        height_change_m : estimated height change in metres
        volume_change_m3 : estimated volume change in cubic metres
        mass_change_t : estimated mass change in metric tonnes
        error_bounds_t : approximate ±1σ error in metric tonnes
        assumptions : dict of all assumptions used
    """
    height_change_m = log_ratio_db / sensitivity_db_per_m
    volume_change_m3 = height_change_m * aoi_area_m2
    mass_change_t = volume_change_m3 * bulk_density_t_m3

    # Error propagation: combine uncertainty from resolution (~10m pixel),
    # bulk density (±10%), and sensitivity parameter (±30%)
    resolution_error_m = 10.0 / np.sqrt(aoi_area_m2 / 100)  # sub-pixel averaging
    sensitivity_error = 0.30  # 30% uncertainty on sensitivity parameter
    density_error = 0.10      # 10% uncertainty on bulk density

    total_relative_error = np.sqrt(
        sensitivity_error**2 + density_error**2 + (resolution_error_m / max(abs(height_change_m), 0.01))**2
    )
    error_bounds_t = abs(mass_change_t) * total_relative_error

    return {
        "height_change_m": round(height_change_m, 3),
        "volume_change_m3": round(volume_change_m3, 1),
        "mass_change_t": round(mass_change_t, 0),
        "error_bounds_t": round(error_bounds_t, 0),
        "assumptions": {
            "bulk_density_t_m3": bulk_density_t_m3,
            "bulk_density_source": "Typical copper concentrate bulk density (ICSG, 2019)",
            "sensitivity_db_per_m": sensitivity_db_per_m,
            "sensitivity_source": "Conservative estimate; calibrate against known stockpile surveys if available",
            "aoi_area_m2": aoi_area_m2,
            "note": "Absolute tonnage figures are order-of-magnitude. Use directional signal as primary indicator.",
        },
    }
