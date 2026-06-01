# Methodology Documentation
## Sentinel-1 SAR Change Detection — Concentrate Stockpile Monitoring

**Version:** 1.0  
**Study Site:** Chuquicamata Copper Mine, Chile  
**Author:** Axel Franke

---

## 1. Data Source

**Sentinel-1 IW GRD** (Interferometric Wide Swath, Ground Range Detected)

- Free of charge via Copernicus Open Access Hub (scihub.copernicus.eu) or Alaska Satellite Facility (ASF)
- Acquisition mode: IW (Interferometric Wide)
- Product type: GRDH (High Resolution)
- Polarisation: VV (primary for surface change; VH optionally included)
- Nominal ground resolution: ~10x10m (after multi-looking)
- Revisit time: 6-12 days over most land areas (combined Sentinel-1A/1B)
- Radiometric accuracy: ~1 dB absolute, better for relative comparisons

---

## 2. Preprocessing Chain

All preprocessing performed in ESA SNAP 9.0 via GPT (Graph Processing Tool).

### Step 1: Apply Orbit File
Corrects satellite position and velocity metadata using precise orbit determination (POD) files. Precise orbits are released ~20 days after acquisition and are preferred over restituted orbits (available within hours but less accurate). Accurate geolocation is critical for co-registration across the image stack.

**Parameter:** Sentinel Precise (Auto Download)

### Step 2: Thermal Noise Removal
Removes the additive thermal noise floor that varies across the swath in range direction. Without this step, the noise pattern introduces a systematic range-dependent bias into the calibrated backscatter values, causing false gradients across the AOI.

**Parameter:** VV polarisation; removeeThermalNoise = true

### Step 3: Radiometric Calibration
Converts raw digital number (DN) pixel values to physical backscatter coefficient sigma-naught (σ⁰). This step makes images acquired at different times, dates, and orbital geometries directly comparable in physical units. Output is linear sigma-naught.

**Parameter:** outputSigmaBand = true

### Step 4: Speckle Filtering
SAR images are inherently affected by speckle — a multiplicative noise caused by coherent interference of backscatter from multiple scatterers within a resolution cell. Speckle filtering reduces this noise before terrain correction.

**Filter chosen:** Refined Lee (5x5 window)  
**Rationale:** Refined Lee adapts to local image statistics and preserves edges better than standard Lee or Boxcar filters. Edge preservation is important for stockpile boundary definition — the transition between stockpile surface and surrounding ground must remain detectable.

Applied **after calibration** and **before terrain correction** to avoid resampling artefacts affecting the filtered values.

### Step 5: Terrain Correction (Range Doppler)
Orthorectifies the SAR image to a map projection, correcting for geometric distortions caused by topographic relief (foreshortening, layover, shadow). Without terrain correction, pixels in mountainous terrain are mislocated, and co-registration between dates fails.

**DEM used:** Copernicus GLO-30 (30m global DEM)  
**Rationale:** The Chuquicamata site sits in the Atacama at ~2,850m elevation with significant topographic relief around the pit. GLO-30 provides better accuracy than SRTM 3Sec in this terrain type (SRTM has known data voids and lower vertical accuracy in steep areas).

**Output projection:** UTM Zone 19S (EPSG:32719)  
**Output pixel spacing:** 10m

---

## 3. AOI Definition

The Area of Interest (AOI) was defined **independently of SAR imagery** using:

1. Sentinel-2 RGB + SWIR composite to identify the concentrate stockpile spectral signature (high SWIR reflectance, distinct from tailings which have lower reflectance and finer texture)
2. Google Earth historical imagery to confirm stockpile location, access roads, and boundary evolution over the analysis period
3. Conservative polygon that excludes: processing plant infrastructure, haul road corridors, adjacent tailings facility, pit rim

**Why this matters:** Defining the AOI using SAR imagery risks confirmation bias — the analyst may unconsciously draw the polygon to include areas that already show interesting backscatter variation. Using optical imagery as the independent reference and documenting the rationale makes the AOI selection defensible to a technically sophisticated audience.

---

## 4. Change Detection Methodology

### Time Series Construction
Mean sigma-naught backscatter is extracted within the AOI for each acquisition date and converted to dB scale:

```
σ⁰_dB = 10 × log₁₀(σ⁰_linear)
```

### Log-Ratio Change Detection
Change between consecutive image pairs is computed as the log-ratio (equivalent to dB difference):

```
ΔdB = σ⁰_dB(t₂) - σ⁰_dB(t₁)
```

**Why log-ratio:** Under the no-change hypothesis, the log-ratio of SAR backscatter values is approximately normally distributed. This makes statistical thresholding well-defined. The metric is symmetric (equal sensitivity to increases and decreases) and physically interpretable.

### Significance Thresholding
The significance threshold is derived empirically from stable reference pixels — areas within the scene expected to show no change (paved roads, bare rock outcrops away from the active site). The standard deviation of log-ratio values in these reference areas characterises the expected no-change distribution. A 2σ threshold is applied.

**Where no reference polygon is available**, a conservative fixed threshold of 2.0 dB is used. This is documented explicitly.

---

## 5. Volumetric Estimation

### Conversion Approach
Backscatter change is converted to estimated height change using an empirical sensitivity parameter, then to volume and mass:

```
Δh = ΔdB / sensitivity (dB/m)
ΔV = Δh × AOI_area
Δm = ΔV × bulk_density
```

### Assumptions and Sources

| Parameter | Value | Source |
|---|---|---|
| Bulk density (copper concentrate) | 1.8 t/m³ | ICSG (2019), typical range 1.6-2.0 |
| Backscatter sensitivity | 0.5 dB/m | Conservative estimate; calibrate against survey data if available |
| AOI area | Calculated from validated polygon | |

### Error Bounds
Uncertainty is propagated from three sources:
- **Spatial resolution** (~10m pixel, sub-pixel averaging within AOI)
- **Bulk density** (±10% relative uncertainty)
- **Sensitivity parameter** (±30% relative uncertainty — dominant term)

Combined relative uncertainty: approximately 40-50% on absolute mass estimates.

### Interpretation Guidance
- **Directional signal** (accumulation vs. drawdown) is reliable
- **Absolute tonnage figures** are order-of-magnitude estimates only
- Results should be expressed as ranges, not point estimates
- Conditions that can cause unreliable estimates:
  - Rainfall events (surface moisture changes backscatter independently of volume)
  - Wind-blown surface disturbance
  - Changes in surface material composition
  - Acquisition geometry changes (different orbital track)

---

## 6. Validation Against Public Data

SAR change events are cross-referenced against publicly available quarterly operational data (Codelco operational reviews for Chuquicamata).

For each significant SAR-detected event:
- Date of detection
- Direction (accumulation / drawdown)
- Magnitude of estimated change
- Subsequent quarterly production figure
- Whether the quarterly figure directionally confirms the SAR signal
- Lead time in days between SAR detection and public reporting date

**Accurate negative results are documented.** If the SAR signal does not precede or confirm the subsequent public data, this is reported explicitly. The methodology is only useful if its failure modes are understood.

---

## 7. Reproducibility

This analysis is fully reproducible from raw Sentinel-1 GRD data downloaded from Copernicus Open Access Hub. All processing parameters are documented in this file and in the notebook code. A second qualified SAR engineer should be able to produce the same outputs given the same input data.

**Software versions:**
- ESA SNAP 9.0
- Python 3.10
- rasterio 1.3+
- GDAL 3.4+
- geopandas 0.12+

---

## 8. Known Limitations and Failure Modes

| Limitation | Impact | Mitigation |
|---|---|---|
| 10m spatial resolution | Limits sub-stockpile precision; edge pixels mixed | AOI set conservatively inward from stockpile boundary |
| Surface moisture sensitivity | Rainfall can produce false change signal | Flag acquisitions within 48hrs of significant precipitation |
| Layover/foreshortening | Affects pixels at steep terrain near pit rim | AOI excludes terrain-affected areas |
| Single polarisation (VV) | Reduced sensitivity vs. dual-pol analysis | VH band included where available for cross-check |
| Volumetric conversion uncertainty | ~40-50% relative error on mass estimates | Report as ranges; directional signal is primary output |
| No independent ground survey | Cannot calibrate sensitivity parameter | Document assumption; flag for future calibration opportunity |
