"""
preprocessing.py
----------------
Sentinel-1 GRD preprocessing utilities.

The full preprocessing chain (orbit file, thermal noise removal, radiometric
calibration, speckle filtering, terrain correction) is best run in ESA SNAP
using the snappy Python interface or the SNAP Graph Builder (GPT).

This module provides:
  - A snappy-based automated preprocessing function
  - A GPT graph XML generator for batch processing without snappy
  - Helper utilities for stack management and file naming

Author: Axel Franke
"""

from pathlib import Path
from typing import Union, List, Optional
import subprocess
import xml.etree.ElementTree as ET
import json


# ---------------------------------------------------------------------------
# SNAP GPT graph-based preprocessing (no snappy dependency)
# ---------------------------------------------------------------------------

def build_gpt_graph(
    input_path: str,
    output_path: str,
    dem_name: str = "Copernicus 30m Global DEM",
    speckle_filter: str = "Refined Lee",
    polarisation: str = "VV",
) -> str:
    """Generate a SNAP GPT XML processing graph for Sentinel-1 GRD preprocessing.

    The graph applies the standard preprocessing chain:
      Apply-Orbit-File -> ThermalNoiseRemoval -> Calibration ->
      Speckle-Filter -> Terrain-Correction

    Parameters
    ----------
    input_path : str
        Path to the input Sentinel-1 GRD product (.zip or .SAFE directory).
    output_path : str
        Path for the output GeoTIFF.
    dem_name : str
        DEM to use for terrain correction.
        Options: "Copernicus 30m Global DEM", "SRTM 3Sec", "SRTM 1Sec HGT"
    speckle_filter : str
        Speckle filter type. Options: "Refined Lee", "Lee", "Boxcar"
    polarisation : str
        Polarisation to process. "VV", "VH", or "VV VH"

    Returns
    -------
    str
        XML string of the processing graph. Write to file and pass to gpt.
    """
    graph_template = f"""<graph id="Graph">
  <version>1.0</version>

  <node id="Read">
    <operator>Read</operator>
    <sources/>
    <parameters>
      <file>{input_path}</file>
    </parameters>
  </node>

  <node id="Apply-Orbit-File">
    <operator>Apply-Orbit-File</operator>
    <sources>
      <sourceProduct refid="Read"/>
    </sources>
    <parameters>
      <orbitType>Sentinel Precise (Auto Download)</orbitType>
      <polyDegree>3</polyDegree>
      <continueOnFail>false</continueOnFail>
    </parameters>
  </node>

  <node id="ThermalNoiseRemoval">
    <operator>ThermalNoiseRemoval</operator>
    <sources>
      <sourceProduct refid="Apply-Orbit-File"/>
    </sources>
    <parameters>
      <selectedPolarisations>{polarisation}</selectedPolarisations>
      <removeThermalNoise>true</removeThermalNoise>
      <reIntroduceThermalNoise>false</reIntroduceThermalNoise>
    </parameters>
  </node>

  <node id="Calibration">
    <operator>Calibration</operator>
    <sources>
      <sourceProduct refid="ThermalNoiseRemoval"/>
    </sources>
    <parameters>
      <selectedPolarisations>{polarisation}</selectedPolarisations>
      <outputSigmaBand>true</outputSigmaBand>
      <outputGammaBand>false</outputGammaBand>
      <outputBetaBand>false</outputBetaBand>
    </parameters>
  </node>

  <node id="Speckle-Filter">
    <operator>Speckle-Filter</operator>
    <sources>
      <sourceProduct refid="Calibration"/>
    </sources>
    <parameters>
      <filter>{speckle_filter}</filter>
      <filterSizeX>5</filterSizeX>
      <filterSizeY>5</filterSizeY>
      <dampingFactor>2</dampingFactor>
      <estimateENL>true</estimateENL>
      <enl>1.0</enl>
      <numLooksStr>1</numLooksStr>
      <targetWindowSizeStr>3x3</targetWindowSizeStr>
      <sigmaStr>0.9</sigmaStr>
      <anSize>50</anSize>
    </parameters>
  </node>

  <node id="Terrain-Correction">
    <operator>Terrain-Correction</operator>
    <sources>
      <sourceProduct refid="Speckle-Filter"/>
    </sources>
    <parameters>
      <demName>{dem_name}</demName>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
      <demResamplingMethod>BILINEAR_INTERPOLATION</demResamplingMethod>
      <imgResamplingMethod>BILINEAR_INTERPOLATION</imgResamplingMethod>
      <pixelSpacingInMeter>10.0</pixelSpacingInMeter>
      <mapProjection>AUTO:42001</mapProjection>
      <nodataValueAtSea>false</nodataValueAtSea>
      <saveDEM>false</saveDEM>
      <saveLatLon>false</saveLatLon>
      <saveIncidenceAngleFromEllipsoid>false</saveIncidenceAngleFromEllipsoid>
      <saveLocalIncidenceAngle>false</saveLocalIncidenceAngle>
      <saveProjectedLocalIncidenceAngle>false</saveProjectedLocalIncidenceAngle>
      <saveSelectedSourceBand>true</saveSelectedSourceBand>
    </parameters>
  </node>

  <node id="Write">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="Terrain-Correction"/>
    </sources>
    <parameters>
      <file>{output_path}</file>
      <formatName>GeoTIFF</formatName>
    </parameters>
  </node>

</graph>"""
    return graph_template


def run_gpt_preprocessing(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    gpt_executable: str = "gpt",
    dem_name: str = "Copernicus 30m Global DEM",
    speckle_filter: str = "Refined Lee",
    polarisation: str = "VV",
    graph_dir: Optional[Union[str, Path]] = None,
) -> Path:
    """Run SNAP GPT preprocessing on a single Sentinel-1 GRD product.

    Parameters
    ----------
    input_path : str or Path
        Input .zip or .SAFE product path.
    output_path : str or Path
        Output GeoTIFF path.
    gpt_executable : str
        Path to the SNAP gpt executable. If SNAP is on PATH, "gpt" works.
    dem_name : str
        DEM name for terrain correction.
    speckle_filter : str
        Speckle filter type.
    polarisation : str
        Polarisation band.
    graph_dir : str or Path, optional
        Directory to write temporary graph XML files. Defaults to output directory.

    Returns
    -------
    Path
        Path to the output GeoTIFF.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    graph_dir = Path(graph_dir) if graph_dir else output_path.parent
    graph_path = graph_dir / f"{output_path.stem}_graph.xml"

    xml_content = build_gpt_graph(
        input_path=str(input_path),
        output_path=str(output_path),
        dem_name=dem_name,
        speckle_filter=speckle_filter,
        polarisation=polarisation,
    )

    with open(graph_path, "w") as f:
        f.write(xml_content)

    cmd = [gpt_executable, str(graph_path), "-e"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"SNAP GPT preprocessing failed for {input_path}.\n"
            f"STDERR: {result.stderr}"
        )

    # Clean up graph file
    graph_path.unlink(missing_ok=True)

    return output_path


def batch_preprocess(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    gpt_executable: str = "gpt",
    pattern: str = "*.zip",
    **kwargs,
) -> List[Path]:
    """Batch preprocess all Sentinel-1 GRD products in a directory.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing Sentinel-1 GRD .zip files.
    output_dir : str or Path
        Output directory for preprocessed GeoTIFFs.
    gpt_executable : str
        Path to gpt.
    pattern : str
        Glob pattern for input files.
    **kwargs
        Additional keyword arguments passed to run_gpt_preprocessing.

    Returns
    -------
    list of Path
        Paths to successfully preprocessed output files.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted(input_dir.glob(pattern))
    if not inputs:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {input_dir}")

    outputs = []
    for i, input_path in enumerate(inputs, 1):
        output_path = output_dir / f"{input_path.stem}_preprocessed.tif"
        print(f"[{i}/{len(inputs)}] Processing: {input_path.name}")
        try:
            out = run_gpt_preprocessing(input_path, output_path, gpt_executable, **kwargs)
            outputs.append(out)
            print(f"  -> Done: {out.name}")
        except RuntimeError as e:
            print(f"  -> FAILED: {e}")

    return outputs


# ---------------------------------------------------------------------------
# File naming and stack management utilities
# ---------------------------------------------------------------------------

def parse_s1_filename(filename: str) -> dict:
    """Parse acquisition metadata from a Sentinel-1 filename.

    Standard format:
    S1A_IW_GRDH_1SDV_20240115T091234_20240115T091259_051234_063ABC_1234.SAFE

    Parameters
    ----------
    filename : str
        Sentinel-1 product filename (with or without .SAFE/.zip extension).

    Returns
    -------
    dict
        Parsed metadata: mission, mode, product_type, date, time, orbit.
    """
    name = Path(filename).stem.replace(".SAFE", "")
    parts = name.split("_")

    if len(parts) < 6:
        return {"raw": filename, "error": "Could not parse filename"}

    return {
        "mission": parts[0],           # S1A or S1B
        "mode": parts[1],              # IW, EW, SM, WV
        "product_type": parts[2],      # GRDH, GRDM, SLC
        "polarisation": parts[3],      # 1SDV, 1SSH, etc.
        "start_datetime": parts[4],    # YYYYMMDDTHHmmss
        "stop_datetime": parts[5],     # YYYYMMDDTHHmmss
        "acquisition_date": parts[4][:8],  # YYYYMMDD
    }


def sort_stack_by_date(image_paths: List[Union[str, Path]]) -> List[Path]:
    """Sort a list of preprocessed Sentinel-1 GeoTIFFs by acquisition date.

    Expects filenames derived from standard Sentinel-1 naming convention.
    """
    def _date_key(p):
        meta = parse_s1_filename(Path(p).name)
        return meta.get("acquisition_date", "00000000")

    return sorted([Path(p) for p in image_paths], key=_date_key)
