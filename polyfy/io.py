import os
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import shapely.geometry as sgeom

from .creation import Feature
from .util import wrap_longitudes, unwrap_longitudes


def _find_schema(records: Iterable[dict]):
    """
    Identify a schema for a homogeneous collection of features

    Argument:
        records: List of features.

    Returns:
        Schema suitable for eg GeoJSON or Shapefile file formats
    """
    # Get types
    geomtypes = set()
    proptypes = {}
    for record in records:
        geomtypes.add(record["geometry"]["type"])
        for key, val in record["properties"].items():
            proptypes[key] = val.__class__.__name__

    # Determine the single geometry type, if possible
    if geomtypes <= {"Polygon", "MultiPolygon"}:
        # Note that the empty set (ie no features at all) is included here,
        # which is fine
        geomtype = "Polygon"
    elif geomtypes <= {"LineString", "MultiLineString"}:
        geomtype = "LineString"
    elif geomtypes == {"Point"} or geomtypes == {"MultiPoint"}:
        # Single types that cannot be mixed
        geomtype = geomtypes.pop()
    else:
        raise ValueError(f"Cannot mix geometry types: {sorted(geomtypes)}")

    return {"geometry": geomtype, "properties": proptypes}


def _to_records(driver: str, features: Iterable[Feature], filename: Path):
    import fiona

    records = [
        {"geometry": f.geometry.__geo_interface__, "properties": f.properties}
        for f in features
    ]
    schema = _find_schema(records)

    with fiona.open(str(filename), "w", driver=driver, schema=schema) as file:
        file.writerecords(records)


def _from_records(filename: Path) -> Iterable[Feature]:
    import fiona

    with fiona.open(str(filename), "r") as file:
        features = getattr(file, "__geo_interface__", file)
        if isinstance(features, dict) and features.get("type") == "FeatureCollection":
            features = features["features"]
        for feature in features:
            feature = getattr(feature, "__geo_interface__", feature)
            if not feature["geometry"]:
                continue
            yield Feature(feature["geometry"], dict(feature["properties"]))


def to_shapefile(features: Iterable[Feature], filename: Path):
    """
    Write a list of polygon features to ESRI shapefile format

    Arguments:
        features: List of features.
        filename: Path to shapefile.
    """
    _to_records("ESRI Shapefile", features, filename)


def from_shapefile(filename: Path) -> Iterable[Feature]:
    """
    Load features from a shapefile

    Arguments:
        filename: Path to shapefile.

    Returns:
        List of features contained in the shapefile
    """
    return _from_records(filename)


def to_geojson(features: Iterable[Feature], filename: Path):
    """
    Write a list of polygon features to GeoJSON format

    Arguments:
        features: List of features.
        filename: Path to JSON file.
    """
    _to_records("GeoJSON", features, filename)


def from_geojson(filename: Path) -> Iterable[Feature]:
    """
    Load features from GeoJSON

    Arguments:
        filename: Path to JSON file.

    Returns:
        List of features contained in the file
    """
    return _from_records(filename)


def to_iwxxm(
    features: Iterable[Feature],
    out_dir: Path,
    model_time: datetime,
    valid_time: datetime,
    fragments=False,
):
    """
    Write a list of polygon features to IWXXM format

    Arguments:
        features: List of features.
        out_dir: Directory in which to save files.  Filenames are generated
            based on other parameters.
        model_time: Model forecast reference time.
        valid_time: Data validity time.
        fragments: Whether to save one file per feature.
    """
    from sigwx_objects.file_formats import iwxxm
    from sigwx_objects.common.classes import QvaObject

    # Convert times to expected types
    day = datetime(*model_time.timetuple()[:3])
    time = str(model_time.hour)
    lead_time = str(int((valid_time - model_time).total_seconds()) // 3600)

    # Convert features to expected type
    def to_sigwx_object(geometry, properties):
        lons, lats = map(np.array, geometry.exterior.xy)
        sigwx_object = QvaObject(wrap_longitudes(lons, -180), lats, 1)
        sigwx_object.severity = properties["severity"].replace(" ", "")
        sigwx_object.base = properties.get("base", 0)
        sigwx_object.top = properties.get("top", 600)
        return sigwx_object

    features = [to_sigwx_object(f.geometry, f.properties) for f in features]

    iwxxm.write_iwxxm(
        features,
        "QVA",
        "EGRR",
        "London",
        day,
        time,
        lead_time,
        out_dir,
        fragments=fragments,
    )


def from_iwxxm(filename: Path) -> Iterable[Feature]:
    """
    Load features from IWXXM files

    Arguments:
        filename: Path to IWXXM files - either a complete xml file or a
            directory containing xml fragments.

    Returns:
        List of features contained in the IWXXM files
    """
    from sigwx_objects.file_formats import iwxxm

    # Split filename so that sigwx can reconstruct it
    datadir, filename = os.path.split(filename)
    datadir += "/"
    filename, ext = os.path.splitext(filename)
    parts = str(filename).split("_")
    if (
        len(parts) != 7
        or parts[:2] != ["sigwx", "iwxxm"]
        or not parts[-1].startswith("T+")
    ):
        raise ValueError(f"{filename} does not appear to match sigwx naming convention")
    model_time = datetime.strptime(parts[4] + parts[5], "%Y%m%d%HZ")
    lead = parts[-1][2:]  # Skip initial 'T+'
    hour = str(model_time.hour)

    # Read features
    read = iwxxm.read_iwxxm if ext == ".xml" else iwxxm.read_iwxxm_fragments
    features = read("EGRR", ["QVA"], datadir, model_time, hour, lead, 0).get("QVA", [])

    # Convert to own feature class
    for feature in features:
        geometry = sgeom.Polygon([*zip(unwrap_longitudes(feature.lons), feature.lats)])
        properties = {
            "severity": feature.severity,
            "base": feature.base,
            "top": feature.top,
        }
        yield Feature(geometry, properties)
