import os
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

from .creation import Feature


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
        sigwx_object = QvaObject(lons, lats, 1)
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
