import os
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

from .creation import Feature


def to_shapefile(features: Iterable[Feature], filename: Path):
    """
    Write a list of polygon features to ESRI shapefile format

    Arguments:
        features: List of features.
        filename: Path to shapefile.
    """
    import fiona

    # Convert to shapefile records, checking that there is enough homogeneity
    # to be saved to a shapefile
    geomtypes = set()
    proptypes = {}
    records = []
    for feature in features:
        geometry = feature.geometry.__geo_interface__
        geomtypes.add(geometry["type"])

        properties = feature.properties
        for key, val in properties.items():
            proptypes[key] = val.__class__.__name__

        records.append({"geometry": geometry, "properties": properties})

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

    schema = {"geometry": geomtype, "properties": proptypes}
    with fiona.open(filename, "w", driver="ESRI Shapefile", schema=schema) as file:
        file.writerecords(records)


def from_shapefile(filename: Path) -> Iterable[Feature]:
    """
    Load features from a shapefile

    Arguments:
        filename: Path to shapefile.

    Returns:
        List of features contained in the shapefile
    """
    import fiona

    with fiona.open(str(filename), "r") as file:
        features = getattr(file, "__geo_interface__", file)
        if isinstance(features, dict) and features.get("type") == "FeatureCollection":
            features = features["features"]
        for feature in features:
            feature = getattr(feature, "__geo_interface__", feature)
            if not feature["geometry"]:
                continue
            yield Feature(feature["geometry"], feature["properties"])


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
