from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

from .creation import Feature


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
