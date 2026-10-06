"""
Tests for ``polyfy.io``
"""

import shapely.geometry as sgeom
from polyfy import Feature, io


def test_save_load_geojson(tmp_path):
    "Round trip save -> load GeoJSON"
    geometry = sgeom.box(1, 2, 3, 4)
    features = [Feature(geometry, {"test_data": 2})]
    filename = tmp_path / "test.json"
    io.to_geojson(features, filename)
    result = list(io.from_geojson(filename))
    assert features == result
