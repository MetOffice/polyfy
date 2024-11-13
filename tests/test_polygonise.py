"""
Tests for polyfy.creation.polygonise_region
"""

import numpy as np
import pytest
import shapely.geometry as sgeom

from polyfy.creation import polygonise_region


def test_dateline_simple():
    "Simple polygon that crosses the dateline"
    width, height = 360, 180
    extent = (width * 7 // 8, height * 3 // 8, width * 9 // 8, height * 5 // 8)
    data = np.zeros((height, width), dtype=bool)
    # Use closed ranges instead of half-open, to match shapely box
    data[extent[1] : extent[3] + 1, extent[0] % width :] = 1
    data[extent[1] : extent[3] + 1, : (extent[2] + 1) % width] = 1

    polygon = polygonise_region(data)
    expected = sgeom.box(*extent)
    assert polygon.equals(expected)


@pytest.mark.parametrize(
    "hole",
    [(7, 3, 9, 5), (9, 3, 10, 5), (6, 3, 7, 5)],
)
def test_dateline_hole(hole):
    "Polygon crossing the dateline and containing a hole"
    width, height = 360, 180
    dims = (width, height, width, height)
    shell = (width * 5 // 8, height * 2 // 8, width * 11 // 8, height * 6 // 8)
    hole = [dim * t // 8 for dim, t in zip(dims, hole)]
    data = np.zeros((height, width), dtype=bool)
    data[shell[1] : shell[3] + 1, shell[0] % width :] = 1
    data[shell[1] : shell[3] + 1, : (shell[2] + 1) % width] = 1
    if hole[0] % width > hole[2] % width:
        # Hole is split across the edge: set in two steps
        data[hole[1] : hole[3] + 1, hole[0] % width :] = 0
        data[hole[1] : hole[3] + 1, : (hole[2] + 1) % width] = 0
    else:
        data[hole[1] : hole[3] + 1, hole[0] % width : (hole[2] + 1) % width] = 0

    polygon = polygonise_region(data)
    expected = sgeom.box(*shell) - sgeom.box(*hole)
    # Can't use strict equality (not in a general way anyway) because the
    # polygonisation algorithm can and will slightly cut concave corners, ie
    # the corners of the hole.  Instead check >= 99% overlap.
    assert (polygon - expected).area / expected.area <= 0.01
    assert (expected - polygon).area / expected.area <= 0.01


def test_north_pole():
    "Simple polygon near the north pole"
    width, height = 360, 180
    extent = (width * 3 // 8, height * 6 // 8, width * 5 // 8, height - 1)
    data = np.zeros((height, width), dtype=bool)
    # Use closed ranges instead of half-open, to match shapely box
    data[extent[1] : extent[3] + 1, extent[0] % width : (extent[2] + 1) % width] = 1

    polygon = polygonise_region(data)
    expected = sgeom.box(*extent)
    assert polygon.equals(expected)
