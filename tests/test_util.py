"""
Tests for `polyfy.util`
"""

import shapely.geometry as sgeom

from polyfy.util import wrap_polygon


def test_split_left():
    geom = sgeom.box(-181, 0, -179, 1)
    result = wrap_polygon(geom)
    expected = sgeom.box(179, 0, 180, 1).union(sgeom.box(-180, 0, -179, 1))
    print("input:", geom)
    print("output:", result)
    print("expected:", expected)
    assert result.equals(expected)


def test_split_right():
    geom = sgeom.box(179, 0, 181, 1)
    result = wrap_polygon(geom)
    expected = sgeom.box(179, 0, 180, 1).union(sgeom.box(-180, 0, -179, 1))
    print("input:", geom)
    print("output:", result)
    print("expected:", expected)
    assert result.equals(expected)


def test_no_split():
    geom = sgeom.box(0, 0, 1, 1)
    result = wrap_polygon(geom)
    expected = geom
    assert result.equals(expected)
