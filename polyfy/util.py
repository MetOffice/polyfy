import iris.analysis
import numpy as np
import shapely.affinity
import shapely.geometry as sgeom
from iris.cube import Cube


def wrap_longitudes(longitudes: np.ndarray, x0: float = -180) -> np.ndarray:
    """
    Wrap longitudes modulo 360 degrees

    Arguments:
        longitudes: Array of longitude values.
        x0: Start of the required range: longitudes will be wrapped into the
            half-open interval [x0, x0 + 360). Default -180.

    Returns:
        New array wrapped into the required range.
    """
    return (longitudes - x0) % 360 + x0


def unwrap_longitudes(longitudes: np.ndarray) -> np.ndarray:
    """
    Resolve discontiguities in longitudes arising as a result of wrapping

    Arguments:
        longitudes: Array of longitudes values.

    Returns:
        New array with discontiguities of magnitude approximately 360 degrees
        removed.
    """
    # Take all diffs, prepending a zero to preserve the first longitude
    diffs = np.diff(longitudes, prepend=0.0)

    # Wrap the diffs into the range [-180, 180], to convert the discontiguities
    # (~360° diffs) into more continuous (~0°) diffs
    diffs = wrap_longitudes(diffs, -180)

    # Sum the diffs to obtain contiguous longitudes
    return np.cumsum(diffs)


def wrap_polygon(polygon: sgeom.Polygon, x0: float = -180) -> sgeom.MultiPolygon:
    """
    Wrap a polygon modulo 360 degrees, splitting it where necessary

    Arguments:
        polygon: Polygon to split and wrap.
        x0: Start of the required range: longitudes will be wrapped into the
            half-open interval [x0, x0 + 360). Default -180.

    Returns:
        Multipolygon whose components all lie in the given interval
    """
    # Quick check for any intersection at all
    xmin, _, xmax, _ = polygon.bounds
    if xmin >= x0 and xmax < x0 + 360:
        return polygon

    # Split into the part "in range" and the rest
    part1 = sgeom.box(x0, -90, x0 + 360, 90).intersection(polygon)
    part2 = polygon - part1

    # Decide which way to shift the rest to keep it in range
    if xmin < x0:
        part2 = shapely.affinity.translate(part2, 360)
    else:
        part2 = shapely.affinity.translate(part2, -360)

    return part1.union(part2)


def flatten_cube(cube: Cube, method=iris.analysis.MAX) -> Cube:
    """
    Collapse a cube to 2 dimensions

    Arguments:
        cube: Cube containing at least latitude and longitude coordinates.
        method: Aggregation method to apply across each axis. Default maximum.

    Returns:
        2D cube, containing only latitude and longitude coordinates
    """
    flat_coords = [
        coord.name()
        for coord in cube.coords(dim_coords=True)
        if coord.name() not in ["latitude", "longitude"]
    ]
    return cube.collapsed(flat_coords, method)
