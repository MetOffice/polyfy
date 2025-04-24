from typing import Iterable

import numpy as np
import shapely.affinity
import shapely.geometry as sgeom
from iris.cube import Cube
from scipy import ndimage
from scipy.spatial import KDTree

from .util import flatten_cube


class Feature:
    def __init__(self, geometry, properties=None):
        if properties is None:
            properties = {}

        self.geometry = geometry
        self.properties = properties


def concave_hull(data: np.ndarray, k: int) -> sgeom.Polygon:
    """
    Find a concave hull of binary data

    Employs the algorithm outlined in Moreira, A. and Yasmina Santos, M.
    (2007), "Concave Hull: a k-nearest neighbours approach for the computation
    of the region occupied by a set of points" (DOI: 10.5220/0002080800610068).

    That is, starting from an arbitrary point (though minimal (y, x) is chosen
    for convenience), walk around the edge of the points, choosing the next
    point out of k nearest neighbours that creates the most concave angle.

    Arguments:
        data: 2D gridded data.
        k: Initial number of nearest neighbours to consider. Will be
            automatically incremented if a polygon cannot be found, up to a
            maximum of 2k.

    Raises:
        RuntimeError: if a polygon cannot be found. This is most likely to
            happen if the data contains too many collinear points.

    Returns:
        A single non-convex polygon covering nonzero data values. It is assumed
        that there is one connected region; if this is not the case then a
        single polygon is still returned, but it may not cover all the data.
        The polygon be empty if there are too few nonzero values.
    """
    pi = np.pi
    tau = 2 * pi

    # Get indices, in (x, y) order rather than (y, x)
    coords = np.stack(np.nonzero(data)[::-1], axis=1)
    n = coords.shape[0]
    if n < 3:
        # Not enough vertices to define a polygon
        return sgeom.Polygon()

    # Clamp k from 4 to n, 4 being the current point + 3 neighbours
    k = min(max(k, 4), n)
    # Set a maximum k to try before giving up
    k_max = min(k * 2, n)

    # Get a starting point, which will have minimal y (and minimal x in case of
    # a tie)
    point = first_point = coords[0]
    # Maintain current angle of travel, which should therefore start pointing
    # right, as we will walk anticlockwise from the bottom
    angle = initial_angle = 0
    # Maintain hull so far to check for self-intersections
    hull = sgeom.LineString()
    # Also maintain a list of the coordinates themselves, because converting a
    # multi line string to a polygon is surprisingly awkward
    points = [first_point]

    tree = KDTree(coords)
    while hull.is_empty or not np.all(first_point == point):
        can_close = False

        # Find k nearest neighbours, with associated distances and angles
        distances, indices = tree.query(point, k)
        neighbours = coords[indices]
        diffs = neighbours - point
        angles = np.arctan2(diffs[:, 1], diffs[:, 0])

        # Sort candidate neighours:
        # - First by new angle, in range [angle-pi, angle+pi] to prioritise
        #   sharpest right-hand turn
        # - Second by negative distance, to prioritise more distant points
        #   making the same angle
        angles = (angles - (angle - pi)) % tau + (angle - pi)
        candidates = sorted(zip(angles, -distances, neighbours))
        for new_angle, distance, new_point in candidates:
            if np.isclose(distance, 0):
                # Repeat of the previous point
                continue
            if np.all(new_point == first_point):
                can_close = True
            new_line = sgeom.LineString([point, new_point])
            if not hull.is_empty:
                if not (hull.intersection(new_line) - hull.boundary).is_empty:
                    # Invalid: intersects the hull so far in a point that is
                    # not at either end of the hull (its boundary)
                    continue

            # Add the new point
            hull = hull.union(new_line)
            point = new_point
            if np.isclose(new_angle, angle) and len(points) > 1:
                # Angle is not changing (and not because we only just started),
                # so replace the latest point instead of adding lots of
                # collinear points
                points[-1] = point
            else:
                points.append(point)
                angle = new_angle
            break

        else:
            # No valid points found (`break` not encountered)

            if can_close:
                # Done whether the resulting hull is valid or not.  Should
                # only occur when the shape was in fact small and/or degenerate
                # enough that we would rather an empty hull than an error.
                break

            if k >= k_max:
                # Bail out if tried too many times
                raise RuntimeError("unable to determine a concave hull")

            # Otherwise reset and try with a larger k
            k += 1
            hull = sgeom.LineString()
            point = first_point
            points = [first_point]
            angle = initial_angle

    if len(points) < 3:
        # Although there were at least 3 points available, they must have
        # described a degenerate shape
        return sgeom.Polygon()

    # Convert multi line string to polygon
    hull = sgeom.Polygon(points)

    if not hull.is_valid:
        # Should only happen if forcibly closing the loop with `can_close`:
        # validity in all other cases is incrementally enforced
        return sgeom.Polygon()

    return hull


def tidy_data(
    data: np.ndarray,
    sigma: float = 1.4,
    box: int = 3,
    threshold: float = 0.1,
    scale: int = 2,
    **kwargs,
) -> np.ndarray:
    """
    Smooth binary (thresholded) data

    Carries out the following steps:

    - Gaussian (default sigma = 1.4)
    - Box (default size = 3)
    - Threshold (default 10%)
    - Downscale (default 2x)

    Arguments:
        data: 2D boolean array.
        sigma: Gaussian blur parameter.
        box: Box filter parameter.
        threshold: Binarisation threshold, in the range 0-1.
        scale: Downscale factor.

    Returns:
        A new 2D array with filters and scaling applied.
    """
    # Low-pass filter to remove noise: Gaussian blur followed by box filter
    wrap_mode = ["constant", "wrap"]
    data = data * 100
    if sigma:
        data = ndimage.gaussian_filter(data, sigma, mode=wrap_mode)
    if box:
        data = ndimage.uniform_filter(data, box, mode=wrap_mode)

    # Apply a threshold to re-binarise the data
    data = data >= threshold * 100

    # Downscale to speed up subsequent processing.  Note that this is only
    # appropriate following the removal of noise.
    if scale > 1:
        data = data[::scale, ::scale]

    return data


def polygonise_region(data: np.ndarray, k: int = 11, **kwargs) -> sgeom.Polygon:
    """
    Convert gridded data to a polygon

    Arguments:
        data: 2D gridded data. It is assumed that appropriate filters and
            thresholds have already been applied. It is also assumed that the
            data represents a single connected region.
        k: Initial number of nearest neighbours to consider when finding a
            concave hull.

    Returns:
        Polygon covering the data.
    """
    # Apply Laplacian filter to detect edges, both speeding up polygonisation
    # and making it easier to identify interior edges.  The filter on its own
    # highlights both the "filled side" of an edge and the "non-filled side" of
    # an edge, so intersect with the original area to get only the filled side.
    wrap_mode = ["constant", "wrap"]
    data = ndimage.laplace(data, mode=wrap_mode) & data

    # Detect where wrapping has occurred
    first_empty = np.argmin(np.any(data, axis=0))
    if first_empty > 0:
        # This data is (most likely) split, so wrap the left part to the right.
        # May genuinely only touch the left edge without crossing it, in which
        # case this operation is harmless.
        data = np.concatenate([data, data[:, :first_empty]], axis=1)
        data[:, :first_empty] = 0

    # Label connected components.  Labelling starts from the corner, so the
    # first component will be the exterior boundary, and subsequent components
    # will be interiors
    data, n_components = ndimage.label(data, np.ones((3, 3)))

    # Convert exterior to a polygon
    polygon = concave_hull(data == 1, k=k)
    # and interiors to holes
    for i in range(2, n_components + 1):
        polygon -= concave_hull(data == i, k=k)

    # Simplify shape, to reduce the number of vertices and reduce how
    # noticeably pixellated it is (which manifests as zig-zags).  Tolerances
    # less than 1 leave a significant amount of zig-zagging, while tolerances
    # greater than 1 carry too much risk of cutting off significant detail.
    polygon = polygon.simplify(1)

    return polygon


def find_regions(data: np.ndarray) -> Iterable[np.ndarray]:
    """
    Yield each connected component of binary data

    Arguments:
        data: 2D gridded data.

    Yields:
        New 2D arrays, each representing a single connected component.
    """
    # Detect connected components.  There is no option to wrap around edges, so
    # we concatenate two copies side-by-side then overlay them.
    width = data.shape[1]
    data, n_components = ndimage.label(
        np.concatenate([data, data], axis=1),
        np.ones((3, 3)),
    )
    data = np.maximum(data[:, :width], data[:, width:])

    for i in range(1, n_components + 1):
        component = data == i
        if not np.any(component):
            # Was a duplicate introduced by our handling of wrapping
            continue

        yield component


def find_objects(cube: Cube, thresholds: dict, **kwargs) -> Iterable[Feature]:
    """
    Find polygons describing where thresholds are exceeded

    Arguments:
        cube: 3D field (Z, Y, X) of gridded data.
        thresholds: Mapping of threshold names to threshold values.

    Yields:
        One feature per polygon identified at each of the requested thresholds.
    """
    # Get coordinate info
    xcoord = cube.coord(axis="x", dim_coords=True)
    ycoord = cube.coord(axis="y", dim_coords=True)
    zcoord = cube.coord(axis="z")
    if cube.coord_dims(xcoord)[0] < cube.coord_dims(ycoord)[0]:
        raise RuntimeError("expected shape (y, x)")

    # Get transformation from grid indices to coordinate space
    scale = kwargs.pop("scale", 2)
    dx = np.diff(xcoord.points).mean()
    dy = np.diff(ycoord.points).mean()
    x0 = xcoord.points[0]
    y0 = ycoord.points[0]
    tf = [dx * scale, 0, 0, dy * scale, x0, y0]

    flat = flatten_cube(cube)
    for label, threshold in thresholds.items():
        # Generate fields of base heights and top heights for this threshold
        data = cube.data >= threshold
        bases = np.ma.masked_invalid(
            np.where(data, zcoord.bounds[:, 0, np.newaxis, np.newaxis], np.nan)
        ).min(axis=0)
        tops = np.ma.masked_invalid(
            np.where(data, zcoord.bounds[:, 1, np.newaxis, np.newaxis], np.nan)
        ).max(axis=0)

        # Scale height fields to match the data
        for i in range(scale):
            for j in range(scale):
                bases[::scale, ::scale] = np.minimum(
                    bases[::scale, ::scale], bases[i::scale, j::scale]
                )
                tops[::scale, ::scale] = np.maximum(
                    tops[::scale, ::scale], tops[i::scale, j::scale]
                )
        bases = bases[::scale, ::scale]
        tops = tops[::scale, ::scale]

        # Remove noise from 2D max field
        data = tidy_data(flat.data >= threshold, **kwargs)

        # Loop through connected components
        for region in find_regions(data):
            # Attempt to convert to a polygon. Very small areas may yield
            # empty polygons, to be skipped.
            polygon = polygonise_region(region, **kwargs)
            if polygon.is_empty:
                continue

            # Convert indices to latitude/longitude
            polygon = shapely.affinity.affine_transform(polygon, tf)

            # Assign additional information
            properties = {
                "severity": label,
                "base": int(np.ma.filled(np.min(bases[region]), 0)),
                "top": int(np.ma.filled(np.max(tops[region]), 600)),
            }

            yield Feature(polygon, properties)
