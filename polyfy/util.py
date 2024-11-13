import iris.analysis
from iris.cube import Cube


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
