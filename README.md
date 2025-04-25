# Polyfy

Convert gridded data to polygons using a concave hull algorithm.


## Installation

Currently, this package can only be installed from source.

First, setup an appropriate environment, such as a conda environment or a python venv, then install the package and its dependencies.

### Conda

```bash
conda env create --file environment.yml
conda activate polyfy
pip install .
```

### Venv

```bash
python -m venv .venv
source .venv/bin/activate
pip install .
```


## Example usage

```python
import iris
import polyfy

data = iris.load("data.nc")
polygons = polyfy.find_objects(data)
polyfy.io.to_shapefile(polygons, "polygons.shp")
```

Other serialisation methods are available, namely `polyfy.io.to_geojson`.
[IWXXM] format via `to_iwxxm` is currently only available via a closed-source dependency.


## Algorithm overview

This code aims to generate prisms - polygons with the added information of a base height and a top height - that represent 3-dimensional gridded environmental data.
The general approach taken is to determine polygon shapes from a 2D field of maximum values, then determine which vertical levels each shape should cover.
Specifically, this is accomplished by the following steps:

1. [Process gridded data](#process-gridded-data) - 3D (z, y, x) fields are preprocessed in order to facilitate identifying distinct objects, including collapsing to 2D (y, x) fields.
1. [Find objects](#find-objects) - initial objects are found by a concave hull algorithm.
1. [Refine objects](#refine-objects) - vertices are removed from the objects in order to simplify them and remove artefacts arising from the grid structure.
1. [Find object characteristics](#find-object-characteristics) - top and base values are found for each object.

Once created, the polygon objects can be saved to a file and plotted.

Further details on some steps are expanded on in the following sections.


## Details

### Process gridded data

First, each 3D field is collapsed to a 2D field by taking the maximum over the vertical dimension.
Polygon objects will be identified from these 2D fields, with [height information added later](#find-object-characteristics).

Each threshold is then applied in turn, and the resulting binary data is smoothed through the use of low-pass filters.
Default filters and their parameters were chosen such that data is spread out by up to two grid cells, and that noisy areas containing few isolated grid cells are filtered away.
These are:

- Gaussian (sigma = 1.4)
- Box (size = 3x3)
- Binarise (10% threshold)

These can, however, be configured.

Further, having removed noise, the gridded data is then downscaled (2x) to speed up subsequent processing.


### Find objects

Objects are found from gridded data as a concave hull of the grid cells of interest.

In contrast to the [convex hull], the unique smallest convex shape that covers a given set of points, a concave hull relaxes the convex requirement, and consequently is not unique.

The algorithm employed to find a concave hull is a modification of a [convex hull algorithm]: finding the locally sharpest right-hand turn instead of the globally shallowest left-hand turn.
This is specifically carried out by choosing one pixel as a starting point, then iteratively identifing the *k* closest pixels and choosing the one that creates the most acute angle to be the next point.

Smaller values of *k* follow the boundary more precisely, which may take longer due to identifying a longer perimeter, and will take longer due to taking smaller steps along it.
Larger values of *k* require less processing time due to taking longer steps, but the allowed step length may be large enough to skip over entire strongly concave regions.
On balance, a value of 11 has been chosen to be the default, though this is configurable.


### Refine objects

The polygon resulting from the concave hull algorithm shows some obvious artefacts from the fact that all points were on a grid.
Further adjustments are therefore made to help obtain a more natural shape.

This is achieved by applying the Douglas-Peucker simplification algorithm (via [`shapely.simplify`]) with a tolerance of 1.
This tolerance was chosen as a balance of being high enough that the gridded artefacts become less noticeable, yet still small enough that the shapes are not changed too much.


### Find object characteristics

Having identified each object from the 2D fields, top and base heights are assigned according to the original 3D fields.

Considering all horizontal grid cells contributing to the object, the lowest and highest vertical levels containing values above the relevant threshold are identified.
The result is a collection of 3D prisms enclosing the gridded data.


[IWXXM]: https://en.wikipedia.org/wiki/IWXXM
[convex hull]: https://en.wikipedia.org/wiki/Convex_hull
[convex hull algorithm]: https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
[`shapely.simplify`]: https://shapely.readthedocs.io/en/latest/reference/shapely.simplify.html
