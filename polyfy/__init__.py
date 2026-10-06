from importlib.metadata import version, PackageNotFoundError

from . import io
from .creation import Feature, find_objects

__all__ = ["Feature", "io", "find_objects"]

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    pass
