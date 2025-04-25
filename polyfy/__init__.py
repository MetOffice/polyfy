from importlib.metadata import version, PackageNotFoundError

from . import io
from .creation import find_objects

__all__ = ["io", "find_objects"]

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    pass
