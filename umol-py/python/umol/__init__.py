"""umol — molecular representation (Python bindings over the Rust core)."""

from importlib.metadata import PackageNotFoundError, version

from . import _native

try:
    __version__ = version("umol")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = ["_native", "__version__"]
