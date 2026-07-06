"""umol — molecular representation (Python bindings over the Rust core)."""

from importlib.metadata import PackageNotFoundError, version

from ._native import (
    AtomAst,
    AtomId,
    AtomView,
    AtomViews,
    Element,
    ElementAst,
    IsotopeMassAst,
    MemOp,
    MoleculeAst,
    RelOp,
    SpinStateAst,
    ValueAst,
    ValuePredicate,
    ValueTerm,
)

try:
    __version__ = version("umol")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "AtomAst",
    "AtomId",
    "AtomView",
    "AtomViews",
    "Element",
    "ElementAst",
    "IsotopeMassAst",
    "MemOp",
    "MoleculeAst",
    "RelOp",
    "SpinStateAst",
    "ValueAst",
    "ValuePredicate",
    "ValueTerm",
    "__version__",
]
