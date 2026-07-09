"""umol — molecular representation (Python bindings over the Rust core)."""

from importlib.metadata import PackageNotFoundError, version

from ._native import (
    AromaticValenceAst,
    AtomAst,
    AtomConstraintAst,
    AtomConstraintKey,
    AtomConstraintsAst,
    AtomId,
    AtomView,
    AtomViews,
    Element,
    ElementAst,
    IsotopeMassAst,
    MemOp,
    MoleculeAst,
    MulticenterValenceAst,
    ParseError,
    Permutation,
    RelOp,
    RingMembershipAst,
    RingScope,
    SpinStateAst,
    StereoCosetAst,
    StereoTerm,
    TetrahedralStereoAst,
    ValueAst,
    ValuePredicate,
    ValueTerm,
)

try:
    __version__ = version("umol")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "AromaticValenceAst",
    "AtomAst",
    "AtomConstraintAst",
    "AtomConstraintKey",
    "AtomConstraintsAst",
    "AtomId",
    "AtomView",
    "AtomViews",
    "Element",
    "ElementAst",
    "IsotopeMassAst",
    "MemOp",
    "MoleculeAst",
    "MulticenterValenceAst",
    "ParseError",
    "Permutation",
    "RelOp",
    "RingMembershipAst",
    "RingScope",
    "SpinStateAst",
    "StereoCosetAst",
    "StereoTerm",
    "TetrahedralStereoAst",
    "ValueAst",
    "ValuePredicate",
    "ValueTerm",
    "__version__",
]
