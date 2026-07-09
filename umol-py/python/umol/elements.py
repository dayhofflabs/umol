"""Terse, unquoted element access: ``E.H``, ``E.Cl``, ``E.As`` -> ``Element(...)``.

The only Python spelling that keeps the bare, unquoted symbol -- ``Element("H")``
and ``E["H"]`` need quotes, and ``Element(H)`` is a ``NameError``. Mirrors the Rust
``e!(H)`` macro.
"""

from ._native import Element


class _Elements:
    """Namespace whose attribute access resolves an element symbol to an ``Element``."""

    __slots__ = ()

    def __getattr__(self, symbol: str) -> Element:
        if symbol.startswith("_"):
            raise AttributeError(symbol)
        try:
            return Element(symbol)
        except ValueError as exc:
            raise AttributeError(f"no element with symbol {symbol!r}") from exc

    def __getitem__(self, symbol: str) -> Element:
        return Element(symbol)

    def __repr__(self) -> str:
        return "E"


E = _Elements()
