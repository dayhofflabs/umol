import umol


def test_import_umol():
    assert isinstance(umol.__version__, str)
    assert umol.__version__ != "0.0.0"


def test_native_present():
    assert hasattr(umol, "_native")


def test_error_imports():
    from umol import (
        ContradictionError,
        InvalidStructureError,
        ModelConversionError,
        ParseError,
        UnderdeterminedError,
    )

    for error_type in (
        ContradictionError,
        InvalidStructureError,
        ModelConversionError,
        ParseError,
        UnderdeterminedError,
    ):
        error = error_type("diagnostic")
        assert getattr(umol, error_type.__name__) is error_type
        assert isinstance(error, Exception)
        assert str(error) == "diagnostic"


def test_smiles_syntax_flags_import():
    from umol import SmilesIoConfig, SmilesSyntaxFlags

    assert umol.SmilesIoConfig is SmilesIoConfig
    assert umol.SmilesSyntaxFlags is SmilesSyntaxFlags
    assert not hasattr(umol, "SmilesParseFlags")


def test_atom_type_registry_import():
    from umol import AtomTypeRegistry, ValenceEntry, ValenceTable

    assert umol.AtomTypeRegistry is AtomTypeRegistry
    assert umol.ValenceEntry is ValenceEntry
    assert umol.ValenceTable is ValenceTable
