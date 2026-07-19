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
    from umol import AtomTypeRegistry, ValenceEntry, ValenceModel, ValenceTable

    assert umol.AtomTypeRegistry is AtomTypeRegistry
    assert umol.ValenceEntry is ValenceEntry
    assert umol.ValenceModel is ValenceModel
    assert umol.ValenceTable is ValenceTable


def test_element_scope_import():
    from umol import ElementScope

    assert umol.ElementScope is ElementScope


def test_ring_limits_import():
    from umol import RingLimits

    assert umol.RingLimits is RingLimits


def test_aromaticity_model_import():
    from umol import AromaticityModel

    assert umol.AromaticityModel is AromaticityModel


def test_inconsistency_policy_import():
    from umol import InconsistencyPolicy

    assert umol.InconsistencyPolicy is InconsistencyPolicy


def test_stereo_kind_model_import():
    from umol import StereoKindModel

    assert umol.StereoKindModel is StereoKindModel


def test_stereo_model_import():
    from umol import StereoModel

    assert umol.StereoModel is StereoModel


def test_chemistry_model_import():
    from umol import ChemistryModel

    assert umol.ChemistryModel is ChemistryModel
