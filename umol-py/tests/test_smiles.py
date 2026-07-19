import pytest

from umol import SmilesIoConfig, SmilesParseFlags


def test_smiles_parse_flags():
    extended_syntax = (
        SmilesParseFlags.EXTENDED_AROMATICS | SmilesParseFlags.EXTENDED_BONDS
    )

    assert SmilesParseFlags(0) == SmilesParseFlags.OPENSMILES
    assert SmilesParseFlags(2) == SmilesParseFlags.EXTENDED_AROMATICS
    assert SmilesParseFlags.CHEMAXON == SmilesParseFlags.CHEMAXON_EXTENSIONS
    assert extended_syntax.bits == 6
    assert repr(extended_syntax) == (
        "SmilesParseFlags(EXTENDED_AROMATICS | EXTENDED_BONDS)"
    )
    assert extended_syntax != SmilesParseFlags.LENIENT


def test_smiles_parse_flags_error():
    with pytest.raises(ValueError, match="^unknown SMILES parse flag bits: 1$"):
        SmilesParseFlags(1)


def test_smiles_parse_flags_immutable():
    with pytest.raises(AttributeError):
        SmilesParseFlags.OPENSMILES.bits = 2


def test_smiles_io_config():
    extended_syntax = (
        SmilesParseFlags.EXTENDED_AROMATICS | SmilesParseFlags.EXTENDED_BONDS
    )
    config = SmilesIoConfig.with_parse_flags(extended_syntax)

    assert SmilesIoConfig.opensmiles().parse_flags == SmilesParseFlags.OPENSMILES
    assert SmilesIoConfig.lenient().parse_flags == SmilesParseFlags.LENIENT
    assert SmilesIoConfig.chemaxon().parse_flags == SmilesParseFlags.CHEMAXON
    assert config.parse_flags == extended_syntax
    assert repr(config) == (
        "SmilesIoConfig.with_parse_flags("
        "SmilesParseFlags(EXTENDED_AROMATICS | EXTENDED_BONDS))"
    )
    assert config == SmilesIoConfig.with_parse_flags(SmilesParseFlags(6))
    assert config != SmilesIoConfig.opensmiles()
    assert not hasattr(config, "lint_flags")
    assert not hasattr(config, "lint_config")
    assert not hasattr(config, "chemistry_model")


def test_smiles_io_config_immutable():
    with pytest.raises(AttributeError):
        SmilesIoConfig.opensmiles().parse_flags = SmilesParseFlags.LENIENT
