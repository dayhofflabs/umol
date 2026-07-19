import pytest

from umol import SmilesParseFlags


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
