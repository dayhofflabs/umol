import pytest

from umol import SmilesIoConfig, SmilesSyntaxFlags


@pytest.mark.parametrize(
    ("bits", "expected"),
    [
        (0, SmilesSyntaxFlags.OPENSMILES),
        (2, SmilesSyntaxFlags.EXTENDED_AROMATICS),
        (4, SmilesSyntaxFlags.EXTENDED_BONDS),
        (6, SmilesSyntaxFlags.LENIENT),
    ],
)
def test_smiles_syntax_flags(bits, expected):
    assert SmilesSyntaxFlags(bits) == expected
    assert SmilesSyntaxFlags(bits).bits == bits


def test_smiles_syntax_flags_or():
    assert (
        SmilesSyntaxFlags.EXTENDED_AROMATICS | SmilesSyntaxFlags.EXTENDED_BONDS
    ) == SmilesSyntaxFlags.LENIENT


@pytest.mark.parametrize("bits", [1, 8, 1024, 1032, 1 << 31])
def test_smiles_syntax_flags_error(bits):
    with pytest.raises(
        ValueError,
        match=rf"^unknown SMILES syntax flag bits: {bits}$",
    ):
        SmilesSyntaxFlags(bits)


@pytest.mark.parametrize(
    "name",
    ["CHEMAXON_EXTENSIONS", "SKIP_UNKNOWN_CHEMAXON_TAGS", "CHEMAXON"],
)
def test_smiles_syntax_flags_cx(name):
    assert not hasattr(SmilesSyntaxFlags, name)


@pytest.mark.parametrize(
    ("flags", "expected"),
    [
        (SmilesSyntaxFlags.OPENSMILES, "SmilesSyntaxFlags.OPENSMILES"),
        (
            SmilesSyntaxFlags.EXTENDED_AROMATICS,
            "SmilesSyntaxFlags.EXTENDED_AROMATICS",
        ),
        (SmilesSyntaxFlags.EXTENDED_BONDS, "SmilesSyntaxFlags.EXTENDED_BONDS"),
        (SmilesSyntaxFlags.LENIENT, "SmilesSyntaxFlags.LENIENT"),
    ],
)
def test_smiles_syntax_flags_repr(flags, expected):
    assert repr(flags) == expected


def test_smiles_syntax_flags_immutable():
    with pytest.raises(AttributeError):
        SmilesSyntaxFlags.OPENSMILES.bits = 2


def test_smiles_io_config():
    extended_syntax = (
        SmilesSyntaxFlags.EXTENDED_AROMATICS | SmilesSyntaxFlags.EXTENDED_BONDS
    )
    config = SmilesIoConfig.with_parse_flags(extended_syntax)

    assert SmilesIoConfig.opensmiles().parse_flags == SmilesSyntaxFlags.OPENSMILES
    assert SmilesIoConfig.lenient().parse_flags == SmilesSyntaxFlags.LENIENT
    assert config.parse_flags == extended_syntax
    assert repr(config) == "SmilesIoConfig.lenient()"
    assert config == SmilesIoConfig.with_parse_flags(SmilesSyntaxFlags(6))
    assert config == SmilesIoConfig.lenient()
    assert config != SmilesIoConfig.opensmiles()
    assert not hasattr(config, "lint_flags")
    assert not hasattr(config, "lint_config")
    assert not hasattr(config, "chemistry_model")


def test_smiles_io_config_immutable():
    with pytest.raises(AttributeError):
        SmilesIoConfig.opensmiles().parse_flags = SmilesSyntaxFlags.LENIENT
