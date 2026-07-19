import inspect

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


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (SmilesIoConfig.opensmiles(), SmilesSyntaxFlags.OPENSMILES),
        (SmilesIoConfig.lenient(), SmilesSyntaxFlags.LENIENT),
        (
            SmilesIoConfig.with_syntax_flags(
                syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS
            ),
            SmilesSyntaxFlags.EXTENDED_AROMATICS,
        ),
    ],
)
def test_smiles_io_config_syntax_flags(config, expected):
    assert config.syntax_flags == expected
    assert config.syntax_flags is not config.syntax_flags


def test_smiles_io_config_with_syntax_flags():
    syntax_flags = (
        SmilesSyntaxFlags.EXTENDED_AROMATICS | SmilesSyntaxFlags.EXTENDED_BONDS
    )
    config = SmilesIoConfig.with_syntax_flags(syntax_flags=syntax_flags)

    assert config == SmilesIoConfig.with_syntax_flags(
        syntax_flags=SmilesSyntaxFlags(6)
    )
    assert config == SmilesIoConfig.lenient()
    assert config != SmilesIoConfig.opensmiles()


def test_smiles_io_config_with_syntax_flags_error():
    assert str(inspect.signature(SmilesIoConfig.with_syntax_flags)) == "(*, syntax_flags)"
    with pytest.raises(
        TypeError,
        match=r"^SmilesIoConfig\.with_syntax_flags\(\) takes 0 positional "
        r"arguments but 1 was given$",
    ):
        SmilesIoConfig.with_syntax_flags(SmilesSyntaxFlags.LENIENT)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (SmilesIoConfig.opensmiles(), "SmilesIoConfig.opensmiles()"),
        (SmilesIoConfig.lenient(), "SmilesIoConfig.lenient()"),
        (
            SmilesIoConfig.with_syntax_flags(
                syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS
            ),
            "SmilesIoConfig.with_syntax_flags("
            "syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS)",
        ),
    ],
)
def test_smiles_io_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize("name", ["chemaxon", "with_parse_flags"])
def test_smiles_io_config_surface(name):
    assert not hasattr(SmilesIoConfig, name)


@pytest.mark.parametrize(
    "name",
    ["parse_flags", "lint_flags", "lint_config", "chemistry_model"],
)
def test_smiles_io_config_value_surface(name):
    assert not hasattr(SmilesIoConfig.opensmiles(), name)


def test_smiles_io_config_immutable():
    with pytest.raises(AttributeError):
        SmilesIoConfig.opensmiles().syntax_flags = SmilesSyntaxFlags.LENIENT
