import pytest

from umol import AromaticityResolveConfig, StereoResolveConfig


def test_aromaticity_resolve_config_default():
    config = AromaticityResolveConfig()

    assert config.delocalize_charge is True
    assert config.reset_aromatic_valence is False
    assert config == AromaticityResolveConfig()


@pytest.mark.parametrize(
    ("delocalize_charge", "reset_aromatic_valence"),
    [
        (True, False),
        (False, False),
        (True, True),
        (False, True),
    ],
)
def test_aromaticity_resolve_config_new(
    delocalize_charge, reset_aromatic_valence
):
    config = AromaticityResolveConfig(
        delocalize_charge=delocalize_charge,
        reset_aromatic_valence=reset_aromatic_valence,
    )

    assert config.delocalize_charge is delocalize_charge
    assert config.reset_aromatic_valence is reset_aromatic_valence


def test_aromaticity_resolve_config_new_error():
    with pytest.raises(TypeError):
        AromaticityResolveConfig(False, True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            AromaticityResolveConfig(),
            "AromaticityResolveConfig(delocalize_charge=True, "
            "reset_aromatic_valence=False)",
        ),
        (
            AromaticityResolveConfig(
                delocalize_charge=False,
                reset_aromatic_valence=True,
            ),
            "AromaticityResolveConfig(delocalize_charge=False, "
            "reset_aromatic_valence=True)",
        ),
    ],
)
def test_aromaticity_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("delocalize_charge", False),
        ("reset_aromatic_valence", True),
    ],
)
def test_aromaticity_resolve_config_mutation(field, value):
    config = AromaticityResolveConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


def test_stereo_resolve_config_default():
    config = StereoResolveConfig()

    assert config.reset_stereo_constraints is False
    assert config == StereoResolveConfig()


@pytest.mark.parametrize("reset_stereo_constraints", [False, True])
def test_stereo_resolve_config_new(reset_stereo_constraints):
    config = StereoResolveConfig(
        reset_stereo_constraints=reset_stereo_constraints
    )

    assert config.reset_stereo_constraints is reset_stereo_constraints


def test_stereo_resolve_config_new_error():
    with pytest.raises(TypeError):
        StereoResolveConfig(True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            StereoResolveConfig(),
            "StereoResolveConfig(reset_stereo_constraints=False)",
        ),
        (
            StereoResolveConfig(reset_stereo_constraints=True),
            "StereoResolveConfig(reset_stereo_constraints=True)",
        ),
    ],
)
def test_stereo_resolve_config_repr(config, expected):
    assert repr(config) == expected


def test_stereo_resolve_config_mutation():
    config = StereoResolveConfig()

    with pytest.raises(AttributeError):
        config.reset_stereo_constraints = True
