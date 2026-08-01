import pytest

from umol import (
    AromaticityConfig,
    AromaticityInconsistencyPolicy,
    AromaticityResolveConfig,
    ResolveConfig,
    StereoInconsistencyPolicy,
    StereoResolveConfig,
)


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (
            AromaticityInconsistencyPolicy.Keep,
            AromaticityInconsistencyPolicy.Keep,
            True,
        ),
        (
            AromaticityInconsistencyPolicy.Error,
            AromaticityInconsistencyPolicy.Error,
            True,
        ),
        (
            AromaticityInconsistencyPolicy.Keep,
            AromaticityInconsistencyPolicy.Error,
            False,
        ),
    ],
)
def test_aromaticity_inconsistency_policy_equality(left, right, expected):
    assert (left == right) is expected


def test_aromaticity_inconsistency_policy_hash():
    policies = {
        AromaticityInconsistencyPolicy.Keep: "keep",
        AromaticityInconsistencyPolicy.Error: "error",
    }

    assert policies[AromaticityInconsistencyPolicy.Keep] == "keep"
    assert policies[AromaticityInconsistencyPolicy.Error] == "error"


@pytest.mark.parametrize(
    ("policy", "expected"),
    [
        (
            AromaticityInconsistencyPolicy.Keep,
            "AromaticityInconsistencyPolicy.Keep",
        ),
        (
            AromaticityInconsistencyPolicy.Error,
            "AromaticityInconsistencyPolicy.Error",
        ),
    ],
)
def test_aromaticity_inconsistency_policy_repr(policy, expected):
    assert repr(policy) == expected


@pytest.mark.parametrize(
    "policy",
    [
        AromaticityInconsistencyPolicy.Keep,
        AromaticityInconsistencyPolicy.Error,
    ],
)
def test_aromaticity_inconsistency_policy_mutation(policy):
    with pytest.raises(AttributeError):
        policy.value = "changed"


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (StereoInconsistencyPolicy.Keep, StereoInconsistencyPolicy.Keep, True),
        (StereoInconsistencyPolicy.Strip, StereoInconsistencyPolicy.Strip, True),
        (StereoInconsistencyPolicy.Error, StereoInconsistencyPolicy.Error, True),
        (StereoInconsistencyPolicy.Keep, StereoInconsistencyPolicy.Strip, False),
        (StereoInconsistencyPolicy.Strip, StereoInconsistencyPolicy.Error, False),
        (StereoInconsistencyPolicy.Error, StereoInconsistencyPolicy.Keep, False),
    ],
)
def test_stereo_inconsistency_policy_equality(left, right, expected):
    assert (left == right) is expected


def test_stereo_inconsistency_policy_hash():
    policies = {
        StereoInconsistencyPolicy.Keep: "keep",
        StereoInconsistencyPolicy.Strip: "strip",
        StereoInconsistencyPolicy.Error: "error",
    }

    assert policies[StereoInconsistencyPolicy.Keep] == "keep"
    assert policies[StereoInconsistencyPolicy.Strip] == "strip"
    assert policies[StereoInconsistencyPolicy.Error] == "error"


@pytest.mark.parametrize(
    ("policy", "expected"),
    [
        (StereoInconsistencyPolicy.Keep, "StereoInconsistencyPolicy.Keep"),
        (StereoInconsistencyPolicy.Strip, "StereoInconsistencyPolicy.Strip"),
        (StereoInconsistencyPolicy.Error, "StereoInconsistencyPolicy.Error"),
    ],
)
def test_stereo_inconsistency_policy_repr(policy, expected):
    assert repr(policy) == expected


@pytest.mark.parametrize(
    "policy",
    [
        StereoInconsistencyPolicy.Keep,
        StereoInconsistencyPolicy.Strip,
        StereoInconsistencyPolicy.Error,
    ],
)
def test_stereo_inconsistency_policy_mutation(policy):
    with pytest.raises(AttributeError):
        policy.value = "changed"


def test_aromaticity_resolve_config_default():
    config = AromaticityResolveConfig()

    assert config.perception == AromaticityConfig()
    assert config.inconsistency == AromaticityInconsistencyPolicy.Error
    assert config.reset_aromatic_valence is False
    assert config == AromaticityResolveConfig()


@pytest.mark.parametrize(
    ("inconsistency", "reset_aromatic_valence"),
    [
        (AromaticityInconsistencyPolicy.Keep, False),
        (AromaticityInconsistencyPolicy.Error, False),
        (AromaticityInconsistencyPolicy.Keep, True),
        (AromaticityInconsistencyPolicy.Error, True),
    ],
)
def test_aromaticity_resolve_config_new(
    inconsistency, reset_aromatic_valence
):
    perception = AromaticityConfig()
    config = AromaticityResolveConfig(
        perception=perception,
        inconsistency=inconsistency,
        reset_aromatic_valence=reset_aromatic_valence,
    )

    assert config.perception == perception
    assert config.perception is not perception
    assert config.inconsistency == inconsistency
    assert config.reset_aromatic_valence is reset_aromatic_valence


def test_aromaticity_resolve_config_new_error():
    with pytest.raises(TypeError):
        AromaticityResolveConfig(False, True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            AromaticityResolveConfig(),
            "AromaticityResolveConfig(perception=AromaticityConfig("
            "ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "inconsistency=AromaticityInconsistencyPolicy.Error, "
            "reset_aromatic_valence=False)",
        ),
        (
            AromaticityResolveConfig(
                inconsistency=AromaticityInconsistencyPolicy.Keep,
                reset_aromatic_valence=True,
            ),
            "AromaticityResolveConfig(perception=AromaticityConfig("
            "ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "inconsistency=AromaticityInconsistencyPolicy.Keep, "
            "reset_aromatic_valence=True)",
        ),
    ],
)
def test_aromaticity_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("perception", AromaticityConfig()),
        ("inconsistency", AromaticityInconsistencyPolicy.Keep),
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
    assert config.inconsistency == StereoInconsistencyPolicy.Error
    assert config == StereoResolveConfig()


@pytest.mark.parametrize(
    ("reset_stereo_constraints", "inconsistency"),
    [
        (False, StereoInconsistencyPolicy.Keep),
        (True, StereoInconsistencyPolicy.Strip),
        (False, StereoInconsistencyPolicy.Error),
    ],
)
def test_stereo_resolve_config_new(
    reset_stereo_constraints, inconsistency
):
    config = StereoResolveConfig(
        reset_stereo_constraints=reset_stereo_constraints,
        inconsistency=inconsistency,
    )

    assert config.reset_stereo_constraints is reset_stereo_constraints
    assert config.inconsistency == inconsistency


def test_stereo_resolve_config_new_error():
    with pytest.raises(TypeError):
        StereoResolveConfig(True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            StereoResolveConfig(),
            "StereoResolveConfig(reset_stereo_constraints=False, "
            "inconsistency=StereoInconsistencyPolicy.Error)",
        ),
        (
            StereoResolveConfig(
                reset_stereo_constraints=True,
                inconsistency=StereoInconsistencyPolicy.Strip,
            ),
            "StereoResolveConfig(reset_stereo_constraints=True, "
            "inconsistency=StereoInconsistencyPolicy.Strip)",
        ),
    ],
)
def test_stereo_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("reset_stereo_constraints", True),
        ("inconsistency", StereoInconsistencyPolicy.Keep),
    ],
)
def test_stereo_resolve_config_mutation(field, value):
    config = StereoResolveConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


def test_resolve_config_default():
    config = ResolveConfig.default()

    assert config.aromaticity == AromaticityResolveConfig()
    assert config.stereo == StereoResolveConfig()
    assert config == ResolveConfig.default()


@pytest.mark.parametrize(
    ("aromaticity", "stereo"),
    [
        (
            AromaticityResolveConfig(
                inconsistency=AromaticityInconsistencyPolicy.Keep,
                reset_aromatic_valence=True,
            ),
            StereoResolveConfig(),
        ),
        (
            AromaticityResolveConfig(),
            StereoResolveConfig(reset_stereo_constraints=True),
        ),
    ],
)
def test_resolve_config_new(aromaticity, stereo):
    config = ResolveConfig(aromaticity=aromaticity, stereo=stereo)

    assert config.aromaticity == aromaticity
    assert config.aromaticity is not aromaticity
    assert config.stereo == stereo
    assert config.stereo is not stereo


def test_resolve_config_new_error():
    with pytest.raises(TypeError):
        ResolveConfig(AromaticityResolveConfig(), StereoResolveConfig())


@pytest.mark.parametrize(
    "other",
    [
        ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                inconsistency=AromaticityInconsistencyPolicy.Keep,
                reset_aromatic_valence=False,
            ),
            stereo=StereoResolveConfig(),
        ),
        ResolveConfig(
            aromaticity=AromaticityResolveConfig(),
            stereo=StereoResolveConfig(
                reset_stereo_constraints=True,
                inconsistency=StereoInconsistencyPolicy.Strip,
            ),
        ),
    ],
)
def test_resolve_config_equality(other):
    assert ResolveConfig.default() != other


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (ResolveConfig.default(), "ResolveConfig.default()"),
        (
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    inconsistency=AromaticityInconsistencyPolicy.Keep,
                    reset_aromatic_valence=True,
                ),
                stereo=StereoResolveConfig(
                    reset_stereo_constraints=True,
                    inconsistency=StereoInconsistencyPolicy.Strip,
                ),
            ),
            "ResolveConfig(aromaticity=AromaticityResolveConfig(perception="
            "AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "inconsistency=AromaticityInconsistencyPolicy.Keep, "
            "reset_aromatic_valence=True), "
            "stereo=StereoResolveConfig(reset_stereo_constraints=True, "
            "inconsistency=StereoInconsistencyPolicy.Strip))",
        ),
    ],
)
def test_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("aromaticity", AromaticityResolveConfig()),
        ("stereo", StereoResolveConfig()),
    ],
)
def test_resolve_config_mutation(field, value):
    config = ResolveConfig.default()

    with pytest.raises(AttributeError):
        setattr(config, field, value)
