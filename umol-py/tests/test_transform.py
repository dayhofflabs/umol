import pytest

from umol import ConnectedComponentsAlgorithm, KekulizeConfig, MaximumMatchingAlgorithm


@pytest.mark.parametrize(
    "kwargs,expected_algorithm",
    [
        pytest.param({}, MaximumMatchingAlgorithm.Edmonds(), id="default"),
        pytest.param(
            {"algorithm": MaximumMatchingAlgorithm.Edmonds()},
            MaximumMatchingAlgorithm.Edmonds(),
            id="edmonds",
        ),
        pytest.param(
            {"algorithm": MaximumMatchingAlgorithm.HopcroftKarp()},
            MaximumMatchingAlgorithm.HopcroftKarp(),
            id="hopcroft-karp",
        ),
    ],
)
def test_kekulize_config_new(kwargs, expected_algorithm):
    config = KekulizeConfig(**kwargs)

    assert config.algorithm == expected_algorithm
    assert isinstance(config.algorithm, type(expected_algorithm))


@pytest.mark.parametrize(
    "args,kwargs",
    [
        pytest.param((MaximumMatchingAlgorithm.Edmonds(),), {}, id="positional"),
        pytest.param((), {"algorithm": "Edmonds"}, id="string"),
        pytest.param(
            (), {"algorithm": ConnectedComponentsAlgorithm.Bfs()}, id="wrong-enum"
        ),
    ],
)
def test_kekulize_config_new_error(args, kwargs):
    with pytest.raises(TypeError):
        KekulizeConfig(*args, **kwargs)


@pytest.mark.parametrize(
    "algorithm",
    [MaximumMatchingAlgorithm.Edmonds(), MaximumMatchingAlgorithm.HopcroftKarp()],
)
def test_kekulize_config_algorithm(algorithm):
    config = KekulizeConfig(algorithm=algorithm)

    with pytest.raises(AttributeError):
        config.algorithm = MaximumMatchingAlgorithm.HopcroftKarp()
    with pytest.raises(AttributeError):
        del config.algorithm
    assert config.algorithm == algorithm


@pytest.mark.parametrize(
    "algorithm,other,expected",
    [
        (MaximumMatchingAlgorithm.Edmonds(), MaximumMatchingAlgorithm.Edmonds(), True),
        (
            MaximumMatchingAlgorithm.HopcroftKarp(),
            MaximumMatchingAlgorithm.HopcroftKarp(),
            True,
        ),
        (
            MaximumMatchingAlgorithm.Edmonds(),
            MaximumMatchingAlgorithm.HopcroftKarp(),
            False,
        ),
    ],
)
def test_kekulize_config_eq(algorithm, other, expected):
    config = KekulizeConfig(algorithm=algorithm)
    other_config = KekulizeConfig(algorithm=other)

    assert (config == other_config) is expected
    assert (config != other_config) is not expected


@pytest.mark.parametrize(
    "algorithm,expected",
    [
        (
            MaximumMatchingAlgorithm.Edmonds(),
            "KekulizeConfig(algorithm=MaximumMatchingAlgorithm.Edmonds())",
        ),
        (
            MaximumMatchingAlgorithm.HopcroftKarp(),
            "KekulizeConfig(algorithm=MaximumMatchingAlgorithm.HopcroftKarp())",
        ),
    ],
)
def test_kekulize_config_repr(algorithm, expected):
    config = KekulizeConfig(algorithm=algorithm)

    assert repr(config) == expected
    assert eval(
        repr(config),
        {
            "KekulizeConfig": KekulizeConfig,
            "MaximumMatchingAlgorithm": MaximumMatchingAlgorithm,
        },
    ) == config
