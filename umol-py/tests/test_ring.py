import pytest
import umol

from umol import (
    RelevantCycleEnumerationAlgorithm,
    RingConfig,
    SimpleCycleEnumerationAlgorithm,
)


def test_ring_exports():
    exports = {
        "RelevantCycleEnumerationAlgorithm": RelevantCycleEnumerationAlgorithm,
        "RingConfig": RingConfig,
        "SimpleCycleEnumerationAlgorithm": SimpleCycleEnumerationAlgorithm,
    }

    assert exports.keys() <= set(umol.__all__)
    assert {name: getattr(umol, name) for name in exports} == exports
    assert not hasattr(umol, "CycleEnumerationAlgorithm")


def test_ring_config_default():
    config = RingConfig()

    assert (
        config.simple_cycle_algorithm
        == SimpleCycleEnumerationAlgorithm.ReadTarjan()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert config == RingConfig()


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "simple_cycle_algorithm": (
                SimpleCycleEnumerationAlgorithm.ReadTarjan()
            ),
        },
        {
            "relevant_cycle_algorithm": (
                RelevantCycleEnumerationAlgorithm.Vismara()
            ),
        },
        {
            "simple_cycle_algorithm": (
                SimpleCycleEnumerationAlgorithm.ReadTarjan()
            ),
            "relevant_cycle_algorithm": (
                RelevantCycleEnumerationAlgorithm.Vismara()
            ),
        },
    ],
)
def test_ring_config_new(kwargs):
    config = RingConfig(**kwargs)

    assert (
        config.simple_cycle_algorithm
        == SimpleCycleEnumerationAlgorithm.ReadTarjan()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )


def test_ring_config_new_error():
    with pytest.raises(TypeError):
        RingConfig(
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
            RelevantCycleEnumerationAlgorithm.Vismara(),
        )


def test_ring_config_repr():
    assert repr(RingConfig()) == (
        "RingConfig("
        "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara())"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        (
            "simple_cycle_algorithm",
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
        ),
        (
            "relevant_cycle_algorithm",
            RelevantCycleEnumerationAlgorithm.Vismara(),
        ),
    ],
)
def test_ring_config_mutation(field, value):
    config = RingConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)
