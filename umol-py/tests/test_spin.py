import pytest

from umol import SpinState, UnpairedElectrons


@pytest.mark.parametrize(
    ("count", "multiplicity"),
    [
        pytest.param(0, 1, id="closed_shell"),
        pytest.param(2, 3, id="open_shell"),
        pytest.param(-1, 0, id="physics_invalid"),
    ],
)
def test_unpaired_electrons(count, multiplicity):
    unpaired_electrons = UnpairedElectrons(count, multiplicity)

    assert unpaired_electrons.count == count
    assert unpaired_electrons.multiplicity == multiplicity


def test_unpaired_electrons_value_semantics():
    first = UnpairedElectrons(2, 3)
    same = UnpairedElectrons(2, 3)
    different = UnpairedElectrons(2, 1)

    assert first == same
    assert first != different
    assert hash(first) == hash(same)
    assert repr(first) == "UnpairedElectrons(count=2, multiplicity=3)"


@pytest.mark.parametrize("attribute", ["count", "multiplicity"])
def test_unpaired_electrons_assignment_error(attribute):
    unpaired_electrons = UnpairedElectrons(2, 3)

    with pytest.raises(AttributeError):
        setattr(unpaired_electrons, attribute, 1)


@pytest.mark.parametrize(
    ("unpaired_electrons", "multiplicity"),
    [
        pytest.param(0, 1, id="closed_shell"),
        pytest.param(1, 2, id="doublet"),
        pytest.param(2, 1, id="open_shell_singlet"),
        pytest.param(2, 3, id="triplet"),
    ],
)
def test_spin_state(unpaired_electrons, multiplicity):
    spin_state = SpinState(
        unpaired_electrons=unpaired_electrons,
        multiplicity=multiplicity,
    )

    assert spin_state.unpaired_electrons == unpaired_electrons
    assert spin_state.multiplicity == multiplicity


@pytest.mark.parametrize(
    ("unpaired_electrons", "multiplicity", "message"),
    [
        pytest.param(-1, 1, "unpaired electrons -1 out of range", id="negative_count"),
        pytest.param(256, 1, "unpaired electrons 256 out of range", id="large_count"),
        pytest.param(0, 0, "multiplicity 0 out of range", id="zero_multiplicity"),
        pytest.param(0, 256, "multiplicity 256 out of range", id="large_multiplicity"),
        pytest.param(
            2,
            2,
            "2 unpaired electrons, 2 multiplicity incompatible",
            id="incompatible",
        ),
    ],
)
def test_spin_state_error(unpaired_electrons, multiplicity, message):
    with pytest.raises(ValueError, match=message):
        SpinState(
            unpaired_electrons=unpaired_electrons,
            multiplicity=multiplicity,
        )


def test_spin_state_signature_error():
    with pytest.raises(TypeError):
        SpinState(2, 3)


def test_spin_state_value_semantics():
    first = SpinState(unpaired_electrons=2, multiplicity=3)
    same = SpinState(unpaired_electrons=2, multiplicity=3)
    different = SpinState(unpaired_electrons=2, multiplicity=1)

    assert first == same
    assert first != different
    assert hash(first) == hash(same)
    assert repr(first) == "SpinState(unpaired_electrons=2, multiplicity=3)"


@pytest.mark.parametrize("attribute", ["unpaired_electrons", "multiplicity"])
def test_spin_state_assignment_error(attribute):
    spin_state = SpinState(unpaired_electrons=2, multiplicity=3)

    with pytest.raises(AttributeError):
        setattr(spin_state, attribute, 1)
