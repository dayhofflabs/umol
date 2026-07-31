import pytest

from umol import UnpairedElectrons


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
