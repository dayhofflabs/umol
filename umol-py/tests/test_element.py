import pytest

from umol import Element


def test_element_from_symbol():
    carbon = Element("C")
    assert carbon.symbol == "C"
    assert carbon.atomic_number == 6


def test_element_from_atomic_number():
    carbon = Element.from_atomic_number(6)
    assert carbon.symbol == "C"
    assert carbon.atomic_number == 6


def test_element_invalid_symbol():
    with pytest.raises(ValueError):
        Element("X")


def test_element_invalid_atomic_number():
    with pytest.raises(ValueError):
        Element.from_atomic_number(0)


def test_element_eq_hash():
    assert Element("C") == Element.from_atomic_number(6)
    assert Element("C") != Element("N")
    assert hash(Element("C")) == hash(Element.from_atomic_number(6))
    assert len({Element("C"), Element.from_atomic_number(6)}) == 1


def test_element_repr():
    assert repr(Element("C")) == "Element('C')"
