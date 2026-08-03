import pytest

from umol import New


def test_new():
    handle = New(3)

    assert handle.index == 3
    assert handle == New(3)
    assert handle != New(4)
    assert repr(handle) == "New(3)"


def test_new_immutability():
    handle = New(3)

    with pytest.raises(AttributeError):
        handle.index = 4
    with pytest.raises(AttributeError):
        handle.extra = 4


@pytest.mark.parametrize("index", [-1, 2**100])
def test_new_error(index):
    with pytest.raises(OverflowError):
        New(index)
