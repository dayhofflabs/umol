import pytest
from umol import StereoInconsistency


@pytest.mark.parametrize(
    "value,expected,fields",
    [
        (
            StereoInconsistency.TetrahedralStereoFailure(atom=2),
            "StereoInconsistency.TetrahedralStereoFailure(atom=2)",
            {"atom": 2},
        ),
        (
            StereoInconsistency.StereoAtomFailure(stereo_atom=1),
            "StereoInconsistency.StereoAtomFailure(stereo_atom=1)",
            {"stereo_atom": 1},
        ),
        (
            StereoInconsistency.TetrahedralStereoMismatch(atom=3, stereo_atom=1),
            "StereoInconsistency.TetrahedralStereoMismatch(atom=3, stereo_atom=1)",
            {"atom": 3, "stereo_atom": 1},
        ),
        (
            StereoInconsistency.CisTransStereoFailure(bond=2),
            "StereoInconsistency.CisTransStereoFailure(bond=2)",
            {"bond": 2},
        ),
        (
            StereoInconsistency.StereoBondFailure(stereo_bond=1),
            "StereoInconsistency.StereoBondFailure(stereo_bond=1)",
            {"stereo_bond": 1},
        ),
        (
            StereoInconsistency.CisTransStereoMismatch(bond=3, stereo_bond=1),
            "StereoInconsistency.CisTransStereoMismatch(bond=3, stereo_bond=1)",
            {"bond": 3, "stereo_bond": 1},
        ),
    ],
)
def test_stereo_inconsistency(value, expected, fields):
    assert repr(value) == expected
    assert eval(repr(value), globals()) == value
    assert hash(value) == hash(eval(repr(value), globals()))
    for name, field in fields.items():
        assert getattr(value, name) == field
        with pytest.raises(AttributeError):
            setattr(value, name, field + 1)


def test_stereo_inconsistency_order():
    first = StereoInconsistency.TetrahedralStereoFailure(atom=0)
    second = StereoInconsistency.TetrahedralStereoFailure(atom=2)
    third = StereoInconsistency.StereoAtomFailure(stereo_atom=0)
    assert sorted([third, second, first]) == [first, second, third]
    assert first <= second
    assert third > second
    assert third >= StereoInconsistency.StereoAtomFailure(stereo_atom=0)
    assert first != 0
