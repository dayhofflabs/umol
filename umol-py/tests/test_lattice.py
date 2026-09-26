import pytest

from umol import AtomConstraintForm, NoJoinError, NumForm, RingMembershipForm, RingScope


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (NumForm.Lit(1), NumForm.Lit(2), NumForm.LitSet({1, 2})),
        (
            AtomConstraintForm.Degree(NumForm.Lit(1)),
            AtomConstraintForm.Degree(NumForm.Lit(2)),
            AtomConstraintForm.Degree(NumForm.LitSet({1, 2})),
        ),
    ],
    ids=["numeric", "constraint"],
)
def test_lattice_join(left, right, expected):
    assert left.join(right) == expected
    assert right.join(left) == expected


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (
            AtomConstraintForm.Degree(NumForm.Lit(1)),
            AtomConstraintForm.Valence(NumForm.Lit(1)),
        ),
        (
            RingMembershipForm(RingScope.All(), 1),
            RingMembershipForm(RingScope.Size(6), 1),
        ),
    ],
    ids=["different_kinds", "different_scopes"],
)
def test_lattice_join_error(left, right):
    with pytest.raises(
        NoJoinError, match="^no join: elements have no least upper bound$"
    ):
        left.join(right)
    with pytest.raises(
        NoJoinError, match="^no join: elements have no least upper bound$"
    ):
        right.join(left)


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (NumForm.LitSet({1, 2}), NumForm.Lit(1), NumForm.Lit(1)),
        (NumForm.Lit(1), NumForm.Lit(2), None),
        (
            AtomConstraintForm.Degree(NumForm.Lit(1)),
            AtomConstraintForm.Degree(NumForm.Lit(2)),
            None,
        ),
    ],
    ids=["intersection", "numeric_bottom", "constraint_bottom"],
)
def test_lattice_meet(left, right, expected):
    assert left.meet(right) == expected
    assert right.meet(left) == expected
