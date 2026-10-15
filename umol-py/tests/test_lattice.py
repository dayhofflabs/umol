import pytest

from umol import (
    AromaticSystemForm,
    AtomConstraintForm,
    AtomForm,
    BondForm,
    BooleanForm,
    ContradictionError,
    DativeBondForm,
    ElectronCountsForm,
    Element,
    ElementForm,
    MulticenterBondForm,
    NoJoinError,
    NoncovalentBondConstraintForm,
    NoncovalentBondConstraintsForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NoncovalentBondKindForm,
    NumForm,
    RingMembershipForm,
    RingScope,
    StereoAtomForm,
    StereoBondForm,
    StereoConfigurationForm,
    StereoCoset,
    StereoKind,
)


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


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (NumForm.LitSet({1}), NumForm.Lit(1)),
        (
            AtomForm(Element("C"), charge=NumForm.LitSet({0})),
            AtomForm(Element("C"), charge=0),
        ),
    ],
    ids=["numeric", "atom"],
)
def test_lattice_normalize(source, expected):
    assert source.normalize() == expected
    assert source.normalized_eq(expected)
    assert expected.normalized_eq(source)
    assert source != expected


@pytest.mark.parametrize(
    "source",
    [
        NumForm.LitSet(set()),
        AtomForm(Element("C"), charge=NumForm.LitSet(set())),
        BondForm(NumForm.LitSet(set())),
        StereoAtomForm(
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.LitSet(set())
            )
        ),
        StereoBondForm(
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.LitSet(set())
            )
        ),
    ],
    ids=["numeric", "atom", "bond", "stereo-atom", "stereo-bond"],
)
def test_lattice_normalize_error(source):
    before = repr(source)
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        source.normalize()
    assert repr(source) == before
    assert source.normalized_eq(source)


@pytest.mark.parametrize(
    ("pattern", "target", "disjoint", "ground", "field", "replacement"),
    [
        pytest.param(
            AtomForm(Element("C"), charge=NumForm.LitSet({0, 1})),
            AtomForm(Element("C"), charge=0),
            AtomForm(Element("C"), charge=2),
            False,
            "charge",
            2,
            id="atom",
        ),
        pytest.param(
            BondForm(NumForm.LitSet({1, 2})),
            BondForm(1),
            BondForm(3),
            False,
            "order",
            3,
            id="bond",
        ),
        pytest.param(
            DativeBondForm(NumForm.LitSet({1, 2})),
            DativeBondForm(1),
            DativeBondForm(3),
            True,
            "order",
            3,
            id="dative-bond",
        ),
        pytest.param(
            AromaticSystemForm([1, 1], charge=NumForm.LitSet({0, 1})),
            AromaticSystemForm([1, 1], charge=0),
            AromaticSystemForm([1, 1], charge=2),
            False,
            "charge",
            2,
            id="aromatic-system",
        ),
        pytest.param(
            MulticenterBondForm([1, 1], charge=NumForm.LitSet({0, 1})),
            MulticenterBondForm([1, 1], charge=0),
            MulticenterBondForm([1, 1], charge=2),
            False,
            "charge",
            2,
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondForm(
                NoncovalentBondKindForm.Undetermined(),
                constraints=NoncovalentBondConstraintsForm(
                    [
                        NoncovalentBondConstraintForm.Intramolecular(
                            BooleanForm.Lit(True)
                        )
                    ]
                ),
            ),
            NoncovalentBondForm(
                NoncovalentBondKind.HydrogenBond,
                constraints=NoncovalentBondConstraintsForm(
                    [
                        NoncovalentBondConstraintForm.Intramolecular(
                            BooleanForm.Lit(True)
                        )
                    ]
                ),
            ),
            NoncovalentBondForm(NoncovalentBondKind.Ionic),
            True,
            "kind",
            NoncovalentBondKind.Ionic,
            id="noncovalent-bond",
        ),
        pytest.param(
            StereoAtomForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Undetermined()
                )
            ),
            StereoAtomForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(0)
                )
            ),
            StereoAtomForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                )
            ),
            True,
            "configuration",
            StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(1)),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Undetermined()
                )
            ),
            StereoBondForm(
                StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(0))
            ),
            StereoBondForm(
                StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(1))
            ),
            True,
            "configuration",
            StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(1)),
            id="stereo-bond",
        ),
    ],
)
def test_entity_form_lattice(pattern, target, disjoint, ground, field, replacement):
    before_pattern = pattern.asdict()
    before_target = target.asdict()
    before_disjoint = disjoint.asdict()

    assert pattern.is_undetermined() is False
    assert pattern.is_ground() is False
    assert target.is_ground() is ground
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False
    assert pattern.is_compatible(target) is True
    assert target.is_compatible(disjoint) is False
    assert pattern.normalized_eq(pattern.copy()) is True
    assert pattern.normalized_eq(target) is False

    meet = pattern.meet(target)
    assert meet.asdict() == before_target
    assert target.meet(disjoint) is None
    join = pattern.join(target)
    assert join.asdict() == before_pattern
    normalized = target.normalize()
    assert normalized.asdict() == before_target

    for result in (meet, join, normalized):
        assert result is not pattern and result is not target
        assert result.readonly is False
        before_result = result.asdict()
        setattr(result, field, replacement)
        assert result.asdict() != before_result

    assert pattern.asdict() == before_pattern
    assert target.asdict() == before_target
    assert disjoint.asdict() == before_disjoint


@pytest.mark.parametrize(
    "source",
    [
        AtomForm(ElementForm.Undetermined()),
        BondForm(NumForm.Undetermined()),
        DativeBondForm(NumForm.Undetermined()),
        AromaticSystemForm(ElectronCountsForm.Undetermined()),
        MulticenterBondForm(ElectronCountsForm.Undetermined()),
        NoncovalentBondForm(NoncovalentBondKindForm.Undetermined()),
        StereoAtomForm(StereoConfigurationForm.Undetermined()),
        StereoBondForm(StereoConfigurationForm.Undetermined()),
    ],
    ids=[
        "atom",
        "bond",
        "dative-bond",
        "aromatic-system",
        "multicenter-bond",
        "noncovalent-bond",
        "stereo-atom",
        "stereo-bond",
    ],
)
def test_entity_form_is_undetermined(source):
    before = source.asdict()
    assert source.is_undetermined() is True
    assert source.is_ground() is False
    assert source.asdict() == before
