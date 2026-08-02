import pytest
from umol import (
    AromaticSystemAst,
    AromaticSystemConstraintAst,
    AromaticSystemConstraintsAst,
    AromaticValenceAst,
    AtomAst,
    AtomConstraintAst,
    AtomConstraintsAst,
    BondAst,
    BondConstraintAst,
    BondConstraintsAst,
    BooleanAst,
    CisTransStereoAst,
    ContradictionError,
    DativeBondAst,
    DativeBondConstraintAst,
    DativeBondConstraintsAst,
    ElectronCountsAst,
    Element,
    ElementAst,
    FluxionalityAst,
    IsotopeMassAst,
    LigandPermutation,
    LigandSymmetryAst,
    MulticenterBondAst,
    MulticenterBondConstraintAst,
    MulticenterBondConstraintsAst,
    MulticenterValenceAst,
    NoncovalentBondAst,
    NoncovalentBondConstraintAst,
    NoncovalentBondConstraintsAst,
    NoncovalentBondKind,
    NoncovalentBondKindAst,
    Orientation,
    OrientedLigandPermutation,
    Permutation,
    RelOp,
    RingMembershipAst,
    RingScope,
    StereoAtomAst,
    StereoBondAst,
    StereoConfigurationAst,
    StereoCoset,
    Stereogenicity,
    StereogenicityAst,
    StereoKind,
    StereoLigandPair,
    StereoTerm,
    TetrahedralStereoAst,
    Topicity,
    TopicityAst,
    TopicityRelationAst,
    UnpairedElectronsAst,
    ValueAst,
    ValuePredicate,
    ValueTerm,
)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(BooleanAst.Undetermined(), True, id="undetermined"),
        pytest.param(BooleanAst.Lit(True), False, id="literal"),
    ],
)
def test_boolean_ast_is_undetermined(value, expected):
    assert value.is_undetermined() is expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(BooleanAst.Undetermined(), False, id="undetermined"),
        pytest.param(BooleanAst.Lit(False), True, id="literal"),
    ],
)
def test_boolean_ast_is_ground(value, expected):
    assert value.is_ground() is expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            BooleanAst.Undetermined(),
            BooleanAst.Lit(True),
            BooleanAst.Lit(True),
            id="top-left",
        ),
        pytest.param(
            BooleanAst.Lit(False),
            BooleanAst.Undetermined(),
            BooleanAst.Lit(False),
            id="top-right",
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(True),
            BooleanAst.Lit(True),
            id="same",
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(False),
            None,
            id="incompatible",
        ),
    ],
)
def test_boolean_ast_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(True),
            BooleanAst.Lit(True),
            id="same",
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(False),
            BooleanAst.Undetermined(),
            id="different",
        ),
        pytest.param(
            BooleanAst.Undetermined(),
            BooleanAst.Lit(True),
            BooleanAst.Undetermined(),
            id="top",
        ),
    ],
)
def test_boolean_ast_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target", "expected"),
    [
        pytest.param(
            BooleanAst.Undetermined(), BooleanAst.Lit(True), True, id="top-pattern"
        ),
        pytest.param(
            BooleanAst.Lit(True), BooleanAst.Lit(True), True, id="same-literal"
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Undetermined(),
            False,
            id="top-target",
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(False),
            False,
            id="different-literal",
        ),
    ],
)
def test_boolean_ast_matches(pattern, target, expected):
    assert pattern.matches(target) is expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            BooleanAst.Undetermined(), BooleanAst.Lit(True), True, id="top-literal"
        ),
        pytest.param(
            BooleanAst.Lit(True), BooleanAst.Lit(True), True, id="same-literal"
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(False),
            False,
            id="different-literal",
        ),
    ],
)
def test_boolean_ast_is_compatible(lhs, rhs, expected):
    assert lhs.is_compatible(rhs) is expected


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(BooleanAst.Undetermined(), id="undetermined"),
        pytest.param(BooleanAst.Lit(True), id="literal"),
    ],
)
def test_boolean_ast_canonicalize(value):
    canonical = value.canonicalize()

    assert canonical == value


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            BooleanAst.Undetermined(),
            BooleanAst.Undetermined(),
            True,
            id="undetermined",
        ),
        pytest.param(
            BooleanAst.Lit(True), BooleanAst.Lit(True), True, id="same-literal"
        ),
        pytest.param(
            BooleanAst.Lit(True),
            BooleanAst.Lit(False),
            False,
            id="different-literal",
        ),
    ],
)
def test_boolean_ast_canonical_eq(lhs, rhs, expected):
    assert lhs.canonical_eq(rhs) is expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            ValueTerm.Sum([ValueTerm.Var("x"), ValueTerm.Lit(0)]),
            ValueTerm.Var("x"),
            id="sum-identity",
        ),
        pytest.param(
            ValueTerm.Product([ValueTerm.Lit(3), ValueTerm.Lit(2)]),
            ValueTerm.Lit(6),
            id="product",
        ),
    ],
)
def test_value_term_canonicalize(value, expected):
    assert value.canonicalize() == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            ValueTerm.Sum([ValueTerm.Var("x"), ValueTerm.Lit(0)]),
            ValueTerm.Var("x"),
            True,
            id="canonical",
        ),
        pytest.param(ValueTerm.Var("x"), ValueTerm.Var("y"), False, id="different"),
    ],
)
def test_value_term_canonical_eq(lhs, rhs, expected):
    assert lhs.canonical_eq(rhs) is expected


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(ValueAst.Undetermined(), True, False, id="undetermined"),
        pytest.param(ValueAst.Lit(1), False, True, id="literal"),
        pytest.param(ValueAst.LitSet({1, 2}), False, False, id="set"),
    ],
)
def test_value_ast_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            ValueAst.LitSet({1, 2}), ValueAst.Lit(2), ValueAst.Lit(2), id="set-literal"
        ),
        pytest.param(ValueAst.Lit(1), ValueAst.Lit(2), None, id="incompatible"),
    ],
)
def test_value_ast_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            ValueAst.Lit(1),
            ValueAst.Lit(2),
            ValueAst.LitSet({1, 2}),
            id="different-literals",
        ),
        pytest.param(
            ValueAst.RangeFrom(2),
            ValueAst.RangeFrom(5),
            ValueAst.RangeFrom(2),
            id="ranges",
        ),
    ],
)
def test_value_ast_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target", "expected"),
    [
        pytest.param(ValueAst.LitSet({1, 2}), ValueAst.Lit(2), True, id="refinement"),
        pytest.param(
            ValueAst.Lit(2), ValueAst.LitSet({1, 2}), False, id="generalization"
        ),
    ],
)
def test_value_ast_matches(pattern, target, expected):
    assert pattern.matches(target) is expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(ValueAst.LitSet({1, 2}), ValueAst.Lit(2), True, id="overlap"),
        pytest.param(ValueAst.Lit(1), ValueAst.Lit(2), False, id="disjoint"),
    ],
)
def test_value_ast_is_compatible(lhs, rhs, expected):
    assert lhs.is_compatible(rhs) is expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(ValueAst.LitSet({3}), ValueAst.Lit(3), id="singleton-set"),
        pytest.param(
            ValueAst.Term(ValueTerm.Sum([ValueTerm.Lit(2), ValueTerm.Lit(3)])),
            ValueAst.Lit(5),
            id="ground-term",
        ),
    ],
)
def test_value_ast_canonicalize(value, expected):
    assert value.canonicalize() == expected


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(ValueAst.LitSet(set()), id="empty-set"),
        pytest.param(
            ValueAst.Predicate(
                ValuePredicate.Rel(ValueTerm.Lit(2), RelOp.Lt, ValueTerm.Lit(1))
            ),
            id="false-predicate",
        ),
    ],
)
def test_value_ast_canonicalize_error(value):
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        value.canonicalize()


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(ValueAst.LitSet({3}), ValueAst.Lit(3), True, id="canonical"),
        pytest.param(ValueAst.Lit(3), ValueAst.Lit(4), False, id="different"),
        pytest.param(
            ValueAst.LitSet(set()),
            ValueAst.Predicate(
                ValuePredicate.Rel(ValueTerm.Lit(2), RelOp.Lt, ValueTerm.Lit(1))
            ),
            True,
            id="contradictions",
        ),
    ],
)
def test_value_ast_canonical_eq(lhs, rhs, expected):
    assert lhs.canonical_eq(rhs) is expected


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(ElementAst.Undetermined(), True, False, id="element-undetermined"),
        pytest.param(ElementAst.Lit(Element("C")), False, True, id="element-literal"),
        pytest.param(
            ElementAst.LitSet({Element("C"), Element("N")}),
            False,
            False,
            id="element-set",
        ),
        pytest.param(
            IsotopeMassAst.Undetermined(), True, False, id="isotope-undetermined"
        ),
        pytest.param(IsotopeMassAst.Natural(), False, True, id="isotope-natural"),
        pytest.param(IsotopeMassAst.LitSet({12, 13}), False, False, id="isotope-set"),
        pytest.param(
            UnpairedElectronsAst(ValueAst.Undetermined(), ValueAst.Undetermined()),
            True,
            False,
            id="unpaired-undetermined",
        ),
        pytest.param(
            UnpairedElectronsAst(ValueAst.LitSet({0, 2}), 3),
            False,
            False,
            id="unpaired-partial",
        ),
        pytest.param(UnpairedElectronsAst(2, 3), False, True, id="unpaired-ground"),
        pytest.param(
            ElectronCountsAst.Undetermined(), True, False, id="electron-undetermined"
        ),
        pytest.param(ElectronCountsAst.Lit([1, 1]), False, True, id="electron-literal"),
        pytest.param(
            AromaticValenceAst.Undetermined(), True, False, id="aromatic-undetermined"
        ),
        pytest.param(AromaticValenceAst.NotAromatic(), False, True, id="not-aromatic"),
        pytest.param(
            AromaticValenceAst.Aromatic(ValueAst.LitSet({1, 2})),
            False,
            False,
            id="aromatic-set",
        ),
        pytest.param(
            MulticenterValenceAst.Undetermined(),
            True,
            False,
            id="multicenter-undetermined",
        ),
        pytest.param(
            MulticenterValenceAst.NotMulticenter(),
            False,
            True,
            id="not-multicenter",
        ),
        pytest.param(
            MulticenterValenceAst.Multicenter(ValueAst.LitSet({2, 3})),
            False,
            False,
            id="multicenter-set",
        ),
    ],
)
def test_lattice_leaf_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            ElementAst.LitSet({Element("C"), Element("N")}),
            ElementAst.Lit(Element("C")),
            ElementAst.Lit(Element("C")),
            id="element-compatible",
        ),
        pytest.param(
            ElementAst.Lit(Element("C")),
            ElementAst.Lit(Element("N")),
            None,
            id="element-incompatible",
        ),
        pytest.param(
            IsotopeMassAst.LitSet({12, 13}),
            IsotopeMassAst.Lit(13),
            IsotopeMassAst.Lit(13),
            id="isotope-compatible",
        ),
        pytest.param(
            IsotopeMassAst.Natural(),
            IsotopeMassAst.Lit(13),
            None,
            id="isotope-incompatible",
        ),
        pytest.param(
            UnpairedElectronsAst(ValueAst.LitSet({0, 2}), ValueAst.LitSet({1, 3})),
            UnpairedElectronsAst(2, 3),
            UnpairedElectronsAst(2, 3),
            id="unpaired-compatible",
        ),
        pytest.param(
            UnpairedElectronsAst(2, 3),
            UnpairedElectronsAst(0, 1),
            None,
            id="unpaired-incompatible",
        ),
        pytest.param(
            ElectronCountsAst.Undetermined(),
            ElectronCountsAst.Lit([1, 1]),
            ElectronCountsAst.Lit([1, 1]),
            id="electron-compatible",
        ),
        pytest.param(
            ElectronCountsAst.Lit([1, 1]),
            ElectronCountsAst.Lit([2, 0]),
            None,
            id="electron-incompatible",
        ),
        pytest.param(
            AromaticValenceAst.Aromatic(ValueAst.LitSet({1, 2})),
            AromaticValenceAst.Aromatic(1),
            AromaticValenceAst.Aromatic(1),
            id="aromatic-compatible",
        ),
        pytest.param(
            AromaticValenceAst.NotAromatic(),
            AromaticValenceAst.Aromatic(1),
            None,
            id="aromatic-incompatible",
        ),
        pytest.param(
            MulticenterValenceAst.Multicenter(ValueAst.LitSet({2, 3})),
            MulticenterValenceAst.Multicenter(2),
            MulticenterValenceAst.Multicenter(2),
            id="multicenter-compatible",
        ),
        pytest.param(
            MulticenterValenceAst.NotMulticenter(),
            MulticenterValenceAst.Multicenter(2),
            None,
            id="multicenter-incompatible",
        ),
    ],
)
def test_lattice_leaf_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            ElementAst.Lit(Element("C")),
            ElementAst.Lit(Element("N")),
            ElementAst.LitSet({Element("C"), Element("N")}),
            id="element",
        ),
        pytest.param(
            IsotopeMassAst.Lit(12),
            IsotopeMassAst.Lit(13),
            IsotopeMassAst.LitSet({12, 13}),
            id="isotope",
        ),
        pytest.param(
            UnpairedElectronsAst(2, 3),
            UnpairedElectronsAst(0, 1),
            UnpairedElectronsAst(ValueAst.LitSet({0, 2}), ValueAst.LitSet({1, 3})),
            id="unpaired",
        ),
        pytest.param(
            ElectronCountsAst.Lit([1, 1]),
            ElectronCountsAst.Lit([2, 0]),
            ElectronCountsAst.Undetermined(),
            id="electron",
        ),
        pytest.param(
            AromaticValenceAst.NotAromatic(),
            AromaticValenceAst.Aromatic(1),
            AromaticValenceAst.Undetermined(),
            id="aromatic",
        ),
        pytest.param(
            MulticenterValenceAst.NotMulticenter(),
            MulticenterValenceAst.Multicenter(2),
            MulticenterValenceAst.Undetermined(),
            id="multicenter",
        ),
    ],
)
def test_lattice_leaf_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            ElementAst.LitSet({Element("C"), Element("N")}),
            ElementAst.Lit(Element("C")),
            id="element",
        ),
        pytest.param(
            IsotopeMassAst.LitSet({12, 13}), IsotopeMassAst.Lit(13), id="isotope"
        ),
        pytest.param(
            UnpairedElectronsAst(ValueAst.LitSet({0, 2}), ValueAst.LitSet({1, 3})),
            UnpairedElectronsAst(2, 3),
            id="unpaired",
        ),
        pytest.param(
            ElectronCountsAst.Undetermined(),
            ElectronCountsAst.Lit([1, 1]),
            id="electron",
        ),
        pytest.param(
            AromaticValenceAst.Aromatic(ValueAst.LitSet({1, 2})),
            AromaticValenceAst.Aromatic(1),
            id="aromatic",
        ),
        pytest.param(
            MulticenterValenceAst.Multicenter(ValueAst.LitSet({2, 3})),
            MulticenterValenceAst.Multicenter(2),
            id="multicenter",
        ),
    ],
)
def test_lattice_leaf_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(
            ElementAst.Lit(Element("C")), ElementAst.Lit(Element("N")), id="element"
        ),
        pytest.param(IsotopeMassAst.Natural(), IsotopeMassAst.Lit(13), id="isotope"),
        pytest.param(
            UnpairedElectronsAst(2, 3), UnpairedElectronsAst(0, 1), id="unpaired"
        ),
        pytest.param(
            ElectronCountsAst.Lit([1, 1]),
            ElectronCountsAst.Lit([2, 0]),
            id="electron",
        ),
        pytest.param(
            AromaticValenceAst.NotAromatic(),
            AromaticValenceAst.Aromatic(1),
            id="aromatic",
        ),
        pytest.param(
            MulticenterValenceAst.NotMulticenter(),
            MulticenterValenceAst.Multicenter(2),
            id="multicenter",
        ),
    ],
)
def test_lattice_leaf_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            ElementAst.LitSet({Element("C")}),
            ElementAst.Lit(Element("C")),
            id="element",
        ),
        pytest.param(IsotopeMassAst.LitSet({13}), IsotopeMassAst.Lit(13), id="isotope"),
        pytest.param(
            UnpairedElectronsAst(
                ValueAst.LitSet({2}),
                ValueAst.Term(ValueTerm.Sum([ValueTerm.Lit(1), ValueTerm.Lit(2)])),
            ),
            UnpairedElectronsAst(2, 3),
            id="unpaired",
        ),
        pytest.param(
            ElectronCountsAst.Lit([1, 1]),
            ElectronCountsAst.Lit([1, 1]),
            id="electron",
        ),
        pytest.param(
            AromaticValenceAst.Aromatic(ValueAst.LitSet({1})),
            AromaticValenceAst.Aromatic(1),
            id="aromatic",
        ),
        pytest.param(
            MulticenterValenceAst.Multicenter(ValueAst.LitSet({2})),
            MulticenterValenceAst.Multicenter(2),
            id="multicenter",
        ),
    ],
)
def test_lattice_leaf_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(ElementAst.LitSet(set()), id="element"),
        pytest.param(IsotopeMassAst.LitSet(set()), id="isotope"),
        pytest.param(UnpairedElectronsAst(ValueAst.LitSet(set()), 1), id="unpaired"),
        pytest.param(
            AromaticValenceAst.Aromatic(ValueAst.LitSet(set())), id="aromatic"
        ),
        pytest.param(
            MulticenterValenceAst.Multicenter(ValueAst.LitSet(set())),
            id="multicenter",
        ),
    ],
)
def test_lattice_leaf_canonicalize_error(value):
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        value.canonicalize()


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            StereoConfigurationAst.Undetermined(), True, False, id="configuration-top"
        ),
        pytest.param(
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
            False,
            False,
            id="configuration-open",
        ),
        pytest.param(
            TetrahedralStereoAst.Undetermined(), True, False, id="tetrahedral-top"
        ),
        pytest.param(
            TetrahedralStereoAst.NotStereo(), False, True, id="tetrahedral-ground"
        ),
        pytest.param(CisTransStereoAst.Undetermined(), True, False, id="cis-trans-top"),
        pytest.param(
            CisTransStereoAst.Stereo(StereoCoset.Lit(0)),
            False,
            True,
            id="cis-trans-ground",
        ),
        pytest.param(
            StereogenicityAst.Undetermined(), True, False, id="stereogenicity-top"
        ),
        pytest.param(
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            False,
            True,
            id="stereogenicity-ground",
        ),
        pytest.param(
            TopicityRelationAst.Undetermined(), True, False, id="topicity-relation-top"
        ),
        pytest.param(
            TopicityRelationAst.Lit(Topicity.Homotopic),
            False,
            True,
            id="topicity-relation-ground",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Undetermined(),
            ),
            True,
            False,
            id="ligand-symmetry-top",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(True)
            ),
            False,
            True,
            id="fluxionality-ground",
        ),
        pytest.param(
            TopicityAst(StereoLigandPair(0, 1), TopicityRelationAst.Undetermined()),
            True,
            False,
            id="topicity-top",
        ),
        pytest.param(
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            False,
            True,
            id="topicity-ground",
        ),
    ],
)
def test_stereo_leaf_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            id="configuration-compatible",
        ),
        pytest.param(
            TetrahedralStereoAst.Stereo(StereoCoset.Undetermined()),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)),
            id="tetrahedral-compatible",
        ),
        pytest.param(
            CisTransStereoAst.NotStereo(),
            CisTransStereoAst.Stereo(StereoCoset.Lit(0)),
            None,
            id="cis-trans-incompatible",
        ),
        pytest.param(
            StereogenicityAst.LitSet(
                {Stereogenicity.Prochiral, Stereogenicity.Stereogenic}
            ),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            id="stereogenicity-compatible",
        ),
        pytest.param(
            TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            TopicityRelationAst.Lit(Topicity.Homotopic),
            TopicityRelationAst.Lit(Topicity.Homotopic),
            id="topicity-relation-compatible",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(False),
            ),
            None,
            id="ligand-symmetry-incompatible",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])),
                BooleanAst.Undetermined(),
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            id="fluxionality-compatible",
        ),
        pytest.param(
            TopicityAst(
                StereoLigandPair(0, 1),
                TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            ),
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            id="topicity-compatible",
        ),
    ],
)
def test_stereo_leaf_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(1)),
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCoset.LitSet({0, 1})
            ),
            id="configuration",
        ),
        pytest.param(
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(0)),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)),
            TetrahedralStereoAst.Stereo(StereoCoset.LitSet({0, 1})),
            id="tetrahedral",
        ),
        pytest.param(
            CisTransStereoAst.NotStereo(),
            CisTransStereoAst.Stereo(StereoCoset.Lit(0)),
            CisTransStereoAst.Undetermined(),
            id="cis-trans",
        ),
        pytest.param(
            StereogenicityAst.Lit(Stereogenicity.Prochiral),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            StereogenicityAst.NotSet({Stereogenicity.Symmetric}),
            id="stereogenicity",
        ),
        pytest.param(
            TopicityRelationAst.Lit(Topicity.Homotopic),
            TopicityRelationAst.Lit(Topicity.Enantiotopic),
            TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            id="topicity-relation",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([1, 0])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            None,
            id="ligand-symmetry-different-fiber",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(True)
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])),
                BooleanAst.Undetermined(),
            ),
            id="fluxionality",
        ),
        pytest.param(
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            TopicityAst(StereoLigandPair(0, 1), Topicity.Enantiotopic),
            TopicityAst(
                StereoLigandPair(0, 1),
                TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            ),
            id="topicity",
        ),
    ],
)
def test_stereo_leaf_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            id="configuration",
        ),
        pytest.param(
            TetrahedralStereoAst.Stereo(StereoCoset.Undetermined()),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(0)),
            id="tetrahedral",
        ),
        pytest.param(
            CisTransStereoAst.Undetermined(),
            CisTransStereoAst.NotStereo(),
            id="cis-trans",
        ),
        pytest.param(
            StereogenicityAst.NotSet({Stereogenicity.Symmetric}),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            id="stereogenicity",
        ),
        pytest.param(
            TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            TopicityRelationAst.Lit(Topicity.Homotopic),
            id="topicity-relation",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Undetermined(),
            ),
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            id="ligand-symmetry",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])),
                BooleanAst.Undetermined(),
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            id="fluxionality",
        ),
        pytest.param(
            TopicityAst(
                StereoLigandPair(0, 1),
                TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            ),
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            id="topicity",
        ),
    ],
)
def test_stereo_leaf_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCoset.Lit(0)),
            id="configuration",
        ),
        pytest.param(
            TetrahedralStereoAst.NotStereo(),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(0)),
            id="tetrahedral",
        ),
        pytest.param(
            CisTransStereoAst.Stereo(StereoCoset.Lit(0)),
            CisTransStereoAst.Stereo(StereoCoset.Lit(1)),
            id="cis-trans",
        ),
        pytest.param(
            StereogenicityAst.Lit(Stereogenicity.Symmetric),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            id="stereogenicity",
        ),
        pytest.param(
            TopicityRelationAst.Lit(Topicity.Homotopic),
            TopicityRelationAst.Lit(Topicity.Diastereotopic),
            id="topicity-relation",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([1, 0])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            id="ligand-symmetry",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([0, 1])), BooleanAst.Lit(True)
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(True)
            ),
            id="fluxionality",
        ),
        pytest.param(
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            TopicityAst(StereoLigandPair(0, 2), Topicity.Homotopic),
            id="topicity",
        ),
    ],
)
def test_stereo_leaf_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral,
                StereoCoset.Term(StereoTerm.Swap(StereoTerm.Lit(0))),
            ),
            StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(1)),
            id="configuration",
        ),
        pytest.param(
            TetrahedralStereoAst.Stereo(
                StereoCoset.Term(StereoTerm.Swap(StereoTerm.Lit(0)))
            ),
            TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)),
            id="tetrahedral",
        ),
        pytest.param(
            CisTransStereoAst.Stereo(
                StereoCoset.Term(StereoTerm.Swap(StereoTerm.Lit(0)))
            ),
            CisTransStereoAst.Stereo(StereoCoset.Lit(1)),
            id="cis-trans",
        ),
        pytest.param(
            StereogenicityAst.LitSet({Stereogenicity.Stereogenic}),
            StereogenicityAst.Lit(Stereogenicity.Stereogenic),
            id="stereogenicity",
        ),
        pytest.param(
            TopicityRelationAst.LitSet({Topicity.Homotopic, Topicity.Enantiotopic}),
            TopicityRelationAst.NotSet({Topicity.Diastereotopic}),
            id="topicity-relation",
        ),
        pytest.param(
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            LigandSymmetryAst(
                OrientedLigandPermutation(
                    LigandPermutation(Permutation([0, 1])), Orientation.Proper
                ),
                BooleanAst.Lit(True),
            ),
            id="ligand-symmetry",
        ),
        pytest.param(
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            FluxionalityAst(
                LigandPermutation(Permutation([1, 0])), BooleanAst.Lit(False)
            ),
            id="fluxionality",
        ),
        pytest.param(
            TopicityAst(
                StereoLigandPair(0, 1),
                TopicityRelationAst.LitSet({Topicity.Homotopic}),
            ),
            TopicityAst(StereoLigandPair(0, 1), Topicity.Homotopic),
            id="topicity",
        ),
    ],
)
def test_stereo_leaf_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            NoncovalentBondKindAst.Undetermined(),
            True,
            False,
            id="noncovalent-top",
        ),
        pytest.param(
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            False,
            True,
            id="noncovalent-ground",
        ),
        pytest.param(
            RingMembershipAst(RingScope.All(), ValueAst.Undetermined()),
            True,
            False,
            id="ring-top",
        ),
        pytest.param(
            RingMembershipAst(RingScope.Size(6), 2),
            False,
            True,
            id="ring-ground",
        ),
    ],
)
def test_remaining_leaf_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            NoncovalentBondKindAst.Undetermined(),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            id="noncovalent-compatible",
        ),
        pytest.param(
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic),
            None,
            id="noncovalent-incompatible",
        ),
        pytest.param(
            RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2})),
            RingMembershipAst(RingScope.Size(6), 2),
            RingMembershipAst(RingScope.Size(6), 2),
            id="ring-compatible",
        ),
        pytest.param(
            RingMembershipAst(RingScope.All(), 2),
            RingMembershipAst(RingScope.Size(6), 2),
            None,
            id="ring-different-fiber",
        ),
    ],
)
def test_remaining_leaf_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic),
            NoncovalentBondKindAst.Undetermined(),
            id="noncovalent",
        ),
        pytest.param(
            RingMembershipAst(RingScope.Size(6), 1),
            RingMembershipAst(RingScope.Size(6), 2),
            RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2})),
            id="ring",
        ),
        pytest.param(
            RingMembershipAst(RingScope.All(), 2),
            RingMembershipAst(RingScope.Size(6), 2),
            None,
            id="ring-different-fiber",
        ),
    ],
)
def test_remaining_leaf_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            NoncovalentBondKindAst.Undetermined(),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HalogenBond),
            id="noncovalent",
        ),
        pytest.param(
            RingMembershipAst(RingScope.All(), ValueAst.LitSet({1, 2})),
            RingMembershipAst(RingScope.All(), 2),
            id="ring",
        ),
    ],
)
def test_remaining_leaf_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic),
            id="noncovalent",
        ),
        pytest.param(
            RingMembershipAst(RingScope.All(), 2),
            RingMembershipAst(RingScope.Size(6), 2),
            id="ring",
        ),
    ],
)
def test_remaining_leaf_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.VanDerWaals),
            NoncovalentBondKindAst.Lit(NoncovalentBondKind.VanDerWaals),
            id="noncovalent",
        ),
        pytest.param(
            RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({2})),
            RingMembershipAst(RingScope.Size(6), 2),
            id="ring",
        ),
    ],
)
def test_remaining_leaf_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(
            RingMembershipAst(RingScope.All(), ValueAst.LitSet(set())), id="ring"
        ),
    ],
)
def test_remaining_leaf_canonicalize_error(value):
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        value.canonicalize()


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            AtomAst(ElementAst.Undetermined()), True, False, id="atom-undetermined"
        ),
        pytest.param(
            AtomAst(
                Element("C"),
                isotope_mass=IsotopeMassAst.Natural(),
                charge=0,
                implicit_hydrogens=4,
                lone_pairs=0,
                unpaired_electrons=UnpairedElectronsAst(0, 1),
            ),
            False,
            True,
            id="atom-ground",
        ),
        pytest.param(
            BondAst(ValueAst.Undetermined()), True, False, id="bond-undetermined"
        ),
        pytest.param(
            BondAst(
                1,
                charge=0,
                unpaired_electrons=UnpairedElectronsAst(0, 1),
            ),
            False,
            True,
            id="bond-ground",
        ),
    ],
)
def test_entity_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            AtomAst(
                ElementAst.LitSet({Element("C"), Element("N")}),
                charge=ValueAst.LitSet({0, 1}),
            ),
            AtomAst(Element("C"), charge=0),
            AtomAst(Element("C"), charge=0),
            id="atom-compatible",
        ),
        pytest.param(
            AtomAst(Element("C")),
            AtomAst(Element("N")),
            None,
            id="atom-incompatible",
        ),
        pytest.param(
            BondAst(ValueAst.LitSet({1, 2}), charge=ValueAst.LitSet({0, 1})),
            BondAst(1, charge=0),
            BondAst(1, charge=0),
            id="bond-compatible",
        ),
        pytest.param(
            BondAst(1),
            BondAst(2),
            None,
            id="bond-incompatible",
        ),
    ],
)
def test_entity_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            AtomAst(Element("C"), charge=0),
            AtomAst(Element("N"), charge=1),
            AtomAst(
                ElementAst.LitSet({Element("C"), Element("N")}),
                charge=ValueAst.LitSet({0, 1}),
            ),
            id="atom",
        ),
        pytest.param(
            BondAst(1, charge=0),
            BondAst(2, charge=1),
            BondAst(ValueAst.LitSet({1, 2}), charge=ValueAst.LitSet({0, 1})),
            id="bond",
        ),
    ],
)
def test_entity_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            AtomAst(
                ElementAst.LitSet({Element("C"), Element("N")}),
                charge=ValueAst.LitSet({0, 1}),
            ),
            AtomAst(Element("C"), charge=0),
            id="atom",
        ),
        pytest.param(
            BondAst(ValueAst.LitSet({1, 2}), charge=ValueAst.LitSet({0, 1})),
            BondAst(1, charge=0),
            id="bond",
        ),
    ],
)
def test_entity_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(AtomAst(Element("C")), AtomAst(Element("N")), id="atom"),
        pytest.param(BondAst(1), BondAst(2), id="bond"),
    ],
)
def test_entity_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            AtomAst(ElementAst.LitSet({Element("C")}), charge=ValueAst.LitSet({0})),
            AtomAst(Element("C"), charge=0),
            id="atom",
        ),
        pytest.param(
            BondAst(ValueAst.LitSet({1}), charge=ValueAst.LitSet({0})),
            BondAst(1, charge=0),
            id="bond",
        ),
    ],
)
def test_entity_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(AtomAst(Element("C"), charge=ValueAst.LitSet(set())), id="atom"),
        pytest.param(BondAst(ValueAst.LitSet(set())), id="bond"),
    ],
)
def test_entity_canonicalize_error(value):
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        value.canonicalize()


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            DativeBondAst(ValueAst.Undetermined()),
            True,
            False,
            id="dative-undetermined",
        ),
        pytest.param(DativeBondAst(1), False, True, id="dative-ground"),
        pytest.param(
            AromaticSystemAst(ElectronCountsAst.Undetermined()),
            True,
            False,
            id="aromatic-undetermined",
        ),
        pytest.param(
            AromaticSystemAst(
                [1, 1],
                charge=0,
                unpaired_electrons=UnpairedElectronsAst(0, 1),
            ),
            False,
            True,
            id="aromatic-ground",
        ),
        pytest.param(
            MulticenterBondAst(ElectronCountsAst.Undetermined()),
            True,
            False,
            id="multicenter-undetermined",
        ),
        pytest.param(
            MulticenterBondAst(
                [1, 0, 1],
                charge=0,
                unpaired_electrons=UnpairedElectronsAst(0, 1),
            ),
            False,
            True,
            id="multicenter-ground",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKindAst.Undetermined()),
            True,
            False,
            id="noncovalent-undetermined",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            False,
            True,
            id="noncovalent-ground",
        ),
    ],
)
def test_overlay_entity_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            DativeBondAst(ValueAst.LitSet({1, 2})),
            DativeBondAst(1),
            DativeBondAst(1),
            id="dative-compatible",
        ),
        pytest.param(
            DativeBondAst(1),
            DativeBondAst(2),
            None,
            id="dative-incompatible",
        ),
        pytest.param(
            AromaticSystemAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            AromaticSystemAst([1, 1], charge=0),
            AromaticSystemAst([1, 1], charge=0),
            id="aromatic-compatible",
        ),
        pytest.param(
            AromaticSystemAst([1, 1]),
            AromaticSystemAst([2, 0]),
            None,
            id="aromatic-incompatible",
        ),
        pytest.param(
            MulticenterBondAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            MulticenterBondAst([1, 0, 1], charge=0),
            MulticenterBondAst([1, 0, 1], charge=0),
            id="multicenter-compatible",
        ),
        pytest.param(
            MulticenterBondAst([1, 0, 1]),
            MulticenterBondAst([2, 0, 0]),
            None,
            id="multicenter-incompatible",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKindAst.Undetermined()),
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            id="noncovalent-compatible",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondAst(NoncovalentBondKind.Ionic),
            None,
            id="noncovalent-incompatible",
        ),
    ],
)
def test_overlay_entity_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            DativeBondAst(1),
            DativeBondAst(2),
            DativeBondAst(ValueAst.LitSet({1, 2})),
            id="dative",
        ),
        pytest.param(
            AromaticSystemAst([1, 1], charge=0),
            AromaticSystemAst([2, 0], charge=1),
            AromaticSystemAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            id="aromatic",
        ),
        pytest.param(
            MulticenterBondAst([1, 0, 1], charge=0),
            MulticenterBondAst([2, 0, 0], charge=1),
            MulticenterBondAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            id="multicenter",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondAst(NoncovalentBondKind.Ionic),
            NoncovalentBondAst(NoncovalentBondKindAst.Undetermined()),
            id="noncovalent",
        ),
    ],
)
def test_overlay_entity_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            DativeBondAst(ValueAst.LitSet({1, 2})),
            DativeBondAst(1),
            id="dative",
        ),
        pytest.param(
            AromaticSystemAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            AromaticSystemAst([1, 1], charge=0),
            id="aromatic",
        ),
        pytest.param(
            MulticenterBondAst(
                ElectronCountsAst.Undetermined(), charge=ValueAst.LitSet({0, 1})
            ),
            MulticenterBondAst([1, 0, 1], charge=0),
            id="multicenter",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKindAst.Undetermined()),
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            id="noncovalent",
        ),
    ],
)
def test_overlay_entity_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(DativeBondAst(1), DativeBondAst(2), id="dative"),
        pytest.param(
            AromaticSystemAst([1, 1]),
            AromaticSystemAst([2, 0]),
            id="aromatic",
        ),
        pytest.param(
            MulticenterBondAst([1, 0, 1]),
            MulticenterBondAst([2, 0, 0]),
            id="multicenter",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            NoncovalentBondAst(NoncovalentBondKind.Ionic),
            id="noncovalent",
        ),
    ],
)
def test_overlay_entity_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            DativeBondAst(ValueAst.LitSet({1})),
            DativeBondAst(1),
            id="dative",
        ),
        pytest.param(
            AromaticSystemAst([1, 1], charge=ValueAst.LitSet({0})),
            AromaticSystemAst([1, 1], charge=0),
            id="aromatic",
        ),
        pytest.param(
            MulticenterBondAst([1, 0, 1], charge=ValueAst.LitSet({0})),
            MulticenterBondAst([1, 0, 1], charge=0),
            id="multicenter",
        ),
        pytest.param(
            NoncovalentBondAst(NoncovalentBondKind.VanDerWaals),
            NoncovalentBondAst(NoncovalentBondKind.VanDerWaals),
            id="noncovalent",
        ),
    ],
)
def test_overlay_entity_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            StereoAtomAst(StereoConfigurationAst.Undetermined()),
            True,
            False,
            id="stereo-atom-undetermined",
        ),
        pytest.param(
            StereoAtomAst.parse("Th*"),
            False,
            False,
            id="stereo-atom-kinded",
        ),
        pytest.param(
            StereoAtomAst.parse("Th0"),
            False,
            True,
            id="stereo-atom-ground",
        ),
        pytest.param(
            StereoBondAst(StereoConfigurationAst.Undetermined()),
            True,
            False,
            id="stereo-bond-undetermined",
        ),
        pytest.param(
            StereoBondAst.parse("Ct*"),
            False,
            False,
            id="stereo-bond-kinded",
        ),
        pytest.param(
            StereoBondAst.parse("Ct0"),
            False,
            True,
            id="stereo-bond-ground",
        ),
    ],
)
def test_stereo_entity_classification(value, expected_undetermined, expected_ground):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            StereoAtomAst.parse("Th*"),
            StereoAtomAst.parse("Th0"),
            StereoAtomAst.parse("Th0"),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondAst.parse("Ct*"),
            StereoBondAst.parse("Ct0"),
            StereoBondAst.parse("Ct0"),
            id="stereo-bond",
        ),
    ],
)
def test_stereo_entity_meet(lhs, rhs, expected):
    result = lhs.meet(rhs)
    assert result is not None
    assert result.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            StereoAtomAst.parse("Th0"),
            StereoAtomAst.parse("Th1"),
            StereoAtomAst.parse("Th{0,1}"),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondAst.parse("Ct0"),
            StereoBondAst.parse("Ct1"),
            StereoBondAst.parse("Ct{0,1}"),
            id="stereo-bond",
        ),
    ],
)
def test_stereo_entity_join(lhs, rhs, expected):
    result = lhs.join(rhs)
    assert result is not None
    assert result.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            StereoAtomAst.parse("Th*"),
            StereoAtomAst.parse("Th0"),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondAst.parse("Ct*"),
            StereoBondAst.parse("Ct0"),
            id="stereo-bond",
        ),
    ],
)
def test_stereo_entity_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(
            StereoAtomAst.parse("Th0"),
            StereoAtomAst.parse("Th1"),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondAst.parse("Ct0"),
            StereoBondAst.parse("Ct1"),
            id="stereo-bond",
        ),
    ],
)
def test_stereo_entity_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            StereoAtomAst.parse("Th{0}"),
            StereoAtomAst.parse("Th0"),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondAst.parse("Ct{0}"),
            StereoBondAst.parse("Ct0"),
            id="stereo-bond",
        ),
    ],
)
def test_stereo_entity_canonicalize(value, expected):
    canonical = value.canonicalize()
    assert str(canonical) == str(expected)
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("value", "expected_undetermined", "expected_ground"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.Undetermined()),
            True,
            False,
            id="atom-constraint-undetermined",
        ),
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            False,
            True,
            id="atom-constraint-ground",
        ),
        pytest.param(
            BondConstraintAst.Aromatic(BooleanAst.Undetermined()),
            True,
            False,
            id="bond-constraint-undetermined",
        ),
        pytest.param(
            BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            False,
            True,
            id="bond-constraint-ground",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Undetermined())]),
            True,
            False,
            id="atom-constraints-undetermined",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
            False,
            True,
            id="atom-constraints-ground",
        ),
        pytest.param(
            BondConstraintsAst([BondConstraintAst.Aromatic(BooleanAst.Undetermined())]),
            True,
            False,
            id="bond-constraints-undetermined",
        ),
        pytest.param(
            BondConstraintsAst([BondConstraintAst.Aromatic(BooleanAst.Lit(True))]),
            False,
            True,
            id="bond-constraints-ground",
        ),
    ],
)
def test_entity_constraint_classification(
    value, expected_undetermined, expected_ground
):
    assert value.is_undetermined() is expected_undetermined
    assert value.is_ground() is expected_ground


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.LitSet({3, 4})),
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            id="atom-constraint",
        ),
        pytest.param(
            BondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 1)),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 1)),
            id="bond-constraint",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.LitSet({3, 4}))]),
            AtomConstraintsAst(
                [
                    AtomConstraintAst.Valence(ValueAst.Lit(4)),
                    AtomConstraintAst.Degree(ValueAst.Lit(3)),
                ]
            ),
            AtomConstraintsAst(
                [
                    AtomConstraintAst.Valence(ValueAst.Lit(4)),
                    AtomConstraintAst.Degree(ValueAst.Lit(3)),
                ]
            ),
            id="atom-constraints",
        ),
        pytest.param(
            BondConstraintsAst(
                [
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    ),
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    ),
                ]
            ),
            id="bond-constraints",
        ),
    ],
)
def test_entity_constraint_meet(lhs, rhs, expected):
    assert lhs.meet(rhs) == expected


@pytest.mark.parametrize(
    ("lhs", "rhs", "expected"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.Lit(3)),
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Valence(ValueAst.LitSet({3, 4})),
            id="atom-constraint-same-fiber",
        ),
        pytest.param(
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 1)),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 2)),
            BondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            id="bond-constraint-same-fiber",
        ),
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(4)),
            None,
            id="atom-constraint-different-fiber",
        ),
        pytest.param(
            BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.All(), 1)),
            None,
            id="bond-constraint-different-fiber",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(3))]),
            AtomConstraintsAst(
                [
                    AtomConstraintAst.Valence(ValueAst.Lit(4)),
                    AtomConstraintAst.Degree(ValueAst.Lit(2)),
                ]
            ),
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.LitSet({3, 4}))]),
            id="atom-constraints",
        ),
        pytest.param(
            BondConstraintsAst(
                [
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 2)
                    ),
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            id="bond-constraints",
        ),
    ],
)
def test_entity_constraint_join(lhs, rhs, expected):
    assert lhs.join(rhs) == expected


@pytest.mark.parametrize(
    ("pattern", "target"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.LitSet({3, 4})),
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            id="atom-constraint",
        ),
        pytest.param(
            BondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 1)),
            id="bond-constraint",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.LitSet({3, 4}))]),
            AtomConstraintsAst(
                [
                    AtomConstraintAst.Valence(ValueAst.Lit(4)),
                    AtomConstraintAst.Degree(ValueAst.Lit(3)),
                ]
            ),
            id="atom-constraints",
        ),
        pytest.param(
            BondConstraintsAst(
                [
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    ),
                ]
            ),
            id="bond-constraints",
        ),
    ],
)
def test_entity_constraint_matches(pattern, target):
    assert pattern.matches(target) is True
    assert target.matches(pattern) is False


@pytest.mark.parametrize(
    ("lhs", "rhs"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(4)),
            id="atom-constraint",
        ),
        pytest.param(
            BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.All(), 1)),
            id="bond-constraint",
        ),
        pytest.param(
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(3))]),
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
            id="atom-constraints",
        ),
        pytest.param(
            BondConstraintsAst([BondConstraintAst.Aromatic(BooleanAst.Lit(True))]),
            BondConstraintsAst([BondConstraintAst.Aromatic(BooleanAst.Lit(False))]),
            id="bond-constraints",
        ),
    ],
)
def test_entity_constraint_is_compatible(lhs, rhs):
    assert lhs.is_compatible(rhs) is False
    assert rhs.is_compatible(lhs) is False


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(
            AtomConstraintAst.Valence(ValueAst.LitSet({4})),
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            id="atom-constraint",
        ),
        pytest.param(
            BondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1}))
            ),
            BondConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), 1)),
            id="bond-constraint",
        ),
        pytest.param(
            AtomConstraintsAst(
                [
                    AtomConstraintAst.Valence(ValueAst.LitSet({4})),
                    AtomConstraintAst.Degree(ValueAst.Undetermined()),
                ]
            ),
            AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
            id="atom-constraints",
        ),
        pytest.param(
            BondConstraintsAst(
                [
                    BondConstraintAst.Aromatic(BooleanAst.Undetermined()),
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1}))
                    ),
                ]
            ),
            BondConstraintsAst(
                [
                    BondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            id="bond-constraints",
        ),
    ],
)
def test_entity_constraint_canonicalize(value, expected):
    assert value.canonicalize() == expected
    assert value.canonical_eq(expected) is True


@pytest.mark.parametrize(
    ("open_constraint", "ground_constraint", "open_container", "ground_container"),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Undetermined()),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Undetermined())]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.Aromatic(BooleanAst.Undetermined()),
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            DativeBondConstraintsAst(
                [DativeBondConstraintAst.Aromatic(BooleanAst.Undetermined())]
            ),
            DativeBondConstraintsAst(
                [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Undetermined()),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Undetermined())]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined()),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined())]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
            ),
            id="noncovalent-bond",
        ),
    ],
)
def test_overlay_constraint_classification(
    open_constraint, ground_constraint, open_container, ground_container
):
    assert open_constraint.is_undetermined() is True
    assert open_constraint.is_ground() is False
    assert ground_constraint.is_undetermined() is False
    assert ground_constraint.is_ground() is True
    assert open_container.is_undetermined() is True
    assert open_container.is_ground() is False
    assert ground_container.is_undetermined() is False
    assert ground_container.is_ground() is True


@pytest.mark.parametrize(
    (
        "constraint_lhs",
        "constraint_rhs",
        "expected_constraint",
        "container_lhs",
        "container_rhs",
        "expected_container",
    ),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8})),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8}))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 1)
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 1)
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4})),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4}))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined()),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False)),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False)),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined())]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
            ),
            id="noncovalent-bond",
        ),
    ],
)
def test_overlay_constraint_meet(
    constraint_lhs,
    constraint_rhs,
    expected_constraint,
    container_lhs,
    container_rhs,
    expected_container,
):
    assert constraint_lhs.meet(constraint_rhs) == expected_constraint
    assert container_lhs.meet(container_rhs) == expected_container


@pytest.mark.parametrize(
    (
        "constraint_lhs",
        "constraint_rhs",
        "expected_constraint",
        "container_lhs",
        "container_rhs",
        "expected_container",
    ),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(8)),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8})),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(8))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8}))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 1)
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 2)
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 2)
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(2)),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4})),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(2))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4}))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False)),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined()),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
            ),
            NoncovalentBondConstraintsAst([]),
            id="noncovalent-bond",
        ),
        pytest.param(
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.All(), 1)
            ),
            None,
            DativeBondConstraintsAst(
                [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.All(), 1)
                    )
                ]
            ),
            DativeBondConstraintsAst([]),
            id="dative-bond-different-fiber",
        ),
    ],
)
def test_overlay_constraint_join(
    constraint_lhs,
    constraint_rhs,
    expected_constraint,
    container_lhs,
    container_rhs,
    expected_container,
):
    assert constraint_lhs.join(constraint_rhs) == expected_constraint
    assert container_lhs.join(container_rhs) == expected_container


@pytest.mark.parametrize(
    (
        "constraint_pattern",
        "constraint_target",
        "container_pattern",
        "container_target",
    ),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8})),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6, 8}))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 1)
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1, 2}))
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4})),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({2, 4}))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined()),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined())]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
            ),
            id="noncovalent-bond",
        ),
    ],
)
def test_overlay_constraint_matches(
    constraint_pattern, constraint_target, container_pattern, container_target
):
    assert constraint_pattern.matches(constraint_target) is True
    assert constraint_target.matches(constraint_pattern) is False
    assert container_pattern.matches(container_target) is True
    assert container_target.matches(container_pattern) is False


@pytest.mark.parametrize(
    ("constraint_lhs", "constraint_rhs", "container_lhs", "container_rhs"),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(8)),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(8))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(False)),
            DativeBondConstraintsAst(
                [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
            ),
            DativeBondConstraintsAst(
                [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(False))]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(2)),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(2))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False)),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
            ),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
            ),
            id="noncovalent-bond",
        ),
    ],
)
def test_overlay_constraint_is_compatible(
    constraint_lhs, constraint_rhs, container_lhs, container_rhs
):
    assert constraint_lhs.is_compatible(constraint_rhs) is False
    assert constraint_rhs.is_compatible(constraint_lhs) is False
    assert container_lhs.is_compatible(container_rhs) is False
    assert container_rhs.is_compatible(container_lhs) is False


@pytest.mark.parametrize(
    (
        "constraint",
        "expected_constraint",
        "container",
        "expected_container",
    ),
    [
        pytest.param(
            AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6})),
            AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.LitSet({6}))]
            ),
            AromaticSystemConstraintsAst(
                [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
            ),
            id="aromatic-system",
        ),
        pytest.param(
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1}))
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), 1)
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), ValueAst.LitSet({1}))
                    )
                ]
            ),
            DativeBondConstraintsAst(
                [
                    DativeBondConstraintAst.RingMembership(
                        RingMembershipAst(RingScope.Size(6), 1)
                    )
                ]
            ),
            id="dative-bond",
        ),
        pytest.param(
            MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({4})),
            MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4)),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.LitSet({4}))]
            ),
            MulticenterBondConstraintsAst(
                [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(4))]
            ),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
            NoncovalentBondConstraintsAst(
                [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Undetermined())]
            ),
            NoncovalentBondConstraintsAst([]),
            id="noncovalent-bond",
        ),
    ],
)
def test_overlay_constraint_canonicalize(
    constraint, expected_constraint, container, expected_container
):
    assert constraint.canonicalize() == expected_constraint
    assert constraint.canonical_eq(expected_constraint) is True
    assert container.canonicalize() == expected_container
    assert container.canonical_eq(expected_container) is True
