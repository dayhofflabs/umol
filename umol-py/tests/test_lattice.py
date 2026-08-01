import pytest

from umol import (
    AromaticValenceAst,
    BooleanAst,
    ContradictionError,
    ElectronCountsAst,
    Element,
    ElementAst,
    IsotopeMassAst,
    MulticenterValenceAst,
    RelOp,
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
