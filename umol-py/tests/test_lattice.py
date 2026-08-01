import pytest

from umol import BooleanAst


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
