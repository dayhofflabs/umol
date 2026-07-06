from umol import MemOp, RelOp, ValueAst, ValuePredicate, ValueTerm


def test_relop_eq():
    assert RelOp.Lt == RelOp.Lt
    assert RelOp.Lt != RelOp.Ge


def test_relop_match():
    match RelOp.Ne:
        case RelOp.Ne:
            pass
        case _:
            raise AssertionError


def test_memop_eq():
    assert MemOp.In == MemOp.In
    assert MemOp.In != MemOp.NotIn


def test_valueterm_lit():
    assert ValueTerm.Lit(5)._0 == 5


def test_valueterm_match():
    match ValueTerm.Lit(7):
        case ValueTerm.Lit(n):
            assert n == 7
        case _:
            raise AssertionError


def test_valueterm_recursive_neg():
    match ValueTerm.Neg(ValueTerm.Lit(3))._0:
        case ValueTerm.Lit(n):
            assert n == 3
        case _:
            raise AssertionError


def test_valueterm_sum():
    terms = ValueTerm.Sum([ValueTerm.Lit(1), ValueTerm.Lit(2)])._0
    assert len(terms) == 2
    assert terms[0]._0 == 1
    assert terms[1]._0 == 2


def test_valueterm_div():
    div = ValueTerm.Div(ValueTerm.Lit(6), ValueTerm.Lit(2))
    assert div._0._0 == 6
    assert div._1._0 == 2


def test_valuepredicate_rel():
    pred = ValuePredicate.Rel(ValueTerm.Var("h"), RelOp.Le, ValueTerm.Lit(3))
    match pred:
        case ValuePredicate.Rel(lhs, op, rhs):
            assert op == RelOp.Le
            assert rhs._0 == 3
        case _:
            raise AssertionError


def test_valuepredicate_mem():
    pred = ValuePredicate.Mem(ValueTerm.Var("x"), MemOp.In, {1, 2, 3})
    assert pred._1 == MemOp.In
    assert pred._2 == {1, 2, 3}


def test_valuepredicate_and():
    pred = ValuePredicate.And(
        [ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Lt, ValueTerm.Lit(2))]
    )
    assert len(pred._0) == 1


def test_valuepredicate_recursive_not():
    inner = ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Eq, ValueTerm.Lit(1))
    match ValuePredicate.Not(inner)._0:
        case ValuePredicate.Rel(_, op, _):
            assert op == RelOp.Eq
        case _:
            raise AssertionError


def test_valueast_lit():
    assert ValueAst.Lit(0)._0 == 0


def test_valueast_undetermined_match():
    match ValueAst.Undetermined():
        case ValueAst.Undetermined():
            pass
        case _:
            raise AssertionError


def test_valueast_litset():
    assert ValueAst.LitSet({1, 2, 3})._0 == {1, 2, 3}


def test_valueast_term_wraps_valueterm():
    match ValueAst.Term(ValueTerm.Var("h"))._0:
        case ValueTerm.Var(name):
            assert name == "h"
        case _:
            raise AssertionError


def test_valueast_predicate_wraps_valuepredicate():
    pred = ValuePredicate.Rel(ValueTerm.Var("h"), RelOp.Le, ValueTerm.Lit(3))
    match ValueAst.Predicate(pred)._0:
        case ValuePredicate.Rel(_, op, _):
            assert op == RelOp.Le
        case _:
            raise AssertionError
