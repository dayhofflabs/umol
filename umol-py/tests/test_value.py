from umol import MemOp, RelOp, ValueAst, ValuePredicate, ValueTerm


def test_relop_eq():
    assert RelOp.Lt == RelOp.Lt
    assert RelOp.Lt != RelOp.Ge


def test_relop_hashable():
    assert len({RelOp.Lt, RelOp.Lt, RelOp.Ge}) == 2


def test_relop_match():
    match RelOp.Ne:
        case RelOp.Ne:
            pass
        case _:
            raise AssertionError


def test_memop_eq():
    assert MemOp.In == MemOp.In
    assert MemOp.In != MemOp.NotIn


def test_memop_hashable():
    assert len({MemOp.In, MemOp.In, MemOp.NotIn}) == 2


def test_valueterm_lit():
    assert ValueTerm.Lit(5)._0 == 5


def test_valueterm_match():
    assert ValueTerm.Lit(7) == ValueTerm.Lit(7)


def test_valueterm_recursive_neg():
    assert ValueTerm.Neg(ValueTerm.Lit(3))._0 == ValueTerm.Lit(3)


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
    assert pred == ValuePredicate.Rel(ValueTerm.Var("h"), RelOp.Le, ValueTerm.Lit(3))


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
    assert ValuePredicate.Not(inner)._0 == inner


def test_valueast_lit():
    assert ValueAst.Lit(0)._0 == 0


def test_valueast_undetermined_match():
    assert ValueAst.Undetermined() == ValueAst.Undetermined()


def test_valueast_litset():
    assert ValueAst.LitSet({1, 2, 3})._0 == {1, 2, 3}


def test_valueast_term_wraps_valueterm():
    assert ValueAst.Term(ValueTerm.Var("h"))._0 == ValueTerm.Var("h")


def test_valueast_predicate_wraps_valuepredicate():
    pred = ValuePredicate.Rel(ValueTerm.Var("h"), RelOp.Le, ValueTerm.Lit(3))
    assert ValueAst.Predicate(pred)._0 == pred


def test_valueast_eq():
    assert ValueAst.Lit(1) == ValueAst.Lit(1)
    assert ValueAst.Lit(1) != ValueAst.Lit(2)
    assert ValueAst.Lit(1) != 5


def test_valueast_hash():
    assert len({ValueAst.Lit(1), ValueAst.Lit(1)}) == 1
    d = {ValueAst.Lit(1): "a"}
    assert d[ValueAst.Lit(1)] == "a"


def test_valueast_repr():
    assert repr(ValueAst.Lit(1)) == "ValueAst.Lit(1)"
    assert repr(ValueAst.Undetermined()) == "ValueAst.Undetermined()"
    x = ValueAst.Term(ValueTerm.Var("h"))
    assert eval(repr(x), {"ValueAst": ValueAst, "ValueTerm": ValueTerm}) == x


def test_valueterm_eq_repr():
    assert ValueTerm.Lit(1) == ValueTerm.Lit(1)
    assert ValueTerm.Lit(1) != ValueTerm.Var("x")
    assert repr(ValueTerm.Sum([ValueTerm.Lit(1), ValueTerm.Lit(2)])) == (
        "ValueTerm.Sum([ValueTerm.Lit(1), ValueTerm.Lit(2)])"
    )


def test_valuepredicate_eq_repr():
    a = ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Le, ValueTerm.Lit(2))
    b = ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Le, ValueTerm.Lit(2))
    assert a == b
    assert a != ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Ge, ValueTerm.Lit(2))
    assert repr(a) == "ValuePredicate.Rel(ValueTerm.Lit(1), RelOp.Le, ValueTerm.Lit(2))"
