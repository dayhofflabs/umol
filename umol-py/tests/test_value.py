from umol import MemOp, RelOp, ValueAst, PredExpr, ArithExpr


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


def test_arith_expr_lit():
    assert ArithExpr.Lit(5)._0 == 5


def test_arith_expr_match():
    assert ArithExpr.Lit(7) == ArithExpr.Lit(7)


def test_arith_expr_recursive_neg():
    assert ArithExpr.Neg(ArithExpr.Lit(3))._0 == ArithExpr.Lit(3)


def test_arith_expr_sum():
    terms = ArithExpr.Sum([ArithExpr.Lit(1), ArithExpr.Lit(2)])._0
    assert len(terms) == 2
    assert terms[0]._0 == 1
    assert terms[1]._0 == 2


def test_arith_expr_div():
    div = ArithExpr.Div(ArithExpr.Lit(6), ArithExpr.Lit(2))
    assert div._0._0 == 6
    assert div._1._0 == 2


def test_pred_expr_rel():
    pred = PredExpr.Rel(ArithExpr.Var("h"), RelOp.Le, ArithExpr.Lit(3))
    assert pred == PredExpr.Rel(ArithExpr.Var("h"), RelOp.Le, ArithExpr.Lit(3))


def test_pred_expr_mem():
    pred = PredExpr.Mem(ArithExpr.Var("x"), MemOp.In, {1, 2, 3})
    assert pred._1 == MemOp.In
    assert pred._2 == {1, 2, 3}


def test_pred_expr_and():
    pred = PredExpr.And(
        [PredExpr.Rel(ArithExpr.Lit(1), RelOp.Lt, ArithExpr.Lit(2))]
    )
    assert len(pred._0) == 1


def test_pred_expr_recursive_not():
    inner = PredExpr.Rel(ArithExpr.Lit(1), RelOp.Eq, ArithExpr.Lit(1))
    assert PredExpr.Not(inner)._0 == inner


def test_valueast_lit():
    assert ValueAst.Lit(0)._0 == 0


def test_valueast_undetermined_match():
    assert ValueAst.Undetermined() == ValueAst.Undetermined()


def test_valueast_litset():
    assert ValueAst.LitSet({1, 2, 3})._0 == {1, 2, 3}


def test_valueast_arith_expr():
    assert ValueAst.ArithExpr(ArithExpr.Var("h"))._0 == ArithExpr.Var("h")


def test_valueast_pred_expr():
    pred = PredExpr.Rel(ArithExpr.Var("h"), RelOp.Le, ArithExpr.Lit(3))
    assert ValueAst.PredExpr(pred)._0 == pred


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
    x = ValueAst.ArithExpr(ArithExpr.Var("h"))
    assert eval(repr(x), {"ValueAst": ValueAst, "ArithExpr": ArithExpr}) == x


def test_arith_expr_eq_repr():
    assert ArithExpr.Lit(1) == ArithExpr.Lit(1)
    assert ArithExpr.Lit(1) != ArithExpr.Var("x")
    assert repr(ArithExpr.Sum([ArithExpr.Lit(1), ArithExpr.Lit(2)])) == (
        "ArithExpr.Sum([ArithExpr.Lit(1), ArithExpr.Lit(2)])"
    )


def test_pred_expr_eq_repr():
    a = PredExpr.Rel(ArithExpr.Lit(1), RelOp.Le, ArithExpr.Lit(2))
    b = PredExpr.Rel(ArithExpr.Lit(1), RelOp.Le, ArithExpr.Lit(2))
    assert a == b
    assert a != PredExpr.Rel(ArithExpr.Lit(1), RelOp.Ge, ArithExpr.Lit(2))
    assert repr(a) == "PredExpr.Rel(ArithExpr.Lit(1), RelOp.Le, ArithExpr.Lit(2))"
