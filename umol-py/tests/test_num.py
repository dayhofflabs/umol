from umol import MemOp, RelOp, NumForm, PredExpr, ArithExpr


def test_relop_eq():
    assert RelOp.Lt == RelOp.Lt
    assert RelOp.Lt != RelOp.Ge


def test_relop_hashable():
    assert len({RelOp.Lt, RelOp.Lt, RelOp.Ge}) == 2


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


def test_num_form_lit():
    assert NumForm.Lit(0)._0 == 0


def test_num_form_undetermined_match():
    assert NumForm.Undetermined() == NumForm.Undetermined()


def test_num_form_litset():
    assert NumForm.LitSet({1, 2, 3})._0 == {1, 2, 3}


def test_num_form_arith_expr():
    assert NumForm.ArithExpr(ArithExpr.Var("h"))._0 == ArithExpr.Var("h")


def test_num_form_pred_expr():
    pred = PredExpr.Rel(ArithExpr.Var("h"), RelOp.Le, ArithExpr.Lit(3))
    assert NumForm.PredExpr(pred)._0 == pred


def test_num_form_eq():
    assert NumForm.Lit(1) == NumForm.Lit(1)
    assert NumForm.Lit(1) != NumForm.Lit(2)
    assert NumForm.Lit(1) != 5


def test_num_form_hash():
    assert len({NumForm.Lit(1), NumForm.Lit(1)}) == 1
    d = {NumForm.Lit(1): "a"}
    assert d[NumForm.Lit(1)] == "a"


def test_num_form_repr():
    assert repr(NumForm.Lit(1)) == "NumForm.Lit(1)"
    assert repr(NumForm.Undetermined()) == "NumForm.Undetermined()"
    x = NumForm.ArithExpr(ArithExpr.Var("h"))
    assert eval(repr(x), {"NumForm": NumForm, "ArithExpr": ArithExpr}) == x


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
