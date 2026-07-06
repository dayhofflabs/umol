from umol import MemOp, RelOp, ValueTerm


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
