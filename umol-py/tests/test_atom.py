from umol import Element, ElementAst, IsotopeMassAst, MemOp, SpinStateAst, ValueAst


def test_elementast_lit():
    match ElementAst.Lit(Element("C")):
        case ElementAst.Lit(e):
            assert e.symbol == "C"
        case _:
            raise AssertionError


def test_elementast_undetermined_match():
    match ElementAst.Undetermined():
        case ElementAst.Undetermined():
            pass
        case _:
            raise AssertionError


def test_elementast_litset():
    members = ElementAst.LitSet({Element("C"), Element("N")})._0
    assert members == {Element("C"), Element("N")}


def test_elementast_notset():
    match ElementAst.NotSet({Element("O")}):
        case ElementAst.NotSet(members):
            assert members == {Element("O")}
        case _:
            raise AssertionError


def test_elementast_var_free():
    var = ElementAst.Var("x", None)
    assert var._0 == "x"
    assert var._1 is None


def test_elementast_var_restricted():
    var = ElementAst.Var("y", (MemOp.In, {Element("C")}))
    op, members = var._1
    assert op == MemOp.In
    assert members == {Element("C")}


def test_isotopemassast_natural_match():
    match IsotopeMassAst.Natural():
        case IsotopeMassAst.Natural():
            pass
        case _:
            raise AssertionError


def test_isotopemassast_lit():
    assert IsotopeMassAst.Lit(13)._0 == 13


def test_isotopemassast_litset():
    assert IsotopeMassAst.LitSet({12, 13, 14})._0 == {12, 13, 14}


def test_isotopemassast_var_free():
    var = IsotopeMassAst.Var("x", None)
    assert var._0 == "x"
    assert var._1 is None


def test_isotopemassast_var_restricted():
    assert IsotopeMassAst.Var("y", {12, 13})._1 == {12, 13}


def test_spinstateast_fields():
    spin = SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))
    assert spin.unpaired._0 == 1
    assert spin.multiplicity._0 == 2


def test_spinstateast_undetermined():
    spin = SpinStateAst(ValueAst.Undetermined(), ValueAst.Undetermined())
    match spin.unpaired:
        case ValueAst.Undetermined():
            pass
        case _:
            raise AssertionError
