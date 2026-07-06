from umol import (
    AtomAst,
    Element,
    ElementAst,
    IsotopeMassAst,
    MemOp,
    SpinStateAst,
    ValueAst,
    ValueTerm,
)


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


def test_atomast_from_element():
    match AtomAst.from_element(Element("C")).element:
        case ElementAst.Lit(e):
            assert e.symbol == "C"
        case _:
            raise AssertionError


def test_atomast_new_from_elementast():
    match AtomAst(ElementAst.Lit(Element("N"))).element:
        case ElementAst.Lit(e):
            assert e.symbol == "N"
        case _:
            raise AssertionError


def test_atomast_default_charge_undetermined():
    match AtomAst.from_element(Element("C")).charge:
        case ValueAst.Undetermined():
            pass
        case _:
            raise AssertionError


def test_atomast_with_charge():
    atom = AtomAst.from_element(Element("C")).with_charge(ValueAst.Lit(1))
    match atom.charge:
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_atomast_with_spin():
    spin = SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))
    atom = AtomAst.from_element(Element("C")).with_spin(spin)
    assert atom.spin.unpaired._0 == 1
    assert atom.spin.multiplicity._0 == 2


def test_atomast_eq():
    assert AtomAst.from_element(Element("C")) == AtomAst.from_element(Element("C"))
    assert AtomAst.from_element(Element("C")) != AtomAst.from_element(Element("N"))


def test_atomast_charge_nested_variant_readable():
    # A from_ast-built nested child (ValueTerm inside ValueAst.Term) must read back
    # as a proper variant from Python, not a base instance — regression for the
    # Py::new-vs-IntoPyObject bug.
    atom = AtomAst.from_element(Element("C")).with_charge(ValueAst.Term(ValueTerm.Var("h")))
    match atom.charge:
        case ValueAst.Term(term):
            match term:
                case ValueTerm.Var(name):
                    assert name == "h"
                case _:
                    raise AssertionError
        case _:
            raise AssertionError
