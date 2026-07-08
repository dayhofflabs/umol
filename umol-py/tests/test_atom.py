import pytest

from umol import (
    AtomAst,
    AtomConstraint,
    AtomConstraintKey,
    AtomConstraints,
    AtomId,
    Element,
    ElementAst,
    IsotopeMassAst,
    MemOp,
    MoleculeAst,
    SpinStateAst,
    ValueAst,
    ValueTerm,
)


def carbon_oxygen():
    return MoleculeAst.from_atoms(
        [AtomAst(Element("C")), AtomAst(Element("O"))]
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


def test_spinstateast_int_literals():
    spin = SpinStateAst(1, 2)
    assert spin.unpaired._0 == 1
    assert spin.multiplicity._0 == 2


def test_spinstateast_undetermined():
    spin = SpinStateAst(ValueAst.Undetermined(), ValueAst.Undetermined())
    match spin.unpaired:
        case ValueAst.Undetermined():
            pass
        case _:
            raise AssertionError


def test_atomast_new_from_element():
    match AtomAst(Element("C")).element:
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
    match AtomAst(Element("C")).charge:
        case ValueAst.Undetermined():
            pass
        case _:
            raise AssertionError


def test_atomast_replace_charge():
    atom = AtomAst(Element("C")).replace(charge=ValueAst.Lit(1))
    match atom.charge:
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_atomast_replace_spin():
    spin = SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))
    atom = AtomAst(Element("C")).replace(spin=spin)
    assert atom.spin.unpaired._0 == 1
    assert atom.spin.multiplicity._0 == 2


def test_atomast_new_from_element_kwargs():
    atom = AtomAst(Element("C"), charge=ValueAst.Lit(-1))
    match atom.charge:
        case ValueAst.Lit(n):
            assert n == -1
        case _:
            raise AssertionError


def test_atomast_new_bad_element_type():
    with pytest.raises(TypeError):
        AtomAst("C")


def test_atomast_new_kwargs():
    atom = AtomAst(ElementAst.Lit(Element("N")), charge=ValueAst.Lit(1))
    match atom.element:
        case ElementAst.Lit(e):
            assert e.symbol == "N"
        case _:
            raise AssertionError
    match atom.charge:
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_atomast_charge_int_literal():
    match AtomAst(Element("C"), charge=-1).charge:
        case ValueAst.Lit(n):
            assert n == -1
        case _:
            raise AssertionError


def test_atomast_isotope_mass_int_literal():
    match AtomAst(Element("C"), isotope_mass=13).isotope_mass:
        case IsotopeMassAst.Lit(n):
            assert n == 13
        case _:
            raise AssertionError


def test_atomast_replace_charge_int_literal():
    match AtomAst(Element("C")).replace(charge=1).charge:
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_atomast_constraints_empty():
    assert len(AtomAst(Element("C")).constraints) == 0


def test_atomast_constraints_kwarg():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraints([AtomConstraint.Valence(ValueAst.Lit(4))]),
    )
    assert len(atom.constraints) == 1
    match atom.constraints.get(AtomConstraintKey.Valence()):
        case AtomConstraint.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError


def test_atomview_constraints():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraints([AtomConstraint.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_atoms([atom])
    assert len(mol.atoms[AtomId(0)].constraints) == 1


def test_atomast_asdict():
    d = AtomAst(Element("C"), charge=ValueAst.Lit(-1)).asdict()
    assert set(d.keys()) == {
        "element",
        "isotope_mass",
        "charge",
        "implicit_hydrogens",
        "lone_pairs",
        "spin",
        "constraints",
    }
    match d["element"]:
        case ElementAst.Lit(e):
            assert e.symbol == "C"
        case _:
            raise AssertionError
    match d["charge"]:
        case ValueAst.Lit(n):
            assert n == -1
        case _:
            raise AssertionError


def test_atomast_eq():
    assert AtomAst(Element("C")) == AtomAst(Element("C"))
    assert AtomAst(Element("C")) != AtomAst(Element("N"))


def test_molecule_atoms_len():
    mol = carbon_oxygen()
    assert len(mol.atoms) == 2
    assert mol.atom_count == 2


def test_molecule_atoms_getitem():
    view = carbon_oxygen().atoms[AtomId(1)]
    assert view.id == AtomId(1)
    match view.element:
        case ElementAst.Lit(e):
            assert e.symbol == "O"
        case _:
            raise AssertionError


def test_molecule_atoms_getitem_out_of_range():
    with pytest.raises(IndexError):
        carbon_oxygen().atoms[AtomId(5)]


def test_molecule_atoms_iter():
    symbols = []
    for view in carbon_oxygen().atoms:
        match view.element:
            case ElementAst.Lit(e):
                symbols.append(e.symbol)
    assert symbols == ["C", "O"]


def test_atomview_charge_through_handle():
    atom = AtomAst(Element("C")).replace(charge=ValueAst.Lit(-1))
    mol = MoleculeAst.from_atoms([atom])
    match mol.atoms[AtomId(0)].charge:
        case ValueAst.Lit(n):
            assert n == -1
        case _:
            raise AssertionError


def test_atomid_index_and_repr():
    aid = AtomId(3)
    assert aid.index == 3
    assert repr(aid) == "AtomId(3)"


def test_atomid_eq_hash():
    assert AtomId(3) == AtomId(3)
    assert AtomId(3) != AtomId(4)
    assert len({AtomId(3), AtomId(3)}) == 1


def test_atomast_charge_nested_variant_readable():
    # A from_ast-built nested child (ValueTerm inside ValueAst.Term) must read back
    # as a proper variant from Python, not a base instance — regression for the
    # Py::new-vs-IntoPyObject bug.
    atom = AtomAst(Element("C")).replace(charge=ValueAst.Term(ValueTerm.Var("h")))
    match atom.charge:
        case ValueAst.Term(term):
            match term:
                case ValueTerm.Var(name):
                    assert name == "h"
                case _:
                    raise AssertionError
        case _:
            raise AssertionError
