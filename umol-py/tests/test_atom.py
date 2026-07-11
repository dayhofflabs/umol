import pytest

from umol import (
    AtomAst,
    AtomConstraintAst,
    AtomConstraintKey,
    AtomConstraintsAst,
    Element,
    ElementAst,
    IsotopeMassAst,
    MemOp,
    MoleculeAst,
    ParseError,
    SpinStateAst,
    ValueAst,
    ValueTerm,
)


def carbon_oxygen():
    return MoleculeAst.from_atoms_and_bonds(
        [AtomAst(Element("C")), AtomAst(Element("O"))]
    )


def test_elementast_lit():
    assert ElementAst.Lit(Element("C")) == ElementAst.Lit(Element("C"))


def test_elementast_undetermined_match():
    assert ElementAst.Undetermined() == ElementAst.Undetermined()


def test_elementast_litset():
    members = ElementAst.LitSet({Element("C"), Element("N")})._0
    assert members == {Element("C"), Element("N")}


def test_elementast_notset():
    assert ElementAst.NotSet({Element("O")}) == ElementAst.NotSet({Element("O")})


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
    assert IsotopeMassAst.Natural() == IsotopeMassAst.Natural()


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
    assert spin.unpaired == ValueAst.Undetermined()


def test_atomast_new_from_element():
    assert AtomAst(Element("C")).element == ElementAst.Lit(Element("C"))


def test_atomast_new_from_elementast():
    assert AtomAst(ElementAst.Lit(Element("N"))).element == ElementAst.Lit(Element("N"))


def test_atomast_default_charge_undetermined():
    assert AtomAst(Element("C")).charge == ValueAst.Undetermined()


def test_atomast_set_charge():
    atom = AtomAst(Element("C"))
    atom.charge = ValueAst.Lit(1)
    assert atom.charge == ValueAst.Lit(1)


def test_atomast_set_spin():
    spin = SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))
    atom = AtomAst(Element("C"))
    atom.spin = spin
    assert atom.spin.unpaired._0 == 1
    assert atom.spin.multiplicity._0 == 2


def test_atomast_new_from_element_kwargs():
    atom = AtomAst(Element("C"), charge=ValueAst.Lit(-1))
    assert atom.charge == ValueAst.Lit(-1)


def test_atomast_new_bad_element_type():
    with pytest.raises(TypeError):
        AtomAst("C")


def test_atomast_new_kwargs():
    atom = AtomAst(ElementAst.Lit(Element("N")), charge=ValueAst.Lit(1))
    assert atom.element == ElementAst.Lit(Element("N"))
    assert atom.charge == ValueAst.Lit(1)


def test_atomast_charge_int_literal():
    assert AtomAst(Element("C"), charge=-1).charge == ValueAst.Lit(-1)


def test_atomast_isotope_mass_int_literal():
    assert AtomAst(Element("C"), isotope_mass=13).isotope_mass == IsotopeMassAst.Lit(13)


def test_atomast_set_charge_int_literal():
    atom = AtomAst(Element("C"))
    atom.charge = 1
    assert atom.charge == ValueAst.Lit(1)


def test_atomast_constraints_empty():
    assert len(AtomAst(Element("C")).constraints) == 0


def test_atomast_constraints_kwarg():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    assert len(atom.constraints) == 1
    assert atom.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_atomview_constraints():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_atoms_and_bonds([atom])
    assert len(mol.atoms[0].constraints) == 1


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
    assert d["element"] == ElementAst.Lit(Element("C"))
    assert d["charge"] == ValueAst.Lit(-1)


def test_atomast_asdict_constraints():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    constraints = atom.asdict()["constraints"]
    assert isinstance(constraints, dict)
    assert set(constraints.keys()) == {"valence"}
    assert constraints["valence"] == ValueAst.Lit(4)


def test_atomast_eq():
    assert AtomAst(Element("C")) == AtomAst(Element("C"))
    assert AtomAst(Element("C")) != AtomAst(Element("N"))


def test_molecule_atoms_len():
    assert len(carbon_oxygen().atoms) == 2


def test_molecule_atoms_getitem():
    view = carbon_oxygen().atoms[1]
    assert view.id == 1
    assert view.element == ElementAst.Lit(Element("O"))


def test_molecule_atoms_getitem_out_of_range():
    with pytest.raises(IndexError):
        carbon_oxygen().atoms[5]


def test_molecule_atoms_negative_index():
    mol = carbon_oxygen()
    assert mol.atoms[-1].id == 1
    assert mol.atoms[-2].id == 0
    mol.atoms[-1] = AtomAst(Element("N"))
    assert mol.atoms[1].element == ElementAst.Lit(Element("N"))
    with pytest.raises(IndexError):
        mol.atoms[-3]


def test_molecule_atoms_setitem():
    mol = carbon_oxygen()
    mol.atoms[0] = AtomAst(Element("N"))
    assert mol.atoms[0].element == ElementAst.Lit(Element("N"))


def test_molecule_atoms_setitem_out_of_range():
    with pytest.raises(IndexError):
        carbon_oxygen().atoms[5] = AtomAst(Element("N"))


def test_molecule_atoms_iter():
    symbols = []
    for view in carbon_oxygen().atoms:
        match view.element:
            case ElementAst.Lit(e):
                symbols.append(e.symbol)
    assert symbols == ["C", "O"]


def test_atomview_charge_through_handle():
    atom = AtomAst(Element("C"))
    atom.charge = ValueAst.Lit(-1)
    mol = MoleculeAst.from_atoms_and_bonds([atom])
    assert mol.atoms[0].charge == ValueAst.Lit(-1)


def test_atomast_charge_nested_variant_readable():
    # A from_ast-built nested child (ValueTerm inside ValueAst.Term) must read back
    # as a proper variant from Python, not a base instance — regression for the
    # Py::new-vs-IntoPyObject bug.
    atom = AtomAst(Element("C"))
    atom.charge = ValueAst.Term(ValueTerm.Var("h"))
    match atom.charge:
        case ValueAst.Term(term):
            match term:
                case ValueTerm.Var(name):
                    assert name == "h"
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


@pytest.mark.parametrize("dsl", ["C", "N#c+", "C#v4", "O#n2", "C#R(6)"])
def test_atomast_parse(dsl):
    atom = AtomAst.parse(dsl)
    assert str(atom) == dsl
    assert repr(atom) == f"AtomAst.parse('{dsl}')"


def test_atomast_parse_error():
    with pytest.raises(ParseError):
        AtomAst.parse("Zz##")


def test_atomview_set_charge():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].charge = ValueAst.Lit(-1)
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.atoms[0].charge == ValueAst.Lit(-1)


def test_atomview_set_charge_int_literal():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].charge = -1
    assert mol.atoms[0].charge == ValueAst.Lit(-1)


def test_atomview_set_element():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].element = Element("N")
    assert mol.atoms[0].element == ElementAst.Lit(Element("N"))


def test_atomview_set_isotope_mass():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].isotope_mass = 13
    assert mol.atoms[0].isotope_mass == IsotopeMassAst.Lit(13)


def test_atomview_set_implicit_hydrogens():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].implicit_hydrogens = ValueAst.Lit(3)
    assert mol.atoms[0].implicit_hydrogens == ValueAst.Lit(3)


def test_atomview_set_lone_pairs():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("O"))])
    mol.atoms[0].lone_pairs = ValueAst.Lit(2)
    assert mol.atoms[0].lone_pairs == ValueAst.Lit(2)


def test_atomview_set_spin():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].spin = SpinStateAst(1, 2)
    spin = mol.atoms[0].spin
    assert spin.unpaired._0 == 1
    assert spin.multiplicity._0 == 2


def test_elementast_as_lit():
    assert ElementAst.Lit(Element("C")).as_lit().symbol == "C"
    assert ElementAst.Undetermined().as_lit() is None


def test_isotopemassast_as_lit():
    assert IsotopeMassAst.Lit(13).as_lit() == 13
    assert IsotopeMassAst.Natural().as_lit() is None


def test_valueast_as_lit():
    assert ValueAst.Lit(4).as_lit() == 4
    assert ValueAst.Undetermined().as_lit() is None


def test_elementast_eq_hash_repr():
    assert ElementAst.Lit(Element("C")) == ElementAst.Lit(Element("C"))
    assert ElementAst.Lit(Element("C")) != ElementAst.Lit(Element("N"))
    assert len({ElementAst.Lit(Element("C")), ElementAst.Lit(Element("C"))}) == 1
    assert repr(ElementAst.Lit(Element("C"))) == "ElementAst.Lit(Element('C'))"


def test_isotopemassast_eq_repr():
    assert IsotopeMassAst.Lit(13) == IsotopeMassAst.Lit(13)
    assert IsotopeMassAst.Lit(13) != IsotopeMassAst.Natural()
    assert repr(IsotopeMassAst.Lit(13)) == "IsotopeMassAst.Lit(13)"


def test_spinstateast_eq_repr():
    assert SpinStateAst(1, 2) == SpinStateAst(1, 2)
    assert SpinStateAst(1, 2) != SpinStateAst(1, 3)
    assert repr(SpinStateAst(1, 2)) == "SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))"


def test_atomview_repr():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C")), AtomAst(Element("O"))])
    assert repr(mol.atoms[0]) == "AtomView(id=0)"
    assert repr(mol.atoms) == "AtomViews(len=2)"


def test_atomview_asdict():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"), charge=ValueAst.Lit(-1))])
    d = mol.atoms[0].asdict()
    assert set(d.keys()) == {
        "element",
        "isotope_mass",
        "charge",
        "implicit_hydrogens",
        "lone_pairs",
        "spin",
        "constraints",
    }
    assert d["element"] == ElementAst.Lit(Element("C"))
    assert d["charge"] == ValueAst.Lit(-1)
    assert isinstance(d["constraints"], dict)


def test_atomview_set_constraints():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    mol.atoms[0].constraints = AtomConstraintsAst([AtomConstraintAst.Degree(ValueAst.Lit(2))])
    assert len(mol.atoms[0].constraints) == 1
    assert mol.atoms[0].constraints.get(AtomConstraintKey.Degree()) == AtomConstraintAst.Degree(
        ValueAst.Lit(2)
    )


def test_atomast_set_constraints():
    atom = AtomAst(Element("C"))
    atom.constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert len(atom.constraints) == 1
    atom.constraints = AtomConstraintsAst([])
    assert not atom.constraints
