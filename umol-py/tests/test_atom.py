import pytest

from umol import (
    AtomAst,
    AtomConstraintAst,
    AtomConstraintKey,
    AtomConstraintsAst,
    AtomUpdate,
    Element,
    ElementForm,
    IsotopeMass,
    IsotopeMassForm,
    MemOp,
    MoleculeAst,
    ParseError,
    UnpairedElectronsForm,
    UnpairedElectronsUpdate,
    NumForm,
    ArithExpr,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            AtomUpdate(),
            (None, None, None, None, None, UnpairedElectronsUpdate(), AtomConstraintsAst([])),
        ),
        (
            AtomUpdate(
                element=Element("N"),
                isotope_mass=15,
                charge=1,
                implicit_hydrogens=2,
                lone_pairs=1,
                unpaired_electrons=UnpairedElectronsUpdate(count=1),
                constraints=AtomConstraintsAst(
                    [AtomConstraintAst.Valence(NumForm.Lit(3))]
                ),
            ),
            (
                ElementForm.Lit(Element("N")),
                IsotopeMassForm.Lit(15),
                NumForm.Lit(1),
                NumForm.Lit(2),
                NumForm.Lit(1),
                UnpairedElectronsUpdate(count=1),
                AtomConstraintsAst([AtomConstraintAst.Valence(NumForm.Lit(3))]),
            ),
        ),
        (
            AtomUpdate(
                constraints=AtomConstraintsAst(
                    [AtomConstraintAst.Valence(NumForm.Undetermined())]
                )
            ),
            (
                None,
                None,
                None,
                None,
                None,
                UnpairedElectronsUpdate(),
                AtomConstraintsAst(
                    [AtomConstraintAst.Valence(NumForm.Undetermined())]
                ),
            ),
        ),
    ],
)
def test_atom_update(update, expected):
    assert (
        update.element,
        update.isotope_mass,
        update.charge,
        update.implicit_hydrogens,
        update.lone_pairs,
        update.unpaired_electrons,
        update.constraints,
    ) == expected


@pytest.mark.parametrize(
    ("dsl", "canonical"),
    [
        ("", ""),
        ("#c-1", "#c-"),
        ("*#i*#c*#h*#n*#u*#s*", "*#i*#c*#h*#n*#u*#s*"),
        ("#v*", "#v*"),
    ],
)
def test_atom_update_parse(dsl, canonical):
    update = AtomUpdate.parse(dsl)
    assert str(update) == canonical
    assert repr(update) == f"AtomUpdate.parse('{canonical}')"
    assert eval(repr(update)) == update


def test_atom_update_parse_error():
    with pytest.raises(ParseError):
        AtomUpdate.parse("#h1#h2")


def carbon_oxygen():
    return MoleculeAst.from_entries(
        [AtomAst(Element("C")), AtomAst(Element("O"))]
    )


def test_elementast_lit():
    assert ElementForm.Lit(Element("C")) == ElementForm.Lit(Element("C"))


def test_elementast_undetermined_match():
    assert ElementForm.Undetermined() == ElementForm.Undetermined()


def test_elementast_litset():
    members = ElementForm.LitSet({Element("C"), Element("N")})._0
    assert members == {Element("C"), Element("N")}


def test_elementast_notset():
    assert ElementForm.NotSet({Element("O")}) == ElementForm.NotSet({Element("O")})


def test_elementast_var_free():
    var = ElementForm.Var("x", None)
    assert var._0 == "x"
    assert var._1 is None


def test_elementast_var_restricted():
    var = ElementForm.Var("y", (MemOp.In, {Element("C")}))
    op, members = var._1
    assert op == MemOp.In
    assert members == {Element("C")}


def test_isotopemassast_natural_match():
    assert IsotopeMassForm.Natural() == IsotopeMassForm.Natural()


def test_isotopemassast_lit():
    assert IsotopeMassForm.Lit(13)._0 == 13


def test_isotopemassast_litset():
    assert IsotopeMassForm.LitSet({12, 13, 14})._0 == {12, 13, 14}


def test_isotopemassast_var_free():
    var = IsotopeMassForm.Var("x", None)
    assert var._0 == "x"
    assert var._1 is None


def test_isotopemassast_var_restricted():
    assert IsotopeMassForm.Var("y", {12, 13})._1 == {12, 13}


def test_unpairedelectronsast_fields():
    unpaired_electrons = UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2))
    assert unpaired_electrons.count._0 == 1
    assert unpaired_electrons.multiplicity._0 == 2


def test_unpairedelectronsast_int_literals():
    unpaired_electrons = UnpairedElectronsForm(count=2, multiplicity=2)
    assert unpaired_electrons.count._0 == 2
    assert unpaired_electrons.multiplicity._0 == 2


def test_unpairedelectronsast_undetermined():
    unpaired_electrons = UnpairedElectronsForm(
        NumForm.Undetermined(), NumForm.Undetermined()
    )
    assert unpaired_electrons.count == NumForm.Undetermined()


def test_atomast_new_from_element():
    assert AtomAst(Element("C")).element == ElementForm.Lit(Element("C"))


def test_atomast_new_from_elementast():
    assert AtomAst(ElementForm.Lit(Element("N"))).element == ElementForm.Lit(Element("N"))


def test_atomast_default_charge_undetermined():
    assert AtomAst(Element("C")).charge == NumForm.Undetermined()


def test_atomast_set_charge():
    atom = AtomAst(Element("C"))
    atom.charge = NumForm.Lit(1)
    assert atom.charge == NumForm.Lit(1)


def test_atomast_set_unpaired_electrons():
    unpaired_electrons = UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2))
    atom = AtomAst(Element("C"))
    atom.unpaired_electrons = unpaired_electrons
    assert atom.unpaired_electrons.count._0 == 1
    assert atom.unpaired_electrons.multiplicity._0 == 2


def test_atomast_new_from_element_kwargs():
    atom = AtomAst(Element("C"), charge=NumForm.Lit(-1))
    assert atom.charge == NumForm.Lit(-1)


def test_atomast_new_bad_element_type():
    with pytest.raises(TypeError):
        AtomAst("C")


def test_atomast_new_kwargs():
    atom = AtomAst(
        ElementForm.Lit(Element("N")),
        charge=NumForm.Lit(1),
        unpaired_electrons=UnpairedElectronsForm(1, 2),
    )
    assert atom.element == ElementForm.Lit(Element("N"))
    assert atom.charge == NumForm.Lit(1)
    assert atom.unpaired_electrons == UnpairedElectronsForm(1, 2)


def test_atomast_charge_int_literal():
    assert AtomAst(Element("C"), charge=-1).charge == NumForm.Lit(-1)


def test_atomast_isotope_mass_int_literal():
    assert AtomAst(Element("C"), isotope_mass=13).isotope_mass == IsotopeMassForm.Lit(13)


def test_atomast_set_charge_int_literal():
    atom = AtomAst(Element("C"))
    atom.charge = 1
    assert atom.charge == NumForm.Lit(1)


def test_atomast_constraints_empty():
    assert len(AtomAst(Element("C")).constraints) == 0


def test_atomast_constraints_kwarg():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(NumForm.Lit(4))]),
    )
    assert len(atom.constraints) == 1
    assert atom.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        NumForm.Lit(4)
    )


def test_atomview_constraints():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(NumForm.Lit(4))]),
    )
    mol = MoleculeAst.from_entries([atom])
    assert len(mol.atoms[0].constraints) == 1


def test_atomast_asdict():
    d = AtomAst(Element("C"), charge=NumForm.Lit(-1)).asdict()
    assert set(d.keys()) == {
        "element",
        "isotope_mass",
        "charge",
        "implicit_hydrogens",
        "lone_pairs",
        "unpaired_electrons",
        "constraints",
    }
    assert d["element"] == ElementForm.Lit(Element("C"))
    assert d["charge"] == NumForm.Lit(-1)


def test_atomast_asdict_constraints():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(NumForm.Lit(4))]),
    )
    constraints = atom.asdict()["constraints"]
    assert isinstance(constraints, dict)
    assert set(constraints.keys()) == {"valence"}
    assert constraints["valence"] == NumForm.Lit(4)


def test_atomast_eq():
    assert AtomAst(Element("C")) == AtomAst(Element("C"))
    assert AtomAst(Element("C")) != AtomAst(Element("N"))


def test_molecule_atoms_len():
    assert len(carbon_oxygen().atoms) == 2


def test_molecule_atoms_getitem():
    view = carbon_oxygen().atoms[1]
    assert view.id == 1
    assert view.element == ElementForm.Lit(Element("O"))


def test_molecule_atoms_getitem_out_of_range():
    with pytest.raises(IndexError):
        carbon_oxygen().atoms[5]


def test_molecule_atoms_negative_index():
    mol = carbon_oxygen()
    assert mol.atoms[-1].id == 1
    assert mol.atoms[-2].id == 0
    mol.atoms[-1] = AtomAst(Element("N"))
    assert mol.atoms[1].element == ElementForm.Lit(Element("N"))
    with pytest.raises(IndexError):
        mol.atoms[-3]


def test_molecule_atoms_setitem():
    mol = carbon_oxygen()
    mol.atoms[0] = AtomAst(Element("N"))
    assert mol.atoms[0].element == ElementForm.Lit(Element("N"))


def test_molecule_atoms_setitem_out_of_range():
    with pytest.raises(IndexError):
        carbon_oxygen().atoms[5] = AtomAst(Element("N"))


def test_molecule_atoms_iter():
    symbols = []
    for view in carbon_oxygen().atoms:
        match view.element:
            case ElementForm.Lit(e):
                symbols.append(e.symbol)
    assert symbols == ["C", "O"]


def test_atomview_charge_through_handle():
    atom = AtomAst(Element("C"))
    atom.charge = NumForm.Lit(-1)
    mol = MoleculeAst.from_entries([atom])
    assert mol.atoms[0].charge == NumForm.Lit(-1)


def test_atomast_charge_nested_variant_readable():
    # A conversion-built nested child (ArithExpr inside NumForm.ArithExpr) must read back
    # as a proper variant from Python, not a base instance — regression for the
    # Py::new-vs-IntoPyObject bug.
    atom = AtomAst(Element("C"))
    atom.charge = NumForm.ArithExpr(ArithExpr.Var("h"))
    match atom.charge:
        case NumForm.ArithExpr(term):
            match term:
                case ArithExpr.Var(name):
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
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].charge = NumForm.Lit(-1)
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.atoms[0].charge == NumForm.Lit(-1)


def test_atomview_set_charge_int_literal():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].charge = -1
    assert mol.atoms[0].charge == NumForm.Lit(-1)


def test_atomview_set_element():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].element = Element("N")
    assert mol.atoms[0].element == ElementForm.Lit(Element("N"))


def test_atomview_set_isotope_mass():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].isotope_mass = 13
    assert mol.atoms[0].isotope_mass == IsotopeMassForm.Lit(13)


def test_atomview_set_implicit_hydrogens():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].implicit_hydrogens = NumForm.Lit(3)
    assert mol.atoms[0].implicit_hydrogens == NumForm.Lit(3)


def test_atomview_set_lone_pairs():
    mol = MoleculeAst.from_entries([AtomAst(Element("O"))])
    mol.atoms[0].lone_pairs = NumForm.Lit(2)
    assert mol.atoms[0].lone_pairs == NumForm.Lit(2)


def test_atomview_set_unpaired_electrons():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].unpaired_electrons = UnpairedElectronsForm(1, 2)
    unpaired_electrons = mol.atoms[0].unpaired_electrons
    assert unpaired_electrons.count._0 == 1
    assert unpaired_electrons.multiplicity._0 == 2


def test_elementast_as_lit():
    assert ElementForm.Lit(Element("C")).as_lit().symbol == "C"
    assert ElementForm.Undetermined().as_lit() is None


def test_isotopemassast_as_lit():
    assert IsotopeMassForm.Lit(13).as_lit() == IsotopeMass.MassNumber(13)
    assert IsotopeMassForm.Natural().as_lit() == IsotopeMass.Natural()
    assert IsotopeMassForm.Undetermined().as_lit() is None


def test_valueast_as_lit():
    assert NumForm.Lit(4).as_lit() == 4
    assert NumForm.Undetermined().as_lit() is None


def test_elementast_eq_hash_repr():
    assert ElementForm.Lit(Element("C")) == ElementForm.Lit(Element("C"))
    assert ElementForm.Lit(Element("C")) != ElementForm.Lit(Element("N"))
    assert len({ElementForm.Lit(Element("C")), ElementForm.Lit(Element("C"))}) == 1
    assert repr(ElementForm.Lit(Element("C"))) == "ElementForm.Lit(Element('C'))"


def test_isotopemassast_eq_repr():
    assert IsotopeMassForm.Lit(13) == IsotopeMassForm.Lit(13)
    assert IsotopeMassForm.Lit(13) != IsotopeMassForm.Natural()
    assert repr(IsotopeMassForm.Lit(13)) == "IsotopeMassForm.Lit(13)"


def test_unpairedelectronsast_eq_repr():
    assert UnpairedElectronsForm(1, 2) == UnpairedElectronsForm(1, 2)
    assert UnpairedElectronsForm(1, 2) != UnpairedElectronsForm(1, 3)
    assert repr(UnpairedElectronsForm(1, 2)) == (
        "UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2))"
    )


def test_atomview_repr():
    mol = MoleculeAst.from_entries([AtomAst(Element("C")), AtomAst(Element("O"))])
    assert repr(mol.atoms[0]) == "AtomView(id=0)"
    assert repr(mol.atoms) == "AtomViews(len=2)"


def test_atomview_asdict():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"), charge=NumForm.Lit(-1))])
    d = mol.atoms[0].asdict()
    assert set(d.keys()) == {
        "element",
        "isotope_mass",
        "charge",
        "implicit_hydrogens",
        "lone_pairs",
        "unpaired_electrons",
        "constraints",
    }
    assert d["element"] == ElementForm.Lit(Element("C"))
    assert d["charge"] == NumForm.Lit(-1)
    assert isinstance(d["constraints"], dict)


def test_atomview_set_constraints():
    mol = MoleculeAst.from_entries([AtomAst(Element("C"))])
    mol.atoms[0].constraints = AtomConstraintsAst([AtomConstraintAst.Degree(NumForm.Lit(2))])
    assert len(mol.atoms[0].constraints) == 1
    assert mol.atoms[0].constraints.get(AtomConstraintKey.Degree()) == AtomConstraintAst.Degree(
        NumForm.Lit(2)
    )


def test_atomast_set_constraints():
    atom = AtomAst(Element("C"))
    atom.constraints = AtomConstraintsAst([AtomConstraintAst.Valence(NumForm.Lit(4))])
    assert len(atom.constraints) == 1
    atom.constraints = AtomConstraintsAst([])
    assert not atom.constraints
