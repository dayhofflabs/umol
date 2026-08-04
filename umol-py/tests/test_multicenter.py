import pytest

from umol import (
    AtomAst,
    ElectronCountsAst,
    Element,
    MoleculeAst,
    MulticenterBondAst,
    MulticenterBondConstraintAst,
    MulticenterBondConstraintKey,
    MulticenterBondConstraintsAst,
    MulticenterBondUpdate,
    ParseError,
    UnpairedElectronsAst,
    UnpairedElectronsUpdate,
    ValueAst,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            MulticenterBondUpdate(),
            (None, None, UnpairedElectronsUpdate(), MulticenterBondConstraintsAst([])),
        ),
        (
            MulticenterBondUpdate(
                electrons=[1, 1, 1],
                charge=1,
                unpaired_electrons=UnpairedElectronsUpdate(count=1),
                constraints=MulticenterBondConstraintsAst(
                    [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(3))]
                ),
            ),
            (
                ElectronCountsAst.Lit([1, 1, 1]),
                ValueAst.Lit(1),
                UnpairedElectronsUpdate(count=1),
                MulticenterBondConstraintsAst(
                    [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(3))]
                ),
            ),
        ),
        (
            MulticenterBondUpdate(
                unpaired_electrons=UnpairedElectronsUpdate(multiplicity=2),
                constraints=MulticenterBondConstraintsAst(
                    [
                        MulticenterBondConstraintAst.ElectronCount(
                            ValueAst.Undetermined()
                        )
                    ]
                ),
            ),
            (
                None,
                None,
                UnpairedElectronsUpdate(multiplicity=2),
                MulticenterBondConstraintsAst(
                    [
                        MulticenterBondConstraintAst.ElectronCount(
                            ValueAst.Undetermined()
                        )
                    ]
                ),
            ),
        ),
    ],
)
def test_multicenter_bond_update(update, expected):
    assert (
        update.electrons,
        update.charge,
        update.unpaired_electrons,
        update.constraints,
    ) == expected


@pytest.mark.parametrize(
    ("dsl", "canonical"),
    [
        ("", ""),
        ("[2,2,2]#c-1#s1", "[2,2,2]#c-#s"),
        ("*#c*#u*#s*#e*", "*#c*#u*#s*#e*"),
        ("#e*", "#e*"),
    ],
)
def test_multicenter_bond_update_parse(dsl, canonical):
    update = MulticenterBondUpdate.parse(dsl)
    assert str(update) == canonical
    assert repr(update) == f"MulticenterBondUpdate.parse('{canonical}')"
    assert eval(repr(update)) == update


def test_multicenter_bond_update_parse_error():
    with pytest.raises(ParseError):
        MulticenterBondUpdate.parse("#c+#c-")


def three_center_bond():
    # three borons (atom ids 0-2), one 3-center multicenter bond over all three
    return MoleculeAst.from_parts(
        [AtomAst(Element("B")) for _ in range(3)],
        multicenter_bonds=[([0, 1, 2], MulticenterBondAst([1, 1, 1]))],
    )


def test_multicenterbondast_new():
    bond = MulticenterBondAst([1, 1, 1])
    assert bond.electrons == ElectronCountsAst.Lit([1, 1, 1])
    assert bond.charge == ValueAst.Undetermined()
    assert len(bond.constraints) == 0


def test_multicenterbondast_new_kwargs():
    bond = MulticenterBondAst(
        [1, 1, 1],
        charge=-1,
        unpaired_electrons=UnpairedElectronsAst(0, 1),
    )
    assert bond.charge == ValueAst.Lit(-1)
    assert bond.unpaired_electrons == UnpairedElectronsAst(0, 1)


def test_multicenterbondast_constraints_kwarg():
    bond = MulticenterBondAst(
        [1, 1, 1],
        constraints=MulticenterBondConstraintsAst(
            [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))]
        ),
    )
    assert len(bond.constraints) == 1
    assert bond.constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondast_electrons_setter():
    bond = MulticenterBondAst([1, 1, 1])
    bond.electrons = [2, 2]
    assert bond.electrons == ElectronCountsAst.Lit([2, 2])


def test_multicenterbondast_charge_setter():
    bond = MulticenterBondAst([1, 1, 1])
    bond.charge = -1
    assert bond.charge == ValueAst.Lit(-1)


def test_multicenterbondast_unpaired_electrons_setter():
    bond = MulticenterBondAst([1, 1, 1])
    bond.unpaired_electrons = UnpairedElectronsAst(0, 1)
    assert bond.unpaired_electrons == UnpairedElectronsAst(0, 1)


@pytest.mark.parametrize("dsl", ["*", "[1,1,1]#e6", "[1,1,1]#c-2"])
def test_multicenterbondast_parse_roundtrip(dsl):
    bond = MulticenterBondAst.parse(dsl)
    assert str(bond) == dsl
    assert repr(bond) == f"MulticenterBondAst.parse('{dsl}')"


def test_multicenterbondast_parse_error():
    with pytest.raises(ParseError):
        MulticenterBondAst.parse("z")


def test_multicenterbondast_asdict():
    bond = MulticenterBondAst(
        [1, 1, 1],
        constraints=MulticenterBondConstraintsAst(
            [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))]
        ),
    )
    d = bond.asdict()
    assert set(d.keys()) == {
        "electrons",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["electrons"] == ElectronCountsAst.Lit([1, 1, 1])
    assert d["constraints"]["electron_count"] == ValueAst.Lit(6)


def test_multicenterbondast_set_constraints():
    bond = MulticenterBondAst([1, 1, 1])
    bond.constraints = MulticenterBondConstraintsAst(
        [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    assert bond.constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondconstraints_electron_count():
    constraints = MulticenterBondConstraintsAst([])
    assert constraints.electron_count == ValueAst.Undetermined()
    constraints.electron_count = 6
    assert constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondconstraints_mapping_ops():
    constraints = MulticenterBondConstraintsAst([])
    constraints.set(MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6)))
    assert len(constraints) == 1
    assert MulticenterBondConstraintKey.ElectronCount() in constraints
    assert constraints[MulticenterBondConstraintKey.ElectronCount()] == (
        MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))
    )
    assert [key for key in constraints] == [MulticenterBondConstraintKey.ElectronCount()]
    del constraints[MulticenterBondConstraintKey.ElectronCount()]
    assert len(constraints) == 0


def test_multicenterbondconstraints_getitem_missing():
    constraints = MulticenterBondConstraintsAst([])
    with pytest.raises(KeyError):
        constraints[MulticenterBondConstraintKey.ElectronCount()]


def test_multicenterbondconstraints_delitem_missing():
    constraints = MulticenterBondConstraintsAst([])
    with pytest.raises(KeyError):
        del constraints[MulticenterBondConstraintKey.ElectronCount()]


def test_multicenterbondconstraintkey_electron_count():
    key = MulticenterBondConstraintKey.ElectronCount()
    assert key == MulticenterBondConstraintKey.ElectronCount()
    assert key.__repr__().startswith("MulticenterBondConstraintKey.ElectronCount")


def test_multicenterbondconstraints_asdict():
    constraints = MulticenterBondConstraintsAst(
        [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"electron_count"}
    assert d["electron_count"] == ValueAst.Lit(6)


def test_multicenterbondview_fields():
    view = three_center_bond().multicenter_bonds[0]
    assert view.id == 0
    assert view.atom_ids == (0, 1, 2)
    assert view.electrons == ElectronCountsAst.Lit([1, 1, 1])
    assert repr(view) == "MulticenterBondView(id=0)"


def test_multicenterbondview_set_electrons():
    mol = three_center_bond()
    mol.multicenter_bonds[0].electrons = [2, 2, 2]
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.multicenter_bonds[0].electrons == ElectronCountsAst.Lit([2, 2, 2])


def test_multicenterbondview_set_charge():
    mol = three_center_bond()
    mol.multicenter_bonds[0].charge = -1
    assert mol.multicenter_bonds[0].charge == ValueAst.Lit(-1)


def test_multicenterbondview_set_unpaired_electrons():
    mol = three_center_bond()
    mol.multicenter_bonds[0].unpaired_electrons = UnpairedElectronsAst(0, 1)
    assert mol.multicenter_bonds[0].unpaired_electrons == UnpairedElectronsAst(0, 1)


def test_multicenterbondview_asdict():
    view = three_center_bond().multicenter_bonds[0]
    d = view.asdict()
    assert set(d.keys()) == {
        "electrons",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["electrons"] == ElectronCountsAst.Lit([1, 1, 1])


def test_multicenterbondview_constraints_write_through():
    mol = three_center_bond()
    mol.multicenter_bonds[0].constraints.set(
        MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))
    )
    constraints = mol.multicenter_bonds[0].constraints
    assert len(constraints) == 1
    assert constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondview_constraints_electron_count_property():
    mol = three_center_bond()
    mol.multicenter_bonds[0].constraints.electron_count = 6
    assert mol.multicenter_bonds[0].constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondview_set_constraints():
    mol = three_center_bond()
    mol.multicenter_bonds[0].constraints = MulticenterBondConstraintsAst(
        [MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    assert mol.multicenter_bonds[0].constraints.electron_count == ValueAst.Lit(6)


def test_multicenterbondviews_len_getitem():
    bonds = three_center_bond().multicenter_bonds
    assert len(bonds) == 1
    assert bonds[0].id == 0
    assert bonds[-1].id == 0
    with pytest.raises(IndexError):
        bonds[5]
    with pytest.raises(IndexError):
        bonds[-2]


def test_multicenterbondviews_setitem():
    mol = three_center_bond()
    mol.multicenter_bonds[0] = MulticenterBondAst([2, 2, 2])
    view = mol.multicenter_bonds[0]
    # value replaced, members preserved
    assert view.electrons == ElectronCountsAst.Lit([2, 2, 2])
    assert view.atom_ids == (0, 1, 2)


def test_multicenterbondviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        three_center_bond().multicenter_bonds[5] = MulticenterBondAst([1, 1, 1])


def test_multicenterbondviews_iter():
    ids = [view.id for view in three_center_bond().multicenter_bonds]
    assert ids == [0]


def test_multicenterbondviews_of():
    mol = three_center_bond()
    assert mol.multicenter_bonds.of([0, 1, 2]).id == 0
    # a subset is not the bond's exact atom set
    assert mol.multicenter_bonds.of([0, 1]) is None


def test_multicenterbondviews_incident():
    # three bonded borons plus one isolated boron (atom id 3)
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("B")) for _ in range(4)],
        multicenter_bonds=[([0, 1, 2], MulticenterBondAst([1, 1, 1]))],
    )
    assert [view.id for view in mol.multicenter_bonds.incident(0)] == [0]
    assert mol.multicenter_bonds.incident(3) == []


def test_multicenterbondviews_repr():
    assert repr(three_center_bond().multicenter_bonds) == "MulticenterBondViews(len=1)"
