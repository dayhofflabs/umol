import pytest

from umol import (
    AromaticSystemForm,
    AromaticSystemConstraintForm,
    AromaticSystemConstraintKey,
    AromaticSystemConstraintsForm,
    AromaticSystemUpdate,
    AtomForm,
    ElectronCountsForm,
    Element,
    Molecule,
    ParseError,
    UnpairedElectronsForm,
    UnpairedElectronsUpdate,
    NumForm,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            AromaticSystemUpdate(),
            (None, None, UnpairedElectronsUpdate(), AromaticSystemConstraintsForm([])),
        ),
        (
            AromaticSystemUpdate(
                electrons=[1, 1, 1],
                charge=-1,
                unpaired_electrons=UnpairedElectronsUpdate(count=2),
                constraints=AromaticSystemConstraintsForm(
                    [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
                ),
            ),
            (
                ElectronCountsForm.Lit([1, 1, 1]),
                NumForm.Lit(-1),
                UnpairedElectronsUpdate(count=2),
                AromaticSystemConstraintsForm(
                    [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
                ),
            ),
        ),
        (
            AromaticSystemUpdate(
                unpaired_electrons=UnpairedElectronsUpdate(multiplicity=1),
                constraints=AromaticSystemConstraintsForm(
                    [
                        AromaticSystemConstraintForm.ElectronCount(
                            NumForm.Undetermined()
                        )
                    ]
                ),
            ),
            (
                None,
                None,
                UnpairedElectronsUpdate(multiplicity=1),
                AromaticSystemConstraintsForm(
                    [
                        AromaticSystemConstraintForm.ElectronCount(
                            NumForm.Undetermined()
                        )
                    ]
                ),
            ),
        ),
    ],
)
def test_aromatic_system_update(update, expected):
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
def test_aromatic_system_update_parse(dsl, canonical):
    update = AromaticSystemUpdate.parse(dsl)
    assert str(update) == canonical
    assert repr(update) == f"AromaticSystemUpdate.parse('{canonical}')"
    assert eval(repr(update)) == update


def test_aromatic_system_update_parse_error():
    with pytest.raises(ParseError):
        AromaticSystemUpdate.parse("#c+#c-")


def benzene():
    # six aromatic carbons (atom ids 0-5), one aromatic system over all six
    return Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(6)],
        aromatic_systems=[([0, 1, 2, 3, 4, 5], AromaticSystemForm([1, 1, 1, 1, 1, 1]))],
    )


def test_aromaticsystemast_new():
    system = AromaticSystemForm([1, 1, 1])
    assert system.electrons == ElectronCountsForm.Lit([1, 1, 1])
    assert system.charge == NumForm.Undetermined()
    assert len(system.constraints) == 0


def test_aromaticsystemast_new_kwargs():
    system = AromaticSystemForm(
        [1, 1, 1],
        charge=-1,
        unpaired_electrons=UnpairedElectronsForm(0, 1),
    )
    assert system.charge == NumForm.Lit(-1)
    assert system.unpaired_electrons == UnpairedElectronsForm(0, 1)


def test_aromaticsystemast_constraints_kwarg():
    system = AromaticSystemForm(
        [1, 1, 1],
        constraints=AromaticSystemConstraintsForm(
            [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
        ),
    )
    assert len(system.constraints) == 1
    assert system.constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemast_electrons_setter():
    system = AromaticSystemForm([1, 1, 1])
    system.electrons = [2, 2]
    assert system.electrons == ElectronCountsForm.Lit([2, 2])


def test_aromaticsystemast_charge_setter():
    system = AromaticSystemForm([1, 1, 1])
    system.charge = -1
    assert system.charge == NumForm.Lit(-1)


def test_aromaticsystemast_unpaired_electrons_setter():
    system = AromaticSystemForm([1, 1, 1])
    system.unpaired_electrons = UnpairedElectronsForm(0, 1)
    assert system.unpaired_electrons == UnpairedElectronsForm(0, 1)


@pytest.mark.parametrize("dsl", ["*", "[1,1,1]#e6", "[1,1,1]#c-2"])
def test_aromaticsystemast_parse_roundtrip(dsl):
    system = AromaticSystemForm.parse(dsl)
    assert str(system) == dsl
    assert repr(system) == f"AromaticSystemForm.parse('{dsl}')"


def test_aromaticsystemast_parse_error():
    with pytest.raises(ParseError):
        AromaticSystemForm.parse("z")


def test_aromaticsystemast_asdict():
    system = AromaticSystemForm(
        [1, 1, 1],
        constraints=AromaticSystemConstraintsForm(
            [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
        ),
    )
    d = system.asdict()
    assert set(d.keys()) == {
        "electrons",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["electrons"] == ElectronCountsForm.Lit([1, 1, 1])
    assert d["constraints"]["electron_count"] == NumForm.Lit(6)


def test_aromaticsystemast_set_constraints():
    system = AromaticSystemForm([1, 1, 1])
    system.constraints = AromaticSystemConstraintsForm(
        [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
    )
    assert system.constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemconstraints_electron_count():
    constraints = AromaticSystemConstraintsForm([])
    assert constraints.electron_count == NumForm.Undetermined()
    constraints.electron_count = 6
    assert constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemconstraints_mapping_ops():
    constraints = AromaticSystemConstraintsForm([])
    constraints.set(AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6)))
    assert len(constraints) == 1
    assert AromaticSystemConstraintKey.ElectronCount() in constraints
    assert constraints[AromaticSystemConstraintKey.ElectronCount()] == (
        AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))
    )
    assert [key for key in constraints] == [AromaticSystemConstraintKey.ElectronCount()]
    del constraints[AromaticSystemConstraintKey.ElectronCount()]
    assert len(constraints) == 0


def test_aromaticsystemconstraints_getitem_missing():
    constraints = AromaticSystemConstraintsForm([])
    with pytest.raises(KeyError):
        constraints[AromaticSystemConstraintKey.ElectronCount()]


def test_aromaticsystemconstraints_delitem_missing():
    constraints = AromaticSystemConstraintsForm([])
    with pytest.raises(KeyError):
        del constraints[AromaticSystemConstraintKey.ElectronCount()]


def test_aromaticsystemconstraintkey_electron_count():
    key = AromaticSystemConstraintKey.ElectronCount()
    assert key == AromaticSystemConstraintKey.ElectronCount()
    assert key.__repr__().startswith("AromaticSystemConstraintKey.ElectronCount")


def test_aromaticsystemconstraints_asdict():
    constraints = AromaticSystemConstraintsForm(
        [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"electron_count"}
    assert d["electron_count"] == NumForm.Lit(6)


def test_aromaticsystemview_fields():
    view = benzene().aromatic_systems[0]
    assert view.id == 0
    assert view.atom_ids == (0, 1, 2, 3, 4, 5)
    assert view.electrons == ElectronCountsForm.Lit([1, 1, 1, 1, 1, 1])
    assert repr(view) == "AromaticSystemView(id=0)"


def test_aromaticsystemview_set_electrons():
    mol = benzene()
    mol.aromatic_systems[0].electrons = [2, 2, 2, 2, 2, 2]
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.aromatic_systems[0].electrons == ElectronCountsForm.Lit([2, 2, 2, 2, 2, 2])


def test_aromaticsystemview_set_charge():
    mol = benzene()
    mol.aromatic_systems[0].charge = -1
    assert mol.aromatic_systems[0].charge == NumForm.Lit(-1)


def test_aromaticsystemview_set_unpaired_electrons():
    mol = benzene()
    mol.aromatic_systems[0].unpaired_electrons = UnpairedElectronsForm(0, 1)
    assert mol.aromatic_systems[0].unpaired_electrons == UnpairedElectronsForm(0, 1)


def test_aromaticsystemview_asdict():
    view = benzene().aromatic_systems[0]
    d = view.asdict()
    assert set(d.keys()) == {
        "electrons",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["electrons"] == ElectronCountsForm.Lit([1, 1, 1, 1, 1, 1])


def test_aromaticsystemview_constraints_write_through():
    mol = benzene()
    mol.aromatic_systems[0].constraints.set(
        AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))
    )
    constraints = mol.aromatic_systems[0].constraints
    assert len(constraints) == 1
    assert constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemview_constraints_electron_count_property():
    mol = benzene()
    mol.aromatic_systems[0].constraints.electron_count = 6
    assert mol.aromatic_systems[0].constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemview_set_constraints():
    mol = benzene()
    mol.aromatic_systems[0].constraints = AromaticSystemConstraintsForm(
        [AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))]
    )
    assert mol.aromatic_systems[0].constraints.electron_count == NumForm.Lit(6)


def test_aromaticsystemviews_len_getitem():
    systems = benzene().aromatic_systems
    assert len(systems) == 1
    assert systems[0].id == 0
    assert systems[-1].id == 0
    with pytest.raises(IndexError):
        systems[5]
    with pytest.raises(IndexError):
        systems[-2]


def test_aromaticsystemviews_setitem():
    mol = benzene()
    mol.aromatic_systems[0] = AromaticSystemForm([2, 2, 2, 2, 2, 2])
    view = mol.aromatic_systems[0]
    # value replaced, members preserved
    assert view.electrons == ElectronCountsForm.Lit([2, 2, 2, 2, 2, 2])
    assert view.atom_ids == (0, 1, 2, 3, 4, 5)


def test_aromaticsystemviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        benzene().aromatic_systems[5] = AromaticSystemForm([1, 1, 1])


def test_aromaticsystemviews_iter():
    ids = [view.id for view in benzene().aromatic_systems]
    assert ids == [0]


def test_aromaticsystemviews_of():
    mol = benzene()
    assert mol.aromatic_systems.of([0, 1, 2, 3, 4, 5]).id == 0
    # a subset is not the system's exact atom set
    assert mol.aromatic_systems.of([0, 1, 2]) is None


def test_aromaticsystemviews_incident():
    # benzene's six carbons plus one isolated carbon (atom id 6)
    mol = Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(7)],
        aromatic_systems=[([0, 1, 2, 3, 4, 5], AromaticSystemForm([1, 1, 1, 1, 1, 1]))],
    )
    assert [view.id for view in mol.aromatic_systems.incident(0)] == [0]
    assert mol.aromatic_systems.incident(6) == []


def test_aromaticsystemviews_repr():
    assert repr(benzene().aromatic_systems) == "AromaticSystemViews(len=1)"
