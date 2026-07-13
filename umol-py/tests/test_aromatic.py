import pytest

from umol import (
    AromaticSystemAst,
    AromaticSystemConstraintAst,
    AromaticSystemConstraintKey,
    AromaticSystemConstraintsAst,
    AtomAst,
    ElectronCountsAst,
    Element,
    MoleculeAst,
    ParseError,
    SpinStateAst,
    ValueAst,
)


def benzene():
    # six aromatic carbons (atom ids 0-5), one aromatic system over all six
    return MoleculeAst.from_parts(
        [AtomAst(Element("C")) for _ in range(6)],
        aromatic_systems=[([0, 1, 2, 3, 4, 5], AromaticSystemAst([1, 1, 1, 1, 1, 1]))],
    )


def test_aromaticsystemast_new():
    system = AromaticSystemAst([1, 1, 1])
    assert system.electrons == ElectronCountsAst.Lit([1, 1, 1])
    assert system.charge == ValueAst.Undetermined()
    assert len(system.constraints) == 0


def test_aromaticsystemast_new_kwargs():
    system = AromaticSystemAst([1, 1, 1], charge=-1, spin=SpinStateAst(0, 1))
    assert system.charge == ValueAst.Lit(-1)
    assert system.spin == SpinStateAst(0, 1)


def test_aromaticsystemast_constraints_kwarg():
    system = AromaticSystemAst(
        [1, 1, 1],
        constraints=AromaticSystemConstraintsAst(
            [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
        ),
    )
    assert len(system.constraints) == 1
    assert system.constraints.electron_count == ValueAst.Lit(6)


def test_aromaticsystemast_electrons_setter():
    system = AromaticSystemAst([1, 1, 1])
    system.electrons = [2, 2]
    assert system.electrons == ElectronCountsAst.Lit([2, 2])


def test_aromaticsystemast_charge_setter():
    system = AromaticSystemAst([1, 1, 1])
    system.charge = -1
    assert system.charge == ValueAst.Lit(-1)


def test_aromaticsystemast_spin_setter():
    system = AromaticSystemAst([1, 1, 1])
    system.spin = SpinStateAst(0, 1)
    assert system.spin == SpinStateAst(0, 1)


@pytest.mark.parametrize("dsl", ["*", "[1,1,1]#e6", "[1,1,1]#c-2"])
def test_aromaticsystemast_parse_roundtrip(dsl):
    system = AromaticSystemAst.parse(dsl)
    assert str(system) == dsl
    assert repr(system) == f"AromaticSystemAst.parse('{dsl}')"


def test_aromaticsystemast_parse_error():
    with pytest.raises(ParseError):
        AromaticSystemAst.parse("z")


def test_aromaticsystemast_asdict():
    system = AromaticSystemAst(
        [1, 1, 1],
        constraints=AromaticSystemConstraintsAst(
            [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
        ),
    )
    d = system.asdict()
    assert set(d.keys()) == {"electrons", "charge", "spin", "constraints"}
    assert d["electrons"] == ElectronCountsAst.Lit([1, 1, 1])
    assert d["constraints"]["electron_count"] == ValueAst.Lit(6)


def test_aromaticsystemast_set_constraints():
    system = AromaticSystemAst([1, 1, 1])
    system.constraints = AromaticSystemConstraintsAst(
        [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    assert system.constraints.electron_count == ValueAst.Lit(6)


def test_aromaticsystemconstraints_electron_count():
    constraints = AromaticSystemConstraintsAst([])
    assert constraints.electron_count == ValueAst.Undetermined()
    constraints.electron_count = 6
    assert constraints.electron_count == ValueAst.Lit(6)


def test_aromaticsystemconstraints_mapping_ops():
    constraints = AromaticSystemConstraintsAst([])
    constraints.set(AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)))
    assert len(constraints) == 1
    assert AromaticSystemConstraintKey.ElectronCount() in constraints
    assert constraints[AromaticSystemConstraintKey.ElectronCount()] == (
        AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))
    )
    assert [key for key in constraints] == [AromaticSystemConstraintKey.ElectronCount()]
    del constraints[AromaticSystemConstraintKey.ElectronCount()]
    assert len(constraints) == 0


def test_aromaticsystemconstraints_getitem_missing():
    constraints = AromaticSystemConstraintsAst([])
    with pytest.raises(KeyError):
        constraints[AromaticSystemConstraintKey.ElectronCount()]


def test_aromaticsystemconstraints_delitem_missing():
    constraints = AromaticSystemConstraintsAst([])
    with pytest.raises(KeyError):
        del constraints[AromaticSystemConstraintKey.ElectronCount()]


def test_aromaticsystemconstraintkey_electron_count():
    key = AromaticSystemConstraintKey.ElectronCount()
    assert key == AromaticSystemConstraintKey.ElectronCount()
    assert key.__repr__().startswith("AromaticSystemConstraintKey.ElectronCount")


def test_aromaticsystemconstraints_asdict():
    constraints = AromaticSystemConstraintsAst(
        [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"electron_count"}
    assert d["electron_count"] == ValueAst.Lit(6)


def test_aromaticsystemview_fields():
    view = benzene().aromatic_systems[0]
    assert view.id == 0
    assert view.atom_ids == (0, 1, 2, 3, 4, 5)
    assert view.electrons == ElectronCountsAst.Lit([1, 1, 1, 1, 1, 1])
    assert repr(view) == "AromaticSystemView(id=0)"


def test_aromaticsystemview_set_electrons():
    mol = benzene()
    mol.aromatic_systems[0].electrons = [2, 2, 2, 2, 2, 2]
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.aromatic_systems[0].electrons == ElectronCountsAst.Lit([2, 2, 2, 2, 2, 2])


def test_aromaticsystemview_set_charge():
    mol = benzene()
    mol.aromatic_systems[0].charge = -1
    assert mol.aromatic_systems[0].charge == ValueAst.Lit(-1)


def test_aromaticsystemview_set_spin():
    mol = benzene()
    mol.aromatic_systems[0].spin = SpinStateAst(0, 1)
    assert mol.aromatic_systems[0].spin == SpinStateAst(0, 1)


def test_aromaticsystemview_asdict():
    view = benzene().aromatic_systems[0]
    d = view.asdict()
    assert set(d.keys()) == {"electrons", "charge", "spin", "constraints"}
    assert d["electrons"] == ElectronCountsAst.Lit([1, 1, 1, 1, 1, 1])


def test_aromaticsystemview_constraints_write_through():
    mol = benzene()
    mol.aromatic_systems[0].constraints.set(
        AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))
    )
    constraints = mol.aromatic_systems[0].constraints
    assert len(constraints) == 1
    assert constraints.electron_count == ValueAst.Lit(6)


def test_aromaticsystemview_constraints_electron_count_property():
    mol = benzene()
    mol.aromatic_systems[0].constraints.electron_count = 6
    assert mol.aromatic_systems[0].constraints.electron_count == ValueAst.Lit(6)


def test_aromaticsystemview_set_constraints():
    mol = benzene()
    mol.aromatic_systems[0].constraints = AromaticSystemConstraintsAst(
        [AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))]
    )
    assert mol.aromatic_systems[0].constraints.electron_count == ValueAst.Lit(6)


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
    mol.aromatic_systems[0] = AromaticSystemAst([2, 2, 2, 2, 2, 2])
    view = mol.aromatic_systems[0]
    # value replaced, members preserved
    assert view.electrons == ElectronCountsAst.Lit([2, 2, 2, 2, 2, 2])
    assert view.atom_ids == (0, 1, 2, 3, 4, 5)


def test_aromaticsystemviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        benzene().aromatic_systems[5] = AromaticSystemAst([1, 1, 1])


def test_aromaticsystemviews_iter():
    ids = [view.id for view in benzene().aromatic_systems]
    assert ids == [0]


def test_aromaticsystemviews_connecting():
    mol = benzene()
    assert mol.aromatic_systems.connecting([0, 1, 2, 3, 4, 5]).id == 0
    # a subset is not the system's exact atom set
    assert mol.aromatic_systems.connecting([0, 1, 2]) is None


def test_aromaticsystemviews_incident():
    # benzene's six carbons plus one isolated carbon (atom id 6)
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("C")) for _ in range(7)],
        aromatic_systems=[([0, 1, 2, 3, 4, 5], AromaticSystemAst([1, 1, 1, 1, 1, 1]))],
    )
    assert [view.id for view in mol.aromatic_systems.incident(0)] == [0]
    assert mol.aromatic_systems.incident(6) == []


def test_aromaticsystemviews_repr():
    assert repr(benzene().aromatic_systems) == "AromaticSystemViews(len=1)"
