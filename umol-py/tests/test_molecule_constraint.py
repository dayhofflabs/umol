from umol import (
    Constraint,
    Constraints,
    ConstraintsView,
    Molecule,
    MoleculeConstraint,
    RelationalConstraint,
    UnpairedElectronsForm,
)


def connected_constraint():
    return Constraint.Molecule(MoleculeConstraint.Connected(None))


def test_moleculeconstraint_unpaired_electron_coupling():
    constraint = MoleculeConstraint.UnpairedElectronCoupling(
        atoms=[2, 3],
        unpaired_electrons=UnpairedElectronsForm(1, 2),
    )

    assert constraint.atoms == [2, 3]
    assert constraint.unpaired_electrons == UnpairedElectronsForm(1, 2)
    assert repr(constraint) == (
        "MoleculeConstraint.UnpairedElectronCoupling(atoms=[2, 3], "
        "unpaired_electrons=UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2)))"
    )


def test_constraints_sequence():
    entry = connected_constraint()
    constraints = Constraints([entry, entry])

    assert len(constraints) == 2
    assert constraints[-1] == entry
    assert list(constraints) == [entry, entry]
    assert repr(constraints) == (
        "Constraints([Constraint.Molecule(MoleculeConstraint.Connected(None)), "
        "Constraint.Molecule(MoleculeConstraint.Connected(None))])"
    )


def test_molecule_from_entries_constraints():
    entry = connected_constraint()
    molecule = Molecule.from_entries([], constraints=[entry])

    assert isinstance(molecule.constraints, ConstraintsView)
    assert list(molecule.constraints) == [entry]


def test_molecule_constraints_live_view():
    molecule = Molecule()
    view = molecule.constraints
    entry = connected_constraint()

    view.append(entry)

    assert list(molecule.constraints) == [entry]
    molecule.constraints.clear()
    assert len(view) == 0


def test_molecule_constraints_set_container():
    molecule = Molecule()
    entry = connected_constraint()

    molecule.constraints = Constraints([entry, entry])

    assert list(molecule.constraints) == [entry, entry]


def test_molecule_constraints_set_view():
    source = Molecule.from_entries([], constraints=[connected_constraint()])
    target = Molecule()

    target.constraints = source.constraints

    assert list(target.constraints) == list(source.constraints)


def test_molecule_constraints_set_self():
    molecule = Molecule.from_entries([], constraints=[connected_constraint()])

    molecule.constraints = molecule.constraints

    assert list(molecule.constraints) == [connected_constraint()]
