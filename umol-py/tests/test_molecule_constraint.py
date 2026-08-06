from umol import (
    Constraint,
    Constraints,
    ConstraintsView,
    MoleculeAst,
    MoleculeConstraint,
    RelationalConstraint,
    SubPatternAnchor,
    UnpairedElectronsAst,
)


def connected_constraint():
    return Constraint.Molecule(MoleculeConstraint.Connected(None))


def test_moleculeconstraint_unpaired_electron_coupling():
    constraint = MoleculeConstraint.UnpairedElectronCoupling(
        atoms=[2, 3],
        unpaired_electrons=UnpairedElectronsAst(1, 2),
    )

    assert constraint.atoms == [2, 3]
    assert constraint.unpaired_electrons == UnpairedElectronsAst(1, 2)
    assert repr(constraint) == (
        "MoleculeConstraint.UnpairedElectronCoupling(atoms=[2, 3], "
        "unpaired_electrons=UnpairedElectronsAst(ValueAst.Lit(1), ValueAst.Lit(2)))"
    )


def test_subpatternanchor_fields():
    anchor = SubPatternAnchor(
        atoms=[(2, 0)],
        bonds=[(3, 1)],
        dative_bonds=[(4, 2)],
        aromatic_systems=[(5, 3)],
        multicenter_bonds=[(6, 4)],
        noncovalent_bonds=[(7, 5)],
        stereo_atoms=[(8, 6)],
        stereo_bonds=[(9, 7)],
    )

    assert anchor.atoms == [(2, 0)]
    assert anchor.bonds == [(3, 1)]
    assert anchor.dative_bonds == [(4, 2)]
    assert anchor.aromatic_systems == [(5, 3)]
    assert anchor.multicenter_bonds == [(6, 4)]
    assert anchor.noncovalent_bonds == [(7, 5)]
    assert anchor.stereo_atoms == [(8, 6)]
    assert anchor.stereo_bonds == [(9, 7)]


def test_constraint_pattern_match():
    constraint = Constraint.And(
        [
            Constraint.Relational(RelationalConstraint.DativeBondDonor(3, 5)),
            Constraint.Not(connected_constraint()),
        ]
    )

    match constraint:
        case Constraint.And(
            [
                Constraint.Relational(RelationalConstraint.DativeBondDonor(bond, atom)),
                Constraint.Not(Constraint.Molecule(MoleculeConstraint.Connected(atoms))),
            ]
        ):
            assert (bond, atom, atoms) == (3, 5, None)
        case _:
            raise AssertionError("constraint tree did not match its structural variants")


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


def test_moleculeast_from_entries_constraints():
    entry = connected_constraint()
    molecule = MoleculeAst.from_entries([], constraints=[entry])

    assert isinstance(molecule.constraints, ConstraintsView)
    assert list(molecule.constraints) == [entry]


def test_moleculeast_constraints_live_view():
    molecule = MoleculeAst()
    view = molecule.constraints
    entry = connected_constraint()

    view.append(entry)

    assert list(molecule.constraints) == [entry]
    molecule.constraints.clear()
    assert len(view) == 0


def test_moleculeast_constraints_set_container():
    molecule = MoleculeAst()
    entry = connected_constraint()

    molecule.constraints = Constraints([entry, entry])

    assert list(molecule.constraints) == [entry, entry]


def test_moleculeast_constraints_set_view():
    source = MoleculeAst.from_entries([], constraints=[connected_constraint()])
    target = MoleculeAst()

    target.constraints = source.constraints

    assert list(target.constraints) == list(source.constraints)


def test_moleculeast_constraints_set_self():
    molecule = MoleculeAst.from_entries([], constraints=[connected_constraint()])

    molecule.constraints = molecule.constraints

    assert list(molecule.constraints) == [connected_constraint()]
