import pytest

from umol import (
    AromaticSystemAst,
    AromaticSystemFieldChange,
    AromaticSystemUpdate,
    AtomAst,
    AtomConstraintAst,
    AtomConstraintsAst,
    AtomFieldChange,
    AtomUpdate,
    BondAst,
    BondFieldChange,
    BondUpdate,
    Constraint,
    ConstraintEdit,
    DativeBondAst,
    DativeBondFieldChange,
    DativeBondUpdate,
    Edit,
    Edits,
    ElectronCountsAst,
    Entity,
    MoleculeDefaults,
    MoleculeConstraint,
    MulticenterBondAst,
    MulticenterBondFieldChange,
    MulticenterBondUpdate,
    New,
    NoncovalentBondAst,
    NoncovalentBondFieldChange,
    NoncovalentBondKind,
    NoncovalentBondKindAst,
    NoncovalentBondUpdate,
    ParseError,
    StereoAtomAst,
    StereoAtomFieldChange,
    StereoAtomUpdate,
    StereoBondAst,
    StereoBondFieldChange,
    StereoBondUpdate,
    StereoConfigurationAst,
    StereoConfigurationUpdate,
    StereoCoset,
    StereoKind,
    StereoLigandKind,
    ValueAst,
)


def test_new():
    handle = New(3)

    assert handle.index == 3
    assert handle == New(3)
    assert handle != New(4)
    assert repr(handle) == "New(3)"


def test_new_immutability():
    handle = New(3)

    with pytest.raises(AttributeError):
        handle.index = 4
    with pytest.raises(AttributeError):
        handle.extra = 4


@pytest.mark.parametrize("index", [-1, 2**100])
def test_new_error(index):
    with pytest.raises(OverflowError):
        New(index)


def test_constraint_edit():
    constraint = Constraint.Atom(
        0, AtomConstraintAst.Valence(ValueAst.Lit(4))
    )

    identity = ConstraintEdit(constraint)
    created = ConstraintEdit(
        constraint,
        handles={Entity.Atom(0): New(0)},
    )

    assert identity == ConstraintEdit(constraint)
    assert created == ConstraintEdit(
        constraint,
        handles={Entity.Atom(0): New(0)},
    )
    assert created != identity
    assert repr(created) == "ConstraintEdit(...)"


@pytest.mark.parametrize(
    "edit",
    [
        Edit.AddAtoms(atoms=[AtomAst.parse("C")]),
        Edit.AddBonds(bonds=[((0, New(0)), BondAst.parse("1"))]),
        Edit.RemoveTopology(atoms=[New(0)], bonds=[1]),
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.AddMoleculeConstraint(
            constraint=ConstraintEdit(
                Constraint.Molecule(MoleculeConstraint.Connected(None))
            )
        ),
    ],
)
def test_edit(edit):
    assert edit == edit
    assert repr(edit).startswith(f"Edit.{type(edit).__name__}(")


def test_edits():
    first = Edit.AddAtoms(atoms=[AtomAst.parse("C")])
    second = Edit.AddAtoms(atoms=[AtomAst.parse("N")])
    edits = Edits([first, first])
    iterator = iter(edits)

    edits.append(second)

    assert edits == Edits([first, first, second])
    assert len(edits) == 3
    assert edits[0] == first
    assert edits[-1] == second
    assert list(iterator) == [first, first]
    assert list(edits) == [first, first, second]


@pytest.mark.parametrize("index", [-4, 3])
def test_edits_getitem_error(index):
    edits = Edits(
        [
            Edit.AddAtoms(atoms=[AtomAst.parse("C")]),
            Edit.AddAtoms(atoms=[AtomAst.parse("N")]),
            Edit.AddAtoms(atoms=[AtomAst.parse("O")]),
        ]
    )

    with pytest.raises(IndexError, match="edit index out of range"):
        edits[index]


def test_edits_add():
    edits = Edits()
    atom = AtomAst.parse("C")
    bond = BondAst.parse("1")
    dative = DativeBondAst(1)
    aromatic = AromaticSystemAst([1, 1])
    multicenter = MulticenterBondAst([1, 1])
    noncovalent = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomAst.parse("Th0")
    stereo_bond = StereoBondAst.parse("Ct0")
    ligands = [(0, StereoLigandKind.Atom)]

    handles = (
        edits.add_atom(atom),
        edits.add_bond(0, New(0), bond),
        edits.add_dative_bond([0, New(0)], dative),
        edits.add_aromatic_system([0, New(0)], aromatic),
        edits.add_multicenter_bond([0, New(0)], multicenter),
        edits.add_noncovalent_bond((0, New(0)), noncovalent),
        edits.add_stereo_atom(New(0), ligands, stereo_atom),
        edits.add_stereo_bond(New(0), ligands, stereo_bond),
    )

    assert handles == (New(0),) * 8
    assert list(edits) == [
        Edit.AddAtoms(atoms=[atom]),
        Edit.AddBonds(bonds=[((0, New(0)), bond)]),
        Edit.AddDativeBond(atoms=[0, New(0)], ast=dative),
        Edit.AddAromaticSystem(atoms=[0, New(0)], ast=aromatic),
        Edit.AddMulticenterBond(atoms=[0, New(0)], ast=multicenter),
        Edit.AddNoncovalentBond(atoms=(0, New(0)), ast=noncovalent),
        Edit.AddStereoAtom(site=New(0), ligands=ligands, ast=stereo_atom),
        Edit.AddStereoBond(site=New(0), ligands=ligands, ast=stereo_bond),
    ]


def test_edits_add_many():
    edits = Edits()
    atoms = [AtomAst.parse("C"), AtomAst.parse("N")]
    bond = BondAst.parse("1")
    dative = DativeBondAst(1)
    aromatic = AromaticSystemAst([1, 1])
    multicenter = MulticenterBondAst([1, 1])
    noncovalent = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomAst.parse("Th0")
    stereo_bond = StereoBondAst.parse("Ct0")
    ligands = [(0, StereoLigandKind.Atom)]
    atom_handles = edits.add_atoms(atoms)
    bond_handles = edits.add_bonds(
        [((0, New(0)), bond), ((New(0), New(1)), bond)]
    )
    dative_handles = edits.add_dative_bonds(
        [([0, New(0)], dative), ([New(0), New(1)], dative)]
    )
    aromatic_handles = edits.add_aromatic_systems(
        [([0, New(0)], aromatic), ([New(0), New(1)], aromatic)]
    )
    multicenter_handles = edits.add_multicenter_bonds(
        [([0, New(0)], multicenter), ([New(0), New(1)], multicenter)]
    )
    noncovalent_handles = edits.add_noncovalent_bonds(
        [((0, New(0)), noncovalent), ((New(0), New(1)), noncovalent)]
    )
    stereo_atom_handles = edits.add_stereo_atoms(
        [
            (New(0), ligands, stereo_atom),
            (New(1), ligands, stereo_atom),
        ]
    )
    stereo_bond_handles = edits.add_stereo_bonds(
        [
            (New(0), ligands, stereo_bond),
            (New(1), ligands, stereo_bond),
        ]
    )

    assert atom_handles == [New(0), New(1)]
    assert bond_handles == [New(0), New(1)]
    assert dative_handles == [New(0), New(1)]
    assert aromatic_handles == [New(0), New(1)]
    assert multicenter_handles == [New(0), New(1)]
    assert noncovalent_handles == [New(0), New(1)]
    assert stereo_atom_handles == [New(0), New(1)]
    assert stereo_bond_handles == [New(0), New(1)]
    assert list(edits) == [
        Edit.AddAtoms(atoms=atoms),
        Edit.AddBonds(
            bonds=[((0, New(0)), bond), ((New(0), New(1)), bond)]
        ),
        Edit.AddDativeBond(atoms=[0, New(0)], ast=dative),
        Edit.AddDativeBond(atoms=[New(0), New(1)], ast=dative),
        Edit.AddAromaticSystem(atoms=[0, New(0)], ast=aromatic),
        Edit.AddAromaticSystem(atoms=[New(0), New(1)], ast=aromatic),
        Edit.AddMulticenterBond(atoms=[0, New(0)], ast=multicenter),
        Edit.AddMulticenterBond(atoms=[New(0), New(1)], ast=multicenter),
        Edit.AddNoncovalentBond(atoms=(0, New(0)), ast=noncovalent),
        Edit.AddNoncovalentBond(
            atoms=(New(0), New(1)), ast=noncovalent
        ),
        Edit.AddStereoAtom(site=New(0), ligands=ligands, ast=stereo_atom),
        Edit.AddStereoAtom(site=New(1), ligands=ligands, ast=stereo_atom),
        Edit.AddStereoBond(site=New(0), ligands=ligands, ast=stereo_bond),
        Edit.AddStereoBond(site=New(1), ligands=ligands, ast=stereo_bond),
    ]


def test_edits_constructor_counters():
    atom = AtomAst.parse("C")
    bond = BondAst.parse("1")
    dative = DativeBondAst(1)
    aromatic = AromaticSystemAst([1, 1])
    multicenter = MulticenterBondAst([1, 1])
    noncovalent = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomAst.parse("Th0")
    stereo_bond = StereoBondAst.parse("Ct0")
    ligands = [(0, StereoLigandKind.Atom)]
    edits = Edits(
        [
            Edit.AddAtoms(atoms=[atom, atom]),
            Edit.AddBonds(
                bonds=[((0, 1), bond), ((1, 0), bond)]
            ),
            Edit.AddDativeBond(atoms=[0, 1], ast=dative),
            Edit.AddAromaticSystem(atoms=[0, 1], ast=aromatic),
            Edit.AddMulticenterBond(atoms=[0, 1], ast=multicenter),
            Edit.AddNoncovalentBond(atoms=(0, 1), ast=noncovalent),
            Edit.AddStereoAtom(site=0, ligands=ligands, ast=stereo_atom),
            Edit.AddStereoBond(site=0, ligands=ligands, ast=stereo_bond),
        ]
    )

    assert edits.add_atom(atom) == New(2)
    assert edits.add_bond(0, 1, bond) == New(2)
    assert edits.add_dative_bond([0, 1], dative) == New(1)
    assert edits.add_aromatic_system([0, 1], aromatic) == New(1)
    assert edits.add_multicenter_bond([0, 1], multicenter) == New(1)
    assert edits.add_noncovalent_bond((0, 1), noncovalent) == New(1)
    assert edits.add_stereo_atom(0, ligands, stereo_atom) == New(1)
    assert edits.add_stereo_bond(0, ligands, stereo_bond) == New(1)


def test_edits_remove():
    edits = Edits()
    dative = DativeBondAst(1)
    aromatic = AromaticSystemAst([1, 1])
    multicenter = MulticenterBondAst([1, 1])
    noncovalent = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomAst.parse("Th0")
    stereo_bond = StereoBondAst.parse("Ct0")
    ligands = [(New(0), StereoLigandKind.Atom)]

    edits.remove_topology([0, New(0)], [1, New(0)])
    edits.remove_dative_bonds([(New(0), [0, New(0)], dative)])
    edits.remove_aromatic_systems([(0, [0, New(0)], aromatic)])
    edits.remove_multicenter_bonds([(New(0), [0, New(0)], multicenter)])
    edits.remove_noncovalent_bonds([(0, (0, New(0)), noncovalent)])
    edits.remove_stereo_atoms(
        [(New(0), New(0), ligands, stereo_atom)]
    )
    edits.remove_stereo_bonds([(0, New(0), ligands, stereo_bond)])

    assert list(edits) == [
        Edit.RemoveTopology(atoms=[0, New(0)], bonds=[1, New(0)]),
        Edit.RemoveDativeBonds(
            removes=[(New(0), [0, New(0)], dative)]
        ),
        Edit.RemoveAromaticSystems(
            removes=[(0, [0, New(0)], aromatic)]
        ),
        Edit.RemoveMulticenterBonds(
            removes=[(New(0), [0, New(0)], multicenter)]
        ),
        Edit.RemoveNoncovalentBonds(
            removes=[(0, (0, New(0)), noncovalent)]
        ),
        Edit.RemoveStereoAtoms(
            removes=[(New(0), New(0), ligands, stereo_atom)]
        ),
        Edit.RemoveStereoBonds(
            removes=[(0, New(0), ligands, stereo_bond)]
        ),
    ]


def test_edits_molecule_constraint():
    constraint = ConstraintEdit(
        Constraint.Molecule(MoleculeConstraint.Connected(None))
    )
    edits = Edits()

    edits.add_molecule_constraint(constraint)
    edits.remove_molecule_constraint(constraint)

    assert list(edits) == [
        Edit.AddMoleculeConstraint(constraint=constraint),
        Edit.RemoveMoleculeConstraint(constraint=constraint),
    ]


def test_edits_update_empty():
    edits = Edits()

    edits.update_atom(0, AtomAst.parse("C"), AtomUpdate())
    edits.update_bond(0, BondAst.parse("1"), BondUpdate())
    edits.update_dative_bond(0, DativeBondAst(1), DativeBondUpdate())
    edits.update_aromatic_system(
        0, AromaticSystemAst([1, 1]), AromaticSystemUpdate()
    )
    edits.update_multicenter_bond(
        0, MulticenterBondAst([1, 1]), MulticenterBondUpdate()
    )
    edits.update_noncovalent_bond(
        0,
        NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
        NoncovalentBondUpdate(),
    )
    edits.update_stereo_atom(
        0, StereoAtomAst.parse("Th0"), StereoAtomUpdate()
    )
    edits.update_stereo_bond(
        0, StereoBondAst.parse("Ct0"), StereoBondUpdate()
    )

    assert list(edits) == []


def test_edits_update():
    edits = Edits()
    atom_constraints = AtomConstraintsAst(
        [AtomConstraintAst.Valence(ValueAst.Lit(4))]
    )

    edits.update_atom(
        New(0),
        AtomAst.parse("C#c0#h3#v3"),
        AtomUpdate(
            charge=1,
            implicit_hydrogens=2,
            constraints=atom_constraints,
        ),
    )
    edits.update_bond(
        0,
        BondAst.parse("1#c0"),
        BondUpdate(order=2, charge=1),
    )
    edits.update_dative_bond(
        New(0), DativeBondAst(1), DativeBondUpdate(order=2)
    )
    edits.update_aromatic_system(
        0,
        AromaticSystemAst([1, 1], charge=0),
        AromaticSystemUpdate(electrons=[2, 0], charge=1),
    )
    edits.update_multicenter_bond(
        New(0),
        MulticenterBondAst([1, 1], charge=0),
        MulticenterBondUpdate(electrons=[2, 0], charge=1),
    )
    edits.update_noncovalent_bond(
        0,
        NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
        NoncovalentBondUpdate(kind=NoncovalentBondKind.Ionic),
    )
    edits.update_stereo_atom(
        New(0),
        StereoAtomAst.parse("Th0"),
        StereoAtomUpdate(
            configuration=StereoConfigurationUpdate.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(1)
            )
        ),
    )
    edits.update_stereo_bond(
        0,
        StereoBondAst.parse("Ct0"),
        StereoBondUpdate(
            configuration=StereoConfigurationUpdate.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(1)
            )
        ),
    )

    assert list(edits) == [
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.ImplicitHydrogens(
                old=ValueAst.Lit(3), new=ValueAst.Lit(2)
            ),
        ),
        Edit.ModifyAtomConstraint(
            id=New(0),
            old=AtomConstraintAst.Valence(ValueAst.Lit(3)),
            new=AtomConstraintAst.Valence(ValueAst.Lit(4)),
        ),
        Edit.ModifyBondField(
            id=0,
            change=BondFieldChange.Order(
                old=ValueAst.Lit(1), new=ValueAst.Lit(2)
            ),
        ),
        Edit.ModifyBondField(
            id=0,
            change=BondFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.ModifyDativeBondField(
            id=New(0),
            change=DativeBondFieldChange.Order(
                old=ValueAst.Lit(1), new=ValueAst.Lit(2)
            ),
        ),
        Edit.ModifyAromaticSystemField(
            id=0,
            change=AromaticSystemFieldChange.Electrons(
                old=ElectronCountsAst.Lit([1, 1]),
                new=ElectronCountsAst.Lit([2, 0]),
            ),
        ),
        Edit.ModifyAromaticSystemField(
            id=0,
            change=AromaticSystemFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.ModifyMulticenterBondField(
            id=New(0),
            change=MulticenterBondFieldChange.Electrons(
                old=ElectronCountsAst.Lit([1, 1]),
                new=ElectronCountsAst.Lit([2, 0]),
            ),
        ),
        Edit.ModifyMulticenterBondField(
            id=New(0),
            change=MulticenterBondFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.ModifyNoncovalentBondField(
            id=0,
            change=NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindAst.Lit(
                    NoncovalentBondKind.HydrogenBond
                ),
                new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic),
            ),
        ),
        Edit.ModifyStereoAtomField(
            id=New(0),
            change=StereoAtomFieldChange.Configuration(
                old=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(0)
                ),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                ),
            ),
        ),
        Edit.ModifyStereoBondField(
            id=0,
            change=StereoBondFieldChange.Configuration(
                old=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                ),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(1)
                ),
            ),
        ),
    ]


def test_edits_parse_render():
    edits = Edits()
    atom = edits.add_atom(AtomAst.parse("C#h3"))
    bond = edits.add_bond(0, atom, BondAst(1))
    edits.add_dative_bond([0, atom], DativeBondAst(1))
    edits.add_aromatic_system([0, atom], AromaticSystemAst([1, 1]))
    edits.add_multicenter_bond([0, atom], MulticenterBondAst([1, 1]))
    edits.add_noncovalent_bond(
        (0, atom),
        NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
    )
    edits.add_stereo_atom(
        atom,
        [(0, StereoLigandKind.Atom)],
        StereoAtomAst.parse("Th0"),
    )
    edits.add_stereo_bond(
        bond,
        [(atom, StereoLigandKind.Atom)],
        StereoBondAst.parse("Ct0"),
    )
    edits.add_molecule_constraint(
        ConstraintEdit(
            Constraint.Molecule(MoleculeConstraint.Connected(None))
        )
    )

    rendered = edits.render()

    assert Edits.parse(rendered) == edits
    assert Edits.parse(rendered).render() == rendered


def test_edits_parse_render_defaults():
    source = '[{:atom {:add "O"}}]'
    defaults = MoleculeDefaults.ground()
    expected = Edits(
        [
            Edit.AddAtoms(
                atoms=[AtomAst.parse("O#i=#c0#h0#n0#u0#s")]
            )
        ]
    )

    parsed = Edits.parse(source, defaults=defaults)

    assert parsed == expected
    assert parsed.render(defaults=defaults) == source


def test_edits_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: trailing content at byte 4$",
    ):
        Edits.parse("not edn")
