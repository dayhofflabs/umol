import pytest

from umol import (
    AromaticSystemForm,
    AromaticSystemFieldChange,
    AromaticSystemUpdate,
    AtomForm,
    AtomConstraintForm,
    AtomConstraintsForm,
    AtomFieldChange,
    AtomUpdate,
    BondForm,
    BondFieldChange,
    BondUpdate,
    Constraint,
    ConstraintEdit,
    DativeBondForm,
    DativeBondFieldChange,
    DativeBondUpdate,
    Edit,
    Edits,
    ElectronCountsForm,
    Entity,
    MoleculeDefaults,
    MoleculeConstraint,
    MulticenterBondForm,
    MulticenterBondFieldChange,
    MulticenterBondUpdate,
    New,
    NoncovalentBondForm,
    NoncovalentBondFieldChange,
    NoncovalentBondKind,
    NoncovalentBondKindForm,
    NoncovalentBondUpdate,
    ParseError,
    StereoAtomForm,
    StereoAtomFieldChange,
    StereoAtomUpdate,
    StereoBondForm,
    StereoBondFieldChange,
    StereoBondUpdate,
    StereoConfigurationForm,
    StereoConfigurationUpdate,
    StereoCoset,
    StereoKind,
    StereoLigandKind,
    NumForm,
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
        0, AtomConstraintForm.Valence(NumForm.Lit(4))
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
        Edit.AddAtoms(atoms=[AtomForm.parse("C")]),
        Edit.AddBonds(bonds=[((0, New(0)), BondForm.parse("1"))]),
        Edit.RemoveTopology(atoms=[New(0)], bonds=[1]),
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(1)
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
    with pytest.raises(AttributeError):
        edit.extra = None


def test_edit_form_fields():
    atom = AtomForm.parse("C")
    bond = BondForm.parse("1")
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")

    edits = [
        Edit.AddAtoms(atoms=[atom]),
        Edit.AddBonds(bonds=[((0, New(0)), bond)]),
        Edit.AddDativeBond(atoms=[0, New(0)], attributes=dative),
        Edit.AddAromaticSystem(atoms=[0, New(0)], attributes=aromatic),
        Edit.AddMulticenterBond(atoms=[0, New(0)], attributes=multicenter),
        Edit.AddNoncovalentBond(
            atoms=(0, New(0)), attributes=noncovalent
        ),
        Edit.AddStereoAtom(
            site=New(0), ligands=[], attributes=stereo_atom
        ),
        Edit.AddStereoBond(
            site=New(0), ligands=[], attributes=stereo_bond
        ),
    ]

    atom.charge = -1
    bond.order = 2
    dative.order = 2
    aromatic.charge = -1
    multicenter.charge = -1
    noncovalent.kind = NoncovalentBondKind.Ionic
    stereo_atom.configuration = StereoConfigurationForm.Undetermined()
    stereo_bond.configuration = StereoConfigurationForm.Undetermined()

    retained = [
        edits[0].atoms[0],
        edits[1].bonds[0][1],
        edits[2].attributes,
        edits[3].attributes,
        edits[4].attributes,
        edits[5].attributes,
        edits[6].attributes,
        edits[7].attributes,
    ]
    assert retained[:6] == [
        AtomForm.parse("C"),
        BondForm.parse("1"),
        DativeBondForm(1),
        AromaticSystemForm([1, 1]),
        MulticenterBondForm([1, 1]),
        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
    ]
    assert retained[6].configuration == StereoAtomForm.parse(
        "Th0"
    ).configuration
    assert retained[7].configuration == StereoBondForm.parse(
        "Ct0"
    ).configuration
    assert [attributes.readonly for attributes in retained] == [True] * 8
    with pytest.raises(TypeError):
        retained[0].charge = 1
    returned_atoms = edits[0].atoms
    returned_atoms.clear()
    assert edits[0].atoms == [retained[0]]
    with pytest.raises(AttributeError):
        edits[0].atoms = []


def test_edit_removal_form_fields():
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")

    edits = [
        Edit.RemoveDativeBonds(removes=[(0, [0, 1], dative)]),
        Edit.RemoveAromaticSystems(removes=[(0, [0, 1], aromatic)]),
        Edit.RemoveMulticenterBonds(removes=[(0, [0, 1], multicenter)]),
        Edit.RemoveNoncovalentBonds(removes=[(0, (0, 1), noncovalent)]),
        Edit.RemoveStereoAtoms(removes=[(0, 0, [], stereo_atom)]),
        Edit.RemoveStereoBonds(removes=[(0, 0, [], stereo_bond)]),
    ]

    dative.order = 2
    aromatic.charge = -1
    multicenter.charge = -1
    noncovalent.kind = NoncovalentBondKind.Ionic
    stereo_atom.configuration = StereoConfigurationForm.Undetermined()
    stereo_bond.configuration = StereoConfigurationForm.Undetermined()

    retained = [edit.removes[0][-1] for edit in edits]
    assert retained[:4] == [
        DativeBondForm(1),
        AromaticSystemForm([1, 1]),
        MulticenterBondForm([1, 1]),
        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
    ]
    assert retained[4].configuration == StereoAtomForm.parse(
        "Th0"
    ).configuration
    assert retained[5].configuration == StereoBondForm.parse(
        "Ct0"
    ).configuration
    assert [attributes.readonly for attributes in retained] == [True] * 6


def test_edits():
    first = Edit.AddAtoms(atoms=[AtomForm.parse("C")])
    second = Edit.AddAtoms(atoms=[AtomForm.parse("N")])
    edits = Edits([first, first])
    iterator = iter(edits)

    edits.append(second)

    assert edits == Edits([first, first, second])
    assert len(edits) == 3
    assert edits[0] == first
    assert edits[-1] == second
    assert list(iterator) == [first, first]
    assert list(edits) == [first, first, second]


def test_edits_extend():
    first = Edit.AddAtoms(atoms=[AtomForm.parse("C")])
    second = Edit.RemoveTopology(atoms=[New(0)], bonds=[])
    source = Edits([first, second])
    target = Edits([second])

    assert target.extend(source) is None
    assert target.extend([first]) is None
    source.append(first)

    assert list(target) == [second, first, second, first]

    target.extend(target)

    assert list(target) == [second, first, second, first] * 2


@pytest.mark.parametrize("index", [-4, 3])
def test_edits_getitem_error(index):
    edits = Edits(
        [
            Edit.AddAtoms(atoms=[AtomForm.parse("C")]),
            Edit.AddAtoms(atoms=[AtomForm.parse("N")]),
            Edit.AddAtoms(atoms=[AtomForm.parse("O")]),
        ]
    )

    with pytest.raises(IndexError, match="edit index out of range"):
        edits[index]


def test_edits_add():
    edits = Edits()
    atom = AtomForm.parse("C")
    bond = BondForm.parse("1")
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")
    ligands = [(0, StereoLigandKind.Atom)]

    handles = (
        edits.add_atom(attributes=atom),
        edits.add_bond(0, New(0), attributes=bond),
        edits.add_dative_bond([0, New(0)], attributes=dative),
        edits.add_aromatic_system([0, New(0)], attributes=aromatic),
        edits.add_multicenter_bond([0, New(0)], attributes=multicenter),
        edits.add_noncovalent_bond((0, New(0)), attributes=noncovalent),
        edits.add_stereo_atom(New(0), ligands, attributes=stereo_atom),
        edits.add_stereo_bond(New(0), ligands, attributes=stereo_bond),
    )

    assert handles == (New(0),) * 8
    assert list(edits) == [
        Edit.AddAtoms(atoms=[atom]),
        Edit.AddBonds(bonds=[((0, New(0)), bond)]),
        Edit.AddDativeBond(atoms=[0, New(0)], attributes=dative),
        Edit.AddAromaticSystem(atoms=[0, New(0)], attributes=aromatic),
        Edit.AddMulticenterBond(atoms=[0, New(0)], attributes=multicenter),
        Edit.AddNoncovalentBond(atoms=(0, New(0)), attributes=noncovalent),
        Edit.AddStereoAtom(site=New(0), ligands=ligands, attributes=stereo_atom),
        Edit.AddStereoBond(site=New(0), ligands=ligands, attributes=stereo_bond),
    ]


def test_edits_add_many():
    edits = Edits()
    atoms = [AtomForm.parse("C"), AtomForm.parse("N")]
    bond = BondForm.parse("1")
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")
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
        Edit.AddDativeBond(atoms=[0, New(0)], attributes=dative),
        Edit.AddDativeBond(atoms=[New(0), New(1)], attributes=dative),
        Edit.AddAromaticSystem(atoms=[0, New(0)], attributes=aromatic),
        Edit.AddAromaticSystem(atoms=[New(0), New(1)], attributes=aromatic),
        Edit.AddMulticenterBond(atoms=[0, New(0)], attributes=multicenter),
        Edit.AddMulticenterBond(atoms=[New(0), New(1)], attributes=multicenter),
        Edit.AddNoncovalentBond(atoms=(0, New(0)), attributes=noncovalent),
        Edit.AddNoncovalentBond(
            atoms=(New(0), New(1)), attributes=noncovalent
        ),
        Edit.AddStereoAtom(site=New(0), ligands=ligands, attributes=stereo_atom),
        Edit.AddStereoAtom(site=New(1), ligands=ligands, attributes=stereo_atom),
        Edit.AddStereoBond(site=New(0), ligands=ligands, attributes=stereo_bond),
        Edit.AddStereoBond(site=New(1), ligands=ligands, attributes=stereo_bond),
    ]


def test_edits_constructor_counters():
    atom = AtomForm.parse("C")
    bond = BondForm.parse("1")
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")
    ligands = [(0, StereoLigandKind.Atom)]
    edits = Edits(
        [
            Edit.AddAtoms(atoms=[atom, atom]),
            Edit.AddBonds(
                bonds=[((0, 1), bond), ((1, 0), bond)]
            ),
            Edit.AddDativeBond(atoms=[0, 1], attributes=dative),
            Edit.AddAromaticSystem(atoms=[0, 1], attributes=aromatic),
            Edit.AddMulticenterBond(atoms=[0, 1], attributes=multicenter),
            Edit.AddNoncovalentBond(atoms=(0, 1), attributes=noncovalent),
            Edit.AddStereoAtom(site=0, ligands=ligands, attributes=stereo_atom),
            Edit.AddStereoBond(site=0, ligands=ligands, attributes=stereo_bond),
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
    dative = DativeBondForm(1)
    aromatic = AromaticSystemForm([1, 1])
    multicenter = MulticenterBondForm([1, 1])
    noncovalent = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    stereo_atom = StereoAtomForm.parse("Th0")
    stereo_bond = StereoBondForm.parse("Ct0")
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

    edits.update_atom(0, AtomForm.parse("C"), AtomUpdate())
    edits.update_bond(0, BondForm.parse("1"), BondUpdate())
    edits.update_dative_bond(0, DativeBondForm(1), DativeBondUpdate())
    edits.update_aromatic_system(
        0, AromaticSystemForm([1, 1]), AromaticSystemUpdate()
    )
    edits.update_multicenter_bond(
        0, MulticenterBondForm([1, 1]), MulticenterBondUpdate()
    )
    edits.update_noncovalent_bond(
        0,
        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
        NoncovalentBondUpdate(),
    )
    edits.update_stereo_atom(
        0, StereoAtomForm.parse("Th0"), StereoAtomUpdate()
    )
    edits.update_stereo_bond(
        0, StereoBondForm.parse("Ct0"), StereoBondUpdate()
    )

    assert list(edits) == []


def test_edits_update():
    edits = Edits()
    atom_constraints = AtomConstraintsForm(
        [AtomConstraintForm.Valence(NumForm.Lit(4))]
    )

    edits.update_atom(
        New(0),
        AtomForm.parse("C#c0#h3#v3"),
        AtomUpdate(
            charge=1,
            implicit_hydrogens=2,
            constraints=atom_constraints,
        ),
    )
    edits.update_bond(
        0,
        BondForm.parse("1#c0"),
        BondUpdate(order=2, charge=1),
    )
    edits.update_dative_bond(
        New(0), DativeBondForm(1), DativeBondUpdate(order=2)
    )
    edits.update_aromatic_system(
        0,
        AromaticSystemForm([1, 1], charge=0),
        AromaticSystemUpdate(electrons=[2, 0], charge=1),
    )
    edits.update_multicenter_bond(
        New(0),
        MulticenterBondForm([1, 1], charge=0),
        MulticenterBondUpdate(electrons=[2, 0], charge=1),
    )
    edits.update_noncovalent_bond(
        0,
        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
        NoncovalentBondUpdate(kind=NoncovalentBondKind.Ionic),
    )
    edits.update_stereo_atom(
        New(0),
        StereoAtomForm.parse("Th0"),
        StereoAtomUpdate(
            configuration=StereoConfigurationUpdate.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(1)
            )
        ),
    )
    edits.update_stereo_bond(
        0,
        StereoBondForm.parse("Ct0"),
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
                old=NumForm.Lit(0), new=NumForm.Lit(1)
            ),
        ),
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.ImplicitHydrogens(
                old=NumForm.Lit(3), new=NumForm.Lit(2)
            ),
        ),
        Edit.ModifyAtomConstraint(
            id=New(0),
            old=AtomConstraintForm.Valence(NumForm.Lit(3)),
            new=AtomConstraintForm.Valence(NumForm.Lit(4)),
        ),
        Edit.ModifyBondField(
            id=0,
            change=BondFieldChange.Order(
                old=NumForm.Lit(1), new=NumForm.Lit(2)
            ),
        ),
        Edit.ModifyBondField(
            id=0,
            change=BondFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(1)
            ),
        ),
        Edit.ModifyDativeBondField(
            id=New(0),
            change=DativeBondFieldChange.Order(
                old=NumForm.Lit(1), new=NumForm.Lit(2)
            ),
        ),
        Edit.ModifyAromaticSystemField(
            id=0,
            change=AromaticSystemFieldChange.Electrons(
                old=ElectronCountsForm.Lit([1, 1]),
                new=ElectronCountsForm.Lit([2, 0]),
            ),
        ),
        Edit.ModifyAromaticSystemField(
            id=0,
            change=AromaticSystemFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(1)
            ),
        ),
        Edit.ModifyMulticenterBondField(
            id=New(0),
            change=MulticenterBondFieldChange.Electrons(
                old=ElectronCountsForm.Lit([1, 1]),
                new=ElectronCountsForm.Lit([2, 0]),
            ),
        ),
        Edit.ModifyMulticenterBondField(
            id=New(0),
            change=MulticenterBondFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(1)
            ),
        ),
        Edit.ModifyNoncovalentBondField(
            id=0,
            change=NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindForm.Lit(
                    NoncovalentBondKind.HydrogenBond
                ),
                new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.Ionic),
            ),
        ),
        Edit.ModifyStereoAtomField(
            id=New(0),
            change=StereoAtomFieldChange.Configuration(
                old=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(0)
                ),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                ),
            ),
        ),
        Edit.ModifyStereoBondField(
            id=0,
            change=StereoBondFieldChange.Configuration(
                old=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                ),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(1)
                ),
            ),
        ),
    ]


def test_edits_parse_render():
    edits = Edits()
    atom = edits.add_atom(AtomForm.parse("C#h3"))
    bond = edits.add_bond(0, atom, BondForm(1))
    edits.add_dative_bond([0, atom], DativeBondForm(1))
    edits.add_aromatic_system([0, atom], AromaticSystemForm([1, 1]))
    edits.add_multicenter_bond([0, atom], MulticenterBondForm([1, 1]))
    edits.add_noncovalent_bond(
        (0, atom),
        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
    )
    edits.add_stereo_atom(
        atom,
        [(0, StereoLigandKind.Atom)],
        StereoAtomForm.parse("Th0"),
    )
    edits.add_stereo_bond(
        bond,
        [(atom, StereoLigandKind.Atom)],
        StereoBondForm.parse("Ct0"),
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
                atoms=[AtomForm.parse("O#i=#c0#h0#n0#u0#s")]
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
