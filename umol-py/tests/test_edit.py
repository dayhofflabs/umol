import pytest

from umol import (
    AromaticSystemForm,
    AromaticSystemConstraintForm,
    AromaticSystemConstraintKey,
    AromaticSystemFieldChange,
    AromaticSystemUpdate,
    AtomForm,
    AtomConstraintForm,
    AtomConstraintKey,
    AtomConstraintsForm,
    AtomFieldChange,
    AtomUpdate,
    BooleanForm,
    BondConstraintForm,
    BondConstraintKey,
    BondForm,
    BondFieldChange,
    BondUpdate,
    Constraint,
    ConstraintEdit,
    ConsumedError,
    DativeBondForm,
    DativeBondConstraintForm,
    DativeBondConstraintKey,
    DativeBondFieldChange,
    DativeBondUpdate,
    Edit,
    Edits,
    ElectronCountsForm,
    Element,
    ElementForm,
    IsotopeMassForm,
    Entity,
    InvalidatedViewError,
    Molecule,
    MoleculeDefaults,
    MoleculeConstraint,
    MulticenterBondForm,
    MulticenterBondConstraintForm,
    MulticenterBondConstraintKey,
    MulticenterBondFieldChange,
    MulticenterBondUpdate,
    New,
    NoncovalentBondForm,
    NoncovalentBondConstraintForm,
    NoncovalentBondConstraintKey,
    NoncovalentBondFieldChange,
    NoncovalentBondKind,
    NoncovalentBondKindForm,
    NoncovalentBondUpdate,
    ParseError,
    RingScope,
    RingMembershipForm,
    UnpairedElectronsForm,
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
    TransactionError,
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


@pytest.mark.parametrize(
    ("variant", "old", "new"),
    [
        (Edit.ReplaceDativeBondDonors, [0], [New(0), 2]),
        (Edit.ReplaceDativeBondAcceptor, 0, New(0)),
        (Edit.ReplaceAromaticSystemAtoms, [0, 1], [New(0), 1, 2]),
        (Edit.ReplaceMulticenterBondAtoms, [0, 1], [New(0), 1, 2]),
        (Edit.ReplaceNoncovalentBondAtoms, (0, 1), (New(0), 2)),
        (Edit.ReplaceStereoAtomSite, 0, New(0)),
        (
            Edit.ReplaceStereoAtomLigands,
            [(0, StereoLigandKind.Atom)],
            [(New(0), StereoLigandKind.ImplicitHydrogen),
             (1, StereoLigandKind.LonePair)],
        ),
        (Edit.ReplaceStereoBondSite, 0, New(0)),
        (
            Edit.ReplaceStereoBondLigands,
            [(0, StereoLigandKind.Atom)],
            [(New(0), StereoLigandKind.ImplicitHydrogen),
             (1, StereoLigandKind.LonePair)],
        ),
    ],
)
def test_edit_replace(variant, old, new):
    edit = variant(id=New(0), old=old, new=new)
    edits = Edits([edit])

    assert (edit.id, edit.old, edit.new) == (New(0), old, new)
    assert repr(edit) == (
        f"Edit.{variant.__name__}(id=New(0), old={old!r}, new={new!r})"
    )
    assert list(edits) == [edit]
    assert Edits.parse(edits.render()) == edits
    for field in ("id", "old", "new"):
        with pytest.raises(AttributeError):
            setattr(edit, field, None)


@pytest.mark.parametrize("atoms", [(), (0,), (0, 1, 2)])
def test_edit_replace_noncovalent_bond_atoms_error(atoms):
    with pytest.raises(ValueError):
        Edit.ReplaceNoncovalentBondAtoms(id=0, old=(0, 1), new=atoms)


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
        Edit.AddDativeBond(donors=[0], acceptor=New(0), attributes=dative),
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
        Edit.RemoveDativeBonds(removes=[(0, [0], 1, dative)]),
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


@pytest.mark.parametrize("count", [0, 1, 2])
def test_edits_iter(count):
    expected = [
        Edit.AddAtoms(atoms=[AtomForm.parse("C")]),
        Edit.RemoveTopology(atoms=[0], bonds=[]),
    ][:count]
    edits = Edits(expected)
    first = iter(edits)
    second = iter(edits)
    assert iter(first) is first

    edits.append(Edit.AddAtoms(atoms=[AtomForm.parse("N")]))
    second_values = list(second)
    first_values = list(first)

    assert first_values == expected
    assert second_values == expected
    for left, right, original in zip(first_values, second_values, expected):
        assert type(left) is type(original)
        assert type(right) is type(original)
        assert left is not right
    if count:
        assert first_values[0].atoms[0] is not second_values[0].atoms[0]
        with pytest.raises(TypeError):
            first_values[0].atoms[0].charge = -1
        assert edits[0] == expected[0]

    edits.append(Edit.AddAtoms(atoms=[AtomForm.parse("O")]))
    assert list(first) == []
    assert list(second) == []
    assert list(first) == []


def test_edits_iter_owner():
    expected = [
        Edit.AddAtoms(atoms=[AtomForm.parse("C")]),
        Edit.AddAtoms(atoms=[AtomForm.parse("N")]),
    ]
    edits = Edits(expected)
    first = iter(edits)
    second = iter(edits)
    assert next(first) == expected[0]
    del edits

    assert list(second) == expected
    assert list(first) == expected[1:]


@pytest.mark.parametrize("receiver", ["molecule", "editor"])
@pytest.mark.parametrize("method", ["apply", "tracked_apply"])
def test_edits_consumption(receiver, method):
    atom = AtomForm.parse("C")
    edit = Edit.AddAtoms(atoms=[atom])
    edits = Edits([edit])
    alias = edits
    copied = Edits(list(edits))
    indexed = edits[0]
    yielded = next(iter(edits))
    molecule = Molecule()
    target = molecule if receiver == "molecule" else molecule.edit()

    result = getattr(target, method)(edits)
    if receiver == "editor":
        result = target.finish()
    elif method == "tracked_apply":
        result = result[0]

    assert result == Molecule.parse('{:atoms ["C"]}')
    assert molecule == Molecule()
    with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
        len(alias)
    assert copied == Edits([edit])
    assert indexed == yielded == edit
    assert indexed.atoms == yielded.atoms == [atom]
    assert repr(copied) == f"Edits([{edit!r}])"
    with pytest.raises(TypeError):
        hash(copied)


@pytest.mark.parametrize("receiver", ["molecule", "editor"])
@pytest.mark.parametrize("method", ["apply", "tracked_apply"])
def test_edits_consumption_error(receiver, method):
    edits = Edits([
        Edit.AddAtoms(atoms=[AtomForm.parse("N")]),
        Edit.ModifyAtomField(
            id=0,
            change=AtomFieldChange.Charge(old=NumForm.Lit(1), new=NumForm.Lit(2)),
        ),
    ])
    molecule = Molecule.parse('{:atoms ["C#c0"]}')
    target = molecule if receiver == "molecule" else molecule.edit()

    with pytest.raises(TransactionError, match="^precondition failed: old state does not match current$"):
        getattr(target, method)(edits)
    with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
        edits.render()
    assert molecule == Molecule.parse('{:atoms ["C#c0"]}')


@pytest.mark.parametrize(
    ("method", "args"),
    [
        ("__len__", ()),
        ("__repr__", ()),
        ("__getitem__", (0,)),
        ("__getitem__", (-1,)),
        ("__iter__", ()),
        ("render", ()),
        ("append", (Edit.AddAtoms(atoms=[]),)),
        ("add_atom", (AtomForm.parse("C"),)),
        ("add_atoms", ([],)),
        ("add_bond", (0, 1, BondForm(1))),
        ("add_bonds", ([],)),
        ("add_dative_bond", ([0], 1, DativeBondForm(1))),
        ("add_dative_bonds", ([],)),
        ("add_aromatic_system", ([0, 1], AromaticSystemForm([1, 1]))),
        ("add_aromatic_systems", ([],)),
        ("add_multicenter_bond", ([0, 1], MulticenterBondForm([1, 1]))),
        ("add_multicenter_bonds", ([],)),
        ("add_noncovalent_bond", ((0, 1), NoncovalentBondForm(NoncovalentBondKind.Ionic))),
        ("add_noncovalent_bonds", ([],)),
        ("add_stereo_atom", (0, [], StereoAtomForm.parse("Th0"))),
        ("add_stereo_atoms", ([],)),
        ("add_stereo_bond", (0, [], StereoBondForm.parse("Ct0"))),
        ("add_stereo_bonds", ([],)),
        ("remove_topology", ([], [])),
        ("remove_dative_bonds", ([],)),
        ("remove_aromatic_systems", ([],)),
        ("remove_multicenter_bonds", ([],)),
        ("remove_noncovalent_bonds", ([],)),
        ("remove_stereo_atoms", ([],)),
        ("remove_stereo_bonds", ([],)),
        ("add_molecule_constraint", (ConstraintEdit(Constraint.Molecule(MoleculeConstraint.Connected(None))),)),
        ("remove_molecule_constraint", (ConstraintEdit(Constraint.Molecule(MoleculeConstraint.Connected(None))),)),
        ("update_atom", (0, AtomForm.parse("C"), AtomUpdate())),
        ("update_bond", (0, BondForm(1), BondUpdate())),
        ("update_dative_bond", (0, DativeBondForm(1), DativeBondUpdate())),
        ("update_aromatic_system", (0, AromaticSystemForm([1, 1]), AromaticSystemUpdate())),
        ("update_multicenter_bond", (0, MulticenterBondForm([1, 1]), MulticenterBondUpdate())),
        ("update_noncovalent_bond", (0, NoncovalentBondForm(NoncovalentBondKind.Ionic), NoncovalentBondUpdate())),
        ("update_stereo_atom", (0, StereoAtomForm.parse("Th0"), StereoAtomUpdate())),
        ("update_stereo_bond", (0, StereoBondForm.parse("Ct0"), StereoBondUpdate())),
    ],
)
def test_edits_consumed(method, args):
    edits = Edits()
    Molecule().apply(edits)

    with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
        getattr(edits, method)(*args)


@pytest.mark.parametrize("lhs_consumed", [False, True])
@pytest.mark.parametrize("rhs_consumed", [False, True])
def test_edits_equality(lhs_consumed, rhs_consumed):
    lhs = Edits()
    rhs = Edits()
    if lhs_consumed:
        Molecule().apply(lhs)
    if rhs_consumed:
        Molecule().apply(rhs)

    if lhs_consumed or rhs_consumed:
        with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
            lhs == rhs
        with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
            lhs != rhs
    else:
        assert lhs == rhs
        assert not lhs != rhs


def test_edits_equality_alias():
    edits = Edits()
    alias = edits
    Molecule().apply(edits)

    with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
        edits == alias


@pytest.mark.parametrize(
    ("count", "consumed"),
    [(0, 0), (2, 0), (2, 1), (2, 2)],
    ids=["empty", "live", "partial", "exhausted"],
)
def test_edit_iter_consumed(count, consumed):
    edits = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("C")])] * count)
    iterator = iter(edits)
    for _ in range(consumed):
        next(iterator)
    if count == consumed:
        with pytest.raises(StopIteration):
            next(iterator)
    Molecule().apply(edits)

    with pytest.raises(
        InvalidatedViewError,
        match="^Edit iterator is invalid because Edits has been consumed$",
    ):
        next(iterator)


def test_edits_append():
    edits = Edits()
    carbon = edits.add_atom(AtomForm.parse("C"))
    nitrogen = edits.add_atom(AtomForm.parse("N"))
    edits.append(Edit.RemoveTopology(atoms=[nitrogen], bonds=[]))

    assert carbon == New(0)
    assert nitrogen == New(1)
    assert Molecule().apply(edits) == Molecule.parse('{:atoms ["C"]}')


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
        edits.add_dative_bond([0], New(0), attributes=dative),
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
        Edit.AddDativeBond(donors=[0], acceptor=New(0), attributes=dative),
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
        [([0], New(0), dative), ([New(0)], New(1), dative)]
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
        Edit.AddDativeBond(donors=[0], acceptor=New(0), attributes=dative),
        Edit.AddDativeBond(donors=[New(0)], acceptor=New(1), attributes=dative),
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
            Edit.AddDativeBond(donors=[0], acceptor=1, attributes=dative),
            Edit.AddAromaticSystem(atoms=[0, 1], attributes=aromatic),
            Edit.AddMulticenterBond(atoms=[0, 1], attributes=multicenter),
            Edit.AddNoncovalentBond(atoms=(0, 1), attributes=noncovalent),
            Edit.AddStereoAtom(site=0, ligands=ligands, attributes=stereo_atom),
            Edit.AddStereoBond(site=0, ligands=ligands, attributes=stereo_bond),
        ]
    )

    assert edits.add_atom(atom) == New(2)
    assert edits.add_bond(0, 1, bond) == New(2)
    assert edits.add_dative_bond([0], 1, dative) == New(1)
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
    edits.remove_dative_bonds([(New(0), [0], New(0), dative)])
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
            removes=[(New(0), [0], New(0), dative)]
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
    edits.add_dative_bond([0], atom, DativeBondForm(1))
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
    defaults = MoleculeDefaults.concrete()
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


@pytest.mark.parametrize(
    "variant,method,field,old,new",
    [
        (Edit.ModifyAtomElements, "modify_atom_elements", "element",
         ElementForm.Lit(Element("C")), ElementForm.Lit(Element("N"))),
        (Edit.ModifyAtomIsotopeMasses, "modify_atom_isotope_masses", "isotope_mass",
         IsotopeMassForm.Lit(12), IsotopeMassForm.Lit(13)),
        (Edit.ModifyAtomCharges, "modify_atom_charges", "charge",
         NumForm.Lit(0), NumForm.LitSet({-1, 1})),
        (Edit.ModifyAtomImplicitHydrogens, "modify_atom_implicit_hydrogens", "implicit_hydrogens",
         NumForm.Lit(1), NumForm.Lit(2)),
        (Edit.ModifyAtomLonePairs, "modify_atom_lone_pairs", "lone_pairs",
         NumForm.Lit(0), NumForm.Lit(1)),
        (Edit.ModifyAtomUnpairedElectrons, "modify_atom_unpaired_electrons", "unpaired_electrons",
         UnpairedElectronsForm(0, 1), UnpairedElectronsForm(1, 2)),
    ],
)
@pytest.mark.parametrize("count", [0, 1, 2])
def test_edits_atom_batches(variant, method, field, old, new, count):
    rows = [(0, old, new), (0, new, old)][:count]
    edit = variant(changes=rows)
    expected_rows = tuple(rows)
    rows.clear()
    assert edit.changes == expected_rows
    assert isinstance(edit.changes, tuple)
    assert all(isinstance(row, tuple) for row in edit.changes)
    with pytest.raises(AttributeError):
        edit.changes = ()
    edits = Edits()
    getattr(edits, method)(list(expected_rows))
    assert list(edits) == [edit]
    assert list(Edits.parse(edits.render())) == [edit]
    assert repr(edit).startswith("Edit." + variant.__name__ + "(changes=")
    atom = AtomForm.parse("C")
    setattr(atom, field, old)
    source = Molecule.from_entries([atom])
    result = source.apply(edits)
    expected_atom = AtomForm.parse("C")
    setattr(expected_atom, field, new if count == 1 else old)
    assert result == Molecule.from_entries([expected_atom])
    assert source == Molecule.from_entries([atom])
    with pytest.raises(ConsumedError):
        getattr(edits, method)([])


@pytest.mark.parametrize(
    "key",
    [
        AtomConstraintKey.Valence(), AtomConstraintKey.DonatedPairs(),
        AtomConstraintKey.AcceptedPairs(), AtomConstraintKey.AromaticValence(),
        AtomConstraintKey.MulticenterValence(), AtomConstraintKey.TetrahedralStereo(),
        AtomConstraintKey.Degree(), AtomConstraintKey.TotalDegree(),
        AtomConstraintKey.TotalValence(), AtomConstraintKey.RingDegree(),
        AtomConstraintKey.RingValence(), AtomConstraintKey.TotalHydrogens(),
        AtomConstraintKey.RingMembership(RingScope.All()),
        AtomConstraintKey.RingMembership(RingScope.Size(5)),
    ],
)
@pytest.mark.parametrize("empty", [False, True])
def test_edits_atom_constraint_keys(key, empty):
    rows = [] if empty else [(0, None, None)]
    edit = Edit.ModifyAtomConstraints(key, rows)
    edits = Edits()
    edits.modify_atom_constraints(key, rows)
    rows.clear()
    assert list(edits) == [edit]
    assert edit.key == key
    assert edit.changes == (() if empty else ((0, None, None),))
    assert list(Edits.parse(edits.render())) == [edit]
    source = Molecule.from_entries([AtomForm.parse("C")])
    assert source.apply(edits) == source
    with pytest.raises(ConsumedError):
        edits.modify_atom_constraints(key, [])


def test_edits_atom_constraint_batches():
    first = AtomConstraintForm.Valence(NumForm.Lit(3))
    second = AtomConstraintForm.Valence(NumForm.LitSet({3, 4}))
    rows = [(0, None, first), (New(0), None, second), (0, first, second)]
    edit = Edit.ModifyAtomConstraints(AtomConstraintKey.Valence(), rows)
    edits = Edits()
    edits.add_atom(AtomForm.parse("N"))
    edits.modify_atom_constraints(AtomConstraintKey.Valence(), rows)
    rows.clear()
    assert edits[1] == edit
    assert list(Edits.parse(edits.render())) == list(edits)
    source = Molecule.from_entries([AtomForm.parse("C")])
    result = source.apply(edits)
    carbon = AtomForm.parse("C")
    nitrogen = AtomForm.parse("N")
    carbon.constraints.set(second)
    nitrogen.constraints.set(second)
    assert result == Molecule.from_entries([carbon, nitrogen])


def test_edits_atom_batches_new():
    edits = Edits()
    created = edits.add_atom(AtomForm.parse("N"))
    edits.modify_atom_charges([(0, NumForm.Undetermined(), NumForm.Lit(1)),
                               (created, NumForm.Undetermined(), NumForm.Lit(-1))])
    edits.modify_atom_lone_pairs([(created, NumForm.Undetermined(), NumForm.Lit(1))])
    source = Molecule.from_entries([AtomForm.parse("C")])
    carbon = AtomForm.parse("C")
    nitrogen = AtomForm.parse("N")
    carbon.charge = 1
    nitrogen.charge = -1
    nitrogen.lone_pairs = 1
    assert source.apply(edits) == Molecule.from_entries([carbon, nitrogen])


@pytest.mark.parametrize("method", ["modify_atom_charges", "modify_atom_lone_pairs"])
def test_edits_atom_batches_error(method):
    edits = Edits()
    getattr(edits, method)([(0, NumForm.Undetermined(), NumForm.Lit(1)),
                            (0, NumForm.Undetermined(), NumForm.Lit(2))])
    source = Molecule.from_entries([AtomForm.parse("C")])
    with pytest.raises(TransactionError):
        source.apply(edits)
    assert source == Molecule.from_entries([AtomForm.parse("C")])


@pytest.mark.parametrize(
    "variant,method,field,old,new",
    [
        (Edit.ModifyBondOrders, "modify_bond_orders", "order", NumForm.Lit(1), NumForm.LitSet({1, 2})),
        (Edit.ModifyBondCharges, "modify_bond_charges", "charge", NumForm.Lit(0), NumForm.Lit(-1)),
        (Edit.ModifyBondUnpairedElectrons, "modify_bond_unpaired_electrons", "unpaired_electrons",
         UnpairedElectronsForm(0, 1), UnpairedElectronsForm(1, 2)),
    ],
)
@pytest.mark.parametrize("count", [0, 1, 2])
def test_edits_bond_batches(variant, method, field, old, new, count):
    rows = [(0, old, new), (0, new, old)][:count]
    edit = variant(changes=rows)
    expected_rows = tuple(rows)
    rows.clear()
    assert edit.changes == expected_rows
    assert isinstance(edit.changes, tuple)
    assert all(isinstance(row, tuple) for row in edit.changes)
    with pytest.raises(AttributeError):
        edit.changes = ()
    edits = Edits()
    getattr(edits, method)(list(expected_rows))
    assert list(edits) == [edit]
    assert list(Edits.parse(edits.render())) == [edit]
    assert repr(edit).startswith("Edit." + variant.__name__ + "(changes=")
    atoms = [AtomForm.parse("C"), AtomForm.parse("N")]
    bond = BondForm.parse("1")
    setattr(bond, field, old)
    source = Molecule.from_entries(atoms, bonds=[(0, 1, bond)])
    result = source.apply(edits)
    expected = BondForm.parse("1")
    setattr(expected, field, new if count == 1 else old)
    assert result == Molecule.from_entries(atoms, bonds=[(0, 1, expected)])
    assert source == Molecule.from_entries(atoms, bonds=[(0, 1, bond)])
    with pytest.raises(ConsumedError):
        getattr(edits, method)([])


@pytest.mark.parametrize("key", [BondConstraintKey.Aromatic(), BondConstraintKey.CisTransStereo(),
                                  BondConstraintKey.RingMembership(RingScope.All()),
                                  BondConstraintKey.RingMembership(RingScope.Size(5))])
@pytest.mark.parametrize("empty", [False, True])
def test_edits_bond_constraint_keys(key, empty):
    rows = [] if empty else [(0, None, None)]
    edit = Edit.ModifyBondConstraints(key, rows)
    edits = Edits()
    edits.modify_bond_constraints(key, rows)
    rows.clear()
    assert list(edits) == [edit]
    assert edit.key == key
    assert edit.changes == (() if empty else ((0, None, None),))
    assert list(Edits.parse(edits.render())) == [edit]
    source = Molecule.from_entries([AtomForm.parse("C")] * 2, bonds=[(0, 1, BondForm.parse("1"))])
    assert source.apply(edits) == source
    with pytest.raises(ConsumedError):
        edits.modify_bond_constraints(key, [])


def test_edits_bond_batches_new():
    atoms = [AtomForm.parse("C")] * 3
    first = BondForm.parse("1")
    source = Molecule.from_entries(atoms, bonds=[(0, 1, first)])
    edits = Edits()
    created = edits.add_bond(1, 2, BondForm.parse("1"))
    edits.modify_bond_orders([(0, NumForm.Lit(1), NumForm.Lit(2)), (created, NumForm.Lit(1), NumForm.Lit(3))])
    edits.modify_bond_charges([(created, NumForm.Undetermined(), NumForm.Lit(-1))])
    edits.modify_bond_unpaired_electrons([(created, UnpairedElectronsForm(NumForm.Undetermined(), NumForm.Undetermined()), UnpairedElectronsForm(0, 1))])
    assert list(Edits.parse(edits.render())) == list(edits)
    first.order = 2
    second = BondForm.parse("3")
    second.charge = -1
    second.unpaired_electrons = UnpairedElectronsForm(0, 1)
    assert source.apply(edits) == Molecule.from_entries(atoms, bonds=[(0, 1, first), (1, 2, second)])


def test_edits_bond_constraint_batches():
    first = BondConstraintForm.Aromatic(BooleanForm.Lit(False))
    second = BondConstraintForm.Aromatic(BooleanForm.Undetermined())
    rows = [(0, None, first), (New(0), None, second), (0, first, second)]
    edit = Edit.ModifyBondConstraints(BondConstraintKey.Aromatic(), rows)
    edits = Edits()
    edits.add_bond(1, 2, BondForm.parse("1"))
    edits.modify_bond_constraints(BondConstraintKey.Aromatic(), rows)
    rows.clear()
    assert edits[1] == edit
    assert list(Edits.parse(edits.render())) == list(edits)
    atoms = [AtomForm.parse("C")] * 3
    source = Molecule.from_entries(atoms, bonds=[(0, 1, BondForm.parse("1"))])
    expected = BondForm.parse("1")
    expected.constraints.set(second)
    assert source.apply(edits) == Molecule.from_entries(atoms, bonds=[(0, 1, expected), (1, 2, expected)])


@pytest.mark.parametrize("key,old,new", [
    (BondConstraintKey.Aromatic(), BondConstraintForm.Aromatic(BooleanForm.Lit(False)), None),
    (BondConstraintKey.CisTransStereo(), None, BondConstraintForm.Aromatic(BooleanForm.Lit(True))),
])
def test_edits_bond_constraint_batches_error(key, old, new):
    edits = Edits()
    edits.modify_bond_orders([(0, NumForm.Lit(1), NumForm.Lit(2))])
    edits.modify_bond_constraints(key, [(0, old, new)])
    source = Molecule.from_entries([AtomForm.parse("C")] * 2, bonds=[(0, 1, BondForm.parse("1"))])
    with pytest.raises(TransactionError):
        source.apply(edits)
    assert source == Molecule.from_entries([AtomForm.parse("C")] * 2, bonds=[(0, 1, BondForm.parse("1"))])


@pytest.mark.parametrize(
    "variant,method,field,form_type,text,old,new,collection,participants",
    [
        (Edit.ModifyDativeBondOrders, "modify_dative_bond_orders", "order", DativeBondForm, "1", NumForm.Lit(1), NumForm.LitSet({1, 2}), "dative_bonds", ([2, 0], 1)),
        (Edit.ModifyAromaticSystemElectrons, "modify_aromatic_system_electrons", "electrons", AromaticSystemForm, "[1,2]", ElectronCountsForm.Lit([1, 2]), ElectronCountsForm.Lit([]), "aromatic_systems", ([2, 0, 1],)),
        (Edit.ModifyAromaticSystemCharges, "modify_aromatic_system_charges", "charge", AromaticSystemForm, "[1,2]", NumForm.Lit(1), NumForm.LitSet({1, 2}), "aromatic_systems", ([2, 0, 1],)),
        (Edit.ModifyAromaticSystemUnpairedElectrons, "modify_aromatic_system_unpaired_electrons", "unpaired_electrons", AromaticSystemForm, "[1,2]", UnpairedElectronsForm(0, 1), UnpairedElectronsForm(1, 2), "aromatic_systems", ([2, 0, 1],)),
        (Edit.ModifyMulticenterBondElectrons, "modify_multicenter_bond_electrons", "electrons", MulticenterBondForm, "[1,2]", ElectronCountsForm.Lit([1, 2]), ElectronCountsForm.Lit([]), "multicenter_bonds", ([2, 0, 1],)),
        (Edit.ModifyMulticenterBondCharges, "modify_multicenter_bond_charges", "charge", MulticenterBondForm, "[1,2]", NumForm.Lit(1), NumForm.LitSet({1, 2}), "multicenter_bonds", ([2, 0, 1],)),
        (Edit.ModifyMulticenterBondUnpairedElectrons, "modify_multicenter_bond_unpaired_electrons", "unpaired_electrons", MulticenterBondForm, "[1,2]", UnpairedElectronsForm(0, 1), UnpairedElectronsForm(1, 2), "multicenter_bonds", ([2, 0, 1],)),
        (Edit.ModifyNoncovalentBondKinds, "modify_noncovalent_bond_kinds", "kind", NoncovalentBondForm, "Hbd", NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond), NoncovalentBondKindForm.Undetermined(), "noncovalent_bonds", ((2, 0),)),
    ],
)
@pytest.mark.parametrize("count", [0, 1, 2])
def test_edits_relation_field_batches(variant, method, field, form_type, text, old, new, collection, participants, count):
    rows = [(0, old, new), (0, new, old)][:count]
    expected_rows = tuple(rows)
    edit = variant(changes=rows)
    rows.clear()
    assert edit.changes == expected_rows
    assert isinstance(edit.changes, tuple)
    assert all(isinstance(row, tuple) for row in edit.changes)
    with pytest.raises(AttributeError):
        edit.changes = ()
    edits = Edits()
    getattr(edits, method)(list(expected_rows))
    assert list(edits) == [edit]
    assert list(Edits.parse(edits.render())) == [edit]
    assert repr(edit).startswith("Edit." + variant.__name__ + "(changes=")
    atoms = [AtomForm.parse("C")] * 3
    form = form_type.parse(text)
    setattr(form, field, old)
    source = Molecule.from_entries(atoms, **{collection: [(*participants, form)]})
    expected = form_type.parse(text)
    setattr(expected, field, new if count == 1 else old)
    assert source.apply(edits) == Molecule.from_entries(atoms, **{collection: [(*participants, expected)]})
    assert source == Molecule.from_entries(atoms, **{collection: [(*participants, form)]})
    with pytest.raises(ConsumedError):
        getattr(edits, method)([])


@pytest.mark.parametrize(
    "variant,method,key,form_type,text,collection,participants",
    [
        (Edit.ModifyDativeBondConstraints, "modify_dative_bond_constraints", DativeBondConstraintKey.Aromatic(), DativeBondForm, "1", "dative_bonds", ([2, 0], 1)),
        (Edit.ModifyDativeBondConstraints, "modify_dative_bond_constraints", DativeBondConstraintKey.RingMembership(RingScope.All()), DativeBondForm, "1", "dative_bonds", ([2, 0], 1)),
        (Edit.ModifyDativeBondConstraints, "modify_dative_bond_constraints", DativeBondConstraintKey.RingMembership(RingScope.Size(5)), DativeBondForm, "1", "dative_bonds", ([2, 0], 1)),
        (Edit.ModifyAromaticSystemConstraints, "modify_aromatic_system_constraints", AromaticSystemConstraintKey.ElectronCount(), AromaticSystemForm, "[1,2]", "aromatic_systems", ([2, 0, 1],)),
        (Edit.ModifyMulticenterBondConstraints, "modify_multicenter_bond_constraints", MulticenterBondConstraintKey.ElectronCount(), MulticenterBondForm, "[1,2]", "multicenter_bonds", ([2, 0, 1],)),
        (Edit.ModifyNoncovalentBondConstraints, "modify_noncovalent_bond_constraints", NoncovalentBondConstraintKey.Intramolecular(), NoncovalentBondForm, "Hbd", "noncovalent_bonds", ((2, 0),)),
    ],
)
@pytest.mark.parametrize("empty", [False, True])
def test_edits_relation_constraint_keys(variant, method, key, form_type, text, collection, participants, empty):
    rows = [] if empty else [(0, None, None)]
    edit = variant(key, rows)
    edits = Edits()
    getattr(edits, method)(key, rows)
    rows.clear()
    assert list(edits) == [edit]
    assert edit.key == key
    assert edit.changes == (() if empty else ((0, None, None),))
    assert list(Edits.parse(edits.render())) == [edit]
    source = Molecule.from_entries([AtomForm.parse("C")] * 3, **{collection: [(*participants, form_type.parse(text))]})
    assert source.apply(edits) == source
    with pytest.raises(ConsumedError):
        getattr(edits, method)(key, [])


@pytest.mark.parametrize(
    "variant,method,add_method,key,first,second,form_type,text,collection,participants,created_participants",
    [
        (Edit.ModifyDativeBondConstraints, "modify_dative_bond_constraints", "add_dative_bond", DativeBondConstraintKey.Aromatic(), DativeBondConstraintForm.Aromatic(BooleanForm.Lit(False)), DativeBondConstraintForm.Aromatic(BooleanForm.Undetermined()), DativeBondForm, "1", "dative_bonds", ([2, 0], 1), ([5, 3], 4)),
        (Edit.ModifyDativeBondConstraints, "modify_dative_bond_constraints", "add_dative_bond", DativeBondConstraintKey.RingMembership(RingScope.Size(5)), DativeBondConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(5), NumForm.Lit(1))), DativeBondConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(5), NumForm.Undetermined())), DativeBondForm, "1", "dative_bonds", ([2, 0], 1), ([5, 3], 4)),
        (Edit.ModifyAromaticSystemConstraints, "modify_aromatic_system_constraints", "add_aromatic_system", AromaticSystemConstraintKey.ElectronCount(), AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(2)), AromaticSystemConstraintForm.ElectronCount(NumForm.LitSet({2, 3})), AromaticSystemForm, "[1,2]", "aromatic_systems", ([2, 0, 1],), ([5, 3, 4],)),
        (Edit.ModifyMulticenterBondConstraints, "modify_multicenter_bond_constraints", "add_multicenter_bond", MulticenterBondConstraintKey.ElectronCount(), MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(2)), MulticenterBondConstraintForm.ElectronCount(NumForm.LitSet({2, 3})), MulticenterBondForm, "[1,2]", "multicenter_bonds", ([2, 0, 1],), ([5, 3, 4],)),
        (Edit.ModifyNoncovalentBondConstraints, "modify_noncovalent_bond_constraints", "add_noncovalent_bond", NoncovalentBondConstraintKey.Intramolecular(), NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(False)), NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Undetermined()), NoncovalentBondForm, "Hbd", "noncovalent_bonds", ((2, 0),), ((5, 3),)),
    ],
)
def test_edits_relation_constraint_batches(variant, method, add_method, key, first, second, form_type, text, collection, participants, created_participants):
    rows = [(0, None, first), (New(0), None, second), (0, first, second), (0, None, None), (0, second, None)]
    expected_rows = tuple(rows)
    edit = variant(key, rows)
    edits = Edits()
    getattr(edits, add_method)(*created_participants, form_type.parse(text))
    getattr(edits, method)(key, rows)
    rows.clear()
    assert edits[1] == edit
    assert edit.key == key
    assert edit.changes == expected_rows
    assert isinstance(edit.changes, tuple)
    assert all(isinstance(row, tuple) for row in edit.changes)
    assert list(Edits.parse(edits.render())) == list(edits)
    atoms = [AtomForm.parse("C")] * 6
    source = Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text))]})
    expected = form_type.parse(text)
    expected.constraints.set(second)
    assert source.apply(edits) == Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text)), (*created_participants, expected)]})
    assert source == Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text))]})
    rejected = Edits()
    getattr(rejected, method)(key, [(0, None, first), (0, None, second)])
    with pytest.raises(TransactionError):
        source.apply(rejected)
    assert source == Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text))]})


@pytest.mark.parametrize(
    "method,add_method,field,old,new,form_type,text,collection,participants,created_participants",
    [
        ("modify_dative_bond_orders", "add_dative_bond", "order", NumForm.Lit(1), NumForm.Lit(2), DativeBondForm, "1", "dative_bonds", ([2, 0], 1), ([5, 3], 4)),
        ("modify_aromatic_system_electrons", "add_aromatic_system", "electrons", ElectronCountsForm.Lit([1, 2]), ElectronCountsForm.Lit([]), AromaticSystemForm, "[1,2]", "aromatic_systems", ([2, 0, 1],), ([5, 3, 4],)),
        ("modify_multicenter_bond_electrons", "add_multicenter_bond", "electrons", ElectronCountsForm.Lit([1, 2]), ElectronCountsForm.Lit([]), MulticenterBondForm, "[1,2]", "multicenter_bonds", ([2, 0, 1],), ([5, 3, 4],)),
        ("modify_noncovalent_bond_kinds", "add_noncovalent_bond", "kind", NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond), NoncovalentBondKindForm.Undetermined(), NoncovalentBondForm, "Hbd", "noncovalent_bonds", ((2, 0),), ((5, 3),)),
    ],
)
def test_edits_relation_field_batches_new(method, add_method, field, old, new, form_type, text, collection, participants, created_participants):
    atoms = [AtomForm.parse("C")] * 6
    source = Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text))]})
    edits = Edits()
    created = getattr(edits, add_method)(*created_participants, form_type.parse(text))
    getattr(edits, method)([(0, old, new), (created, old, new), (created, new, old)])
    assert list(Edits.parse(edits.render())) == list(edits)
    expected = form_type.parse(text)
    setattr(expected, field, new)
    assert source.apply(edits) == Molecule.from_entries(atoms, **{collection: [(*participants, expected), (*created_participants, form_type.parse(text))]})
    rejected = Edits()
    getattr(rejected, method)([(0, old, new), (0, old, new)])
    with pytest.raises(TransactionError):
        source.apply(rejected)
    assert source == Molecule.from_entries(atoms, **{collection: [(*participants, form_type.parse(text))]})
