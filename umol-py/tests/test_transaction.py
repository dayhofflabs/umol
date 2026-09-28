import pytest

from umol import (
    AromaticSystemForm,
    AtomForm,
    BondForm,
    ConsumedError,
    Correspondence,
    DativeBondForm,
    Edit,
    Edits,
    Element,
    Molecule,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    StereoAtomForm,
    StereoBondForm,
    StereoConfigurationForm,
    StereoLigand,
    StereoLigandKind,
    TetrahedralConfiguration,
    TransactionError,
)


def rich_molecule():
    return Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(5)],
        bonds=[
            (0, 1, BondForm(2)),
            (0, 2, BondForm(1)),
            (0, 3, BondForm(1)),
            (0, 4, BondForm(1)),
            (1, 3, BondForm(1)),
        ],
        dative_bonds=[([2], 1, DativeBondForm(1))],
        aromatic_systems=[([0, 1, 2], AromaticSystemForm([1, 1, 1]))],
        multicenter_bonds=[([0, 1, 2], MulticenterBondForm([1, 1, 1]))],
        noncovalent_bonds=[
            ([0, 2], NoncovalentBondForm(NoncovalentBondKind.HydrogenBond))
        ],
        stereo_atoms=[
            (
                0,
                [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
                StereoAtomForm(TetrahedralConfiguration.Ccw),
            )
        ],
        stereo_bonds=[
            (
                0,
                [
                    StereoLigand(2, StereoLigandKind.Atom),
                    StereoLigand(0, StereoLigandKind.ImplicitHydrogen),
                    StereoLigand(3, StereoLigandKind.Atom),
                    StereoLigand(1, StereoLigandKind.ImplicitHydrogen),
                ],
                StereoBondForm.parse("Ct0"),
            )
        ],
    )


def add_carbon_edits():
    return Edits.parse('[{:atom {:add "C"}}]')


@pytest.fixture(params=[
    pytest.param((
        Edit.ReplaceDativeBondDonors(id=0, old=[0], new=[2, 3]),
        {"dative_bonds": [([0], 1, DativeBondForm(1))]},
        {"dative_bonds": [([2, 3], 1, DativeBondForm(1))]},
    ), id="dative_bond_donors"),
    pytest.param((
        Edit.ReplaceDativeBondAcceptor(id=0, old=1, new=2),
        {"dative_bonds": [([0], 1, DativeBondForm(1))]},
        {"dative_bonds": [([0], 2, DativeBondForm(1))]},
    ), id="dative_bond_acceptor"),
    pytest.param((
        Edit.ReplaceAromaticSystemAtoms(id=0, old=[0, 1], new=[2, 3]),
        {"aromatic_systems": [([0, 1], AromaticSystemForm([1, 2]))]},
        {"aromatic_systems": [([2, 3], AromaticSystemForm([1, 2]))]},
    ), id="aromatic_system_atoms"),
    pytest.param((
        Edit.ReplaceMulticenterBondAtoms(id=0, old=[0, 1], new=[2, 3]),
        {"multicenter_bonds": [([0, 1], MulticenterBondForm([1, 2]))]},
        {"multicenter_bonds": [([2, 3], MulticenterBondForm([1, 2]))]},
    ), id="multicenter_bond_atoms"),
    pytest.param((
        Edit.ReplaceNoncovalentBondAtoms(id=0, old=(0, 1), new=(2, 3)),
        {"noncovalent_bonds": [
            ((0, 1), NoncovalentBondForm(NoncovalentBondKind.HydrogenBond))
        ]},
        {"noncovalent_bonds": [
            ((2, 3), NoncovalentBondForm(NoncovalentBondKind.HydrogenBond))
        ]},
    ), id="noncovalent_bond_atoms"),
    pytest.param((
        Edit.ReplaceStereoAtomSite(id=0, old=0, new=2),
        {"stereo_atoms": [(
            0, [StereoLigand(1, StereoLigandKind.Atom)],
            StereoAtomForm(StereoConfigurationForm.Undetermined()),
        )]},
        {"stereo_atoms": [(
            2, [StereoLigand(1, StereoLigandKind.Atom)],
            StereoAtomForm(StereoConfigurationForm.Undetermined()),
        )]},
    ), id="stereo_atom_site"),
    pytest.param((
        Edit.ReplaceStereoAtomLigands(
            id=0, old=[(3, StereoLigandKind.Atom)], new=[(4, StereoLigandKind.Atom)],
        ),
        {"stereo_atoms": [(
            0, [StereoLigand(3, StereoLigandKind.Atom)],
            StereoAtomForm(StereoConfigurationForm.Undetermined()),
        )]},
        {"stereo_atoms": [(
            0, [StereoLigand(4, StereoLigandKind.Atom)],
            StereoAtomForm(StereoConfigurationForm.Undetermined()),
        )]},
    ), id="stereo_atom_ligands"),
    pytest.param((
        Edit.ReplaceStereoBondSite(id=0, old=0, new=1),
        {"stereo_bonds": [(
            0, [StereoLigand(i, StereoLigandKind.Atom) for i in [3, 4, 5, 6]],
            StereoBondForm.parse("Ct0"),
        )]},
        {"stereo_bonds": [(
            1, [StereoLigand(i, StereoLigandKind.Atom) for i in [3, 4, 5, 6]],
            StereoBondForm.parse("Ct0"),
        )]},
    ), id="stereo_bond_site"),
    pytest.param((
        Edit.ReplaceStereoBondLigands(
            id=0,
            old=[(i, StereoLigandKind.Atom) for i in [3, 4, 5, 6]],
            new=[(0, StereoLigandKind.ImplicitHydrogen)]
                + [(i, StereoLigandKind.Atom) for i in [4, 5, 6]],
        ),
        {"stereo_bonds": [(
            0, [StereoLigand(i, StereoLigandKind.Atom) for i in [3, 4, 5, 6]],
            StereoBondForm.parse("Ct0"),
        )]},
        {"stereo_bonds": [(
            0, [StereoLigand(0, StereoLigandKind.ImplicitHydrogen)]
                + [StereoLigand(i, StereoLigandKind.Atom) for i in [4, 5, 6]],
            StereoBondForm.parse("Ct0"),
        )]},
    ), id="stereo_bond_ligands"),
])
def replacement_case(request):
    edit, original, changed = request.param
    atoms = [AtomForm.parse("C") for _ in range(7)]
    bonds = [(a, b, BondForm(1)) for a, b in [
        (0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (2, 5), (2, 6), (1, 2),
    ]]
    return (
        Molecule.from_entries(atoms, bonds=bonds, **original),
        Molecule.from_entries(atoms, bonds=bonds, **changed),
        edit,
    )


def test_molecule_editor_tracked_snapshot_and_build():
    molecule = Molecule.parse('{:atoms ["N#h3"]}')
    editor = molecule.edit()

    snapshot, correspondence = editor.tracked_snapshot()

    assert snapshot == editor.snapshot() == molecule
    assert correspondence.atoms == Correspondence([(0, 0)], 1, 1)

    plain = molecule.edit().build()
    tracked, correspondence = molecule.edit().tracked_build()

    assert tracked == plain
    assert correspondence.atoms == Correspondence([(0, 0)], 1, 1)


@pytest.mark.parametrize("method", ["build", "tracked_build", "apply", "tracked_apply"])
def test_molecule_editor_consumed(method):
    molecule = Molecule.parse('{:atoms ["C"]}')
    editor = molecule.edit()
    alias = editor
    if method in ("apply", "tracked_apply"):
        result = getattr(editor, method)(Edits())
        if method == "tracked_apply":
            result = result[0]
        assert result.build() == molecule
    else:
        result = getattr(editor, method)()
        assert (result[0] if method == "tracked_build" else result) == molecule

    with pytest.raises(ConsumedError, match="^MoleculeEditor has been consumed$") as error:
        alias.snapshot()
    assert type(error.value) is ConsumedError


def test_molecule_editor_tracked_apply():
    molecule = Molecule.parse('{:atoms ["N#h3"]}')

    plain = molecule.edit().apply(add_carbon_edits())
    tracked, correspondence = molecule.edit().tracked_apply(add_carbon_edits())

    assert tracked.build() == plain.build()
    assert correspondence.atoms == Correspondence([(0, 0)], 1, 2)


def test_molecule_editor_tracked_transact_and_rollback():
    molecule = Molecule.parse('{:atoms ["N#h3"]}')
    plain_editor = molecule.edit()
    tracked_editor = molecule.edit()

    plain_transaction = plain_editor.transact(add_carbon_edits())
    tracked_transaction, forward = tracked_editor.tracked_transact(add_carbon_edits())

    assert tracked_editor.snapshot() == plain_editor.snapshot()
    assert forward.atoms == Correspondence([(0, 0)], 1, 2)

    plain_transaction.rollback(plain_editor)
    reverse = tracked_transaction.tracked_rollback(tracked_editor)

    assert tracked_editor.snapshot() == plain_editor.snapshot() == molecule
    assert reverse.atoms == Correspondence([(0, 0)], 2, 1)


@pytest.mark.parametrize("method", ["rollback", "tracked_rollback"])
def test_transaction_consumed(method):
    molecule = Molecule.parse('{:atoms ["N"]}')
    editor = molecule.edit()
    edits = Edits()
    edits.add_atom(AtomForm.parse("C"))
    transaction = editor.transact(edits)
    alias = transaction

    getattr(transaction, method)(editor)

    assert editor.snapshot() == molecule
    with pytest.raises(ConsumedError, match="^Transaction has been consumed$") as error:
        getattr(alias, method)(editor)
    assert type(error.value) is ConsumedError


def test_molecule_apply_replace(replacement_case):
    original, expected, edit = replacement_case

    assert original.apply(Edits([edit])) == expected


def test_molecule_editor_transact_replace(replacement_case):
    original, expected, edit = replacement_case
    editor = original.edit()

    transaction = editor.transact(Edits([edit]))

    assert editor.snapshot() == expected
    transaction.rollback(editor)
    assert editor.build() == original


def test_molecule_editor_transact_replace_new_handles():
    original = Molecule.parse('{:atoms ["C" "N"]}')
    editor = original.edit()
    edits = Edits()
    atom = edits.add_atom(AtomForm.parse("O"))
    dative = edits.add_dative_bond([0], 1, DativeBondForm(1))
    edits.append(Edit.ReplaceDativeBondDonors(id=dative, old=[0], new=[atom]))
    edits.append(Edit.ReplaceDativeBondAcceptor(id=dative, old=1, new=0))
    expected = Molecule.from_entries(
        [AtomForm.parse("C"), AtomForm.parse("N"), AtomForm.parse("O")],
        dative_bonds=[([2], 0, DativeBondForm(1))],
    )

    transaction = editor.transact(edits)

    assert editor.snapshot() == expected
    transaction.rollback(editor)
    assert editor.build() == original


@pytest.mark.parametrize("method", ["apply", "transact"])
def test_molecule_editor_replace_error(replacement_case, method):
    original, _, edit = replacement_case
    editor = original.edit()
    mismatch = type(edit)(id=edit.id, old=edit.new, new=edit.old)
    edits = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("N")]), mismatch])

    with pytest.raises(
        TransactionError,
        match="^precondition failed: old state does not match current$",
    ):
        getattr(editor, method)(edits)
    if method == "transact":
        assert editor.build() == original


def test_molecule_editor_remove_topology():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(1)), (1, 2, BondForm(1))],
    )
    plain = molecule.edit()
    tracked = molecule.edit()

    plain.remove_topology([1], [])
    tracked.remove_topology([1], [])
    result, correspondence = tracked.tracked_build()

    assert result == plain.build() == Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("N"))]
    )
    assert correspondence.atoms == Correspondence([(0, 0), (2, 1)], 3, 2)
    assert correspondence.bonds == Correspondence([], 2, 0)


@pytest.mark.parametrize(
    ("method", "field"),
    [
        ("remove_dative_bonds", "dative_bonds"),
        ("remove_aromatic_systems", "aromatic_systems"),
        ("remove_multicenter_bonds", "multicenter_bonds"),
        ("remove_noncovalent_bonds", "noncovalent_bonds"),
        ("remove_stereo_atoms", "stereo_atoms"),
        ("remove_stereo_bonds", "stereo_bonds"),
    ],
)
def test_molecule_editor_remove_entity_family(method, field):
    plain = rich_molecule().edit()
    tracked = rich_molecule().edit()

    getattr(plain, method)([0])
    getattr(tracked, method)([0])
    result, correspondence = tracked.tracked_build()

    assert result == plain.build()
    assert getattr(correspondence, field) == Correspondence([], 1, 0)


@pytest.mark.parametrize(
    ("method", "arguments", "message"),
    [
        ("remove_topology", ([5], []), "atom id out of range"),
        ("remove_topology", ([], [5]), "bond id out of range"),
        ("remove_dative_bonds", ([1],), "dative bond id out of range"),
        ("remove_aromatic_systems", ([1],), "aromatic system id out of range"),
        ("remove_multicenter_bonds", ([1],), "multicenter bond id out of range"),
        ("remove_noncovalent_bonds", ([1],), "noncovalent bond id out of range"),
        ("remove_stereo_atoms", ([1],), "stereo atom id out of range"),
        ("remove_stereo_bonds", ([1],), "stereo bond id out of range"),
    ],
)
def test_molecule_editor_remove_error(method, arguments, message):
    editor = rich_molecule().edit()

    with pytest.raises(IndexError, match=f"^{message}$"):
        getattr(editor, method)(*arguments)

    assert editor.snapshot() == rich_molecule()
