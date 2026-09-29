import pytest

from umol import (
    AromaticSystemForm,
    AtomForm,
    AtomFieldChange,
    BondForm,
    ConsumedError,
    Correspondence,
    DativeBondForm,
    Edit,
    Edits,
    Element,
    InvalidatedViewError,
    InvalidStructureError,
    Molecule,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NumForm,
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


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
@pytest.mark.parametrize("container", [list, tuple, iter], ids=["list", "tuple", "iterator"])
def test_molecule_transact(method, container):
    molecule = Molecule.parse('{:atoms ["C"]}')
    batches = []
    for element in ("N", "O"):
        edits = Edits()
        atom = edits.add_atom(AtomForm.parse(element))
        edits.add_bond(0, atom, BondForm(1))
        batches.append(edits)
    aliases = list(batches)
    iterator = iter(batches[0])
    saved_edit = next(iterator)

    result = getattr(molecule, method)(container(batches))

    assert molecule == Molecule.parse(
        '{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [0 2 "1"]]}'
    )
    if method == "tracked_transact":
        assert result.atoms == Correspondence([(0, 0)], 1, 3)
        assert result.bonds == Correspondence([], 0, 2)
    else:
        assert result is None
    for batch in aliases:
        with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
            len(batch)
    with pytest.raises(InvalidatedViewError, match="Edits has been consumed"):
        next(iterator)
    assert saved_edit == Edit.AddAtoms(atoms=[AtomForm.parse("N")])


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
@pytest.mark.parametrize("failure", ["consumed", "alias"])
def test_molecule_transact_consumption_error(method, failure):
    molecule = Molecule.parse('{:atoms ["C"]}')
    view = molecule.atoms[0]
    first = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("N")])])
    untouched = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("O")])])
    iterator = iter(first)
    if failure == "consumed":
        second = Edits()
        Molecule().transact([second])
    else:
        second = first

    with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
        getattr(molecule, method)([first, second, untouched])

    assert molecule == Molecule.parse('{:atoms ["C"]}')
    assert view.id == 0
    assert list(untouched) == [Edit.AddAtoms(atoms=[AtomForm.parse("O")])]
    for batch in (first, second):
        with pytest.raises(ConsumedError, match="^Edits has been consumed$"):
            len(batch)
    with pytest.raises(InvalidatedViewError, match="Edits has been consumed"):
        next(iterator)


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_argument_error(method):
    molecule = Molecule.parse('{:atoms ["C"]}')
    view = molecule.atoms[0]
    edit = Edit.AddAtoms(atoms=[AtomForm.parse("N")])
    first = Edits([edit])

    with pytest.raises(TypeError):
        getattr(molecule, method)(iter([first, None]))

    assert list(first) == [edit]
    assert view.id == 0
    assert molecule == Molecule.parse('{:atoms ["C"]}')


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
@pytest.mark.parametrize("fail", [False, True], ids=["complete", "iteration_error"])
def test_molecule_transact_generator(method, fail):
    molecule = Molecule.parse('{:atoms ["C"]}')
    view = molecule.atoms[0]
    edit = Edit.AddAtoms(atoms=[AtomForm.parse("N")])
    first = Edits([edit])

    def batches():
        yield first
        assert list(first) == [edit]
        assert view.id == 0
        if fail:
            raise ValueError("batch preparation failed")
        yield Edits([Edit.AddAtoms(atoms=[AtomForm.parse("O")])])

    if fail:
        with pytest.raises(ValueError, match="^batch preparation failed$"):
            getattr(molecule, method)(batches())
        assert list(first) == [edit]
        assert view.id == 0
        assert molecule == Molecule.parse('{:atoms ["C"]}')
    else:
        getattr(molecule, method)(batches())
        assert molecule == Molecule.parse('{:atoms ["C" "N" "O"]}')
        with pytest.raises(ConsumedError):
            len(first)
        with pytest.raises(InvalidatedViewError):
            view.id


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_fields(method):
    molecule = Molecule.from_entries([AtomForm(Element("C"), charge=0)])
    view = molecule.atoms[0]
    constraints = view.constraints
    bonds = molecule.bonds
    first = Edits([Edit.ModifyAtomField(
        id=0, change=AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(1)),
    )])
    second = Edits([Edit.ModifyAtomField(
        id=0, change=AtomFieldChange.Charge(old=NumForm.Lit(1), new=NumForm.Lit(2)),
    )])

    result = getattr(molecule, method)([first, second])

    assert molecule == Molecule.from_entries([AtomForm(Element("C"), charge=2)])
    assert molecule.atoms[0].charge == NumForm.Lit(2)
    if method == "tracked_transact":
        assert result.atoms == Correspondence([(0, 0)], 1, 1)
    for access in (lambda: view.id, lambda: len(constraints), lambda: len(bonds)):
        with pytest.raises(InvalidatedViewError):
            access()


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_error(method):
    molecule = Molecule.parse('{:atoms ["C"]}')
    view = molecule.atoms[0]
    first = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("N")])])
    second = Edits([Edit.ModifyAtomField(
        id=7, change=AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(1)),
    )])
    third = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("O")])])

    with pytest.raises(TransactionError, match="^atom handle 7 is out of range for 2 entries$"):
        getattr(molecule, method)([first, second, third])

    assert molecule == Molecule.parse('{:atoms ["C"]}')
    assert molecule.atoms[0].id == 0
    with pytest.raises(InvalidatedViewError):
        view.id
    for edits in (first, second, third):
        with pytest.raises(ConsumedError):
            len(edits)


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_integrity(method):
    molecule = Molecule.parse('{:atoms ["C" "N"] :bonds [[0 1 "1"]]}')
    first = Edits()
    first.add_bond(0, 1, BondForm(2))
    second = Edits([Edit.RemoveTopology(atoms=[], bonds=[0])])

    result = getattr(molecule, method)([first, second])

    assert molecule == Molecule.parse('{:atoms ["C" "N"] :bonds [[0 1 "2"]]}')
    if method == "tracked_transact":
        assert result.atoms == Correspondence([(0, 0), (1, 1)], 2, 2)
        assert result.bonds == Correspondence([], 1, 1)


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_integrity_error(method):
    molecule = Molecule.parse('{:atoms ["C" "N"] :bonds [[0 1 "1"]]}')
    view = molecule.bonds[0]
    edits = Edits()
    edits.add_bond(0, 1, BondForm(1))

    with pytest.raises(
        InvalidStructureError,
        match=r"^bond: parallel bonds on atoms \[AtomId\(0\), AtomId\(1\)\]$",
    ):
        getattr(molecule, method)([edits])

    assert molecule == Molecule.parse('{:atoms ["C" "N"] :bonds [[0 1 "1"]]}')
    assert molecule.bonds[0].id == 0
    with pytest.raises(InvalidatedViewError):
        view.id
    with pytest.raises(ConsumedError):
        len(edits)


def test_molecule_tracked_transact_compaction():
    molecule = Molecule.parse('{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')
    first = Edits([Edit.RemoveTopology(atoms=[0], bonds=[])])
    second = Edits()
    atom = second.add_atom(AtomForm.parse("F"))
    second.add_bond(1, atom, BondForm(1))

    correspondence = molecule.tracked_transact([first, second])

    assert molecule == Molecule.parse(
        '{:atoms ["N" "O" "F"] :bonds [[0 1 "1"] [1 2 "1"]]}'
    )
    assert correspondence.atoms == Correspondence([(1, 0), (2, 1)], 3, 3)
    assert correspondence.bonds == Correspondence([(1, 0)], 2, 2)


def test_molecule_apply_replace(replacement_case):
    original, expected, edit = replacement_case

    assert original.apply(Edits([edit])) == expected


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_replace(replacement_case, method):
    molecule, expected, edit = replacement_case

    getattr(molecule, method)([Edits([edit])])

    assert molecule == expected


@pytest.mark.parametrize("method", ["transact", "tracked_transact"])
def test_molecule_transact_replace_error(replacement_case, method):
    molecule, _, edit = replacement_case
    original = Molecule.parse(str(molecule))
    first = Edits([edit])
    second = Edits([edit])

    with pytest.raises(
        TransactionError, match="^precondition failed: old state does not match current$",
    ):
        getattr(molecule, method)([first, second])

    assert molecule == original


@pytest.mark.parametrize("rollback", [False, True], ids=["commit", "rollback"])
def test_molecule_transact_replace_new_handles(rollback):
    molecule = Molecule.parse('{:atoms ["C" "N"]}')
    edits = Edits()
    atom = edits.add_atom(AtomForm.parse("O"))
    dative = edits.add_dative_bond([0], 1, DativeBondForm(1))
    edits.append(Edit.ReplaceDativeBondDonors(id=dative, old=[0], new=[atom]))
    edits.append(Edit.ReplaceDativeBondAcceptor(id=dative, old=1, new=0))
    expected = Molecule.from_entries(
        [AtomForm.parse("C"), AtomForm.parse("N"), AtomForm.parse("O")],
        dative_bonds=[([2], 0, DativeBondForm(1))],
    )

    if rollback:
        with pytest.raises(
            TransactionError, match="^atom handle 7 is out of range for 3 entries$",
        ):
            molecule.transact([edits, Edits([Edit.RemoveTopology(atoms=[7], bonds=[])])])
        assert molecule == Molecule.parse('{:atoms ["C" "N"]}')
    else:
        molecule.transact([edits])
        assert molecule == expected


def test_molecule_editor_apply_replace_error(replacement_case):
    original, _, edit = replacement_case
    editor = original.edit()
    mismatch = type(edit)(id=edit.id, old=edit.new, new=edit.old)
    edits = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("N")]), mismatch])

    with pytest.raises(
        TransactionError,
        match="^precondition failed: old state does not match current$",
    ):
        editor.apply(edits)


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
