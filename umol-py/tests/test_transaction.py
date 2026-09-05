import pytest

from umol import (
    AromaticSystemForm,
    AtomForm,
    BondForm,
    Correspondence,
    DativeBondForm,
    Edits,
    Element,
    Molecule,
    MoleculeCompaction,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    StereoAtomForm,
    StereoBondForm,
    StereoLigand,
    StereoLigandKind,
    TetrahedralConfiguration,
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


def test_molecule_editor_tracked_remove():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(1)), (1, 2, BondForm(1))],
    )
    plain = molecule.edit()
    tracked = molecule.edit()

    plain.remove([1], [])
    compaction = tracked.tracked_remove([1], [])

    assert tracked.build() == plain.build()
    assert isinstance(compaction, MoleculeCompaction)
    assert compaction.atoms.removed == [1]
    assert compaction.bonds.removed == [0, 1]


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
def test_molecule_editor_tracked_remove_entity_family(method, field):
    plain = rich_molecule().edit()
    tracked = rich_molecule().edit()

    getattr(plain, method)([0])
    compaction = getattr(tracked, f"tracked_{method}")([0])

    assert tracked.build() == plain.build()
    assert getattr(compaction, field).removed == [0]


@pytest.mark.parametrize(
    ("method", "arguments", "message"),
    [
        ("remove", ([5], []), "atom id out of range"),
        ("remove", ([], [5]), "bond id out of range"),
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
