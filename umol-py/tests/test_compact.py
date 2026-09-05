import pytest

from umol import Compaction, Correspondence, MoleculeCompaction


def test_compaction_empty():
    assert Compaction.empty() == Compaction(0, [])


@pytest.mark.parametrize("source_count", [0, 1, 4])
def test_compaction_identity(source_count):
    assert Compaction.identity(source_count) == Compaction(source_count, [])


def test_compaction_constructor():
    compaction = Compaction(5, [3, 1, 3])

    assert compaction.source_count == 5
    assert compaction.result_count == 3
    assert compaction.removed == [1, 3]
    assert repr(compaction) == "Compaction(source_count=5, removed=[1, 3])"


def test_compaction_constructor_error():
    with pytest.raises(
        ValueError,
        match=r"^removed id NodeId\(2\) is out of range for 2 entries$",
    ):
        Compaction(2, [2])


@pytest.mark.parametrize(
    ("source", "result"),
    [(0, 0), (1, None), (2, 1), (3, None), (4, 2), (5, None)],
)
def test_compaction_compact(source, result):
    assert Compaction(5, [1, 3]).compact(source) == result


@pytest.mark.parametrize(("result", "source"), [(0, 0), (1, 2), (2, 4)])
def test_compaction_uncompact(result, source):
    compaction = Compaction(5, [1, 3])

    assert compaction.try_uncompact(result) == source
    assert compaction.uncompact(result) == source


def test_compaction_uncompact_error():
    compaction = Compaction(5, [1, 3])

    assert compaction.try_uncompact(3) is None
    with pytest.raises(ValueError, match="^id outside compaction result domain$"):
        compaction.uncompact(3)


def test_compaction_to_correspondence():
    assert Compaction(5, [1, 3]).to_correspondence() == Correspondence(
        [(0, 0), (2, 1), (4, 2)], 5, 3
    )


def test_compaction_frozen():
    compaction = Compaction(3, [1])
    removed = compaction.removed
    removed.append(2)

    assert compaction.removed == [1]
    with pytest.raises(AttributeError):
        compaction.removed = []


def molecule_compaction():
    return MoleculeCompaction(
        Compaction(4, [1, 3]),
        Compaction(2, [0]),
        Compaction.identity(1),
        Compaction.empty(),
        Compaction(2, [1]),
        Compaction.identity(1),
        Compaction(2, [0]),
        Compaction.empty(),
    )


@pytest.mark.parametrize(
    ("field", "source_count", "removed"),
    [
        ("atoms", 4, [1, 3]),
        ("bonds", 2, [0]),
        ("dative_bonds", 1, []),
        ("aromatic_systems", 0, []),
        ("multicenter_bonds", 2, [1]),
        ("noncovalent_bonds", 1, []),
        ("stereo_atoms", 2, [0]),
        ("stereo_bonds", 0, []),
    ],
)
def test_molecule_compaction_components(field, source_count, removed):
    component = getattr(molecule_compaction(), field)

    assert component == Compaction(source_count, removed)


@pytest.mark.parametrize(
    ("method", "source", "result"),
    [
        ("compact_atom", 2, 1),
        ("compact_bond", 0, None),
        ("compact_dative_bond", 0, 0),
        ("compact_aromatic_system", 0, None),
        ("compact_multicenter_bond", 1, None),
        ("compact_noncovalent_bond", 0, 0),
        ("compact_stereo_atom", 1, 0),
        ("compact_stereo_bond", 0, None),
    ],
)
def test_molecule_compaction_compact(method, source, result):
    assert getattr(molecule_compaction(), method)(source) == result


def test_molecule_compaction_to_correspondence():
    compaction = molecule_compaction()
    correspondence = compaction.to_correspondence()

    assert correspondence.atoms == compaction.atoms.to_correspondence()
    assert correspondence.bonds == compaction.bonds.to_correspondence()
    assert correspondence.stereo_atoms == compaction.stereo_atoms.to_correspondence()


@pytest.mark.parametrize(
    "field",
    [
        "atoms",
        "bonds",
        "dative_bonds",
        "aromatic_systems",
        "multicenter_bonds",
        "noncovalent_bonds",
        "stereo_atoms",
        "stereo_bonds",
    ],
)
def test_molecule_compaction_empty(field):
    assert getattr(MoleculeCompaction.empty(), field) == Compaction.empty()


def test_molecule_compaction_frozen():
    with pytest.raises(AttributeError):
        molecule_compaction().atoms = Compaction.identity(4)
