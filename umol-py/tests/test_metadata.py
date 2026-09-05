import pytest

from umol import (
    AtomForm,
    BondForm,
    Element,
    Entity,
    MetadataError,
    Molecule,
    MoleculeMetadata,
    ReactionMetadata,
)


@pytest.mark.parametrize(
    ("entity", "keyword"),
    [
        pytest.param(Entity.Atom(0), "atom", id="atom"),
        pytest.param(Entity.Bond(0), "bond", id="bond"),
        pytest.param(Entity.DativeBond(0), "dative", id="dative-bond"),
        pytest.param(
            Entity.AromaticSystem(0),
            "aromatic",
            id="aromatic-system",
        ),
        pytest.param(
            Entity.MulticenterBond(0),
            "multicenter",
            id="multicenter-bond",
        ),
        pytest.param(
            Entity.NoncovalentBond(0),
            "noncovalent",
            id="noncovalent-bond",
        ),
        pytest.param(Entity.StereoAtom(0), "stereo_atom", id="stereo-atom"),
        pytest.param(Entity.StereoBond(0), "stereo_bond", id="stereo-bond"),
    ],
)
def test_molecule_metadata_keyword(entity, keyword):
    metadata = MoleculeMetadata()

    metadata.set_keyword(entity, keyword)
    metadata.set_keyword(entity, keyword)

    assert metadata.keyword(entity) == keyword
    assert metadata.entity(keyword) == entity
    assert metadata.keyword(Entity.Atom(99)) is None
    assert metadata.entity("missing") is None


def test_molecule_metadata_set_keyword_error():
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(0), "site")

    with pytest.raises(MetadataError, match="^duplicate keyword: site$"):
        metadata.set_keyword(Entity.Bond(0), "site")

    assert metadata.entity("site") == Entity.Atom(0)
    assert metadata.keyword(Entity.Bond(0)) is None


def test_molecule_metadata_remap():
    source = Molecule.from_entries(
        [AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(2))],
    )
    combined = Molecule.from_entries(
        [AtomForm(Element("C"))]
    ).combine(source)
    _, source_to_component = combined.tracked_split()[1]
    correspondence = source_to_component.reverse()
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(0), "oxygen")
    metadata.set_keyword(Entity.Atom(1), "nitrogen")
    metadata.set_keyword(Entity.Bond(0), "bond")

    remapped = metadata.remap(correspondence)

    assert remapped.keyword(Entity.Atom(1)) == "oxygen"
    assert remapped.keyword(Entity.Atom(2)) == "nitrogen"
    assert remapped.keyword(Entity.Bond(0)) == "bond"
    assert remapped.entity("oxygen") == Entity.Atom(1)
    assert metadata.entity("oxygen") == Entity.Atom(0)


def test_molecule_metadata_repr():
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(0), "carbon")

    assert repr(MoleculeMetadata()) == (
        "MoleculeMetadata(keywords=[], atom_alias_count=0)"
    )
    assert repr(metadata) == (
        'MoleculeMetadata(keywords=[(Entity.Atom(0), "carbon")], '
        "atom_alias_count=0)"
    )


@pytest.mark.parametrize(
    ("entity", "keyword"),
    [
        pytest.param(Entity.Atom(0), "atom", id="atom"),
        pytest.param(Entity.Bond(0), "bond", id="bond"),
        pytest.param(Entity.DativeBond(0), "dative", id="dative-bond"),
        pytest.param(
            Entity.AromaticSystem(0),
            "aromatic",
            id="aromatic-system",
        ),
        pytest.param(
            Entity.MulticenterBond(0),
            "multicenter",
            id="multicenter-bond",
        ),
        pytest.param(
            Entity.NoncovalentBond(0),
            "noncovalent",
            id="noncovalent-bond",
        ),
        pytest.param(Entity.StereoAtom(0), "stereo_atom", id="stereo-atom"),
        pytest.param(Entity.StereoBond(0), "stereo_bond", id="stereo-bond"),
    ],
)
def test_reaction_metadata_delta_keyword(entity, keyword):
    metadata = ReactionMetadata()

    metadata.set_delta_keyword(entity, keyword)
    metadata.set_delta_keyword(entity, keyword)

    assert metadata.keyword(entity) == keyword
    assert metadata.entity(keyword) == entity
    assert metadata.delta_keyword(entity) == keyword
    assert metadata.delta_entity(keyword) == entity


def test_reaction_metadata_scope():
    lhs = MoleculeMetadata()
    lhs.set_keyword(Entity.Atom(0), "lhs")
    metadata = ReactionMetadata(lhs)
    metadata.set_delta_keyword(Entity.Bond(0), "delta")

    detached_lhs = metadata.lhs
    detached_lhs.set_keyword(Entity.Atom(1), "detached")

    assert metadata.keyword(Entity.Atom(0)) == "lhs"
    assert metadata.entity("lhs") == Entity.Atom(0)
    assert metadata.delta_keyword(Entity.Atom(0)) is None
    assert metadata.delta_entity("lhs") is None
    assert metadata.keyword(Entity.Bond(0)) == "delta"
    assert metadata.entity("delta") == Entity.Bond(0)
    assert metadata.delta_keyword(Entity.Bond(0)) == "delta"
    assert metadata.delta_entity("delta") == Entity.Bond(0)
    assert metadata.entity("detached") is None


def test_reaction_metadata_set_delta_keyword_error():
    lhs = MoleculeMetadata()
    lhs.set_keyword(Entity.Atom(0), "site")
    metadata = ReactionMetadata(lhs)

    with pytest.raises(MetadataError, match="^duplicate keyword: site$"):
        metadata.set_delta_keyword(Entity.Bond(0), "site")

    assert metadata.entity("site") == Entity.Atom(0)
    assert metadata.delta_entity("site") is None
    assert metadata.delta_keyword(Entity.Bond(0)) is None


def test_reaction_metadata_repr():
    lhs = MoleculeMetadata()
    lhs.set_keyword(Entity.Atom(0), "lhs")
    metadata = ReactionMetadata(lhs)
    metadata.set_delta_keyword(Entity.Bond(1), "delta")

    assert repr(ReactionMetadata()) == (
        "ReactionMetadata(lhs=MoleculeMetadata(keywords=[], "
        "atom_alias_count=0), delta_keywords=[], "
        "reaction_atom_alias_count=0)"
    )
    assert repr(metadata) == (
        "ReactionMetadata(lhs=MoleculeMetadata("
        'keywords=[(Entity.Atom(0), "lhs")], atom_alias_count=0), '
        'delta_keywords=[(Entity.Bond(1), "delta")], '
        "reaction_atom_alias_count=0)"
    )
