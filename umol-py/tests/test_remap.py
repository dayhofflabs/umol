import pytest

from umol import MoleculeRemapping, Remapping


@pytest.mark.parametrize("images", [[], [0], [2, 0, 1]])
def test_remapping_constructor(images):
    remapping = Remapping(images)
    assert remapping.images == images
    assert len(remapping) == len(images)
    assert remapping.is_empty() == (images == [])
    assert remapping.to_correspondence().matched_pairs == list(enumerate(images))
    assert remapping.to_correspondence().left_count == len(images)
    assert remapping.to_correspondence().right_count == len(images)


@pytest.mark.parametrize(
    ("images", "message"),
    [
        ([1], r"^image NodeId\(1\) is out of range for 1 entries$"),
        ([1, 1], r"^image NodeId\(1\) occurs more than once$"),
    ],
)
def test_remapping_constructor_error(images, message):
    with pytest.raises(ValueError, match=message):
        Remapping(images)


@pytest.mark.parametrize(("source", "expected"), [(0, 2), (1, 0), (2, 1)])
def test_remapping_map(source, expected):
    remapping = Remapping([2, 0, 1])
    assert remapping.map(source) == expected
    assert remapping.try_map(source) == expected


def test_remapping_map_error():
    remapping = Remapping([0])
    assert remapping.try_map(1) is None
    with pytest.raises(ValueError, match="^id outside remapping source domain$"):
        remapping.map(1)


def test_remapping_frozen():
    remapping = Remapping([1, 0])
    images = remapping.images
    images[0] = 0
    assert remapping.images == [1, 0]
    with pytest.raises(AttributeError):
        remapping.images = [0, 1]


@pytest.fixture
def molecule_remapping():
    return MoleculeRemapping(
        Remapping([2, 0, 1]), Remapping([1, 0]), Remapping([0]), Remapping([]),
        Remapping([1, 0]), Remapping([0]), Remapping([0, 1]), Remapping([]),
    )


@pytest.mark.parametrize(
    ("field", "images"),
    [
        ("atoms", [2, 0, 1]), ("bonds", [1, 0]), ("dative_bonds", [0]),
        ("aromatic_systems", []), ("multicenter_bonds", [1, 0]),
        ("noncovalent_bonds", [0]), ("stereo_atoms", [0, 1]), ("stereo_bonds", []),
    ],
)
def test_molecule_remapping_components(molecule_remapping, field, images):
    assert getattr(molecule_remapping, field) == Remapping(images)
    correspondence = getattr(molecule_remapping.to_correspondence(), field)
    assert correspondence.matched_pairs == list(enumerate(images))
    assert correspondence.left_count == len(images)
    assert correspondence.right_count == len(images)


def test_molecule_remapping_frozen(molecule_remapping):
    with pytest.raises(AttributeError):
        molecule_remapping.atoms = Remapping([0, 1, 2])
