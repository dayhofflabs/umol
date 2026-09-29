import pytest

from umol import (
    AromaticSystemForm,
    AtomForm,
    BondForm,
    DativeBondForm,
    Edits,
    Element,
    InvalidatedViewError,
    Molecule,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NumForm,
    StereoAtomForm,
    StereoBondForm,
    StereoLigand,
    StereoLigandKind,
    TetrahedralConfiguration,
)


@pytest.fixture
def molecule():
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


@pytest.fixture(params=[
    ("atoms", "charge", AtomForm(Element("N"))),
    ("bonds", "order", BondForm(1)),
    ("dative_bonds", "order", DativeBondForm(2)),
    ("aromatic_systems", "electrons", AromaticSystemForm([2, 0, 1])),
    ("multicenter_bonds", "electrons", MulticenterBondForm([2, 0, 1])),
    ("noncovalent_bonds", "kind", NoncovalentBondForm(NoncovalentBondKind.Ionic)),
    ("stereo_atoms", "configuration", StereoAtomForm(TetrahedralConfiguration.Cw)),
    ("stereo_bonds", "configuration", StereoBondForm.parse("Ct1")),
])
def entity(request):
    return request.param


@pytest.mark.parametrize("operation", ["combine_empty", "combine_nonempty", "transact", "tracked_transact"])
def test_molecule_view_invalidation(molecule, entity, operation):
    name, field, form = entity
    collection = getattr(molecule, name)
    view = collection[0]
    constraints = view.constraints
    iterator = iter(collection)
    first = next(iterator)
    exhausted = iter(collection)
    list(exhausted)
    copied = getattr(view, field)
    saved = view.asdict()
    count = len(collection)

    if operation.startswith("combine"):
        molecule.combine_from(molecule if operation == "combine_nonempty" else Molecule())
    else:
        getattr(molecule, operation)([Edits()])

    for access in (
        lambda: view.id,
        lambda: first.id,
        lambda: repr(view),
        lambda: getattr(view, field),
        lambda: setattr(view, field, getattr(form, field)),
        lambda: view.asdict(),
        lambda: view.constraints,
        lambda: setattr(view, "constraints", constraints),
        lambda: len(constraints),
        lambda: repr(constraints),
        lambda: iter(constraints),
        lambda: constraints.update([]),
        lambda: setattr(form, "constraints", constraints),
        lambda: len(collection),
        lambda: repr(collection),
        lambda: collection[0],
        lambda: collection[-1],
        lambda: collection.__setitem__(0, form),
        lambda: iter(collection),
        lambda: iter(iterator),
        lambda: next(iterator),
        lambda: next(exhausted),
    ):
        with pytest.raises(InvalidatedViewError, match="was invalidated by Molecule mutation"):
            access()

    if not name.startswith("stereo_"):
        with pytest.raises(InvalidatedViewError):
            constraints.asdict()

    fresh = getattr(molecule, name)
    assert len(fresh) == count * (2 if operation == "combine_nonempty" else 1)
    assert fresh[0].id == 0
    assert fresh[0].asdict() == saved
    assert copied == saved[field]
    assert next(iter(fresh)).id == 0


def test_molecule_view_setters(molecule, entity):
    name, field, form = entity
    collection = getattr(molecule, name)
    view = collection[0]
    constraints = view.constraints
    iterator = iter(collection)
    expected = getattr(form, field)

    setattr(view, field, expected)

    assert getattr(view, field) == expected
    assert getattr(getattr(molecule, name)[0], field) == expected
    assert next(iterator).id == 0
    assert list(constraints) == list(view.constraints)

    with pytest.raises(TypeError):
        setattr(view, field, object())
    assert getattr(view, field) == expected

    collection[0] = form
    assert view.asdict() == form.asdict()
    assert list(constraints) == list(view.constraints)
    view.constraints = constraints
    assert view.asdict() == form.asdict()


def test_molecule_combine_from_argument_error(molecule, entity):
    name, field, _ = entity
    view = getattr(molecule, name)[0]
    value = getattr(view, field)

    with pytest.raises(TypeError):
        molecule.combine_from(None)

    assert view.id == 0
    assert getattr(view, field) == value


@pytest.mark.parametrize("name", [
    "atoms", "bonds", "dative_bonds", "aromatic_systems",
    "multicenter_bonds", "noncovalent_bonds", "stereo_atoms", "stereo_bonds",
])
@pytest.mark.parametrize("method", ["combine_from", "transact", "tracked_transact"])
def test_molecule_view_empty_iterator(name, method):
    molecule = Molecule()
    collection = getattr(molecule, name)
    iterator = iter(collection)
    assert next(iterator, None) is None

    getattr(molecule, method)(Molecule() if method == "combine_from" else [])

    with pytest.raises(InvalidatedViewError):
        next(iterator, None)
    with pytest.raises(InvalidatedViewError):
        iter(iterator)
    assert list(getattr(molecule, name)) == []


@pytest.mark.parametrize("name", ["atoms", "bonds", "dative_bonds"])
@pytest.mark.parametrize("method", ["combine_from", "transact", "tracked_transact"])
def test_molecule_view_ring_size_invalidation(molecule, name, method):
    view = getattr(molecule, name)[0]
    constraints = view.constraints
    sizes = constraints.ring_size_count
    sizes[6] = 1
    copied_iterator = iter(sizes)
    copied_items = iter(constraints)
    copied = list(constraints)
    assert sizes[6] == NumForm.Lit(1)
    assert view.constraints.ring_size_count[6] == NumForm.Lit(1)

    getattr(molecule, method)(Molecule() if method == "combine_from" else [])

    for access in (
        lambda: constraints.ring_size_count,
        lambda: sizes[6],
        lambda: sizes.__setitem__(6, 2),
        lambda: sizes.__delitem__(6),
        lambda: 6 in sizes,
        lambda: len(sizes),
        lambda: repr(sizes),
        lambda: iter(sizes),
    ):
        with pytest.raises(InvalidatedViewError):
            access()

    assert list(copied_iterator) == [6]
    assert list(copied_items) == copied
    fresh = getattr(molecule, name)[0].constraints.ring_size_count
    assert fresh[6] == NumForm.Lit(1)
    fresh[6] = 2
    assert fresh[6] == NumForm.Lit(2)
