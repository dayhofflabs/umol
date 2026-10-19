import pytest

from umol import (
    AromaticSystemForm,
    AromaticityModel,
    AtomForm,
    BondForm,
    DativeBondForm,
    Edit,
    Edits,
    Element,
    InvalidatedViewError,
    InvalidStructureError,
    TransactionError,
    Molecule,
    MoleculeDefaults,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NumForm,
    StereoAtomForm,
    StereoBondForm,
    StereoLigand,
    StereoLigandKind,
    TetrahedralConfiguration,
    UnderdeterminedError,
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


@pytest.mark.parametrize("operation", [
    "combine_empty", "combine_nonempty", "transact", "tracked_transact",
])
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

    assert copied == saved[field]
    fresh = getattr(molecule, name)
    assert len(fresh) == count * (2 if operation == "combine_nonempty" else 1)
    assert fresh[0].id == 0
    assert fresh[0].asdict() == saved
    assert copied == saved[field]
    assert next(iter(fresh)).id == 0


@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply", "apply_error", "tracked_apply_error"])
def test_molecule_view_edit(molecule, entity, method):
    name, field, form = entity
    original = molecule.copy()
    collection = getattr(molecule, name)
    view = collection[0]
    constraints = view.constraints
    saved_constraints = list(constraints)
    iterator = iter(collection)
    first = next(iterator)
    exhausted = iter(collection)
    list(exhausted)
    result = None
    if method == "edit":
        result = molecule.edit().finish()
    elif method.endswith("_error"):
        edits = Edits([Edit.AddAtoms(atoms=[AtomForm.parse("O")]),
                       Edit.RemoveTopology(atoms=[8], bonds=[])])
        with pytest.raises(TransactionError, match="^atom handle 8 is out of range for 5 entries$"):
            getattr(molecule, method.removesuffix("_error"))(edits)
    else:
        applied = getattr(molecule, method)(Edits())
        result = applied[0] if method == "tracked_apply" else applied

    assert molecule == original
    assert view.id == first.id == 0
    assert view.asdict() == getattr(original, name)[0].asdict()
    assert list(constraints) == saved_constraints
    assert [entry.id for entry in iterator] == list(range(1, len(collection)))
    assert next(exhausted, None) is None
    assert iter(exhausted) is exhausted
    if result is not None:
        assert result == original
        getattr(result, name)[0] = form
        assert molecule == original
        assert getattr(view, field) == getattr(getattr(original, name)[0], field)


@pytest.mark.parametrize("in_place", [False, True])
@pytest.mark.parametrize(
    "method,args",
    [
        ("aromatize", (AromaticityModel.daylight(),)),
        ("kekulize", ([],)),
        ("delocalize_charge", ()),
        ("fold_hydrogens", ()),
        ("unfold_hydrogens", ()),
    ],
)
@pytest.mark.parametrize(
    "source",
    [
        '{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}',
        '{:atoms ["C"]}',
    ],
    ids=["explicit-hydrogen", "no-hydrogen"],
)
def test_molecule_view_transform(source, method, args, in_place):
    molecule = Molecule.parse(source, defaults=MoleculeDefaults.concrete())
    original = molecule.copy()
    alias = molecule
    collection = molecule.atoms
    view = collection[0]
    constraints = view.constraints
    sizes = constraints.ring_size_count
    saved_constraints = list(constraints)
    saved_sizes = list(sizes)
    iterator = iter(collection)
    first = next(iterator)
    exhausted = iter(collection)
    list(exhausted)

    result = getattr(alias, method + ("_into" if in_place else ""))(*args)

    if in_place:
        assert result is None
        for access in (
            lambda: view.id,
            lambda: first.charge,
            lambda: len(constraints),
            lambda: len(sizes),
            lambda: len(collection),
            lambda: next(iterator, None),
            lambda: next(exhausted, None),
        ):
            with pytest.raises(InvalidatedViewError):
                access()
        assert alias.atoms[0].id == 0
    else:
        assert molecule == original
        assert alias == original
        assert view.id == first.id == 0
        assert view.charge == original.atoms[0].charge
        assert list(constraints) == saved_constraints
        assert list(sizes) == saved_sizes
        assert [atom.id for atom in iterator] == list(range(1, len(original.atoms)))
        assert next(exhausted, None) is None
        result.atoms[0].charge = 1
        assert molecule == original
        molecule.atoms[0].charge = -1
        assert view.charge == NumForm.Lit(-1)
        assert result.atoms[0].charge == NumForm.Lit(1)


@pytest.mark.parametrize("in_place", [False, True])
@pytest.mark.parametrize(
    "source,method,args,error",
    [
        pytest.param(
            '{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] '
            ':bonds [[0 1 "*"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]}',
            "aromatize", (AromaticityModel.daylight(),), UnderdeterminedError,
            id="aromatize",
        ),
        pytest.param(
            '{:atoms ["C#a" "C#a" "C#a"] :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]] '
            ':aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]}',
            "kekulize", ([],), ValueError,
            id="kekulize",
        ),
        pytest.param(
            '{:atoms ["C#h-1" "H"] :bonds [[0 1 "1"]]}',
            "fold_hydrogens", (), InvalidStructureError,
            id="fold",
        ),
        pytest.param(
            '{:atoms ["C#h3" "O#h-1"]}',
            "unfold_hydrogens", (), InvalidStructureError,
            id="unfold",
        ),
    ],
)
def test_molecule_view_transform_error(source, method, args, error, in_place):
    molecule = Molecule.parse(source, defaults=MoleculeDefaults.concrete())
    original = molecule.copy()
    alias = molecule
    collection = molecule.atoms
    view = collection[0]
    constraints = view.constraints
    saved_constraints = list(constraints)
    iterator = iter(collection)
    exhausted = iter(collection)
    list(exhausted)

    with pytest.raises(error):
        getattr(alias, method + ("_into" if in_place else ""))(*args)

    assert molecule == original
    assert alias == original
    if in_place:
        for access in (
            lambda: view.charge,
            lambda: len(constraints),
            lambda: len(collection),
            lambda: next(iterator, None),
            lambda: next(exhausted, None),
        ):
            with pytest.raises(InvalidatedViewError):
                access()
    else:
        assert view.charge == original.atoms[0].charge
        assert list(constraints) == saved_constraints
        assert [atom.id for atom in iterator] == list(range(len(original.atoms)))
        assert next(exhausted, None) is None
    assert molecule.atoms[0].asdict() == original.atoms[0].asdict()


@pytest.mark.parametrize(
    "method,args,kwargs,error",
    [
        ("aromatize_into", (None,), {}, TypeError),
        ("aromatize_into", (AromaticityModel.daylight(),), {"config": object()}, TypeError),
        ("kekulize_into", ([None],), {}, TypeError),
        ("kekulize_into", ([-1],), {}, OverflowError),
        ("kekulize_into", ([],), {"config": object()}, TypeError),
        ("delocalize_charge_into", (None,), {}, TypeError),
        ("fold_hydrogens_into", (None,), {}, TypeError),
        ("unfold_hydrogens_into", (None,), {}, TypeError),
    ],
)
def test_molecule_view_transform_argument_error(method, args, kwargs, error):
    molecule = Molecule.parse('{:atoms ["C#h4"]}', defaults=MoleculeDefaults.concrete())
    original = molecule.copy()
    collection = molecule.atoms
    view = collection[0]
    constraints = view.constraints
    saved_constraints = list(constraints)
    iterator = iter(collection)

    with pytest.raises(error):
        getattr(molecule, method)(*args, **kwargs)

    assert molecule == original
    assert view.asdict() == original.atoms[0].asdict()
    assert list(constraints) == saved_constraints
    assert len(collection) == 1
    assert next(iterator).id == 0


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
@pytest.mark.parametrize("method", ["combine_from", "transact", "tracked_transact", "edit", "apply", "tracked_apply"])
def test_molecule_view_empty_iterator(name, method):
    molecule = Molecule()
    collection = getattr(molecule, name)
    iterator = iter(collection)
    assert next(iterator, None) is None

    if method == "edit":
        assert molecule.edit().finish() == molecule
    elif method in ("apply", "tracked_apply"):
        getattr(molecule, method)(Edits())
    else:
        getattr(molecule, method)(Molecule() if method == "combine_from" else [])

    if method in ("edit", "apply", "tracked_apply"):
        assert next(iterator, None) is None
        assert iter(iterator) is iterator
    else:
        with pytest.raises(InvalidatedViewError):
            next(iterator, None)
        with pytest.raises(InvalidatedViewError):
            iter(iterator)
    assert list(getattr(molecule, name)) == []


@pytest.mark.parametrize("name", ["atoms", "bonds", "dative_bonds"])
@pytest.mark.parametrize("method", ["combine_from", "transact", "tracked_transact", "edit", "apply", "tracked_apply"])
def test_molecule_view_ring_size(molecule, name, method):
    view = getattr(molecule, name)[0]
    constraints = view.constraints
    sizes = constraints.ring_size_count
    sizes[6] = 1
    copied_iterator = iter(sizes)
    copied_items = iter(constraints)
    copied = list(constraints)
    assert sizes[6] == NumForm.Lit(1)
    assert view.constraints.ring_size_count[6] == NumForm.Lit(1)

    if method == "edit":
        assert molecule.edit().finish() == molecule
    elif method in ("apply", "tracked_apply"):
        getattr(molecule, method)(Edits())
    else:
        getattr(molecule, method)(Molecule() if method == "combine_from" else [])

    if method in ("edit", "apply", "tracked_apply"):
        assert sizes[6] == NumForm.Lit(1)
        sizes[6] = 2
        assert view.constraints.ring_size_count[6] == NumForm.Lit(2)
        sizes[6] = 1
    else:
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
