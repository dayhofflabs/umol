import pytest

from umol import (
    AtomConstraintForm,
    AtomForm,
    Constraint,
    ConstraintDelta,
    ConstraintEdit,
    ContradictionError,
    Edits,
    Element,
    InvalidatedViewError,
    NumForm,
    Constraints,
    ConstraintsView,
    Molecule,
    MoleculeConstraint,
    RelationalConstraint,
    UnpairedElectronsForm,
)


@pytest.mark.parametrize(
    ("variant", "fields", "expected"),
    [
        (RelationalConstraint.DativeBondDonors, ('bond', 'atoms'), (2, (3, 1, 3))),
        (RelationalConstraint.DativeBondDonor, ('bond', 'atom'), (2, 3)),
        (RelationalConstraint.DativeBondContainsAllDonors, ('bond', 'atoms'), (2, (3, 1, 3))),
        (
            RelationalConstraint.DativeBondAllDonors,
            ('bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (
            RelationalConstraint.DativeBondAnyDonor,
            ('bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (RelationalConstraint.DativeBondAcceptor, ('bond', 'atom'), (2, 3)),
        (
            RelationalConstraint.DativeBondAcceptorSatisfies,
            ('bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (RelationalConstraint.DativeBondParallels, ('dative', 'parallel'), (2, 3)),
        (RelationalConstraint.AromaticSystemAtoms, ('system', 'atoms'), (2, (3, 1, 3))),
        (RelationalConstraint.AromaticSystemContains, ('system', 'atom'), (2, 3)),
        (RelationalConstraint.AromaticSystemContainsAll, ('system', 'atoms'), (2, (3, 1, 3))),
        (
            RelationalConstraint.AromaticSystemAllAtoms,
            ('system', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (
            RelationalConstraint.AromaticSystemAnyAtom,
            ('system', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (RelationalConstraint.MulticenterBondAtoms, ('bond', 'atoms'), (2, (3, 1, 3))),
        (RelationalConstraint.MulticenterBondContains, ('bond', 'atom'), (2, 3)),
        (RelationalConstraint.MulticenterBondContainsAll, ('bond', 'atoms'), (2, (3, 1, 3))),
        (
            RelationalConstraint.MulticenterBondAllAtoms,
            ('bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (
            RelationalConstraint.MulticenterBondAnyAtom,
            ('bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (RelationalConstraint.NoncovalentBondEnds, ('bond', 'atoms'), (2, (3, 1))),
        (RelationalConstraint.NoncovalentBondContains, ('bond', 'atom'), (2, 3)),
        (
            RelationalConstraint.NoncovalentBondEndsSatisfy,
            ('bond', 'predicates'),
            (
                2,
                (
                    AtomConstraintForm.Valence(NumForm.Lit(3)),
                    AtomConstraintForm.Valence(NumForm.Lit(1)),
                ),
            ),
        ),
        (RelationalConstraint.StereoAtomSite, ('stereo_atom', 'atom'), (2, 3)),
        (RelationalConstraint.StereoAtomContains, ('stereo_atom', 'atom'), (2, 3)),
        (RelationalConstraint.StereoAtomLigands, ('stereo_atom', 'atoms'), (2, (3, 1, 3))),
        (
            RelationalConstraint.StereoAtomAllLigands,
            ('stereo_atom', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (
            RelationalConstraint.StereoAtomAnyLigand,
            ('stereo_atom', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (RelationalConstraint.StereoBondSite, ('stereo_bond', 'bond'), (2, 3)),
        (RelationalConstraint.StereoBondContains, ('stereo_bond', 'atom'), (2, 3)),
        (RelationalConstraint.StereoBondLigands, ('stereo_bond', 'atoms'), (2, (3, 1, 3))),
        (
            RelationalConstraint.StereoBondAllLigands,
            ('stereo_bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
        (
            RelationalConstraint.StereoBondAnyLigand,
            ('stereo_bond', 'predicate'),
            (2, AtomConstraintForm.Valence(NumForm.Lit(3))),
        ),
    ],
)
def test_relational_constraint_fields(variant, fields, expected):
    arguments = [list(value) if isinstance(value, tuple) else value for value in expected]
    value = variant(*arguments)
    for argument in arguments:
        if isinstance(argument, list):
            argument.clear()
    assert value.__match_args__ == fields
    assert tuple(getattr(value, field) for field in fields) == expected
    assert variant(**dict(zip(fields, expected))) == value
    assert variant(**{field: getattr(value, field) for field in fields}) == value
    for field, payload in zip(fields, expected):
        if isinstance(payload, tuple):
            assert type(getattr(value, field)) is tuple
        with pytest.raises(AttributeError):
            setattr(value, field, None)
    assert not hasattr(value, "_0")
    expected_fields = ", ".join(
        f"{field}={payload!r}" for field, payload in zip(fields, expected)
    )
    assert repr(value) == f"RelationalConstraint.{variant.__name__}({expected_fields})"


@pytest.mark.parametrize(
    ("variant", "fields", "expected"),
    [
        (MoleculeConstraint.ChargeSum, ("atoms", "sum"), ((3, 1, 3), NumForm.Lit(1))),
        (MoleculeConstraint.ChargeSum, ("atoms", "sum"), (None, NumForm.Lit(0))),
        (MoleculeConstraint.BondOrderSum, ("bonds", "sum"), ((2, 0, 2), NumForm.Lit(3))),
        (MoleculeConstraint.Connected, ("atoms",), ((3, 1, 3),)),
        (MoleculeConstraint.Connected, ("atoms",), (None,)),
        (
            MoleculeConstraint.UnpairedElectronCoupling,
            ("atoms", "unpaired_electrons"),
            ((3, 1, 3), UnpairedElectronsForm(1, 2)),
        ),
    ],
)
def test_molecule_constraint_fields(variant, fields, expected):
    arguments = [list(value) if isinstance(value, tuple) else value for value in expected]
    value = variant(*arguments)
    for argument in arguments:
        if isinstance(argument, list):
            argument.clear()
    assert value.__match_args__ == fields
    assert tuple(getattr(value, field) for field in fields) == expected
    assert variant(**dict(zip(fields, expected))) == value
    assert variant(**{field: getattr(value, field) for field in fields}) == value
    if expected[0] is not None:
        assert type(getattr(value, fields[0])) is tuple
    for field in fields:
        with pytest.raises(AttributeError):
            setattr(value, field, None)
    expected_fields = ", ".join(
        f"{field}={payload!r}" for field, payload in zip(fields, expected)
    )
    assert repr(value) == f"MoleculeConstraint.{variant.__name__}({expected_fields})"


def connected_constraint():
    return Constraint.Molecule(MoleculeConstraint.Connected(None))


def test_moleculeconstraint_unpaired_electron_coupling():
    constraint = MoleculeConstraint.UnpairedElectronCoupling(
        atoms=[2, 3],
        unpaired_electrons=UnpairedElectronsForm(1, 2),
    )

    assert constraint.atoms == (2, 3)
    assert constraint.unpaired_electrons == UnpairedElectronsForm(1, 2)
    assert repr(constraint) == (
        "MoleculeConstraint.UnpairedElectronCoupling(atoms=(2, 3), "
        "unpaired_electrons=UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2)))"
    )


def test_constraints_sequence():
    entry = connected_constraint()
    constraints = Constraints([entry, entry])

    assert len(constraints) == 2
    assert constraints[-1] == entry
    assert list(constraints) == [entry, entry]
    assert repr(constraints) == (
        "Constraints([Constraint.Molecule(MoleculeConstraint.Connected(atoms=None)), "
        "Constraint.Molecule(MoleculeConstraint.Connected(atoms=None))])"
    )


def test_molecule_from_entries_constraints():
    entry = connected_constraint()
    molecule = Molecule.from_entries([], constraints=[entry])

    assert isinstance(molecule.constraints, ConstraintsView)
    assert list(molecule.constraints) == [entry]


@pytest.mark.parametrize("method", ["append", "clear", "update"])
def test_molecule_constraints_readonly(method):
    molecule = Molecule.from_entries([], constraints=[connected_constraint()])
    with pytest.raises(AttributeError):
        getattr(molecule.constraints, method)
    assert list(molecule.constraints) == [connected_constraint()]


def test_molecule_constraints_assignment():
    molecule = Molecule.from_entries([], constraints=[connected_constraint()])
    with pytest.raises(AttributeError):
        molecule.constraints = Constraints([])
    assert list(molecule.constraints) == [connected_constraint()]


def test_constraints_mutation():
    entry = connected_constraint()
    constraints = Constraints([])
    constraints.append(entry)
    constraints.update(Constraints([entry]))
    assert list(constraints) == [entry, entry]
    constraints.clear()
    assert list(constraints) == []


@pytest.mark.parametrize(
    ("entries", "expected"),
    [
        (
            [
                Constraint.Molecule(
                    MoleculeConstraint.ChargeSum([3, 1, 3], NumForm.LitSet({0}))
                ),
                Constraint.Molecule(
                    MoleculeConstraint.ChargeSum([1, 3], NumForm.Lit(0))
                ),
            ],
            [Constraint.Molecule(MoleculeConstraint.ChargeSum([1, 3], NumForm.Lit(0)))],
        ),
    ],
)
def test_constraints_normalize(entries, expected):
    source = Constraints(entries)
    normalized = source.normalize()
    assert normalized == Constraints(expected)
    assert source.normalized_eq(normalized) is True
    assert normalized.normalized_eq(source) is True
    assert source != normalized
    normalized.clear()
    assert normalized == Constraints([])
    assert source == Constraints(entries)


@pytest.mark.parametrize(
    "entries",
    [
        [],
        [
            Constraint.Molecule(
                MoleculeConstraint.ChargeSum(list(range(16)), NumForm.Lit(0))
            )
        ],
    ],
    ids=["empty", "charge-sum"],
)
def test_constraints_normalize_identity(entries):
    source = Constraints(entries)
    equal = Constraints(entries)
    assert source.normalized_eq(equal) is True
    normalized = source.normalize()
    assert normalized == source
    normalized.append(Constraint.Molecule(MoleculeConstraint.Connected(None)))
    assert normalized == Constraints(
        entries + [Constraint.Molecule(MoleculeConstraint.Connected(None))]
    )
    assert source == Constraints(entries)
    assert equal == Constraints(entries)


@pytest.mark.parametrize(
    "other",
    [
        [],
        [
            Constraint.Molecule(
                MoleculeConstraint.ChargeSum(list(range(16)), NumForm.Lit(1))
            )
        ],
    ],
    ids=["empty", "different-sum"],
)
def test_constraints_normalized_eq(other):
    entries = [
        Constraint.Molecule(
            MoleculeConstraint.ChargeSum(list(range(16)), NumForm.Lit(0))
        )
    ]
    source = Constraints(entries)
    different = Constraints(other)
    assert source.normalized_eq(different) is False
    assert different.normalized_eq(source) is False
    assert source == Constraints(entries)
    assert different == Constraints(other)


def test_constraints_normalize_error():
    entries = [
        Constraint.Molecule(
            MoleculeConstraint.ChargeSum([3, 1, 3], NumForm.LitSet(set()))
        )
    ]
    source = Constraints(entries)
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        source.normalize()
    assert source == Constraints(entries)


def test_constraints_view_composition():
    first = Constraint.And([
        Constraint.Atom(0, AtomConstraintForm.Degree(NumForm.Lit(2))),
        Constraint.Not(Constraint.Or([])),
    ])
    second = connected_constraint()
    molecule = Molecule.from_entries([AtomForm(Element("C"))], constraints=[first, second])
    iterator = iter(molecule.constraints)
    entry = next(iterator)

    assert isinstance(entry, Constraint.And)
    assert entry == first
    assert list(entry._0) == list(first._0)
    assert next(iterator) == second
    with pytest.raises(StopIteration):
        next(iterator)
    match entry:
        case Constraint.And(children):
            assert isinstance(children[0], Constraint.Atom)
            assert children[0]._0 == 0
            assert children[0]._1 == AtomConstraintForm.Degree(NumForm.Lit(2))
            assert isinstance(children[1]._0, Constraint.Or)
        case _:
            pytest.fail("constraint variant matching failed")

    with pytest.raises(AttributeError):
        entry._0 = []
    with pytest.raises(TypeError):
        entry._0[0] = second
    with pytest.raises(AttributeError):
        entry._0[0]._0 = 1


def test_constraint_variant_protocol():
    leaf = Constraint.Atom(_0=0, _1=AtomConstraintForm.Degree(NumForm.Lit(2)))
    entry = Constraint.And(_0=[leaf, Constraint.Not(_0=Constraint.Or(_0=[]))])

    assert Constraint.Atom.__qualname__ == "Constraint.Atom"
    assert len(leaf) == 2
    assert leaf[0] == 0
    assert leaf[1] == AtomConstraintForm.Degree(NumForm.Lit(2))
    with pytest.raises(IndexError):
        leaf[2]
    match entry:
        case Constraint.And([Constraint.Atom(0, _), Constraint.Not(Constraint.Or([]))]):
            pass
        case _:
            pytest.fail("nested constraint pattern did not match")


@pytest.mark.parametrize("index", [1, -2])
def test_constraints_view_index_error(index):
    molecule = Molecule.from_entries([], constraints=[connected_constraint()])
    with pytest.raises(IndexError):
        molecule.constraints[index]


@pytest.mark.parametrize("operation", ["combine_from", "transact"])
def test_constraints_view_invalidation(operation):
    first = Constraint.And([Constraint.Not(Constraint.Or([]))])
    molecule = Molecule.from_entries([], constraints=[first])
    collection = molecule.constraints
    iterator = iter(collection)
    entry = next(iterator)
    children = entry._0
    child = children[0]._0
    copied = entry.copy()
    delta = ConstraintDelta.Add(constraint=entry)
    if operation == "transact":
        molecule.transact([])
    else:
        molecule.combine_from(Molecule())

    for access in (
        lambda: len(collection),
        lambda: next(iterator),
        lambda: iter(iterator),
        lambda: entry._0,
        lambda: children[0],
        lambda: repr(child),
        lambda: entry.copy(),
        lambda: ConstraintEdit(entry),
        lambda: delta.inverse(),
        lambda: entry == copied,
    ):
        with pytest.raises(InvalidatedViewError):
            access()

    assert copied == first
    assert list(molecule.constraints) == [first]


@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply"])
def test_constraints_view_edit(method):
    first = Constraint.And([Constraint.Not(Constraint.Or([]))])
    molecule = Molecule.from_entries([], constraints=[first])
    collection = molecule.constraints
    iterator = iter(collection)
    entry = next(iterator)
    children = entry._0
    child = children[0]._0
    delta = ConstraintDelta.Add(constraint=entry)
    if method == "edit":
        product = molecule.edit().finish()
    else:
        result = getattr(molecule, method)(Edits())
        product = result[0] if method == "tracked_apply" else result

    assert product == molecule
    assert list(collection) == [first]
    assert next(iterator, None) is None
    assert entry == entry.copy() == first
    assert list(children) == [Constraint.Not(Constraint.Or([]))]
    assert child == Constraint.Or([])
    assert ConstraintEdit(entry) == ConstraintEdit(first)
    assert delta.inverse() == ConstraintDelta.Remove(constraint=first)


def test_constraints_view_attribute_mutation():
    entry = Constraint.Atom(0, AtomConstraintForm.Degree(NumForm.Lit(2)))
    molecule = Molecule.from_entries([AtomForm(Element("C"))], constraints=[entry])
    view = molecule.constraints[0]
    molecule.atoms[0].charge = 1
    assert view == entry
    assert view._1 == AtomConstraintForm.Degree(NumForm.Lit(2))
