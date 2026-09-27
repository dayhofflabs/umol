import pytest

from umol import (
    AtomConstraintForm,
    AtomForm,
    Constraint,
    ConstraintDelta,
    ConstraintEdit,
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


def connected_constraint():
    return Constraint.Molecule(MoleculeConstraint.Connected(None))


def test_moleculeconstraint_unpaired_electron_coupling():
    constraint = MoleculeConstraint.UnpairedElectronCoupling(
        atoms=[2, 3],
        unpaired_electrons=UnpairedElectronsForm(1, 2),
    )

    assert constraint.atoms == [2, 3]
    assert constraint.unpaired_electrons == UnpairedElectronsForm(1, 2)
    assert repr(constraint) == (
        "MoleculeConstraint.UnpairedElectronCoupling(atoms=[2, 3], "
        "unpaired_electrons=UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2)))"
    )


def test_constraints_sequence():
    entry = connected_constraint()
    constraints = Constraints([entry, entry])

    assert len(constraints) == 2
    assert constraints[-1] == entry
    assert list(constraints) == [entry, entry]
    assert repr(constraints) == (
        "Constraints([Constraint.Molecule(MoleculeConstraint.Connected(None)), "
        "Constraint.Molecule(MoleculeConstraint.Connected(None))])"
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


def test_constraints_view_invalidation():
    first = Constraint.And([Constraint.Not(Constraint.Or([]))])
    molecule = Molecule.from_entries([], constraints=[first])
    collection = molecule.constraints
    iterator = iter(collection)
    entry = next(iterator)
    children = entry._0
    child = children[0]._0
    copied = entry.copy()
    delta = ConstraintDelta.Add(constraint=entry)
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


def test_constraints_view_attribute_mutation():
    entry = Constraint.Atom(0, AtomConstraintForm.Degree(NumForm.Lit(2)))
    molecule = Molecule.from_entries([AtomForm(Element("C"))], constraints=[entry])
    view = molecule.constraints[0]
    molecule.atoms[0].charge = 1
    assert view == entry
    assert view._1 == AtomConstraintForm.Degree(NumForm.Lit(2))
