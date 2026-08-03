import pytest

from umol import (
    AtomAst,
    AtomConstraintAst,
    AtomFieldChange,
    BondAst,
    Constraint,
    ConstraintEdit,
    Edit,
    Entity,
    MoleculeConstraint,
    New,
    ValueAst,
)


def test_new():
    handle = New(3)

    assert handle.index == 3
    assert handle == New(3)
    assert handle != New(4)
    assert repr(handle) == "New(3)"


def test_new_immutability():
    handle = New(3)

    with pytest.raises(AttributeError):
        handle.index = 4
    with pytest.raises(AttributeError):
        handle.extra = 4


@pytest.mark.parametrize("index", [-1, 2**100])
def test_new_error(index):
    with pytest.raises(OverflowError):
        New(index)


def test_constraint_edit():
    constraint = Constraint.Atom(
        0, AtomConstraintAst.Valence(ValueAst.Lit(4))
    )

    identity = ConstraintEdit(constraint)
    created = ConstraintEdit(
        constraint,
        handles={Entity.Atom(0): New(0)},
    )

    assert identity == ConstraintEdit(constraint)
    assert created == ConstraintEdit(
        constraint,
        handles={Entity.Atom(0): New(0)},
    )
    assert created != identity
    assert repr(created) == "ConstraintEdit(...)"


@pytest.mark.parametrize(
    "edit",
    [
        Edit.AddAtoms(atoms=[AtomAst.parse("C")]),
        Edit.AddBonds(bonds=[((0, New(0)), BondAst.parse("1"))]),
        Edit.RemoveTopology(atoms=[New(0)], bonds=[1]),
        Edit.ModifyAtomField(
            id=New(0),
            change=AtomFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(1)
            ),
        ),
        Edit.AddMoleculeConstraint(
            constraint=ConstraintEdit(
                Constraint.Molecule(MoleculeConstraint.Connected(None))
            )
        ),
    ],
)
def test_edit(edit):
    assert edit == edit
    assert repr(edit).startswith(f"Edit.{type(edit).__name__}(")
