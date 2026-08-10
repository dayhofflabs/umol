import pytest

from umol import (
    AtomForm,
    BooleanForm,
    DativeBondForm,
    DativeBondConstraintForm,
    DativeBondConstraintKey,
    DativeBondConstraintsForm,
    DativeBondUpdate,
    Element,
    Molecule,
    ParseError,
    RingMembershipForm,
    RingScope,
    NumForm,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (DativeBondUpdate(), (None, DativeBondConstraintsForm([]))),
        (DativeBondUpdate(order=2), (NumForm.Lit(2), DativeBondConstraintsForm([]))),
        (
            DativeBondUpdate(
                constraints=DativeBondConstraintsForm(
                    [DativeBondConstraintForm.Aromatic(BooleanForm.Undetermined())]
                )
            ),
            (
                None,
                DativeBondConstraintsForm(
                    [DativeBondConstraintForm.Aromatic(BooleanForm.Undetermined())]
                ),
            ),
        ),
    ],
)
def test_dative_bond_update(update, expected):
    assert (update.order, update.constraints) == expected


@pytest.mark.parametrize(
    ("dsl", "canonical"),
    [
        ("", ""),
        ("2", "2"),
        ("*", "*"),
        ("#R(6)*", "#R(6)*"),
    ],
)
def test_dative_bond_update_parse(dsl, canonical):
    update = DativeBondUpdate.parse(dsl)
    assert str(update) == canonical
    assert repr(update) == f"DativeBondUpdate.parse('{canonical}')"
    assert eval(repr(update)) == update


def test_dative_bond_update_parse_error():
    with pytest.raises(ParseError):
        DativeBondUpdate.parse("#a#a")


def ammonia_borane():
    # borane B (id 0) accepts from ammonia N (id 1); dative bond id 0
    return Molecule.from_entries(
        [AtomForm(Element("B")), AtomForm(Element("N"))],
        dative_bonds=[([1], 0, DativeBondForm(1))],
    )


def test_dativebond_form_new():
    bond = DativeBondForm(1)
    assert bond.order == NumForm.Lit(1)
    assert len(bond.constraints) == 0


def test_dativebond_form_constraints_kwarg():
    bond = DativeBondForm(
        1,
        constraints=DativeBondConstraintsForm(
            [DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))]
        ),
    )
    assert len(bond.constraints) == 1
    assert bond.constraints.aromatic == BooleanForm.Lit(True)


def test_dativebond_form_order_setter():
    bond = DativeBondForm(1)
    bond.order = 2
    assert bond.order == NumForm.Lit(2)


@pytest.mark.parametrize("dsl", ["1", "1#a", "1#R(6)"])
def test_dativebond_form_parse_roundtrip(dsl):
    bond = DativeBondForm.parse(dsl)
    assert str(bond) == dsl
    assert repr(bond) == f"DativeBondForm.parse('{dsl}')"


def test_dativebond_form_parse_error():
    with pytest.raises(ParseError):
        DativeBondForm.parse("x#")


def test_dativebond_form_asdict():
    bond = DativeBondForm(
        1,
        constraints=DativeBondConstraintsForm(
            [DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))]
        ),
    )
    d = bond.asdict()
    assert set(d.keys()) == {"order", "constraints"}
    assert d["order"] == NumForm.Lit(1)
    assert d["constraints"]["aromatic"] == BooleanForm.Lit(True)


def test_dativebond_form_set_constraints():
    bond = DativeBondForm(1)
    bond.constraints = DativeBondConstraintsForm(
        [DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))]
    )
    assert bond.constraints.aromatic == BooleanForm.Lit(True)


def test_dativebondconstraints_aromatic_and_ring():
    constraints = DativeBondConstraintsForm(
        [
            DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            DativeBondConstraintForm.RingMembership(
                RingMembershipForm(RingScope.All(), NumForm.Lit(2))
            ),
        ]
    )
    assert len(constraints) == 2
    assert constraints.aromatic == BooleanForm.Lit(True)
    assert constraints.ring_count == NumForm.Lit(2)


def test_dativebondconstraints_ring_size_count():
    constraints = DativeBondConstraintsForm(
        [
            DativeBondConstraintForm.RingMembership(
                RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
            )
        ]
    )
    assert constraints.ring_size_count[6] == NumForm.Lit(1)
    assert constraints.ring_size_count[5] is None
    assert constraints.ring_count is None


def test_dativebondconstraints_mapping_ops():
    constraints = DativeBondConstraintsForm([])
    constraints.set(DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)))
    assert len(constraints) == 1
    assert DativeBondConstraintKey.Aromatic() in constraints
    assert constraints[DativeBondConstraintKey.Aromatic()] == (
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))
    )
    assert [key for key in constraints] == [DativeBondConstraintKey.Aromatic()]
    del constraints[DativeBondConstraintKey.Aromatic()]
    assert len(constraints) == 0


def test_dativebondconstraints_getitem_missing():
    constraints = DativeBondConstraintsForm([])
    with pytest.raises(KeyError):
        constraints[DativeBondConstraintKey.Aromatic()]


def test_dativebondconstraintkey_ring_membership():
    key = DativeBondConstraintKey.RingMembership(RingScope.Size(6))
    assert key == DativeBondConstraintKey.RingMembership(RingScope.Size(6))
    assert key != DativeBondConstraintKey.RingMembership(RingScope.Size(5))


def test_dativebondconstraints_asdict():
    constraints = DativeBondConstraintsForm(
        [
            DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            DativeBondConstraintForm.RingMembership(
                RingMembershipForm(RingScope.All(), NumForm.Lit(2))
            ),
            DativeBondConstraintForm.RingMembership(
                RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
            ),
        ]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"aromatic", "ring_count", "ring_size_count_6"}
    assert d["ring_count"] == NumForm.Lit(2)
    assert d["ring_size_count_6"] == NumForm.Lit(1)


def test_dativebondview_fields():
    view = ammonia_borane().dative_bonds[0]
    assert view.id == 0
    assert view.order == NumForm.Lit(1)
    assert view.acceptor == 0
    assert view.donors == [1]
    # atom_ids is the donors followed by the acceptor
    assert view.atom_ids == (1, 0)
    assert repr(view) == "DativeBondView(id=0)"


def test_dativebondview_set_order():
    mol = ammonia_borane()
    mol.dative_bonds[0].order = 2
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.dative_bonds[0].order == NumForm.Lit(2)


def test_dativebondview_asdict():
    view = ammonia_borane().dative_bonds[0]
    d = view.asdict()
    assert set(d.keys()) == {"order", "constraints"}
    assert d["order"] == NumForm.Lit(1)


def test_dativebondview_constraints_write_through():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.set(
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))
    )
    constraints = mol.dative_bonds[0].constraints
    assert len(constraints) == 1
    assert constraints.aromatic == BooleanForm.Lit(True)


def test_dativebondview_constraints_aromatic_property():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.aromatic = True
    assert mol.dative_bonds[0].constraints.aromatic == BooleanForm.Lit(True)


def test_dativebondview_constraints_ring_size_count():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.ring_size_count[6] = 3
    assert mol.dative_bonds[0].constraints.ring_size_count[6].as_lit() == 3
    del mol.dative_bonds[0].constraints.ring_size_count[6]
    assert mol.dative_bonds[0].constraints.ring_size_count[6] is None


def test_dativebondview_set_constraints():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints = DativeBondConstraintsForm(
        [DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))]
    )
    assert mol.dative_bonds[0].constraints.aromatic == BooleanForm.Lit(True)


def test_dativebondviews_len_getitem():
    dative = ammonia_borane().dative_bonds
    assert len(dative) == 1
    assert dative[0].id == 0
    assert dative[-1].id == 0
    with pytest.raises(IndexError):
        dative[5]
    with pytest.raises(IndexError):
        dative[-2]


def test_dativebondviews_setitem():
    mol = ammonia_borane()
    mol.dative_bonds[0] = DativeBondForm(2)
    view = mol.dative_bonds[0]
    # value replaced, participants preserved
    assert view.order == NumForm.Lit(2)
    assert view.acceptor == 0
    assert view.donors == [1]


def test_dativebondviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        ammonia_borane().dative_bonds[5] = DativeBondForm(2)


def test_dativebondviews_iter():
    orders = [view.order for view in ammonia_borane().dative_bonds]
    assert orders == [NumForm.Lit(1)]


def test_dativebondviews_of():
    mol = ammonia_borane()
    assert mol.dative_bonds.of([1], 0).id == 0
    # roles swapped: no such dative bond
    assert mol.dative_bonds.of([0], 1) is None


def test_dativebondviews_incident():
    # B(0) accepts from N(1); C(2) isolated
    mol = Molecule.from_entries(
        [AtomForm(Element("B")), AtomForm(Element("N")), AtomForm(Element("C"))],
        dative_bonds=[([1], 0, DativeBondForm(1))],
    )
    assert [view.id for view in mol.dative_bonds.incident(0)] == [0]
    assert [view.id for view in mol.dative_bonds.incident(1)] == [0]
    assert mol.dative_bonds.incident(2) == []


def test_dativebondviews_repr():
    assert repr(ammonia_borane().dative_bonds) == "DativeBondViews(len=1)"
