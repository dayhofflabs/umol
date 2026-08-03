import pytest

from umol import (
    AtomAst,
    BooleanAst,
    DativeBondAst,
    DativeBondConstraintAst,
    DativeBondConstraintKey,
    DativeBondConstraintsAst,
    DativeBondUpdate,
    Element,
    MoleculeAst,
    ParseError,
    RingMembershipAst,
    RingScope,
    ValueAst,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (DativeBondUpdate(), (None, DativeBondConstraintsAst([]))),
        (DativeBondUpdate(order=2), (ValueAst.Lit(2), DativeBondConstraintsAst([]))),
        (
            DativeBondUpdate(
                constraints=DativeBondConstraintsAst(
                    [DativeBondConstraintAst.Aromatic(BooleanAst.Undetermined())]
                )
            ),
            (
                None,
                DativeBondConstraintsAst(
                    [DativeBondConstraintAst.Aromatic(BooleanAst.Undetermined())]
                ),
            ),
        ),
    ],
)
def test_dative_bond_update(update, expected):
    assert (update.order, update.constraints) == expected


def ammonia_borane():
    # borane B (id 0) accepts from ammonia N (id 1); dative bond id 0
    return MoleculeAst.from_parts(
        [AtomAst(Element("B")), AtomAst(Element("N"))],
        dative_bonds=[([1], 0, DativeBondAst(1))],
    )


def test_dativebondast_new():
    bond = DativeBondAst(1)
    assert bond.order == ValueAst.Lit(1)
    assert len(bond.constraints) == 0


def test_dativebondast_constraints_kwarg():
    bond = DativeBondAst(
        1,
        constraints=DativeBondConstraintsAst(
            [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
        ),
    )
    assert len(bond.constraints) == 1
    assert bond.constraints.aromatic == BooleanAst.Lit(True)


def test_dativebondast_order_setter():
    bond = DativeBondAst(1)
    bond.order = 2
    assert bond.order == ValueAst.Lit(2)


@pytest.mark.parametrize("dsl", ["1", "1#a", "1#R(6)"])
def test_dativebondast_parse_roundtrip(dsl):
    bond = DativeBondAst.parse(dsl)
    assert str(bond) == dsl
    assert repr(bond) == f"DativeBondAst.parse('{dsl}')"


def test_dativebondast_parse_error():
    with pytest.raises(ParseError):
        DativeBondAst.parse("x#")


def test_dativebondast_asdict():
    bond = DativeBondAst(
        1,
        constraints=DativeBondConstraintsAst(
            [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
        ),
    )
    d = bond.asdict()
    assert set(d.keys()) == {"order", "constraints"}
    assert d["order"] == ValueAst.Lit(1)
    assert d["constraints"]["aromatic"] == BooleanAst.Lit(True)


def test_dativebondast_set_constraints():
    bond = DativeBondAst(1)
    bond.constraints = DativeBondConstraintsAst(
        [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
    )
    assert bond.constraints.aromatic == BooleanAst.Lit(True)


def test_dativebondconstraints_aromatic_and_ring():
    constraints = DativeBondConstraintsAst(
        [
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.All(), ValueAst.Lit(2))
            ),
        ]
    )
    assert len(constraints) == 2
    assert constraints.aromatic == BooleanAst.Lit(True)
    assert constraints.ring_count == ValueAst.Lit(2)


def test_dativebondconstraints_ring_size_count():
    constraints = DativeBondConstraintsAst(
        [
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
            )
        ]
    )
    assert constraints.ring_size_count[6] == ValueAst.Lit(1)
    assert constraints.ring_size_count[5] is None
    assert constraints.ring_count is None


def test_dativebondconstraints_mapping_ops():
    constraints = DativeBondConstraintsAst([])
    constraints.set(DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)))
    assert len(constraints) == 1
    assert DativeBondConstraintKey.Aromatic() in constraints
    assert constraints[DativeBondConstraintKey.Aromatic()] == (
        DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))
    )
    assert [key for key in constraints] == [DativeBondConstraintKey.Aromatic()]
    del constraints[DativeBondConstraintKey.Aromatic()]
    assert len(constraints) == 0


def test_dativebondconstraints_getitem_missing():
    constraints = DativeBondConstraintsAst([])
    with pytest.raises(KeyError):
        constraints[DativeBondConstraintKey.Aromatic()]


def test_dativebondconstraintkey_ring_membership():
    key = DativeBondConstraintKey.RingMembership(RingScope.Size(6))
    assert key == DativeBondConstraintKey.RingMembership(RingScope.Size(6))
    assert key != DativeBondConstraintKey.RingMembership(RingScope.Size(5))


def test_dativebondconstraints_asdict():
    constraints = DativeBondConstraintsAst(
        [
            DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.All(), ValueAst.Lit(2))
            ),
            DativeBondConstraintAst.RingMembership(
                RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
            ),
        ]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"aromatic", "ring_count", "ring_size_count_6"}
    assert d["ring_count"] == ValueAst.Lit(2)
    assert d["ring_size_count_6"] == ValueAst.Lit(1)


def test_dativebondview_fields():
    view = ammonia_borane().dative_bonds[0]
    assert view.id == 0
    assert view.order == ValueAst.Lit(1)
    assert view.acceptor == 0
    assert view.donors == [1]
    # atom_ids is the donors followed by the acceptor
    assert view.atom_ids == (1, 0)
    assert repr(view) == "DativeBondView(id=0)"


def test_dativebondview_set_order():
    mol = ammonia_borane()
    mol.dative_bonds[0].order = 2
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.dative_bonds[0].order == ValueAst.Lit(2)


def test_dativebondview_asdict():
    view = ammonia_borane().dative_bonds[0]
    d = view.asdict()
    assert set(d.keys()) == {"order", "constraints"}
    assert d["order"] == ValueAst.Lit(1)


def test_dativebondview_constraints_write_through():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.set(
        DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))
    )
    constraints = mol.dative_bonds[0].constraints
    assert len(constraints) == 1
    assert constraints.aromatic == BooleanAst.Lit(True)


def test_dativebondview_constraints_aromatic_property():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.aromatic = True
    assert mol.dative_bonds[0].constraints.aromatic == BooleanAst.Lit(True)


def test_dativebondview_constraints_ring_size_count():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints.ring_size_count[6] = 3
    assert mol.dative_bonds[0].constraints.ring_size_count[6].as_lit() == 3
    del mol.dative_bonds[0].constraints.ring_size_count[6]
    assert mol.dative_bonds[0].constraints.ring_size_count[6] is None


def test_dativebondview_set_constraints():
    mol = ammonia_borane()
    mol.dative_bonds[0].constraints = DativeBondConstraintsAst(
        [DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))]
    )
    assert mol.dative_bonds[0].constraints.aromatic == BooleanAst.Lit(True)


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
    mol.dative_bonds[0] = DativeBondAst(2)
    view = mol.dative_bonds[0]
    # value replaced, participants preserved
    assert view.order == ValueAst.Lit(2)
    assert view.acceptor == 0
    assert view.donors == [1]


def test_dativebondviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        ammonia_borane().dative_bonds[5] = DativeBondAst(2)


def test_dativebondviews_iter():
    orders = [view.order for view in ammonia_borane().dative_bonds]
    assert orders == [ValueAst.Lit(1)]


def test_dativebondviews_of():
    mol = ammonia_borane()
    assert mol.dative_bonds.of([1], 0).id == 0
    # roles swapped: no such dative bond
    assert mol.dative_bonds.of([0], 1) is None


def test_dativebondviews_incident():
    # B(0) accepts from N(1); C(2) isolated
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("B")), AtomAst(Element("N")), AtomAst(Element("C"))],
        dative_bonds=[([1], 0, DativeBondAst(1))],
    )
    assert [view.id for view in mol.dative_bonds.incident(0)] == [0]
    assert [view.id for view in mol.dative_bonds.incident(1)] == [0]
    assert mol.dative_bonds.incident(2) == []


def test_dativebondviews_repr():
    assert repr(ammonia_borane().dative_bonds) == "DativeBondViews(len=1)"
