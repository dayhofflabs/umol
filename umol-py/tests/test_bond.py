import pytest

from umol import (
    AtomForm,
    BondForm,
    BondConstraintForm,
    BondConstraintKey,
    BondConstraintsForm,
    BondUpdate,
    BooleanForm,
    CisTransConfiguration,
    CisTransStereoForm,
    Element,
    Molecule,
    ParseError,
    RingMembershipForm,
    RingScope,
    StereoCoset,
    UnpairedElectronsForm,
    UnpairedElectronsUpdate,
    NumForm,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            BondUpdate(),
            (None, None, UnpairedElectronsUpdate(), BondConstraintsForm([])),
        ),
        (
            BondUpdate(order=2, charge=-1, unpaired_electrons=UnpairedElectronsUpdate(count=2)),
            (
                NumForm.Lit(2),
                NumForm.Lit(-1),
                UnpairedElectronsUpdate(count=2),
                BondConstraintsForm([]),
            ),
        ),
        (
            BondUpdate(
                unpaired_electrons=UnpairedElectronsUpdate(multiplicity=1),
                constraints=BondConstraintsForm(
                    [BondConstraintForm.Aromatic(BooleanForm.Undetermined())]
                ),
            ),
            (
                None,
                None,
                UnpairedElectronsUpdate(multiplicity=1),
                BondConstraintsForm(
                    [BondConstraintForm.Aromatic(BooleanForm.Undetermined())]
                ),
            ),
        ),
    ],
)
def test_bond_update(update, expected):
    assert (
        update.order,
        update.charge,
        update.unpaired_electrons,
        update.constraints,
    ) == expected


@pytest.mark.parametrize(
    ("dsl", "canonical"),
    [
        ("", ""),
        ("#c-1", "#c-"),
        ("*#c*#u*#s*", "*#c*#u*#s*"),
        ("#a*", "#a*"),
    ],
)
def test_bond_update_parse(dsl, canonical):
    update = BondUpdate.parse(dsl)
    assert str(update) == canonical
    assert repr(update) == f"BondUpdate.parse('{canonical}')"
    assert eval(repr(update)) == update


def test_bond_update_parse_error():
    with pytest.raises(ParseError):
        BondUpdate.parse("#c+#c-")


def ethene():
    # two carbons joined by one double bond (bond id 0, atoms 0-1)
    return Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("C"))],
        bonds=[(0, 1, BondForm(2))],
    )


def test_bondast_new():
    bond = BondForm(2)
    assert bond.order == NumForm.Lit(2)
    assert bond.charge == NumForm.Undetermined()


@pytest.mark.parametrize(
    ("actual", "expected"),
    [
        pytest.param(BondForm.single(), BondForm(1), id="single"),
        pytest.param(BondForm.double(), BondForm(2), id="double"),
        pytest.param(BondForm.triple(), BondForm(3), id="triple"),
        pytest.param(BondForm.quadruple(), BondForm(4), id="quadruple"),
        pytest.param(
            BondForm.aromatic(),
            BondForm(
                1,
                constraints=BondConstraintsForm(
                    [BondConstraintForm.Aromatic(BooleanForm.Lit(True))]
                ),
            ),
            id="aromatic",
        ),
    ],
)
def test_bondast_keyword_constructors(actual, expected):
    assert actual == expected


def test_bondast_new_kwargs():
    bond = BondForm(
        1,
        charge=NumForm.Lit(-1),
        unpaired_electrons=UnpairedElectronsForm(0, 1),
    )
    assert bond.order == NumForm.Lit(1)
    assert bond.charge == NumForm.Lit(-1)
    assert bond.unpaired_electrons == UnpairedElectronsForm(0, 1)


def test_bondast_constraints_kwarg():
    bond = BondForm(
        1,
        constraints=BondConstraintsForm(
            [BondConstraintForm.Aromatic(BooleanForm.Lit(True))]
        ),
    )
    assert len(bond.constraints) == 1
    assert bond.constraints.aromatic == BooleanForm.Lit(True)


def test_bondast_order_setter():
    bond = BondForm(1)
    bond.order = 2
    assert bond.order == NumForm.Lit(2)


def test_bondast_charge_setter():
    bond = BondForm(1)
    bond.charge = -1
    assert bond.charge == NumForm.Lit(-1)


def test_bondast_unpaired_electrons_setter():
    bond = BondForm(1)
    bond.unpaired_electrons = UnpairedElectronsForm(0, 1)
    assert bond.unpaired_electrons == UnpairedElectronsForm(0, 1)


def test_bondast_asdict():
    d = BondForm(2, charge=NumForm.Lit(-1)).asdict()
    assert set(d.keys()) == {
        "order",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["order"] == NumForm.Lit(2)
    assert d["charge"] == NumForm.Lit(-1)


def test_bondast_asdict_constraints():
    bond = BondForm(
        1,
        constraints=BondConstraintsForm(
            [BondConstraintForm.Aromatic(BooleanForm.Lit(True))]
        ),
    )
    constraints = bond.asdict()["constraints"]
    assert isinstance(constraints, dict)
    assert set(constraints.keys()) == {"aromatic"}
    assert constraints["aromatic"] == BooleanForm.Lit(True)


def test_bondast_eq():
    assert BondForm(1) == BondForm(1)
    assert BondForm(1) != BondForm(2)


@pytest.mark.parametrize("dsl", ["1", "2#c-", "1#a", "1#R(6)"])
def test_bondast_parse(dsl):
    bond = BondForm.parse(dsl)
    assert str(bond) == dsl
    assert repr(bond) == f"BondForm.parse('{dsl}')"


def test_bondast_parse_error():
    with pytest.raises(ParseError):
        BondForm.parse("x#")


def test_bondconstraint_key_aromatic():
    assert BondConstraintForm.Aromatic(BooleanForm.Lit(True)).key == BondConstraintKey.Aromatic()


def test_bondconstraint_key_cis_trans_stereo():
    constraint = BondConstraintForm.CisTransStereo(CisTransStereoForm.NotStereo())
    assert constraint.key == BondConstraintKey.CisTransStereo()


def test_bondconstraint_key_ring_membership():
    constraint = BondConstraintForm.RingMembership(
        RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
    )
    assert constraint.key == BondConstraintKey.RingMembership(RingScope.Size(6))


def test_bondconstraints_iter():
    constraints = BondConstraintsForm(
        [
            BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(2))),
        ]
    )
    assert len(constraints) == 2
    assert list(constraints) == list(constraints.keys())
    assert list(constraints) == [
        BondConstraintKey.Aromatic(),
        BondConstraintKey.RingMembership(RingScope.All()),
    ]
    assert list(constraints.values()) == [
        BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(2))),
    ]
    assert list(constraints.items()) == [
        (BondConstraintKey.Aromatic(), BondConstraintForm.Aromatic(BooleanForm.Lit(True))),
        (
            BondConstraintKey.RingMembership(RingScope.All()),
            BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(2))),
        ),
    ]


def test_bondconstraints_get():
    constraints = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert BondConstraintKey.Aromatic() in constraints
    assert BondConstraintKey.CisTransStereo() not in constraints
    assert constraints.get(BondConstraintKey.CisTransStereo()) is None
    assert constraints.get(BondConstraintKey.CisTransStereo(), 0) == 0
    assert constraints.get(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )


def test_bondconstraints_aromatic():
    empty = BondConstraintsForm([])
    # aromatic is non-optional: unset reads back as Undetermined
    assert empty.aromatic == BooleanForm.Undetermined()
    assert empty.cis_trans_stereo is None
    assert empty.ring_count is None
    constraints = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert constraints.aromatic == BooleanForm.Lit(True)


def test_bondconstraints_asdict():
    constraints = BondConstraintsForm(
        [
            BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(2))),
            BondConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))),
        ]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"aromatic", "ring_count", "ring_size_count_6"}
    assert d["aromatic"] == BooleanForm.Lit(True)
    assert d["ring_count"] == NumForm.Lit(2)
    assert d["ring_size_count_6"] == NumForm.Lit(1)


def test_bondconstraints_ring_size_count():
    constraints = BondConstraintsForm(
        [BondConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))]
    )
    assert constraints.ring_size_count[6] == NumForm.Lit(1)
    assert constraints.ring_size_count[5] is None
    assert constraints.ring_count is None


def test_bondconstraintsast_set():
    constraints = BondConstraintsForm([])
    constraints.set(BondConstraintForm.Aromatic(BooleanForm.Lit(True)))
    assert len(constraints) == 1
    assert constraints.get(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )


def test_bondconstraintsast_pop():
    constraints = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert constraints.pop(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )
    assert len(constraints) == 0
    assert constraints.pop(BondConstraintKey.Aromatic()) is None


def test_bondconstraintsast_update():
    constraints = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    constraints.update(
        BondConstraintsForm(
            [
                BondConstraintForm.Aromatic(BooleanForm.Lit(False)),
                BondConstraintForm.RingMembership(
                    RingMembershipForm(RingScope.All(), NumForm.Lit(2))
                ),
            ]
        )
    )
    assert len(constraints) == 2
    assert constraints.aromatic == BooleanForm.Lit(False)
    assert constraints.ring_count == NumForm.Lit(2)


def test_bondconstraints_update_iterable():
    cs = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    cs.update(
        [
            BondConstraintForm.Aromatic(BooleanForm.Lit(False)),
            BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(1))),
        ]
    )
    assert len(cs) == 2
    assert cs.aromatic == BooleanForm.Lit(False)


def test_bondconstraints_aromatic_property():
    cs = BondConstraintsForm([])
    cs.aromatic = True
    assert cs.aromatic == BooleanForm.Lit(True)


def test_bondconstraints_cis_trans_stereo_config():
    cs = BondConstraintsForm([])
    cs.cis_trans_stereo = CisTransConfiguration.E
    assert cs.cis_trans_stereo == CisTransStereoForm.Stereo(StereoCoset.Lit(1))


def test_bondconstraints_cis_trans_stereo_false():
    cs = BondConstraintsForm([])
    cs.cis_trans_stereo = False
    assert cs.cis_trans_stereo == CisTransStereoForm.NotStereo()


def test_bondconstraints_cis_trans_stereo_true_error():
    cs = BondConstraintsForm([])
    with pytest.raises(ValueError):
        cs.cis_trans_stereo = True


def test_bondconstraints_ring_count_property():
    cs = BondConstraintsForm([])
    cs.ring_count = 2
    assert cs.ring_count.as_lit() == 2


def test_bondconstraints_ring_size_count_subscript():
    cs = BondConstraintsForm([])
    cs.ring_size_count[6] = 3
    assert cs.ring_size_count[6].as_lit() == 3
    del cs.ring_size_count[6]
    assert cs.ring_size_count[6] is None


def test_bondringsizecounts_len_iter_contains():
    cs = BondConstraintsForm([])
    cs.ring_size_count[6] = 3
    cs.ring_size_count[5] = 1
    rsc = cs.ring_size_count
    assert len(rsc) == 2
    assert sorted(rsc) == [5, 6]
    assert 6 in rsc
    assert 4 not in rsc


def test_bondconstraints_getitem_delitem():
    cs = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert BondConstraintKey.Aromatic() in cs
    assert cs[BondConstraintKey.Aromatic()] == BondConstraintForm.Aromatic(BooleanForm.Lit(True))
    with pytest.raises(KeyError):
        cs[BondConstraintKey.CisTransStereo()]
    del cs[BondConstraintKey.Aromatic()]
    assert not cs
    with pytest.raises(KeyError):
        del cs[BondConstraintKey.Aromatic()]


def test_bondconstraintsview_set():
    bond = BondForm(1)
    bond.constraints.set(BondConstraintForm.Aromatic(BooleanForm.Lit(True)))
    # a fresh view proves the write mutated the standalone bond in place
    assert len(bond.constraints) == 1
    assert bond.constraints.get(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )


def test_bondconstraintsview_pop():
    bond = BondForm(
        1,
        constraints=BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))]),
    )
    assert bond.constraints.pop(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )
    assert len(bond.constraints) == 0


def test_bondconstraintsview_update():
    bond = BondForm(1)
    bond.constraints.update(
        BondConstraintsForm(
            [
                BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
                BondConstraintForm.RingMembership(
                    RingMembershipForm(RingScope.All(), NumForm.Lit(2))
                ),
            ]
        )
    )
    assert len(bond.constraints) == 2
    assert bond.constraints.aromatic == BooleanForm.Lit(True)


def test_bondconstraintsview_aromatic_property():
    bond = BondForm(1)
    bond.constraints.aromatic = True
    # a fresh view proves the write hit the bond
    assert bond.constraints.aromatic == BooleanForm.Lit(True)


def test_bondconstraintsview_ring_size_count():
    bond = BondForm(1)
    bond.constraints.ring_size_count[6] = 3
    assert bond.constraints.ring_size_count[6].as_lit() == 3
    del bond.constraints.ring_size_count[6]
    assert bond.constraints.ring_size_count[6] is None


def test_bondconstraintsview_reads():
    bond = BondForm(
        1,
        constraints=BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))]),
    )
    constraints = bond.constraints
    assert constraints
    assert BondConstraintKey.Aromatic() in constraints
    assert constraints.get(BondConstraintKey.CisTransStereo()) is None
    assert constraints.aromatic == BooleanForm.Lit(True)
    assert set(constraints.asdict().keys()) == {"aromatic"}
    assert list(constraints) == [BondConstraintKey.Aromatic()]
    assert list(constraints.values()) == [BondConstraintForm.Aromatic(BooleanForm.Lit(True))]


def test_bondconstraintsview_getitem_delitem():
    bond = BondForm(
        1,
        constraints=BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))]),
    )
    cs = bond.constraints
    assert BondConstraintKey.Aromatic() in cs
    assert cs[BondConstraintKey.Aromatic()] == BondConstraintForm.Aromatic(BooleanForm.Lit(True))
    with pytest.raises(KeyError):
        cs[BondConstraintKey.CisTransStereo()]
    del bond.constraints[BondConstraintKey.Aromatic()]
    assert not bond.constraints
    with pytest.raises(KeyError):
        del bond.constraints[BondConstraintKey.Aromatic()]


def test_bondconstraintsview_update_from_view():
    src = BondForm(
        1,
        constraints=BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))]),
    )
    dst = BondForm(2)
    dst.constraints.update(src.constraints)
    assert BondConstraintKey.Aromatic() in dst.constraints


def test_bondast_set_constraints_from_value():
    dst = BondForm(2)
    dst.constraints = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert dst.constraints.get(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )


def test_bondast_set_constraints_from_view():
    src = BondForm(
        1,
        constraints=BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))]),
    )
    dst = BondForm(2)
    dst.constraints = src.constraints  # RHS is a live view, not a value container
    assert dst.constraints.get(BondConstraintKey.Aromatic()) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )


def test_bondconstraint_eq_hash():
    assert BondConstraintForm.Aromatic(BooleanForm.Lit(True)) == BondConstraintForm.Aromatic(
        BooleanForm.Lit(True)
    )
    assert BondConstraintForm.Aromatic(BooleanForm.Lit(True)) != BondConstraintForm.Aromatic(
        BooleanForm.Lit(False)
    )
    assert (
        len(
            {
                BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
                BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            }
        )
        == 1
    )


def test_bondconstraint_repr():
    x = BondConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))
    env = {
        "BondConstraintForm": BondConstraintForm,
        "RingMembershipForm": RingMembershipForm,
        "RingScope": RingScope,
        "NumForm": NumForm,
    }
    assert eval(repr(x), env) == x


def test_bondconstraintkey_eq_hash():
    assert BondConstraintKey.Aromatic() == BondConstraintKey.Aromatic()
    assert BondConstraintKey.Aromatic() != BondConstraintKey.CisTransStereo()
    assert BondConstraintKey.RingMembership(RingScope.Size(6)) == BondConstraintKey.RingMembership(
        RingScope.Size(6)
    )
    assert BondConstraintKey.RingMembership(RingScope.Size(6)) != BondConstraintKey.RingMembership(
        RingScope.Size(5)
    )
    assert len({BondConstraintKey.Aromatic(), BondConstraintKey.Aromatic()}) == 1


def test_bondconstraintsast_eq_repr():
    a = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    b = BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])
    assert a == b
    assert a != BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(False))])
    assert repr(a) == "BondConstraintsForm([BondConstraintForm.Aromatic(BooleanForm.Lit(True))])"


def test_bondconstraintsast_unhashable():
    # mutable container: value-equal but unhashable, like BondForm
    with pytest.raises(TypeError):
        hash(BondConstraintsForm([]))


def test_bondconstraintsview_repr():
    bond = BondForm(1)
    assert repr(bond.constraints) == "BondConstraintsView(0 entries)"


def test_bondview_fields():
    view = ethene().bonds[0]
    assert view.id == 0
    assert view.order == NumForm.Lit(2)
    assert view.atom_ids == (0, 1)
    assert repr(view) == "BondView(id=0)"


def test_bondview_set_order():
    mol = ethene()
    mol.bonds[0].order = 1
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.bonds[0].order == NumForm.Lit(1)


def test_bondview_set_charge():
    mol = ethene()
    mol.bonds[0].charge = -1
    assert mol.bonds[0].charge == NumForm.Lit(-1)


def test_bondview_set_unpaired_electrons():
    mol = ethene()
    mol.bonds[0].unpaired_electrons = UnpairedElectronsForm(0, 1)
    assert mol.bonds[0].unpaired_electrons == UnpairedElectronsForm(0, 1)


def test_bondview_asdict():
    view = ethene().bonds[0]
    d = view.asdict()
    assert set(d.keys()) == {
        "order",
        "charge",
        "unpaired_electrons",
        "constraints",
    }
    assert d["order"] == NumForm.Lit(2)


def test_bondview_constraints_write_through():
    mol = ethene()
    mol.bonds[0].constraints.set(BondConstraintForm.Aromatic(BooleanForm.Lit(True)))
    # a fresh view proves the write hit the molecule, not a transient copy
    constraints = mol.bonds[0].constraints
    assert len(constraints) == 1
    assert constraints.aromatic == BooleanForm.Lit(True)


def test_bondview_constraints_aromatic_property():
    mol = ethene()
    mol.bonds[0].constraints.aromatic = True
    assert mol.bonds[0].constraints.aromatic == BooleanForm.Lit(True)


def test_bondview_constraints_ring_size_count():
    mol = ethene()
    mol.bonds[0].constraints.ring_size_count[6] = 3
    assert mol.bonds[0].constraints.ring_size_count[6].as_lit() == 3
    del mol.bonds[0].constraints.ring_size_count[6]
    assert mol.bonds[0].constraints.ring_size_count[6] is None


def test_bondview_set_constraints():
    mol = ethene()
    mol.bonds[0].constraints = BondConstraintsForm(
        [BondConstraintForm.Aromatic(BooleanForm.Lit(True))]
    )
    assert mol.bonds[0].constraints.aromatic == BooleanForm.Lit(True)


def test_bondviews_len_getitem():
    bonds = ethene().bonds
    assert len(bonds) == 1
    assert bonds[0].id == 0
    assert bonds[-1].id == 0
    with pytest.raises(IndexError):
        bonds[5]
    with pytest.raises(IndexError):
        bonds[-2]


def test_bondviews_setitem():
    mol = ethene()
    mol.bonds[0] = BondForm(1)
    view = mol.bonds[0]
    # value replaced, endpoints preserved
    assert view.order == NumForm.Lit(1)
    assert view.atom_ids == (0, 1)


def test_bondviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        ethene().bonds[5] = BondForm(1)


def test_bondviews_iter():
    orders = [view.order for view in ethene().bonds]
    assert orders == [NumForm.Lit(2)]


def test_bondviews_of():
    mol = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("C")), AtomForm(Element("C"))],
        bonds=[(0, 1, BondForm(1))],
    )
    assert mol.bonds.of(0, 1).id == 0
    assert mol.bonds.of(1, 0).id == 0
    assert mol.bonds.of(1, 2) is None


def test_bondviews_repr():
    assert repr(ethene().bonds) == "BondViews(len=1)"
