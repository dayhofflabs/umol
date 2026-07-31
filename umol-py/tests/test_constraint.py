import pytest

from umol import (
    AromaticValence,
    AromaticValenceAst,
    AtomAst,
    AtomConstraintAst,
    AtomConstraintKey,
    AtomConstraintsAst,
    Element,
    MoleculeAst,
    MulticenterValence,
    MulticenterValenceAst,
    RingMembershipAst,
    RingScope,
    StereoCoset,
    TetrahedralConfiguration,
    TetrahedralStereoAst,
    ValueAst,
)


def test_aromaticvalenceast_aromatic():
    assert AromaticValenceAst.Aromatic(ValueAst.Lit(1)) == AromaticValenceAst.Aromatic(
        ValueAst.Lit(1)
    )


def test_aromaticvalenceast_not_aromatic():
    assert AromaticValenceAst.NotAromatic() == AromaticValenceAst.NotAromatic()


def test_aromaticvalenceast_as_lit():
    assert AromaticValenceAst.NotAromatic().as_lit() == AromaticValence.NotAromatic()
    assert AromaticValenceAst.Aromatic(2).as_lit() == AromaticValence.Aromatic(2)
    assert AromaticValenceAst.Undetermined().as_lit() is None
    assert AromaticValence.NotAromatic().valence_count() == 0
    assert AromaticValence.Aromatic(2).valence_count() == 2


def test_multicentervalenceast_multicenter():
    assert MulticenterValenceAst.Multicenter(ValueAst.Lit(2)) == MulticenterValenceAst.Multicenter(
        ValueAst.Lit(2)
    )


def test_multicentervalenceast_as_lit():
    assert MulticenterValenceAst.NotMulticenter().as_lit() == MulticenterValence.NotMulticenter()
    assert MulticenterValenceAst.Multicenter(3).as_lit() == MulticenterValence.Multicenter(3)
    assert MulticenterValenceAst.Undetermined().as_lit() is None
    assert MulticenterValence.NotMulticenter().valence_count() == 0
    assert MulticenterValence.Multicenter(3).valence_count() == 3


def test_ringscope_size():
    assert RingScope.Size(6) == RingScope.Size(6)


def test_ringmembershipast_fields():
    rm = RingMembershipAst(RingScope.All(), ValueAst.Lit(2))
    assert rm.scope == RingScope.All()
    assert rm.count == ValueAst.Lit(2)


def test_ringmembershipast_int_count():
    assert RingMembershipAst(RingScope.All(), 2).count == ValueAst.Lit(2)


def test_atomconstraint_key_valence():
    assert AtomConstraintAst.Valence(ValueAst.Lit(4)).key == AtomConstraintKey.Valence()


def test_atomconstraint_key_tetrahedral_stereo():
    constraint = AtomConstraintAst.TetrahedralStereo(TetrahedralStereoAst.NotStereo())
    assert constraint.key == AtomConstraintKey.TetrahedralStereo()


def test_atomconstraint_key_ring_membership():
    constraint = AtomConstraintAst.RingMembership(
        RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
    )
    assert constraint.key == AtomConstraintKey.RingMembership(RingScope.Size(6))


def test_atomconstraints_iter():
    constraints = AtomConstraintsAst(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
        ]
    )
    assert len(constraints) == 2
    # mapping-style: iteration and keys() yield keys in canonical order
    assert list(constraints) == list(constraints.keys())
    assert list(constraints) == [AtomConstraintKey.Valence(), AtomConstraintKey.Degree()]
    assert list(constraints.values()) == [
        AtomConstraintAst.Valence(ValueAst.Lit(4)),
        AtomConstraintAst.Degree(ValueAst.Lit(3)),
    ]
    assert list(constraints.items()) == [
        (AtomConstraintKey.Valence(), AtomConstraintAst.Valence(ValueAst.Lit(4))),
        (AtomConstraintKey.Degree(), AtomConstraintAst.Degree(ValueAst.Lit(3))),
    ]


def test_atomconstraints_get():
    constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert AtomConstraintKey.Valence() in constraints
    assert AtomConstraintKey.Degree() not in constraints
    assert constraints.get(AtomConstraintKey.Degree()) is None
    assert constraints.get(AtomConstraintKey.Degree(), 0) == 0
    assert constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_atomconstraints_get_ring_membership():
    constraints = AtomConstraintsAst(
        [AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))]
    )
    assert AtomConstraintKey.RingMembership(RingScope.Size(6)) in constraints
    assert AtomConstraintKey.RingMembership(RingScope.All()) not in constraints
    assert constraints.get(
        AtomConstraintKey.RingMembership(RingScope.Size(6))
    ) == AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))


def test_atomconstraints_valence():
    constraints = AtomConstraintsAst(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
        ]
    )
    assert constraints.valence == ValueAst.Lit(4)
    assert constraints.degree == ValueAst.Lit(3)
    assert constraints.total_valence is None
    assert constraints.aromatic_valence is None


def test_atomconstraints_asdict():
    constraints = AtomConstraintsAst(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
            AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.All(), ValueAst.Lit(2))),
            AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))),
        ]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"valence", "degree", "ring_count", "ring_size_count_6"}
    assert d["valence"] == ValueAst.Lit(4)
    assert d["degree"] == ValueAst.Lit(3)
    assert d["ring_count"] == ValueAst.Lit(2)
    assert d["ring_size_count_6"] == ValueAst.Lit(1)


def test_atomconstraints_ring_size_count():
    constraints = AtomConstraintsAst(
        [AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))]
    )
    assert constraints.ring_size_count[6] == ValueAst.Lit(1)
    assert constraints.ring_size_count[5] is None
    assert constraints.ring_count is None


def test_atomconstraintsast_set():
    constraints = AtomConstraintsAst([])
    constraints.set(AtomConstraintAst.Valence(ValueAst.Lit(4)))
    assert len(constraints) == 1
    assert constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_atomconstraintsast_pop():
    constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert constraints.pop(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )
    assert len(constraints) == 0
    assert constraints.pop(AtomConstraintKey.Valence()) is None


def test_atomconstraintsast_update():
    constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    constraints.update(
        AtomConstraintsAst(
            [
                AtomConstraintAst.Valence(ValueAst.Lit(3)),
                AtomConstraintAst.Degree(ValueAst.Lit(2)),
            ]
        )
    )
    assert len(constraints) == 2
    assert constraints.valence == ValueAst.Lit(3)
    assert constraints.degree == ValueAst.Lit(2)


def test_atomconstraintsview_set():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    mol.atoms[0].constraints.set(
        AtomConstraintAst.AromaticValence(AromaticValenceAst.Aromatic(ValueAst.Lit(1)))
    )
    # a fresh view proves the write hit the molecule, not a transient copy
    constraints = mol.atoms[0].constraints
    assert len(constraints) == 1
    assert constraints.get(
        AtomConstraintKey.AromaticValence()
    ) == AtomConstraintAst.AromaticValence(AromaticValenceAst.Aromatic(ValueAst.Lit(1)))


def test_atomconstraintsview_pop():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_parts([atom])
    assert mol.atoms[0].constraints.pop(
        AtomConstraintKey.Valence()
    ) == AtomConstraintAst.Valence(ValueAst.Lit(4))
    assert len(mol.atoms[0].constraints) == 0


def test_atomconstraintsview_update():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    mol.atoms[0].constraints.update(
        AtomConstraintsAst(
            [
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
                AtomConstraintAst.Degree(ValueAst.Lit(3)),
            ]
        )
    )
    constraints = mol.atoms[0].constraints
    assert len(constraints) == 2
    assert constraints.valence == ValueAst.Lit(4)


def test_atomconstraintsview_atom_backed_set():
    atom = AtomAst(Element("C"))
    atom.constraints.set(AtomConstraintAst.Valence(ValueAst.Lit(4)))
    # a fresh view proves the write mutated the standalone atom in place
    assert len(atom.constraints) == 1
    assert atom.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_atomconstraintsview_atom_backed_pop():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    assert atom.constraints.pop(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )
    assert len(atom.constraints) == 0


def test_atomconstraintsview_atom_backed_update():
    atom = AtomAst(Element("C"))
    atom.constraints.update(
        AtomConstraintsAst(
            [
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
                AtomConstraintAst.Degree(ValueAst.Lit(3)),
            ]
        )
    )
    assert len(atom.constraints) == 2
    assert atom.constraints.valence == ValueAst.Lit(4)


def test_atomconstraintsview_reads():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_parts([atom])
    constraints = mol.atoms[0].constraints
    assert constraints
    assert AtomConstraintKey.Valence() in constraints
    assert constraints.get(AtomConstraintKey.Degree()) is None
    assert constraints.valence == ValueAst.Lit(4)
    assert set(constraints.asdict().keys()) == {"valence"}
    assert list(constraints) == [AtomConstraintKey.Valence()]
    assert list(constraints.values()) == [AtomConstraintAst.Valence(ValueAst.Lit(4))]


def test_atomconstraints_valence_property():
    cs = AtomConstraintsAst([])
    cs.valence = 4
    assert cs.valence.as_lit() == 4


def test_atomconstraints_aromatic_valence_int():
    cs = AtomConstraintsAst([])
    cs.aromatic_valence = 1
    assert cs.aromatic_valence == AromaticValenceAst.Aromatic(ValueAst.Lit(1))


def test_atomconstraints_aromatic_valence_false():
    cs = AtomConstraintsAst([])
    cs.aromatic_valence = False
    assert cs.aromatic_valence == AromaticValenceAst.NotAromatic()


def test_atomconstraints_aromatic_valence_true_error():
    cs = AtomConstraintsAst([])
    with pytest.raises(ValueError):
        cs.aromatic_valence = True


def test_atomconstraints_multicenter_valence_int():
    cs = AtomConstraintsAst([])
    cs.multicenter_valence = 2
    assert cs.multicenter_valence == MulticenterValenceAst.Multicenter(ValueAst.Lit(2))


def test_atomconstraints_tetrahedral_stereo_config():
    cs = AtomConstraintsAst([])
    cs.tetrahedral_stereo = TetrahedralConfiguration.Cw
    assert cs.tetrahedral_stereo == TetrahedralStereoAst.Stereo(StereoCoset.Lit(1))


def test_atomconstraints_tetrahedral_stereo_false():
    cs = AtomConstraintsAst([])
    cs.tetrahedral_stereo = False
    assert cs.tetrahedral_stereo == TetrahedralStereoAst.NotStereo()


def test_atomconstraints_ring_count_property():
    cs = AtomConstraintsAst([])
    cs.ring_count = 2
    assert cs.ring_count.as_lit() == 2


def test_atomconstraints_ring_size_count_subscript():
    cs = AtomConstraintsAst([])
    cs.ring_size_count[6] = 3
    assert cs.ring_size_count[6].as_lit() == 3
    del cs.ring_size_count[6]
    assert cs.ring_size_count[6] is None


def test_atomconstraintsview_property_on_molecule():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    mol.atoms[0].constraints.aromatic_valence = 1
    # a fresh view proves the write hit the molecule
    assert mol.atoms[0].constraints.aromatic_valence == AromaticValenceAst.Aromatic(ValueAst.Lit(1))


def test_atomconstraintsview_ring_size_count_on_molecule():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    mol.atoms[0].constraints.ring_size_count[6] = 3
    assert mol.atoms[0].constraints.ring_size_count[6].as_lit() == 3
    del mol.atoms[0].constraints.ring_size_count[6]
    assert mol.atoms[0].constraints.ring_size_count[6] is None


def test_aromaticvalenceast_aromatic_int():
    assert AromaticValenceAst.Aromatic(1) == AromaticValenceAst.Aromatic(ValueAst.Lit(1))


def test_multicentervalenceast_multicenter_int():
    assert MulticenterValenceAst.Multicenter(2) == MulticenterValenceAst.Multicenter(ValueAst.Lit(2))


def test_tetrahedralconfiguration_enum():
    assert TetrahedralConfiguration.Ccw == TetrahedralConfiguration.Ccw
    assert TetrahedralConfiguration.Ccw != TetrahedralConfiguration.Cw
    assert len(
        {
            TetrahedralConfiguration.Cw,
            TetrahedralConfiguration.Cw,
            TetrahedralConfiguration.Ccw,
        }
    ) == 2


def test_atomconstraint_eq_hash():
    assert AtomConstraintAst.Valence(ValueAst.Lit(4)) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )
    assert AtomConstraintAst.Valence(ValueAst.Lit(4)) != AtomConstraintAst.Valence(
        ValueAst.Lit(5)
    )
    assert AtomConstraintAst.Valence(ValueAst.Lit(4)) != AtomConstraintAst.Degree(
        ValueAst.Lit(4)
    )
    assert (
        len(
            {
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
            }
        )
        == 1
    )


def test_atomconstraint_repr():
    x = AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))
    env = {
        "AtomConstraintAst": AtomConstraintAst,
        "RingMembershipAst": RingMembershipAst,
        "RingScope": RingScope,
        "ValueAst": ValueAst,
    }
    assert eval(repr(x), env) == x


def test_atomconstraintkey_eq_hash():
    assert AtomConstraintKey.Valence() == AtomConstraintKey.Valence()
    assert AtomConstraintKey.Valence() != AtomConstraintKey.Degree()
    assert AtomConstraintKey.RingMembership(RingScope.Size(6)) == AtomConstraintKey.RingMembership(
        RingScope.Size(6)
    )
    assert AtomConstraintKey.RingMembership(RingScope.Size(6)) != AtomConstraintKey.RingMembership(
        RingScope.Size(5)
    )
    assert len({AtomConstraintKey.Valence(), AtomConstraintKey.Valence()}) == 1


def test_atomconstraintsast_eq_repr():
    a = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    b = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert a == b
    assert a != AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(5))])
    assert repr(a) == "AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])"


def test_atomconstraintsast_unhashable():
    # mutable container: value-equal but unhashable, like AtomAst
    with pytest.raises(TypeError):
        hash(AtomConstraintsAst([]))


def test_aromaticvalenceast_eq():
    assert AromaticValenceAst.Aromatic(1) == AromaticValenceAst.Aromatic(1)
    assert AromaticValenceAst.Aromatic(1) != AromaticValenceAst.NotAromatic()


def test_ringmembershipast_eq_repr():
    a = RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
    assert a == RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
    assert a != RingMembershipAst(RingScope.All(), ValueAst.Lit(1))
    assert repr(a) == "RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))"


def test_atomconstraintsview_repr():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    assert repr(mol.atoms[0].constraints) == "AtomConstraintsView(0 entries)"


def test_atomconstraints_getitem_delitem():
    cs = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert AtomConstraintKey.Valence() in cs
    assert cs[AtomConstraintKey.Valence()] == AtomConstraintAst.Valence(ValueAst.Lit(4))
    with pytest.raises(KeyError):
        cs[AtomConstraintKey.Degree()]
    del cs[AtomConstraintKey.Valence()]
    assert not cs
    with pytest.raises(KeyError):
        del cs[AtomConstraintKey.Valence()]


def test_atomconstraints_update_iterable():
    cs = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    cs.update(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(2)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
        ]
    )
    assert len(cs) == 2
    assert cs.valence == ValueAst.Lit(2)


def test_atomconstraints_update_container():
    cs = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    cs.update(AtomConstraintsAst([AtomConstraintAst.Degree(ValueAst.Lit(3))]))
    assert len(cs) == 2


def test_atomconstraintsview_getitem_delitem():
    mol = MoleculeAst.from_parts(
        [
            AtomAst(
                Element("C"),
                constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
            )
        ]
    )
    cs = mol.atoms[0].constraints
    assert AtomConstraintKey.Valence() in cs
    assert cs[AtomConstraintKey.Valence()] == AtomConstraintAst.Valence(ValueAst.Lit(4))
    with pytest.raises(KeyError):
        cs[AtomConstraintKey.Degree()]
    del mol.atoms[0].constraints[AtomConstraintKey.Valence()]
    assert not mol.atoms[0].constraints
    with pytest.raises(KeyError):
        del mol.atoms[0].constraints[AtomConstraintKey.Valence()]


def test_atomconstraintsview_update_from_view():
    src = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    mol.atoms[0].constraints.update(src.constraints)
    assert AtomConstraintKey.Valence() in mol.atoms[0].constraints


def test_atomast_set_constraints_from_value():
    dst = AtomAst(Element("N"))
    dst.constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert dst.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_atomast_set_constraints_from_view():
    src = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    dst = AtomAst(Element("N"))
    dst.constraints = src.constraints  # RHS is a live view, not a value container
    assert dst.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintAst.Valence(
        ValueAst.Lit(4)
    )


def test_ringsizecounts_len_iter_contains():
    cs = AtomConstraintsAst([])
    cs.ring_size_count[6] = 3
    cs.ring_size_count[5] = 1
    rsc = cs.ring_size_count
    assert len(rsc) == 2
    assert sorted(rsc) == [5, 6]
    assert 6 in rsc
    assert 4 not in rsc
