import pytest

from umol import (
    AromaticValence,
    AromaticValenceForm,
    AtomForm,
    AtomConstraintForm,
    AtomConstraintKey,
    AtomConstraintsForm,
    Element,
    Molecule,
    MulticenterValence,
    MulticenterValenceForm,
    RingMembershipForm,
    RingScope,
    StereoCoset,
    TetrahedralConfiguration,
    TetrahedralStereoForm,
    NumForm,
)


def test_aromatic_valence_form_aromatic():
    assert AromaticValenceForm.Aromatic(NumForm.Lit(1)) == AromaticValenceForm.Aromatic(
        NumForm.Lit(1)
    )


def test_aromatic_valence_form_not_aromatic():
    assert AromaticValenceForm.NotAromatic() == AromaticValenceForm.NotAromatic()


def test_aromatic_valence_form_as_lit():
    assert AromaticValenceForm.NotAromatic().as_lit() == AromaticValence.NotAromatic()
    assert AromaticValenceForm.Aromatic(2).as_lit() == AromaticValence.Aromatic(2)
    assert AromaticValenceForm.Undetermined().as_lit() is None
    assert AromaticValence.NotAromatic().valence_count() == 0
    assert AromaticValence.Aromatic(2).valence_count() == 2


def test_multicenter_valence_form_multicenter():
    assert MulticenterValenceForm.Multicenter(NumForm.Lit(2)) == MulticenterValenceForm.Multicenter(
        NumForm.Lit(2)
    )


def test_multicenter_valence_form_as_lit():
    assert MulticenterValenceForm.NotMulticenter().as_lit() == MulticenterValence.NotMulticenter()
    assert MulticenterValenceForm.Multicenter(3).as_lit() == MulticenterValence.Multicenter(3)
    assert MulticenterValenceForm.Undetermined().as_lit() is None
    assert MulticenterValence.NotMulticenter().valence_count() == 0
    assert MulticenterValence.Multicenter(3).valence_count() == 3


def test_ringscope_size():
    assert RingScope.Size(6) == RingScope.Size(6)


def test_ring_membership_form_fields():
    rm = RingMembershipForm(RingScope.All(), NumForm.Lit(2))
    assert rm.scope == RingScope.All()
    assert rm.count == NumForm.Lit(2)


def test_ring_membership_form_int_count():
    assert RingMembershipForm(RingScope.All(), 2).count == NumForm.Lit(2)


def test_atomconstraint_key_valence():
    assert AtomConstraintForm.Valence(NumForm.Lit(4)).key == AtomConstraintKey.Valence()


def test_atomconstraint_key_tetrahedral_stereo():
    constraint = AtomConstraintForm.TetrahedralStereo(TetrahedralStereoForm.NotStereo())
    assert constraint.key == AtomConstraintKey.TetrahedralStereo()


def test_atomconstraint_key_ring_membership():
    constraint = AtomConstraintForm.RingMembership(
        RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
    )
    assert constraint.key == AtomConstraintKey.RingMembership(RingScope.Size(6))


def test_atomconstraints_iter():
    constraints = AtomConstraintsForm(
        [
            AtomConstraintForm.Valence(NumForm.Lit(4)),
            AtomConstraintForm.Degree(NumForm.Lit(3)),
        ]
    )
    assert len(constraints) == 2
    # mapping-style: iteration and keys() yield keys in canonical order
    assert list(constraints) == list(constraints.keys())
    assert list(constraints) == [AtomConstraintKey.Valence(), AtomConstraintKey.Degree()]
    assert list(constraints.values()) == [
        AtomConstraintForm.Valence(NumForm.Lit(4)),
        AtomConstraintForm.Degree(NumForm.Lit(3)),
    ]
    assert list(constraints.items()) == [
        (AtomConstraintKey.Valence(), AtomConstraintForm.Valence(NumForm.Lit(4))),
        (AtomConstraintKey.Degree(), AtomConstraintForm.Degree(NumForm.Lit(3))),
    ]


def test_atomconstraints_get():
    constraints = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    assert AtomConstraintKey.Valence() in constraints
    assert AtomConstraintKey.Degree() not in constraints
    assert constraints.get(AtomConstraintKey.Degree()) is None
    assert constraints.get(AtomConstraintKey.Degree(), 0) == 0
    assert constraints.get(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )


def test_atomconstraints_get_ring_membership():
    constraints = AtomConstraintsForm(
        [AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))]
    )
    assert AtomConstraintKey.RingMembership(RingScope.Size(6)) in constraints
    assert AtomConstraintKey.RingMembership(RingScope.All()) not in constraints
    assert constraints.get(
        AtomConstraintKey.RingMembership(RingScope.Size(6))
    ) == AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))


def test_atomconstraints_valence():
    constraints = AtomConstraintsForm(
        [
            AtomConstraintForm.Valence(NumForm.Lit(4)),
            AtomConstraintForm.Degree(NumForm.Lit(3)),
        ]
    )
    assert constraints.valence == NumForm.Lit(4)
    assert constraints.degree == NumForm.Lit(3)
    assert constraints.total_valence is None
    assert constraints.aromatic_valence is None


def test_atomconstraints_asdict():
    constraints = AtomConstraintsForm(
        [
            AtomConstraintForm.Valence(NumForm.Lit(4)),
            AtomConstraintForm.Degree(NumForm.Lit(3)),
            AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(2))),
            AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))),
        ]
    )
    d = constraints.asdict()
    assert set(d.keys()) == {"valence", "degree", "ring_count", "ring_size_count_6"}
    assert d["valence"] == NumForm.Lit(4)
    assert d["degree"] == NumForm.Lit(3)
    assert d["ring_count"] == NumForm.Lit(2)
    assert d["ring_size_count_6"] == NumForm.Lit(1)


def test_atomconstraints_ring_size_count():
    constraints = AtomConstraintsForm(
        [AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))]
    )
    assert constraints.ring_size_count[6] == NumForm.Lit(1)
    assert constraints.ring_size_count[5] is None
    assert constraints.ring_count is None


def test_atom_constraints_form_set():
    constraints = AtomConstraintsForm([])
    constraints.set(AtomConstraintForm.Valence(NumForm.Lit(4)))
    assert len(constraints) == 1
    assert constraints.get(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )


def test_atom_constraints_form_pop():
    constraints = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    assert constraints.pop(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )
    assert len(constraints) == 0
    assert constraints.pop(AtomConstraintKey.Valence()) is None


def test_atom_constraints_form_update():
    constraints = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    constraints.update(
        AtomConstraintsForm(
            [
                AtomConstraintForm.Valence(NumForm.Lit(3)),
                AtomConstraintForm.Degree(NumForm.Lit(2)),
            ]
        )
    )
    assert len(constraints) == 2
    assert constraints.valence == NumForm.Lit(3)
    assert constraints.degree == NumForm.Lit(2)


def test_atomconstraintsview_set():
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    mol.atoms[0].constraints.set(
        AtomConstraintForm.AromaticValence(AromaticValenceForm.Aromatic(NumForm.Lit(1)))
    )
    # a fresh view proves the write hit the molecule, not a transient copy
    constraints = mol.atoms[0].constraints
    assert len(constraints) == 1
    assert constraints.get(
        AtomConstraintKey.AromaticValence()
    ) == AtomConstraintForm.AromaticValence(AromaticValenceForm.Aromatic(NumForm.Lit(1)))


def test_atomconstraintsview_pop():
    atom = AtomForm(
        Element("C"),
        constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
    )
    mol = Molecule.from_entries([atom])
    assert mol.atoms[0].constraints.pop(
        AtomConstraintKey.Valence()
    ) == AtomConstraintForm.Valence(NumForm.Lit(4))
    assert len(mol.atoms[0].constraints) == 0


def test_atomconstraintsview_update():
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    mol.atoms[0].constraints.update(
        AtomConstraintsForm(
            [
                AtomConstraintForm.Valence(NumForm.Lit(4)),
                AtomConstraintForm.Degree(NumForm.Lit(3)),
            ]
        )
    )
    constraints = mol.atoms[0].constraints
    assert len(constraints) == 2
    assert constraints.valence == NumForm.Lit(4)


def test_atomconstraintsview_atom_backed_set():
    atom = AtomForm(Element("C"))
    atom.constraints.set(AtomConstraintForm.Valence(NumForm.Lit(4)))
    # a fresh view proves the write mutated the standalone atom in place
    assert len(atom.constraints) == 1
    assert atom.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )


def test_atomconstraintsview_atom_backed_pop():
    atom = AtomForm(
        Element("C"),
        constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
    )
    assert atom.constraints.pop(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )
    assert len(atom.constraints) == 0


def test_atomconstraintsview_atom_backed_update():
    atom = AtomForm(Element("C"))
    atom.constraints.update(
        AtomConstraintsForm(
            [
                AtomConstraintForm.Valence(NumForm.Lit(4)),
                AtomConstraintForm.Degree(NumForm.Lit(3)),
            ]
        )
    )
    assert len(atom.constraints) == 2
    assert atom.constraints.valence == NumForm.Lit(4)


def test_atomconstraintsview_reads():
    atom = AtomForm(
        Element("C"),
        constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
    )
    mol = Molecule.from_entries([atom])
    constraints = mol.atoms[0].constraints
    assert constraints
    assert AtomConstraintKey.Valence() in constraints
    assert constraints.get(AtomConstraintKey.Degree()) is None
    assert constraints.valence == NumForm.Lit(4)
    assert set(constraints.asdict().keys()) == {"valence"}
    assert list(constraints) == [AtomConstraintKey.Valence()]
    assert list(constraints.values()) == [AtomConstraintForm.Valence(NumForm.Lit(4))]


def test_atomconstraints_valence_property():
    cs = AtomConstraintsForm([])
    cs.valence = 4
    assert cs.valence.as_lit() == 4


def test_atomconstraints_aromatic_valence_int():
    cs = AtomConstraintsForm([])
    cs.aromatic_valence = 1
    assert cs.aromatic_valence == AromaticValenceForm.Aromatic(NumForm.Lit(1))


def test_atomconstraints_aromatic_valence_false():
    cs = AtomConstraintsForm([])
    cs.aromatic_valence = False
    assert cs.aromatic_valence == AromaticValenceForm.NotAromatic()


def test_atomconstraints_aromatic_valence_true_error():
    cs = AtomConstraintsForm([])
    with pytest.raises(ValueError):
        cs.aromatic_valence = True


def test_atomconstraints_multicenter_valence_int():
    cs = AtomConstraintsForm([])
    cs.multicenter_valence = 2
    assert cs.multicenter_valence == MulticenterValenceForm.Multicenter(NumForm.Lit(2))


def test_atomconstraints_tetrahedral_stereo_config():
    cs = AtomConstraintsForm([])
    cs.tetrahedral_stereo = TetrahedralConfiguration.Cw
    assert cs.tetrahedral_stereo == TetrahedralStereoForm.Stereo(StereoCoset.Lit(1))


def test_atomconstraints_tetrahedral_stereo_false():
    cs = AtomConstraintsForm([])
    cs.tetrahedral_stereo = False
    assert cs.tetrahedral_stereo == TetrahedralStereoForm.NotStereo()


def test_atomconstraints_ring_count_property():
    cs = AtomConstraintsForm([])
    cs.ring_count = 2
    assert cs.ring_count.as_lit() == 2


def test_atomconstraints_ring_size_count_subscript():
    cs = AtomConstraintsForm([])
    cs.ring_size_count[6] = 3
    assert cs.ring_size_count[6].as_lit() == 3
    del cs.ring_size_count[6]
    assert cs.ring_size_count[6] is None


def test_atomconstraintsview_property_on_molecule():
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    mol.atoms[0].constraints.aromatic_valence = 1
    # a fresh view proves the write hit the molecule
    assert mol.atoms[0].constraints.aromatic_valence == AromaticValenceForm.Aromatic(NumForm.Lit(1))


def test_atomconstraintsview_ring_size_count_on_molecule():
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    mol.atoms[0].constraints.ring_size_count[6] = 3
    assert mol.atoms[0].constraints.ring_size_count[6].as_lit() == 3
    del mol.atoms[0].constraints.ring_size_count[6]
    assert mol.atoms[0].constraints.ring_size_count[6] is None


def test_aromatic_valence_form_aromatic_int():
    assert AromaticValenceForm.Aromatic(1) == AromaticValenceForm.Aromatic(NumForm.Lit(1))


def test_multicenter_valence_form_multicenter_int():
    assert MulticenterValenceForm.Multicenter(2) == MulticenterValenceForm.Multicenter(NumForm.Lit(2))


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
    assert AtomConstraintForm.Valence(NumForm.Lit(4)) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )
    assert AtomConstraintForm.Valence(NumForm.Lit(4)) != AtomConstraintForm.Valence(
        NumForm.Lit(5)
    )
    assert AtomConstraintForm.Valence(NumForm.Lit(4)) != AtomConstraintForm.Degree(
        NumForm.Lit(4)
    )
    assert (
        len(
            {
                AtomConstraintForm.Valence(NumForm.Lit(4)),
                AtomConstraintForm.Valence(NumForm.Lit(4)),
            }
        )
        == 1
    )


def test_atomconstraint_repr():
    x = AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.Size(6), NumForm.Lit(1)))
    env = {
        "AtomConstraintForm": AtomConstraintForm,
        "RingMembershipForm": RingMembershipForm,
        "RingScope": RingScope,
        "NumForm": NumForm,
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


def test_atom_constraints_form_eq_repr():
    a = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    b = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    assert a == b
    assert a != AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(5))])
    assert repr(a) == "AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])"


def test_atom_constraints_form_unhashable():
    # mutable container: value-equal but unhashable, like AtomForm
    with pytest.raises(TypeError):
        hash(AtomConstraintsForm([]))


def test_aromatic_valence_form_eq():
    assert AromaticValenceForm.Aromatic(1) == AromaticValenceForm.Aromatic(1)
    assert AromaticValenceForm.Aromatic(1) != AromaticValenceForm.NotAromatic()


def test_ring_membership_form_eq_repr():
    a = RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
    assert a == RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))
    assert a != RingMembershipForm(RingScope.All(), NumForm.Lit(1))
    assert repr(a) == "RingMembershipForm(RingScope.Size(6), NumForm.Lit(1))"


def test_atomconstraintsview_repr():
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    assert repr(mol.atoms[0].constraints) == "AtomConstraintsView(0 entries)"


def test_atomconstraints_getitem_delitem():
    cs = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    assert AtomConstraintKey.Valence() in cs
    assert cs[AtomConstraintKey.Valence()] == AtomConstraintForm.Valence(NumForm.Lit(4))
    with pytest.raises(KeyError):
        cs[AtomConstraintKey.Degree()]
    del cs[AtomConstraintKey.Valence()]
    assert not cs
    with pytest.raises(KeyError):
        del cs[AtomConstraintKey.Valence()]


def test_atomconstraints_update_iterable():
    cs = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    cs.update(
        [
            AtomConstraintForm.Valence(NumForm.Lit(2)),
            AtomConstraintForm.Degree(NumForm.Lit(3)),
        ]
    )
    assert len(cs) == 2
    assert cs.valence == NumForm.Lit(2)


def test_atomconstraints_update_container():
    cs = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    cs.update(AtomConstraintsForm([AtomConstraintForm.Degree(NumForm.Lit(3))]))
    assert len(cs) == 2


def test_atomconstraintsview_getitem_delitem():
    mol = Molecule.from_entries(
        [
            AtomForm(
                Element("C"),
                constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
            )
        ]
    )
    cs = mol.atoms[0].constraints
    assert AtomConstraintKey.Valence() in cs
    assert cs[AtomConstraintKey.Valence()] == AtomConstraintForm.Valence(NumForm.Lit(4))
    with pytest.raises(KeyError):
        cs[AtomConstraintKey.Degree()]
    del mol.atoms[0].constraints[AtomConstraintKey.Valence()]
    assert not mol.atoms[0].constraints
    with pytest.raises(KeyError):
        del mol.atoms[0].constraints[AtomConstraintKey.Valence()]


def test_atomconstraintsview_update_from_view():
    src = AtomForm(
        Element("C"),
        constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
    )
    mol = Molecule.from_entries([AtomForm(Element("C"))])
    mol.atoms[0].constraints.update(src.constraints)
    assert AtomConstraintKey.Valence() in mol.atoms[0].constraints


def test_atom_form_set_constraints_from_value():
    dst = AtomForm(Element("N"))
    dst.constraints = AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))])
    assert dst.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )


def test_atom_form_set_constraints_from_view():
    src = AtomForm(
        Element("C"),
        constraints=AtomConstraintsForm([AtomConstraintForm.Valence(NumForm.Lit(4))]),
    )
    dst = AtomForm(Element("N"))
    dst.constraints = src.constraints  # RHS is a live view, not a value container
    assert dst.constraints.get(AtomConstraintKey.Valence()) == AtomConstraintForm.Valence(
        NumForm.Lit(4)
    )


def test_ringsizecounts_len_iter_contains():
    cs = AtomConstraintsForm([])
    cs.ring_size_count[6] = 3
    cs.ring_size_count[5] = 1
    rsc = cs.ring_size_count
    assert len(rsc) == 2
    assert sorted(rsc) == [5, 6]
    assert 6 in rsc
    assert 4 not in rsc
