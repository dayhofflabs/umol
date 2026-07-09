from umol import (
    AromaticValenceAst,
    AtomAst,
    AtomConstraintAst,
    AtomConstraintKey,
    AtomConstraintsAst,
    AtomId,
    Element,
    MoleculeAst,
    MulticenterValenceAst,
    RingMembershipAst,
    RingScope,
    TetrahedralStereoAst,
    ValueAst,
)


def test_aromaticvalenceast_aromatic():
    match AromaticValenceAst.Aromatic(ValueAst.Lit(1)):
        case AromaticValenceAst.Aromatic(v):
            match v:
                case ValueAst.Lit(n):
                    assert n == 1
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_aromaticvalenceast_not_aromatic():
    match AromaticValenceAst.NotAromatic():
        case AromaticValenceAst.NotAromatic():
            pass
        case _:
            raise AssertionError


def test_multicentervalenceast_multicenter():
    match MulticenterValenceAst.Multicenter(ValueAst.Lit(2)):
        case MulticenterValenceAst.Multicenter(v):
            match v:
                case ValueAst.Lit(n):
                    assert n == 2
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_ringscope_size():
    match RingScope.Size(6):
        case RingScope.Size(s):
            assert s == 6
        case _:
            raise AssertionError


def test_ringmembershipast_fields():
    rm = RingMembershipAst(RingScope.All(), ValueAst.Lit(2))
    match rm.scope:
        case RingScope.All():
            pass
        case _:
            raise AssertionError
    match rm.count:
        case ValueAst.Lit(n):
            assert n == 2
        case _:
            raise AssertionError


def test_ringmembershipast_int_count():
    match RingMembershipAst(RingScope.All(), 2).count:
        case ValueAst.Lit(n):
            assert n == 2
        case _:
            raise AssertionError


def test_atomconstraint_key_valence():
    match AtomConstraintAst.Valence(ValueAst.Lit(4)).key:
        case AtomConstraintKey.Valence():
            pass
        case _:
            raise AssertionError


def test_atomconstraint_key_tetrahedral_stereo():
    constraint = AtomConstraintAst.TetrahedralStereo(TetrahedralStereoAst.NotStereo())
    match constraint.key:
        case AtomConstraintKey.TetrahedralStereo():
            pass
        case _:
            raise AssertionError


def test_atomconstraint_key_ring_membership():
    constraint = AtomConstraintAst.RingMembership(
        RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1))
    )
    match constraint.key:
        case AtomConstraintKey.RingMembership(RingScope.Size(s)):
            assert s == 6
        case _:
            raise AssertionError


def test_atomconstraints_iter():
    constraints = AtomConstraintsAst(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
        ]
    )
    assert len(constraints) == 2
    keys = set()
    for constraint in constraints:
        match constraint.key:
            case AtomConstraintKey.Valence():
                keys.add("valence")
            case AtomConstraintKey.Degree():
                keys.add("degree")
            case _:
                raise AssertionError
    assert keys == {"valence", "degree"}


def test_atomconstraints_get():
    constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    assert constraints.contains(AtomConstraintKey.Valence())
    assert not constraints.contains(AtomConstraintKey.Degree())
    assert constraints.get(AtomConstraintKey.Degree()) is None
    match constraints.get(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError


def test_atomconstraints_get_ring_membership():
    constraints = AtomConstraintsAst(
        [AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))]
    )
    assert constraints.contains(AtomConstraintKey.RingMembership(RingScope.Size(6)))
    assert not constraints.contains(AtomConstraintKey.RingMembership(RingScope.All()))
    match constraints.get(AtomConstraintKey.RingMembership(RingScope.Size(6))):
        case AtomConstraintAst.RingMembership(rm):
            match rm.count:
                case ValueAst.Lit(n):
                    assert n == 1
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_atomconstraints_valence():
    constraints = AtomConstraintsAst(
        [
            AtomConstraintAst.Valence(ValueAst.Lit(4)),
            AtomConstraintAst.Degree(ValueAst.Lit(3)),
        ]
    )
    match constraints.valence():
        case ValueAst.Lit(n):
            assert n == 4
        case _:
            raise AssertionError
    match constraints.degree():
        case ValueAst.Lit(n):
            assert n == 3
        case _:
            raise AssertionError
    assert constraints.total_valence() is None
    assert constraints.aromatic_valence() is None


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
    match d["valence"]:
        case ValueAst.Lit(n):
            assert n == 4
        case _:
            raise AssertionError
    match d["degree"]:
        case ValueAst.Lit(n):
            assert n == 3
        case _:
            raise AssertionError
    match d["ring_count"]:
        case ValueAst.Lit(n):
            assert n == 2
        case _:
            raise AssertionError
    match d["ring_size_count_6"]:
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_atomconstraints_ring_size_count():
    constraints = AtomConstraintsAst(
        [AtomConstraintAst.RingMembership(RingMembershipAst(RingScope.Size(6), ValueAst.Lit(1)))]
    )
    match constraints.ring_size_count(6):
        case ValueAst.Lit(n):
            assert n == 1
        case _:
            raise AssertionError
    assert constraints.ring_size_count(5) is None
    assert constraints.ring_count() is None


def test_atomconstraintsast_set():
    constraints = AtomConstraintsAst([])
    constraints.set(AtomConstraintAst.Valence(ValueAst.Lit(4)))
    assert len(constraints) == 1
    match constraints.get(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError


def test_atomconstraintsast_remove():
    constraints = AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))])
    match constraints.remove(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError
    assert len(constraints) == 0
    assert constraints.remove(AtomConstraintKey.Valence()) is None


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
    match constraints.valence():
        case ValueAst.Lit(n):
            assert n == 3
        case _:
            raise AssertionError
    match constraints.degree():
        case ValueAst.Lit(n):
            assert n == 2
        case _:
            raise AssertionError


def test_atomconstraintsview_set():
    mol = MoleculeAst.from_atoms([AtomAst(Element("C"))])
    mol.atoms[AtomId(0)].constraints.set(
        AtomConstraintAst.AromaticValence(AromaticValenceAst.Aromatic(ValueAst.Lit(1)))
    )
    # a fresh view proves the write hit the molecule, not a transient copy
    constraints = mol.atoms[AtomId(0)].constraints
    assert len(constraints) == 1
    match constraints.get(AtomConstraintKey.AromaticValence()):
        case AtomConstraintAst.AromaticValence(AromaticValenceAst.Aromatic(ValueAst.Lit(n))):
            assert n == 1
        case _:
            raise AssertionError


def test_atomconstraintsview_remove():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_atoms([atom])
    match mol.atoms[AtomId(0)].constraints.remove(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError
    assert len(mol.atoms[AtomId(0)].constraints) == 0


def test_atomconstraintsview_update():
    mol = MoleculeAst.from_atoms([AtomAst(Element("C"))])
    mol.atoms[AtomId(0)].constraints.update(
        AtomConstraintsAst(
            [
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
                AtomConstraintAst.Degree(ValueAst.Lit(3)),
            ]
        )
    )
    constraints = mol.atoms[AtomId(0)].constraints
    assert len(constraints) == 2
    match constraints.valence():
        case ValueAst.Lit(n):
            assert n == 4
        case _:
            raise AssertionError


def test_atomconstraintsview_atom_backed_set():
    atom = AtomAst(Element("C"))
    atom.constraints.set(AtomConstraintAst.Valence(ValueAst.Lit(4)))
    # a fresh view proves the write mutated the standalone atom in place
    assert len(atom.constraints) == 1
    match atom.constraints.get(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError


def test_atomconstraintsview_atom_backed_remove():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    match atom.constraints.remove(AtomConstraintKey.Valence()):
        case AtomConstraintAst.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError
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
    match atom.constraints.valence():
        case ValueAst.Lit(n):
            assert n == 4
        case _:
            raise AssertionError


def test_atomconstraintsview_reads():
    atom = AtomAst(
        Element("C"),
        constraints=AtomConstraintsAst([AtomConstraintAst.Valence(ValueAst.Lit(4))]),
    )
    mol = MoleculeAst.from_atoms([atom])
    constraints = mol.atoms[AtomId(0)].constraints
    assert not constraints.is_empty()
    assert constraints.contains(AtomConstraintKey.Valence())
    assert constraints.get(AtomConstraintKey.Degree()) is None
    match constraints.valence():
        case ValueAst.Lit(n):
            assert n == 4
        case _:
            raise AssertionError
    assert set(constraints.asdict().keys()) == {"valence"}
    keys = set()
    for constraint in constraints:
        match constraint.key:
            case AtomConstraintKey.Valence():
                keys.add("valence")
            case _:
                raise AssertionError
    assert keys == {"valence"}
