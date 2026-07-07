from umol import (
    AromaticValenceAst,
    AtomConstraint,
    AtomConstraintKind,
    AtomConstraints,
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


def test_atomconstraint_valence():
    constraint = AtomConstraint.Valence(ValueAst.Lit(4))
    assert constraint.kind == AtomConstraintKind.Valence
    match constraint:
        case AtomConstraint.Valence(v):
            match v:
                case ValueAst.Lit(n):
                    assert n == 4
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_atomconstraint_tetrahedral_stereo():
    constraint = AtomConstraint.TetrahedralStereo(TetrahedralStereoAst.NotStereo())
    assert constraint.kind == AtomConstraintKind.TetrahedralStereo
    match constraint:
        case AtomConstraint.TetrahedralStereo(TetrahedralStereoAst.NotStereo()):
            pass
        case _:
            raise AssertionError


def test_atomconstraints_len_iter():
    constraints = AtomConstraints(
        [
            AtomConstraint.Valence(ValueAst.Lit(4)),
            AtomConstraint.Degree(ValueAst.Lit(3)),
        ]
    )
    assert len(constraints) == 2
    kinds = [constraint.kind for constraint in constraints]
    assert AtomConstraintKind.Valence in kinds
    assert AtomConstraintKind.Degree in kinds


def test_atomconstraints_get():
    constraints = AtomConstraints([AtomConstraint.Valence(ValueAst.Lit(4))])
    assert constraints.contains(AtomConstraintKind.Valence)
    assert not constraints.contains(AtomConstraintKind.Degree)
    assert constraints.get(AtomConstraintKind.Degree) is None
    match constraints.get(AtomConstraintKind.Valence):
        case AtomConstraint.Valence(ValueAst.Lit(n)):
            assert n == 4
        case _:
            raise AssertionError
