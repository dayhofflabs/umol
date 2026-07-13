import pytest

from umol import (
    AtomAst,
    BooleanAst,
    Element,
    MoleculeAst,
    NoncovalentBondAst,
    NoncovalentBondConstraintAst,
    NoncovalentBondConstraintKey,
    NoncovalentBondConstraintsAst,
    NoncovalentBondKind,
    NoncovalentBondKindAst,
    ParseError,
)


def hbond_molecule():
    # two oxygens (atom ids 0-1), one hydrogen bond over them (noncovalent id 0)
    return MoleculeAst.from_parts(
        [AtomAst(Element("O")) for _ in range(2)],
        noncovalent_bonds=[([0, 1], NoncovalentBondAst(NoncovalentBondKind.HydrogenBond))],
    )


def test_noncovalentbondkind_members():
    kinds = [
        NoncovalentBondKind.HydrogenBond,
        NoncovalentBondKind.HalogenBond,
        NoncovalentBondKind.ChalcogenBond,
        NoncovalentBondKind.Ionic,
        NoncovalentBondKind.VanDerWaals,
    ]
    assert len(set(kinds)) == 5
    assert NoncovalentBondKind.HydrogenBond == NoncovalentBondKind.HydrogenBond
    assert NoncovalentBondKind.HydrogenBond != NoncovalentBondKind.Ionic


def test_noncovalentbondkindast_as_lit():
    assert NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic).as_lit() == (
        NoncovalentBondKind.Ionic
    )
    assert NoncovalentBondKindAst.Undetermined().as_lit() is None


def test_noncovalentbondast_new():
    bond = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    assert bond.kind == NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)
    assert len(bond.constraints) == 0


def test_noncovalentbondast_new_constraints_kwarg():
    bond = NoncovalentBondAst(
        NoncovalentBondKind.HalogenBond,
        constraints=NoncovalentBondConstraintsAst(
            [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
        ),
    )
    assert bond.constraints.intramolecular == BooleanAst.Lit(True)


def test_noncovalentbondast_kind_setter():
    bond = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    bond.kind = NoncovalentBondKind.Ionic
    assert bond.kind == NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic)


@pytest.mark.parametrize("dsl", ["Hbd", "Hbd#I", "Hbd#I!", "*"])
def test_noncovalentbondast_parse_roundtrip(dsl):
    bond = NoncovalentBondAst.parse(dsl)
    assert str(bond) == dsl
    assert repr(bond) == f"NoncovalentBondAst.parse('{dsl}')"


def test_noncovalentbondast_parse_error():
    with pytest.raises(ParseError):
        NoncovalentBondAst.parse("z")


def test_noncovalentbondast_asdict():
    bond = NoncovalentBondAst.parse("Hbd#I")
    d = bond.asdict()
    assert set(d.keys()) == {"kind", "constraints"}
    assert d["kind"] == NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)
    assert d["constraints"] == {"intramolecular": BooleanAst.Lit(True)}


def test_noncovalentbondast_set_constraints():
    bond = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    bond.constraints = NoncovalentBondConstraintsAst(
        [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
    )
    assert bond.constraints.intramolecular == BooleanAst.Lit(False)


def test_noncovalentbondast_constraints_self_assign():
    # regression: assigning the bond's own constraints view back to it is a no-op, not a panic
    bond = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    bond.constraints.intramolecular = True
    bond.constraints = bond.constraints
    bond.constraints.update(bond.constraints)
    assert bond.constraints.intramolecular == BooleanAst.Lit(True)


def test_noncovalentbondconstraints_intramolecular():
    constraints = NoncovalentBondConstraintsAst([])
    assert constraints.intramolecular == BooleanAst.Undetermined()
    constraints.intramolecular = True
    assert constraints.intramolecular == BooleanAst.Lit(True)


def test_noncovalentbondconstraints_mapping_ops():
    constraints = NoncovalentBondConstraintsAst([])
    constraints.set(NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)))
    assert len(constraints) == 1
    assert NoncovalentBondConstraintKey.Intramolecular() in constraints
    assert constraints[NoncovalentBondConstraintKey.Intramolecular()] == (
        NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))
    )
    assert [key for key in constraints] == [NoncovalentBondConstraintKey.Intramolecular()]
    del constraints[NoncovalentBondConstraintKey.Intramolecular()]
    assert len(constraints) == 0


def test_noncovalentbondconstraints_getitem_missing():
    constraints = NoncovalentBondConstraintsAst([])
    with pytest.raises(KeyError):
        constraints[NoncovalentBondConstraintKey.Intramolecular()]


def test_noncovalentbondconstraints_delitem_missing():
    constraints = NoncovalentBondConstraintsAst([])
    with pytest.raises(KeyError):
        del constraints[NoncovalentBondConstraintKey.Intramolecular()]


def test_noncovalentbondconstraintkey_intramolecular():
    key = NoncovalentBondConstraintKey.Intramolecular()
    assert key == NoncovalentBondConstraintKey.Intramolecular()
    assert key.__repr__().startswith("NoncovalentBondConstraintKey.Intramolecular")


def test_noncovalentbondconstraints_asdict():
    constraints = NoncovalentBondConstraintsAst(
        [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))]
    )
    assert constraints.asdict() == {"intramolecular": BooleanAst.Lit(True)}


def test_noncovalentbondview_fields():
    view = hbond_molecule().noncovalent_bonds[0]
    assert view.id == 0
    assert view.atom_ids == (0, 1)
    assert view.kind == NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)
    assert repr(view) == "NoncovalentBondView(id=0)"


def test_noncovalentbondview_set_kind():
    mol = hbond_molecule()
    mol.noncovalent_bonds[0].kind = NoncovalentBondKind.Ionic
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.noncovalent_bonds[0].kind == NoncovalentBondKindAst.Lit(
        NoncovalentBondKind.Ionic
    )


def test_noncovalentbondview_asdict():
    view = hbond_molecule().noncovalent_bonds[0]
    d = view.asdict()
    assert set(d.keys()) == {"kind", "constraints"}
    assert d["kind"] == NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)


def test_noncovalentbondview_constraints_write_through():
    mol = hbond_molecule()
    mol.noncovalent_bonds[0].constraints.intramolecular = True
    assert mol.noncovalent_bonds[0].constraints.intramolecular == BooleanAst.Lit(True)


def test_noncovalentbondview_set_constraints():
    mol = hbond_molecule()
    mol.noncovalent_bonds[0].constraints = NoncovalentBondConstraintsAst(
        [NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(False))]
    )
    assert mol.noncovalent_bonds[0].constraints.intramolecular == BooleanAst.Lit(False)


def test_noncovalentbondviews_len_getitem():
    bonds = hbond_molecule().noncovalent_bonds
    assert len(bonds) == 1
    assert bonds[0].id == 0
    assert bonds[-1].id == 0
    with pytest.raises(IndexError):
        bonds[5]
    with pytest.raises(IndexError):
        bonds[-2]


def test_noncovalentbondviews_setitem():
    mol = hbond_molecule()
    mol.noncovalent_bonds[0] = NoncovalentBondAst(NoncovalentBondKind.Ionic)
    view = mol.noncovalent_bonds[0]
    # value replaced, endpoints preserved
    assert view.kind == NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic)
    assert view.atom_ids == (0, 1)


def test_noncovalentbondviews_setitem_out_of_range():
    with pytest.raises(IndexError):
        hbond_molecule().noncovalent_bonds[5] = NoncovalentBondAst(
            NoncovalentBondKind.HydrogenBond
        )


def test_noncovalentbondviews_iter():
    ids = [view.id for view in hbond_molecule().noncovalent_bonds]
    assert ids == [0]


def test_noncovalentbondviews_of():
    # three oxygens, one hydrogen bond over (0, 1); atom 2 isolated
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("O")) for _ in range(3)],
        noncovalent_bonds=[([0, 1], NoncovalentBondAst(NoncovalentBondKind.HydrogenBond))],
    )
    # unordered pair — both orders find the same bond
    assert mol.noncovalent_bonds.of(0, 1).id == 0
    assert mol.noncovalent_bonds.of(1, 0).id == 0
    # no bond between 0 and the isolated atom 2
    assert mol.noncovalent_bonds.of(0, 2) is None


def test_noncovalentbondviews_incident():
    # three oxygens, one hydrogen bond over (0, 1); atom 2 isolated
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("O")) for _ in range(3)],
        noncovalent_bonds=[([0, 1], NoncovalentBondAst(NoncovalentBondKind.HydrogenBond))],
    )
    assert [view.id for view in mol.noncovalent_bonds.incident(0)] == [0]
    assert mol.noncovalent_bonds.incident(2) == []


def test_noncovalentbondviews_repr():
    assert repr(hbond_molecule().noncovalent_bonds) == "NoncovalentBondViews(len=1)"
