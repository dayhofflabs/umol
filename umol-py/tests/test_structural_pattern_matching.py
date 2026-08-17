from umol import (
    AromaticSystemConstraintForm,
    AromaticSystemDelta,
    AromaticSystemFieldChange,
    AromaticSystemForm,
    AtomConstraintForm,
    AtomDelta,
    AtomFieldChange,
    AtomForm,
    BondConstraintForm,
    BondDelta,
    BondFieldChange,
    BondForm,
    BooleanForm,
    Constraint,
    ConstraintDelta,
    DativeBondConstraintForm,
    DativeBondDelta,
    DativeBondFieldChange,
    DativeBondForm,
    Delta,
    ElectronCountsForm,
    Element,
    ElementForm,
    MoleculeConstraint,
    MulticenterBondConstraintForm,
    MulticenterBondDelta,
    MulticenterBondFieldChange,
    MulticenterBondForm,
    NoncovalentBondConstraintForm,
    NoncovalentBondDelta,
    NoncovalentBondFieldChange,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NoncovalentBondKindForm,
    NumForm,
    RelOp,
    RelationalConstraint,
    StereoAtomConstraintForm,
    StereoAtomDelta,
    StereoAtomFieldChange,
    StereoAtomForm,
    StereoBondConstraintForm,
    StereoBondDelta,
    StereoBondFieldChange,
    StereoBondForm,
    StereoConfigurationForm,
    StereoCoset,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereogenicityForm,
    UnpairedElectronsForm,
)


def test_relop_match():
    match RelOp.Ne:
        case RelOp.Ne:
            pass
        case _:
            raise AssertionError


def test_constraint_pattern_match():
    constraint = Constraint.And(
        [
            Constraint.Relational(RelationalConstraint.DativeBondDonor(3, 5)),
            Constraint.Not(Constraint.Molecule(MoleculeConstraint.Connected(None))),
        ]
    )

    match constraint:
        case Constraint.And(
            [
                Constraint.Relational(RelationalConstraint.DativeBondDonor(bond, atom)),
                Constraint.Not(Constraint.Molecule(MoleculeConstraint.Connected(atoms))),
            ]
        ):
            assert (bond, atom, atoms) == (3, 5, None)
        case _:
            raise AssertionError("constraint tree did not match its structural variants")


def test_atomfieldchange_match_scalar():
    change = AtomFieldChange.ImplicitHydrogens(
        old=NumForm.Lit(3), new=NumForm.Lit(2)
    )

    match change:
        case AtomFieldChange.ImplicitHydrogens(old, new):
            assert (old, new) == (NumForm.Lit(3), NumForm.Lit(2))
        case _:
            raise AssertionError("atom field change did not match its scalar variant")


def test_atomfieldchange_match_structured():
    change = AtomFieldChange.Element(
        old=ElementForm.Lit(Element("C")),
        new=ElementForm.Lit(Element("N")),
    )

    match change:
        case AtomFieldChange.Element(old=old, new=new):
            assert (old, new) == (
                ElementForm.Lit(Element("C")),
                ElementForm.Lit(Element("N")),
            )
        case _:
            raise AssertionError("atom field change did not match its structured variant")


def test_bondfieldchange_match():
    change = BondFieldChange.UnpairedElectrons(
        old=UnpairedElectronsForm(0, 1), new=UnpairedElectronsForm(1, 2)
    )

    match change:
        case BondFieldChange.UnpairedElectrons(old, new):
            assert (old, new) == (
                UnpairedElectronsForm(0, 1),
                UnpairedElectronsForm(1, 2),
            )
        case _:
            raise AssertionError("bond field change did not match its variant")


def test_dativebondfieldchange_match():
    change = DativeBondFieldChange.Order(
        old=NumForm.Lit(1), new=NumForm.Lit(2)
    )

    match change:
        case DativeBondFieldChange.Order(old=old, new=new):
            assert (old, new) == (NumForm.Lit(1), NumForm.Lit(2))
        case _:
            raise AssertionError("dative bond field change did not match its variant")

    inverse = change.inverse()
    assert isinstance(inverse, DativeBondFieldChange.Order)
    assert inverse == DativeBondFieldChange.Order(
        old=NumForm.Lit(2), new=NumForm.Lit(1)
    )


def test_aromaticsystemfieldchange_match():
    change = AromaticSystemFieldChange.Electrons(
        old=ElectronCountsForm.Lit([2, 0, 2]),
        new=ElectronCountsForm.Lit([1, 1, 1]),
    )

    match change:
        case AromaticSystemFieldChange.Electrons(old=old, new=new):
            assert (old, new) == (
                ElectronCountsForm.Lit([2, 0, 2]),
                ElectronCountsForm.Lit([1, 1, 1]),
            )
        case _:
            raise AssertionError("aromatic field change did not match its variant")


def test_multicenterbondfieldchange_match():
    change = MulticenterBondFieldChange.Electrons(
        old=ElectronCountsForm.Lit([1, 0, 1]),
        new=ElectronCountsForm.Lit([2, 0, 1]),
    )

    match change:
        case MulticenterBondFieldChange.Electrons(old, new):
            assert (old, new) == (
                ElectronCountsForm.Lit([1, 0, 1]),
                ElectronCountsForm.Lit([2, 0, 1]),
            )
        case _:
            raise AssertionError("multicenter field change did not match its variant")

    inverse = change.inverse()
    assert isinstance(inverse, MulticenterBondFieldChange.Electrons)
    assert inverse.old == ElectronCountsForm.Lit([2, 0, 1])
    assert inverse.new == ElectronCountsForm.Lit([1, 0, 1])


def test_noncovalentbondfieldchange_match():
    change = NoncovalentBondFieldChange.Kind(
        old=NoncovalentBondKindForm.Undetermined(),
        new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond),
    )

    match change:
        case NoncovalentBondFieldChange.Kind(old=old, new=new):
            assert old == NoncovalentBondKindForm.Undetermined()
            assert new == NoncovalentBondKindForm.Lit(
                NoncovalentBondKind.HydrogenBond
            )
        case _:
            raise AssertionError("noncovalent field change did not match its variant")


def test_stereoatomfieldchange_match_positional():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Undetermined()
        ),
        new=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(1)
        ),
    )

    match change:
        case StereoAtomFieldChange.Configuration(old, new):
            assert old == StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            )
            assert new == StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(1)
            )
        case _:
            raise AssertionError("stereo atom field change did not match its variant")


def test_stereoatomfieldchange_match_named():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationForm.Undetermined(),
        new=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Undetermined()
        ),
    )

    match change:
        case StereoAtomFieldChange.Configuration(old=old, new=new):
            assert old == StereoConfigurationForm.Undetermined()
            assert new == StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            )
        case _:
            raise AssertionError("stereo atom field change did not match its variant")


def test_stereobondfieldchange_match_positional():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Undetermined()
        ),
        new=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Lit(1)
        ),
    )

    match change:
        case StereoBondFieldChange.Configuration(old, new):
            assert old == StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Undetermined()
            )
            assert new == StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(1)
            )
        case _:
            raise AssertionError("stereo bond field change did not match its variant")


def test_stereobondfieldchange_match_named():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationForm.Undetermined(),
        new=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Undetermined()
        ),
    )

    match change:
        case StereoBondFieldChange.Configuration(old=old, new=new):
            assert old == StereoConfigurationForm.Undetermined()
            assert new == StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Undetermined()
            )
        case _:
            raise AssertionError("stereo bond field change did not match its variant")


def test_atomdelta_add_match():
    delta = AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))

    match delta:
        case AtomDelta.Add(id=id, attributes=attributes):
            assert id == 3
            assert attributes == AtomForm(Element("C"))
        case _:
            raise AssertionError("atom delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AtomDelta.Remove)
    assert inverse.id == 3
    assert inverse.attributes == AtomForm(Element("C"))
    assert inverse.inverse() == delta


def test_atomdelta_modifyfield_match():
    delta = AtomDelta.ModifyField(
        id=3,
        change=AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1)),
    )

    match delta:
        case AtomDelta.ModifyField(id, change):
            assert id == 3
            assert change == AtomFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(-1)
            )
        case _:
            raise AssertionError("atom delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AtomDelta.ModifyField)
    assert inverse.inverse() == delta


def test_atomdelta_modifyconstraint_match():
    delta = AtomDelta.ModifyConstraint(
        id=3,
        old=None,
        new=AtomConstraintForm.Valence(NumForm.Lit(4)),
    )

    match delta:
        case AtomDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 3
            assert old is None
            assert new == AtomConstraintForm.Valence(NumForm.Lit(4))
        case _:
            raise AssertionError("atom delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AtomDelta.ModifyConstraint)
    assert inverse.old == AtomConstraintForm.Valence(NumForm.Lit(4))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_bonddelta_add_match():
    delta = BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm(1))

    match delta:
        case BondDelta.Add(id=id, atoms=atoms, attributes=attributes):
            assert id == 2
            assert atoms == (5, 1)
            assert attributes == BondForm(1)
        case _:
            raise AssertionError("bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, BondDelta.Remove)
    assert inverse.id == 2
    assert inverse.atoms == (5, 1)
    assert inverse.attributes == BondForm(1)
    assert inverse.inverse() == delta


def test_bonddelta_modifyfield_match():
    delta = BondDelta.ModifyField(
        id=2,
        change=BondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2)),
    )

    match delta:
        case BondDelta.ModifyField(id, change):
            assert id == 2
            assert change == BondFieldChange.Order(
                old=NumForm.Lit(1), new=NumForm.Lit(2)
            )
        case _:
            raise AssertionError("bond delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, BondDelta.ModifyField)
    assert inverse.inverse() == delta


def test_bonddelta_modifyconstraint_match():
    delta = BondDelta.ModifyConstraint(
        id=2,
        old=None,
        new=BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
    )

    match delta:
        case BondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 2
            assert old is None
            assert new == BondConstraintForm.Aromatic(BooleanForm.Lit(True))
        case _:
            raise AssertionError("bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, BondDelta.ModifyConstraint)
    assert inverse.old == BondConstraintForm.Aromatic(BooleanForm.Lit(True))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_dativebonddelta_add_match():
    delta = DativeBondDelta.Add(
        id=1,
        donors=[4, 2, 4],
        acceptor=3,
        attributes=DativeBondForm(1),
    )

    match delta:
        case DativeBondDelta.Add(id=id, donors=donors, acceptor=acceptor, attributes=attributes):
            assert id == 1
            assert donors == [4, 2, 4]
            assert acceptor == 3
            assert attributes == DativeBondForm(1)
        case _:
            raise AssertionError("dative bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, DativeBondDelta.Remove)
    assert inverse.id == 1
    assert inverse.donors == [4, 2, 4]
    assert inverse.acceptor == 3
    assert inverse.attributes == DativeBondForm(1)
    assert inverse.inverse() == delta


def test_dativebonddelta_modifyfield_match():
    delta = DativeBondDelta.ModifyField(
        id=1,
        change=DativeBondFieldChange.Order(
            old=NumForm.Lit(1), new=NumForm.Lit(2)
        ),
    )

    match delta:
        case DativeBondDelta.ModifyField(id, change):
            assert id == 1
            assert change == DativeBondFieldChange.Order(
                old=NumForm.Lit(1), new=NumForm.Lit(2)
            )
        case _:
            raise AssertionError("dative bond delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, DativeBondDelta.ModifyField)
    assert inverse.inverse() == delta


def test_dativebonddelta_modifyconstraint_match():
    delta = DativeBondDelta.ModifyConstraint(
        id=1,
        old=None,
        new=DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
    )

    match delta:
        case DativeBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 1
            assert old is None
            assert new == DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))
        case _:
            raise AssertionError("dative bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, DativeBondDelta.ModifyConstraint)
    assert inverse.old == DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_aromaticsystemdelta_add_match():
    delta = AromaticSystemDelta.Add(
        id=2,
        atoms=[4, 2, 4],
        attributes=AromaticSystemForm([1, 1, 1]),
    )

    match delta:
        case AromaticSystemDelta.Add(id=id, atoms=atoms, attributes=attributes):
            assert id == 2
            assert atoms == [4, 2, 4]
            assert attributes == AromaticSystemForm([1, 1, 1])
        case _:
            raise AssertionError("aromatic system delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AromaticSystemDelta.Remove)
    assert inverse.id == 2
    assert inverse.atoms == [4, 2, 4]
    assert inverse.attributes == AromaticSystemForm([1, 1, 1])
    assert inverse.inverse() == delta


def test_aromaticsystemdelta_modifyfield_match():
    delta = AromaticSystemDelta.ModifyField(
        id=2,
        change=AromaticSystemFieldChange.Charge(
            old=NumForm.Lit(0), new=NumForm.Lit(-1)
        ),
    )

    match delta:
        case AromaticSystemDelta.ModifyField(id, change):
            assert id == 2
            assert change == AromaticSystemFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(-1)
            )
        case _:
            raise AssertionError("aromatic system delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AromaticSystemDelta.ModifyField)
    assert inverse.inverse() == delta


def test_aromaticsystemdelta_modifyconstraint_match():
    delta = AromaticSystemDelta.ModifyConstraint(
        id=2,
        old=None,
        new=AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6)),
    )

    match delta:
        case AromaticSystemDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 2
            assert old is None
            assert new == AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))
        case _:
            raise AssertionError(
                "aromatic system delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, AromaticSystemDelta.ModifyConstraint)
    assert inverse.old == AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_multicenterbonddelta_add_match():
    delta = MulticenterBondDelta.Add(
        id=3,
        atoms=[4, 2, 4],
        attributes=MulticenterBondForm([1, 1, 1]),
    )

    match delta:
        case MulticenterBondDelta.Add(id=id, atoms=atoms, attributes=attributes):
            assert id == 3
            assert atoms == [4, 2, 4]
            assert attributes == MulticenterBondForm([1, 1, 1])
        case _:
            raise AssertionError("multicenter bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, MulticenterBondDelta.Remove)
    assert inverse.id == 3
    assert inverse.atoms == [4, 2, 4]
    assert inverse.attributes == MulticenterBondForm([1, 1, 1])
    assert inverse.inverse() == delta


def test_multicenterbonddelta_modifyfield_match():
    delta = MulticenterBondDelta.ModifyField(
        id=3,
        change=MulticenterBondFieldChange.Charge(
            old=NumForm.Lit(0), new=NumForm.Lit(-1)
        ),
    )

    match delta:
        case MulticenterBondDelta.ModifyField(id, change):
            assert id == 3
            assert change == MulticenterBondFieldChange.Charge(
                old=NumForm.Lit(0), new=NumForm.Lit(-1)
            )
        case _:
            raise AssertionError("multicenter bond delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, MulticenterBondDelta.ModifyField)
    assert inverse.inverse() == delta


def test_multicenterbonddelta_modifyconstraint_match():
    delta = MulticenterBondDelta.ModifyConstraint(
        id=3,
        old=None,
        new=MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6)),
    )

    match delta:
        case MulticenterBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 3
            assert old is None
            assert new == MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6))
        case _:
            raise AssertionError(
                "multicenter bond delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, MulticenterBondDelta.ModifyConstraint)
    assert inverse.old == MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_noncovalentbonddelta_add_match():
    delta = NoncovalentBondDelta.Add(
        id=4,
        atoms=(5, 2),
        attributes=NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
    )

    match delta:
        case NoncovalentBondDelta.Add(id=id, atoms=atoms, attributes=attributes):
            assert id == 4
            assert atoms == (5, 2)
            assert attributes == NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
        case _:
            raise AssertionError("noncovalent bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, NoncovalentBondDelta.Remove)
    assert inverse.id == 4
    assert inverse.atoms == (5, 2)
    assert inverse.attributes == NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    assert inverse.inverse() == delta


def test_noncovalentbonddelta_modifyfield_match():
    delta = NoncovalentBondDelta.ModifyField(
        id=4,
        change=NoncovalentBondFieldChange.Kind(
            old=NoncovalentBondKindForm.Undetermined(),
            new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond),
        ),
    )

    match delta:
        case NoncovalentBondDelta.ModifyField(id, change):
            assert id == 4
            assert change == NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindForm.Undetermined(),
                new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond),
            )
        case _:
            raise AssertionError("noncovalent bond delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, NoncovalentBondDelta.ModifyField)
    assert inverse.inverse() == delta


def test_noncovalentbonddelta_modifyconstraint_match():
    delta = NoncovalentBondDelta.ModifyConstraint(
        id=4,
        old=None,
        new=NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)),
    )

    match delta:
        case NoncovalentBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 4
            assert old is None
            assert new == NoncovalentBondConstraintForm.Intramolecular(
                BooleanForm.Lit(True)
            )
        case _:
            raise AssertionError(
                "noncovalent bond delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, NoncovalentBondDelta.ModifyConstraint)
    assert inverse.old == NoncovalentBondConstraintForm.Intramolecular(
        BooleanForm.Lit(True)
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_stereoatomdelta_add_match():
    delta = StereoAtomDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
        ],
        attributes=StereoAtomForm(
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(0)
            )
        ),
    )

    match delta:
        case StereoAtomDelta.Add(id=id, site=site, ligands=ligands, attributes=attributes):
            assert id == 5
            assert site == 3
            assert ligands == [
                StereoLigand(4, StereoLigandKind.Atom),
                StereoLigand(2, StereoLigandKind.LonePair),
            ]
            assert attributes.configuration == StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(0)
            )
        case _:
            raise AssertionError("stereo atom delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoAtomDelta.Remove)
    assert inverse.site == 3
    assert inverse.ligands == delta.ligands
    assert inverse.inverse() == delta


def test_stereoatomdelta_modifyfield_match():
    delta = StereoAtomDelta.ModifyField(
        id=5,
        change=StereoAtomFieldChange.Configuration(
            old=StereoConfigurationForm.Undetermined(),
            new=StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(0)
            ),
        ),
    )

    match delta:
        case StereoAtomDelta.ModifyField(id, change):
            assert id == 5
            assert change == StereoAtomFieldChange.Configuration(
                old=StereoConfigurationForm.Undetermined(),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(0)
                ),
            )
        case _:
            raise AssertionError("stereo atom delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoAtomDelta.ModifyField)
    assert inverse != delta
    assert inverse.inverse() == delta


def test_stereoatomdelta_modifyconstraint_match():
    delta = StereoAtomDelta.ModifyConstraint(
        id=5,
        kind=StereoKind.Tetrahedral,
        old=None,
        new=StereoAtomConstraintForm.Stereogenicity(
            StereogenicityForm.Undetermined()
        ),
    )

    match delta:
        case StereoAtomDelta.ModifyConstraint(
            id=id, kind=kind, old=old, new=new
        ):
            assert id == 5
            assert kind is StereoKind.Tetrahedral
            assert old is None
            assert new == StereoAtomConstraintForm.Stereogenicity(
                StereogenicityForm.Undetermined()
            )
        case _:
            raise AssertionError("stereo atom delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoAtomDelta.ModifyConstraint)
    assert inverse.kind is StereoKind.Tetrahedral
    assert inverse.old == StereoAtomConstraintForm.Stereogenicity(
        StereogenicityForm.Undetermined()
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_stereobonddelta_add_match():
    delta = StereoBondDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
        ],
        attributes=StereoBondForm(
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(0)
            )
        ),
    )

    match delta:
        case StereoBondDelta.Add(id=id, site=site, ligands=ligands, attributes=attributes):
            assert id == 5
            assert site == 3
            assert ligands == [
                StereoLigand(4, StereoLigandKind.Atom),
                StereoLigand(2, StereoLigandKind.LonePair),
            ]
            assert attributes.configuration == StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(0)
            )
        case _:
            raise AssertionError("stereo bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoBondDelta.Remove)
    assert inverse.site == 3
    assert inverse.ligands == delta.ligands
    assert inverse.inverse() == delta


def test_stereobonddelta_modifyfield_match():
    delta = StereoBondDelta.ModifyField(
        id=5,
        change=StereoBondFieldChange.Configuration(
            old=StereoConfigurationForm.Undetermined(),
            new=StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(0)
            ),
        ),
    )

    match delta:
        case StereoBondDelta.ModifyField(id, change):
            assert id == 5
            assert change == StereoBondFieldChange.Configuration(
                old=StereoConfigurationForm.Undetermined(),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                ),
            )
        case _:
            raise AssertionError("stereo bond delta did not match its field variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoBondDelta.ModifyField)
    assert inverse != delta
    assert inverse.inverse() == delta


def test_stereobonddelta_modifyconstraint_match():
    delta = StereoBondDelta.ModifyConstraint(
        id=5,
        kind=StereoKind.CisTrans,
        old=None,
        new=StereoBondConstraintForm.Stereogenicity(
            StereogenicityForm.Undetermined()
        ),
    )

    match delta:
        case StereoBondDelta.ModifyConstraint(
            id=id, kind=kind, old=old, new=new
        ):
            assert id == 5
            assert kind is StereoKind.CisTrans
            assert old is None
            assert new == StereoBondConstraintForm.Stereogenicity(
                StereogenicityForm.Undetermined()
            )
        case _:
            raise AssertionError("stereo bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoBondDelta.ModifyConstraint)
    assert inverse.kind is StereoKind.CisTrans
    assert inverse.old == StereoBondConstraintForm.Stereogenicity(
        StereogenicityForm.Undetermined()
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_constraintdelta_add_match():
    delta = ConstraintDelta.Add(
        constraint=Constraint.Atom(
            3,
            AtomConstraintForm.Degree(NumForm.Lit(2)),
        )
    )

    match delta:
        case ConstraintDelta.Add(Constraint.Atom(atom_id, constraint)):
            assert (atom_id, constraint) == (
                3,
                AtomConstraintForm.Degree(NumForm.Lit(2)),
            )
        case _:
            raise AssertionError("constraint delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, ConstraintDelta.Remove)
    assert inverse.constraint == delta.constraint
    assert inverse.inverse() == delta


def test_constraintdelta_remove_match():
    delta = ConstraintDelta.Remove(
        constraint=Constraint.And(
            [
                Constraint.Atom(
                    7,
                    AtomConstraintForm.Valence(NumForm.Lit(4)),
                ),
                Constraint.Not(Constraint.Or([])),
            ]
        )
    )

    match delta:
        case ConstraintDelta.Remove(
            constraint=Constraint.And(
                [
                    Constraint.Atom(atom_id, constraint),
                    Constraint.Not(Constraint.Or([])),
                ]
            )
        ):
            assert (atom_id, constraint) == (
                7,
                AtomConstraintForm.Valence(NumForm.Lit(4)),
            )
        case _:
            raise AssertionError("constraint delta did not match its remove variant")

    inverse = delta.inverse()
    assert isinstance(inverse, ConstraintDelta.Add)
    assert inverse.constraint == delta.constraint
    assert inverse.inverse() == delta


def test_delta_atom_match():
    delta = Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))))

    match delta:
        case Delta.Atom(AtomDelta.Add(id=atom_id, attributes=attributes)):
            assert atom_id == 3
            assert attributes == AtomForm(Element("C"))
        case _:
            raise AssertionError("delta did not match its atom family and add operation")

    inverse = delta.inverse()
    assert type(inverse) is Delta.Atom
    assert type(inverse._0) is AtomDelta.Remove
    assert inverse._0 is not delta._0
    assert inverse != delta
    assert inverse.inverse() == delta


def test_delta_match():
    seen = []
    for delta in [
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Bond(BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm(1))),
        Delta.DativeBond(
            DativeBondDelta.Add(
                id=1,
                donors=[4, 2],
                acceptor=3,
                attributes=DativeBondForm(1),
            )
        ),
        Delta.AromaticSystem(
            AromaticSystemDelta.Add(
                id=2,
                atoms=[4, 2],
                attributes=AromaticSystemForm([1, 1]),
            )
        ),
        Delta.MulticenterBond(
            MulticenterBondDelta.Add(
                id=3,
                atoms=[4, 2],
                attributes=MulticenterBondForm([1, 1]),
            )
        ),
        Delta.NoncovalentBond(
            NoncovalentBondDelta.Add(
                id=4,
                atoms=(5, 2),
                attributes=NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
            )
        ),
        Delta.StereoAtom(
            StereoAtomDelta.Add(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoAtomForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.Tetrahedral,
                        StereoCoset.Lit(0),
                    )
                ),
            )
        ),
        Delta.StereoBond(
            StereoBondDelta.Add(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoBondForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.CisTrans,
                        StereoCoset.Lit(0),
                    )
                ),
            )
        ),
        Delta.Constraint(
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
            )
        ),
    ]:
        match delta:
            case Delta.Atom(AtomDelta.Add(id=3, attributes=attributes)):
                assert attributes == AtomForm(Element("C"))
                seen.append("atom")
            case Delta.Bond(BondDelta.Add(id=2, atoms=atoms, attributes=attributes)):
                assert (atoms, attributes) == ((5, 1), BondForm(1))
                seen.append("bond")
            case Delta.DativeBond(
                DativeBondDelta.Add(id=1, donors=donors, acceptor=3, attributes=attributes)
            ):
                assert (donors, attributes) == ([4, 2], DativeBondForm(1))
                seen.append("dative-bond")
            case Delta.AromaticSystem(
                AromaticSystemDelta.Add(id=2, atoms=atoms, attributes=attributes)
            ):
                assert (atoms, attributes) == ([4, 2], AromaticSystemForm([1, 1]))
                seen.append("aromatic-system")
            case Delta.MulticenterBond(
                MulticenterBondDelta.Add(id=3, atoms=atoms, attributes=attributes)
            ):
                assert (atoms, attributes) == ([4, 2], MulticenterBondForm([1, 1]))
                seen.append("multicenter-bond")
            case Delta.NoncovalentBond(
                NoncovalentBondDelta.Add(id=4, atoms=atoms, attributes=attributes)
            ):
                assert (atoms, attributes) == (
                    (5, 2),
                    NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
                )
                seen.append("noncovalent-bond")
            case Delta.StereoAtom(
                StereoAtomDelta.Add(id=5, site=3, ligands=ligands, attributes=attributes)
            ):
                assert ligands == [StereoLigand(4, StereoLigandKind.Atom)]
                assert attributes.configuration == StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral,
                    StereoCoset.Lit(0),
                )
                seen.append("stereo-atom")
            case Delta.StereoBond(
                StereoBondDelta.Add(id=5, site=3, ligands=ligands, attributes=attributes)
            ):
                assert ligands == [StereoLigand(4, StereoLigandKind.Atom)]
                assert attributes.configuration == StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans,
                    StereoCoset.Lit(0),
                )
                seen.append("stereo-bond")
            case Delta.Constraint(ConstraintDelta.Add(constraint=constraint)):
                assert constraint == Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
                seen.append("constraint")
            case _:
                raise AssertionError("delta did not match a complete family arm")

    assert seen == [
        "atom",
        "bond",
        "dative-bond",
        "aromatic-system",
        "multicenter-bond",
        "noncovalent-bond",
        "stereo-atom",
        "stereo-bond",
        "constraint",
    ]

