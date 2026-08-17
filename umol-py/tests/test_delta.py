import pytest

from umol import (
    AromaticSystemForm,
    AromaticSystemConstraintForm,
    AromaticSystemDelta,
    AromaticSystemFieldChange,
    AtomForm,
    AtomConstraintForm,
    AtomDelta,
    AtomFieldChange,
    BondForm,
    BondConstraintForm,
    BondDelta,
    BondFieldChange,
    BooleanForm,
    Constraint,
    ConstraintDelta,
    ContradictionError,
    DativeBondForm,
    DativeBondConstraintForm,
    DativeBondDelta,
    DativeBondFieldChange,
    Delta,
    Deltas,
    Element,
    ElementForm,
    ElectronCountsForm,
    IsotopeMassForm,
    Molecule,
    MoleculeConstraint,
    MulticenterBondForm,
    MulticenterBondConstraintForm,
    MulticenterBondDelta,
    MulticenterBondFieldChange,
    NoncovalentBondForm,
    NoncovalentBondConstraintForm,
    NoncovalentBondDelta,
    NoncovalentBondFieldChange,
    NoncovalentBondKind,
    NoncovalentBondKindForm,
    StereoAtomForm,
    StereoAtomConstraintForm,
    StereoAtomDelta,
    StereoAtomFieldChange,
    StereoBondForm,
    StereoBondConstraintForm,
    StereoBondDelta,
    StereoBondFieldChange,
    StereoConfigurationForm,
    StereoCoset,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereogenicityForm,
    UnpairedElectronsForm,
    NumForm,
)


@pytest.mark.parametrize(
    (
        "source",
        "make_entity_delta",
        "wrap_delta",
        "field",
        "source_update",
        "copy_update",
        "constraint",
    ),
    [
        pytest.param(
            AtomForm(Element("C")),
            lambda value: AtomDelta.Add(id=0, attributes=value),
            Delta.Atom,
            "charge",
            -1,
            1,
            AtomConstraintForm.Valence(NumForm.Lit(4)),
            id="atom",
        ),
        pytest.param(
            BondForm(1),
            lambda value: BondDelta.Add(id=0, atoms=(0, 1), attributes=value),
            Delta.Bond,
            "order",
            2,
            3,
            BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            id="bond",
        ),
        pytest.param(
            DativeBondForm(1),
            lambda value: DativeBondDelta.Add(
                id=0, donors=[0], acceptor=1, attributes=value
            ),
            Delta.DativeBond,
            "order",
            2,
            3,
            DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            id="dative-bond",
        ),
        pytest.param(
            AromaticSystemForm([1, 1]),
            lambda value: AromaticSystemDelta.Add(
                id=0, atoms=[0, 1], attributes=value
            ),
            Delta.AromaticSystem,
            "charge",
            -1,
            1,
            AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(2)),
            id="aromatic-system",
        ),
        pytest.param(
            MulticenterBondForm([1, 1]),
            lambda value: MulticenterBondDelta.Add(
                id=0, atoms=[0, 1], attributes=value
            ),
            Delta.MulticenterBond,
            "charge",
            -1,
            1,
            MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(2)),
            id="multicenter-bond",
        ),
        pytest.param(
            NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
            lambda value: NoncovalentBondDelta.Add(
                id=0, atoms=(0, 1), attributes=value
            ),
            Delta.NoncovalentBond,
            "kind",
            NoncovalentBondKind.Ionic,
            NoncovalentBondKind.HydrogenBond,
            NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)),
            id="noncovalent-bond",
        ),
        pytest.param(
            StereoAtomForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(0)
                )
            ),
            lambda value: StereoAtomDelta.Add(
                id=0,
                site=0,
                ligands=[StereoLigand(1, StereoLigandKind.Atom)],
                attributes=value,
            ),
            Delta.StereoAtom,
            "configuration",
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(1)
            ),
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Lit(0)
            ),
            StereoAtomConstraintForm.Stereogenicity(
                StereogenicityForm.Undetermined()
            ),
            id="stereo-atom",
        ),
        pytest.param(
            StereoBondForm(
                StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                )
            ),
            lambda value: StereoBondDelta.Add(
                id=0,
                site=0,
                ligands=[StereoLigand(1, StereoLigandKind.Atom)],
                attributes=value,
            ),
            Delta.StereoBond,
            "configuration",
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(1)
            ),
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(0)
            ),
            StereoBondConstraintForm.Stereogenicity(
                StereogenicityForm.Undetermined()
            ),
            id="stereo-bond",
        ),
    ],
)
def test_entity_delta_readonly_contract(
    source,
    make_entity_delta,
    wrap_delta,
    field,
    source_update,
    copy_update,
    constraint,
):
    expected = str(source)
    entity_delta = make_entity_delta(source)
    attributes = entity_delta.attributes

    assert source.readonly is False
    assert type(source).parse(str(source)).readonly is False
    assert attributes.readonly is True
    assert attributes is entity_delta.attributes
    with pytest.raises(AttributeError):
        source.readonly = True

    setattr(source, field, source_update)
    assert str(attributes) == expected

    with pytest.raises(TypeError):
        setattr(attributes, field, copy_update)
    with pytest.raises(TypeError):
        attributes.constraints.set(constraint)
    with pytest.raises(TypeError):
        attributes.constraints = attributes.constraints

    writable = attributes.copy()
    assert writable.readonly is False
    setattr(writable, field, copy_update)

    normalized = attributes.normalize()
    meet = attributes.meet(attributes)
    join = attributes.join(attributes)
    assert normalized.readonly is False
    assert meet is not None and meet.readonly is False
    assert join is not None and join.readonly is False

    assert entity_delta == make_entity_delta(attributes.copy())
    inverse = entity_delta.inverse()
    assert inverse is not entity_delta
    assert inverse.attributes.readonly is True
    assert inverse.inverse() == entity_delta

    wrapped = wrap_delta(entity_delta)
    deltas = Deltas([wrapped])
    deltas.append(wrapped)
    deltas.extend([wrapped])
    assert list(deltas) == [wrapped, wrapped, wrapped]


def test_atomfieldchange_fields():
    change = AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1))

    assert change.old == NumForm.Lit(0)
    assert change.new == NumForm.Lit(-1)
    assert repr(change) == (
        "AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1))"
    )
    with pytest.raises(AttributeError):
        change.old = NumForm.Lit(1)
    with pytest.raises(TypeError):
        hash(change)


def test_atomfieldchange_inverse():
    change = AtomFieldChange.UnpairedElectrons(
        old=UnpairedElectronsForm(0, 1), new=UnpairedElectronsForm(1, 2)
    )

    inverse = change.inverse()

    assert isinstance(inverse, AtomFieldChange.UnpairedElectrons)
    assert inverse.old == UnpairedElectronsForm(1, 2)
    assert inverse.new == UnpairedElectronsForm(0, 1)
    assert inverse.inverse() == change


def test_bondfieldchange_fields():
    change = BondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2))

    assert change.old == NumForm.Lit(1)
    assert change.new == NumForm.Lit(2)
    assert repr(change) == (
        "BondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2))"
    )


def test_dativebondfieldchange_fields():
    change = DativeBondFieldChange.Order(
        old=NumForm.Lit(1), new=NumForm.Lit(2)
    )

    assert change.old == NumForm.Lit(1)
    assert change.new == NumForm.Lit(2)
    assert repr(change) == (
        "DativeBondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2))"
    )


def test_aromaticsystemfieldchange_fields():
    change = AromaticSystemFieldChange.Electrons(
        old=ElectronCountsForm.Undetermined(),
        new=ElectronCountsForm.Lit([1, 1, 1]),
    )

    assert change.old == ElectronCountsForm.Undetermined()
    assert change.new == ElectronCountsForm.Lit([1, 1, 1])
    assert repr(change) == (
        "AromaticSystemFieldChange.Electrons("
        "old=ElectronCountsForm.Undetermined(), "
        "new=ElectronCountsForm.Lit([1, 1, 1]))"
    )


def test_multicenterbondfieldchange_fields():
    change = MulticenterBondFieldChange.Charge(
        old=NumForm.Lit(0), new=NumForm.Lit(1)
    )

    assert change.old == NumForm.Lit(0)
    assert change.new == NumForm.Lit(1)
    assert repr(change) == (
        "MulticenterBondFieldChange.Charge("
        "old=NumForm.Lit(0), new=NumForm.Lit(1))"
    )


def test_stereoatomfieldchange_fields():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationForm.Undetermined(),
        new=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Undetermined()
        ),
    )

    assert change.old == StereoConfigurationForm.Undetermined()
    assert change.new == StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Undetermined()
    )
    assert repr(change) == (
        "StereoAtomFieldChange.Configuration("
        "old=StereoConfigurationForm.Undetermined(), "
        "new=StereoConfigurationForm.Kinded("
        "StereoKind.Tetrahedral, StereoCoset.Undetermined()))"
    )
    with pytest.raises(AttributeError):
        change.old = StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(0)
        )
    with pytest.raises(TypeError):
        hash(change)


def test_stereoatomfieldchange_inverse():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Undetermined()
        ),
        new=StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(1)
        ),
    )

    inverse = change.inverse()

    assert isinstance(inverse, StereoAtomFieldChange.Configuration)
    assert inverse.old == StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(1)
    )
    assert inverse.new == StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Undetermined()
    )
    assert inverse != change
    assert inverse.inverse() == change


def test_stereobondfieldchange_fields():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationForm.Undetermined(),
        new=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Undetermined()
        ),
    )

    assert change.old == StereoConfigurationForm.Undetermined()
    assert change.new == StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Undetermined()
    )
    assert repr(change) == (
        "StereoBondFieldChange.Configuration("
        "old=StereoConfigurationForm.Undetermined(), "
        "new=StereoConfigurationForm.Kinded("
        "StereoKind.CisTrans, StereoCoset.Undetermined()))"
    )
    with pytest.raises(AttributeError):
        change.old = StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Lit(0)
        )
    with pytest.raises(TypeError):
        hash(change)


def test_stereobondfieldchange_inverse():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Undetermined()
        ),
        new=StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Lit(1)
        ),
    )

    inverse = change.inverse()

    assert isinstance(inverse, StereoBondFieldChange.Configuration)
    assert inverse.old == StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(1)
    )
    assert inverse.new == StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Undetermined()
    )
    assert inverse != change
    assert inverse.inverse() == change


@pytest.mark.parametrize(
    ("change", "expected_repr"),
    [
        (
            StereoAtomFieldChange.Configuration(
                old=StereoConfigurationForm.Undetermined(),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Undetermined()
                ),
            ),
            "StereoAtomFieldChange.Configuration("
            "old=StereoConfigurationForm.Undetermined(), "
            "new=StereoConfigurationForm.Kinded("
            "StereoKind.Tetrahedral, StereoCoset.Undetermined()))",
        ),
        (
            StereoAtomFieldChange.Configuration(
                old=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Undetermined()
                ),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                ),
            ),
            "StereoAtomFieldChange.Configuration("
            "old=StereoConfigurationForm.Kinded("
            "StereoKind.Tetrahedral, StereoCoset.Undetermined()), "
            "new=StereoConfigurationForm.Kinded("
            "StereoKind.Tetrahedral, StereoCoset.Lit(1)))",
        ),
        (
            StereoBondFieldChange.Configuration(
                old=StereoConfigurationForm.Undetermined(),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Undetermined()
                ),
            ),
            "StereoBondFieldChange.Configuration("
            "old=StereoConfigurationForm.Undetermined(), "
            "new=StereoConfigurationForm.Kinded("
            "StereoKind.CisTrans, StereoCoset.Undetermined()))",
        ),
        (
            StereoBondFieldChange.Configuration(
                old=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Undetermined()
                ),
                new=StereoConfigurationForm.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(1)
                ),
            ),
            "StereoBondFieldChange.Configuration("
            "old=StereoConfigurationForm.Kinded("
            "StereoKind.CisTrans, StereoCoset.Undetermined()), "
            "new=StereoConfigurationForm.Kinded("
            "StereoKind.CisTrans, StereoCoset.Lit(1)))",
        ),
    ],
    ids=[
        "atom-geometry-unknown",
        "atom-coset-resolved",
        "bond-geometry-unknown",
        "bond-coset-resolved",
    ],
)
def test_stereofieldchange_closure(change, expected_repr):
    assert repr(change) == expected_repr
    inverse = change.inverse()
    assert type(inverse) is type(change)
    assert inverse != change
    assert inverse.inverse() == change


@pytest.mark.parametrize(
    ("change", "expected_repr"),
    [
        (
            AtomFieldChange.Element(
                old=ElementForm.Lit(Element("C")),
                new=ElementForm.Lit(Element("N")),
            ),
            "AtomFieldChange.Element(old=ElementForm.Lit(Element('C')), "
            "new=ElementForm.Lit(Element('N')))",
        ),
        (
            AtomFieldChange.IsotopeMass(
                old=IsotopeMassForm.Lit(12),
                new=IsotopeMassForm.Lit(13),
            ),
            "AtomFieldChange.IsotopeMass(old=IsotopeMassForm.Lit(12), "
            "new=IsotopeMassForm.Lit(13))",
        ),
        (
            AtomFieldChange.Charge(
                old=NumForm.Lit(0),
                new=NumForm.Lit(-1),
            ),
            "AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1))",
        ),
        (
            AtomFieldChange.ImplicitHydrogens(
                old=NumForm.Lit(3),
                new=NumForm.Lit(2),
            ),
            "AtomFieldChange.ImplicitHydrogens(old=NumForm.Lit(3), "
            "new=NumForm.Lit(2))",
        ),
        (
            AtomFieldChange.LonePairs(
                old=NumForm.Lit(1),
                new=NumForm.Lit(2),
            ),
            "AtomFieldChange.LonePairs(old=NumForm.Lit(1), new=NumForm.Lit(2))",
        ),
        (
            AtomFieldChange.UnpairedElectrons(
                old=UnpairedElectronsForm(0, 1),
                new=UnpairedElectronsForm(1, 2),
            ),
            "AtomFieldChange.UnpairedElectrons(old=UnpairedElectronsForm(NumForm.Lit(0), "
            "NumForm.Lit(1)), new=UnpairedElectronsForm(NumForm.Lit(1), "
            "NumForm.Lit(2)))",
        ),
        (
            BondFieldChange.Order(
                old=NumForm.Lit(1),
                new=NumForm.Lit(2),
            ),
            "BondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2))",
        ),
        (
            BondFieldChange.Charge(
                old=NumForm.Lit(0),
                new=NumForm.Lit(1),
            ),
            "BondFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(1))",
        ),
        (
            BondFieldChange.UnpairedElectrons(
                old=UnpairedElectronsForm(0, 1),
                new=UnpairedElectronsForm(1, 2),
            ),
            "BondFieldChange.UnpairedElectrons(old=UnpairedElectronsForm(NumForm.Lit(0), "
            "NumForm.Lit(1)), new=UnpairedElectronsForm(NumForm.Lit(1), "
            "NumForm.Lit(2)))",
        ),
        (
            DativeBondFieldChange.Order(
                old=NumForm.Lit(1),
                new=NumForm.Lit(2),
            ),
            "DativeBondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2))",
        ),
        (
            AromaticSystemFieldChange.Electrons(
                old=ElectronCountsForm.Undetermined(),
                new=ElectronCountsForm.Lit([1, 1, 1]),
            ),
            "AromaticSystemFieldChange.Electrons("
            "old=ElectronCountsForm.Undetermined(), "
            "new=ElectronCountsForm.Lit([1, 1, 1]))",
        ),
        (
            AromaticSystemFieldChange.Charge(
                old=NumForm.Lit(0),
                new=NumForm.Lit(-1),
            ),
            "AromaticSystemFieldChange.Charge(old=NumForm.Lit(0), "
            "new=NumForm.Lit(-1))",
        ),
        (
            AromaticSystemFieldChange.UnpairedElectrons(
                old=UnpairedElectronsForm(0, 1),
                new=UnpairedElectronsForm(1, 2),
            ),
            "AromaticSystemFieldChange.UnpairedElectrons(old=UnpairedElectronsForm("
            "NumForm.Lit(0), NumForm.Lit(1)), new=UnpairedElectronsForm("
            "NumForm.Lit(1), NumForm.Lit(2)))",
        ),
        (
            MulticenterBondFieldChange.Electrons(
                old=ElectronCountsForm.Lit([1, 0, 1]),
                new=ElectronCountsForm.Lit([2, 0, 1]),
            ),
            "MulticenterBondFieldChange.Electrons("
            "old=ElectronCountsForm.Lit([1, 0, 1]), "
            "new=ElectronCountsForm.Lit([2, 0, 1]))",
        ),
        (
            MulticenterBondFieldChange.Charge(
                old=NumForm.Lit(0),
                new=NumForm.Lit(1),
            ),
            "MulticenterBondFieldChange.Charge(old=NumForm.Lit(0), "
            "new=NumForm.Lit(1))",
        ),
        (
            MulticenterBondFieldChange.UnpairedElectrons(
                old=UnpairedElectronsForm(0, 1),
                new=UnpairedElectronsForm(2, 3),
            ),
            "MulticenterBondFieldChange.UnpairedElectrons(old=UnpairedElectronsForm("
            "NumForm.Lit(0), NumForm.Lit(1)), new=UnpairedElectronsForm("
            "NumForm.Lit(2), NumForm.Lit(3)))",
        ),
        (
            NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindForm.Undetermined(),
                new=NoncovalentBondKindForm.Lit(
                    NoncovalentBondKind.HydrogenBond
                ),
            ),
            "NoncovalentBondFieldChange.Kind("
            "old=NoncovalentBondKindForm.Undetermined(), "
            "new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond))",
        ),
    ],
    ids=[
        "atom-element",
        "atom-isotope-mass",
        "atom-charge",
        "atom-implicit-hydrogens",
        "atom-lone-pairs",
        "atom-unpaired-electrons",
        "bond-order",
        "bond-charge",
        "bond-unpaired-electrons",
        "dative-order",
        "aromatic-electrons",
        "aromatic-charge",
        "aromatic-unpaired-electrons",
        "multicenter-electrons",
        "multicenter-charge",
        "multicenter-unpaired-electrons",
        "noncovalent-kind",
    ],
)
def test_fieldchange_closure(change, expected_repr):
    assert repr(change) == expected_repr
    inverse = change.inverse()
    assert type(inverse) is type(change)
    assert inverse != change
    assert inverse.inverse() == change


def test_atomdelta_fields():
    source = AtomForm(Element("C"))
    delta = AtomDelta.Add(id=3, attributes=source)

    source.charge = -1

    assert delta.id == 3
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.attributes.charge == NumForm.Undetermined()
    assert repr(delta) == "AtomDelta.Add(id=3, attributes=AtomForm.parse('C'))"
    with pytest.raises(TypeError):
        delta.attributes.charge = 1
    with pytest.raises(AttributeError):
        delta.id = 4
    with pytest.raises(TypeError):
        hash(delta)


def test_bonddelta_fields():
    source = BondForm(1)
    delta = BondDelta.Add(id=2, atoms=(5, 1), attributes=source)

    source.order = 2

    assert delta.id == 2
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.atoms == (5, 1)
    assert isinstance(delta.atoms, tuple)
    assert delta.attributes.order == NumForm.Lit(1)
    assert repr(delta) == (
        "BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm.parse('1'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.order = 3
    with pytest.raises(AttributeError):
        delta.atoms = (1, 5)
    with pytest.raises(TypeError):
        hash(delta)


def test_dativebonddelta_fields():
    source = DativeBondForm(1)
    delta = DativeBondDelta.Add(
        id=1,
        donors=[4, 2, 4],
        acceptor=3,
        attributes=source,
    )

    source.order = 2

    assert delta.id == 1
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.donors == [4, 2, 4]
    assert isinstance(delta.donors, list)
    assert delta.acceptor == 3
    assert delta.attributes.order == NumForm.Lit(1)
    assert repr(delta) == (
        "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, "
        "attributes=DativeBondForm.parse('1'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.order = 3
    with pytest.raises(AttributeError):
        delta.donors = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_aromaticsystemdelta_fields():
    source = AromaticSystemForm([1, 1, 1])
    delta = AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], attributes=source)

    source.electrons = [2, 0, 1]

    assert delta.id == 2
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.atoms == [4, 2, 4]
    assert isinstance(delta.atoms, list)
    assert delta.attributes.electrons == ElectronCountsForm.Lit([1, 1, 1])
    assert repr(delta) == (
        "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], "
        "attributes=AromaticSystemForm.parse('[1,1,1]'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.charge = -1
    with pytest.raises(AttributeError):
        delta.atoms = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_multicenterbonddelta_fields():
    source = MulticenterBondForm([1, 1, 1])
    delta = MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], attributes=source)

    source.electrons = [2, 0, 1]

    assert delta.id == 3
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.atoms == [4, 2, 4]
    assert isinstance(delta.atoms, list)
    assert delta.attributes.electrons == ElectronCountsForm.Lit([1, 1, 1])
    assert repr(delta) == (
        "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], "
        "attributes=MulticenterBondForm.parse('[1,1,1]'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.charge = -1
    with pytest.raises(AttributeError):
        delta.atoms = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_noncovalentbonddelta_fields():
    source = NoncovalentBondForm(NoncovalentBondKind.HydrogenBond)
    delta = NoncovalentBondDelta.Add(id=4, atoms=(5, 2), attributes=source)

    source.kind = NoncovalentBondKind.Ionic

    assert delta.id == 4
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.atoms == (5, 2)
    assert isinstance(delta.atoms, tuple)
    assert delta.attributes.kind == NoncovalentBondKindForm.Lit(
        NoncovalentBondKind.HydrogenBond
    )
    assert repr(delta) == (
        "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), "
        "attributes=NoncovalentBondForm.parse('Hbd'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.kind = NoncovalentBondKind.Ionic
    with pytest.raises(AttributeError):
        delta.atoms = (2, 5)
    with pytest.raises(TypeError):
        hash(delta)


def test_stereoatomdelta_fields():
    source = StereoAtomForm(
        StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(0)
        )
    )
    delta = StereoAtomDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
            StereoLigand(4, StereoLigandKind.Atom),
        ],
        attributes=source,
    )

    source.configuration = StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(1)
    )

    assert delta.id == 5
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.site == 3
    assert delta.ligands == [
        StereoLigand(4, StereoLigandKind.Atom),
        StereoLigand(2, StereoLigandKind.LonePair),
        StereoLigand(4, StereoLigandKind.Atom),
    ]
    assert isinstance(delta.ligands, list)
    assert delta.attributes.configuration == StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(0)
    )
    assert repr(delta) == (
        "StereoAtomDelta.Add(id=5, site=3, ligands=["
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), "
        "StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair), "
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], "
        "attributes=StereoAtomForm.parse('Th0'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.configuration = StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(1)
        )
    with pytest.raises(AttributeError):
        delta.ligands = []
    with pytest.raises(TypeError):
        hash(delta)


def test_stereoatomdelta_modifyconstraint_kind_none():
    delta = StereoAtomDelta.ModifyConstraint(
        id=5,
        kind=None,
        old=StereoAtomConstraintForm.Stereogenicity(
            StereogenicityForm.Undetermined()
        ),
        new=None,
    )

    assert delta.kind is None
    inverse = delta.inverse()
    assert inverse.kind is None
    assert inverse.old is None
    assert inverse.new == StereoAtomConstraintForm.Stereogenicity(
        StereogenicityForm.Undetermined()
    )
    assert inverse.inverse() == delta


def test_stereobonddelta_fields():
    source = StereoBondForm(
        StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Lit(0)
        )
    )
    delta = StereoBondDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
            StereoLigand(4, StereoLigandKind.Atom),
        ],
        attributes=source,
    )

    source.configuration = StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(1)
    )

    assert delta.id == 5
    assert source.readonly is False
    assert delta.attributes.readonly is True
    assert delta.site == 3
    assert delta.ligands == [
        StereoLigand(4, StereoLigandKind.Atom),
        StereoLigand(2, StereoLigandKind.LonePair),
        StereoLigand(4, StereoLigandKind.Atom),
    ]
    assert isinstance(delta.ligands, list)
    assert delta.attributes.configuration == StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(0)
    )
    assert repr(delta) == (
        "StereoBondDelta.Add(id=5, site=3, ligands=["
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), "
        "StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair), "
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], "
        "attributes=StereoBondForm.parse('Ct0'))"
    )
    with pytest.raises(TypeError):
        delta.attributes.configuration = StereoConfigurationForm.Kinded(
            StereoKind.CisTrans, StereoCoset.Lit(1)
        )
    with pytest.raises(AttributeError):
        delta.ligands = []
    with pytest.raises(TypeError):
        hash(delta)


def test_stereobonddelta_modifyconstraint_kind_none():
    delta = StereoBondDelta.ModifyConstraint(
        id=5,
        kind=None,
        old=StereoBondConstraintForm.Stereogenicity(
            StereogenicityForm.Undetermined()
        ),
        new=None,
    )

    assert delta.kind is None
    inverse = delta.inverse()
    assert inverse.kind is None
    assert inverse.old is None
    assert inverse.new == StereoBondConstraintForm.Stereogenicity(
        StereogenicityForm.Undetermined()
    )
    assert inverse.inverse() == delta



@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type"),
    [
        (
            AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))),
            "AtomDelta.Add(id=3, attributes=AtomForm.parse('C'))",
            AtomDelta.Remove,
        ),
        (
            AtomDelta.Remove(id=3, attributes=AtomForm(Element("N"))),
            "AtomDelta.Remove(id=3, attributes=AtomForm.parse('N'))",
            AtomDelta.Add,
        ),
        (
            AtomDelta.ModifyField(
                id=3,
                change=AtomFieldChange.Charge(
                    old=NumForm.Lit(0), new=NumForm.Lit(-1)
                ),
            ),
            "AtomDelta.ModifyField(id=3, "
            "change=AtomFieldChange.Charge("
            "old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
            AtomDelta.ModifyField,
        ),
        (
            AtomDelta.ModifyConstraint(
                id=3,
                old=None,
                new=AtomConstraintForm.Valence(NumForm.Lit(4)),
            ),
            "AtomDelta.ModifyConstraint(id=3, old=None, "
            "new=AtomConstraintForm.Valence(NumForm.Lit(4)))",
            AtomDelta.ModifyConstraint,
        ),
        (
            BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm(1)),
            "BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm.parse('1'))",
            BondDelta.Remove,
        ),
        (
            BondDelta.Remove(id=2, atoms=(1, 5), attributes=BondForm(2)),
            "BondDelta.Remove(id=2, atoms=(1, 5), attributes=BondForm.parse('2'))",
            BondDelta.Add,
        ),
        (
            BondDelta.ModifyField(
                id=2,
                change=BondFieldChange.Order(
                    old=NumForm.Lit(1), new=NumForm.Lit(2)
                ),
            ),
            "BondDelta.ModifyField(id=2, "
            "change=BondFieldChange.Order("
            "old=NumForm.Lit(1), new=NumForm.Lit(2)))",
            BondDelta.ModifyField,
        ),
        (
            BondDelta.ModifyConstraint(
                id=2,
                old=None,
                new=BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            ),
            "BondDelta.ModifyConstraint(id=2, old=None, "
            "new=BondConstraintForm.Aromatic(BooleanForm.Lit(True)))",
            BondDelta.ModifyConstraint,
        ),
    ],
    ids=[
        "atom-add",
        "atom-remove",
        "atom-modify-field",
        "atom-modify-constraint",
        "bond-add",
        "bond-remove",
        "bond-modify-field",
        "bond-modify-constraint",
    ],
)
def test_entitydelta_closure(delta, expected_repr, inverse_type):
    assert repr(delta) == expected_repr
    inverse = delta.inverse()
    assert type(inverse) is inverse_type
    assert inverse != delta
    assert inverse.inverse() == delta


@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type"),
    [
        (
            DativeBondDelta.Add(
                id=1, donors=[4, 2, 4], acceptor=3, attributes=DativeBondForm(1)
            ),
            "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, "
            "attributes=DativeBondForm.parse('1'))",
            DativeBondDelta.Remove,
        ),
        (
            DativeBondDelta.Remove(
                id=1, donors=[2, 4, 2], acceptor=3, attributes=DativeBondForm(2)
            ),
            "DativeBondDelta.Remove(id=1, donors=[2, 4, 2], acceptor=3, "
            "attributes=DativeBondForm.parse('2'))",
            DativeBondDelta.Add,
        ),
        (
            DativeBondDelta.ModifyField(
                id=1,
                change=DativeBondFieldChange.Order(
                    old=NumForm.Lit(1), new=NumForm.Lit(2)
                ),
            ),
            "DativeBondDelta.ModifyField(id=1, "
            "change=DativeBondFieldChange.Order("
            "old=NumForm.Lit(1), new=NumForm.Lit(2)))",
            DativeBondDelta.ModifyField,
        ),
        (
            DativeBondDelta.ModifyConstraint(
                id=1,
                old=None,
                new=DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
            ),
            "DativeBondDelta.ModifyConstraint(id=1, old=None, "
            "new=DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)))",
            DativeBondDelta.ModifyConstraint,
        ),
        (
            AromaticSystemDelta.Add(
                id=2, atoms=[4, 2, 4], attributes=AromaticSystemForm([1, 1, 1])
            ),
            "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], "
            "attributes=AromaticSystemForm.parse('[1,1,1]'))",
            AromaticSystemDelta.Remove,
        ),
        (
            AromaticSystemDelta.Remove(
                id=2, atoms=[2, 4, 2], attributes=AromaticSystemForm([2, 0, 1])
            ),
            "AromaticSystemDelta.Remove(id=2, atoms=[2, 4, 2], "
            "attributes=AromaticSystemForm.parse('[2,0,1]'))",
            AromaticSystemDelta.Add,
        ),
        (
            AromaticSystemDelta.ModifyField(
                id=2,
                change=AromaticSystemFieldChange.Charge(
                    old=NumForm.Lit(0), new=NumForm.Lit(-1)
                ),
            ),
            "AromaticSystemDelta.ModifyField(id=2, "
            "change=AromaticSystemFieldChange.Charge("
            "old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
            AromaticSystemDelta.ModifyField,
        ),
        (
            AromaticSystemDelta.ModifyConstraint(
                id=2,
                old=None,
                new=AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6)),
            ),
            "AromaticSystemDelta.ModifyConstraint(id=2, old=None, "
            "new=AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6)))",
            AromaticSystemDelta.ModifyConstraint,
        ),
        (
            MulticenterBondDelta.Add(
                id=3, atoms=[4, 2, 4], attributes=MulticenterBondForm([1, 1, 1])
            ),
            "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], "
            "attributes=MulticenterBondForm.parse('[1,1,1]'))",
            MulticenterBondDelta.Remove,
        ),
        (
            MulticenterBondDelta.Remove(
                id=3, atoms=[2, 4, 2], attributes=MulticenterBondForm([2, 0, 1])
            ),
            "MulticenterBondDelta.Remove(id=3, atoms=[2, 4, 2], "
            "attributes=MulticenterBondForm.parse('[2,0,1]'))",
            MulticenterBondDelta.Add,
        ),
        (
            MulticenterBondDelta.ModifyField(
                id=3,
                change=MulticenterBondFieldChange.Charge(
                    old=NumForm.Lit(0), new=NumForm.Lit(-1)
                ),
            ),
            "MulticenterBondDelta.ModifyField(id=3, "
            "change=MulticenterBondFieldChange.Charge("
            "old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
            MulticenterBondDelta.ModifyField,
        ),
        (
            MulticenterBondDelta.ModifyConstraint(
                id=3,
                old=None,
                new=MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6)),
            ),
            "MulticenterBondDelta.ModifyConstraint(id=3, old=None, "
            "new=MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6)))",
            MulticenterBondDelta.ModifyConstraint,
        ),
        (
            NoncovalentBondDelta.Add(
                id=4,
                atoms=(5, 2),
                attributes=NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
            ),
            "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), "
            "attributes=NoncovalentBondForm.parse('Hbd'))",
            NoncovalentBondDelta.Remove,
        ),
        (
            NoncovalentBondDelta.Remove(
                id=4,
                atoms=(2, 5),
                attributes=NoncovalentBondForm(NoncovalentBondKind.Ionic),
            ),
            "NoncovalentBondDelta.Remove(id=4, atoms=(2, 5), "
            "attributes=NoncovalentBondForm.parse('Ion'))",
            NoncovalentBondDelta.Add,
        ),
        (
            NoncovalentBondDelta.ModifyField(
                id=4,
                change=NoncovalentBondFieldChange.Kind(
                    old=NoncovalentBondKindForm.Undetermined(),
                    new=NoncovalentBondKindForm.Lit(
                        NoncovalentBondKind.HydrogenBond
                    ),
                ),
            ),
            "NoncovalentBondDelta.ModifyField(id=4, "
            "change=NoncovalentBondFieldChange.Kind("
            "old=NoncovalentBondKindForm.Undetermined(), "
            "new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond)))",
            NoncovalentBondDelta.ModifyField,
        ),
        (
            NoncovalentBondDelta.ModifyConstraint(
                id=4,
                old=None,
                new=NoncovalentBondConstraintForm.Intramolecular(
                    BooleanForm.Lit(True)
                ),
            ),
            "NoncovalentBondDelta.ModifyConstraint(id=4, old=None, "
            "new=NoncovalentBondConstraintForm.Intramolecular("
            "BooleanForm.Lit(True)))",
            NoncovalentBondDelta.ModifyConstraint,
        ),
    ],
    ids=[
        "dative-add",
        "dative-remove",
        "dative-modify-field",
        "dative-modify-constraint",
        "aromatic-add",
        "aromatic-remove",
        "aromatic-modify-field",
        "aromatic-modify-constraint",
        "multicenter-add",
        "multicenter-remove",
        "multicenter-modify-field",
        "multicenter-modify-constraint",
        "noncovalent-add",
        "noncovalent-remove",
        "noncovalent-modify-field",
        "noncovalent-modify-constraint",
    ],
)
def test_overlaydelta_closure(delta, expected_repr, inverse_type):
    assert repr(delta) == expected_repr
    inverse = delta.inverse()
    assert type(inverse) is inverse_type
    assert inverse != delta
    assert inverse.inverse() == delta


@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type", "self_inverse"),
    [
        (
            StereoAtomDelta.Add(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoAtomForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.Tetrahedral, StereoCoset.Lit(0)
                    )
                ),
            ),
            "StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoAtomForm.parse('Th0'))",
            StereoAtomDelta.Remove,
            False,
        ),
        (
            StereoAtomDelta.Remove(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoAtomForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.Tetrahedral, StereoCoset.Lit(0)
                    )
                ),
            ),
            "StereoAtomDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoAtomForm.parse('Th0'))",
            StereoAtomDelta.Add,
            False,
        ),
        (
            StereoAtomDelta.ModifyField(
                id=5,
                change=StereoAtomFieldChange.Configuration(
                    old=StereoConfigurationForm.Undetermined(),
                    new=StereoConfigurationForm.Kinded(
                        StereoKind.Tetrahedral, StereoCoset.Lit(0)
                    ),
                ),
            ),
            "StereoAtomDelta.ModifyField(id=5, change=StereoAtomFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0))))",
            StereoAtomDelta.ModifyField,
            False,
        ),
        (
            StereoAtomDelta.ModifyConstraint(
                id=5,
                kind=StereoKind.Tetrahedral,
                old=None,
                new=StereoAtomConstraintForm.Stereogenicity(
                    StereogenicityForm.Undetermined()
                ),
            ),
            "StereoAtomDelta.ModifyConstraint(id=5, kind=StereoKind.Tetrahedral, old=None, new=StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()))",
            StereoAtomDelta.ModifyConstraint,
            False,
        ),
        (
            StereoBondDelta.Add(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoBondForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.CisTrans, StereoCoset.Lit(0)
                    )
                ),
            ),
            "StereoBondDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoBondForm.parse('Ct0'))",
            StereoBondDelta.Remove,
            False,
        ),
        (
            StereoBondDelta.Remove(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                attributes=StereoBondForm(
                    StereoConfigurationForm.Kinded(
                        StereoKind.CisTrans, StereoCoset.Lit(0)
                    )
                ),
            ),
            "StereoBondDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoBondForm.parse('Ct0'))",
            StereoBondDelta.Add,
            False,
        ),
        (
            StereoBondDelta.ModifyField(
                id=5,
                change=StereoBondFieldChange.Configuration(
                    old=StereoConfigurationForm.Undetermined(),
                    new=StereoConfigurationForm.Kinded(
                        StereoKind.CisTrans, StereoCoset.Lit(0)
                    ),
                ),
            ),
            "StereoBondDelta.ModifyField(id=5, change=StereoBondFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(0))))",
            StereoBondDelta.ModifyField,
            False,
        ),
        (
            StereoBondDelta.ModifyConstraint(
                id=5,
                kind=StereoKind.CisTrans,
                old=None,
                new=StereoBondConstraintForm.Stereogenicity(
                    StereogenicityForm.Undetermined()
                ),
            ),
            "StereoBondDelta.ModifyConstraint(id=5, kind=StereoKind.CisTrans, old=None, new=StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()))",
            StereoBondDelta.ModifyConstraint,
            False,
        ),
    ],
    ids=[
        "atom-add",
        "atom-remove",
        "atom-modify-field",
        "atom-modify-constraint",
        "bond-add",
        "bond-remove",
        "bond-modify-field",
        "bond-modify-constraint",
    ],
)
def test_stereodelta_closure(delta, expected_repr, inverse_type, self_inverse):
    assert repr(delta) == expected_repr
    inverse = delta.inverse()
    assert type(inverse) is inverse_type
    assert (inverse == delta) is self_inverse
    assert inverse.inverse() == delta


def test_constraintdelta_fields():
    source = Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2)))
    delta = ConstraintDelta.Add(constraint=source)

    assert delta.constraint == source
    assert delta.constraint is source
    assert delta.constraint is delta.constraint
    assert repr(delta) == (
        "ConstraintDelta.Add(constraint=Constraint.Atom(3, "
        "AtomConstraintForm.Degree(NumForm.Lit(2))))"
    )
    with pytest.raises(AttributeError):
        delta.constraint = Constraint.Or([])
    with pytest.raises(TypeError):
        hash(delta)


@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type"),
    [
        (
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
            ),
            "ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2))))",
            ConstraintDelta.Remove,
        ),
        (
            ConstraintDelta.Remove(
                constraint=Constraint.And(
                    [
                        Constraint.Atom(
                            7,
                            AtomConstraintForm.Valence(NumForm.Lit(4)),
                        ),
                        Constraint.Not(Constraint.Or([])),
                    ]
                )
            ),
            "ConstraintDelta.Remove(constraint=Constraint.And([Constraint.Atom(7, AtomConstraintForm.Valence(NumForm.Lit(4))), Constraint.Not(Constraint.Or([]))]))",
            ConstraintDelta.Add,
        ),
    ],
    ids=["add-leaf", "remove-recursive"],
)
def test_constraintdelta_closure(delta, expected_repr, inverse_type):
    assert repr(delta) == expected_repr
    inverse = delta.inverse()
    assert type(inverse) is inverse_type
    assert inverse != delta
    assert inverse.inverse() == delta


def test_delta_fields():
    child = AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))
    delta = Delta.Atom(child)

    assert delta._0 is child
    assert delta == Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))))
    assert repr(delta) == "Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C')))"

    with pytest.raises(TypeError):
        child.attributes.charge = 1
    assert delta._0.attributes.charge == NumForm.Undetermined()

    with pytest.raises(AttributeError):
        delta._0 = AtomDelta.Remove(id=3, attributes=AtomForm(Element("C")))
    with pytest.raises(TypeError):
        hash(delta)


@pytest.mark.parametrize(
    ("delta", "expected_repr", "outer_type", "inverse_type"),
    [
        (
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
            "Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C')))",
            Delta.Atom,
            AtomDelta.Remove,
        ),
        (
            Delta.Bond(BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm(1))),
            "Delta.Bond(BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm.parse('1')))",
            Delta.Bond,
            BondDelta.Remove,
        ),
        (
            Delta.DativeBond(
                DativeBondDelta.Add(
                    id=1,
                    donors=[4, 2],
                    acceptor=3,
                    attributes=DativeBondForm(1),
                )
            ),
            "Delta.DativeBond(DativeBondDelta.Add(id=1, donors=[4, 2], acceptor=3, attributes=DativeBondForm.parse('1')))",
            Delta.DativeBond,
            DativeBondDelta.Remove,
        ),
        (
            Delta.AromaticSystem(
                AromaticSystemDelta.Add(
                    id=2,
                    atoms=[4, 2],
                    attributes=AromaticSystemForm([1, 1]),
                )
            ),
            "Delta.AromaticSystem(AromaticSystemDelta.Add(id=2, atoms=[4, 2], attributes=AromaticSystemForm.parse('[1,1]')))",
            Delta.AromaticSystem,
            AromaticSystemDelta.Remove,
        ),
        (
            Delta.MulticenterBond(
                MulticenterBondDelta.Add(
                    id=3,
                    atoms=[4, 2],
                    attributes=MulticenterBondForm([1, 1]),
                )
            ),
            "Delta.MulticenterBond(MulticenterBondDelta.Add(id=3, atoms=[4, 2], attributes=MulticenterBondForm.parse('[1,1]')))",
            Delta.MulticenterBond,
            MulticenterBondDelta.Remove,
        ),
        (
            Delta.NoncovalentBond(
                NoncovalentBondDelta.Add(
                    id=4,
                    atoms=(5, 2),
                    attributes=NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
                )
            ),
            "Delta.NoncovalentBond(NoncovalentBondDelta.Add(id=4, atoms=(5, 2), attributes=NoncovalentBondForm.parse('Hbd')))",
            Delta.NoncovalentBond,
            NoncovalentBondDelta.Remove,
        ),
        (
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
            "Delta.StereoAtom(StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoAtomForm.parse('Th0')))",
            Delta.StereoAtom,
            StereoAtomDelta.Remove,
        ),
        (
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
            "Delta.StereoBond(StereoBondDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoBondForm.parse('Ct0')))",
            Delta.StereoBond,
            StereoBondDelta.Remove,
        ),
        (
            Delta.Constraint(
                ConstraintDelta.Add(
                    constraint=Constraint.Atom(
                        3,
                        AtomConstraintForm.Degree(NumForm.Lit(2)),
                    )
                )
            ),
            "Delta.Constraint(ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2)))))",
            Delta.Constraint,
            ConstraintDelta.Remove,
        ),
    ],
    ids=[
        "atom",
        "bond",
        "dative-bond",
        "aromatic-system",
        "multicenter-bond",
        "noncovalent-bond",
        "stereo-atom",
        "stereo-bond",
        "constraint",
    ],
)
def test_delta_closure(delta, expected_repr, outer_type, inverse_type):
    assert repr(delta) == expected_repr
    assert type(delta) is outer_type
    inverse = delta.inverse()
    assert type(inverse) is outer_type
    assert type(inverse._0) is inverse_type
    assert inverse != delta
    assert inverse.inverse() == delta


def test_deltas_sequence():
    source = Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))))
    deltas = Deltas([source, source])

    with pytest.raises(TypeError):
        source._0.attributes.charge = 1

    assert Deltas() == Deltas([])
    assert len(deltas) == 2
    assert bool(deltas)
    assert not Deltas()
    assert deltas == Deltas(
        [
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        ]
    )
    assert repr(deltas) == (
        "Deltas([Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C'))), "
        "Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C')))])"
    )
    assert deltas[0]._0.attributes.charge == NumForm.Undetermined()
    with pytest.raises(TypeError):
        hash(deltas)


def test_deltas_append():
    source = Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))))
    deltas = Deltas()

    assert deltas.append(source) is None
    deltas.append(source)
    with pytest.raises(TypeError):
        source._0.attributes.charge = 1

    assert list(deltas) == [
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
    ]


def test_deltas_extend():
    target = Deltas(
        [Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C"))))]
    )
    container_source = Deltas(
        [
            Delta.Constraint(
                ConstraintDelta.Add(
                    constraint=Constraint.Atom(
                        3,
                        AtomConstraintForm.Degree(NumForm.Lit(2)),
                    )
                )
            )
        ]
    )
    iterable_source = [
        Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
        Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
    ]

    assert target.extend(container_source) is None
    target.extend(iterable_source)
    container_source.append(
        Delta.Atom(AtomDelta.Add(id=9, attributes=AtomForm(Element("O"))))
    )
    with pytest.raises(TypeError):
        iterable_source[0]._0.attributes.charge = 1

    assert list(target) == [
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Constraint(
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
            )
        ),
        Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
        Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
    ]


def test_deltas_extend_self():
    deltas = Deltas(
        [
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
            Delta.Constraint(
                ConstraintDelta.Add(
                    constraint=Constraint.Atom(
                        3,
                        AtomConstraintForm.Degree(NumForm.Lit(2)),
                    )
                )
            ),
        ]
    )

    deltas.extend(deltas)

    assert list(deltas) == [
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Constraint(
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
            )
        ),
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Constraint(
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintForm.Degree(NumForm.Lit(2)),
                )
            )
        ),
    ]


def test_deltas_getitem():
    deltas = Deltas(
        [
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
            Delta.Constraint(
                ConstraintDelta.Add(
                    constraint=Constraint.Atom(
                        3,
                        AtomConstraintForm.Degree(NumForm.Lit(2)),
                    )
                )
            ),
        ]
    )

    first = deltas[0]
    assert type(first) is Delta.Atom
    assert first._0 == AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))
    assert type(deltas[-1]) is Delta.Constraint

    with pytest.raises(TypeError):
        first._0.attributes.charge = 1
    assert deltas[0]._0.attributes.charge == NumForm.Undetermined()

    with pytest.raises(IndexError, match="delta index out of range"):
        deltas[2]
    with pytest.raises(IndexError, match="delta index out of range"):
        deltas[-3]


def test_deltas_iter():
    deltas = Deltas(
        [
            Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
            Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
        ]
    )

    entries = list(deltas)

    assert entries == [
        Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("C")))),
        Delta.Atom(AtomDelta.Add(id=4, attributes=AtomForm(Element("N")))),
    ]
    with pytest.raises(TypeError):
        entries[0]._0.attributes.charge = -1
    assert deltas[0]._0.attributes.charge == NumForm.Undetermined()


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        pytest.param(
            Deltas(
                [
                    Delta.Atom(
                        AtomDelta.ModifyField(
                            id=0,
                            change=AtomFieldChange.Charge(
                                old=NumForm.Lit(0), new=NumForm.Lit(1)
                            ),
                        )
                    ),
                    Delta.Atom(
                        AtomDelta.ModifyField(
                            id=0,
                            change=AtomFieldChange.Charge(
                                old=NumForm.Lit(1), new=NumForm.Lit(2)
                            ),
                        )
                    ),
                ]
            ),
            Deltas(
                [
                    Delta.Atom(
                        AtomDelta.ModifyField(
                            id=0,
                            change=AtomFieldChange.Charge(
                                old=NumForm.Lit(0), new=NumForm.Lit(2)
                            ),
                        )
                    )
                ]
            ),
            id="field-fusion",
        ),
        pytest.param(
            Deltas(
                [
                    Delta.Atom(AtomDelta.Add(id=0, attributes=AtomForm(Element("C")))),
                    Delta.Atom(AtomDelta.Remove(id=0, attributes=AtomForm(Element("C")))),
                ]
            ),
            Deltas(),
            id="add-remove-cancellation",
        ),
        pytest.param(
            Deltas(
                [
                    Delta.Bond(
                        BondDelta.ModifyField(
                            id=0,
                            change=BondFieldChange.Order(
                                old=NumForm.Lit(1), new=NumForm.Lit(2)
                            ),
                        )
                    ),
                    Delta.Atom(
                        AtomDelta.ModifyField(
                            id=0,
                            change=AtomFieldChange.Charge(
                                old=NumForm.Lit(0), new=NumForm.Lit(1)
                            ),
                        )
                    ),
                ]
            ),
            Deltas(
                [
                    Delta.Atom(
                        AtomDelta.ModifyField(
                            id=0,
                            change=AtomFieldChange.Charge(
                                old=NumForm.Lit(0), new=NumForm.Lit(1)
                            ),
                        )
                    ),
                    Delta.Bond(
                        BondDelta.ModifyField(
                            id=0,
                            change=BondFieldChange.Order(
                                old=NumForm.Lit(1), new=NumForm.Lit(2)
                            ),
                        )
                    ),
                ]
            ),
            id="family-order",
        ),
        pytest.param(
            Deltas(
                [
                    Delta.Constraint(
                        ConstraintDelta.Add(
                            constraint=Constraint.Atom(
                                3,
                                AtomConstraintForm.Degree(NumForm.Lit(2)),
                            )
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
                    Delta.Constraint(
                        ConstraintDelta.Remove(
                            constraint=Constraint.Atom(
                                3,
                                AtomConstraintForm.Degree(NumForm.Lit(2)),
                            )
                        )
                    ),
                ]
            ),
            Deltas(),
            id="constraint-set-semantics",
        ),
    ],
)
def test_deltas_normalize(source, expected):
    snapshot = Deltas(source)

    normalized = source.normalize()

    assert normalized == expected
    assert normalized is not source
    assert source == snapshot
    assert normalized.normalize() == normalized


def test_deltas_normalize_error():
    source = Deltas(
        [
            Delta.Atom(
                AtomDelta.ModifyField(
                    id=0,
                    change=AtomFieldChange.Charge(
                        old=NumForm.Lit(0), new=NumForm.Lit(1)
                    ),
                )
            ),
            Delta.Atom(
                AtomDelta.ModifyField(
                    id=0,
                    change=AtomFieldChange.Charge(
                        old=NumForm.Lit(2), new=NumForm.Lit(3)
                    ),
                )
            ),
        ]
    )
    snapshot = Deltas(source)

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        source.normalize()

    assert source == snapshot


def test_deltas_equiv():
    lhs = Deltas(
        [
            Delta.Bond(
                BondDelta.ModifyField(
                    id=0,
                    change=BondFieldChange.Order(
                        old=NumForm.Lit(1), new=NumForm.Lit(2)
                    ),
                )
            ),
            Delta.Atom(
                AtomDelta.ModifyField(
                    id=0,
                    change=AtomFieldChange.Charge(
                        old=NumForm.Lit(0), new=NumForm.Lit(1)
                    ),
                )
            ),
        ]
    )
    rhs = Deltas(
        [
            Delta.Atom(
                AtomDelta.ModifyField(
                    id=0,
                    change=AtomFieldChange.Charge(
                        old=NumForm.Lit(0), new=NumForm.Lit(1)
                    ),
                )
            ),
            Delta.Bond(
                BondDelta.ModifyField(
                    id=0,
                    change=BondFieldChange.Order(
                        old=NumForm.Lit(1), new=NumForm.Lit(2)
                    ),
                )
            ),
        ]
    )

    assert lhs != rhs
    assert lhs.equiv(rhs) is True
