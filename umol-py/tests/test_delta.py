import pytest

from umol import (
    AromaticSystemAst,
    AromaticSystemConstraintAst,
    AromaticSystemDelta,
    AromaticSystemFieldChange,
    AtomAst,
    AtomConstraintAst,
    AtomDelta,
    AtomFieldChange,
    BondAst,
    BondConstraintAst,
    BondDelta,
    BondFieldChange,
    BooleanAst,
    Constraint,
    ConstraintDelta,
    DativeBondAst,
    DativeBondConstraintAst,
    DativeBondDelta,
    DativeBondFieldChange,
    Element,
    ElementAst,
    ElectronCountsAst,
    IsotopeMassAst,
    MoleculeAst,
    MoleculeConstraint,
    MulticenterBondAst,
    MulticenterBondConstraintAst,
    MulticenterBondDelta,
    MulticenterBondFieldChange,
    NoncovalentBondAst,
    NoncovalentBondConstraintAst,
    NoncovalentBondDelta,
    NoncovalentBondFieldChange,
    NoncovalentBondKind,
    NoncovalentBondKindAst,
    Permutation,
    SpinStateAst,
    StereoAtomAst,
    StereoAtomConstraintAst,
    StereoAtomDelta,
    StereoAtomFieldChange,
    StereoBondAst,
    StereoBondConstraintAst,
    StereoBondDelta,
    StereoBondFieldChange,
    StereoConfigurationAst,
    StereoCosetAst,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereogenicityAst,
    SubPatternAnchor,
    ValueAst,
)


def test_atomfieldchange_fields():
    change = AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1))

    assert change.old == ValueAst.Lit(0)
    assert change.new == ValueAst.Lit(-1)
    assert repr(change) == (
        "AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1))"
    )
    with pytest.raises(AttributeError):
        change.old = ValueAst.Lit(1)
    with pytest.raises(TypeError):
        hash(change)


def test_atomfieldchange_match_scalar():
    change = AtomFieldChange.ImplicitHydrogens(
        old=ValueAst.Lit(3), new=ValueAst.Lit(2)
    )

    match change:
        case AtomFieldChange.ImplicitHydrogens(old, new):
            assert (old, new) == (ValueAst.Lit(3), ValueAst.Lit(2))
        case _:
            raise AssertionError("atom field change did not match its scalar variant")


def test_atomfieldchange_match_structured():
    change = AtomFieldChange.Element(
        old=ElementAst.Lit(Element("C")),
        new=ElementAst.Lit(Element("N")),
    )

    match change:
        case AtomFieldChange.Element(old=old, new=new):
            assert (old, new) == (
                ElementAst.Lit(Element("C")),
                ElementAst.Lit(Element("N")),
            )
        case _:
            raise AssertionError("atom field change did not match its structured variant")


def test_atomfieldchange_inverse():
    change = AtomFieldChange.Spin(old=SpinStateAst(0, 1), new=SpinStateAst(1, 2))

    inverse = change.inverse()

    assert isinstance(inverse, AtomFieldChange.Spin)
    assert inverse.old == SpinStateAst(1, 2)
    assert inverse.new == SpinStateAst(0, 1)
    assert inverse.inverse() == change


def test_bondfieldchange_fields():
    change = BondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2))

    assert change.old == ValueAst.Lit(1)
    assert change.new == ValueAst.Lit(2)
    assert repr(change) == (
        "BondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2))"
    )


def test_bondfieldchange_match():
    change = BondFieldChange.Spin(old=SpinStateAst(0, 1), new=SpinStateAst(1, 2))

    match change:
        case BondFieldChange.Spin(old, new):
            assert (old, new) == (SpinStateAst(0, 1), SpinStateAst(1, 2))
        case _:
            raise AssertionError("bond field change did not match its variant")


def test_dativebondfieldchange_fields():
    change = DativeBondFieldChange.Order(
        old=ValueAst.Lit(1), new=ValueAst.Lit(2)
    )

    assert change.old == ValueAst.Lit(1)
    assert change.new == ValueAst.Lit(2)
    assert repr(change) == (
        "DativeBondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2))"
    )


def test_dativebondfieldchange_match():
    change = DativeBondFieldChange.Order(
        old=ValueAst.Lit(1), new=ValueAst.Lit(2)
    )

    match change:
        case DativeBondFieldChange.Order(old=old, new=new):
            assert (old, new) == (ValueAst.Lit(1), ValueAst.Lit(2))
        case _:
            raise AssertionError("dative bond field change did not match its variant")

    inverse = change.inverse()
    assert isinstance(inverse, DativeBondFieldChange.Order)
    assert inverse == DativeBondFieldChange.Order(
        old=ValueAst.Lit(2), new=ValueAst.Lit(1)
    )


def test_aromaticsystemfieldchange_fields():
    change = AromaticSystemFieldChange.Electrons(
        old=ElectronCountsAst.Undetermined(),
        new=ElectronCountsAst.Lit([1, 1, 1]),
    )

    assert change.old == ElectronCountsAst.Undetermined()
    assert change.new == ElectronCountsAst.Lit([1, 1, 1])
    assert repr(change) == (
        "AromaticSystemFieldChange.Electrons("
        "old=ElectronCountsAst.Undetermined(), "
        "new=ElectronCountsAst.Lit([1, 1, 1]))"
    )


def test_aromaticsystemfieldchange_match():
    change = AromaticSystemFieldChange.Electrons(
        old=ElectronCountsAst.Lit([2, 0, 2]),
        new=ElectronCountsAst.Lit([1, 1, 1]),
    )

    match change:
        case AromaticSystemFieldChange.Electrons(old=old, new=new):
            assert (old, new) == (
                ElectronCountsAst.Lit([2, 0, 2]),
                ElectronCountsAst.Lit([1, 1, 1]),
            )
        case _:
            raise AssertionError("aromatic field change did not match its variant")


def test_multicenterbondfieldchange_fields():
    change = MulticenterBondFieldChange.Charge(
        old=ValueAst.Lit(0), new=ValueAst.Lit(1)
    )

    assert change.old == ValueAst.Lit(0)
    assert change.new == ValueAst.Lit(1)
    assert repr(change) == (
        "MulticenterBondFieldChange.Charge("
        "old=ValueAst.Lit(0), new=ValueAst.Lit(1))"
    )


def test_multicenterbondfieldchange_match():
    change = MulticenterBondFieldChange.Electrons(
        old=ElectronCountsAst.Lit([1, 0, 1]),
        new=ElectronCountsAst.Lit([2, 0, 1]),
    )

    match change:
        case MulticenterBondFieldChange.Electrons(old, new):
            assert (old, new) == (
                ElectronCountsAst.Lit([1, 0, 1]),
                ElectronCountsAst.Lit([2, 0, 1]),
            )
        case _:
            raise AssertionError("multicenter field change did not match its variant")

    inverse = change.inverse()
    assert isinstance(inverse, MulticenterBondFieldChange.Electrons)
    assert inverse.old == ElectronCountsAst.Lit([2, 0, 1])
    assert inverse.new == ElectronCountsAst.Lit([1, 0, 1])


def test_noncovalentbondfieldchange_match():
    change = NoncovalentBondFieldChange.Kind(
        old=NoncovalentBondKindAst.Undetermined(),
        new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
    )

    match change:
        case NoncovalentBondFieldChange.Kind(old=old, new=new):
            assert old == NoncovalentBondKindAst.Undetermined()
            assert new == NoncovalentBondKindAst.Lit(
                NoncovalentBondKind.HydrogenBond
            )
        case _:
            raise AssertionError("noncovalent field change did not match its variant")


def test_stereoatomfieldchange_fields():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationAst.Undetermined(),
        new=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
        ),
    )

    assert change.old == StereoConfigurationAst.Undetermined()
    assert change.new == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
    )
    assert repr(change) == (
        "StereoAtomFieldChange.Configuration("
        "old=StereoConfigurationAst.Undetermined(), "
        "new=StereoConfigurationAst.Kinded("
        "StereoKind.Tetrahedral, StereoCosetAst.Undetermined()))"
    )
    with pytest.raises(AttributeError):
        change.old = StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
        )
    with pytest.raises(TypeError):
        hash(change)


def test_stereoatomfieldchange_match_positional():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
        ),
        new=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
        ),
    )

    match change:
        case StereoAtomFieldChange.Configuration(old, new):
            assert old == StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
            )
            assert new == StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
            )
        case _:
            raise AssertionError("stereo atom field change did not match its variant")


def test_stereoatomfieldchange_match_named():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationAst.Undetermined(),
        new=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
        ),
    )

    match change:
        case StereoAtomFieldChange.Configuration(old=old, new=new):
            assert old == StereoConfigurationAst.Undetermined()
            assert new == StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
            )
        case _:
            raise AssertionError("stereo atom field change did not match its variant")


def test_stereoatomfieldchange_inverse():
    change = StereoAtomFieldChange.Configuration(
        old=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
        ),
        new=StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
        ),
    )

    inverse = change.inverse()

    assert isinstance(inverse, StereoAtomFieldChange.Configuration)
    assert inverse.old == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
    )
    assert inverse.new == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
    )
    assert inverse != change
    assert inverse.inverse() == change


def test_stereobondfieldchange_fields():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationAst.Undetermined(),
        new=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Undetermined()
        ),
    )

    assert change.old == StereoConfigurationAst.Undetermined()
    assert change.new == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Undetermined()
    )
    assert repr(change) == (
        "StereoBondFieldChange.Configuration("
        "old=StereoConfigurationAst.Undetermined(), "
        "new=StereoConfigurationAst.Kinded("
        "StereoKind.CisTrans, StereoCosetAst.Undetermined()))"
    )
    with pytest.raises(AttributeError):
        change.old = StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Lit(0)
        )
    with pytest.raises(TypeError):
        hash(change)


def test_stereobondfieldchange_match_positional():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Undetermined()
        ),
        new=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Lit(1)
        ),
    )

    match change:
        case StereoBondFieldChange.Configuration(old, new):
            assert old == StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Undetermined()
            )
            assert new == StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Lit(1)
            )
        case _:
            raise AssertionError("stereo bond field change did not match its variant")


def test_stereobondfieldchange_match_named():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationAst.Undetermined(),
        new=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Undetermined()
        ),
    )

    match change:
        case StereoBondFieldChange.Configuration(old=old, new=new):
            assert old == StereoConfigurationAst.Undetermined()
            assert new == StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Undetermined()
            )
        case _:
            raise AssertionError("stereo bond field change did not match its variant")


def test_stereobondfieldchange_inverse():
    change = StereoBondFieldChange.Configuration(
        old=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Undetermined()
        ),
        new=StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Lit(1)
        ),
    )

    inverse = change.inverse()

    assert isinstance(inverse, StereoBondFieldChange.Configuration)
    assert inverse.old == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(1)
    )
    assert inverse.new == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Undetermined()
    )
    assert inverse != change
    assert inverse.inverse() == change


@pytest.mark.parametrize(
    ("change", "expected_repr"),
    [
        (
            StereoAtomFieldChange.Configuration(
                old=StereoConfigurationAst.Undetermined(),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
                ),
            ),
            "StereoAtomFieldChange.Configuration("
            "old=StereoConfigurationAst.Undetermined(), "
            "new=StereoConfigurationAst.Kinded("
            "StereoKind.Tetrahedral, StereoCosetAst.Undetermined()))",
        ),
        (
            StereoAtomFieldChange.Configuration(
                old=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCosetAst.Undetermined()
                ),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
                ),
            ),
            "StereoAtomFieldChange.Configuration("
            "old=StereoConfigurationAst.Kinded("
            "StereoKind.Tetrahedral, StereoCosetAst.Undetermined()), "
            "new=StereoConfigurationAst.Kinded("
            "StereoKind.Tetrahedral, StereoCosetAst.Lit(1)))",
        ),
        (
            StereoBondFieldChange.Configuration(
                old=StereoConfigurationAst.Undetermined(),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCosetAst.Undetermined()
                ),
            ),
            "StereoBondFieldChange.Configuration("
            "old=StereoConfigurationAst.Undetermined(), "
            "new=StereoConfigurationAst.Kinded("
            "StereoKind.CisTrans, StereoCosetAst.Undetermined()))",
        ),
        (
            StereoBondFieldChange.Configuration(
                old=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCosetAst.Undetermined()
                ),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCosetAst.Lit(1)
                ),
            ),
            "StereoBondFieldChange.Configuration("
            "old=StereoConfigurationAst.Kinded("
            "StereoKind.CisTrans, StereoCosetAst.Undetermined()), "
            "new=StereoConfigurationAst.Kinded("
            "StereoKind.CisTrans, StereoCosetAst.Lit(1)))",
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
                old=ElementAst.Lit(Element("C")),
                new=ElementAst.Lit(Element("N")),
            ),
            "AtomFieldChange.Element(old=ElementAst.Lit(Element('C')), "
            "new=ElementAst.Lit(Element('N')))",
        ),
        (
            AtomFieldChange.IsotopeMass(
                old=IsotopeMassAst.Lit(12),
                new=IsotopeMassAst.Lit(13),
            ),
            "AtomFieldChange.IsotopeMass(old=IsotopeMassAst.Lit(12), "
            "new=IsotopeMassAst.Lit(13))",
        ),
        (
            AtomFieldChange.Charge(
                old=ValueAst.Lit(0),
                new=ValueAst.Lit(-1),
            ),
            "AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1))",
        ),
        (
            AtomFieldChange.ImplicitHydrogens(
                old=ValueAst.Lit(3),
                new=ValueAst.Lit(2),
            ),
            "AtomFieldChange.ImplicitHydrogens(old=ValueAst.Lit(3), "
            "new=ValueAst.Lit(2))",
        ),
        (
            AtomFieldChange.LonePairs(
                old=ValueAst.Lit(1),
                new=ValueAst.Lit(2),
            ),
            "AtomFieldChange.LonePairs(old=ValueAst.Lit(1), new=ValueAst.Lit(2))",
        ),
        (
            AtomFieldChange.Spin(
                old=SpinStateAst(0, 1),
                new=SpinStateAst(1, 2),
            ),
            "AtomFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), "
            "ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2)))",
        ),
        (
            BondFieldChange.Order(
                old=ValueAst.Lit(1),
                new=ValueAst.Lit(2),
            ),
            "BondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2))",
        ),
        (
            BondFieldChange.Charge(
                old=ValueAst.Lit(0),
                new=ValueAst.Lit(1),
            ),
            "BondFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(1))",
        ),
        (
            BondFieldChange.Spin(
                old=SpinStateAst(0, 1),
                new=SpinStateAst(1, 2),
            ),
            "BondFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), "
            "ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2)))",
        ),
        (
            DativeBondFieldChange.Order(
                old=ValueAst.Lit(1),
                new=ValueAst.Lit(2),
            ),
            "DativeBondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2))",
        ),
        (
            AromaticSystemFieldChange.Electrons(
                old=ElectronCountsAst.Undetermined(),
                new=ElectronCountsAst.Lit([1, 1, 1]),
            ),
            "AromaticSystemFieldChange.Electrons("
            "old=ElectronCountsAst.Undetermined(), "
            "new=ElectronCountsAst.Lit([1, 1, 1]))",
        ),
        (
            AromaticSystemFieldChange.Charge(
                old=ValueAst.Lit(0),
                new=ValueAst.Lit(-1),
            ),
            "AromaticSystemFieldChange.Charge(old=ValueAst.Lit(0), "
            "new=ValueAst.Lit(-1))",
        ),
        (
            AromaticSystemFieldChange.Spin(
                old=SpinStateAst(0, 1),
                new=SpinStateAst(1, 2),
            ),
            "AromaticSystemFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), "
            "ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2)))",
        ),
        (
            MulticenterBondFieldChange.Electrons(
                old=ElectronCountsAst.Lit([1, 0, 1]),
                new=ElectronCountsAst.Lit([2, 0, 1]),
            ),
            "MulticenterBondFieldChange.Electrons("
            "old=ElectronCountsAst.Lit([1, 0, 1]), "
            "new=ElectronCountsAst.Lit([2, 0, 1]))",
        ),
        (
            MulticenterBondFieldChange.Charge(
                old=ValueAst.Lit(0),
                new=ValueAst.Lit(1),
            ),
            "MulticenterBondFieldChange.Charge(old=ValueAst.Lit(0), "
            "new=ValueAst.Lit(1))",
        ),
        (
            MulticenterBondFieldChange.Spin(
                old=SpinStateAst(0, 1),
                new=SpinStateAst(2, 3),
            ),
            "MulticenterBondFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), "
            "ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(2), ValueAst.Lit(3)))",
        ),
        (
            NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindAst.Undetermined(),
                new=NoncovalentBondKindAst.Lit(
                    NoncovalentBondKind.HydrogenBond
                ),
            ),
            "NoncovalentBondFieldChange.Kind("
            "old=NoncovalentBondKindAst.Undetermined(), "
            "new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond))",
        ),
    ],
    ids=[
        "atom-element",
        "atom-isotope-mass",
        "atom-charge",
        "atom-implicit-hydrogens",
        "atom-lone-pairs",
        "atom-spin",
        "bond-order",
        "bond-charge",
        "bond-spin",
        "dative-order",
        "aromatic-electrons",
        "aromatic-charge",
        "aromatic-spin",
        "multicenter-electrons",
        "multicenter-charge",
        "multicenter-spin",
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
    source = AtomAst(Element("C"))
    delta = AtomDelta.Add(id=3, ast=source)

    source.charge = -1

    assert delta.id == 3
    assert delta.ast.charge == ValueAst.Undetermined()
    assert repr(delta) == "AtomDelta.Add(id=3, ast=AtomAst.parse('C'))"
    delta.ast.charge = 1
    assert delta.ast.charge == ValueAst.Lit(1)
    with pytest.raises(AttributeError):
        delta.id = 4
    with pytest.raises(TypeError):
        hash(delta)


def test_atomdelta_add_match():
    delta = AtomDelta.Add(id=3, ast=AtomAst(Element("C")))

    match delta:
        case AtomDelta.Add(id=id, ast=ast):
            assert id == 3
            assert ast == AtomAst(Element("C"))
        case _:
            raise AssertionError("atom delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AtomDelta.Remove)
    assert inverse.id == 3
    assert inverse.ast == AtomAst(Element("C"))
    assert inverse.inverse() == delta


def test_atomdelta_modifyfield_match():
    delta = AtomDelta.ModifyField(
        id=3,
        change=AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1)),
    )

    match delta:
        case AtomDelta.ModifyField(id, change):
            assert id == 3
            assert change == AtomFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
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
        new=AtomConstraintAst.Valence(ValueAst.Lit(4)),
    )

    match delta:
        case AtomDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 3
            assert old is None
            assert new == AtomConstraintAst.Valence(ValueAst.Lit(4))
        case _:
            raise AssertionError("atom delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AtomDelta.ModifyConstraint)
    assert inverse.old == AtomConstraintAst.Valence(ValueAst.Lit(4))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_bonddelta_fields():
    source = BondAst(1)
    delta = BondDelta.Add(id=2, atoms=(5, 1), ast=source)

    source.order = 2

    assert delta.id == 2
    assert delta.atoms == (5, 1)
    assert isinstance(delta.atoms, tuple)
    assert delta.ast.order == ValueAst.Lit(1)
    assert repr(delta) == (
        "BondDelta.Add(id=2, atoms=(5, 1), ast=BondAst.parse('1'))"
    )
    delta.ast.order = 3
    assert delta.ast.order == ValueAst.Lit(3)
    with pytest.raises(AttributeError):
        delta.atoms = (1, 5)
    with pytest.raises(TypeError):
        hash(delta)


def test_bonddelta_add_match():
    delta = BondDelta.Add(id=2, atoms=(5, 1), ast=BondAst(1))

    match delta:
        case BondDelta.Add(id=id, atoms=atoms, ast=ast):
            assert id == 2
            assert atoms == (5, 1)
            assert ast == BondAst(1)
        case _:
            raise AssertionError("bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, BondDelta.Remove)
    assert inverse.id == 2
    assert inverse.atoms == (5, 1)
    assert inverse.ast == BondAst(1)
    assert inverse.inverse() == delta


def test_bonddelta_modifyfield_match():
    delta = BondDelta.ModifyField(
        id=2,
        change=BondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2)),
    )

    match delta:
        case BondDelta.ModifyField(id, change):
            assert id == 2
            assert change == BondFieldChange.Order(
                old=ValueAst.Lit(1), new=ValueAst.Lit(2)
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
        new=BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
    )

    match delta:
        case BondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 2
            assert old is None
            assert new == BondConstraintAst.Aromatic(BooleanAst.Lit(True))
        case _:
            raise AssertionError("bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, BondDelta.ModifyConstraint)
    assert inverse.old == BondConstraintAst.Aromatic(BooleanAst.Lit(True))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_dativebonddelta_fields():
    source = DativeBondAst(1)
    delta = DativeBondDelta.Add(
        id=1,
        donors=[4, 2, 4],
        acceptor=3,
        ast=source,
    )

    source.order = 2

    assert delta.id == 1
    assert delta.donors == [4, 2, 4]
    assert isinstance(delta.donors, list)
    assert delta.acceptor == 3
    assert delta.ast.order == ValueAst.Lit(1)
    assert repr(delta) == (
        "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, "
        "ast=DativeBondAst.parse('1'))"
    )
    delta.ast.order = 3
    assert delta.ast.order == ValueAst.Lit(3)
    with pytest.raises(AttributeError):
        delta.donors = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_dativebonddelta_add_match():
    delta = DativeBondDelta.Add(
        id=1,
        donors=[4, 2, 4],
        acceptor=3,
        ast=DativeBondAst(1),
    )

    match delta:
        case DativeBondDelta.Add(id=id, donors=donors, acceptor=acceptor, ast=ast):
            assert id == 1
            assert donors == [4, 2, 4]
            assert acceptor == 3
            assert ast == DativeBondAst(1)
        case _:
            raise AssertionError("dative bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, DativeBondDelta.Remove)
    assert inverse.id == 1
    assert inverse.donors == [4, 2, 4]
    assert inverse.acceptor == 3
    assert inverse.ast == DativeBondAst(1)
    assert inverse.inverse() == delta


def test_dativebonddelta_modifyfield_match():
    delta = DativeBondDelta.ModifyField(
        id=1,
        change=DativeBondFieldChange.Order(
            old=ValueAst.Lit(1), new=ValueAst.Lit(2)
        ),
    )

    match delta:
        case DativeBondDelta.ModifyField(id, change):
            assert id == 1
            assert change == DativeBondFieldChange.Order(
                old=ValueAst.Lit(1), new=ValueAst.Lit(2)
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
        new=DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
    )

    match delta:
        case DativeBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 1
            assert old is None
            assert new == DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))
        case _:
            raise AssertionError("dative bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, DativeBondDelta.ModifyConstraint)
    assert inverse.old == DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_aromaticsystemdelta_fields():
    source = AromaticSystemAst([1, 1, 1])
    delta = AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], ast=source)

    source.electrons = [2, 0, 1]

    assert delta.id == 2
    assert delta.atoms == [4, 2, 4]
    assert isinstance(delta.atoms, list)
    assert delta.ast.electrons == ElectronCountsAst.Lit([1, 1, 1])
    assert repr(delta) == (
        "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], "
        "ast=AromaticSystemAst.parse('[1,1,1]'))"
    )
    delta.ast.charge = -1
    assert delta.ast.charge == ValueAst.Lit(-1)
    with pytest.raises(AttributeError):
        delta.atoms = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_aromaticsystemdelta_add_match():
    delta = AromaticSystemDelta.Add(
        id=2,
        atoms=[4, 2, 4],
        ast=AromaticSystemAst([1, 1, 1]),
    )

    match delta:
        case AromaticSystemDelta.Add(id=id, atoms=atoms, ast=ast):
            assert id == 2
            assert atoms == [4, 2, 4]
            assert ast == AromaticSystemAst([1, 1, 1])
        case _:
            raise AssertionError("aromatic system delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, AromaticSystemDelta.Remove)
    assert inverse.id == 2
    assert inverse.atoms == [4, 2, 4]
    assert inverse.ast == AromaticSystemAst([1, 1, 1])
    assert inverse.inverse() == delta


def test_aromaticsystemdelta_modifyfield_match():
    delta = AromaticSystemDelta.ModifyField(
        id=2,
        change=AromaticSystemFieldChange.Charge(
            old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
        ),
    )

    match delta:
        case AromaticSystemDelta.ModifyField(id, change):
            assert id == 2
            assert change == AromaticSystemFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
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
        new=AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
    )

    match delta:
        case AromaticSystemDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 2
            assert old is None
            assert new == AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))
        case _:
            raise AssertionError(
                "aromatic system delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, AromaticSystemDelta.ModifyConstraint)
    assert inverse.old == AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_multicenterbonddelta_fields():
    source = MulticenterBondAst([1, 1, 1])
    delta = MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], ast=source)

    source.electrons = [2, 0, 1]

    assert delta.id == 3
    assert delta.atoms == [4, 2, 4]
    assert isinstance(delta.atoms, list)
    assert delta.ast.electrons == ElectronCountsAst.Lit([1, 1, 1])
    assert repr(delta) == (
        "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], "
        "ast=MulticenterBondAst.parse('[1,1,1]'))"
    )
    delta.ast.charge = -1
    assert delta.ast.charge == ValueAst.Lit(-1)
    with pytest.raises(AttributeError):
        delta.atoms = [2, 4, 4]
    with pytest.raises(TypeError):
        hash(delta)


def test_multicenterbonddelta_add_match():
    delta = MulticenterBondDelta.Add(
        id=3,
        atoms=[4, 2, 4],
        ast=MulticenterBondAst([1, 1, 1]),
    )

    match delta:
        case MulticenterBondDelta.Add(id=id, atoms=atoms, ast=ast):
            assert id == 3
            assert atoms == [4, 2, 4]
            assert ast == MulticenterBondAst([1, 1, 1])
        case _:
            raise AssertionError("multicenter bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, MulticenterBondDelta.Remove)
    assert inverse.id == 3
    assert inverse.atoms == [4, 2, 4]
    assert inverse.ast == MulticenterBondAst([1, 1, 1])
    assert inverse.inverse() == delta


def test_multicenterbonddelta_modifyfield_match():
    delta = MulticenterBondDelta.ModifyField(
        id=3,
        change=MulticenterBondFieldChange.Charge(
            old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
        ),
    )

    match delta:
        case MulticenterBondDelta.ModifyField(id, change):
            assert id == 3
            assert change == MulticenterBondFieldChange.Charge(
                old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
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
        new=MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6)),
    )

    match delta:
        case MulticenterBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 3
            assert old is None
            assert new == MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))
        case _:
            raise AssertionError(
                "multicenter bond delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, MulticenterBondDelta.ModifyConstraint)
    assert inverse.old == MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_noncovalentbonddelta_fields():
    source = NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    delta = NoncovalentBondDelta.Add(id=4, atoms=(5, 2), ast=source)

    source.kind = NoncovalentBondKind.Ionic

    assert delta.id == 4
    assert delta.atoms == (5, 2)
    assert isinstance(delta.atoms, tuple)
    assert delta.ast.kind == NoncovalentBondKindAst.Lit(
        NoncovalentBondKind.HydrogenBond
    )
    assert repr(delta) == (
        "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), "
        "ast=NoncovalentBondAst.parse('Hbd'))"
    )
    delta.ast.kind = NoncovalentBondKind.Ionic
    assert delta.ast.kind == NoncovalentBondKindAst.Lit(NoncovalentBondKind.Ionic)
    with pytest.raises(AttributeError):
        delta.atoms = (2, 5)
    with pytest.raises(TypeError):
        hash(delta)


def test_noncovalentbonddelta_add_match():
    delta = NoncovalentBondDelta.Add(
        id=4,
        atoms=(5, 2),
        ast=NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
    )

    match delta:
        case NoncovalentBondDelta.Add(id=id, atoms=atoms, ast=ast):
            assert id == 4
            assert atoms == (5, 2)
            assert ast == NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
        case _:
            raise AssertionError("noncovalent bond delta did not match its add variant")

    inverse = delta.inverse()
    assert isinstance(inverse, NoncovalentBondDelta.Remove)
    assert inverse.id == 4
    assert inverse.atoms == (5, 2)
    assert inverse.ast == NoncovalentBondAst(NoncovalentBondKind.HydrogenBond)
    assert inverse.inverse() == delta


def test_noncovalentbonddelta_modifyfield_match():
    delta = NoncovalentBondDelta.ModifyField(
        id=4,
        change=NoncovalentBondFieldChange.Kind(
            old=NoncovalentBondKindAst.Undetermined(),
            new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
        ),
    )

    match delta:
        case NoncovalentBondDelta.ModifyField(id, change):
            assert id == 4
            assert change == NoncovalentBondFieldChange.Kind(
                old=NoncovalentBondKindAst.Undetermined(),
                new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond),
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
        new=NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)),
    )

    match delta:
        case NoncovalentBondDelta.ModifyConstraint(id=id, old=old, new=new):
            assert id == 4
            assert old is None
            assert new == NoncovalentBondConstraintAst.Intramolecular(
                BooleanAst.Lit(True)
            )
        case _:
            raise AssertionError(
                "noncovalent bond delta did not match its constraint variant"
            )

    inverse = delta.inverse()
    assert isinstance(inverse, NoncovalentBondDelta.ModifyConstraint)
    assert inverse.old == NoncovalentBondConstraintAst.Intramolecular(
        BooleanAst.Lit(True)
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_stereoatomdelta_fields():
    source = StereoAtomAst(
        StereoConfigurationAst.Kinded(
            StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
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
        ast=source,
    )

    source.configuration = StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
    )

    assert delta.id == 5
    assert delta.site == 3
    assert delta.ligands == [
        StereoLigand(4, StereoLigandKind.Atom),
        StereoLigand(2, StereoLigandKind.LonePair),
        StereoLigand(4, StereoLigandKind.Atom),
    ]
    assert isinstance(delta.ligands, list)
    assert delta.ast.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
    )
    assert repr(delta) == (
        "StereoAtomDelta.Add(id=5, site=3, ligands=["
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), "
        "StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair), "
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], "
        "ast=StereoAtomAst.parse('Th0'))"
    )
    delta.ast.configuration = StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
    )
    assert delta.ast.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(1)
    )
    with pytest.raises(AttributeError):
        delta.ligands = []
    with pytest.raises(TypeError):
        hash(delta)


def test_stereoatomdelta_add_match():
    delta = StereoAtomDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
        ],
        ast=StereoAtomAst(
            StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
            )
        ),
    )

    match delta:
        case StereoAtomDelta.Add(id=id, site=site, ligands=ligands, ast=ast):
            assert id == 5
            assert site == 3
            assert ligands == [
                StereoLigand(4, StereoLigandKind.Atom),
                StereoLigand(2, StereoLigandKind.LonePair),
            ]
            assert ast.configuration == StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
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
            old=StereoConfigurationAst.Undetermined(),
            new=StereoConfigurationAst.Kinded(
                StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
            ),
        ),
    )

    match delta:
        case StereoAtomDelta.ModifyField(id, change):
            assert id == 5
            assert change == StereoAtomFieldChange.Configuration(
                old=StereoConfigurationAst.Undetermined(),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
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
        new=StereoAtomConstraintAst.Stereogenicity(
            StereogenicityAst.Undetermined()
        ),
    )

    match delta:
        case StereoAtomDelta.ModifyConstraint(
            id=id, kind=kind, old=old, new=new
        ):
            assert id == 5
            assert kind is StereoKind.Tetrahedral
            assert old is None
            assert new == StereoAtomConstraintAst.Stereogenicity(
                StereogenicityAst.Undetermined()
            )
        case _:
            raise AssertionError("stereo atom delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoAtomDelta.ModifyConstraint)
    assert inverse.kind is StereoKind.Tetrahedral
    assert inverse.old == StereoAtomConstraintAst.Stereogenicity(
        StereogenicityAst.Undetermined()
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_stereoatomdelta_modifyconstraint_kind_none():
    delta = StereoAtomDelta.ModifyConstraint(
        id=5,
        kind=None,
        old=StereoAtomConstraintAst.Stereogenicity(
            StereogenicityAst.Undetermined()
        ),
        new=None,
    )

    assert delta.kind is None
    inverse = delta.inverse()
    assert inverse.kind is None
    assert inverse.old is None
    assert inverse.new == StereoAtomConstraintAst.Stereogenicity(
        StereogenicityAst.Undetermined()
    )
    assert inverse.inverse() == delta


def test_stereoatomdelta_apply_match():
    delta = StereoAtomDelta.Apply(
        id=5,
        kind=StereoKind.Tetrahedral,
        permutation=Permutation([1, 2, 0, 3]),
    )

    match delta:
        case StereoAtomDelta.Apply(id, kind, permutation):
            assert id == 5
            assert kind is StereoKind.Tetrahedral
            assert permutation.degree == 4
            assert permutation.image() == [1, 2, 0, 3]
        case _:
            raise AssertionError("stereo atom delta did not match its apply variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoAtomDelta.Apply)
    assert inverse.kind is StereoKind.Tetrahedral
    assert inverse.permutation.degree == 4
    assert inverse.permutation.image() == [2, 0, 1, 3]
    assert inverse != delta
    assert inverse.inverse() == delta


def test_stereoatomdelta_involutions():
    swap = StereoAtomDelta.Swap(id=5, kind=StereoKind.Tetrahedral)
    mirror = StereoAtomDelta.Mirror(id=5, kind=StereoKind.Tetrahedral)

    match swap:
        case StereoAtomDelta.Swap(id=id, kind=kind):
            assert id == 5
            assert kind is StereoKind.Tetrahedral
        case _:
            raise AssertionError("stereo atom delta did not match its swap variant")
    match mirror:
        case StereoAtomDelta.Mirror(id=id, kind=kind):
            assert id == 5
            assert kind is StereoKind.Tetrahedral
        case _:
            raise AssertionError("stereo atom delta did not match its mirror variant")

    assert isinstance(swap.inverse(), StereoAtomDelta.Swap)
    assert swap.inverse() == swap
    assert isinstance(mirror.inverse(), StereoAtomDelta.Mirror)
    assert mirror.inverse() == mirror

def test_stereobonddelta_fields():
    source = StereoBondAst(
        StereoConfigurationAst.Kinded(
            StereoKind.CisTrans, StereoCosetAst.Lit(0)
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
        ast=source,
    )

    source.configuration = StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(1)
    )

    assert delta.id == 5
    assert delta.site == 3
    assert delta.ligands == [
        StereoLigand(4, StereoLigandKind.Atom),
        StereoLigand(2, StereoLigandKind.LonePair),
        StereoLigand(4, StereoLigandKind.Atom),
    ]
    assert isinstance(delta.ligands, list)
    assert delta.ast.configuration == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(0)
    )
    assert repr(delta) == (
        "StereoBondDelta.Add(id=5, site=3, ligands=["
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), "
        "StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair), "
        "StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], "
        "ast=StereoBondAst.parse('Ct0'))"
    )
    delta.ast.configuration = StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(1)
    )
    assert delta.ast.configuration == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(1)
    )
    with pytest.raises(AttributeError):
        delta.ligands = []
    with pytest.raises(TypeError):
        hash(delta)


def test_stereobonddelta_add_match():
    delta = StereoBondDelta.Add(
        id=5,
        site=3,
        ligands=[
            StereoLigand(4, StereoLigandKind.Atom),
            StereoLigand(2, StereoLigandKind.LonePair),
        ],
        ast=StereoBondAst(
            StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Lit(0)
            )
        ),
    )

    match delta:
        case StereoBondDelta.Add(id=id, site=site, ligands=ligands, ast=ast):
            assert id == 5
            assert site == 3
            assert ligands == [
                StereoLigand(4, StereoLigandKind.Atom),
                StereoLigand(2, StereoLigandKind.LonePair),
            ]
            assert ast.configuration == StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Lit(0)
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
            old=StereoConfigurationAst.Undetermined(),
            new=StereoConfigurationAst.Kinded(
                StereoKind.CisTrans, StereoCosetAst.Lit(0)
            ),
        ),
    )

    match delta:
        case StereoBondDelta.ModifyField(id, change):
            assert id == 5
            assert change == StereoBondFieldChange.Configuration(
                old=StereoConfigurationAst.Undetermined(),
                new=StereoConfigurationAst.Kinded(
                    StereoKind.CisTrans, StereoCosetAst.Lit(0)
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
        new=StereoBondConstraintAst.Stereogenicity(
            StereogenicityAst.Undetermined()
        ),
    )

    match delta:
        case StereoBondDelta.ModifyConstraint(
            id=id, kind=kind, old=old, new=new
        ):
            assert id == 5
            assert kind is StereoKind.CisTrans
            assert old is None
            assert new == StereoBondConstraintAst.Stereogenicity(
                StereogenicityAst.Undetermined()
            )
        case _:
            raise AssertionError("stereo bond delta did not match its constraint variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoBondDelta.ModifyConstraint)
    assert inverse.kind is StereoKind.CisTrans
    assert inverse.old == StereoBondConstraintAst.Stereogenicity(
        StereogenicityAst.Undetermined()
    )
    assert inverse.new is None
    assert inverse.inverse() == delta


def test_stereobonddelta_modifyconstraint_kind_none():
    delta = StereoBondDelta.ModifyConstraint(
        id=5,
        kind=None,
        old=StereoBondConstraintAst.Stereogenicity(
            StereogenicityAst.Undetermined()
        ),
        new=None,
    )

    assert delta.kind is None
    inverse = delta.inverse()
    assert inverse.kind is None
    assert inverse.old is None
    assert inverse.new == StereoBondConstraintAst.Stereogenicity(
        StereogenicityAst.Undetermined()
    )
    assert inverse.inverse() == delta


def test_stereobonddelta_apply_match():
    delta = StereoBondDelta.Apply(
        id=5,
        kind=StereoKind.CisTrans,
        permutation=Permutation([1, 2, 0, 3]),
    )

    match delta:
        case StereoBondDelta.Apply(id, kind, permutation):
            assert id == 5
            assert kind is StereoKind.CisTrans
            assert permutation.degree == 4
            assert permutation.image() == [1, 2, 0, 3]
        case _:
            raise AssertionError("stereo bond delta did not match its apply variant")

    inverse = delta.inverse()
    assert isinstance(inverse, StereoBondDelta.Apply)
    assert inverse.kind is StereoKind.CisTrans
    assert inverse.permutation.degree == 4
    assert inverse.permutation.image() == [2, 0, 1, 3]
    assert inverse != delta
    assert inverse.inverse() == delta


def test_stereobonddelta_involutions():
    swap = StereoBondDelta.Swap(id=5, kind=StereoKind.CisTrans)
    mirror = StereoBondDelta.Mirror(id=5, kind=StereoKind.CisTrans)

    match swap:
        case StereoBondDelta.Swap(id=id, kind=kind):
            assert id == 5
            assert kind is StereoKind.CisTrans
        case _:
            raise AssertionError("stereo bond delta did not match its swap variant")
    match mirror:
        case StereoBondDelta.Mirror(id=id, kind=kind):
            assert id == 5
            assert kind is StereoKind.CisTrans
        case _:
            raise AssertionError("stereo bond delta did not match its mirror variant")

    assert isinstance(swap.inverse(), StereoBondDelta.Swap)
    assert swap.inverse() == swap
    assert isinstance(mirror.inverse(), StereoBondDelta.Mirror)
    assert mirror.inverse() == mirror



@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type"),
    [
        (
            AtomDelta.Add(id=3, ast=AtomAst(Element("C"))),
            "AtomDelta.Add(id=3, ast=AtomAst.parse('C'))",
            AtomDelta.Remove,
        ),
        (
            AtomDelta.Remove(id=3, ast=AtomAst(Element("N"))),
            "AtomDelta.Remove(id=3, ast=AtomAst.parse('N'))",
            AtomDelta.Add,
        ),
        (
            AtomDelta.ModifyField(
                id=3,
                change=AtomFieldChange.Charge(
                    old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
                ),
            ),
            "AtomDelta.ModifyField(id=3, "
            "change=AtomFieldChange.Charge("
            "old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
            AtomDelta.ModifyField,
        ),
        (
            AtomDelta.ModifyConstraint(
                id=3,
                old=None,
                new=AtomConstraintAst.Valence(ValueAst.Lit(4)),
            ),
            "AtomDelta.ModifyConstraint(id=3, old=None, "
            "new=AtomConstraintAst.Valence(ValueAst.Lit(4)))",
            AtomDelta.ModifyConstraint,
        ),
        (
            BondDelta.Add(id=2, atoms=(5, 1), ast=BondAst(1)),
            "BondDelta.Add(id=2, atoms=(5, 1), ast=BondAst.parse('1'))",
            BondDelta.Remove,
        ),
        (
            BondDelta.Remove(id=2, atoms=(1, 5), ast=BondAst(2)),
            "BondDelta.Remove(id=2, atoms=(1, 5), ast=BondAst.parse('2'))",
            BondDelta.Add,
        ),
        (
            BondDelta.ModifyField(
                id=2,
                change=BondFieldChange.Order(
                    old=ValueAst.Lit(1), new=ValueAst.Lit(2)
                ),
            ),
            "BondDelta.ModifyField(id=2, "
            "change=BondFieldChange.Order("
            "old=ValueAst.Lit(1), new=ValueAst.Lit(2)))",
            BondDelta.ModifyField,
        ),
        (
            BondDelta.ModifyConstraint(
                id=2,
                old=None,
                new=BondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            ),
            "BondDelta.ModifyConstraint(id=2, old=None, "
            "new=BondConstraintAst.Aromatic(BooleanAst.Lit(True)))",
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
                id=1, donors=[4, 2, 4], acceptor=3, ast=DativeBondAst(1)
            ),
            "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, "
            "ast=DativeBondAst.parse('1'))",
            DativeBondDelta.Remove,
        ),
        (
            DativeBondDelta.Remove(
                id=1, donors=[2, 4, 2], acceptor=3, ast=DativeBondAst(2)
            ),
            "DativeBondDelta.Remove(id=1, donors=[2, 4, 2], acceptor=3, "
            "ast=DativeBondAst.parse('2'))",
            DativeBondDelta.Add,
        ),
        (
            DativeBondDelta.ModifyField(
                id=1,
                change=DativeBondFieldChange.Order(
                    old=ValueAst.Lit(1), new=ValueAst.Lit(2)
                ),
            ),
            "DativeBondDelta.ModifyField(id=1, "
            "change=DativeBondFieldChange.Order("
            "old=ValueAst.Lit(1), new=ValueAst.Lit(2)))",
            DativeBondDelta.ModifyField,
        ),
        (
            DativeBondDelta.ModifyConstraint(
                id=1,
                old=None,
                new=DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)),
            ),
            "DativeBondDelta.ModifyConstraint(id=1, old=None, "
            "new=DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)))",
            DativeBondDelta.ModifyConstraint,
        ),
        (
            AromaticSystemDelta.Add(
                id=2, atoms=[4, 2, 4], ast=AromaticSystemAst([1, 1, 1])
            ),
            "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], "
            "ast=AromaticSystemAst.parse('[1,1,1]'))",
            AromaticSystemDelta.Remove,
        ),
        (
            AromaticSystemDelta.Remove(
                id=2, atoms=[2, 4, 2], ast=AromaticSystemAst([2, 0, 1])
            ),
            "AromaticSystemDelta.Remove(id=2, atoms=[2, 4, 2], "
            "ast=AromaticSystemAst.parse('[2,0,1]'))",
            AromaticSystemDelta.Add,
        ),
        (
            AromaticSystemDelta.ModifyField(
                id=2,
                change=AromaticSystemFieldChange.Charge(
                    old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
                ),
            ),
            "AromaticSystemDelta.ModifyField(id=2, "
            "change=AromaticSystemFieldChange.Charge("
            "old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
            AromaticSystemDelta.ModifyField,
        ),
        (
            AromaticSystemDelta.ModifyConstraint(
                id=2,
                old=None,
                new=AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)),
            ),
            "AromaticSystemDelta.ModifyConstraint(id=2, old=None, "
            "new=AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)))",
            AromaticSystemDelta.ModifyConstraint,
        ),
        (
            MulticenterBondDelta.Add(
                id=3, atoms=[4, 2, 4], ast=MulticenterBondAst([1, 1, 1])
            ),
            "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], "
            "ast=MulticenterBondAst.parse('[1,1,1]'))",
            MulticenterBondDelta.Remove,
        ),
        (
            MulticenterBondDelta.Remove(
                id=3, atoms=[2, 4, 2], ast=MulticenterBondAst([2, 0, 1])
            ),
            "MulticenterBondDelta.Remove(id=3, atoms=[2, 4, 2], "
            "ast=MulticenterBondAst.parse('[2,0,1]'))",
            MulticenterBondDelta.Add,
        ),
        (
            MulticenterBondDelta.ModifyField(
                id=3,
                change=MulticenterBondFieldChange.Charge(
                    old=ValueAst.Lit(0), new=ValueAst.Lit(-1)
                ),
            ),
            "MulticenterBondDelta.ModifyField(id=3, "
            "change=MulticenterBondFieldChange.Charge("
            "old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
            MulticenterBondDelta.ModifyField,
        ),
        (
            MulticenterBondDelta.ModifyConstraint(
                id=3,
                old=None,
                new=MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6)),
            ),
            "MulticenterBondDelta.ModifyConstraint(id=3, old=None, "
            "new=MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6)))",
            MulticenterBondDelta.ModifyConstraint,
        ),
        (
            NoncovalentBondDelta.Add(
                id=4,
                atoms=(5, 2),
                ast=NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
            ),
            "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), "
            "ast=NoncovalentBondAst.parse('Hbd'))",
            NoncovalentBondDelta.Remove,
        ),
        (
            NoncovalentBondDelta.Remove(
                id=4,
                atoms=(2, 5),
                ast=NoncovalentBondAst(NoncovalentBondKind.Ionic),
            ),
            "NoncovalentBondDelta.Remove(id=4, atoms=(2, 5), "
            "ast=NoncovalentBondAst.parse('Ion'))",
            NoncovalentBondDelta.Add,
        ),
        (
            NoncovalentBondDelta.ModifyField(
                id=4,
                change=NoncovalentBondFieldChange.Kind(
                    old=NoncovalentBondKindAst.Undetermined(),
                    new=NoncovalentBondKindAst.Lit(
                        NoncovalentBondKind.HydrogenBond
                    ),
                ),
            ),
            "NoncovalentBondDelta.ModifyField(id=4, "
            "change=NoncovalentBondFieldChange.Kind("
            "old=NoncovalentBondKindAst.Undetermined(), "
            "new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)))",
            NoncovalentBondDelta.ModifyField,
        ),
        (
            NoncovalentBondDelta.ModifyConstraint(
                id=4,
                old=None,
                new=NoncovalentBondConstraintAst.Intramolecular(
                    BooleanAst.Lit(True)
                ),
            ),
            "NoncovalentBondDelta.ModifyConstraint(id=4, old=None, "
            "new=NoncovalentBondConstraintAst.Intramolecular("
            "BooleanAst.Lit(True)))",
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
                ast=StereoAtomAst(
                    StereoConfigurationAst.Kinded(
                        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
                    )
                ),
            ),
            "StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], ast=StereoAtomAst.parse('Th0'))",
            StereoAtomDelta.Remove,
            False,
        ),
        (
            StereoAtomDelta.Remove(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                ast=StereoAtomAst(
                    StereoConfigurationAst.Kinded(
                        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
                    )
                ),
            ),
            "StereoAtomDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], ast=StereoAtomAst.parse('Th0'))",
            StereoAtomDelta.Add,
            False,
        ),
        (
            StereoAtomDelta.ModifyField(
                id=5,
                change=StereoAtomFieldChange.Configuration(
                    old=StereoConfigurationAst.Undetermined(),
                    new=StereoConfigurationAst.Kinded(
                        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
                    ),
                ),
            ),
            "StereoAtomDelta.ModifyField(id=5, change=StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(0))))",
            StereoAtomDelta.ModifyField,
            False,
        ),
        (
            StereoAtomDelta.ModifyConstraint(
                id=5,
                kind=StereoKind.Tetrahedral,
                old=None,
                new=StereoAtomConstraintAst.Stereogenicity(
                    StereogenicityAst.Undetermined()
                ),
            ),
            "StereoAtomDelta.ModifyConstraint(id=5, kind=StereoKind.Tetrahedral, old=None, new=StereoAtomConstraintAst.Stereogenicity(StereogenicityAst.Undetermined()))",
            StereoAtomDelta.ModifyConstraint,
            False,
        ),
        (
            StereoAtomDelta.Apply(
                id=5,
                kind=StereoKind.Tetrahedral,
                permutation=Permutation([1, 2, 0, 3]),
            ),
            "StereoAtomDelta.Apply(id=5, kind=StereoKind.Tetrahedral, permutation=Permutation([1, 2, 0, 3]))",
            StereoAtomDelta.Apply,
            False,
        ),
        (
            StereoAtomDelta.Swap(id=5, kind=StereoKind.Tetrahedral),
            "StereoAtomDelta.Swap(id=5, kind=StereoKind.Tetrahedral)",
            StereoAtomDelta.Swap,
            True,
        ),
        (
            StereoAtomDelta.Mirror(id=5, kind=StereoKind.Tetrahedral),
            "StereoAtomDelta.Mirror(id=5, kind=StereoKind.Tetrahedral)",
            StereoAtomDelta.Mirror,
            True,
        ),
        (
            StereoBondDelta.Add(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                ast=StereoBondAst(
                    StereoConfigurationAst.Kinded(
                        StereoKind.CisTrans, StereoCosetAst.Lit(0)
                    )
                ),
            ),
            "StereoBondDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], ast=StereoBondAst.parse('Ct0'))",
            StereoBondDelta.Remove,
            False,
        ),
        (
            StereoBondDelta.Remove(
                id=5,
                site=3,
                ligands=[StereoLigand(4, StereoLigandKind.Atom)],
                ast=StereoBondAst(
                    StereoConfigurationAst.Kinded(
                        StereoKind.CisTrans, StereoCosetAst.Lit(0)
                    )
                ),
            ),
            "StereoBondDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], ast=StereoBondAst.parse('Ct0'))",
            StereoBondDelta.Add,
            False,
        ),
        (
            StereoBondDelta.ModifyField(
                id=5,
                change=StereoBondFieldChange.Configuration(
                    old=StereoConfigurationAst.Undetermined(),
                    new=StereoConfigurationAst.Kinded(
                        StereoKind.CisTrans, StereoCosetAst.Lit(0)
                    ),
                ),
            ),
            "StereoBondDelta.ModifyField(id=5, change=StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(0))))",
            StereoBondDelta.ModifyField,
            False,
        ),
        (
            StereoBondDelta.ModifyConstraint(
                id=5,
                kind=StereoKind.CisTrans,
                old=None,
                new=StereoBondConstraintAst.Stereogenicity(
                    StereogenicityAst.Undetermined()
                ),
            ),
            "StereoBondDelta.ModifyConstraint(id=5, kind=StereoKind.CisTrans, old=None, new=StereoBondConstraintAst.Stereogenicity(StereogenicityAst.Undetermined()))",
            StereoBondDelta.ModifyConstraint,
            False,
        ),
        (
            StereoBondDelta.Apply(
                id=5,
                kind=StereoKind.CisTrans,
                permutation=Permutation([1, 2, 0, 3]),
            ),
            "StereoBondDelta.Apply(id=5, kind=StereoKind.CisTrans, permutation=Permutation([1, 2, 0, 3]))",
            StereoBondDelta.Apply,
            False,
        ),
        (
            StereoBondDelta.Swap(id=5, kind=StereoKind.CisTrans),
            "StereoBondDelta.Swap(id=5, kind=StereoKind.CisTrans)",
            StereoBondDelta.Swap,
            True,
        ),
        (
            StereoBondDelta.Mirror(id=5, kind=StereoKind.CisTrans),
            "StereoBondDelta.Mirror(id=5, kind=StereoKind.CisTrans)",
            StereoBondDelta.Mirror,
            True,
        ),
    ],
    ids=[
        "atom-add",
        "atom-remove",
        "atom-modify-field",
        "atom-modify-constraint",
        "atom-apply",
        "atom-swap",
        "atom-mirror",
        "bond-add",
        "bond-remove",
        "bond-modify-field",
        "bond-modify-constraint",
        "bond-apply",
        "bond-swap",
        "bond-mirror",
    ],
)
def test_stereodelta_closure(delta, expected_repr, inverse_type, self_inverse):
    assert repr(delta) == expected_repr
    inverse = delta.inverse()
    assert type(inverse) is inverse_type
    assert (inverse == delta) is self_inverse
    assert inverse.inverse() == delta


def test_constraintdelta_fields():
    source = Constraint.Atom(3, AtomConstraintAst.Degree(ValueAst.Lit(2)))
    delta = ConstraintDelta.Add(constraint=source)

    assert delta.constraint == source
    assert delta.constraint is not source
    assert delta.constraint is delta.constraint
    assert repr(delta) == (
        "ConstraintDelta.Add(constraint=Constraint.Atom(3, "
        "AtomConstraintAst.Degree(ValueAst.Lit(2))))"
    )
    with pytest.raises(AttributeError):
        delta.constraint = Constraint.Or([])
    with pytest.raises(TypeError):
        hash(delta)


def test_constraintdelta_add_match():
    delta = ConstraintDelta.Add(
        constraint=Constraint.Atom(
            3,
            AtomConstraintAst.Degree(ValueAst.Lit(2)),
        )
    )

    match delta:
        case ConstraintDelta.Add(Constraint.Atom(atom_id, constraint)):
            assert (atom_id, constraint) == (
                3,
                AtomConstraintAst.Degree(ValueAst.Lit(2)),
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
                    AtomConstraintAst.Valence(ValueAst.Lit(4)),
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
                AtomConstraintAst.Valence(ValueAst.Lit(4)),
            )
        case _:
            raise AssertionError("constraint delta did not match its remove variant")

    inverse = delta.inverse()
    assert isinstance(inverse, ConstraintDelta.Add)
    assert inverse.constraint == delta.constraint
    assert inverse.inverse() == delta


def test_constraintdelta_payload_ownership():
    source_molecule = MoleculeAst.from_parts([AtomAst(Element("C"))])
    source = Constraint.Molecule(
        MoleculeConstraint.SubPattern(SubPatternAnchor(), source_molecule)
    )
    delta = ConstraintDelta.Add(constraint=source)

    source_molecule.atoms[0].charge = 1

    match delta.constraint:
        case Constraint.Molecule(MoleculeConstraint.SubPattern(_, stored_molecule)):
            assert stored_molecule.atoms[0].charge == ValueAst.Undetermined()
            stored_molecule.atoms[0].charge = -1
        case _:
            raise AssertionError("constraint delta did not retain its stored subpattern")

    inverse = delta.inverse()
    match inverse.constraint:
        case Constraint.Molecule(MoleculeConstraint.SubPattern(_, stored_molecule)):
            assert stored_molecule.atoms[0].charge == ValueAst.Lit(-1)
        case _:
            raise AssertionError("inverse did not retain the changed stored subpattern")


@pytest.mark.parametrize(
    ("delta", "expected_repr", "inverse_type"),
    [
        (
            ConstraintDelta.Add(
                constraint=Constraint.Atom(
                    3,
                    AtomConstraintAst.Degree(ValueAst.Lit(2)),
                )
            ),
            "ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintAst.Degree(ValueAst.Lit(2))))",
            ConstraintDelta.Remove,
        ),
        (
            ConstraintDelta.Remove(
                constraint=Constraint.And(
                    [
                        Constraint.Atom(
                            7,
                            AtomConstraintAst.Valence(ValueAst.Lit(4)),
                        ),
                        Constraint.Not(Constraint.Or([])),
                    ]
                )
            ),
            "ConstraintDelta.Remove(constraint=Constraint.And([Constraint.Atom(7, AtomConstraintAst.Valence(ValueAst.Lit(4))), Constraint.Not(Constraint.Or([]))]))",
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
