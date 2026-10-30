import pytest
from umol import (
    Molecule,
    MoleculeDefaults,
    StereoAtomConstraintForm,
    StereoAtomForm,
    StereoBondConstraintForm,
    StereoBondForm,
    StereoConfigurationForm,
    StereoCoset,
    StereoDerivation,
    StereogenicityForm,
    StereoInconsistency,
    StereoLigand,
    StereoLigandKind,
    StereoModel,
    StereoPerception,
    StereoSite,
)


@pytest.fixture
def stereo_input():
    return """{
        :atoms ["C#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2#C0"] [5 7 "1"] [6 8 "1"]]
    }"""


def test_stereo_perception_new(stereo_input):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    disabled = StereoModel(
        kind_models={}, para_stereo=False, stereo_bond_minimum_ring_size=0
    )
    perception = StereoPerception(disabled)
    assert perception.derive_stereo_atom(molecule, 0, StereoCoset.Lit(1)) is None
    with pytest.raises(TypeError):
        StereoPerception("default")


def test_stereo_perception_derive(stereo_input):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    result = StereoPerception(StereoModel.default()).derive(molecule)
    expected_atoms = (
        (
            0,
            tuple(StereoLigand(i, StereoLigandKind.Atom) for i in (1, 2, 3, 4)),
            StereoAtomForm.parse("Th1"),
        ),
    )
    expected_bonds = (
        (
            4,
            (
                StereoLigand(7, StereoLigandKind.Atom),
                StereoLigand(5, StereoLigandKind.ImplicitHydrogen),
                StereoLigand(8, StereoLigandKind.Atom),
                StereoLigand(6, StereoLigandKind.ImplicitHydrogen),
            ),
            StereoBondForm.parse("Ct0"),
        ),
    )
    assert result == StereoDerivation(
        stereo_atoms=expected_atoms, stereo_bonds=expected_bonds
    )
    for actual, expected in (
        (result.stereo_atoms, expected_atoms),
        (result.stereo_bonds, expected_bonds),
    ):
        assert tuple(
            (site, ligands, form.asdict()) for site, ligands, form in actual
        ) == tuple((site, ligands, form.asdict()) for site, ligands, form in expected)
    assert result.skipped_stereo_bonds == ()
    assert result.inconsistencies == ()
    assert molecule == Molecule.parse(
        stereo_input, defaults=MoleculeDefaults.concrete()
    )


@pytest.mark.parametrize(
    "value,expected",
    [
        pytest.param(0, StereoCoset.Lit(0), id="zero"),
        pytest.param(1, StereoCoset.Lit(1), id="one"),
    ],
)
def test_stereo_perception_derive_stereo_atom(stereo_input, value, expected):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    perception = StereoPerception(StereoModel.default())
    ligands, form = perception.derive_stereo_atom(molecule, 0, expected)
    assert ligands == tuple(
        StereoLigand(i, StereoLigandKind.Atom) for i in (1, 2, 3, 4)
    )
    assert form.asdict() == StereoAtomForm.parse(f"Th{value}").asdict()
    assert perception.derive_stereo_atom(molecule, 1, expected) is None


def test_stereo_perception_derive_stereo_bond(stereo_input):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    perception = StereoPerception(StereoModel.default())
    coset = StereoCoset.Undetermined()
    ligands, form = perception.derive_stereo_bond(molecule, 4, coset)
    assert ligands == (
        StereoLigand(7, StereoLigandKind.Atom),
        StereoLigand(5, StereoLigandKind.ImplicitHydrogen),
        StereoLigand(8, StereoLigandKind.Atom),
        StereoLigand(6, StereoLigandKind.ImplicitHydrogen),
    )
    assert form.asdict() == StereoBondForm.parse("Ct*").asdict()
    assert perception.derive_stereo_bond(molecule, 0, coset) is None


@pytest.mark.parametrize(
    "method,index,site",
    [
        ("derive_stereo_atom", 9, StereoSite.Atom(9)),
        ("derive_stereo_bond", 7, StereoSite.Bond(7)),
    ],
)
def test_stereo_perception_derive_stereo_error(stereo_input, method, index, site):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    perception = StereoPerception(StereoModel.default())
    with pytest.raises(IndexError) as caught:
        getattr(perception, method)(molecule, index, StereoCoset.Undetermined())
    assert caught.value.site == site
    assert molecule == Molecule.parse(
        stereo_input, defaults=MoleculeDefaults.concrete()
    )


def test_stereo_perception_skipped_stereo_bonds():
    molecule = Molecule.parse(
        """{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "2#C0"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )
    perception = StereoPerception(StereoModel.default())
    skipped = perception.skipped_stereo_bonds(molecule)
    assert skipped == frozenset({0})
    assert type(skipped) is frozenset
    result = perception.derive(molecule)
    assert result.skipped_stereo_bonds == (0,)
    assert result.stereo_bonds == ()
    assert result.inconsistencies == ()


def test_stereo_derivation_new():
    atom_form = StereoAtomForm.parse("Th1")
    bond_form = StereoBondForm.parse("Ct0")
    ligand = StereoLigand(1, StereoLigandKind.Atom)
    inconsistency = StereoInconsistency.TetrahedralStereoFailure(atom=3)
    result = StereoDerivation(
        stereo_atoms=[(0, [ligand], atom_form)],
        stereo_bonds=[(2, [ligand], bond_form)],
        skipped_stereo_bonds=[4, 2, 4],
        inconsistencies=[inconsistency, inconsistency],
    )
    atom_form.configuration = StereoConfigurationForm.Undetermined()
    bond_form.configuration = StereoConfigurationForm.Undetermined()
    assert result.stereo_atoms[0][:2] == (0, (ligand,))
    assert result.stereo_bonds[0][:2] == (2, (ligand,))
    assert result.stereo_atoms[0][2].asdict() == StereoAtomForm.parse("Th1").asdict()
    assert result.stereo_bonds[0][2].asdict() == StereoBondForm.parse("Ct0").asdict()
    assert result.skipped_stereo_bonds == (4, 2, 4)
    assert result.inconsistencies == (inconsistency, inconsistency)
    assert StereoDerivation() == StereoDerivation(
        stereo_atoms=(), stereo_bonds=(), skipped_stereo_bonds=(), inconsistencies=()
    )
    assert eval(repr(result), globals()) == result


def test_stereo_derivation_properties(stereo_input):
    molecule = Molecule.parse(stereo_input, defaults=MoleculeDefaults.concrete())
    result = StereoPerception(StereoModel.default()).derive(molecule)
    with pytest.raises(AttributeError):
        result.stereo_atoms = ()
    with pytest.raises(AttributeError):
        del result.stereo_bonds
    for rows, constraint in (
        (
            result.stereo_atoms,
            StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        ),
        (
            result.stereo_bonds,
            StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        ),
    ):
        assert type(rows) is tuple
        assert type(rows[0][1]) is tuple
        form = rows[0][2]
        assert form.readonly is True
        with pytest.raises(TypeError, match="read-only entity form"):
            form.configuration = StereoConfigurationForm.Undetermined()
        with pytest.raises(TypeError):
            form.constraints.set(constraint)
    assert eval(repr(result), globals()) == result


@pytest.mark.parametrize(
    "args,kwargs,error",
    [
        (([],), {}, TypeError),
        ((), {"skipped_stereo_bonds": [-1]}, OverflowError),
        ((), {"stereo_atoms": [(0, [], "Th0")]}, TypeError),
    ],
)
def test_stereo_derivation_new_error(args, kwargs, error):
    with pytest.raises(error):
        StereoDerivation(*args, **kwargs)


@pytest.mark.parametrize(
    "value,expected,fields",
    [
        (
            StereoInconsistency.TetrahedralStereoFailure(atom=2),
            "StereoInconsistency.TetrahedralStereoFailure(atom=2)",
            {"atom": 2},
        ),
        (
            StereoInconsistency.StereoAtomFailure(stereo_atom=1),
            "StereoInconsistency.StereoAtomFailure(stereo_atom=1)",
            {"stereo_atom": 1},
        ),
        (
            StereoInconsistency.TetrahedralStereoMismatch(atom=3, stereo_atom=1),
            "StereoInconsistency.TetrahedralStereoMismatch(atom=3, stereo_atom=1)",
            {"atom": 3, "stereo_atom": 1},
        ),
        (
            StereoInconsistency.CisTransStereoFailure(bond=2),
            "StereoInconsistency.CisTransStereoFailure(bond=2)",
            {"bond": 2},
        ),
        (
            StereoInconsistency.StereoBondFailure(stereo_bond=1),
            "StereoInconsistency.StereoBondFailure(stereo_bond=1)",
            {"stereo_bond": 1},
        ),
        (
            StereoInconsistency.CisTransStereoMismatch(bond=3, stereo_bond=1),
            "StereoInconsistency.CisTransStereoMismatch(bond=3, stereo_bond=1)",
            {"bond": 3, "stereo_bond": 1},
        ),
    ],
)
def test_stereo_inconsistency(value, expected, fields):
    assert repr(value) == expected
    assert eval(repr(value), globals()) == value
    assert hash(value) == hash(eval(repr(value), globals()))
    for name, field in fields.items():
        assert getattr(value, name) == field
        with pytest.raises(AttributeError):
            setattr(value, name, field + 1)


def test_stereo_inconsistency_order():
    first = StereoInconsistency.TetrahedralStereoFailure(atom=0)
    second = StereoInconsistency.TetrahedralStereoFailure(atom=2)
    third = StereoInconsistency.StereoAtomFailure(stereo_atom=0)
    assert sorted([third, second, first]) == [first, second, third]
    assert first <= second
    assert third > second
    assert third >= StereoInconsistency.StereoAtomFailure(stereo_atom=0)
    assert first != 0
