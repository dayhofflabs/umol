import pytest
from umol import (
    AddStereoConfig,
    ContradictionError,
    InvalidatedViewError,
    ModelConversionError,
    Molecule,
    MoleculeDefaults,
    RemoveStereoConfig,
    RemoveStereoMode,
    StereoAtomForm,
    StereoBondForm,
    Stereogenicity,
    StereoInconsistency,
    StereoKind,
    StereoModel,
    StereoSelection,
    StereoSite,
)


@pytest.fixture
def stereo_source():
    return Molecule.parse(
        """{
        :atoms ["C#D4#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2#C0"] [5 7 "1"] [6 8 "1"]]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )


@pytest.fixture
def stereo_determined():
    return Molecule.parse(
        """{
        :atoms ["C#D4#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2#C0"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct0"}]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )


@pytest.mark.parametrize("method", ["add_stereo", "add_stereo_into"])
@pytest.mark.parametrize("kwargs", [{}, {"config": AddStereoConfig()}])
def test_molecule_add_stereo(stereo_source, stereo_determined, method, kwargs):
    original = stereo_source.copy()
    alias = stereo_source
    atom = stereo_source.atoms[0]
    constraints = atom.constraints
    iterator = iter(stereo_source.atoms)
    result = getattr(stereo_source, method)(StereoModel.default(), **kwargs)
    if method == "add_stereo":
        assert result == stereo_determined
        assert stereo_source == original
        assert atom.element == original.atoms[0].element
        assert next(iterator).element == original.atoms[0].element
        assert len(constraints) == 2
        result.atoms[0].charge = 1
        assert alias == original
    else:
        assert result is None
        assert alias == stereo_determined
        with pytest.raises(InvalidatedViewError):
            _ = atom.element
        with pytest.raises(InvalidatedViewError):
            _ = len(constraints)
        with pytest.raises(InvalidatedViewError):
            next(iterator)


@pytest.mark.parametrize("method", ["add_stereo", "add_stereo_into"])
def test_molecule_add_stereo_identity(stereo_determined, method):
    original = stereo_determined.copy()
    atom = stereo_determined.atoms[0]
    result = getattr(stereo_determined, method)(StereoModel.default())
    if method == "add_stereo":
        assert result == original
        assert atom.element == original.atoms[0].element
    else:
        assert result is None
        with pytest.raises(InvalidatedViewError):
            _ = atom.element
    assert stereo_determined == original


@pytest.mark.parametrize("method", ["add_stereo", "add_stereo_into"])
def test_molecule_add_stereo_selection(stereo_source, method):
    expected = Molecule.parse(
        """{
        :atoms ["C#D4#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2#C0"] [5 7 "1"] [6 8 "1"]]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct0"}]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )
    config = AddStereoConfig(
        selection=StereoSelection(
            sites={StereoSite.Bond(4)}, kinds={StereoKind.CisTrans}
        )
    )
    result = getattr(stereo_source, method)(StereoModel.default(), config=config)
    assert (result if method == "add_stereo" else stereo_source) == expected
    if method.endswith("_into"):
        assert result is None


@pytest.mark.parametrize("method", ["remove_stereo", "remove_stereo_into"])
@pytest.mark.parametrize(
    "kwargs", [{}, {"config": RemoveStereoConfig(mode=RemoveStereoMode.All)}]
)
def test_molecule_remove_stereo(stereo_determined, method, kwargs):
    original = stereo_determined.copy()
    alias = stereo_determined
    atom = stereo_determined.atoms[0]
    result = getattr(stereo_determined, method)(StereoModel.default(), **kwargs)
    expected = Molecule.parse(
        """{
        :atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2"] [5 7 "1"] [6 8 "1"]]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )
    if method == "remove_stereo":
        assert result == expected
        assert alias == original
        assert atom.element == original.atoms[0].element
    else:
        assert result is None
        assert alias == expected
        with pytest.raises(InvalidatedViewError):
            _ = atom.element


@pytest.mark.parametrize("method", ["remove_stereo", "remove_stereo_into"])
def test_molecule_remove_stereo_identity(stereo_determined, method):
    original = stereo_determined.copy()
    atom = stereo_determined.atoms[0]
    result = getattr(stereo_determined, method)(
        StereoModel.default(),
        config=RemoveStereoConfig(mode=RemoveStereoMode.Undetermined),
    )
    if method == "remove_stereo":
        assert result == original
        assert atom.element == original.atoms[0].element
    else:
        assert result is None
        with pytest.raises(InvalidatedViewError):
            _ = atom.element
    assert stereo_determined == original


@pytest.mark.parametrize("method", ["remove_stereo", "remove_stereo_into"])
def test_molecule_remove_stereo_selection(stereo_determined, method):
    expected = Molecule.parse(
        """{
        :atoms ["C#D4#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )
    result = getattr(stereo_determined, method)(
        StereoModel.default(),
        config=RemoveStereoConfig(
            selection=StereoSelection(
                sites=[StereoSite.Atom(0), StereoSite.Bond(4)],
                kinds=[StereoKind.CisTrans],
                stereogenicities=[Stereogenicity.Stereogenic],
            )
        ),
    )
    assert (result if method == "remove_stereo" else stereo_determined) == expected


@pytest.mark.parametrize("method", ["remove_stereo", "remove_stereo_into"])
def test_molecule_remove_stereo_undetermined(stereo_determined, method):
    stereo_determined.stereo_bonds[0] = StereoBondForm.parse("Ct*")
    expected = Molecule.parse(
        """{
        :atoms ["C#D4#T1" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )
    result = getattr(stereo_determined, method)(
        StereoModel.default(),
        config=RemoveStereoConfig(mode=RemoveStereoMode.Undetermined),
    )
    assert (result if method == "remove_stereo" else stereo_determined) == expected


@pytest.mark.parametrize(
    "method,config",
    [
        (
            "add_stereo",
            AddStereoConfig(selection=StereoSelection(sites={StereoSite.Atom(9)})),
        ),
        (
            "add_stereo_into",
            AddStereoConfig(selection=StereoSelection(sites={StereoSite.Atom(9)})),
        ),
        (
            "remove_stereo",
            RemoveStereoConfig(selection=StereoSelection(sites={StereoSite.Atom(9)})),
        ),
        (
            "remove_stereo_into",
            RemoveStereoConfig(selection=StereoSelection(sites={StereoSite.Atom(9)})),
        ),
    ],
)
def test_molecule_stereo_error(stereo_determined, method, config):
    original = stereo_determined.copy()
    atom = stereo_determined.atoms[0]
    with pytest.raises(IndexError) as caught:
        getattr(stereo_determined, method)(StereoModel.default(), config=config)
    assert caught.value.site == StereoSite.Atom(9)
    assert (
        str(caught.value)
        == "stereo selection references unavailable site Atom(AtomId(9))"
    )
    assert stereo_determined == original
    if method.endswith("_into"):
        with pytest.raises(InvalidatedViewError):
            _ = atom.element
    else:
        assert atom.element == original.atoms[0].element


@pytest.mark.parametrize("method", ["add_stereo", "add_stereo_into"])
def test_molecule_add_stereo_evidence(stereo_determined, method):
    stereo_determined.stereo_atoms[0] = StereoAtomForm.parse("Th0")
    original = stereo_determined.copy()
    with pytest.raises(ContradictionError) as caught:
        getattr(stereo_determined, method)(StereoModel.default())
    assert caught.value.inconsistency == StereoInconsistency.TetrahedralStereoMismatch(
        atom=0, stereo_atom=0
    )
    assert stereo_determined == original


@pytest.mark.parametrize("method", ["remove_stereo", "remove_stereo_into"])
def test_molecule_remove_stereo_classification(stereo_determined, method):
    stereo_determined.stereo_atoms[0] = StereoAtomForm.parse("Sp0")
    original = stereo_determined.copy()
    with pytest.raises(ModelConversionError) as caught:
        getattr(stereo_determined, method)(
            StereoModel.default(),
            config=RemoveStereoConfig(
                selection=StereoSelection(stereogenicities={Stereogenicity.Symmetric})
            ),
        )
    assert caught.value.kind == StereoKind.SquarePlanar
    assert stereo_determined == original


@pytest.mark.parametrize("method", ["add_stereo_into", "remove_stereo_into"])
def test_molecule_stereo_argument_error(stereo_source, method):
    atom = stereo_source.atoms[0]
    original = stereo_source.copy()
    with pytest.raises(TypeError):
        getattr(stereo_source, method)(StereoModel.default(), config="default")
    assert atom.element == original.atoms[0].element
    assert stereo_source == original


def test_molecule_remove_stereo_into_sequence(stereo_source, stereo_determined):
    assert stereo_source.add_stereo_into(StereoModel.default()) is None
    assert stereo_source == stereo_determined
    with pytest.raises(IndexError):
        stereo_source.remove_stereo_into(
            StereoModel.default(),
            config=RemoveStereoConfig(
                selection=StereoSelection(sites={StereoSite.Bond(7)})
            ),
        )
    assert stereo_source == stereo_determined
