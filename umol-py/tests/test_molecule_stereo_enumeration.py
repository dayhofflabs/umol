import pytest
import umol
from umol import (
    EnumerateStereoConfig,
    Molecule,
    MoleculeDefaults,
    StereoAtomForm,
    StereoConfigurationForm,
    StereoCoset,
    Stereogenicity,
    StereoKind,
    StereoModel,
    StereoSelection,
    StereoSite,
)


@pytest.fixture
def stereo_source():
    return Molecule.parse(
        """{
        :atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]
    }""",
        defaults=MoleculeDefaults.concrete(),
    )


@pytest.mark.parametrize("kwargs", [{}, {"config": EnumerateStereoConfig()}])
@pytest.mark.parametrize(
    "order,values",
    [
        ([StereoSite.Atom(0), StereoSite.Bond(4)], [(0, 0), (0, 1), (1, 0), (1, 1)]),
        ([StereoSite.Bond(4), StereoSite.Atom(0)], [(0, 0), (1, 0), (0, 1), (1, 1)]),
    ],
)
def test_molecule_enumerate_stereo(stereo_source, kwargs, order, values):
    original = stereo_source.copy()
    alias = stereo_source
    atom = stereo_source.atoms[0]
    constraints = atom.constraints
    iterator = iter(stereo_source.atoms)
    order = list(order)
    run = stereo_source.enumerate_stereo(StereoModel.default(), sites=order, **kwargs)
    order.clear()
    expected = [
        Molecule.parse(
            f"""{{
            :atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                    [5 6 "2"] [5 7 "1"] [6 8 "1"]]
            :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "Th{th}"}}]
            :stereo-bonds [{{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct{ct}"}}]
        }}""",
            defaults=MoleculeDefaults.concrete(),
        )
        for th, ct in values
    ]
    assert list(run) == expected
    assert stereo_source == original
    assert alias == original
    assert atom.element == original.atoms[0].element
    assert len(constraints) == 1
    assert next(iterator).element == original.atoms[0].element
    with pytest.raises(StopIteration):
        next(run)
    with pytest.raises(StopIteration):
        next(run)


@pytest.mark.parametrize(
    "selection,values",
    [
        (StereoSelection(sites=[StereoSite.Atom(0)]), [(0, "*"), (1, "*")]),
        (StereoSelection(kinds=[StereoKind.CisTrans]), [("*", 0), ("*", 1)]),
        (
            StereoSelection(stereogenicities=[Stereogenicity.Stereogenic]),
            [(0, 0), (0, 1), (1, 0), (1, 1)],
        ),
    ],
)
def test_molecule_enumerate_stereo_selection(stereo_source, selection, values):
    run = stereo_source.enumerate_stereo(
        StereoModel.default(),
        [
            StereoSite.Atom(99),
            StereoSite.Atom(0),
            StereoSite.Bond(4),
            StereoSite.Atom(99),
        ],
        config=EnumerateStereoConfig(selection=selection),
    )
    expected = [
        Molecule.parse(
            f"""{{
            :atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                    [5 6 "2"] [5 7 "1"] [6 8 "1"]]
            :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "Th{th}"}}]
            :stereo-bonds [{{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct{ct}"}}]
        }}""",
            defaults=MoleculeDefaults.concrete(),
        )
        for th, ct in values
    ]
    assert list(run) == expected


@pytest.mark.parametrize(
    "coset,values",
    [
        (StereoCoset.Lit(1), [1]),
        (StereoCoset.LitSet({0, 1}), [0, 1]),
        (StereoCoset.LitSet({1}), [1]),
    ],
)
def test_molecule_enumerate_stereo_domain(stereo_source, coset, values):
    stereo_source.stereo_atoms[0].configuration = StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, coset
    )
    run = stereo_source.enumerate_stereo(
        StereoModel.default(),
        [StereoSite.Atom(0)],
        config=EnumerateStereoConfig(
            selection=StereoSelection(sites=[StereoSite.Atom(0)])
        ),
    )
    expected = []
    for value in values:
        molecule = stereo_source.copy()
        molecule.stereo_atoms[0].configuration = StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(value)
        )
        expected.append(molecule)
    assert list(run) == expected


@pytest.mark.parametrize(
    "selection", [StereoSelection(sites=[]), StereoSelection(kinds=[])]
)
def test_molecule_enumerate_stereo_identity(stereo_source, selection):
    original = stereo_source.copy()
    atom = stereo_source.atoms[0]
    assert list(
        stereo_source.enumerate_stereo(
            StereoModel.default(), [], config=EnumerateStereoConfig(selection=selection)
        )
    ) == [original]
    assert stereo_source == original
    assert atom.element == original.atoms[0].element


def test_molecule_enumerate_stereo_multiplicity():
    source = Molecule.parse(
        """{:atoms ["C" "H" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}""",
        defaults=MoleculeDefaults.concrete(),
    )
    expected = [
        Molecule.parse(
            f"""{{:atoms ["C" "H" "H" "H" "H"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "Th{value}"}}]}}""",
            defaults=MoleculeDefaults.concrete(),
        )
        for value in (0, 1)
    ]
    assert (
        list(source.enumerate_stereo(StereoModel.default(), [StereoSite.Atom(0)]))
        == expected
    )


def test_molecule_enumerate_stereo_ownership(stereo_source):
    expected = stereo_source.copy()
    model = StereoModel.default()
    config = EnumerateStereoConfig()
    run = stereo_source.enumerate_stereo(
        model, [StereoSite.Atom(0), StereoSite.Bond(4)], config=config
    )
    assert iter(run) is run
    assert not hasattr(umol, "EnumerateStereoIter")
    assert not hasattr(umol._native, "EnumerateStereoIter")
    with pytest.raises(TypeError):
        type(run)()
    stereo_source.atoms[0].charge = 1
    del stereo_source, model, config
    expected.stereo_atoms[0].configuration = StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(0)
    )
    expected.stereo_bonds[0].configuration = StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(0)
    )
    first = next(run)
    assert first == expected
    first.atoms[0].charge = 2
    first.stereo_atoms[0].configuration = StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(1)
    )
    changed = first.copy()
    second = next(run)
    expected.stereo_bonds[0].configuration = StereoConfigurationForm.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(1)
    )
    assert second == expected
    assert first == changed
    del run
    assert first == changed
    assert second == expected


def test_molecule_enumerate_stereo_prefix():
    atoms = " ".join(['"C" "F" "Cl" "Br" "I"'] * 64)
    bonds = " ".join(
        f'[{site} {site + offset} "1"]'
        for site in range(0, 320, 5)
        for offset in range(1, 5)
    )
    stereo = " ".join(
        f'{{:site {site} :ligands [{site + 1} {site + 2} {site + 3} {site + 4}] :attrs "Th*"}}'
        for site in range(0, 320, 5)
    )
    source = Molecule.parse(
        f"{{:atoms [{atoms}] :bonds [{bonds}] :stereo-atoms [{stereo}]}}",
        defaults=MoleculeDefaults.concrete(),
    )
    run = source.enumerate_stereo(
        StereoModel.default(), [StereoSite.Atom(site) for site in range(0, 320, 5)]
    )
    first = next(run)
    second = next(run)
    expected = source.copy()
    for form in expected.stereo_atoms:
        form.configuration = StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.Lit(0)
        )
    assert first == expected
    expected.stereo_atoms[63].configuration = StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(1)
    )
    assert second == expected
    del run
    assert all(form.coset == StereoCoset.Undetermined() for form in source.stereo_atoms)


@pytest.mark.parametrize(
    "order,kwargs",
    [
        ([], {}),
        ([StereoSite.Atom(0), StereoSite.Atom(0)], {}),
        (
            [],
            {
                "config": EnumerateStereoConfig(
                    selection=StereoSelection(sites=[StereoSite.Atom(9)])
                )
            },
        ),
    ],
)
def test_molecule_enumerate_stereo_rejection(stereo_source, order, kwargs):
    original = stereo_source.copy()
    atom = stereo_source.atoms[0]
    run = stereo_source.enumerate_stereo(StereoModel.default(), order, **kwargs)
    with pytest.raises(StopIteration):
        next(run)
    with pytest.raises(StopIteration):
        next(run)
    assert stereo_source == original
    assert atom.element == original.atoms[0].element


@pytest.mark.parametrize(
    "configuration",
    [
        StereoConfigurationForm.Kinded(
            StereoKind.Tetrahedral, StereoCoset.LitSet(set())
        ),
        StereoConfigurationForm.Undetermined(),
        StereoConfigurationForm.Kinded(StereoKind.SquarePlanar, StereoCoset.Lit(0)),
        StereoAtomForm.parse("Th?configuration").configuration,
    ],
)
def test_molecule_enumerate_stereo_configuration_rejection(
    stereo_source, configuration
):
    stereo_source.stereo_atoms[0].configuration = configuration
    original = stereo_source.copy()
    atom = stereo_source.atoms[0]
    run = stereo_source.enumerate_stereo(
        StereoModel.default(),
        [StereoSite.Atom(0), StereoSite.Bond(4)],
        config=EnumerateStereoConfig(selection=StereoSelection()),
    )
    assert list(run) == []
    with pytest.raises(StopIteration):
        next(run)
    assert stereo_source == original
    assert atom.element == original.atoms[0].element


@pytest.mark.parametrize(
    "args,kwargs",
    [
        ((), {}),
        ((StereoModel.default(),), {}),
        (("default", []), {}),
        ((StereoModel.default(), [0]), {}),
        ((StereoModel.default(), []), {"config": "default"}),
    ],
)
def test_molecule_enumerate_stereo_argument_error(stereo_source, args, kwargs):
    original = stereo_source.copy()
    atom = stereo_source.atoms[0]
    with pytest.raises(TypeError):
        stereo_source.enumerate_stereo(*args, **kwargs)
    assert stereo_source == original
    assert atom.element == original.atoms[0].element
