import sys

import pytest

from umol import (
    AtomConstraintForm,
    AutomorphismAlgorithm,
    ConsumedError,
    Edits,
    Entity,
    GraphSymmetry,
    GraphSymmetryConfig,
    GraphSymmetryView,
    InvalidatedViewError,
    LigandPermutation,
    Molecule,
    MoleculeDefaults,
    NumForm,
    Orientation,
    OrientedLigandPermutation,
    Permutation,
    StereoConfigurationForm,
    StereogenicSiteKind,
    Stereogenicity,
    Topicity,
)


@pytest.fixture
def meso():
    return Molecule.parse(
        '''{:atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                   [1 5 "1"] [1 6 "1"] [1 7 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                          {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]}''',
        defaults=MoleculeDefaults.concrete(),
    )


@pytest.fixture(
    params=[
        (
            ('{:atoms ["C" "H" "F" "Cl" "Br"] '
             ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
             ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}'),
            "stereo_atoms", (0, 1), Topicity.Distinct, [0, 1, 2, 3],
        ),
        (
            ('{:atoms ["C" "C" "F" "Cl" "F" "Cl"] '
             ':bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] '
             ':stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}'),
            "stereo_bonds", (0, 2), Topicity.Homotopic, [2, 3, 0, 1],
        ),
    ],
    ids=["atom", "bond"],
)
def stereo(request):
    source, namespace, pair, topicity, image = request.param
    molecule = Molecule.parse(source, defaults=MoleculeDefaults.concrete())
    return molecule, namespace, pair, topicity, image


@pytest.mark.parametrize("owner", [GraphSymmetry, GraphSymmetryView])
def test_graph_symmetry_new_error(owner):
    with pytest.raises(TypeError):
        owner()


def test_molecule_graph_symmetry(meso):
    references = sys.getrefcount(meso)
    view = meso.graph_symmetry(GraphSymmetryConfig())
    assert isinstance(view, GraphSymmetryView)
    assert sys.getrefcount(meso) == references + 1
    symmetry = view.into_graph_symmetry()
    assert isinstance(symmetry, GraphSymmetry)
    del view
    assert sys.getrefcount(meso) == references
    assert not meso.is_chiral(symmetry)
    del meso
    copied = symmetry.copy()
    assert isinstance(copied, GraphSymmetry)


@pytest.mark.parametrize("include_constraints", [False, True])
def test_graph_symmetry_view_from_graph_symmetry(meso, include_constraints):
    original = meso.copy()
    owned = meso.graph_symmetry(GraphSymmetryConfig(include_constraints=include_constraints))
    symmetry = owned.into_graph_symmetry()
    references = sys.getrefcount(symmetry)
    borrowed = GraphSymmetryView.from_graph_symmetry(symmetry, meso)
    assert sys.getrefcount(symmetry) == references + 1
    copied = borrowed.into_graph_symmetry()
    for data in (symmetry, copied, symmetry.copy()):
        view = GraphSymmetryView.from_graph_symmetry(data, meso)
        for operation, args, expected in [
            ("is_chiral", (), False),
            ("same_proper_orbit", (0, 0), True),
            ("same_proper_orbit", (0, 1), False),
            ("same_star_orbit", (0, 1), True),
            ("proper_orbit_of", (0,), [0]),
            ("star_orbit_of", (0,), [0, 1]),
            ("same_proper_orbit", (0, 8), False),
            ("same_star_orbit", (8, 0), False),
            ("proper_orbit_of", (8,), []),
            ("star_orbit_of", (8,), []),
        ]:
            assert getattr(view, operation)(*args) == expected
            assert getattr(meso, operation)(data, *args) == expected
    assert meso == original


def test_graph_symmetry_view_is_chiral(stereo):
    molecule, namespace, _, _, _ = stereo
    expected = namespace == "stereo_atoms"
    view = molecule.graph_symmetry(GraphSymmetryConfig())
    assert view.is_chiral() == expected
    data = view.into_graph_symmetry()
    assert molecule.is_chiral(data) == expected
    assert GraphSymmetryView.from_graph_symmetry(data, molecule).is_chiral() == expected


@pytest.mark.parametrize("borrowed", [False, True])
@pytest.mark.parametrize(
    "operation,args",
    [
        ("is_chiral", ()),
        ("same_proper_orbit", (0, 1)),
        ("same_star_orbit", (0, 1)),
        ("proper_orbit_of", (0,)),
        ("star_orbit_of", (0,)),
        ("stereogenic_site_kinds", ()),
        ("into_graph_symmetry", ()),
    ],
)
def test_graph_symmetry_view_into_graph_symmetry(meso, borrowed, operation, args):
    view = meso.graph_symmetry(GraphSymmetryConfig())
    if borrowed:
        source = view.into_graph_symmetry()
        view = GraphSymmetryView.from_graph_symmetry(source, meso)
    alias = view
    symmetry = view.into_graph_symmetry()
    with pytest.raises(ConsumedError, match="GraphSymmetryView"):
        getattr(alias, operation)(*args)
    assert not meso.is_chiral(symmetry)
    if borrowed:
        assert not meso.is_chiral(source)


@pytest.mark.parametrize("borrowed", [False, True])
@pytest.mark.parametrize("mutation", ["transact", "tracked_transact", "combine_from"])
def test_graph_symmetry_view_invalidation(meso, borrowed, mutation):
    view = meso.graph_symmetry(GraphSymmetryConfig())
    data = view.into_graph_symmetry()
    view = (GraphSymmetryView.from_graph_symmetry(data, meso) if borrowed else
            meso.graph_symmetry(GraphSymmetryConfig()))
    saved = view.star_orbit_of(0)
    if mutation == "combine_from":
        meso.combine_from(Molecule())
    else:
        getattr(meso, mutation)([Edits()])
    for operation, args in [
        ("is_chiral", ()), ("same_proper_orbit", (0, 1)),
        ("same_star_orbit", (0, 1)), ("proper_orbit_of", (0,)),
        ("star_orbit_of", (0,)), ("stereogenic_site_kinds", ()),
        ("into_graph_symmetry", ()),
    ]:
        with pytest.raises(InvalidatedViewError, match="GraphSymmetryView"):
            getattr(view, operation)(*args)
    assert saved == [0, 1]
    assert not meso.is_chiral(data)


def test_graph_symmetry_view_lifetime(meso):
    view = meso.graph_symmetry(GraphSymmetryConfig())
    del meso
    assert view.star_orbit_of(0) == [0, 1]
    data = view.into_graph_symmetry()
    del view
    assert isinstance(data.copy(), GraphSymmetry)


@pytest.mark.parametrize("include_constraints", [False, True])
def test_graph_symmetry_view_stereogenic_site_kinds(meso, include_constraints):
    view = meso.graph_symmetry(GraphSymmetryConfig(include_constraints=include_constraints))
    expected = [(Entity.StereoAtom(0), StereogenicSiteKind.Constitutional),
                (Entity.StereoAtom(1), StereogenicSiteKind.Constitutional)]
    result = view.stereogenic_site_kinds()
    assert result == expected
    result.clear()
    assert view.stereogenic_site_kinds() == expected
    data = view.into_graph_symmetry()
    assert meso.stereogenic_site_kinds(data) == expected
    assert GraphSymmetryView.from_graph_symmetry(data, meso).stereogenic_site_kinds() == expected


@pytest.mark.parametrize(
    "source,expected",
    [
        ('{:atoms ["C#h4"]}', []),
        ('{:atoms ["C" "H" "H" "F" "Cl"] '
         ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
         ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}',
         [(Entity.StereoAtom(0), StereogenicSiteKind.Prostereogenic)]),
        ('{:atoms ["C" "H" "H" "H" "H"] '
         ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
         ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}',
         [(Entity.StereoAtom(0), StereogenicSiteKind.Nonstereogenic)]),
        ('{:atoms ["C" "C" "F" "Cl" "F" "Cl"] '
         ':bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] '
         ':stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}',
         [(Entity.StereoBond(0), StereogenicSiteKind.Constitutional)]),
    ],
)
def test_molecule_stereogenic_site_kinds(source, expected):
    molecule = Molecule.parse(source, defaults=MoleculeDefaults.concrete())
    view = molecule.graph_symmetry(GraphSymmetryConfig())
    assert view.stereogenic_site_kinds() == expected
    data = view.into_graph_symmetry()
    assert molecule.stereogenic_site_kinds(data) == expected


def test_molecule_stereogenic_site_kinds_empty(stereo):
    molecule, namespace, _, _, _ = stereo
    getattr(molecule, namespace)[0].configuration = StereoConfigurationForm.Undetermined()
    assert molecule.graph_symmetry(GraphSymmetryConfig()).stereogenic_site_kinds() == []
    assert Molecule().graph_symmetry(GraphSymmetryConfig()).stereogenic_site_kinds() == []


@pytest.mark.parametrize("include_constraints", [False, True])
def test_stereo_view_with_symmetry(stereo, include_constraints):
    molecule, namespace, pair, topicity, image = stereo
    for atom in molecule.atoms:
        atom.charge = NumForm.LitSet({0})
    original = molecule.copy()
    data = molecule.graph_symmetry(GraphSymmetryConfig(include_constraints=include_constraints)).into_graph_symmetry()
    source = getattr(molecule, namespace)[0]
    references = sys.getrefcount(data)
    view = source.with_symmetry(data)
    assert sys.getrefcount(data) == references + 1
    permutation = OrientedLigandPermutation(LigandPermutation(Permutation(image)), Orientation.Proper)
    assert view.stereogenicity() == Stereogenicity.Stereogenic
    assert view.topicity(*pair) == topicity
    assert view.ligand_symmetry(permutation) is True
    assert view.ligand_symmetry(OrientedLigandPermutation(
        LigandPermutation(Permutation([1, 0, 2, 3])), Orientation.Proper)) is False
    assert view.ligand_symmetry(OrientedLigandPermutation(
        LigandPermutation(Permutation.identity(2)), Orientation.Proper)) is False
    for method, args in [("stereogenicity", ()), ("topicity", pair), ("ligand_symmetry", (permutation,))]:
        with pytest.raises(ValueError, match="requires symmetry context"):
            getattr(source, method)(*args)
    assert molecule == original
    del molecule, data, source
    assert view.stereogenicity() == Stereogenicity.Stereogenic


def test_stereo_view_derived_kind_error(stereo):
    molecule, namespace, pair, _, image = stereo
    getattr(molecule, namespace)[0].configuration = StereoConfigurationForm.Undetermined()
    data = molecule.graph_symmetry(GraphSymmetryConfig()).into_graph_symmetry()
    view = getattr(molecule, namespace)[0].with_symmetry(data)
    permutation = OrientedLigandPermutation(LigandPermutation(Permutation(image)), Orientation.Proper)
    for method, args in [("stereogenicity", ()), ("topicity", pair), ("ligand_symmetry", (permutation,))]:
        with pytest.raises(ValueError, match="requires a stereo kind"):
            getattr(view, method)(*args)


def test_stereo_view_topicity_error(stereo):
    molecule, namespace, _, _, _ = stereo
    data = molecule.graph_symmetry(GraphSymmetryConfig()).into_graph_symmetry()
    view = getattr(molecule, namespace)[0].with_symmetry(data)
    with pytest.raises(IndexError, match="ligand position"):
        view.topicity(0, 4)


@pytest.mark.parametrize("mutation", ["transact", "tracked_transact"])
def test_stereo_view_with_symmetry_invalidation(stereo, mutation):
    molecule, namespace, _, _, _ = stereo
    data = molecule.graph_symmetry(GraphSymmetryConfig()).into_graph_symmetry()
    view = getattr(molecule, namespace)[0].with_symmetry(data)
    constraints = view.constraints
    getattr(molecule, mutation)([Edits()])
    for method, args in [("with_symmetry", (data,)), ("stereogenicity", ())]:
        with pytest.raises(InvalidatedViewError):
            getattr(view, method)(*args)
    with pytest.raises(InvalidatedViewError):
        constraints.with_symmetry(data)


@pytest.mark.parametrize("include_constraints", [False, True])
def test_molecule_same_proper_orbit_constraints(include_constraints):
    molecule = Molecule.parse(
        '{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1"]]}',
        defaults=MoleculeDefaults.concrete(),
    )
    molecule.atoms[0].constraints.set(AtomConstraintForm.Degree(NumForm.Lit(1)))
    original = molecule.copy()
    view = molecule.graph_symmetry(GraphSymmetryConfig(include_constraints=include_constraints))
    assert view.same_proper_orbit(0, 1) == (not include_constraints)
    data = view.into_graph_symmetry()
    assert molecule.same_proper_orbit(data, 0, 1) == (not include_constraints)
    assert molecule == original


def test_molecule_graph_symmetry_raw():
    molecule = Molecule.parse(
        '{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1"]]}',
        defaults=MoleculeDefaults.concrete(),
    )
    molecule.atoms[0].charge = NumForm.LitSet({0})
    original = molecule.copy()
    view = molecule.graph_symmetry(GraphSymmetryConfig())
    assert view.same_proper_orbit(0, 1) is False
    assert view.proper_orbit_of(0) == [0]
    assert molecule == original


@pytest.mark.parametrize("kindless", [None, "stereo_atoms", "stereo_bonds"])
def test_molecule_stereogenic_site_kinds_mixed(kindless):
    molecule = Molecule.parse(
        '{:atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] '
        '        [5 6 "2"] [5 7 "1"] [6 8 "1"]] '
        ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}] '
        ':stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct0"}]}',
        defaults=MoleculeDefaults.concrete(),
    )
    expected = [(Entity.StereoAtom(0), StereogenicSiteKind.Constitutional),
                (Entity.StereoBond(0), StereogenicSiteKind.Constitutional)]
    if kindless:
        getattr(molecule, kindless)[0].configuration = StereoConfigurationForm.Undetermined()
        expected = [entry for entry in expected if entry[0] != (
            Entity.StereoAtom(0) if kindless == "stereo_atoms" else Entity.StereoBond(0))]
    view = molecule.graph_symmetry(GraphSymmetryConfig())
    assert view.stereogenic_site_kinds() == expected
    data = view.into_graph_symmetry()
    assert molecule.stereogenic_site_kinds(data) == expected


def test_graph_symmetry_config_default():
    config = GraphSymmetryConfig.default()
    assert config == GraphSymmetryConfig()
    assert config.automorphism_algorithm == AutomorphismAlgorithm.Nauty()
    assert config.max_classification_trials == 4096
    assert config.include_constraints is False


@pytest.mark.parametrize("include_constraints", [False, True])
@pytest.mark.parametrize("max_classification_trials", [0, 16])
def test_graph_symmetry_config_new(include_constraints, max_classification_trials):
    config = GraphSymmetryConfig(
        automorphism_algorithm=AutomorphismAlgorithm.Nauty(),
        max_classification_trials=max_classification_trials,
        include_constraints=include_constraints,
    )
    assert config.automorphism_algorithm == AutomorphismAlgorithm.Nauty()
    assert config.max_classification_trials == max_classification_trials
    assert config.include_constraints is include_constraints
    assert eval(repr(config), {
        "GraphSymmetryConfig": GraphSymmetryConfig,
        "AutomorphismAlgorithm": AutomorphismAlgorithm,
    }) == config


def test_graph_symmetry_config_copy():
    config = GraphSymmetryConfig()
    copied = config.copy()
    copied.automorphism_algorithm = AutomorphismAlgorithm.Nauty()
    copied.max_classification_trials = 16
    copied.include_constraints = True
    assert copied == GraphSymmetryConfig(max_classification_trials=16, include_constraints=True)
    assert config == GraphSymmetryConfig()
    with pytest.raises(TypeError):
        hash(config)
