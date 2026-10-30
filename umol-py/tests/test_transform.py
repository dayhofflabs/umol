import pytest
from umol import (
    AddStereoConfig,
    AtomFieldKind,
    AutomorphismAlgorithm,
    ConnectedComponentsAlgorithm,
    EnumerateStereoConfig,
    KekulizeConfig,
    MaximumMatchingAlgorithm,
    RelevantCycleEnumerationAlgorithm,
    RemoveStereoConfig,
    RemoveStereoMode,
    RingConfig,
    SimpleCycleEnumerationAlgorithm,
    StereoConfigurationForm,
    StereoCoset,
    Stereogenicity,
    StereoKind,
    StereoSelection,
    StereoSite,
    StereoTerm,
)


@pytest.mark.parametrize(
    "kwargs,expected_algorithm",
    [
        pytest.param({}, MaximumMatchingAlgorithm.Edmonds(), id="default"),
        pytest.param(
            {"algorithm": MaximumMatchingAlgorithm.Edmonds()},
            MaximumMatchingAlgorithm.Edmonds(),
            id="edmonds",
        ),
        pytest.param(
            {"algorithm": MaximumMatchingAlgorithm.HopcroftKarp()},
            MaximumMatchingAlgorithm.HopcroftKarp(),
            id="hopcroft-karp",
        ),
    ],
)
def test_kekulize_config_new(kwargs, expected_algorithm):
    config = KekulizeConfig(**kwargs)

    assert config.algorithm == expected_algorithm
    assert isinstance(config.algorithm, type(expected_algorithm))


@pytest.mark.parametrize(
    "args,kwargs",
    [
        pytest.param((MaximumMatchingAlgorithm.Edmonds(),), {}, id="positional"),
        pytest.param((), {"algorithm": "Edmonds"}, id="string"),
        pytest.param(
            (), {"algorithm": ConnectedComponentsAlgorithm.Bfs()}, id="wrong-enum"
        ),
    ],
)
def test_kekulize_config_new_error(args, kwargs):
    with pytest.raises(TypeError):
        KekulizeConfig(*args, **kwargs)


@pytest.mark.parametrize(
    "algorithm",
    [MaximumMatchingAlgorithm.Edmonds(), MaximumMatchingAlgorithm.HopcroftKarp()],
)
def test_kekulize_config_algorithm(algorithm):
    config = KekulizeConfig(algorithm=algorithm)

    with pytest.raises(AttributeError):
        config.algorithm = MaximumMatchingAlgorithm.HopcroftKarp()
    with pytest.raises(AttributeError):
        del config.algorithm
    assert config.algorithm == algorithm


@pytest.mark.parametrize(
    "algorithm,other,expected",
    [
        (MaximumMatchingAlgorithm.Edmonds(), MaximumMatchingAlgorithm.Edmonds(), True),
        (
            MaximumMatchingAlgorithm.HopcroftKarp(),
            MaximumMatchingAlgorithm.HopcroftKarp(),
            True,
        ),
        (
            MaximumMatchingAlgorithm.Edmonds(),
            MaximumMatchingAlgorithm.HopcroftKarp(),
            False,
        ),
    ],
)
def test_kekulize_config_eq(algorithm, other, expected):
    config = KekulizeConfig(algorithm=algorithm)
    other_config = KekulizeConfig(algorithm=other)

    assert (config == other_config) is expected
    assert (config != other_config) is not expected


@pytest.mark.parametrize(
    "algorithm,expected",
    [
        (
            MaximumMatchingAlgorithm.Edmonds(),
            "KekulizeConfig(algorithm=MaximumMatchingAlgorithm.Edmonds())",
        ),
        (
            MaximumMatchingAlgorithm.HopcroftKarp(),
            "KekulizeConfig(algorithm=MaximumMatchingAlgorithm.HopcroftKarp())",
        ),
    ],
)
def test_kekulize_config_repr(algorithm, expected):
    config = KekulizeConfig(algorithm=algorithm)

    assert repr(config) == expected
    assert (
        eval(
            repr(config),
            {
                "KekulizeConfig": KekulizeConfig,
                "MaximumMatchingAlgorithm": MaximumMatchingAlgorithm,
            },
        )
        == config
    )


@pytest.mark.parametrize(
    "site,index,expected",
    [
        (StereoSite.Atom(0), 0, "StereoSite.Atom(0)"),
        (StereoSite.Bond(0), 0, "StereoSite.Bond(0)"),
        (StereoSite.Atom(8), 8, "StereoSite.Atom(8)"),
    ],
)
def test_stereo_site(site, index, expected):
    assert repr(site) == expected
    assert site._0 == index
    assert eval(repr(site), globals()) == site
    assert hash(site) == hash(eval(repr(site), globals()))
    with pytest.raises(AttributeError):
        site._0 = 1
    assert {StereoSite.Atom(0), StereoSite.Bond(0)} == frozenset(
        [StereoSite.Atom(0), StereoSite.Bond(0)]
    )


@pytest.mark.parametrize(
    "constructor,value,error",
    [
        (StereoSite.Atom, -1, OverflowError),
        (StereoSite.Bond, 2**32, OverflowError),
        (StereoSite.Atom, "0", TypeError),
    ],
)
def test_stereo_site_new_error(constructor, value, error):
    with pytest.raises(error):
        constructor(value)


@pytest.mark.parametrize(
    "kwargs,sites,kinds,classes",
    [
        ({}, None, None, None),
        (
            {"sites": [], "kinds": [], "stereogenicities": []},
            frozenset(),
            frozenset(),
            frozenset(),
        ),
        (
            {
                "sites": [StereoSite.Atom(1), StereoSite.Bond(1)],
                "kinds": {StereoKind.Tetrahedral, StereoKind.CisTrans},
                "stereogenicities": {Stereogenicity.Prochiral},
            },
            frozenset({StereoSite.Atom(1), StereoSite.Bond(1)}),
            frozenset({StereoKind.Tetrahedral, StereoKind.CisTrans}),
            frozenset({Stereogenicity.Prochiral}),
        ),
    ],
)
def test_stereo_selection_new(kwargs, sites, kinds, classes):
    selection = StereoSelection(**kwargs)
    assert selection.sites == sites
    assert selection.kinds == kinds
    assert selection.stereogenicities == classes
    for field in ("sites", "kinds", "stereogenicities"):
        value = getattr(selection, field)
        assert value is None or type(value) is frozenset
        with pytest.raises(AttributeError):
            setattr(selection, field, None)
        with pytest.raises(AttributeError):
            delattr(selection, field)
    assert eval(repr(selection), globals()) == selection


@pytest.mark.parametrize(
    "kwargs,site,kind,stereogenicity,expected",
    [
        ({}, StereoSite.Atom(0), None, None, True),
        ({"sites": []}, StereoSite.Atom(0), None, None, False),
        ({"kinds": []}, StereoSite.Atom(0), None, None, False),
        ({"stereogenicities": []}, StereoSite.Atom(0), None, None, False),
        ({"kinds": [StereoKind.Tetrahedral]}, StereoSite.Atom(0), None, None, False),
        (
            {"stereogenicities": [Stereogenicity.Symmetric]},
            StereoSite.Atom(0),
            None,
            None,
            False,
        ),
        (
            {
                "sites": [StereoSite.Atom(0)],
                "kinds": [StereoKind.CisTrans],
                "stereogenicities": [Stereogenicity.Prochiral],
            },
            StereoSite.Atom(0),
            StereoKind.CisTrans,
            Stereogenicity.Prochiral,
            True,
        ),
        (
            {
                "sites": [StereoSite.Atom(0)],
                "kinds": [StereoKind.CisTrans],
                "stereogenicities": [Stereogenicity.Prochiral],
            },
            StereoSite.Bond(0),
            StereoKind.CisTrans,
            Stereogenicity.Prochiral,
            False,
        ),
        (
            {
                "sites": [StereoSite.Atom(0)],
                "kinds": [StereoKind.CisTrans],
                "stereogenicities": [Stereogenicity.Prochiral],
            },
            StereoSite.Atom(0),
            StereoKind.Tetrahedral,
            Stereogenicity.Prochiral,
            False,
        ),
        (
            {
                "sites": [StereoSite.Atom(0)],
                "kinds": [StereoKind.CisTrans],
                "stereogenicities": [Stereogenicity.Prochiral],
            },
            StereoSite.Atom(0),
            StereoKind.CisTrans,
            Stereogenicity.Stereogenic,
            False,
        ),
    ],
)
def test_stereo_selection_matches(kwargs, site, kind, stereogenicity, expected):
    selection = StereoSelection(**kwargs)
    assert selection.matches(site, kind, stereogenicity) is expected
    assert (
        StereoSelection(
            sites=selection.sites,
            kinds=selection.kinds,
            stereogenicities=selection.stereogenicities,
        )
        == selection
    )


@pytest.mark.parametrize(
    "kwargs,error",
    [
        ({"sites": [0]}, TypeError),
        ({"kinds": ["Th"]}, TypeError),
        ({"stereogenicities": [StereoKind.Tetrahedral]}, TypeError),
    ],
)
def test_stereo_selection_new_error(kwargs, error):
    with pytest.raises(error):
        StereoSelection(**kwargs)


@pytest.mark.parametrize(
    "configuration,undetermined",
    [
        (StereoConfigurationForm.Undetermined(), False),
        (
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
            True,
        ),
        (
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.Undetermined()
            ),
            True,
        ),
        (
            StereoConfigurationForm.Kinded(
                StereoKind.SquarePlanar, StereoCoset.Undetermined()
            ),
            True,
        ),
        (
            StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0)),
            False,
        ),
        (
            StereoConfigurationForm.Kinded(
                StereoKind.CisTrans, StereoCoset.LitSet({0, 1})
            ),
            False,
        ),
        (
            StereoConfigurationForm.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Term(StereoTerm.Var("c", None))
            ),
            False,
        ),
    ],
)
def test_remove_stereo_mode_matches(configuration, undetermined):
    assert RemoveStereoMode.All.matches(configuration) is True
    assert RemoveStereoMode.Undetermined.matches(configuration) is undetermined
    assert eval(repr(RemoveStereoMode.All), globals()) == RemoveStereoMode.All


@pytest.mark.parametrize(
    "constructor,kinds,classes",
    [
        (
            AddStereoConfig,
            frozenset({StereoKind.Tetrahedral, StereoKind.CisTrans}),
            frozenset({Stereogenicity.Stereogenic}),
        ),
        (RemoveStereoConfig, None, None),
        (
            EnumerateStereoConfig,
            frozenset({StereoKind.Tetrahedral, StereoKind.CisTrans}),
            None,
        ),
    ],
)
def test_stereo_config_new(constructor, kinds, classes):
    config = constructor()
    assert config.selection.sites is None
    assert config.selection.kinds == kinds
    assert config.selection.stereogenicities == classes
    assert config.automorphism_algorithm == AutomorphismAlgorithm.Nauty()
    assert (
        eval(
            repr(config),
            globals()
            | {
                "SimpleCycleEnumerationAlgorithm": SimpleCycleEnumerationAlgorithm,
                "RelevantCycleEnumerationAlgorithm": RelevantCycleEnumerationAlgorithm,
            },
        )
        == config
    )
    assert constructor(selection=None) == config
    unrestricted = constructor(selection=StereoSelection())
    assert unrestricted.selection == StereoSelection()
    restricted = constructor(selection=StereoSelection(kinds=[]))
    assert restricted.selection.kinds == frozenset()
    assert restricted != unrestricted
    if constructor is AddStereoConfig:
        assert config.ring_config == RingConfig()
    if constructor is RemoveStereoConfig:
        assert config.mode == RemoveStereoMode.All
        assert (
            constructor(mode=RemoveStereoMode.Undetermined).mode
            == RemoveStereoMode.Undetermined
        )


@pytest.mark.parametrize(
    "constructor,kwargs,error",
    [
        (AddStereoConfig, {"selection": []}, TypeError),
        (RemoveStereoConfig, {"mode": "All"}, TypeError),
        (
            EnumerateStereoConfig,
            {"automorphism_algorithm": MaximumMatchingAlgorithm.Edmonds()},
            TypeError,
        ),
        (AddStereoConfig, {"ring_config": {}}, TypeError),
        (EnumerateStereoConfig, {"kinds": []}, TypeError),
    ],
)
def test_stereo_config_new_error(constructor, kwargs, error):
    with pytest.raises(error):
        constructor(**kwargs)


@pytest.mark.parametrize(
    "constructor,kwargs",
    [
        (AddStereoConfig, {"ring_config": RingConfig()}),
        (RemoveStereoConfig, {"mode": RemoveStereoMode.Undetermined}),
        (EnumerateStereoConfig, {}),
    ],
)
def test_stereo_config_properties(constructor, kwargs):
    selection = StereoSelection(
        sites={StereoSite.Bond(1)},
        kinds={StereoKind.CisTrans},
        stereogenicities={Stereogenicity.Symmetric},
    )
    config = constructor(selection=selection, **kwargs)
    assert config.selection == selection
    for name in ("selection", "automorphism_algorithm", *kwargs):
        value = getattr(config, name)
        with pytest.raises(AttributeError):
            setattr(config, name, value)
        with pytest.raises(AttributeError):
            delattr(config, name)
        assert getattr(config, name) == value
    assert (
        eval(
            repr(config),
            globals()
            | {
                "SimpleCycleEnumerationAlgorithm": SimpleCycleEnumerationAlgorithm,
                "RelevantCycleEnumerationAlgorithm": RelevantCycleEnumerationAlgorithm,
            },
        )
        == config
    )


@pytest.mark.parametrize(
    "field",
    [
        AtomFieldKind.Element,
        AtomFieldKind.IsotopeMass,
        AtomFieldKind.Charge,
        AtomFieldKind.ImplicitHydrogens,
        AtomFieldKind.LonePairs,
        AtomFieldKind.UnpairedElectrons,
    ],
)
def test_atom_field_kind(field):
    assert eval(repr(field), globals()) == field
    assert hash(field) == hash(eval(repr(field), globals()))


@pytest.mark.parametrize(
    "values,expected",
    [
        (
            [StereoSite.Bond(0), StereoSite.Atom(2), StereoSite.Atom(0)],
            [StereoSite.Atom(0), StereoSite.Atom(2), StereoSite.Bond(0)],
        ),
        (
            [AtomFieldKind.Charge, AtomFieldKind.Element, AtomFieldKind.IsotopeMass],
            [AtomFieldKind.Element, AtomFieldKind.IsotopeMass, AtomFieldKind.Charge],
        ),
    ],
)
def test_stereo_value_order(values, expected):
    assert sorted(values) == expected
    assert expected[0] < expected[1]
    assert expected[1] <= expected[1]
    assert expected[2] > expected[1]
    assert expected[2] >= expected[2]
    assert expected[0] != expected[1]
    assert expected[0] != 0
