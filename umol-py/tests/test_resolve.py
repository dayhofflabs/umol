import pytest

from umol import (
    AromaticityConfig,
    AromaticityFailurePolicy,
    AromaticityModel,
    AromaticityRule,
    AromaticityResolveConfig,
    AtomCompletions,
    AtomForm,
    AtomSelectionAlgorithm,
    AtomTypeRegistry,
    ChemistryModel,
    ElementScope,
    ConsumedError,
    Edits,
    InvalidatedViewError,
    NumForm,
    IsotopePolicy,
    Molecule,
    MoleculeDefaults,
    ResolveConfig,
    ResolveReport,
    Solution,
    StereoFailurePolicy,
    StereoResolveConfig,
    UnderdeterminedError,
    ValenceModel,
    ValenceTable,
)


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (
            AromaticityFailurePolicy.Keep,
            AromaticityFailurePolicy.Keep,
            True,
        ),
    ],
)
def test_aromaticity_policy_equality(left, right, expected):
    assert (left == right) is expected


def test_aromaticity_policy_hash():
    policies = {
        AromaticityFailurePolicy.Keep: "failure",
    }

    assert policies[AromaticityFailurePolicy.Keep] == "failure"


@pytest.mark.parametrize(
    ("policy", "expected"),
    [
        (
            AromaticityFailurePolicy.Keep,
            "AromaticityFailurePolicy.Keep",
        ),
    ],
)
def test_aromaticity_policy_repr(policy, expected):
    assert repr(policy) == expected


@pytest.mark.parametrize(
    "policy",
    [
        AromaticityFailurePolicy.Keep,
    ],
)
def test_aromaticity_policy_mutation(policy):
    with pytest.raises(AttributeError):
        policy.value = "changed"


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (StereoFailurePolicy.Keep, StereoFailurePolicy.Keep, True),
        (StereoFailurePolicy.Remove, StereoFailurePolicy.Error, False),
    ],
)
def test_stereo_policy_equality(left, right, expected):
    assert (left == right) is expected


def test_stereo_policy_hash():
    policies = {
        StereoFailurePolicy.Remove: "failure",
    }

    assert policies[StereoFailurePolicy.Remove] == "failure"


@pytest.mark.parametrize(
    ("policy", "expected"),
    [
        (StereoFailurePolicy.Remove, "StereoFailurePolicy.Remove"),
    ],
)
def test_stereo_policy_repr(policy, expected):
    assert repr(policy) == expected


@pytest.mark.parametrize(
    "policy",
    [
        StereoFailurePolicy.Keep,
    ],
)
def test_stereo_policy_mutation(policy):
    with pytest.raises(AttributeError):
        policy.value = "changed"


@pytest.mark.parametrize(
    ("algorithm", "expected"),
    [
        (AtomSelectionAlgorithm.FullDfs, "AtomSelectionAlgorithm.FullDfs"),
        (AtomSelectionAlgorithm.BacktrackingDfs, "AtomSelectionAlgorithm.BacktrackingDfs"),
    ],
)
def test_atom_selection_algorithm_repr(algorithm, expected):
    assert repr(algorithm) == expected


def test_atom_selection_algorithm_hash():
    algorithms = {
        AtomSelectionAlgorithm.FullDfs: "full",
        AtomSelectionAlgorithm.BacktrackingDfs: "pruned",
    }

    assert algorithms[AtomSelectionAlgorithm.FullDfs] == "full"
    assert algorithms[AtomSelectionAlgorithm.BacktrackingDfs] == "pruned"


def test_aromaticity_resolve_config_default():
    config = AromaticityResolveConfig()

    assert config.perception == AromaticityConfig()
    assert config.atom_selection_algorithm == AtomSelectionAlgorithm.BacktrackingDfs
    assert config.aromatic_valence_failure == AromaticityFailurePolicy.Error
    assert config.aromatic_system_failure == AromaticityFailurePolicy.Error
    assert config == AromaticityResolveConfig()


@pytest.mark.parametrize(
    (
        "atom_selection_algorithm",
        "aromatic_valence_failure",
        "aromatic_system_failure",
    ),
    [
        (
            AtomSelectionAlgorithm.BacktrackingDfs,
            AromaticityFailurePolicy.Error,
            AromaticityFailurePolicy.Error,
        ),
        (
            AtomSelectionAlgorithm.FullDfs,
            AromaticityFailurePolicy.Keep,
            AromaticityFailurePolicy.Keep,
        ),
    ],
)
def test_aromaticity_resolve_config_new(
    atom_selection_algorithm,
    aromatic_valence_failure,
    aromatic_system_failure,
):
    perception = AromaticityConfig()
    config = AromaticityResolveConfig(
        perception=perception,
        atom_selection_algorithm=atom_selection_algorithm,
        aromatic_valence_failure=aromatic_valence_failure,
        aromatic_system_failure=aromatic_system_failure,
    )

    assert config.perception == perception
    assert config.perception is not perception
    assert config.atom_selection_algorithm == atom_selection_algorithm
    assert config.aromatic_valence_failure == aromatic_valence_failure
    assert config.aromatic_system_failure == aromatic_system_failure


def test_aromaticity_resolve_config_new_error():
    with pytest.raises(TypeError):
        AromaticityResolveConfig(False, True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            AromaticityResolveConfig(),
            "AromaticityResolveConfig(perception=AromaticityConfig("
            "ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "atom_selection_algorithm=AtomSelectionAlgorithm.BacktrackingDfs, "
            "aromatic_valence_failure=AromaticityFailurePolicy.Error, "
            "aromatic_system_failure=AromaticityFailurePolicy.Error)",
        ),
        (
            AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep,
                aromatic_system_failure=AromaticityFailurePolicy.Keep,
            ),
            "AromaticityResolveConfig(perception=AromaticityConfig("
            "ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "atom_selection_algorithm=AtomSelectionAlgorithm.BacktrackingDfs, "
            "aromatic_valence_failure=AromaticityFailurePolicy.Keep, "
            "aromatic_system_failure=AromaticityFailurePolicy.Keep)",
        ),
    ],
)
def test_aromaticity_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("perception", AromaticityConfig()),
        ("atom_selection_algorithm", AtomSelectionAlgorithm.FullDfs),
        ("aromatic_valence_failure", AromaticityFailurePolicy.Keep),
        ("aromatic_system_failure", AromaticityFailurePolicy.Keep),
    ],
)
def test_aromaticity_resolve_config_mutation(field, value):
    config = AromaticityResolveConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


def test_stereo_resolve_config_default():
    config = StereoResolveConfig()

    assert config.tetrahedral_stereo_failure == StereoFailurePolicy.Error
    assert config.stereo_atom_failure == StereoFailurePolicy.Error
    assert config.cis_trans_stereo_failure == StereoFailurePolicy.Error
    assert config.stereo_bond_failure == StereoFailurePolicy.Error
    assert config == StereoResolveConfig()


@pytest.mark.parametrize(
    (
        "tetrahedral_stereo_failure",
        "stereo_atom_failure",
        "cis_trans_stereo_failure",
        "stereo_bond_failure",
    ),
    [
        (
            StereoFailurePolicy.Error,
            StereoFailurePolicy.Error,
            StereoFailurePolicy.Error,
            StereoFailurePolicy.Error,
        ),
        (
            StereoFailurePolicy.Keep,
            StereoFailurePolicy.Remove,
            StereoFailurePolicy.Remove,
            StereoFailurePolicy.Keep,
        ),
    ],
)
def test_stereo_resolve_config_new(
    tetrahedral_stereo_failure,
    stereo_atom_failure,
    cis_trans_stereo_failure,
    stereo_bond_failure,
):
    config = StereoResolveConfig(
        tetrahedral_stereo_failure=tetrahedral_stereo_failure,
        stereo_atom_failure=stereo_atom_failure,
        cis_trans_stereo_failure=cis_trans_stereo_failure,
        stereo_bond_failure=stereo_bond_failure,
    )

    assert config.tetrahedral_stereo_failure == tetrahedral_stereo_failure
    assert config.stereo_atom_failure == stereo_atom_failure
    assert config.cis_trans_stereo_failure == cis_trans_stereo_failure
    assert config.stereo_bond_failure == stereo_bond_failure


def test_stereo_resolve_config_new_error():
    with pytest.raises(TypeError):
        StereoResolveConfig(True)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            StereoResolveConfig(),
            "StereoResolveConfig(tetrahedral_stereo_failure="
            "StereoFailurePolicy.Error, stereo_atom_failure="
            "StereoFailurePolicy.Error, cis_trans_stereo_failure="
            "StereoFailurePolicy.Error, stereo_bond_failure="
            "StereoFailurePolicy.Error)",
        ),
        (
            StereoResolveConfig(
                tetrahedral_stereo_failure=StereoFailurePolicy.Keep,
                stereo_atom_failure=StereoFailurePolicy.Remove,
                cis_trans_stereo_failure=StereoFailurePolicy.Remove,
                stereo_bond_failure=StereoFailurePolicy.Keep,
            ),
            "StereoResolveConfig(tetrahedral_stereo_failure="
            "StereoFailurePolicy.Keep, stereo_atom_failure="
            "StereoFailurePolicy.Remove, cis_trans_stereo_failure="
            "StereoFailurePolicy.Remove, stereo_bond_failure="
            "StereoFailurePolicy.Keep)",
        ),
    ],
)
def test_stereo_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("tetrahedral_stereo_failure", StereoFailurePolicy.Keep),
        ("stereo_atom_failure", StereoFailurePolicy.Remove),
        ("cis_trans_stereo_failure", StereoFailurePolicy.Remove),
        ("stereo_bond_failure", StereoFailurePolicy.Keep),
    ],
)
def test_stereo_resolve_config_mutation(field, value):
    config = StereoResolveConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


@pytest.mark.parametrize("isotope", [IsotopePolicy.Strict, IsotopePolicy.Natural])
def test_resolve_config_isotope(isotope):
    config = ResolveConfig(
        isotope=isotope,
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    assert config.isotope == isotope
    assert (config == ResolveConfig.default()) == (isotope == IsotopePolicy.Strict)
    assert {IsotopePolicy.Strict: "strict", IsotopePolicy.Natural: "natural"}[isotope] == (
        "strict" if isotope == IsotopePolicy.Strict else "natural"
    )


def test_resolve_config_default():
    config = ResolveConfig.default()

    assert config.isotope == IsotopePolicy.Strict
    assert config.aromaticity == AromaticityResolveConfig()
    assert config.stereo == StereoResolveConfig()
    assert config == ResolveConfig.default()


@pytest.mark.parametrize(
    ("aromaticity", "stereo"),
    [
        (
            AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep,
            ),
            StereoResolveConfig(),
        ),
        (
            AromaticityResolveConfig(),
            StereoResolveConfig(tetrahedral_stereo_failure=StereoFailurePolicy.Keep),
        ),
    ],
)
def test_resolve_config_new(aromaticity, stereo):
    config = ResolveConfig(aromaticity=aromaticity, stereo=stereo)

    assert config.aromaticity == aromaticity
    assert config.aromaticity is not aromaticity
    assert config.stereo == stereo
    assert config.stereo is not stereo


def test_resolve_config_new_error():
    with pytest.raises(TypeError):
        ResolveConfig(AromaticityResolveConfig(), StereoResolveConfig())


@pytest.mark.parametrize(
    "other",
    [
        ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                aromatic_system_failure=AromaticityFailurePolicy.Keep,
            ),
            stereo=StereoResolveConfig(),
        ),
        ResolveConfig(
            aromaticity=AromaticityResolveConfig(),
            stereo=StereoResolveConfig(
                tetrahedral_stereo_failure=StereoFailurePolicy.Keep,
            ),
        ),
    ],
)
def test_resolve_config_equality(other):
    assert ResolveConfig.default() != other


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (ResolveConfig.default(), "ResolveConfig.default()"),
        (
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    aromatic_valence_failure=AromaticityFailurePolicy.Keep,
                    aromatic_system_failure=AromaticityFailurePolicy.Keep,
                ),
                stereo=StereoResolveConfig(
                    tetrahedral_stereo_failure=StereoFailurePolicy.Keep,
                ),
            ),
            "ResolveConfig(isotope=IsotopePolicy.Strict, aromaticity=AromaticityResolveConfig(perception="
            "AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm="
            "SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()), "
            "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
            "maximum_independent_set_algorithm="
            "MaximumIndependentSetAlgorithm.BranchAndBound()), "
            "atom_selection_algorithm=AtomSelectionAlgorithm.BacktrackingDfs, "
            "aromatic_valence_failure=AromaticityFailurePolicy.Keep, "
            "aromatic_system_failure=AromaticityFailurePolicy.Keep), "
            "stereo=StereoResolveConfig(tetrahedral_stereo_failure="
            "StereoFailurePolicy.Keep, stereo_atom_failure="
            "StereoFailurePolicy.Error, cis_trans_stereo_failure="
            "StereoFailurePolicy.Error, stereo_bond_failure="
            "StereoFailurePolicy.Error))",
        ),
    ],
)
def test_resolve_config_repr(config, expected):
    assert repr(config) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("isotope", IsotopePolicy.Natural),
        ("aromaticity", AromaticityResolveConfig()),
        ("stereo", StereoResolveConfig()),
    ],
)
def test_resolve_config_mutation(field, value):
    config = ResolveConfig.default()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


@pytest.mark.parametrize(
    ("source", "valence"),
    [
        ("C", ValenceModel.counts(ValenceTable.default())),
        ("*", ValenceModel.smiles()),
    ],
)
def test_molecule_from_smiles_underdetermined(source, valence):
    default = ChemistryModel.default()
    model = ChemistryModel(
        connectivity=default.connectivity,
        valence=valence,
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )
    with pytest.raises(UnderdeterminedError) as excinfo:
        Molecule.from_smiles(source, chemistry_model=model)

    assert str(excinfo.value) == "resolution underdetermined"
    assert not hasattr(excinfo.value, "report")


@pytest.mark.parametrize(
    "source,expected",
    [
        pytest.param(
            '{:atoms ["C#T1" "F#n3" "Cl#n3" "Br#n3" "I#n3"] '
            ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
            ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}',
            "stereo inconsistency: tetrahedral stereo constraint at atom AtomId(0) "
            "disagrees with stereo atom StereoAtomId(0)",
            id="tetrahedral",
        ),
        pytest.param(
            '{:atoms ["C" "C" "F#n3" "Cl#n3" "Br#n3" "I#n3"] '
            ':bonds [[0 1 "2#C1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] '
            ':stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}',
            "stereo inconsistency: cis-trans stereo constraint at bond BondId(0) "
            "disagrees with stereo bond StereoBondId(0)",
            id="cis_trans",
        ),
    ],
)
@pytest.mark.parametrize("failure", [StereoFailurePolicy.Error, StereoFailurePolicy.Keep, StereoFailurePolicy.Remove])
def test_molecule_resolve_into_stereo_mismatch(source, expected, failure):
    molecule = Molecule.parse(source, defaults=MoleculeDefaults.concrete())
    original = molecule.copy()
    result = molecule.resolve_into(
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(),
            stereo=StereoResolveConfig(
                tetrahedral_stereo_failure=failure,
                stereo_atom_failure=failure,
                cis_trans_stereo_failure=failure,
                stereo_bond_failure=failure,
            ),
        ),
    )
    assert isinstance(result, Solution.Contradictory)
    assert str(result.contradiction) == expected
    assert molecule == original


@pytest.mark.parametrize("method", ["resolve_with_report", "resolve_into_with_report"])
def test_molecule_resolve_with_report(method):
    source = Molecule.parse('{:atoms ["C#i=#c0#u0#s#v0#a!"]}')
    default = ChemistryModel.default()
    result, report = getattr(source, method)(
        chemistry_model=ChemistryModel(
            connectivity=default.connectivity,
            valence=ValenceModel.counts(ValenceTable.default()),
            aromaticity=default.aromaticity,
            stereo=default.stereo,
        ),
    )
    assert isinstance(result, Solution.Underdetermined)
    assert isinstance(report, ResolveReport)
    assert report.tie_breaks == []
    completions = report.unresolved
    assert isinstance(completions, AtomCompletions)
    assert not completions.is_empty()
    assert len(completions) == 1
    assert completions.get(0) == [
        AtomForm.parse("C#i=#c0#h0#n2#u0#s#v0#a!"),
        AtomForm.parse("C#i=#c0#h2#n#u0#s#v0#a!"),
        AtomForm.parse("C#i=#c0#h4#n0#u0#s#v0#a!"),
    ]
    assert completions.get(1) is None
    assert completions.items() == [
        (
            0,
            [
                AtomForm.parse("C#i=#c0#h0#n2#u0#s#v0#a!"),
                AtomForm.parse("C#i=#c0#h2#n#u0#s#v0#a!"),
                AtomForm.parse("C#i=#c0#h4#n0#u0#s#v0#a!"),
            ],
        )
    ]
    assert repr(report) == (
        "ResolveReport(unresolved=AtomCompletions({0: "
        '["C#i=#c0#h0#n2#u0#s#v0#a!", "C#i=#c0#h2#n#u0#s#v0#a!", '
        '"C#i=#c0#h4#n0#u0#s#v0#a!"]}), tie_breaks=[])'
    )


@pytest.mark.parametrize("method", ["resolve_with_report", "resolve_into_with_report"])
def test_molecule_resolve_with_report_empty(method):
    source = Molecule.parse('{:atoms ["C#c0#h4"]}')
    result, report = getattr(source, method)(
        resolve_config=ResolveConfig(
            isotope=IsotopePolicy.Strict,
            aromaticity=AromaticityResolveConfig(),
            stereo=StereoResolveConfig(),
        )
    )
    assert isinstance(result, Solution.Underdetermined)
    assert report.unresolved.is_empty()
    assert len(report.unresolved) == 0
    assert report.unresolved.items() == []
    assert report.tie_breaks == []
    assert repr(report) == "ResolveReport(unresolved=AtomCompletions({}), tie_breaks=[])"


@pytest.mark.parametrize("valence", [ValenceModel.smiles(), ChemistryModel.default().valence])
@pytest.mark.parametrize("isotope", [IsotopePolicy.Strict, IsotopePolicy.Natural])
@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("[CH4]", '{:atoms ["C#i=#c0#h4#n0#u0#s"]}'),
        ("[13CH4]", '{:atoms ["C#i13#c0#h4#n0#u0#s"]}'),
    ],
)
def test_molecule_from_smiles_isotope(valence, isotope, source, expected):
    defaults = ChemistryModel.default()
    model = ChemistryModel(
        connectivity=defaults.connectivity,
        valence=valence,
        aromaticity=defaults.aromaticity,
        stereo=defaults.stereo,
    )
    config = ResolveConfig(
        isotope=isotope,
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    assert Molecule.from_smiles(
        source, chemistry_model=model, resolve_config=config
    ) == Molecule.parse(expected)


def test_molecule_from_smiles_isotope_default():
    expected = Molecule.parse('{:atoms ["C#i=#c0#h4#n0#u0#s"]}')
    assert Molecule.from_smiles("C", resolve_config=ResolveConfig.default()) == expected
    assert Molecule.from_smiles("C") == expected


@pytest.mark.parametrize(
    ("method", "consuming", "with_report"),
    [
        ("resolve", True, False),
        ("resolve_with_report", True, True),
        ("resolve_into", False, False),
        ("resolve_into_with_report", False, True),
    ],
)
@pytest.mark.parametrize("isotope", [IsotopePolicy.Strict, IsotopePolicy.Natural])
def test_molecule_resolve_isotope(method, consuming, with_report, isotope):
    source = Molecule.parse('{:atoms ["C#c0#h4"]}')
    alias = source
    original = source.copy()
    atoms = source.atoms
    atom = atoms[0]
    constraints = atom.constraints
    atom_iter = iter(atoms)
    molecule_constraints = source.constraints
    constraint_iter = iter(molecule_constraints)
    config = ResolveConfig(
        isotope=isotope,
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    result = getattr(source, method)(resolve_config=config)
    if with_report:
        result, report = result
        assert report.unresolved.items() == []
        assert report.tie_breaks == []
    assert not hasattr(result, "report")
    if isotope == IsotopePolicy.Strict:
        assert result == Solution.Underdetermined()
        if not consuming:
            assert alias == original
    else:
        expected = Molecule.parse('{:atoms ["C#i=#c0#h4#n0#u0#s"]}')
        assert isinstance(result, Solution.Determined)
        if consuming:
            assert result.molecule == expected
            assert result.molecule.edit().finish() == expected
            assert result.molecule == expected
        else:
            assert result == Solution.Determined()
            assert result.molecule is None
            assert alias == expected
            assert source.edit().finish() == expected
            assert source == expected
    if consuming:
        with pytest.raises(ConsumedError, match="^Molecule has been consumed$"):
            repr(alias)
    else:
        assert source.atoms[0].id == 0
    for access in (
        lambda: len(atoms), lambda: atom.id, lambda: atom.charge,
        lambda: len(constraints), lambda: next(atom_iter),
        lambda: len(molecule_constraints), lambda: next(constraint_iter),
    ):
        with pytest.raises(InvalidatedViewError):
            access()
    assert original == Molecule.parse('{:atoms ["C#c0#h4"]}')


@pytest.mark.parametrize(
    ("method", "consuming", "with_report"),
    [
        ("resolve", True, False),
        ("resolve_with_report", True, True),
        ("resolve_into", False, False),
        ("resolve_into_with_report", False, True),
    ],
)
def test_molecule_resolve_identity(method, consuming, with_report):
    source = Molecule.from_smiles("C")
    original = source.copy()
    atom = source.atoms[0]
    result = getattr(source, method)()
    if with_report:
        result, report = result
        assert report.unresolved.items() == []
        assert report.tie_breaks == []
    assert not hasattr(result, "report")
    if consuming:
        assert result.molecule == original
        with pytest.raises(ConsumedError, match="^Molecule has been consumed$"):
            repr(source)
    else:
        assert result == Solution.Determined()
        assert source == original
    with pytest.raises(InvalidatedViewError):
        atom.charge = 1


@pytest.mark.parametrize(
    ("method", "consuming", "with_report"),
    [
        ("resolve", True, False),
        ("resolve_with_report", True, True),
        ("resolve_into", False, False),
        ("resolve_into_with_report", False, True),
    ],
)
def test_molecule_resolve_contradiction(method, consuming, with_report):
    source = Molecule.parse('{:atoms ["C#c0#h5"]}')
    original = source.copy()
    atom = source.atoms[0]
    config = ResolveConfig(
        isotope=IsotopePolicy.Natural,
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    default = ChemistryModel.default()
    model = ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.smiles(),
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )
    result = getattr(source, method)(chemistry_model=model, resolve_config=config)
    if with_report:
        result, report = result
        assert report is None
    assert isinstance(result, Solution.Contradictory)
    assert str(result.contradiction) == "no matching valence state"
    assert not hasattr(result, "report")
    if consuming:
        with pytest.raises(ConsumedError, match="^Molecule has been consumed$"):
            repr(source)
    else:
        assert source == original
        assert source.atoms[0].implicit_hydrogens == NumForm.Lit(5)
    with pytest.raises(InvalidatedViewError):
        atom.id


@pytest.mark.parametrize(
    ("method", "consuming"),
    [
        ("resolve", True), ("resolve_with_report", True),
        ("resolve_into", False), ("resolve_into_with_report", False),
    ],
)
def test_molecule_resolve_error(method, consuming):
    source = Molecule.parse(
        '{:atoms ["C#v2#a2" "C#v2#a2" "C#v2#a2" '
        '"C#v2#a2" "C#v2#a2" "C#v2#a2"] '
        ':bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] '
        '[3 4 "1"] [4 5 "1"] [5 0 "1"]]}'
    )
    original = source.copy()
    atom = source.atoms[0]
    default = ChemistryModel.default()
    model = ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.atom_typing(AtomTypeRegistry.from_atoms([
            AtomForm.parse("C#c0#h0#n0#u0#s#v2#a2"),
        ])),
        aromaticity=AromaticityModel(
            scope=ElementScope.Any(),
            rule=AromaticityRule.Hmo(stabilization_threshold=0.5),
        ),
        stereo=default.stereo,
    )
    config = ResolveConfig(
        isotope=IsotopePolicy.Natural,
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    with pytest.raises(
        RuntimeError,
        match="^hmo: missing parameters: no Van-Catledge parameters for C with 2 pi-electrons$",
    ):
        getattr(source, method)(chemistry_model=model, resolve_config=config)
    if consuming:
        with pytest.raises(ConsumedError, match="^Molecule has been consumed$"):
            repr(source)
    else:
        assert source == original
        assert source.atoms[0].id == 0
    with pytest.raises(InvalidatedViewError):
        atom.id


@pytest.mark.parametrize(
    "method", ["resolve", "resolve_with_report", "resolve_into", "resolve_into_with_report"],
)
@pytest.mark.parametrize("argument", ["chemistry_model", "resolve_config"])
def test_molecule_resolve_arguments(method, argument):
    source = Molecule.from_smiles("C")
    original = source.copy()
    atom = source.atoms[0]
    with pytest.raises(TypeError):
        getattr(source, method)(**{argument: object()})
    assert source == original
    assert atom.id == 0


@pytest.mark.parametrize(
    "method", ["resolve", "resolve_with_report", "resolve_into", "resolve_into_with_report"],
)
def test_molecule_resolve_consumed(method):
    source = Molecule.from_smiles("C")
    source.resolve()
    with pytest.raises(ConsumedError, match="^Molecule has been consumed$"):
        getattr(source, method)()


@pytest.mark.parametrize(
    ("solution", "expected"),
    [
        (Solution.Determined(), "Solution.Determined(molecule=None)"),
        (Solution.Underdetermined(), "Solution.Underdetermined()"),
    ],
)
def test_solution_repr(solution, expected):
    assert repr(solution) == expected
    assert not hasattr(solution, "report")


@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply"])
def test_solution_molecule(method):
    source = Molecule.from_smiles("C")
    original = source.copy()
    result = source.resolve()
    molecule = result.molecule
    copied = molecule.copy()
    same = Solution.Determined(molecule=molecule)
    independent = Solution.Determined(molecule=copied)

    assert repr(result) == "Solution.Determined(molecule=Molecule(atoms=1, bonds=0))"
    assert result.molecule is molecule
    assert same.molecule is molecule
    assert result == independent
    molecule.atoms[0].charge = 1
    assert result.molecule.atoms[0].charge == NumForm.Lit(1)
    assert same == result
    assert independent != result
    assert independent.molecule == original

    view = molecule.atoms[0]
    if method == "edit":
        published = molecule.edit().finish()
    else:
        applied = getattr(molecule, method)(Edits())
        published = applied[0] if method == "tracked_apply" else applied
    assert published.atoms[0].charge == NumForm.Lit(1)
    assert result.molecule is molecule
    assert same == result
    assert result != independent
    assert independent != result
    assert repr(result) == "Solution.Determined(molecule=Molecule(atoms=1, bonds=0))"
    published.atoms[0].charge = -1
    assert view.charge == NumForm.Lit(1)
    assert same.molecule.atoms[0].charge == NumForm.Lit(1)
    assert published.atoms[0].charge == NumForm.Lit(-1)
    assert copied == original
