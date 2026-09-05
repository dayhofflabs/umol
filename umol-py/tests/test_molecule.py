from umol import MoleculeRemapping
import re

import pytest

from umol import (
    AromaticSystemForm,
    AromaticValenceForm,
    AromaticityConfig,
    AromaticityFailurePolicy,
    AromaticityModel,
    AromaticityRule,
    AromaticityResolveConfig,
    AutomorphismAlgorithm,
    AtomForm,
    AtomTypeRegistry,
    BondForm,
    CanonicalizeConfig,
    ChemistryModel,
    ConnectedComponentsAlgorithm,
    Constraint,
    ContradictionError,
    Correspondence,
    DativeBondForm,
    ElectronCountsForm,
    Element,
    ElementForm,
    ElementScope,
    Entity,
    InvalidStructureError,
    MaximumIndependentSetAlgorithm,
    MetadataError,
    ModelConversionError,
    Molecule,
    MoleculeConstraint,
    MoleculeCorrespondence,
    MoleculeDefaults,
    MoleculeMetadata,
    MulticenterBondForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    ParseError,
    ResolveConfig,
    RelevantCycleEnumerationAlgorithm,
    RingConfig,
    RingLimits,
    SimpleCycleEnumerationAlgorithm,
    SmilesIoConfig,
    Solution,
    SmilesSyntaxFlags,
    StereoAtomForm,
    StereoBondForm,
    StereoConfigurationForm,
    StereoCoset,
    StereoFailurePolicy,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereoModel,
    StereoResolveConfig,
    TetrahedralConfiguration,
    UnderdeterminedError,
    ValenceCandidateSource,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
    ValenceTieBreak,
    NumForm,
)


def test_canonicalization_config():
    config = CanonicalizeConfig()

    assert config == CanonicalizeConfig.default()
    assert config.automorphism_algorithm == AutomorphismAlgorithm.Nauty()
    assert repr(config) == "CanonicalizeConfig.default()"
    with pytest.raises(AttributeError):
        config.automorphism_algorithm = AutomorphismAlgorithm.Nauty()


def test_molecule_new():
    assert len(Molecule().atoms) == 0
    assert len(Molecule().bonds) == 0


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (MoleculeDefaults(), MoleculeDefaults(), "MoleculeDefaults()"),
        (
            MoleculeDefaults.concrete(),
            MoleculeDefaults.concrete(),
            "MoleculeDefaults.concrete()",
        ),
    ],
)
def test_molecule_defaults(value, expected, expected_repr):
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize(
    ("source", "defaults", "expected"),
    [
        (
            '{:atoms ["C"]}',
            MoleculeDefaults(),
            Molecule.from_entries([AtomForm.parse("C")]),
        ),
        (
            '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
            MoleculeDefaults.concrete(),
            Molecule.from_entries(
                [AtomForm.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
            ),
        ),
    ],
)
def test_molecule_parse(source, defaults, expected):
    assert Molecule.parse(source, defaults=defaults) == expected


def test_molecule_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        Molecule.parse("not edn")


def test_molecule_parse_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Molecule.parse\\(\\) takes 1 positional arguments but 2 were given$"
        ),
    ):
        Molecule.parse('{:atoms ["C"]}', MoleculeDefaults.concrete())


def test_molecule_parse_with_metadata():
    source = (
        '{:atom-aliases [:x "C"] :atoms [[:carbon :x]] :bonds []}'
    )

    molecule, metadata = Molecule.parse_with_metadata(source)

    assert molecule == Molecule.from_entries([AtomForm(Element("C"))])
    assert metadata.keyword(Entity.Atom(0)) == "carbon"
    assert metadata.entity("carbon") == Entity.Atom(0)
    assert repr(metadata) == (
        'MoleculeMetadata(keywords=[(Entity.Atom(0), "carbon")], '
        "atom_alias_count=1)"
    )
    assert molecule.render_with_metadata(metadata) == source


def test_molecule_parse_with_metadata_defaults():
    molecule, metadata = Molecule.parse_with_metadata(
        '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.concrete(),
    )

    assert molecule == Molecule.from_entries(
        [
            AtomForm.parse(
                "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"
            )
        ]
    )
    assert metadata == MoleculeMetadata()


def test_molecule_parse_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Molecule.parse_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        Molecule.parse_with_metadata(
            '{:atoms ["C"]}',
            MoleculeDefaults.concrete(),
        )


@pytest.mark.parametrize(
    ("molecule", "defaults", "expected"),
    [
        (
            Molecule.parse('{:atoms ["C"]}'),
            MoleculeDefaults(),
            '{:atoms ["C"] :bonds []}',
        ),
        (
            Molecule.parse(
                '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
                defaults=MoleculeDefaults.concrete(),
            ),
            MoleculeDefaults.concrete(),
            '{:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}',
        ),
    ],
)
def test_molecule_render(molecule, defaults, expected):
    assert molecule.render(defaults=defaults) == expected


def test_molecule_render_keyword_error():
    with pytest.raises(
        TypeError,
        match="^Molecule.render\\(\\) takes 0 positional arguments but 1 was given$",
    ):
        Molecule.parse('{:atoms ["C"]}').render(MoleculeDefaults())


def test_molecule_render_with_metadata():
    source = (
        '{:atom-aliases [:x "C"] :atoms [[:carbon :x]] :bonds []}'
    )
    molecule, metadata = Molecule.parse_with_metadata(source)

    assert molecule.render_with_metadata(metadata) == source
    assert molecule.render() == '{:atoms ["C"] :bonds []}'


def test_molecule_render_with_metadata_error():
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(1), "outside")

    with pytest.raises(
        MetadataError,
        match="^metadata entity is out of range: atom 1$",
    ):
        Molecule.parse('{:atoms ["C"]}').render_with_metadata(metadata)


def test_molecule_render_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Molecule.render_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        Molecule.parse('{:atoms ["C"]}').render_with_metadata(
            MoleculeMetadata(),
            MoleculeDefaults(),
        )


def test_molecule_str():
    molecule = Molecule.parse(
        '{:atoms ["C" "O"] :bonds [[0 1 "1"]]}'
    )

    assert str(molecule) == molecule.render()


def test_molecule_from_entries():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(5)],
        bonds=[
            (0, 1, BondForm(2)),
            (0, 2, BondForm(1)),
            (0, 3, BondForm(1)),
            (0, 4, BondForm(1)),
            (1, 3, BondForm(1)),
        ],
        dative_bonds=[([2], 1, DativeBondForm(1))],
        aromatic_systems=[([0, 1, 2], AromaticSystemForm([1, 1, 1]))],
        multicenter_bonds=[([0, 1, 2], MulticenterBondForm([1, 1, 1]))],
        noncovalent_bonds=[
            ([0, 2], NoncovalentBondForm(NoncovalentBondKind.HydrogenBond))
        ],
        stereo_atoms=[
            (
                0,
                [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
                StereoAtomForm(TetrahedralConfiguration.Ccw),
            )
        ],
        stereo_bonds=[
            (
                0,
                [
                    StereoLigand(2, StereoLigandKind.Atom),
                    StereoLigand(0, StereoLigandKind.ImplicitHydrogen),
                    StereoLigand(3, StereoLigandKind.Atom),
                    StereoLigand(1, StereoLigandKind.ImplicitHydrogen),
                ],
                StereoBondForm.parse("Ct0"),
            )
        ],
        constraints=[
            Constraint.Molecule(MoleculeConstraint.Connected([0, 1, 2, 3, 4]))
        ],
    )

    assert len(molecule.atoms) == 5
    assert len(molecule.bonds) == 5
    assert len(molecule.dative_bonds) == 1
    assert len(molecule.aromatic_systems) == 1
    assert len(molecule.multicenter_bonds) == 1
    assert len(molecule.noncovalent_bonds) == 1
    assert len(molecule.stereo_atoms) == 1
    assert len(molecule.stereo_bonds) == 1
    assert list(molecule.constraints) == [
        Constraint.Molecule(MoleculeConstraint.Connected([0, 1, 2, 3, 4]))
    ]


def test_molecule_from_entries_default():
    molecule = Molecule.from_entries([AtomForm(Element("C"))])

    assert len(molecule.atoms) == 1
    assert len(molecule.bonds) == 0


def test_molecule_from_entries_atom_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule references unavailable atom 1$",
    ):
        Molecule.from_entries(
            [AtomForm(Element("C"))],
            bonds=[(0, 1, BondForm(1))],
        )


def test_molecule_from_entries_bond_site_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule references unavailable bond 0$",
    ):
        Molecule.from_entries(
            [AtomForm(Element("C"))],
            stereo_bonds=[
                (
                    0,
                    [StereoLigand(0, StereoLigandKind.Atom)],
                    StereoBondForm.parse("Ct0"),
                )
            ],
        )


def test_molecule_from_entries_ligand_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule references unavailable atom 1$",
    ):
        Molecule.from_entries(
            [AtomForm(Element("C"))],
            stereo_atoms=[
                (
                    0,
                    [StereoLigand(1, StereoLigandKind.Atom)],
                    StereoAtomForm(TetrahedralConfiguration.Ccw),
                )
            ],
        )


def test_molecule_from_entries_constraint_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule references unavailable atom 1$",
    ):
        Molecule.from_entries(
            [AtomForm(Element("C"))],
            constraints=[Constraint.Molecule(MoleculeConstraint.Connected([1]))],
        )


def test_molecule_canonicalize():
    source = Molecule.parse(
        '{:atoms ["C#c+" "C"] :constraints '
        '[{:charge-sum {:atoms [0 0] :sum 0}}]}'
    )
    expected = Molecule.parse(
        '{:atoms ["C" "C#c+"] :constraints '
        '[{:charge-sum {:atoms [1] :sum 0}}]}'
    )

    canonical = source.canonicalize()

    assert canonical is not source
    assert canonical == expected
    assert source != expected
    assert source.canonical_eq(expected)


def test_molecule_canonicalize_with_remapping():
    source = Molecule.parse(
        '{:atoms ["C#c+" "C"] :constraints '
        '[{:charge-sum {:atoms [0 0] :sum 0}}]}'
    )

    canonical, remapping = source.canonicalize_with_remapping()

    assert canonical == source.canonicalize()
    assert isinstance(remapping, MoleculeRemapping)
    assert remapping.to_correspondence().is_total()
    assert remapping.atoms.images == [1, 0]


def test_molecule_canonicalize_error():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C"), charge=NumForm.LitSet(set()))]
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        molecule.canonicalize()
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        molecule.canonicalize_with_remapping()


def test_molecule_stereo_mutation_integrity_error():
    molecule = Molecule.parse(
        '{:atoms ["C" "F" "Cl" "Br" "I"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
        ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}'
    )
    with pytest.raises(InvalidStructureError, match="ligands"):
        molecule.stereo_atoms[0].configuration = (
            StereoConfigurationForm.Kinded(
                StereoKind.Octahedral,
                StereoCoset.Lit(0),
            )
        )

    assert molecule.stereo_atoms[0].configuration == StereoConfigurationForm.Kinded(
        StereoKind.Tetrahedral,
        StereoCoset.Lit(0),
    )
    molecule.canonicalize()
    molecule.canonicalize_with_remapping()


def test_molecule_from_smiles():
    assert Molecule.from_smiles("C") == Molecule.parse('{:atoms ["C#i=#c0#h4#n0#u0#s"]}')


@pytest.mark.parametrize(
    "io_config",
    [
        SmilesIoConfig.with_syntax_flags(
            syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS
        ),
        SmilesIoConfig.lenient(),
    ],
)
def test_molecule_from_smiles_io_config(io_config):
    molecule = Molecule.from_smiles("[se]1cccc1", io_config=io_config)

    assert [atom.element for atom in molecule.atoms] == [
        ElementForm.Lit(Element("Se")),
        *[ElementForm.Lit(Element("C")) for _ in range(4)],
    ]
    assert [(bond.atom_ids, bond.order) for bond in molecule.bonds] == [
        ((0, 4), NumForm.Lit(1)),
        ((0, 1), NumForm.Lit(1)),
        ((1, 2), NumForm.Lit(1)),
        ((2, 3), NumForm.Lit(1)),
        ((3, 4), NumForm.Lit(1)),
    ]
    assert [
        (system.atom_ids, system.electrons, system.charge)
        for system in molecule.aromatic_systems
    ] == [
        (
            (0, 1, 2, 3, 4),
            ElectronCountsForm.Lit([2, 1, 1, 1, 1]),
            NumForm.Lit(0),
        )
    ]


def test_molecule_from_smiles_io_config_dative():
    assert Molecule.from_smiles(
        "C->N", io_config=SmilesIoConfig.lenient()
    ) == Molecule.parse(
        '{:atoms ["C#i=#c0#h4#n0#u0#s" "N#i=#c0#h#n2#u0#s"] '
        ":dative-bonds [{:acceptor 1 :attrs :single :donors [0]}]}"
    )


def test_molecule_from_smiles_io_config_error():
    with pytest.raises(ParseError, match="^Invalid token at position 2$"):
        Molecule.from_smiles("C->N", io_config=SmilesIoConfig.opensmiles())


@pytest.mark.parametrize(
    ("valence_model", "expected"),
    [
        (
            ValenceModel.atom_typing(
                AtomTypeRegistry.from_atoms(
                    [AtomForm.parse("C#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
                )
            ),
            AtomForm.parse("C#i=#c0#h4#n0#u0#s"),
        ),
        (
            ValenceModel(
                candidates=ValenceCandidateSource.Counts(
                    table=ValenceTable.default()
                ),
                tie_break=ValenceTieBreak.MostSaturated,
            ),
            AtomForm.parse("C#i=#c0#h4#n0#u0#s"),
        ),
    ],
)
def test_molecule_from_smiles_chemistry_model_valence(valence_model, expected):
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        connectivity=ChemistryModel.default().connectivity,
        valence=valence_model,
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )

    assert Molecule.from_smiles(
        "C", chemistry_model=chemistry_model
    ) == Molecule.from_entries([expected])


def test_molecule_from_smiles_chemistry_model_aromaticity():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        connectivity=ChemistryModel.default().connectivity,
        valence=default.valence,
        aromaticity=AromaticityModel(
            scope=ElementScope.Any(),
            rule=AromaticityRule.Hmo(stabilization_threshold=0.375),
        ),
        stereo=default.stereo,
    )

    molecule = Molecule.from_smiles(
        "c1ccccc1",
        chemistry_model=chemistry_model,
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep
            ),
            stereo=StereoResolveConfig(),
        ),
    )

    assert [atom.implicit_hydrogens for atom in molecule.atoms] == [NumForm.Lit(1)] * 6
    assert [atom.constraints.aromatic_valence for atom in molecule.atoms] == [
        AromaticValenceForm.Aromatic(NumForm.Undetermined())
    ] * 6
    assert list(molecule.aromatic_systems) == []


def _smiles_valence_model():
    default = ChemistryModel.default()
    return ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.smiles(),
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )


def _counts_strict_model():
    default = ChemistryModel.default()
    return ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.counts(ValenceTable.default()),
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )


def test_molecule_resolve():
    molecule = Molecule.parse('{:atoms ["C#c0"]}')

    solution = molecule.resolve(chemistry_model=_smiles_valence_model())

    assert isinstance(solution, Solution.Determined)
    assert solution.molecule == Molecule.parse(
        '{:atoms ["C#i=#c0#h4#n0#u0#s"]}'
    )
    assert solution.report.tie_breaks == [0]
    assert molecule == Molecule.parse('{:atoms ["C#c0"]}')


def test_molecule_resolve_underdetermined():
    molecule = Molecule.parse('{:atoms ["C#c0"]}')

    solution = molecule.resolve(chemistry_model=_counts_strict_model())

    assert isinstance(solution, Solution.Underdetermined)
    assert len(solution.report.unresolved.get(0)) == 5
    assert molecule == Molecule.parse('{:atoms ["C#c0"]}')


def test_molecule_resolve_contradiction():
    molecule = Molecule.parse('{:atoms ["C#c0#h5"]}')

    solution = molecule.resolve(chemistry_model=_smiles_valence_model())

    assert isinstance(solution, Solution.Contradictory)
    assert str(solution.contradiction) == "no matching valence state"
    assert molecule == Molecule.parse('{:atoms ["C#c0#h5"]}')


def test_molecule_resolve_default_model():
    solution = Molecule.parse('{:atoms ["C"]}').resolve()

    assert isinstance(solution, Solution.Underdetermined)
    assert len(solution.report.unresolved.get(0)) == 9




@pytest.mark.parametrize(
    ("source", "expected"),
    [
        pytest.param(
            "o1cccc1",
            '{:atoms ["O#i=#c0#h0#n#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="furan",
        ),
        pytest.param(
            "s1cccc1",
            '{:atoms ["S#i=#c0#h0#n#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="thiophene",
        ),
        pytest.param(
            "[nH]1cccc1",
            '{:atoms ["N#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+" '
            '"C#i=#c0#h#n0#u0#s#a+"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="pyrrole",
        ),
    ],
)
def test_molecule_from_smiles_aromaticity_policy(source, expected):
    default = ChemistryModel.default()

    assert Molecule.from_smiles(
        source,
        chemistry_model=ChemistryModel(
            connectivity=ChemistryModel.default().connectivity,
            valence=default.valence,
            aromaticity=AromaticityModel.mdl(),
            stereo=default.stereo,
        ),
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep
            ),
            stereo=StereoResolveConfig(),
        ),
    ) == Molecule.parse(expected)


def test_molecule_from_smiles_chemistry_model_stereo():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        connectivity=ChemistryModel.default().connectivity,
        valence=ValenceModel.smiles(),
        aromaticity=default.aromaticity,
        stereo=StereoModel(
            kind_models={},
            para_stereo=False,
        ),
    )

    molecule = Molecule.from_smiles(
        "C[C@H](N)O",
        chemistry_model=chemistry_model,
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(),
            stereo=StereoResolveConfig(
                tetrahedral_stereo_failure=StereoFailurePolicy.Remove
            ),
        ),
    )

    assert [atom.implicit_hydrogens for atom in molecule.atoms] == [
        NumForm.Lit(3),
        NumForm.Lit(1),
        NumForm.Lit(2),
        NumForm.Lit(1),
    ]
    assert [atom.constraints.tetrahedral_stereo for atom in molecule.atoms] == [
        None
    ] * 4
    assert list(molecule.stereo_atoms) == []


@pytest.mark.parametrize(
    ("source", "resolve_config", "expected"),
    [
        (
            "[cH+]1[cH][cH]1",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    perception=AromaticityConfig(
                        ring_config=RingConfig(
                            simple_cycle_algorithm=(
                                SimpleCycleEnumerationAlgorithm.ReadTarjan()
                            ),
                            relevant_cycle_algorithm=(
                                RelevantCycleEnumerationAlgorithm.Vismara()
                            ),
                        ),
                        connected_components_algorithm=(
                            ConnectedComponentsAlgorithm.Bfs()
                        ),
                        maximum_independent_set_algorithm=(
                            MaximumIndependentSetAlgorithm.BranchAndBound()
                        ),
                    ),
                    reset_aromatic_valence=False,
                ),
                stereo=StereoResolveConfig(),
            ),
            (
                [NumForm.Lit(1), NumForm.Lit(0), NumForm.Lit(0)],
                [None] * 3,
                [None] * 3,
                [((0, 1, 2), NumForm.Lit(0))],
                [],
            ),
        ),
        (
            "c1ccccc1",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    reset_aromatic_valence=True,
                ),
                stereo=StereoResolveConfig(),
            ),
            (
                [NumForm.Lit(0)] * 6,
                [None] * 6,
                [None] * 6,
                [((0, 1, 2, 3, 4, 5), NumForm.Lit(0))],
                [],
            ),
        ),
        (
            "C[C@H](N)O",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(),
                stereo=StereoResolveConfig(reset_stereo_constraints=True),
            ),
            (
                [NumForm.Lit(0)] * 4,
                [None] * 4,
                [None] * 4,
                [],
                [(1, StereoKind.Tetrahedral, StereoCoset.Lit(0))],
            ),
        ),
    ],
)
def test_molecule_from_smiles_resolve_config(source, resolve_config, expected):
    molecule = Molecule.from_smiles(source, resolve_config=resolve_config)

    assert (
        [atom.charge for atom in molecule.atoms],
        [atom.constraints.aromatic_valence for atom in molecule.atoms],
        [atom.constraints.tetrahedral_stereo for atom in molecule.atoms],
        [(system.atom_ids, system.charge) for system in molecule.aromatic_systems],
        [
            (stereo.site_id, stereo.kind, stereo.coset)
            for stereo in molecule.stereo_atoms
        ],
    ) == expected


@pytest.mark.parametrize(
    ("source", "kwargs", "error_type", "message"),
    [
        (" C", {}, ParseError, "Leading whitespace"),
        (
            "C[S@]C",
            {},
            ModelConversionError,
            "tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands",
        ),
        (
            "[nH]1cccc1",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel(
                        scope=ElementScope.Any(),
                        rule=AromaticityRule.Clar(),
                    ),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "clar: non-benzenoid input: Clar model requires benzenoid input "
            "but non-carbon aromatic atoms are present",
        ),
        pytest.param(
            "o1cccc1",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "aromaticity inconsistency: aromatic valence at atom AtomId(0) "
            "cannot produce a valid aromatic system",
            id="mdl-furan",
        ),
        pytest.param(
            "s1cccc1",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "aromaticity inconsistency: aromatic valence at atom AtomId(0) "
            "cannot produce a valid aromatic system",
            id="mdl-thiophene",
        ),
        pytest.param(
            "[nH]1cccc1",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "aromaticity inconsistency: aromatic valence at atom AtomId(0) "
            "cannot produce a valid aromatic system",
            id="mdl-pyrrole",
        ),
        ("*", {}, UnderdeterminedError, "resolution underdetermined"),
        (
            "c1ccccc1",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ValenceModel.counts(
                        ValenceTable(
                            entries={
                                Element("C"): ValenceEntry(
                                    target_covalences=[4],
                                    aromatic_valences=[0],
                                )
                            }
                        )
                    ),
                    aromaticity=AromaticityModel(
                        scope=ElementScope.Any(),
                        rule=AromaticityRule.Hmo(stabilization_threshold=0.375),
                    ),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            RuntimeError,
            "hmo: missing parameters: no Van-Catledge parameters for C with "
            "0 pi-electrons",
        ),
    ],
)
def test_molecule_from_smiles_error(source, kwargs, error_type, message):
    with pytest.raises(error_type, match=f"^{re.escape(message)}$"):
        Molecule.from_smiles(source, **kwargs)


def test_molecule_from_smiles_keyword_error():
    with pytest.raises(TypeError):
        Molecule.from_smiles("C", SmilesIoConfig.opensmiles())


def test_molecule_from_smiles_ownership():
    io_config = SmilesIoConfig.opensmiles()
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.smiles(),
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )
    resolve_config = ResolveConfig.default()

    first = Molecule.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    second = Molecule.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    first.atoms[0].charge = 1

    assert second == Molecule.parse('{:atoms ["C#i=#c0#h4#n0#u0#s"]}')
    assert first != second
    assert io_config == SmilesIoConfig.opensmiles()
    assert chemistry_model == ChemistryModel(
        connectivity=default.connectivity,
        valence=ValenceModel.smiles(),
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )
    assert resolve_config == ResolveConfig.default()


def test_molecule_combine():
    left = Molecule.from_entries([AtomForm(Element("C"))])
    right = Molecule.from_entries(
        [AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(2))],
    )
    left_before = Molecule.from_entries([AtomForm(Element("C"))])
    right_before = Molecule.from_entries(
        [AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(2))],
    )

    combined, correspondence = left.combine(right)

    assert combined == Molecule.from_entries(
        [
            AtomForm(Element("C")),
            AtomForm(Element("O")),
            AtomForm(Element("N")),
        ],
        bonds=[(1, 2, BondForm(2))],
    )
    assert isinstance(correspondence, MoleculeCorrespondence)
    assert correspondence.atoms.matched_pairs == [(0, 1), (1, 2)]
    assert correspondence.bonds.matched_pairs == [(0, 0)]
    assert left == left_before
    assert right == right_before


def test_molecule_combine_from():
    recipient = Molecule.from_entries([AtomForm(Element("C"))])
    other = Molecule.from_entries(
        [AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(2))],
    )
    other_before = Molecule.from_entries(
        [AtomForm(Element("O")), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(2))],
    )

    correspondence = recipient.combine_from(other)

    assert recipient == Molecule.from_entries(
        [
            AtomForm(Element("C")),
            AtomForm(Element("O")),
            AtomForm(Element("N")),
        ],
        bonds=[(1, 2, BondForm(2))],
    )
    assert correspondence.atoms.matched_pairs == [(0, 1), (1, 2)]
    assert correspondence.bonds.matched_pairs == [(0, 0)]
    assert other == other_before


def test_molecule_combine_from_alias():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("O"))],
        bonds=[(0, 1, BondForm(1))],
    )

    correspondence = molecule.combine_from(molecule)

    assert molecule == Molecule.from_entries(
        [
            AtomForm(Element("C")),
            AtomForm(Element("O")),
            AtomForm(Element("C")),
            AtomForm(Element("O")),
        ],
        bonds=[(0, 1, BondForm(1)), (2, 3, BondForm(1))],
    )
    assert correspondence.atoms.matched_pairs == [(0, 2), (1, 3)]
    assert correspondence.bonds.matched_pairs == [(0, 1)]


def test_molecule_combine_all():
    molecules = [
        Molecule.from_entries([AtomForm(Element("C"))]),
        Molecule.from_entries(
            [AtomForm(Element("O")), AtomForm(Element("N"))],
            bonds=[(0, 1, BondForm(2))],
        ),
        Molecule.from_entries([AtomForm(Element("F"))]),
    ]
    snapshots = [
        Molecule.from_entries([AtomForm(Element("C"))]),
        Molecule.from_entries(
            [AtomForm(Element("O")), AtomForm(Element("N"))],
            bonds=[(0, 1, BondForm(2))],
        ),
        Molecule.from_entries([AtomForm(Element("F"))]),
    ]

    combined, correspondences = Molecule.combine_all(
        molecule for molecule in molecules
    )

    assert combined == Molecule.from_entries(
        [
            AtomForm(Element("C")),
            AtomForm(Element("O")),
            AtomForm(Element("N")),
            AtomForm(Element("F")),
        ],
        bonds=[(1, 2, BondForm(2))],
    )
    assert [correspondence.atoms.matched_pairs for correspondence in correspondences] == [
        [(0, 0)],
        [(0, 1), (1, 2)],
        [(0, 3)],
    ]
    assert [correspondence.bonds.matched_pairs for correspondence in correspondences] == [
        [],
        [(0, 0)],
        [],
    ]
    assert molecules == snapshots


def test_molecule_combine_all_empty():
    assert Molecule.combine_all([]) == (Molecule(), [])


def test_correspondence_constructor():
    correspondence = Correspondence([(2, 0), (0, 2)], 3, 3)

    assert correspondence.matched_pairs == [(0, 2), (2, 0)]
    assert correspondence.left_count == 3
    assert correspondence.right_count == 3
    assert correspondence.left_unmatched == [1]
    assert correspondence.right_unmatched == [1]
    assert repr(correspondence) == (
        "Correspondence(matched_pairs=[(0, 2), (2, 0)], "
        "left_count=3, right_count=3)"
    )


@pytest.mark.parametrize(
    "matched_pairs,left_count,right_count,message",
    [
        pytest.param(
            [(2, 0)],
            2,
            1,
            "left id 2 is out of range for 2 entries",
            id="left-out-of-range",
        ),
        pytest.param(
            [(0, 1)],
            1,
            1,
            "right id 1 is out of range for 1 entries",
            id="right-out-of-range",
        ),
        pytest.param(
            [(0, 0), (0, 1)],
            1,
            2,
            "left id 0 occurs more than once",
            id="duplicate-left",
        ),
        pytest.param(
            [(0, 0), (1, 0)],
            2,
            1,
            "right id 0 occurs more than once",
            id="duplicate-right",
        ),
    ],
)
def test_correspondence_constructor_error(
    matched_pairs, left_count, right_count, message
):
    with pytest.raises(ValueError, match=rf"^{re.escape(message)}$"):
        Correspondence(matched_pairs, left_count, right_count)


def test_correspondence_value():
    _, molecule_correspondence = Molecule.from_entries(
        [AtomForm(Element("C"))]
    ).combine(
        Molecule.from_entries(
            [AtomForm(Element("O")), AtomForm(Element("N"))],
            bonds=[(0, 1, BondForm(2))],
        )
    )
    correspondence = molecule_correspondence.atoms

    assert isinstance(correspondence, Correspondence)
    assert len(correspondence) == 2
    assert correspondence.matched_pairs == [(0, 1), (1, 2)]
    assert correspondence.left_count == 2
    assert correspondence.right_count == 3
    assert correspondence.left_unmatched == []
    assert correspondence.right_unmatched == [0]
    assert correspondence.right_of(0) == 1
    assert correspondence.right_of(2) is None
    assert correspondence.left_of(2) == 1
    assert correspondence.left_of(0) is None
    assert not correspondence.is_total()

    reverse = correspondence.reverse()
    composite = correspondence.compose(reverse)

    assert reverse.matched_pairs == [(1, 0), (2, 1)]
    assert reverse.left_count == 3
    assert reverse.right_count == 2
    assert composite.matched_pairs == [(0, 0), (1, 1)]
    assert composite.is_total()
    assert Correspondence.compose_all(
        item for item in [correspondence, reverse]
    ) == composite
    assert Correspondence.compose_all(iter(())) is None


@pytest.mark.parametrize("right_count,next_left_count", [(3, 1), (1, 3), (0, 1)])
def test_correspondence_compose_error(right_count, next_left_count):
    left = Correspondence([], 2, right_count)
    right = Correspondence([], next_left_count, 2)
    message = f"intermediate counts differ: {right_count} and {next_left_count}"
    with pytest.raises(ValueError, match=message):
        left.compose(right)
    with pytest.raises(ValueError, match=message):
        Correspondence.compose_all([left.reverse(), left, right])


def test_molecule_correspondence_compose_error():
    molecule = Molecule.from_entries([AtomForm(Element("C"))])
    _, correspondence = molecule.combine(molecule)
    with pytest.raises(ValueError, match="atom: intermediate counts differ: 2 and 1"):
        correspondence.compose(correspondence)
    with pytest.raises(ValueError, match="atom: intermediate counts differ: 2 and 1"):
        MoleculeCorrespondence.compose_all([correspondence, correspondence])


def test_molecule_correspondence_value():
    _, correspondence = Molecule.from_entries(
        [AtomForm(Element("C"))]
    ).combine(
        Molecule.from_entries(
            [AtomForm(Element("O")), AtomForm(Element("N"))],
            bonds=[(0, 1, BondForm(2))],
        )
    )

    assert not correspondence.is_total()

    reverse = correspondence.reverse()
    composite = correspondence.compose(reverse)

    assert reverse.atoms.matched_pairs == [(1, 0), (2, 1)]
    assert reverse.bonds.matched_pairs == [(0, 0)]
    assert not reverse.is_total()
    assert composite.atoms.matched_pairs == [(0, 0), (1, 1)]
    assert composite.bonds.matched_pairs == [(0, 0)]
    assert composite.is_total()
    assert MoleculeCorrespondence.compose_all(
        item for item in [correspondence, reverse]
    ) == composite
    assert MoleculeCorrespondence.compose_all(iter(())) is None


def test_molecule_split():
    molecule = Molecule.from_entries(
        [
            AtomForm(Element("C")),
            AtomForm(Element("O")),
            AtomForm(Element("N")),
        ],
        bonds=[(1, 2, BondForm(2))],
    )

    components = molecule.split()

    assert [component for component, _ in components] == [
        Molecule.from_entries([AtomForm(Element("C"))]),
        Molecule.from_entries(
            [AtomForm(Element("O")), AtomForm(Element("N"))],
            bonds=[(0, 1, BondForm(2))],
        ),
    ]
    assert [
        correspondence.atoms.matched_pairs for _, correspondence in components
    ] == [[(0, 0)], [(0, 1), (1, 2)]]
    assert [
        correspondence.bonds.matched_pairs for _, correspondence in components
    ] == [[], [(0, 0)]]


def test_molecule_split_empty():
    assert Molecule().split() == []


def test_molecule_bonds_error():
    with pytest.raises(IndexError):
        Molecule.from_entries([AtomForm(Element("C"))]).bonds[0]


def test_molecule_eq():
    assert Molecule() == Molecule()


@pytest.mark.parametrize(
    ("molecule", "expected"),
    [
        (Molecule(), "Molecule(atoms=0, bonds=0)"),
        (
            Molecule.from_entries(
                [
                    AtomForm(Element("C")),
                    AtomForm(Element("C")),
                    AtomForm(Element("O")),
                ],
                bonds=[(0, 1, BondForm(1))],
                aromatic_systems=[([0, 1, 2], AromaticSystemForm([1, 1, 1]))],
                noncovalent_bonds=[
                    (
                        [0, 2],
                        NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
                    )
                ],
            ),
            "Molecule(atoms=3, bonds=1, aromatic_systems=1, noncovalent_bonds=1)",
        ),
    ],
)
def test_molecule_repr(molecule, expected):
    assert repr(molecule) == expected
