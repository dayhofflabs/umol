import re

import pytest

from umol import (
    AromaticSystemAst,
    AromaticValenceAst,
    AromaticityConfig,
    AromaticityFailurePolicy,
    AromaticityModel,
    AromaticityResolveConfig,
    AtomAst,
    AtomTypeRegistry,
    BondAst,
    ChemistryModel,
    ConnectedComponentsAlgorithm,
    Constraint,
    ContradictionError,
    Correspondence,
    DativeBondAst,
    ElectronCountsAst,
    Element,
    ElementAst,
    ElementScope,
    Entity,
    MaximumIndependentSetAlgorithm,
    MetadataError,
    ModelConversionError,
    MoleculeAst,
    MoleculeConstraint,
    MoleculeCorrespondence,
    MoleculeDefaults,
    MoleculeMetadata,
    MulticenterBondAst,
    NoncovalentBondAst,
    NoncovalentBondKind,
    ParseError,
    ResolveConfig,
    RelevantCycleEnumerationAlgorithm,
    RingConfig,
    RingLimits,
    SimpleCycleEnumerationAlgorithm,
    SmilesIoConfig,
    SmilesSyntaxFlags,
    StereoAtomAst,
    StereoBondAst,
    StereoCoset,
    StereoFailurePolicy,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereoModel,
    StereoResolveConfig,
    TetrahedralConfiguration,
    UnderdeterminedError,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
    ValueAst,
)


def test_molecule_ast_new():
    assert len(MoleculeAst().atoms) == 0
    assert len(MoleculeAst().bonds) == 0


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (MoleculeDefaults(), MoleculeDefaults(), "MoleculeDefaults()"),
        (
            MoleculeDefaults.ground(),
            MoleculeDefaults.ground(),
            "MoleculeDefaults.ground()",
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
            MoleculeAst.from_entries([AtomAst.parse("C")]),
        ),
        (
            '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
            MoleculeDefaults.ground(),
            MoleculeAst.from_entries(
                [AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
            ),
        ),
    ],
)
def test_molecule_ast_parse(source, defaults, expected):
    assert MoleculeAst.parse(source, defaults=defaults) == expected


def test_molecule_ast_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        MoleculeAst.parse("not edn")


def test_molecule_ast_parse_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^MoleculeAst.parse\\(\\) takes 1 positional arguments but 2 were given$"
        ),
    ):
        MoleculeAst.parse('{:atoms ["C"]}', MoleculeDefaults.ground())


def test_molecule_ast_parse_with_metadata():
    source = (
        '{:atom-aliases [:x "C"] :atoms [[:carbon :x]] :bonds []}'
    )

    molecule, metadata = MoleculeAst.parse_with_metadata(source)

    assert molecule == MoleculeAst.from_entries([AtomAst(Element("C"))])
    assert metadata.keyword(Entity.Atom(0)) == "carbon"
    assert metadata.entity("carbon") == Entity.Atom(0)
    assert repr(metadata) == (
        'MoleculeMetadata(keywords=[(Entity.Atom(0), "carbon")], '
        "atom_alias_count=1)"
    )
    assert molecule.render_with_metadata(metadata) == source


def test_molecule_ast_parse_with_metadata_defaults():
    molecule, metadata = MoleculeAst.parse_with_metadata(
        '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.ground(),
    )

    assert molecule == MoleculeAst.from_entries(
        [
            AtomAst.parse(
                "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"
            )
        ]
    )
    assert metadata == MoleculeMetadata()


def test_molecule_ast_parse_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^MoleculeAst.parse_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        MoleculeAst.parse_with_metadata(
            '{:atoms ["C"]}',
            MoleculeDefaults.ground(),
        )


@pytest.mark.parametrize(
    ("molecule", "defaults", "expected"),
    [
        (
            MoleculeAst.parse('{:atoms ["C"]}'),
            MoleculeDefaults(),
            '{:atoms ["C"] :bonds []}',
        ),
        (
            MoleculeAst.parse(
                '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
                defaults=MoleculeDefaults.ground(),
            ),
            MoleculeDefaults.ground(),
            '{:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}',
        ),
    ],
)
def test_molecule_ast_render(molecule, defaults, expected):
    assert molecule.render(defaults=defaults) == expected


def test_molecule_ast_render_keyword_error():
    with pytest.raises(
        TypeError,
        match="^MoleculeAst.render\\(\\) takes 0 positional arguments but 1 was given$",
    ):
        MoleculeAst.parse('{:atoms ["C"]}').render(MoleculeDefaults())


def test_molecule_ast_render_with_metadata():
    source = (
        '{:atom-aliases [:x "C"] :atoms [[:carbon :x]] :bonds []}'
    )
    molecule, metadata = MoleculeAst.parse_with_metadata(source)

    assert molecule.render_with_metadata(metadata) == source
    assert molecule.render() == '{:atoms ["C"] :bonds []}'


def test_molecule_ast_render_with_metadata_error():
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(1), "outside")

    with pytest.raises(
        MetadataError,
        match="^metadata entity is out of range: atom 1$",
    ):
        MoleculeAst.parse('{:atoms ["C"]}').render_with_metadata(metadata)


def test_molecule_ast_render_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^MoleculeAst.render_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        MoleculeAst.parse('{:atoms ["C"]}').render_with_metadata(
            MoleculeMetadata(),
            MoleculeDefaults(),
        )


def test_molecule_ast_str():
    molecule = MoleculeAst.parse(
        '{:atoms ["C" "O"] :bonds [[0 1 "1"]]}'
    )

    assert str(molecule) == molecule.render()


def test_molecule_ast_from_entries():
    molecule = MoleculeAst.from_entries(
        [AtomAst(Element("C")) for _ in range(5)],
        bonds=[(0, 1, BondAst(2))],
        dative_bonds=[([2], 1, DativeBondAst(1))],
        aromatic_systems=[([0, 1, 2], AromaticSystemAst([1, 1, 1]))],
        multicenter_bonds=[([0, 1, 2], MulticenterBondAst([1, 1, 1]))],
        noncovalent_bonds=[
            ([0, 2], NoncovalentBondAst(NoncovalentBondKind.HydrogenBond))
        ],
        stereo_atoms=[
            (
                0,
                [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
                StereoAtomAst(TetrahedralConfiguration.Ccw),
            )
        ],
        stereo_bonds=[
            (
                0,
                [
                    StereoLigand(2, StereoLigandKind.Atom),
                    StereoLigand(3, StereoLigandKind.Atom),
                ],
                StereoBondAst.parse("Ct0"),
            )
        ],
        constraints=[
            Constraint.Molecule(MoleculeConstraint.Connected([0, 1, 2, 3, 4]))
        ],
    )

    assert len(molecule.atoms) == 5
    assert len(molecule.bonds) == 1
    assert len(molecule.dative_bonds) == 1
    assert len(molecule.aromatic_systems) == 1
    assert len(molecule.multicenter_bonds) == 1
    assert len(molecule.noncovalent_bonds) == 1
    assert len(molecule.stereo_atoms) == 1
    assert len(molecule.stereo_bonds) == 1
    assert list(molecule.constraints) == [
        Constraint.Molecule(MoleculeConstraint.Connected([0, 1, 2, 3, 4]))
    ]


def test_molecule_ast_from_entries_default():
    molecule = MoleculeAst.from_entries([AtomAst(Element("C"))])

    assert len(molecule.atoms) == 1
    assert len(molecule.bonds) == 0


def test_molecule_ast_from_entries_atom_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule entries reference unavailable atom 1$",
    ):
        MoleculeAst.from_entries(
            [AtomAst(Element("C"))],
            bonds=[(0, 1, BondAst(1))],
        )


def test_molecule_ast_from_entries_bond_site_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule entries reference unavailable bond 0$",
    ):
        MoleculeAst.from_entries(
            [AtomAst(Element("C"))],
            stereo_bonds=[
                (
                    0,
                    [StereoLigand(0, StereoLigandKind.Atom)],
                    StereoBondAst.parse("Ct0"),
                )
            ],
        )


def test_molecule_ast_from_entries_ligand_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule entries reference unavailable atom 1$",
    ):
        MoleculeAst.from_entries(
            [AtomAst(Element("C"))],
            stereo_atoms=[
                (
                    0,
                    [StereoLigand(1, StereoLigandKind.Atom)],
                    StereoAtomAst(TetrahedralConfiguration.Ccw),
                )
            ],
        )


def test_molecule_ast_from_entries_constraint_reference_error():
    with pytest.raises(
        ValueError,
        match="^molecule entries reference unavailable atom 1$",
    ):
        MoleculeAst.from_entries(
            [AtomAst(Element("C"))],
            constraints=[Constraint.Molecule(MoleculeConstraint.Connected([1]))],
        )


def test_molecule_ast_from_smiles():
    assert MoleculeAst.from_smiles("C") == MoleculeAst.parse(
        '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.ground(),
    )


@pytest.mark.parametrize(
    "io_config",
    [
        SmilesIoConfig.with_syntax_flags(
            syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS
        ),
        SmilesIoConfig.lenient(),
    ],
)
def test_molecule_ast_from_smiles_io_config(io_config):
    molecule = MoleculeAst.from_smiles("[se]1cccc1", io_config=io_config)

    assert [atom.element for atom in molecule.atoms] == [
        ElementAst.Lit(Element("Se")),
        *[ElementAst.Lit(Element("C")) for _ in range(4)],
    ]
    assert [(bond.atom_ids, bond.order) for bond in molecule.bonds] == [
        ((0, 4), ValueAst.Lit(1)),
        ((0, 1), ValueAst.Lit(1)),
        ((1, 2), ValueAst.Lit(1)),
        ((2, 3), ValueAst.Lit(1)),
        ((3, 4), ValueAst.Lit(1)),
    ]
    assert [
        (system.atom_ids, system.electrons, system.charge)
        for system in molecule.aromatic_systems
    ] == [
        (
            (0, 1, 2, 3, 4),
            ElectronCountsAst.Lit([2, 1, 1, 1, 1]),
            ValueAst.Lit(0),
        )
    ]


def test_molecule_ast_from_smiles_io_config_error():
    with pytest.raises(ParseError, match="^Invalid token at position 2$"):
        MoleculeAst.from_smiles("C->N", io_config=SmilesIoConfig.opensmiles())
    with pytest.raises(
        ContradictionError,
        match="^no atom-typing match for AtomId\\(0\\) "
        "\\(element C, charge Some\\(0\\)\\)$",
    ):
        MoleculeAst.from_smiles("C->N", io_config=SmilesIoConfig.lenient())


@pytest.mark.parametrize(
    ("valence_model", "expected"),
    [
        (
            ValenceModel.AtomTyping(
                registry=AtomTypeRegistry.from_atoms(
                    [AtomAst.parse("C#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
                )
            ),
            AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"),
        ),
        (
            ValenceModel.Counts(table=ValenceTable.default()),
            AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#a!"),
        ),
    ],
)
def test_molecule_ast_from_smiles_chemistry_model_valence(valence_model, expected):
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=valence_model,
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )

    assert MoleculeAst.from_smiles(
        "C", chemistry_model=chemistry_model
    ) == MoleculeAst.from_entries([expected])


def test_molecule_ast_from_smiles_chemistry_model_aromaticity():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=default.valence,
        aromaticity=AromaticityModel.Hmo(
            scope=ElementScope.Any(), stabilization_threshold=0.375
        ),
        stereo=default.stereo,
    )

    molecule = MoleculeAst.from_smiles(
        "c1ccccc1",
        chemistry_model=chemistry_model,
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep
            ),
            stereo=StereoResolveConfig(),
        ),
    )

    assert [atom.implicit_hydrogens for atom in molecule.atoms] == [ValueAst.Lit(1)] * 6
    assert [atom.constraints.aromatic_valence for atom in molecule.atoms] == [
        AromaticValenceAst.Aromatic(ValueAst.Lit(1))
    ] * 6
    assert list(molecule.aromatic_systems) == []


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        pytest.param(
            "o1cccc1",
            '{:atoms ["O#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="furan",
        ),
        pytest.param(
            "s1cccc1",
            '{:atoms ["S#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="thiophene",
        ),
        pytest.param(
            "[nH]1cccc1",
            '{:atoms ["N#i=#c0#h#n0#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}',
            id="pyrrole",
        ),
    ],
)
def test_molecule_ast_from_smiles_aromaticity_policy(source, expected):
    default = ChemistryModel.default()

    assert MoleculeAst.from_smiles(
        source,
        chemistry_model=ChemistryModel(
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
    ) == MoleculeAst.parse(expected)


def test_molecule_ast_from_smiles_chemistry_model_stereo():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=default.valence,
        aromaticity=default.aromaticity,
        stereo=StereoModel(
            kind_models={},
            para_stereo=False,
        ),
    )

    molecule = MoleculeAst.from_smiles(
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
        ValueAst.Lit(3),
        ValueAst.Lit(1),
        ValueAst.Lit(2),
        ValueAst.Lit(1),
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
                [ValueAst.Lit(1), ValueAst.Lit(0), ValueAst.Lit(0)],
                [
                    AromaticValenceAst.Aromatic(ValueAst.Lit(0)),
                    AromaticValenceAst.Aromatic(ValueAst.Lit(1)),
                    AromaticValenceAst.Aromatic(ValueAst.Lit(1)),
                ],
                [None] * 3,
                [((0, 1, 2), ValueAst.Lit(0))],
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
                [ValueAst.Lit(0)] * 6,
                [None] * 6,
                [None] * 6,
                [((0, 1, 2, 3, 4, 5), ValueAst.Lit(0))],
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
                [ValueAst.Lit(0)] * 4,
                [AromaticValenceAst.NotAromatic()] * 4,
                [None] * 4,
                [],
                [(1, StereoKind.Tetrahedral, StereoCoset.Lit(0))],
            ),
        ),
    ],
)
def test_molecule_ast_from_smiles_resolve_config(source, resolve_config, expected):
    molecule = MoleculeAst.from_smiles(source, resolve_config=resolve_config)

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
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.Clar(
                        scope=ElementScope.Any(), ring_limits=RingLimits()
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
        pytest.param(
            "c1cccn1",
            {},
            ContradictionError,
            "aromaticity inconsistency: aromatic valence at atom AtomId(0) "
            "cannot produce a valid aromatic system",
            id="bare-aromatic-nitrogen",
        ),
        ("*", {}, UnderdeterminedError, "resolution underdetermined"),
        (
            "c1ccccc1",
            {
                "chemistry_model": ChemistryModel(
                    valence=ValenceModel.Counts(
                        table=ValenceTable(
                            entries={
                                Element("C"): ValenceEntry(
                                    target_covalences=[4],
                                    aromatic_valences=[0],
                                )
                            }
                        )
                    ),
                    aromaticity=AromaticityModel.Hmo(
                        scope=ElementScope.Any(),
                        stabilization_threshold=0.375,
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
def test_molecule_ast_from_smiles_error(source, kwargs, error_type, message):
    with pytest.raises(error_type, match=f"^{re.escape(message)}$"):
        MoleculeAst.from_smiles(source, **kwargs)


def test_molecule_ast_from_smiles_keyword_error():
    with pytest.raises(TypeError):
        MoleculeAst.from_smiles("C", SmilesIoConfig.opensmiles())


def test_molecule_ast_from_smiles_ownership():
    io_config = SmilesIoConfig.opensmiles()
    chemistry_model = ChemistryModel.default()
    resolve_config = ResolveConfig.default()

    first = MoleculeAst.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    second = MoleculeAst.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    first.atoms[0].charge = 1

    assert second == MoleculeAst.parse(
        '{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.ground(),
    )
    assert first != second
    assert io_config == SmilesIoConfig.opensmiles()
    assert chemistry_model == ChemistryModel.default()
    assert resolve_config == ResolveConfig.default()


def test_molecule_ast_combine():
    left = MoleculeAst.from_entries([AtomAst(Element("C"))])
    right = MoleculeAst.from_entries(
        [AtomAst(Element("O")), AtomAst(Element("N"))],
        bonds=[(0, 1, BondAst(2))],
    )
    left_before = MoleculeAst.from_entries([AtomAst(Element("C"))])
    right_before = MoleculeAst.from_entries(
        [AtomAst(Element("O")), AtomAst(Element("N"))],
        bonds=[(0, 1, BondAst(2))],
    )

    combined, correspondence = left.combine(right)

    assert combined == MoleculeAst.from_entries(
        [
            AtomAst(Element("C")),
            AtomAst(Element("O")),
            AtomAst(Element("N")),
        ],
        bonds=[(1, 2, BondAst(2))],
    )
    assert isinstance(correspondence, MoleculeCorrespondence)
    assert correspondence.atoms.matched_pairs == [(0, 1), (1, 2)]
    assert correspondence.bonds.matched_pairs == [(0, 0)]
    assert left == left_before
    assert right == right_before


def test_molecule_ast_combine_from():
    recipient = MoleculeAst.from_entries([AtomAst(Element("C"))])
    other = MoleculeAst.from_entries(
        [AtomAst(Element("O")), AtomAst(Element("N"))],
        bonds=[(0, 1, BondAst(2))],
    )
    other_before = MoleculeAst.from_entries(
        [AtomAst(Element("O")), AtomAst(Element("N"))],
        bonds=[(0, 1, BondAst(2))],
    )

    correspondence = recipient.combine_from(other)

    assert recipient == MoleculeAst.from_entries(
        [
            AtomAst(Element("C")),
            AtomAst(Element("O")),
            AtomAst(Element("N")),
        ],
        bonds=[(1, 2, BondAst(2))],
    )
    assert correspondence.atoms.matched_pairs == [(0, 1), (1, 2)]
    assert correspondence.bonds.matched_pairs == [(0, 0)]
    assert other == other_before


def test_molecule_ast_combine_from_alias():
    molecule = MoleculeAst.from_entries(
        [AtomAst(Element("C")), AtomAst(Element("O"))],
        bonds=[(0, 1, BondAst(1))],
    )

    correspondence = molecule.combine_from(molecule)

    assert molecule == MoleculeAst.from_entries(
        [
            AtomAst(Element("C")),
            AtomAst(Element("O")),
            AtomAst(Element("C")),
            AtomAst(Element("O")),
        ],
        bonds=[(0, 1, BondAst(1)), (2, 3, BondAst(1))],
    )
    assert correspondence.atoms.matched_pairs == [(0, 2), (1, 3)]
    assert correspondence.bonds.matched_pairs == [(0, 1)]


def test_molecule_ast_combine_all():
    molecules = [
        MoleculeAst.from_entries([AtomAst(Element("C"))]),
        MoleculeAst.from_entries(
            [AtomAst(Element("O")), AtomAst(Element("N"))],
            bonds=[(0, 1, BondAst(2))],
        ),
        MoleculeAst.from_entries([AtomAst(Element("F"))]),
    ]
    snapshots = [
        MoleculeAst.from_entries([AtomAst(Element("C"))]),
        MoleculeAst.from_entries(
            [AtomAst(Element("O")), AtomAst(Element("N"))],
            bonds=[(0, 1, BondAst(2))],
        ),
        MoleculeAst.from_entries([AtomAst(Element("F"))]),
    ]

    combined, correspondences = MoleculeAst.combine_all(
        molecule for molecule in molecules
    )

    assert combined == MoleculeAst.from_entries(
        [
            AtomAst(Element("C")),
            AtomAst(Element("O")),
            AtomAst(Element("N")),
            AtomAst(Element("F")),
        ],
        bonds=[(1, 2, BondAst(2))],
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


def test_molecule_ast_combine_all_empty():
    assert MoleculeAst.combine_all([]) == (MoleculeAst(), [])


def test_correspondence_value():
    _, molecule_correspondence = MoleculeAst.from_entries(
        [AtomAst(Element("C"))]
    ).combine(
        MoleculeAst.from_entries(
            [AtomAst(Element("O")), AtomAst(Element("N"))],
            bonds=[(0, 1, BondAst(2))],
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


def test_molecule_correspondence_value():
    _, correspondence = MoleculeAst.from_entries(
        [AtomAst(Element("C"))]
    ).combine(
        MoleculeAst.from_entries(
            [AtomAst(Element("O")), AtomAst(Element("N"))],
            bonds=[(0, 1, BondAst(2))],
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


def test_molecule_ast_split():
    molecule = MoleculeAst.from_entries(
        [
            AtomAst(Element("C")),
            AtomAst(Element("O")),
            AtomAst(Element("N")),
        ],
        bonds=[(1, 2, BondAst(2))],
    )

    components = molecule.split()

    assert [component for component, _ in components] == [
        MoleculeAst.from_entries([AtomAst(Element("C"))]),
        MoleculeAst.from_entries(
            [AtomAst(Element("O")), AtomAst(Element("N"))],
            bonds=[(0, 1, BondAst(2))],
        ),
    ]
    assert [
        correspondence.atoms.matched_pairs for _, correspondence in components
    ] == [[(0, 0)], [(0, 1), (1, 2)]]
    assert [
        correspondence.bonds.matched_pairs for _, correspondence in components
    ] == [[], [(0, 0)]]


def test_molecule_ast_split_empty():
    assert MoleculeAst().split() == []


def test_molecule_ast_bonds_error():
    with pytest.raises(IndexError):
        MoleculeAst.from_entries([AtomAst(Element("C"))]).bonds[0]


def test_molecule_ast_eq():
    assert MoleculeAst() == MoleculeAst()


@pytest.mark.parametrize(
    ("molecule", "expected"),
    [
        (MoleculeAst(), "MoleculeAst(atoms=0, bonds=0)"),
        (
            MoleculeAst.from_entries(
                [
                    AtomAst(Element("C")),
                    AtomAst(Element("C")),
                    AtomAst(Element("O")),
                ],
                bonds=[(0, 1, BondAst(1))],
                aromatic_systems=[([0, 1, 2], AromaticSystemAst([1, 1, 1]))],
                noncovalent_bonds=[
                    (
                        [0, 2],
                        NoncovalentBondAst(NoncovalentBondKind.HydrogenBond),
                    )
                ],
            ),
            "MoleculeAst(atoms=3, bonds=1, aromatic_systems=1, noncovalent_bonds=1)",
        ),
    ],
)
def test_molecule_ast_repr(molecule, expected):
    assert repr(molecule) == expected
