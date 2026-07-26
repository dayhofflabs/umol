import pytest

from umol import (
    AutomorphismAlgorithm,
    AtomAst,
    BitFp,
    CountedHashedFeatureSet,
    EcfpHashScheme,
    HashedFingerprintConfig,
    HashedFeatureSet,
    MoleculeAst,
    PatternFingerprintConfig,
    ReactionAst,
    ReactionCombinedFingerprint,
    ReactionCombinedFingerprintConfig,
    ReactionDefaults,
    ReactionSide,
    RefinementRounds,
    RelevantCycleEnumerationAlgorithm,
    RingConfig,
    RoleTaggedHashedFeatureSet,
    SignedHashedFeatureSet,
    SimpleCycleEnumerationAlgorithm,
    StructuralFeatureSet,
    StructuralFingerprintConfig,
    SubgraphEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
    UnderdeterminedError,
    WlHashScheme,
)


@pytest.fixture
def ethanol():
    return MoleculeAst.from_smiles("CCO")


@pytest.fixture
def undetermined_molecule():
    return MoleculeAst.from_parts([AtomAst.parse("C")])


@pytest.fixture
def ethanol_deoxygenation():
    return ReactionAst.parse(
        "{:deltas [{:atom {:remove 2}} {:bond {:remove 1}}] "
        ':lhs {:atoms ["C#h3#v#d0#t0#a!#m!" '
        '"C#h2#v2#d0#t0#a!#m!" "O#h#n2#v#d0#t0#a!#m!"] '
        ':bonds [[0 1 "1"] [1 2 "1"]]}}',
        defaults=ReactionDefaults.ground(),
    )


@pytest.mark.parametrize("rounds", [0, 3])
def test_refinement_rounds_fixed(rounds):
    value = RefinementRounds.Fixed(rounds=rounds)

    assert value.rounds == rounds
    assert value == RefinementRounds.Fixed(rounds=rounds)
    assert repr(value) == f"RefinementRounds.Fixed(rounds={rounds})"


def test_refinement_rounds_to_fixpoint():
    value = RefinementRounds.ToFixpoint()

    assert value == RefinementRounds.ToFixpoint()
    assert value != RefinementRounds.Fixed(rounds=0)
    assert repr(value) == "RefinementRounds.ToFixpoint()"
    with pytest.raises(AttributeError):
        value.rounds


def test_wl_hash_scheme():
    value = WlHashScheme.Xxh3SortedWidth64V1()

    assert value == WlHashScheme.Xxh3SortedWidth64V1()
    assert repr(value) == "WlHashScheme.Xxh3SortedWidth64V1()"


def test_ecfp_hash_scheme():
    value = EcfpHashScheme.Xxh3Width64V1()

    assert value == EcfpHashScheme.Xxh3Width64V1()
    assert repr(value) == "EcfpHashScheme.Xxh3Width64V1()"


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (
            HashedFingerprintConfig.Morgan(),
            HashedFingerprintConfig.Morgan(
                radius=2,
                ring_config=RingConfig(),
            ),
            "HashedFingerprintConfig.Morgan(radius=2, "
            "ring_config=RingConfig("
            "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()))",
        ),
        (
            HashedFingerprintConfig.Ecfp(),
            HashedFingerprintConfig.Ecfp(
                radius=2,
                hashing_scheme=EcfpHashScheme.Xxh3Width64V1(),
                ring_config=RingConfig(),
            ),
            "HashedFingerprintConfig.Ecfp(radius=2, "
            "hashing_scheme=EcfpHashScheme.Xxh3Width64V1(), "
            "ring_config=RingConfig("
            "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()))",
        ),
        (
            HashedFingerprintConfig.Wl(rounds=RefinementRounds.ToFixpoint()),
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.ToFixpoint(),
                hashing_scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            "HashedFingerprintConfig.Wl(rounds=RefinementRounds.ToFixpoint(), "
            "hashing_scheme=WlHashScheme.Xxh3SortedWidth64V1())",
        ),
    ],
)
def test_hashed_fingerprint_config_defaults(value, expected, expected_repr):
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (
            HashedFingerprintConfig.Morgan(
                radius=3,
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                ),
            ),
            HashedFingerprintConfig.Morgan(
                radius=3,
                ring_config=RingConfig(),
            ),
            "HashedFingerprintConfig.Morgan(radius=3, "
            "ring_config=RingConfig("
            "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()))",
        ),
        (
            HashedFingerprintConfig.Ecfp(
                radius=3,
                hashing_scheme=EcfpHashScheme.Xxh3Width64V1(),
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                ),
            ),
            HashedFingerprintConfig.Ecfp(
                radius=3,
                hashing_scheme=EcfpHashScheme.Xxh3Width64V1(),
                ring_config=RingConfig(),
            ),
            "HashedFingerprintConfig.Ecfp(radius=3, "
            "hashing_scheme=EcfpHashScheme.Xxh3Width64V1(), "
            "ring_config=RingConfig("
            "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
            "relevant_cycle_algorithm="
            "RelevantCycleEnumerationAlgorithm.Vismara()))",
        ),
        (
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
                hashing_scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
                hashing_scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            "HashedFingerprintConfig.Wl(rounds=RefinementRounds.Fixed(rounds=3), "
            "hashing_scheme=WlHashScheme.Xxh3SortedWidth64V1())",
        ),
    ],
)
def test_hashed_fingerprint_config(value, expected, expected_repr):
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        (
            HashedFingerprintConfig.Morgan(),
            RingConfig(),
        ),
        (
            HashedFingerprintConfig.Ecfp(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            ),
            RingConfig(),
        ),
    ],
)
def test_hashed_fingerprint_config_ring_config(config, expected):
    ring_config = config.ring_config

    assert ring_config == expected
    assert ring_config is not config.ring_config


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (
            PatternFingerprintConfig(),
            PatternFingerprintConfig(
                width=2048,
                match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(),
                subgraph_isomorphism_algorithm=(
                    SubgraphIsomorphismAlgorithm.Vf2Rdkit()
                ),
            ),
            "PatternFingerprintConfig(width=2048, "
            "match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(), "
            "subgraph_isomorphism_algorithm="
            "SubgraphIsomorphismAlgorithm.Vf2Rdkit())",
        ),
        (
            PatternFingerprintConfig(
                width=512,
                match_algorithm=SubstructureMatchAlgorithm.Incidence(),
                subgraph_isomorphism_algorithm=(
                    SubgraphIsomorphismAlgorithm.Ullmann()
                ),
            ),
            PatternFingerprintConfig(
                width=512,
                match_algorithm=SubstructureMatchAlgorithm.Incidence(),
                subgraph_isomorphism_algorithm=(
                    SubgraphIsomorphismAlgorithm.Ullmann()
                ),
            ),
            "PatternFingerprintConfig(width=512, "
            "match_algorithm=SubstructureMatchAlgorithm.Incidence(), "
            "subgraph_isomorphism_algorithm="
            "SubgraphIsomorphismAlgorithm.Ullmann())",
        ),
    ],
)
def test_pattern_fingerprint_config(value, expected, expected_repr):
    assert value.width == expected.width
    assert value.match_algorithm == expected.match_algorithm
    assert (
        value.subgraph_isomorphism_algorithm
        == expected.subgraph_isomorphism_algorithm
    )
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize("width", [0, -1])
def test_pattern_fingerprint_config_error(width):
    with pytest.raises(ValueError, match="width must be positive"):
        PatternFingerprintConfig(width=width)


@pytest.mark.parametrize("max_bonds", [0, 3])
def test_structural_fingerprint_config(max_bonds):
    value = StructuralFingerprintConfig(max_bonds=max_bonds)
    expected = StructuralFingerprintConfig(
        max_bonds=max_bonds,
        subgraph_enumeration_algorithm=SubgraphEnumerationAlgorithm.Esu(),
        automorphism_algorithm=AutomorphismAlgorithm.Nauty(),
    )

    assert value.max_bonds == max_bonds
    assert (
        value.subgraph_enumeration_algorithm
        == expected.subgraph_enumeration_algorithm
    )
    assert value.automorphism_algorithm == expected.automorphism_algorithm
    assert value == expected
    assert repr(value) == (
        f"StructuralFingerprintConfig(max_bonds={max_bonds}, "
        "subgraph_enumeration_algorithm=SubgraphEnumerationAlgorithm.Esu(), "
        "automorphism_algorithm=AutomorphismAlgorithm.Nauty())"
    )


@pytest.mark.parametrize(
    ("variant", "variant_name"),
    [
        (
            ReactionCombinedFingerprintConfig.Difference,
            "Difference",
        ),
        (
            ReactionCombinedFingerprintConfig.DisjointUnion,
            "DisjointUnion",
        ),
    ],
)
@pytest.mark.parametrize(
    "molecule",
    [
        HashedFingerprintConfig.Morgan(),
        HashedFingerprintConfig.Ecfp(),
        HashedFingerprintConfig.Wl(
            rounds=RefinementRounds.Fixed(rounds=3),
        ),
    ],
)
def test_reaction_combined_fingerprint_config(variant, variant_name, molecule):
    value = variant(molecule=molecule)

    assert value.molecule == molecule
    assert value == variant(molecule=molecule)
    assert repr(value) == (
        f"ReactionCombinedFingerprintConfig.{variant_name}(molecule={molecule!r})"
    )


@pytest.mark.parametrize(
    "constructor",
    [
        RefinementRounds.Fixed,
        HashedFingerprintConfig.Wl,
        StructuralFingerprintConfig,
        ReactionCombinedFingerprintConfig.Difference,
        ReactionCombinedFingerprintConfig.DisjointUnion,
    ],
)
def test_fingerprint_config_required_error(constructor):
    with pytest.raises(TypeError):
        constructor()


@pytest.mark.parametrize(
    ("constructor", "argument"),
    [
        (RefinementRounds.Fixed, 3),
        (HashedFingerprintConfig.Morgan, 3),
        (HashedFingerprintConfig.Ecfp, 3),
        (
            HashedFingerprintConfig.Wl,
            RefinementRounds.Fixed(rounds=3),
        ),
        (PatternFingerprintConfig, 512),
        (StructuralFingerprintConfig, 3),
        (
            ReactionCombinedFingerprintConfig.Difference,
            HashedFingerprintConfig.Morgan(),
        ),
        (
            ReactionCombinedFingerprintConfig.DisjointUnion,
            HashedFingerprintConfig.Morgan(),
        ),
    ],
)
def test_fingerprint_config_keyword_error(constructor, argument):
    with pytest.raises(TypeError):
        constructor(argument)


@pytest.mark.parametrize(
    "result_type",
    [
        HashedFeatureSet,
        CountedHashedFeatureSet,
        BitFp,
        StructuralFeatureSet,
        SignedHashedFeatureSet,
        RoleTaggedHashedFeatureSet,
        ReactionCombinedFingerprint,
    ],
)
def test_fingerprint_result_constructor_error(result_type):
    with pytest.raises(TypeError):
        result_type()


@pytest.mark.parametrize(
    ("config", "expected_ids"),
    [
        (
            HashedFingerprintConfig.Morgan(),
            [
                864662311,
                1535166686,
                2245384272,
                2246728737,
                3542456614,
                4018048386,
            ],
        ),
        pytest.param(
            HashedFingerprintConfig.Morgan(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            ),
            [
                864662311,
                1535166686,
                2245384272,
                2246728737,
                3542456614,
                4018048386,
            ],
            id="morgan-explicit-ring-config",
        ),
        (
            HashedFingerprintConfig.Ecfp(),
            [
                63839236075656913,
                1189585227353469813,
                3822471596818936039,
                13652293261850732425,
                15001976065402722634,
                16149328945726899460,
            ],
        ),
        pytest.param(
            HashedFingerprintConfig.Ecfp(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            ),
            [
                63839236075656913,
                1189585227353469813,
                3822471596818936039,
                13652293261850732425,
                15001976065402722634,
                16149328945726899460,
            ],
            id="ecfp-explicit-ring-config",
        ),
        (
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
            ),
            [
                2520347590860685079,
                3352603313223549703,
                4152249898001161146,
                5715207763479934940,
                5807737097854608645,
                7542810387455301591,
                11457795998246593156,
                11986000156817227245,
                12895020514073294021,
                13932567567828606490,
                17305796300852423160,
                17417400371411086222,
            ],
        ),
    ],
)
def test_molecule_ast_hashed_fingerprint(ethanol, config, expected_ids):
    fingerprint = ethanol.hashed_fingerprint(config=config)
    ids = fingerprint.ids

    assert ids == expected_ids
    assert list(fingerprint) == expected_ids
    assert fingerprint.id_width == 64
    assert all(type(identifier) is int for identifier in fingerprint)
    ids.append(9)
    assert fingerprint.ids == expected_ids


@pytest.mark.parametrize(
    ("config", "expected_entries"),
    [
        (
            HashedFingerprintConfig.Morgan(),
            [(2246728737, 2), (3545175291, 1)],
        ),
        pytest.param(
            HashedFingerprintConfig.Morgan(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            ),
            [(2246728737, 2), (3545175291, 1)],
            id="morgan-explicit-ring-config",
        ),
        (
            HashedFingerprintConfig.Ecfp(),
            [(5513743581508886362, 1), (16149328945726899460, 2)],
        ),
        pytest.param(
            HashedFingerprintConfig.Ecfp(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            ),
            [(5513743581508886362, 1), (16149328945726899460, 2)],
            id="ecfp-explicit-ring-config",
        ),
        (
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
            ),
            [
                (2659163409134283895, 2),
                (7542810387455301591, 2),
                (9541344068636876323, 2),
                (12512207080905326651, 2),
            ],
        ),
    ],
)
def test_molecule_ast_counted_hashed_fingerprint(config, expected_entries):
    fingerprint = MoleculeAst.from_smiles("CC").counted_hashed_fingerprint(
        config=config
    )
    entries = fingerprint.entries

    assert entries == expected_entries
    assert list(fingerprint) == expected_entries
    assert fingerprint.id_width == 64
    assert fingerprint.count(expected_entries[0][0]) == expected_entries[0][1]
    assert fingerprint.count(0) == 0
    assert all(
        type(identifier) is int and type(count) is int
        for identifier, count in fingerprint
    )
    entries.append((9, 3))
    assert fingerprint.entries == expected_entries


@pytest.mark.parametrize(
    ("config", "expected_width", "expected_bits"),
    [
        (
            None,
            2048,
            [
                54,
                173,
                217,
                429,
                622,
                759,
                778,
                874,
                946,
                967,
                1022,
                1033,
                1061,
                1236,
                1289,
                1295,
            ],
        ),
        (
            PatternFingerprintConfig(),
            2048,
            [
                54,
                173,
                217,
                429,
                622,
                759,
                778,
                874,
                946,
                967,
                1022,
                1033,
                1061,
                1236,
                1289,
                1295,
            ],
        ),
        (
            PatternFingerprintConfig(
                width=64,
                match_algorithm=SubstructureMatchAlgorithm.Incidence(),
                subgraph_isomorphism_algorithm=(
                    SubgraphIsomorphismAlgorithm.Ullmann()
                ),
            ),
            64,
            [7, 9, 10, 15, 20, 25, 37, 42, 45, 46, 50, 54, 55, 62],
        ),
    ],
)
def test_molecule_ast_pattern_fingerprint(
    ethanol, config, expected_width, expected_bits
):
    if config is None:
        fingerprint = ethanol.pattern_fingerprint()
    else:
        fingerprint = ethanol.pattern_fingerprint(config=config)

    assert fingerprint.width == expected_width
    assert fingerprint.count_ones() == len(expected_bits)
    assert [bit for bit in range(expected_width) if fingerprint[bit]] == expected_bits
    assert type(fingerprint[0]) is bool


@pytest.mark.parametrize(
    ("config", "expected_keys"),
    [
        (
            StructuralFingerprintConfig(max_bonds=0),
            [
                bytes.fromhex("01 00 00 00 05 00 00 00 00 06 00 00 00 00 00 00 00"),
                bytes.fromhex("01 00 00 00 05 00 00 00 00 08 00 00 00 00 00 00 00"),
            ],
        ),
        (
            StructuralFingerprintConfig(
                max_bonds=2,
                subgraph_enumeration_algorithm=SubgraphEnumerationAlgorithm.Esu(),
                automorphism_algorithm=AutomorphismAlgorithm.Nauty(),
            ),
            [
                bytes.fromhex("01 00 00 00 05 00 00 00 00 06 00 00 00 00 00 00 00"),
                bytes.fromhex("01 00 00 00 05 00 00 00 00 08 00 00 00 00 00 00 00"),
                bytes.fromhex(
                    "03 00 00 00 05 00 00 00 00 06 00 00 00 05 00 00 00 00 "
                    "06 00 00 00 03 00 00 00 01 01 00 02 00 00 00 00 00 00 "
                    "00 02 00 00 00 01 00 00 00 02 00 00 00"
                ),
                bytes.fromhex(
                    "03 00 00 00 05 00 00 00 00 06 00 00 00 05 00 00 00 00 "
                    "08 00 00 00 03 00 00 00 01 01 00 02 00 00 00 00 00 00 "
                    "00 02 00 00 00 01 00 00 00 02 00 00 00"
                ),
                bytes.fromhex(
                    "05 00 00 00 05 00 00 00 00 06 00 00 00 05 00 00 00 00 "
                    "06 00 00 00 05 00 00 00 00 08 00 00 00 03 00 00 00 01 "
                    "01 00 03 00 00 00 01 01 00 04 00 00 00 00 00 00 00 03 "
                    "00 00 00 01 00 00 00 03 00 00 00 01 00 00 00 04 00 00 "
                    "00 02 00 00 00 04 00 00 00"
                ),
            ],
        ),
    ],
)
def test_molecule_ast_structural_fingerprint(ethanol, config, expected_keys):
    fingerprint = ethanol.structural_fingerprint(config=config)
    keys = fingerprint.keys

    assert keys == expected_keys
    assert list(fingerprint) == expected_keys
    assert all(type(key) is bytes for key in fingerprint)
    keys.append(b"detached")
    assert fingerprint.keys == expected_keys


def test_hashed_feature_set_operations():
    config = HashedFingerprintConfig.Morgan()
    ethane = MoleculeAst.from_smiles("CC").hashed_fingerprint(config=config)
    propane = MoleculeAst.from_smiles("CCC").hashed_fingerprint(config=config)
    folded = ethane.fold(64)

    assert ethane.tanimoto(propane) == pytest.approx(0.2)
    assert ethane.dice(propane) == pytest.approx(1 / 3)
    assert ethane.is_subset(propane) is False
    assert propane.is_subset(ethane) is False
    assert folded.width == 64
    assert [bit for bit in range(64) if folded[bit]] == [33, 59]


def test_hashed_feature_set_fold_error():
    fingerprint = MoleculeAst.from_smiles("CC").hashed_fingerprint(
        config=HashedFingerprintConfig.Morgan()
    )

    with pytest.raises(ValueError, match="width must be positive"):
        fingerprint.fold(0)


def test_bit_fp_operations():
    config = PatternFingerprintConfig(width=64)
    ethane = MoleculeAst.from_smiles("CC").pattern_fingerprint(config=config)
    propane = MoleculeAst.from_smiles("CCC").pattern_fingerprint(config=config)

    assert [bit for bit in range(64) if ethane[bit]] == [10, 15, 20, 37, 45, 62]
    assert [bit for bit in range(64) if propane[bit]] == [
        7,
        9,
        10,
        15,
        20,
        25,
        37,
        45,
        46,
        62,
    ]
    assert ethane.tanimoto(propane) == pytest.approx(0.6)
    assert ethane.dice(propane) == pytest.approx(0.75)
    assert ethane.is_subset(propane) is True
    assert propane.is_subset(ethane) is False


@pytest.mark.parametrize("index", [64, -65])
def test_bit_fp_getitem_error(index):
    fingerprint = MoleculeAst.from_smiles("CC").pattern_fingerprint(
        config=PatternFingerprintConfig(width=64)
    )

    with pytest.raises(IndexError, match="bit index out of range"):
        fingerprint[index]


@pytest.mark.parametrize("operation", ["tanimoto", "dice", "is_subset"])
def test_bit_fp_operations_error(operation):
    molecule = MoleculeAst.from_smiles("CC")
    narrow = molecule.pattern_fingerprint(config=PatternFingerprintConfig(width=64))
    wide = molecule.pattern_fingerprint(config=PatternFingerprintConfig(width=128))

    with pytest.raises(ValueError, match="fingerprint width mismatch: 64 != 128"):
        getattr(narrow, operation)(wide)


def test_structural_feature_set_is_subset():
    config = StructuralFingerprintConfig(max_bonds=2)
    ethane = MoleculeAst.from_smiles("CC").structural_fingerprint(config=config)
    propane = MoleculeAst.from_smiles("CCC").structural_fingerprint(config=config)

    assert ethane.is_subset(propane) is True
    assert propane.is_subset(ethane) is False


@pytest.mark.parametrize(
    ("method_name", "config"),
    [
        ("hashed_fingerprint", HashedFingerprintConfig.Morgan()),
        ("counted_hashed_fingerprint", HashedFingerprintConfig.Morgan()),
        ("pattern_fingerprint", PatternFingerprintConfig()),
        ("structural_fingerprint", StructuralFingerprintConfig(max_bonds=2)),
    ],
)
def test_molecule_ast_fingerprint_keyword_error(ethanol, method_name, config):
    with pytest.raises(TypeError):
        getattr(ethanol, method_name)(config)


@pytest.mark.parametrize(
    "method_name",
    [
        "hashed_fingerprint",
        "counted_hashed_fingerprint",
        "structural_fingerprint",
    ],
)
def test_molecule_ast_fingerprint_required_error(ethanol, method_name):
    with pytest.raises(TypeError):
        getattr(ethanol, method_name)()


@pytest.mark.parametrize(
    ("method_name", "kwargs"),
    [
        ("hashed_fingerprint", {"config": HashedFingerprintConfig.Morgan()}),
        (
            "counted_hashed_fingerprint",
            {"config": HashedFingerprintConfig.Morgan()},
        ),
        ("pattern_fingerprint", {}),
        (
            "structural_fingerprint",
            {"config": StructuralFingerprintConfig(max_bonds=2)},
        ),
    ],
)
def test_molecule_ast_fingerprint_error(undetermined_molecule, method_name, kwargs):
    with pytest.raises(
        UnderdeterminedError,
        match="fingerprint requires a determined molecule",
    ):
        getattr(undetermined_molecule, method_name)(**kwargs)


def test_reaction_ast_combined_fingerprint_difference(ethanol_deoxygenation):
    expected = [
        (864662311, -1),
        (1535166686, -1),
        (2245384272, -1),
        (2246997334, 1),
        (3542456614, -1),
        (3548082732, 1),
        (4018048386, -1),
    ]

    result = ethanol_deoxygenation.combined_fingerprint(
        config=ReactionCombinedFingerprintConfig.Difference(
            molecule=HashedFingerprintConfig.Morgan(
                ring_config=RingConfig(
                    simple_cycle_algorithm=(
                        SimpleCycleEnumerationAlgorithm.ReadTarjan()
                    ),
                    relevant_cycle_algorithm=(
                        RelevantCycleEnumerationAlgorithm.Vismara()
                    ),
                )
            )
        )
    )
    features = result.features

    assert type(result) is ReactionCombinedFingerprint
    assert type(features) is SignedHashedFeatureSet
    assert features.id_width == 64
    assert features.entries == expected
    assert list(features) == expected
    assert all(type(identifier) is int for identifier, _ in features)
    assert all(type(count) is int for _, count in features)


def test_reaction_ast_combined_fingerprint_disjoint_union(ethanol_deoxygenation):
    expected = [
        (ReactionSide.Reactant, 864662311),
        (ReactionSide.Reactant, 1535166686),
        (ReactionSide.Reactant, 2245384272),
        (ReactionSide.Reactant, 2246728737),
        (ReactionSide.Reactant, 3542456614),
        (ReactionSide.Reactant, 4018048386),
        (ReactionSide.Product, 2246728737),
        (ReactionSide.Product, 2246997334),
        (ReactionSide.Product, 3548082732),
    ]

    result = ethanol_deoxygenation.combined_fingerprint(
        config=ReactionCombinedFingerprintConfig.DisjointUnion(
            molecule=HashedFingerprintConfig.Morgan()
        )
    )
    features = result.features

    assert type(result) is ReactionCombinedFingerprint
    assert type(features) is RoleTaggedHashedFeatureSet
    assert features.id_width == 64
    assert features.ids == expected
    assert list(features) == expected
    assert all(type(side) is ReactionSide for side, _ in features)
    assert all(type(identifier) is int for _, identifier in features)


def test_reaction_ast_combined_fingerprint_feature_types(ethanol_deoxygenation):
    molecule_config = HashedFingerprintConfig.Morgan()
    molecular_binary = ethanol_deoxygenation.lhs.hashed_fingerprint(
        config=molecule_config
    )
    molecular_counted = ethanol_deoxygenation.lhs.counted_hashed_fingerprint(
        config=molecule_config
    )
    reaction_signed = ethanol_deoxygenation.combined_fingerprint(
        config=ReactionCombinedFingerprintConfig.Difference(molecule=molecule_config)
    ).features
    reaction_tagged = ethanol_deoxygenation.combined_fingerprint(
        config=ReactionCombinedFingerprintConfig.DisjointUnion(molecule=molecule_config)
    ).features

    assert type(molecular_binary) is HashedFeatureSet
    assert type(molecular_counted) is CountedHashedFeatureSet
    assert type(reaction_signed) is SignedHashedFeatureSet
    assert type(reaction_tagged) is RoleTaggedHashedFeatureSet
    assert molecular_binary != reaction_tagged
    assert molecular_counted != reaction_signed
    with pytest.raises(
        TypeError,
        match=(
            "^'RoleTaggedHashedFeatureSet' object is not an instance of "
            "'HashedFeatureSet'\\nwhile processing 'other'$"
        ),
    ):
        molecular_binary.is_subset(reaction_tagged)


def test_reaction_ast_combined_fingerprint_required_error(ethanol_deoxygenation):
    with pytest.raises(
        TypeError,
        match=(
            "^ReactionAst.combined_fingerprint\\(\\) missing 1 required keyword "
            "argument: 'config'$"
        ),
    ):
        ethanol_deoxygenation.combined_fingerprint()


def test_reaction_ast_combined_fingerprint_keyword_error(ethanol_deoxygenation):
    config = ReactionCombinedFingerprintConfig.Difference(
        molecule=HashedFingerprintConfig.Morgan()
    )

    with pytest.raises(
        TypeError,
        match=(
            "^ReactionAst.combined_fingerprint\\(\\) takes 0 positional arguments "
            "but 1 was given$"
        ),
    ):
        ethanol_deoxygenation.combined_fingerprint(config)
