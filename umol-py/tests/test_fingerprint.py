import pytest

from umol import (
    BitFp,
    CountedHashedFeatureSet,
    EcfpHashScheme,
    HashedFingerprintConfig,
    HashedFeatureSet,
    PatternFingerprintConfig,
    ReactionCombinedFingerprint,
    ReactionCombinedFingerprintConfig,
    RefinementRounds,
    RoleTaggedHashedFeatureSet,
    SignedHashedFeatureSet,
    StructuralFeatureSet,
    StructuralFingerprintConfig,
    WlHashScheme,
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
            HashedFingerprintConfig.Morgan(radius=2),
            "HashedFingerprintConfig.Morgan(radius=2)",
        ),
        (
            HashedFingerprintConfig.Ecfp(),
            HashedFingerprintConfig.Ecfp(
                radius=2,
                scheme=EcfpHashScheme.Xxh3Width64V1(),
            ),
            "HashedFingerprintConfig.Ecfp(radius=2, "
            "scheme=EcfpHashScheme.Xxh3Width64V1())",
        ),
        (
            HashedFingerprintConfig.Wl(rounds=RefinementRounds.ToFixpoint()),
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.ToFixpoint(),
                scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            "HashedFingerprintConfig.Wl(rounds=RefinementRounds.ToFixpoint(), "
            "scheme=WlHashScheme.Xxh3SortedWidth64V1())",
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
            HashedFingerprintConfig.Morgan(radius=3),
            HashedFingerprintConfig.Morgan(radius=3),
            "HashedFingerprintConfig.Morgan(radius=3)",
        ),
        (
            HashedFingerprintConfig.Ecfp(
                radius=3,
                scheme=EcfpHashScheme.Xxh3Width64V1(),
            ),
            HashedFingerprintConfig.Ecfp(
                radius=3,
                scheme=EcfpHashScheme.Xxh3Width64V1(),
            ),
            "HashedFingerprintConfig.Ecfp(radius=3, "
            "scheme=EcfpHashScheme.Xxh3Width64V1())",
        ),
        (
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
                scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            HashedFingerprintConfig.Wl(
                rounds=RefinementRounds.Fixed(rounds=3),
                scheme=WlHashScheme.Xxh3SortedWidth64V1(),
            ),
            "HashedFingerprintConfig.Wl(rounds=RefinementRounds.Fixed(rounds=3), "
            "scheme=WlHashScheme.Xxh3SortedWidth64V1())",
        ),
    ],
)
def test_hashed_fingerprint_config(value, expected, expected_repr):
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (
            PatternFingerprintConfig(),
            PatternFingerprintConfig(width=2048),
            "PatternFingerprintConfig(width=2048)",
        ),
        (
            PatternFingerprintConfig(width=512),
            PatternFingerprintConfig(width=512),
            "PatternFingerprintConfig(width=512)",
        ),
    ],
)
def test_pattern_fingerprint_config(value, expected, expected_repr):
    assert value.width == expected.width
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize("width", [0, -1])
def test_pattern_fingerprint_config_error(width):
    with pytest.raises(ValueError, match="width must be positive"):
        PatternFingerprintConfig(width=width)


@pytest.mark.parametrize("max_bonds", [0, 3])
def test_structural_fingerprint_config(max_bonds):
    value = StructuralFingerprintConfig(max_bonds=max_bonds)

    assert value.max_bonds == max_bonds
    assert value == StructuralFingerprintConfig(max_bonds=max_bonds)
    assert repr(value) == f"StructuralFingerprintConfig(max_bonds={max_bonds})"


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
