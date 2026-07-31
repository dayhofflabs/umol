import inspect

import pytest
import umol
import umol._native as native


PUBLIC_EXPORTS = frozenset(
    """
    AromaticSystemAst
    AromaticSystemConstraintAst
    AromaticSystemConstraintKey
    AromaticSystemConstraintsAst
    AromaticSystemConstraintsView
    AromaticSystemDelta
    AromaticSystemFieldChange
    AromaticSystemView
    AromaticSystemViews
    AromaticValence
    AromaticValenceAst
    AromaticityConfig
    AromaticityModel
    AromaticityResolveConfig
    AtomAst
    AtomConstraintAst
    AtomConstraintKey
    AtomConstraintsAst
    AtomConstraintsView
    AtomDelta
    AtomFieldChange
    AtomRingSizeCounts
    AtomTypeRegistry
    AtomView
    AtomViews
    AutomorphismAlgorithm
    BitFp
    BondAst
    BondConstraintAst
    BondConstraintKey
    BondConstraintsAst
    BondConstraintsView
    BondDelta
    BondFieldChange
    BondRingSizeCounts
    BondView
    BondViews
    BooleanAst
    ChemistryModel
    CisTransConfiguration
    CisTransStereo
    CisTransStereoAst
    CommonSubgraphEnumerationAlgorithm
    ConnectedComponentsAlgorithm
    Constraint
    ConstraintDelta
    Constraints
    ConstraintsView
    ContradictionError
    Correspondence
    CountedHashedFeatureSet
    DativeBondAst
    DativeBondConstraintAst
    DativeBondConstraintKey
    DativeBondConstraintsAst
    DativeBondConstraintsView
    DativeBondDelta
    DativeBondFieldChange
    DativeBondRingSizeCounts
    DativeBondView
    DativeBondViews
    Delta
    Deltas
    E
    EcfpHashScheme
    ElectronCountsAst
    Element
    ElementAst
    ElementScope
    Entity
    FluxionalityAst
    HashedFeatureSet
    HashedFingerprintConfig
    InconsistencyPolicy
    InvalidStructureError
    IsotopeMass
    IsotopeMassAst
    LigandPermutation
    LigandSymmetryAst
    MaximumIndependentSetAlgorithm
    MemOp
    MetadataError
    ModelConversionError
    MoleculeAst
    MoleculeConstraint
    MoleculeCorrespondence
    MoleculeDefaults
    MoleculeMetadata
    MulticenterBondAst
    MulticenterBondConstraintAst
    MulticenterBondConstraintKey
    MulticenterBondConstraintsAst
    MulticenterBondConstraintsView
    MulticenterBondDelta
    MulticenterBondFieldChange
    MulticenterBondView
    MulticenterBondViews
    MulticenterValence
    MulticenterValenceAst
    NoncovalentBondAst
    NoncovalentBondConstraintAst
    NoncovalentBondConstraintKey
    NoncovalentBondConstraintsAst
    NoncovalentBondConstraintsView
    NoncovalentBondDelta
    NoncovalentBondFieldChange
    NoncovalentBondKind
    NoncovalentBondKindAst
    NoncovalentBondView
    NoncovalentBondViews
    Orientation
    OrientedLigandPermutation
    ParseError
    PatternFingerprintConfig
    Permutation
    ReactionApplicationConfig
    ReactionAst
    ReactionCombinedFingerprint
    ReactionCombinedFingerprintConfig
    ReactionCompositionConfig
    ReactionDefaults
    ReactionDerivation
    ReactionMetadata
    ReactionSide
    RefinementRounds
    RelOp
    RelationalConstraint
    RelevantCycleEnumerationAlgorithm
    ResolveConfig
    RingConfig
    RingLimits
    RingMembershipAst
    RingScope
    RoleTaggedHashedFeatureSet
    SignedHashedFeatureSet
    SimpleCycleEnumerationAlgorithm
    SmilesIoConfig
    SmilesSyntaxFlags
    SpinState
    StereoAtomAst
    StereoAtomConstraintAst
    StereoAtomConstraintKey
    StereoAtomConstraintsAst
    StereoAtomConstraintsView
    StereoAtomDelta
    StereoAtomFieldChange
    StereoAtomView
    StereoAtomViews
    StereoBondAst
    StereoBondConstraintAst
    StereoBondConstraintKey
    StereoBondConstraintsAst
    StereoBondConstraintsView
    StereoBondDelta
    StereoBondFieldChange
    StereoBondView
    StereoBondViews
    StereoConfigurationAst
    StereoCoset
    StereoKind
    StereoKindModel
    StereoLigand
    StereoLigandKind
    StereoLigandPair
    StereoModel
    StereoResolveConfig
    StereoTerm
    Stereogenicity
    StereogenicityAst
    StructuralFeatureSet
    StructuralFingerprintConfig
    SubPatternAnchor
    SubgraphEnumerationAlgorithm
    SubgraphIsomorphismAlgorithm
    SubstructureMatchAlgorithm
    SubstructureSearchConfig
    TetrahedralConfiguration
    TetrahedralStereo
    TetrahedralStereoAst
    Topicity
    TopicityAst
    TopicityRelationAst
    UnpairedElectrons
    UnpairedElectronsAst
    UnderdeterminedError
    ValenceEntry
    ValenceModel
    ValenceTable
    ValueAst
    ValuePredicate
    ValueTerm
    WlHashScheme
    __version__
    """.split()
)


def test_package_metadata():
    assert isinstance(umol.__version__, str)
    assert umol.__version__ != "0.0.0"
    assert umol._native is native


def test_public_exports():
    native_exports = frozenset(
        name for name in vars(native) if not name.startswith("_")
    )
    native_package_exports = PUBLIC_EXPORTS - {"E", "__version__"}

    assert frozenset(umol.__all__) == PUBLIC_EXPORTS
    assert len(umol.__all__) == len(PUBLIC_EXPORTS)
    assert native_exports == native_package_exports
    assert {name: getattr(umol, name) for name in native_package_exports} == {
        name: getattr(native, name) for name in native_package_exports
    }


@pytest.mark.parametrize(
    "name",
    [
        "BridgitConfig",
        "CycleEnumerationAlgorithm",
        "DrfpConfig",
        "ReactionApplicationIter",
        "ReactionCombinator",
        "SmilesLintConfig",
        "SmilesLintFlags",
    ],
)
def test_deferred_export(name):
    assert name not in umol.__all__
    assert not hasattr(umol, name)
    assert not hasattr(native, name)


@pytest.mark.parametrize(
    ("owner", "name"),
    [
        (umol.WlHashScheme, "seed"),
        (umol.WlHashScheme, "from_seed"),
        (umol.EcfpHashScheme, "seed"),
        (umol.EcfpHashScheme, "from_seed"),
        (umol.HashedFeatureSet, "__array__"),
        (umol.HashedFeatureSet, "to_numpy"),
        (umol.CountedHashedFeatureSet, "__array__"),
        (umol.CountedHashedFeatureSet, "to_numpy"),
        (umol.BitFp, "__array__"),
        (umol.BitFp, "to_numpy"),
        (umol.StructuralFeatureSet, "__array__"),
        (umol.StructuralFeatureSet, "to_numpy"),
        (umol.SignedHashedFeatureSet, "__array__"),
        (umol.SignedHashedFeatureSet, "to_numpy"),
        (umol.RoleTaggedHashedFeatureSet, "__array__"),
        (umol.RoleTaggedHashedFeatureSet, "to_numpy"),
    ],
)
def test_deferred_member(owner, name):
    assert not hasattr(owner, name)


@pytest.mark.parametrize(
    ("operation", "expected"),
    [
        (
            umol.MoleculeAst.from_smiles,
            "(source, *, io_config=None, chemistry_model=None, resolve_config=None)",
        ),
        (umol.MoleculeAst.combine, "(self, /, other)"),
        (umol.MoleculeAst.combine_from, "(self, /, other)"),
        (umol.MoleculeAst.combine_all, "(molecules)"),
        (umol.MoleculeAst.split, "(self, /)"),
        (
            umol.MoleculeAst.substructure_matches,
            "(self, /, host, *, config=None)",
        ),
        (umol.MoleculeAst.hashed_fingerprint, "(self, /, *, config)"),
        (umol.MoleculeAst.counted_hashed_fingerprint, "(self, /, *, config)"),
        (umol.MoleculeAst.pattern_fingerprint, "(self, /, *, config=None)"),
        (umol.MoleculeAst.structural_fingerprint, "(self, /, *, config)"),
        (
            umol.ReactionAst.from_reaction_smiles,
            "(source, *, io_config=None, chemistry_model=None, resolve_config=None)",
        ),
        (
            umol.ReactionAst.compose,
            "(self, /, other, *, config=None)",
        ),
        (umol.ReactionAst.apply, "(self, /, host, *, config=None)"),
        (umol.ReactionAst.combined_fingerprint, "(self, /, *, config)"),
    ],
)
def test_public_operation_signature(operation, expected):
    assert str(inspect.signature(operation)) == expected


@pytest.mark.parametrize(
    ("constructor", "expected"),
    [
        (umol.SmilesIoConfig.opensmiles, "()"),
        (umol.SmilesIoConfig.lenient, "()"),
        (umol.SmilesIoConfig.with_syntax_flags, "(*, syntax_flags)"),
        (umol.RefinementRounds.Fixed, "(*, rounds)"),
        (umol.RefinementRounds.ToFixpoint, "()"),
        (
            umol.HashedFingerprintConfig.Morgan,
            "(*, radius=2, ring_config=Ellipsis)",
        ),
        (
            umol.HashedFingerprintConfig.Ecfp,
            "(*, radius, hashing_scheme=Ellipsis, ring_config=Ellipsis)",
        ),
        (
            umol.HashedFingerprintConfig.Wl,
            "(*, rounds, hashing_scheme=Ellipsis)",
        ),
        (
            umol.ReactionCombinedFingerprintConfig.Difference,
            "(*, molecule)",
        ),
        (
            umol.ReactionCombinedFingerprintConfig.DisjointUnion,
            "(*, molecule)",
        ),
        (
            umol.PatternFingerprintConfig,
            "(*, width=2048, match_algorithm=Ellipsis, "
            "subgraph_isomorphism_algorithm=Ellipsis)",
        ),
        (
            umol.StructuralFingerprintConfig,
            "(*, max_bonds, subgraph_enumeration_algorithm=Ellipsis, "
            "automorphism_algorithm=Ellipsis)",
        ),
        (
            umol.RingConfig,
            "(*, simple_cycle_algorithm=None, relevant_cycle_algorithm=None)",
        ),
        (
            umol.AromaticityConfig,
            "(*, ring_config=Ellipsis, connected_components_algorithm=Ellipsis, "
            "maximum_independent_set_algorithm=Ellipsis)",
        ),
        (
            umol.AromaticityResolveConfig,
            "(*, perception=Ellipsis, delocalize_charge=True, "
            "reset_aromatic_valence=False)",
        ),
        (
            umol.StereoResolveConfig,
            "(*, reset_stereo_constraints=False, inconsistency=Ellipsis)",
        ),
        (umol.ResolveConfig, "(*, aromaticity, stereo)"),
        (umol.ChemistryModel, "(*, valence, aromaticity, stereo)"),
        (
            umol.ReactionApplicationConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis)",
        ),
        (
            umol.ReactionCompositionConfig,
            "(*, common_subgraph_enumeration_algorithm=Ellipsis)",
        ),
        (
            umol.SubstructureSearchConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis)",
        ),
    ],
)
def test_public_constructor_signature(constructor, expected):
    assert str(inspect.signature(constructor)) == expected


@pytest.mark.parametrize(
    "error_type",
    [
        umol.ContradictionError,
        umol.InvalidStructureError,
        umol.ModelConversionError,
        umol.ParseError,
        umol.UnderdeterminedError,
    ],
)
def test_error_import(error_type):
    error = error_type("diagnostic")

    assert getattr(umol, error_type.__name__) is error_type
    assert isinstance(error, Exception)
    assert str(error) == "diagnostic"
