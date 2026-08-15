import inspect

import pytest
import umol
import umol._native as native


PUBLIC_EXPORTS = frozenset(
    """
    AromaticBondConstraintMismatchPolicy
    AromaticSystemForm
    AromaticSystemUpdate
    AromaticSystemConstraintForm
    AromaticSystemConstraintKey
    AromaticSystemConstraintsForm
    AromaticSystemConstraintsView
    AromaticSystemDelta
    AromaticSystemFieldChange
    AromaticSystemView
    AromaticSystemViews
    AromaticValence
    AromaticValenceForm
    AromaticityConfig
    AromaticityFailurePolicy
    AromaticityMismatchPolicy
    AromaticityModel
    AromaticityRule
    AromaticityTieBreak
    AromaticityResolveConfig
    AtomForm
    AtomUpdate
    AtomCompletions
    AtomConstraintForm
    AtomConstraintKey
    AtomConstraintsForm
    AtomConstraintsView
    AtomDelta
    AtomFieldChange
    AtomRingSizeCounts
    AtomTypeRegistry
    AtomView
    AtomViews
    AutomorphismAlgorithm
    BitFp
    BondForm
    BondUpdate
    BondConstraintForm
    BondConstraintKey
    BondConstraintsForm
    BondConstraintsView
    BondDelta
    BondFieldChange
    BondRingSizeCounts
    BondView
    BondViews
    BooleanForm
    CanonicalizationConfig
    CanonicalizationLevel
    ChemistryModel
    CisTransConfiguration
    CisTransStereo
    CisTransStereoForm
    CommonSubgraphEnumerationAlgorithm
    ConnectedComponentsAlgorithm
    ConnectivityModel
    Constraint
    ConstraintDelta
    ConstraintEdit
    Constraints
    ConstraintsView
    ContradictionError
    Correspondence
    CountedHashedFeatureSet
    DativeBondForm
    DativeBondUpdate
    DativeBondConstraintForm
    DativeBondConstraintKey
    DativeBondConstraintsForm
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
    ElectronCountsForm
    Edit
    Edits
    Element
    ElementForm
    ElementScope
    Entity
    FluxionalityForm
    HashedFeatureSet
    HashedFingerprintConfig
    InvalidStructureError
    IsotopeMass
    IsotopeMassForm
    LigandPermutation
    LigandSymmetryForm
    MaximumIndependentSetAlgorithm
    MemOp
    MetadataError
    ModelConversionError
    Molecule
    MoleculeConstraint
    MoleculeCorrespondence
    MoleculeDefaults
    MoleculeEditor
    MoleculeMetadata
    MulticenterBondForm
    MulticenterBondUpdate
    MulticenterBondConstraintForm
    MulticenterBondConstraintKey
    MulticenterBondConstraintsForm
    MulticenterBondConstraintsView
    MulticenterBondDelta
    MulticenterBondFieldChange
    MulticenterBondView
    MulticenterBondViews
    MulticenterValence
    MulticenterValenceForm
    New
    NoncovalentBondForm
    NoncovalentBondConstraintForm
    NoncovalentBondConstraintKey
    NoncovalentBondConstraintsForm
    NoncovalentBondConstraintsView
    NoncovalentBondDelta
    NoncovalentBondFieldChange
    NoncovalentBondKind
    NoncovalentBondKindForm
    NoncovalentBondUpdate
    NoncovalentBondView
    NoncovalentBondViews
    Orientation
    OrientedLigandPermutation
    ParseError
    PatternFingerprintConfig
    Permutation
    ReactionApplicationConfig
    Reaction
    ReactionCombinedFingerprint
    ReactionCombinedFingerprintConfig
    ReactionCompositionConfig
    ReactionDefaults
    ReactionDerivation
    ReactionMetadata
    ReactionSide
    ReactionSpan
    RefinementRounds
    RelOp
    RelationalConstraint
    RelevantCycleEnumerationAlgorithm
    ResolveConfig
    ResolveReport
    RingConfig
    RingLimits
    RingMembershipForm
    RingScope
    RoleTaggedHashedFeatureSet
    SignedHashedFeatureSet
    SimpleCycleEnumerationAlgorithm
    SmilesIoConfig
    SmilesSyntaxFlags
    SpinState
    StereoAtomForm
    StereoAtomUpdate
    StereoAtomConstraintForm
    StereoAtomConstraintKey
    StereoAtomConstraintsForm
    StereoAtomConstraintsView
    StereoAtomDelta
    StereoAtomFieldChange
    StereoAtomView
    StereoAtomViews
    StereoBondForm
    StereoBondUpdate
    StereoBondConstraintForm
    StereoBondConstraintKey
    StereoBondConstraintsForm
    StereoBondConstraintsView
    StereoBondDelta
    StereoBondFieldChange
    StereoBondView
    StereoBondViews
    StereoConfigurationForm
    StereoConfigurationUpdate
    StereoCoset
    StereoFailurePolicy
    StereoKind
    StereoKindModel
    StereoLigand
    StereoLigandKind
    StereoLigandPair
    StereoModel
    StereoMismatchPolicy
    StereoResolveConfig
    StereoTerm
    Stereogenicity
    StereogenicityForm
    StructuralFeatureSet
    StructuralFingerprintConfig
    SubgraphEnumerationAlgorithm
    SubgraphIsomorphismAlgorithm
    SubstructureMatchAlgorithm
    SubstructureSearchConfig
    TetrahedralConfiguration
    TetrahedralStereo
    TetrahedralStereoForm
    Topicity
    TopicityForm
    TopicityRelationForm
    Transaction
    TransactionError
    UnpairedElectrons
    UnpairedElectronsForm
    UnpairedElectronsUpdate
    UnderdeterminedError
    ValenceCandidateSource
    ValenceEntry
    ValenceModel
    ValenceTable
    ValenceTieBreak
    NumForm
    PredExpr
    ArithExpr
    WlHashScheme
    __version__
    """.split()
)

EDITING_EXPORTS = frozenset(
    {
        "AromaticSystemUpdate",
        "AtomUpdate",
        "BondUpdate",
        "ConstraintEdit",
        "DativeBondUpdate",
        "Edit",
        "Edits",
        "MulticenterBondUpdate",
        "New",
        "NoncovalentBondUpdate",
        "StereoAtomUpdate",
        "StereoBondUpdate",
        "StereoConfigurationUpdate",
        "UnpairedElectronsUpdate",
    }
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


def test_editing_exports():
    package_exports = frozenset(name for name in EDITING_EXPORTS if name in umol.__all__)
    native_exports = frozenset(name for name in EDITING_EXPORTS if hasattr(native, name))

    assert package_exports == EDITING_EXPORTS
    assert native_exports == EDITING_EXPORTS
    assert {name: getattr(umol, name) for name in EDITING_EXPORTS} == {
        name: getattr(native, name) for name in EDITING_EXPORTS
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
            umol.Molecule.from_smiles,
            "(source, *, io_config=None, chemistry_model=None, resolve_config=None)",
        ),
        (umol.Molecule.edit, "(self, /)"),
        (umol.Molecule.apply, "(self, /, edits)"),
        (umol.Molecule.combine, "(self, /, other)"),
        (umol.Molecule.combine_from, "(self, /, other)"),
        (umol.Molecule.combine_all, "(molecules)"),
        (umol.Molecule.react, "(self, /, reaction, *, config=None)"),
        (umol.Molecule.react_all, "(reactants, reaction, *, config=None)"),
        (umol.Molecule.split, "(self, /)"),
        (
            umol.Molecule.canonicalize,
            "(self, /, *, stereo_model=None, config=None)",
        ),
        (
            umol.Molecule.canonicalize_by,
            "(self, /, level, *, stereo_model=None, config=None)",
        ),
        (
            umol.Molecule.canonical_eq,
            "(self, /, other, *, stereo_model=None, config=None)",
        ),
        (
            umol.Molecule.canonical_eq_by,
            "(self, /, other, level, *, stereo_model=None, config=None)",
        ),
        (
            umol.Molecule.substructure_matches,
            "(self, /, host, *, config=None)",
        ),
        (umol.Molecule.hashed_fingerprint, "(self, /, *, config)"),
        (umol.Molecule.counted_hashed_fingerprint, "(self, /, *, config)"),
        (umol.Molecule.pattern_fingerprint, "(self, /, *, config=None)"),
        (umol.Molecule.structural_fingerprint, "(self, /, *, config)"),
        (umol.MoleculeEditor.snapshot, "(self, /)"),
        (umol.MoleculeEditor.build, "(self, /)"),
        (umol.MoleculeEditor.transact, "(self, /, edits)"),
        (
            umol.Reaction.from_reaction_smiles,
            "(source, *, io_config=None, chemistry_model=None, resolve_config=None)",
        ),
        (
            umol.Reaction.compose,
            "(self, /, other, *, config=None)",
        ),
        (umol.Reaction.apply, "(self, /, host, *, config=None)"),
        (umol.Reaction.combined_fingerprint, "(self, /, *, config)"),
        (
            umol.Reaction.canonicalize,
            "(self, /, *, stereo_model=None, config=None)",
        ),
        (
            umol.Reaction.canonicalize_by,
            "(self, /, level, *, stereo_model=None, config=None)",
        ),
        (
            umol.Reaction.canonical_eq,
            "(self, /, other, *, stereo_model=None, config=None)",
        ),
        (
            umol.Reaction.canonical_eq_by,
            "(self, /, other, level, *, stereo_model=None, config=None)",
        ),
        (
            umol.ReactionSpan.canonicalize,
            "(self, /, *, stereo_model=None, config=None)",
        ),
        (
            umol.ReactionSpan.canonicalize_by,
            "(self, /, level, *, stereo_model=None, config=None)",
        ),
        (
            umol.ReactionSpan.canonical_eq,
            "(self, /, other, *, stereo_model=None, config=None)",
        ),
        (
            umol.ReactionSpan.canonical_eq_by,
            "(self, /, other, level, *, stereo_model=None, config=None)",
        ),
        (umol.Transaction.rollback, "(self, /, editor)"),
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
            "subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
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
            "(*, perception=Ellipsis, aromatic_valence_failure=Ellipsis, "
            "aromatic_system_failure=Ellipsis, "
            "aromatic_valence_mismatch=Ellipsis, "
            "aromatic_bond_constraint_mismatch=Ellipsis, "
            "reset_aromatic_valence=False)",
        ),
        (
            umol.StereoResolveConfig,
            "(*, tetrahedral_stereo_failure=Ellipsis, "
            "stereo_atom_failure=Ellipsis, "
            "tetrahedral_stereo_mismatch=Ellipsis, "
            "cis_trans_stereo_failure=Ellipsis, "
            "stereo_bond_failure=Ellipsis, "
            "cis_trans_stereo_mismatch=Ellipsis, "
            "reset_stereo_constraints=False)",
        ),
        (umol.ResolveConfig, "(*, aromaticity, stereo)"),
        (
            umol.ConnectivityModel,
            "(*, allow_disconnected, allow_disconnected_dative, "
            "allow_disconnected_aromatic, allow_disconnected_multicenter, "
            "allow_disconnected_noncovalent, allow_disconnected_stereo_atom, "
            "allow_disconnected_stereo_bond, allow_disconnected_constraints)",
        ),
        (
            umol.ChemistryModel,
            "(*, connectivity, valence, aromaticity, stereo)",
        ),
        (
            umol.ReactionApplicationConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
        ),
        (
            umol.ReactionCompositionConfig,
            "(*, common_subgraph_enumeration_algorithm=Ellipsis)",
        ),
        (
            umol.CanonicalizationConfig,
            "(*, automorphism_algorithm=Ellipsis)",
        ),
        (
            umol.SubstructureSearchConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
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
        umol.TransactionError,
        umol.UnderdeterminedError,
    ],
)
def test_error_import(error_type):
    error = error_type("diagnostic")

    assert getattr(umol, error_type.__name__) is error_type
    assert isinstance(error, Exception)
    assert str(error) == "diagnostic"
