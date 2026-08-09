import pytest

from umol import (
    AromaticBondConstraintMismatchPolicy,
    AromaticValenceAst,
    AromaticityConfig,
    AromaticityFailurePolicy,
    AromaticityMismatchPolicy,
    AromaticityResolveConfig,
    AtomAst,
    AtomFieldChange,
    AtomUpdate,
    AutomorphismAlgorithm,
    BondAst,
    ChemistryModel,
    ConnectedComponentsAlgorithm,
    Correspondence,
    Edit,
    Edits,
    HashedFingerprintConfig,
    InvalidStructureError,
    MaximumIndependentSetAlgorithm,
    MoleculeAst,
    ParseError,
    PatternFingerprintConfig,
    ReactionApplicationConfig,
    ReactionAst,
    ReactionCombinedFingerprintConfig,
    RefinementRounds,
    RelevantCycleEnumerationAlgorithm,
    ResolveConfig,
    RingConfig,
    SimpleCycleEnumerationAlgorithm,
    SmilesIoConfig,
    StereoFailurePolicy,
    StereoResolveConfig,
    StructuralFingerprintConfig,
    SubgraphEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
    SubstructureSearchConfig,
    TransactionError,
    UnderdeterminedError,
    NumForm,
)


def test_resolved_smiles_workflow():
    io_config = SmilesIoConfig.lenient()
    chemistry_model = ChemistryModel.default()
    resolve_config = ResolveConfig(
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
                connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(),
                maximum_independent_set_algorithm=(
                    MaximumIndependentSetAlgorithm.BranchAndBound()
                ),
            ),
            aromatic_valence_failure=AromaticityFailurePolicy.Error,
            aromatic_system_failure=AromaticityFailurePolicy.Error,
            aromatic_valence_mismatch=AromaticityMismatchPolicy.Error,
            aromatic_bond_constraint_mismatch=(
                AromaticBondConstraintMismatchPolicy.Error
            ),
            reset_aromatic_valence=False,
        ),
        stereo=StereoResolveConfig(reset_stereo_constraints=False),
    )

    molecule = MoleculeAst.from_smiles(
        "[cH+]1[cH][cH]1",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    independent = MoleculeAst.from_smiles(
        "[cH+]1[cH][cH]1",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )

    assert [atom.charge for atom in molecule.atoms] == [
        NumForm.Lit(1),
        NumForm.Lit(0),
        NumForm.Lit(0),
    ]
    assert [atom.constraints.aromatic_valence for atom in molecule.atoms] == [
        AromaticValenceAst.Aromatic(NumForm.Lit(0)),
        AromaticValenceAst.Aromatic(NumForm.Lit(1)),
        AromaticValenceAst.Aromatic(NumForm.Lit(1)),
    ]
    assert [
        (system.atom_ids, system.charge) for system in molecule.aromatic_systems
    ] == [((0, 1, 2), NumForm.Lit(0))]
    assert list(molecule.stereo_atoms) == []
    assert molecule == independent
    assert io_config == SmilesIoConfig.lenient()
    assert chemistry_model == ChemistryModel.default()
    assert (
        resolve_config.aromaticity.aromatic_valence_failure
        == AromaticityFailurePolicy.Error
    )
    assert (
        resolve_config.stereo.tetrahedral_stereo_failure
        == StereoFailurePolicy.Error
    )

    molecule.atoms[0].charge = 3

    assert independent.atoms[0].charge == NumForm.Lit(1)

    with pytest.raises(ParseError, match="^Invalid token at position 2$"):
        MoleculeAst.from_smiles(
            "C->N",
            io_config=SmilesIoConfig.opensmiles(),
            chemistry_model=chemistry_model,
            resolve_config=resolve_config,
        )


def test_molecule_editing_workflow():
    molecule = MoleculeAst.parse('{:atoms ["N#h3"]}')
    original = MoleculeAst.parse('{:atoms ["N#h3"]}')
    expected = MoleculeAst.parse(
        '{:atoms ["N#h2" "C#h3"] :bonds [[0 1 "1"]]}'
    )
    edits = Edits()
    edits.update_atom(
        0,
        AtomAst.parse("N#h3"),
        AtomUpdate(implicit_hydrogens=2),
    )
    methyl = edits.add_atom(AtomAst.parse("C#h3"))
    edits.add_bond(0, methyl, BondAst(1))

    rendered = edits.render()
    parsed = Edits.parse(rendered)
    applied = molecule.apply(parsed)

    assert rendered == (
        '[{:atom {:modify [0 {:expect "#h3" :update "#h2"}]}} '
        '{:atom {:add "C#h3"}} '
        '{:bond {:add [0 {:new 0} :single]}}]'
    )
    assert parsed == edits
    assert parsed.render() == rendered
    assert applied == expected
    assert molecule == original

    editor = molecule.edit()
    transaction = editor.transact(parsed)

    assert editor.snapshot() == expected
    assert molecule == original

    transaction.rollback(editor)

    assert editor.snapshot() == original
    assert editor.build() == original
    assert molecule == original

    failing = Edits()
    failing.add_atom(AtomAst.parse("O"))
    failing.append(
        Edit.ModifyAtomField(
            id=7,
            change=AtomFieldChange.Charge(
                old=NumForm.Lit(0),
                new=NumForm.Lit(1),
            ),
        )
    )

    with pytest.raises(
        TransactionError,
        match="^atom handle 7 is out of range for 1 entries$",
    ):
        molecule.apply(failing)

    assert molecule == original


def test_fingerprint_workflow():
    molecule = MoleculeAst.from_smiles("CO")
    product = MoleculeAst.from_smiles("C")
    molecule_snapshot = MoleculeAst.from_smiles("CO")
    product_snapshot = MoleculeAst.from_smiles("C")
    hashed_config = HashedFingerprintConfig.Morgan(radius=0)
    pattern_config = PatternFingerprintConfig(
        width=16,
        match_algorithm=SubstructureMatchAlgorithm.Incidence(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Ullmann(),
    )
    structural_config = StructuralFingerprintConfig(
        max_bonds=0,
        subgraph_enumeration_algorithm=SubgraphEnumerationAlgorithm.Esu(),
        automorphism_algorithm=AutomorphismAlgorithm.Nauty(),
    )

    hashed = molecule.hashed_fingerprint(config=hashed_config)
    counted = molecule.counted_hashed_fingerprint(config=hashed_config)
    pattern = molecule.pattern_fingerprint(config=pattern_config)
    structural = molecule.structural_fingerprint(config=structural_config)
    reaction = ReactionAst.from_sides(
        molecule,
        product,
        Correspondence([(0, 0)], 2, 1),
    )
    reaction_snapshot = ReactionAst(reaction.lhs, reaction.deltas)
    combined = reaction.combined_fingerprint(
        config=ReactionCombinedFingerprintConfig.Difference(molecule=hashed_config)
    )

    assert hashed.ids == [864662311, 2246728737]
    assert counted.entries == [(864662311, 1), (2246728737, 1)]
    assert [bit for bit in range(pattern.width) if pattern[bit]] == [
        4,
        5,
        6,
        7,
        9,
        13,
        14,
        15,
    ]
    assert structural.keys == [
        bytes.fromhex("01 00 00 00 05 00 00 00 00 06 00 00 00 00 00 00 00"),
        bytes.fromhex("01 00 00 00 05 00 00 00 00 08 00 00 00 00 00 00 00"),
    ]
    assert combined.features.entries == [
        (864662311, -1),
        (2246728737, -1),
        (2246733040, 1),
    ]
    assert molecule == molecule_snapshot
    assert product == product_snapshot
    assert reaction == reaction_snapshot

    hashed.ids.append(7)
    counted.entries.append((7, 2))
    structural.keys.append(b"detached")
    combined_features = combined.features
    combined_features.entries.append((7, 2))

    assert hashed.ids == [864662311, 2246728737]
    assert counted.entries == [(864662311, 1), (2246728737, 1)]
    assert structural.keys == [
        bytes.fromhex("01 00 00 00 05 00 00 00 00 06 00 00 00 00 00 00 00"),
        bytes.fromhex("01 00 00 00 05 00 00 00 00 08 00 00 00 00 00 00 00"),
    ]
    assert combined.features.entries == [
        (864662311, -1),
        (2246728737, -1),
        (2246733040, 1),
    ]

    with pytest.raises(
        UnderdeterminedError,
        match="^fingerprint requires a determined molecule$",
    ):
        MoleculeAst.parse('{:atoms ["C"]}').hashed_fingerprint(
            config=HashedFingerprintConfig.Wl(rounds=RefinementRounds.Fixed(rounds=1))
        )


def test_substructure_workflow():
    pattern = MoleculeAst.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')
    host = MoleculeAst.from_smiles("CCO")
    pattern_snapshot = MoleculeAst.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')
    host_snapshot = MoleculeAst.from_smiles("CCO")
    config = SubstructureSearchConfig(
        match_algorithm=SubstructureMatchAlgorithm.Incidence(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Ullmann(),
    )

    matches = pattern.substructure_matches(host, config=config)

    assert len(matches) == 1
    assert matches[0].atoms.matched_pairs == [(0, 1), (1, 2)]
    assert matches[0].bonds.matched_pairs == [(0, 1)]
    assert matches[0].dative_bonds.matched_pairs == []
    assert matches[0].aromatic_systems.matched_pairs == []
    assert matches[0].multicenter_bonds.matched_pairs == []
    assert matches[0].noncovalent_bonds.matched_pairs == []
    assert matches[0].stereo_atoms.matched_pairs == []
    assert matches[0].stereo_bonds.matched_pairs == []
    assert pattern == pattern_snapshot
    assert host == host_snapshot

    atom_matched_pairs = matches[0].atoms.matched_pairs
    atom_matched_pairs.append((0, 0))
    pattern.atoms[0].charge = 2
    host.atoms[0].charge = 3

    assert matches[0].atoms.matched_pairs == [(0, 1), (1, 2)]
    assert matches[0].bonds.matched_pairs == [(0, 1)]

    with pytest.raises(
        TypeError,
        match=(
            "^MoleculeAst.substructure_matches\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        pattern_snapshot.substructure_matches(host_snapshot, config)


def test_reaction_application_workflow():
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"]]} '
        ':deltas [{:atom {:modify [0 "#h3#v1"]}} '
        "{:atom {:remove 1}} {:bond {:remove 0}}]}"
    )
    host = MoleculeAst.from_smiles("CCO")
    reaction_snapshot = ReactionAst(reaction.lhs, reaction.deltas)
    host_snapshot = MoleculeAst.from_smiles("CCO")
    config = ReactionApplicationConfig(
        match_algorithm=SubstructureMatchAlgorithm.Incidence(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Ullmann(),
    )

    derivations = list(reaction.apply(host, config=config))

    assert len(derivations) == 1
    assert derivations[0].lhs == host_snapshot
    assert derivations[0].rhs == MoleculeAst.from_smiles("CC")
    assert derivations[0].atom_correspondence.matched_pairs == [(0, 0), (1, 1)]
    assert derivations[0].atom_correspondence.left_count == 3
    assert derivations[0].atom_correspondence.right_count == 2
    assert derivations[0].comap.bonds.matched_pairs == [(0, 0)]
    assert reaction == reaction_snapshot
    assert host == host_snapshot

    detached_lhs = derivations[0].lhs
    detached_rhs = derivations[0].rhs
    detached_lhs.atoms[0].charge = 4
    detached_rhs.atoms[0].charge = 5

    assert derivations[0].lhs == host_snapshot
    assert derivations[0].rhs == MoleculeAst.from_smiles("CC")
    assert derivations[0].rhs.hashed_fingerprint(
        config=HashedFingerprintConfig.Morgan(radius=0)
    ).ids == [2246728737]

    invalid_host = ReactionAst.parse(
        '{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"] [0 1 "2"]]} :deltas []}'
    ).lhs

    with pytest.raises(
        InvalidStructureError,
        match=(
            r"^invalid host: bond: parallel bonds on atoms "
            r"\[AtomId\(0\), AtomId\(1\)\]$"
        ),
    ):
        reaction.apply(invalid_host, config=config)
