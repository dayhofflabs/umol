//! Compile-and-run check of the whitepaper's Rust-primer appendix listings.
use std::borrow::Cow;
use std::error::Error;

use umol_edn::FromEdn;
use umol_graph::fingerprint::EcfpFeaturizer;
use umol_graph::ingest;
use umol_graph::ingest::ingest_smiles_with;
use umol_graph::ops::model::{
    ChemistryModel, ValenceCandidateSource, ValenceModel, ValenceTieBreak,
};
use umol_graph::ops::resolve::ResolveConfig;
use umol_graph::ops::valence::AtomTypeRegistry;
use umol_graph_core::algorithms::common_subgraph::CommonSubgraphEnumerationAlgorithm;
use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{
    FromIr, IntoIr, Molecule, React, SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use umol_io::smiles::SmilesIoConfig;

#[test]
fn whitepaper_rust_primer() -> Result<(), Box<dyn Error>> {
    // Listing: reading SMILES and lowering to the notation.
    let mol = ingest::ingest_smiles("CCO")?;
    let text = mol.to_string();
    let compact = MoleculeDsl::from_ir(&mol, &MoleculeDefaults::concrete()).to_string();
    assert!(!text.is_empty());
    assert!(compact.len() <= text.len());

    // Listing: notation parse and the metadata round trip.
    let diborane: Molecule = r#"
  {:atoms ["B" "H" "B" "H" "H" "H" "H" "H"]
   :bonds [[0 4 "1"] [0 5 "1"] [2 6 "1"] [2 7 "1"]]
   :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}
                       {:atoms [0 3 2] :attrs "[0,1,1]"}]}"#
        .parse()?;
    assert_eq!(diborane.multicenter_bonds().count(), 2);
    let (mol_back, metadata) = MoleculeDsl::from_edn_str(&text)?.into_parts();
    let round_tripped = MoleculeDsl::new(mol_back, metadata)?.to_string();
    assert_eq!(round_tripped, text);

    // Listing: the chemistry-model example through the model types.
    let registry =
        AtomTypeRegistry::from_toml_str("[P]\n0 = [\"P #n #v3\"]\n\n[F]\n0 = [\"F #n3 #v\"]")?;
    let strict = ChemistryModel {
        valence: ValenceModel {
            candidates: ValenceCandidateSource::AtomTyping {
                registry: Cow::Owned(registry),
            },
            tie_break: ValenceTieBreak::Strict,
        },
        ..ChemistryModel::default()
    };
    let io = SmilesIoConfig::opensmiles();
    let accepted = ingest_smiles_with("FP(F)F", &io, &strict, &ResolveConfig::default())?;
    assert_eq!(accepted.atoms().count(), 4);
    let refused = ingest_smiles_with("FP(F)(F)(F)F", &io, &strict, &ResolveConfig::default());
    assert!(refused
        .unwrap_err()
        .to_string()
        .contains("no atom-typing match"));

    // Listing: substructure search and reaction application through the
    // gathered configuration.
    let pattern =
        MoleculeDsl::from_edn_str(r##"{:atoms ["N#a2"]}"##)?.into_ir(&MoleculeDefaults::default());
    let rule = ingest::ingest_reaction_smiles(
        "[CH3:1][C:2](=[O:3])[CH2:4][CH2:5][C:6](=[O:7])[CH3:8]\
         >>[CH3:1][c:2]1[cH:4][cH:5][c:6]([CH3:8])[o:3]1.[OH2:7]",
    )?;
    let host = ingest::ingest_smiles("O=C([C@H](Cc1c[nH]cn1)[NH3+])[O-]")?;
    let diketone = ingest::ingest_smiles("CC(=O)CCC(=O)C")?;
    let config = SubstructureMatchConfig {
        match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    };
    let matches = pattern.substructure_matches(&host, config)?;
    assert_eq!(matches.len(), 1);
    let derivations = rule.apply(&diketone, config)?;
    assert_eq!(derivations.count(), 2);

    // Listing: fingerprint featurization, folding, similarity.
    let ecfp = EcfpFeaturizer::new(2);
    let a = ecfp.featurize(&ingest::ingest_smiles("CCO")?)?.fold(2048)?;
    let b = ecfp.featurize(&ingest::ingest_smiles("COC")?)?.fold(2048)?;
    let similarity = a.tanimoto(&b)?;
    assert!((similarity - 0.1111).abs() < 5e-5);

    // Listing: combining and splitting.
    let ammonia = ingest::ingest_smiles("N")?;
    let chloride = ingest::ingest_smiles("[Cl-]")?;
    let (combined, correspondence) = ammonia.combine(&chloride);
    let components = combined.split();
    assert_eq!(correspondence.atoms().matched_pairs().len(), 1);
    assert_eq!(components.len(), 2);

    // Listing: react on a molecule and on a slice of reactants.
    for products in diketone.react(&rule, config)? {
        let products = products?;
        assert_eq!(products.len(), 2);
    }
    let with_reagent = [diketone, ammonia].react(&rule, config)?;
    assert_eq!(with_reagent.count(), 2);

    // Listing: reaction composition through the named algorithm.
    let first = rule.clone();
    let second = rule;
    let composites = first.compose(
        &second,
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
    );
    assert!(!composites.is_empty());

    Ok(())
}
