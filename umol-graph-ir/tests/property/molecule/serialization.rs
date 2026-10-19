//! Molecule graph IR and DSL serialization properties.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use rstest::rstest;
use umol_edn::{DeError, Edn, EdnError, ParseError as EdnParseError};

use crate::strategies::*;

fn hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    #[test]
    fn test_molecule_dsl_to_edn_from_edn_tree_roundtrip(dsl in molecule_dsl_strategy()) {
        let edn = dsl.to_edn();
        let parsed = MoleculeDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_to_edn_from_edn_str_roundtrip(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_streaming_matches_tree(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let via_stream = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming: {e}")))?;
        let tree = read_string(&rendered)
            .map_err(|e| TestCaseError::fail(format!("edn parse: {e}")))?;
        let via_tree = MoleculeDsl::from_edn(&tree)
            .map_err(|e| TestCaseError::fail(format!("tree: {e}")))?;
        prop_assert_eq!(via_stream, via_tree);
    }

    #[test]
    fn test_molecule_dsl_parser_parity(input in any::<String>()) {
        let via_stream = MoleculeDsl::from_edn_str(&input).ok();
        let via_tree = read_string(&input)
            .ok()
            .and_then(|edn| MoleculeDsl::from_edn(&edn).ok());
        prop_assert_eq!(via_stream, via_tree);
    }

    /// Adding an unknown field to a valid structural reference is rejected by both parsers.
    #[test]
    fn test_molecule_dsl_from_edn_unknown_key(
        suffix in "[a-z]{1,12}",
        value in any::<i64>(),
        first in any::<bool>(),
        namespaced in any::<bool>(),
    ) {
        let key = if namespaced { format!("{suffix}/atoms") } else { format!("unknown-{suffix}") };
        let entry = format!(":{key} {value}");
        let fields = if first {
            format!("{entry} :atoms [0 1]")
        } else {
            format!(":atoms [0 1] {entry}")
        };
        let source = format!(
            r#"{{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{{:bond [{{{fields}}} {{:aromatic true}}]}}]}}"#,
        );
        let expected = DeError::UnknownField {
            key,
            path: vec!["bond structural ref".into()],
        };
        let tree = read_string(&source).unwrap();
        prop_assert_eq!(MoleculeDsl::from_edn(&tree), Err(expected.clone()));
        prop_assert_eq!(MoleculeDsl::from_edn_str(&source), Err(expected.into()));
    }

    /// EDN permits non-keyword map keys; the molecule DSL does not.
    #[test]
    fn test_molecule_dsl_from_edn_non_keyword_key(
        values in prop::collection::vec(any::<i64>(), 0..5),
        first in any::<bool>(),
    ) {
        let key = Edn::Vector(values.into_iter().map(Edn::Int).collect::<Vec<_>>().into());
        let fields = if first {
            format!("{key} 0 :atoms [0 1]")
        } else {
            format!(":atoms [0 1] {key} 0")
        };
        let source = format!(
            r#"{{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{{:bond [{{{fields}}} {{:aromatic true}}]}}]}}"#,
        );
        let tree = read_string(&source).unwrap();
        prop_assert_eq!(MoleculeDsl::from_edn(&tree), Err(DeError::TypeMismatch {
            expected: "keyword key",
            got: "vector",
            path: vec!["bond structural ref".into()],
        }));
        prop_assert!(matches!(
            MoleculeDsl::from_edn_str(&source),
            Err(EdnError::Parse(EdnParseError::UnexpectedToken { found: '[', .. }))
        ), "non-keyword map key was not rejected at its opening bracket");
    }

    /// Qualifying a grammar keyword does not retain its unqualified meaning.
    #[test]
    fn test_molecule_dsl_from_edn_namespaced_keyword(
        namespace in "[a-z]{1,12}",
        (source, keyword) in prop::sample::select(vec![
            (r#"{:atoms ["C" "C"] :bonds [[0 1 :single]]}"#, "single"),
            (r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs :single}]}"#, "single"),
            (r#"{:atoms ["C"] :constraints [{:atom [0 {:valence :undetermined}]}]}"#, "undetermined"),
            (r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:aromatic :undetermined}]}]}"#, "undetermined"),
            (r#"{:atoms ["C"] :constraints [{:atom [0 {:aromatic-valence :not-aromatic}]}]}"#, "not-aromatic"),
            (r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs :ccw}]}"#, "ccw"),
        ]),
    ) {
        let valid = MoleculeDsl::from_edn_str(source).unwrap();
        prop_assert_eq!(MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap(), valid);
        let source = source.replace(&format!(":{keyword}"), &format!(":{namespace}/{keyword}"));
        prop_assert_eq!(MoleculeDsl::from_edn(&read_string(&source).unwrap()).ok(), None);
        prop_assert_eq!(MoleculeDsl::from_edn_str(&source).ok(), None);
    }

    /// Repeating a map field is invalid even when the repeated values are identical.
    #[test]
    fn test_molecule_dsl_from_edn_str_duplicate_key(
        dsl in molecule_dsl_strategy(),
        index in 0..10usize,
        first in any::<bool>(),
    ) {
        let edn = dsl.to_edn();
        let Edn::Map(map) = &edn else { unreachable!() };
        let (key, value) = map.iter().nth(index % map.len()).unwrap();
        let rendered = edn.to_string();
        let original_fields = &rendered[1..rendered.len() - 1];
        let source = if first {
            format!("{{{key} {value} {original_fields}}}")
        } else {
            format!("{{{original_fields} {key} {value}}}")
        };
        let expected = read_string(&source).unwrap_err();
        prop_assert!(matches!(expected, EdnError::Parse(EdnParseError::DuplicateKey { .. })),
            "repeated field was not rejected as a duplicate key");
        prop_assert_eq!(MoleculeDsl::from_edn_str(&source), Err(expected));
    }

    /// Direct `Molecule::ToEdn` / `FromEdn` round-trips are the identity.
    /// Refs render as positional integers (no id keywords); the molecule carries
    /// no metadata, so canonical EDN parses back to an equal molecule.
    #[test]
    fn test_molecule_to_edn_from_edn_tree_roundtrip(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let edn = molecule.to_edn();
        let parsed = Molecule::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(molecule, parsed);
    }

    #[test]
    fn test_molecule_to_edn_from_edn_str_roundtrip(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let rendered = molecule.to_edn().to_string();
        let parsed = Molecule::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(molecule, parsed);
    }

    #[test]
    fn test_molecule_structural_equality_implies_equal_hash(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let parsed = Molecule::from_edn(&molecule.to_edn())
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(&molecule, &parsed);
        prop_assert_eq!(hash(&molecule), hash(&parsed));
    }

    #[test]
    fn test_molecule_dsl_double_render_is_stable(dsl in molecule_dsl_strategy()) {
        let s1 = dsl.to_edn().to_string();
        let d1 = MoleculeDsl::from_edn_str(&s1)
            .map_err(|e| TestCaseError::fail(format!("first parse: {e}")))?;
        let s2 = d1.to_edn().to_string();
        prop_assert_eq!(s1, s2);
    }

    #[test]
    fn test_molecule_defaults_roundtrip(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let defaults = MoleculeDefaults::new();
        let rebuilt = MoleculeDsl::from_ir(&molecule, &defaults).into_ir(&defaults);
        prop_assert_eq!(rebuilt, molecule);
    }

    #[test]
    fn test_molecule_defaults_roundtrip_ground(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let required = MoleculeDefaults::new();
        let ground = MoleculeDefaults::concrete();
        let grounded = MoleculeDsl::from_ir(&molecule, &required).into_ir(&ground);
        let rebuilt = MoleculeDsl::from_ir(&grounded, &ground).into_ir(&ground);
        prop_assert_eq!(rebuilt, grounded);
    }
}

/// When `MoleculeMetadata` records an id for an entity, refs in molecule constraints
/// render as the keyword `:id` rather than the positional integer. Rendered
/// EDN must carry the keyword form, never the integer index, and must
/// roundtrip back through both the tree and streaming parsers.
#[rstest]
fn test_constraint_ref_uses_keyword_when_metadata_binding_present() {
    let atoms = vec![AtomForm::default(), AtomForm::default()];
    let mut cs = Constraints::new();
    cs.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms,
        constraints: cs,
        ..Default::default()
    });

    let mut metadata = MoleculeMetadata::new();
    metadata
        .set_keyword(Entity::Atom(AtomId(0)), "carbon".to_string())
        .unwrap();

    let dsl = MoleculeDsl::new(molecule, metadata).unwrap();
    let rendered = dsl.to_edn().to_string();

    assert!(
        rendered.contains(":carbon"),
        "expected :carbon in rendered output: {rendered}",
    );
    assert!(
        !rendered.contains("[0 {:valence"),
        "rendered output must not use positional ref when id is present: {rendered}",
    );

    let via_tree = MoleculeDsl::from_edn(&dsl.to_edn()).expect("tree parse");
    assert_eq!(dsl, via_tree, "tree roundtrip");
    let via_stream = MoleculeDsl::from_edn_str(&rendered).expect("streaming parse");
    assert_eq!(dsl, via_stream, "streaming roundtrip");
}
