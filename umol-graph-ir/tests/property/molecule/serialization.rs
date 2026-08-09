//! Molecule AST and DSL serialization properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use rstest::rstest;

use crate::strategies::*;

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

    /// Direct `MoleculeAst::ToEdn` / `FromEdn` round-trips are the identity.
    /// Refs render as positional integers (no id keywords); the AST carries
    /// no metadata, so canonical EDN parses back to an equal AST.
    #[test]
    fn test_molecule_ast_to_edn_from_edn_tree_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let edn = ast.to_edn();
        let parsed = MoleculeAst::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(ast, parsed);
    }

    #[test]
    fn test_molecule_ast_to_edn_from_edn_str_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let rendered = ast.to_edn().to_string();
        let parsed = MoleculeAst::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(ast, parsed);
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
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let defaults = MoleculeDefaults::new();
        let rebuilt = MoleculeDsl::from_ir(&ast, &defaults).into_ir(&defaults);
        prop_assert_eq!(rebuilt, ast);
    }

    #[test]
    fn test_molecule_defaults_roundtrip_ground(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let required = MoleculeDefaults::new();
        let ground = MoleculeDefaults::ground();
        let grounded = MoleculeDsl::from_ir(&ast, &required).into_ir(&ground);
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
    let atoms = vec![AtomAst::default(), AtomAst::default()];
    let mut cs = Constraints::new();
    cs.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        constraints: cs,
        ..Default::default()
    });

    let mut metadata = MoleculeMetadata::new();
    metadata
        .set_keyword(Entity::Atom(AtomId(0)), "carbon".to_string())
        .unwrap();

    let dsl = MoleculeDsl::new(ast, metadata).unwrap();
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
