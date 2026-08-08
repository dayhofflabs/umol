use proptest::prelude::*;
use umol_graph_ir::ir::Canonicalize;

use crate::strategies::*;

proptest! {
    /// `canonicalize` is idempotent: re-canonicalizing the canonical form is a no-op.
    #[test]
    fn test_value_ast_canonicalize_idempotent(v in any_value_ast_strategy()) {
        let once = v.canonicalize();
        let twice = once.clone().and_then(Canonicalize::canonicalize);
        prop_assert_eq!(once, twice);
    }

    /// `canonicalize()` is the canonical form: for any generated `ValueAst`,
    /// rendering and parsing yields a value that — once canonicalized — equals
    /// `canonicalize()` on the original. The parser is faithful (no folding);
    /// `canonicalize` completes the canonicalization on both sides.
    #[test]
    fn test_value_ast_render_parse_equals_canonicalize(v in any_value_ast_strategy()) {
        let dsl = ValueDsl(v.clone());
        let rendered = dsl.to_string();
        let parsed = parse_value(&rendered).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered:?}"))
        })?;
        prop_assert_eq!(parsed.canonicalize(), v.canonicalize());
    }
}
