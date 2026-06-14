use proptest::prelude::*;
use crate::strategies::*;

proptest! {
    /// `simplify` is idempotent: `x.simplify().simplify() == x.simplify()`.
    #[test]
    fn test_value_ast_simplify_idempotent(v in any_value_ast_strategy()) {
        let once = v.simplify();
        let twice = once.clone().simplify();
        prop_assert_eq!(once, twice);
    }

    /// `simplify()` is the canonical form: for any generated `ValueAst`,
    /// rendering and parsing yields a value that — once simplified —
    /// equals `simplify()` on the original. The parser produces a partly
    /// canonical form (it folds within `ValueExpr` but doesn't always lift
    /// `ValueExpr(Lit(n))` to `ValueAst::Lit(n)`); simplify completes the
    /// canonicalization on both sides.
    #[test]
    fn test_value_ast_render_parse_equals_simplify(v in any_value_ast_strategy()) {
        let dsl = ValueDsl(v.clone());
        let rendered = dsl.to_string();
        let parsed = parse_value(&rendered).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered:?}"))
        })?;
        prop_assert_eq!(parsed.simplify(), v.simplify());
    }
}
