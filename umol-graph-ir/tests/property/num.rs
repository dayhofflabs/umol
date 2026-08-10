use proptest::prelude::*;
use umol_graph_ir::ir::Normalize;

use crate::strategies::*;

proptest! {
    /// `normalize` is idempotent: normalizing the normal form is a no-op.
    #[test]
    fn test_num_form_normalize_idempotent(v in any_num_form_strategy()) {
        let once = v.normalize();
        let twice = once.clone().and_then(Normalize::normalize);
        prop_assert_eq!(once, twice);
    }

    /// `normalize()` is the normal form: for any generated `NumForm`,
    /// rendering and parsing yields a value that — once normalized — equals
    /// `normalize()` on the original. The parser is faithful (no folding);
    /// `normalize` completes normalization on both sides.
    #[test]
    fn test_num_form_render_parse_equals_normalize(v in any_num_form_strategy()) {
        let dsl = NumDsl(v.clone());
        let rendered = dsl.to_string();
        let parsed = parse_num(&rendered).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered:?}"))
        })?;
        prop_assert_eq!(parsed.normalize(), v.normalize());
    }
}
