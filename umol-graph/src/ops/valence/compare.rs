//! Shared tie-break ordering for valence resolution candidates.
//!
//! [`compare_valence_preference`] returns [`Ordering::Greater`] when `a` is
//! preferred over `b` (max `#h`, max `#n`, min `#u`). Use with
//! [`Iterator::max_by`](Iterator::max_by); for ascending sort, reverse the
//! arguments or use [`Iterator::min_by`](Iterator::min_by) on the reversed cmp.

use std::cmp::Ordering;

use umol_ast::ast::{AsLit, AtomAst};

/// Prefer higher implicit hydrogens, then lone pairs, then fewer unpaired electrons.
pub fn compare_valence_preference(a: &AtomAst, b: &AtomAst) -> Ordering {
    let ha = a
        .implicit_hydrogens
        .as_lit_expect("valence preference requires literal implicit hydrogens");
    let hb = b
        .implicit_hydrogens
        .as_lit_expect("valence preference requires literal implicit hydrogens");
    let na = a
        .lone_pairs
        .as_lit_expect("valence preference requires literal lone pairs");
    let nb = b
        .lone_pairs
        .as_lit_expect("valence preference requires literal lone pairs");
    let ua = a
        .spin
        .unpaired
        .as_lit_expect("valence preference requires literal unpaired electrons");
    let ub = b
        .spin
        .unpaired
        .as_lit_expect("valence preference requires literal unpaired electrons");
    ha.cmp(&hb).then(na.cmp(&nb)).then(ub.cmp(&ua))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_ast::ast::{AtomAst, SpinStateAst, ValueAst};
    use umol_chem::element::Element;

    use super::*;

    fn atom(h: i64, n: i64, u: i64) -> AtomAst {
        let mut a = AtomAst::from_element(Element::C);
        a.implicit_hydrogens = ValueAst::Lit(h);
        a.lone_pairs = ValueAst::Lit(n);
        a.spin = SpinStateAst {
            unpaired: ValueAst::Lit(u),
            multiplicity: ValueAst::Lit(1),
        };
        a
    }

    #[rstest]
    #[case::higher_h(atom(3, 0, 0), atom(1, 1, 0), Ordering::Greater)]
    #[case::higher_n(atom(1, 1, 0), atom(3, 0, 0), Ordering::Less)]
    #[case::lower_u(atom(3, 0, 0), atom(3, 0, 2), Ordering::Greater)]
    #[case::equal(atom(2, 1, 0), atom(2, 1, 0), Ordering::Equal)]
    fn test_compare_valence_preference(
        #[case] a: AtomAst,
        #[case] b: AtomAst,
        #[case] expected: Ordering,
    ) {
        assert_eq!(compare_valence_preference(&a, &b), expected);
        assert_eq!(
            compare_valence_preference(&b, &a),
            expected.reverse(),
            "antisymmetric"
        );
    }

    #[rstest]
    fn test_compare_valence_preference_max_by_picks_higher_h() {
        let candidates = [atom(1, 1, 0), atom(3, 0, 0), atom(2, 0, 0)];
        let best = candidates
            .iter()
            .max_by(|a, b| compare_valence_preference(a, b))
            .unwrap();
        assert_eq!(best.implicit_hydrogens.as_lit(), Some(3));
    }
}
