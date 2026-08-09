//! Shared tie-break ordering for valence resolution candidates.
//!
//! [`compare_valence_preference`] returns [`Ordering::Greater`] when `a` is
//! preferred over `b` (max `#h`, max `#n`, min `#u`). Use with
//! [`Iterator::max_by`](Iterator::max_by); for ascending sort, reverse the
//! arguments or use [`Iterator::min_by`](Iterator::min_by) on the reversed cmp.

use std::cmp::Ordering;

use umol_graph_ir::ir::{AsLit, AtomAst};

/// Prefer higher implicit hydrogens, then lone pairs, then fewer unpaired electrons.
pub fn compare_valence_preference(a: &AtomAst, b: &AtomAst) -> Ordering {
    let ha = a
        .implicit_hydrogens
        .as_lit()
        .expect("valence preference requires literal implicit hydrogens");
    let hb = b
        .implicit_hydrogens
        .as_lit()
        .expect("valence preference requires literal implicit hydrogens");
    let na = a
        .lone_pairs
        .as_lit()
        .expect("valence preference requires literal lone pairs");
    let nb = b
        .lone_pairs
        .as_lit()
        .expect("valence preference requires literal lone pairs");
    let ua = a
        .unpaired_electrons
        .count
        .as_lit()
        .expect("valence preference requires literal unpaired electrons");
    let ub = b
        .unpaired_electrons
        .count
        .as_lit()
        .expect("valence preference requires literal unpaired electrons");
    ha.cmp(&hb).then(na.cmp(&nb)).then(ub.cmp(&ua))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomAst, ElementAst, NumForm, UnpairedElectronsAst};

    use super::*;

    #[rstest]
    #[case::higher_h(
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(1), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Greater,
    )]
    #[case::higher_n(
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(1), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Less,
    )]
    #[case::lower_u(
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((2_u8, 1_u8)), ..Default::default() },
        Ordering::Greater,
    )]
    #[case::equal(
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(2), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        AtomAst { element: ElementAst::Lit(Element::C), implicit_hydrogens: NumForm::Lit(2), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Equal,
    )]
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
    fn test_compare_valence_preference_iterator() {
        let candidates = [
            AtomAst {
                element: ElementAst::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(1),
                lone_pairs: NumForm::Lit(1),
                unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
                ..Default::default()
            },
            AtomAst {
                element: ElementAst::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(3),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
                ..Default::default()
            },
            AtomAst {
                element: ElementAst::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
                ..Default::default()
            },
        ];
        let best = candidates
            .iter()
            .max_by(|a, b| compare_valence_preference(a, b))
            .unwrap();
        assert_eq!(best, &candidates[1]);
    }
}
