//! Shared tie-break comparison for valence resolution candidates.
//!
//! [`compare_by_key`] evaluates a policy's lexicographic key
//! (`ValenceTieBreak::key`); the greatest candidate under the composed
//! ordering is the selected one. Use with
//! [`Iterator::max_by`].

use std::cmp::Ordering;

use umol_graph_ir::ir::{AsLit, AtomFieldKind, AtomForm};

use crate::utils::SortingDirection;

/// Lexicographic comparison of `a` and `b` under `key`: each component
/// compares the field's literal value in its direction — `Ascending` is the
/// field's natural order, `Descending` its reverse — and the first non-equal
/// component decides. An empty key yields `Ordering::Equal`.
///
/// # Panics
///
/// Panics when a compared field is non-literal, or on a kind with no
/// tie-break ordering (`Element`, `IsotopeMass`) — policy keys are crate
/// constants and carry only ordered numeric kinds.
pub fn compare_by_key(
    key: &[(AtomFieldKind, SortingDirection)],
    a: &AtomForm,
    b: &AtomForm,
) -> Ordering {
    key.iter()
        .map(|&(kind, direction)| {
            let ordering = field_rank(kind, a).cmp(&field_rank(kind, b));
            match direction {
                SortingDirection::Ascending => ordering,
                SortingDirection::Descending => ordering.reverse(),
            }
        })
        .find(|ordering| *ordering != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

fn field_rank(kind: AtomFieldKind, atom: &AtomForm) -> i64 {
    match kind {
        AtomFieldKind::Charge => atom
            .charge
            .as_lit()
            .expect("tie-break requires a literal charge"),
        AtomFieldKind::ImplicitHydrogens => atom
            .implicit_hydrogens
            .as_lit()
            .expect("tie-break requires literal implicit hydrogens"),
        AtomFieldKind::LonePairs => atom
            .lone_pairs
            .as_lit()
            .expect("tie-break requires literal lone pairs"),
        AtomFieldKind::UnpairedElectrons => atom
            .unpaired_electrons
            .count
            .as_lit()
            .expect("tie-break requires literal unpaired electrons"),
        AtomFieldKind::Element | AtomFieldKind::IsotopeMass => {
            panic!("no tie-break ordering for {kind:?}")
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomForm, ElementForm, NumForm, UnpairedElectronsForm};

    use super::*;

    const MOST_SATURATED: [(AtomFieldKind, SortingDirection); 3] = [
        (
            AtomFieldKind::ImplicitHydrogens,
            SortingDirection::Ascending,
        ),
        (AtomFieldKind::LonePairs, SortingDirection::Ascending),
        (
            AtomFieldKind::UnpairedElectrons,
            SortingDirection::Descending,
        ),
    ];

    #[rstest]
    #[case::higher_h(
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(1), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Greater,
    )]
    #[case::higher_n(
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(1), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Less,
    )]
    #[case::lower_u(
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((2_u8, 1_u8)), ..Default::default() },
        Ordering::Greater,
    )]
    #[case::equal(
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(2), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(2), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        Ordering::Equal,
    )]
    fn test_compare_by_key(#[case] a: AtomForm, #[case] b: AtomForm, #[case] expected: Ordering) {
        assert_eq!(compare_by_key(&MOST_SATURATED, &a, &b), expected);
        assert_eq!(
            compare_by_key(&MOST_SATURATED, &b, &a),
            expected.reverse(),
            "antisymmetric"
        );
    }

    #[rstest]
    #[case::differing_forms(
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
        AtomForm { element: ElementForm::Lit(Element::C), implicit_hydrogens: NumForm::Lit(1), lone_pairs: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() },
    )]
    fn test_compare_by_key_empty(#[case] a: AtomForm, #[case] b: AtomForm) {
        assert_eq!(compare_by_key(&[], &a, &b), Ordering::Equal);
    }

    #[rstest]
    fn test_compare_by_key_iterator() {
        let candidates = [
            AtomForm {
                element: ElementForm::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(1),
                lone_pairs: NumForm::Lit(1),
                unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                ..Default::default()
            },
            AtomForm {
                element: ElementForm::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(3),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                ..Default::default()
            },
            AtomForm {
                element: ElementForm::Lit(Element::C),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                ..Default::default()
            },
        ];
        let best = candidates
            .iter()
            .max_by(|a, b| compare_by_key(&MOST_SATURATED, a, b))
            .unwrap();
        assert_eq!(best, &candidates[1]);
    }
}
