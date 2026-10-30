//! Exhaust every binary assignment, filter by stated domains/admissibility, then sort.
//! This reference does not use a transformer, production traversal, or stereo validator.

use std::collections::BTreeSet;

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_graph_ir::ir::{StereoCoset, StereoTerm};

#[derive(Debug, PartialEq, Eq)]
pub(in crate::ops::transform) enum ReferenceError {
    UnavailableSite(usize),
    InvalidOrder {
        missing: Vec<usize>,
        duplicates: Vec<usize>,
    },
    UnsupportedExpression(usize),
    TooManySites(usize),
}

pub(in crate::ops::transform) fn absolute_assignments(
    source: &[StereoCoset],
    selected: &BTreeSet<usize>,
    order: &[usize],
    admissible: impl Fn(&[StereoCoset]) -> bool,
) -> Result<Vec<Vec<StereoCoset>>, ReferenceError> {
    if let Some(&site) = selected.iter().find(|&&site| site >= source.len()) {
        return Err(ReferenceError::UnavailableSite(site));
    }
    if selected.len() > 8 {
        return Err(ReferenceError::TooManySites(selected.len()));
    }
    let mut seen = BTreeSet::new();
    let mut duplicates = BTreeSet::new();
    let mut processing_order = Vec::new();
    for &site in order.iter().filter(|site| selected.contains(site)) {
        if seen.insert(site) {
            processing_order.push(site);
        } else {
            duplicates.insert(site);
        }
    }
    let missing: Vec<_> = selected.difference(&seen).copied().collect();
    if !missing.is_empty() || !duplicates.is_empty() {
        return Err(ReferenceError::InvalidOrder {
            missing,
            duplicates: duplicates.into_iter().collect(),
        });
    }
    for &site in selected {
        if matches!(&source[site], StereoCoset::Term(_)) {
            return Err(ReferenceError::UnsupportedExpression(site));
        }
    }

    let mut results = Vec::new();
    for bits in 0..(1usize << selected.len()) {
        let values: Vec<_> = (0..selected.len())
            .map(|bit| ((bits >> bit) & 1) as u32)
            .collect();
        let allowed = selected
            .iter()
            .zip(&values)
            .all(|(&site, value)| match &source[site] {
                StereoCoset::Undetermined => true,
                StereoCoset::Lit(fixed) => fixed == value,
                StereoCoset::LitSet(domain) => domain.contains(value),
                StereoCoset::Term(_) => unreachable!("selected expressions were rejected"),
            });
        if !allowed {
            continue;
        }
        let mut candidate = source.to_vec();
        for (&site, &value) in selected.iter().zip(&values) {
            candidate[site] = StereoCoset::Lit(value);
        }
        if admissible(&candidate) {
            let key: Vec<_> = processing_order
                .iter()
                .map(|site| {
                    let bit = selected
                        .iter()
                        .position(|selected_site| selected_site == site)
                        .expect("processing sites belong to the selection");
                    values[bit]
                })
                .collect();
            results.push((key, candidate));
        }
    }
    results.sort_by(|(left, _), (right, _)| left.cmp(right));
    Ok(results
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect())
}

#[rstest]
#[case::binary_product(
    vec![StereoCoset::Undetermined, StereoCoset::Undetermined],
    BTreeSet::from([0, 1]), vec![0, 1], |_: &[StereoCoset]| true,
    vec![
        vec![StereoCoset::Lit(0), StereoCoset::Lit(0)],
        vec![StereoCoset::Lit(0), StereoCoset::Lit(1)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(0)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(1)],
    ],
)]
#[case::reverse_order(
    vec![StereoCoset::Undetermined, StereoCoset::Undetermined],
    BTreeSet::from([0, 1]), vec![1, 0], |_: &[StereoCoset]| true,
    vec![
        vec![StereoCoset::Lit(0), StereoCoset::Lit(0)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(0)],
        vec![StereoCoset::Lit(0), StereoCoset::Lit(1)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(1)],
    ],
)]
#[case::fixed_and_finite(
    vec![StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]), StereoCoset::lit_set([0])],
    BTreeSet::from([0, 1, 2]), vec![0, 1, 2], |_: &[StereoCoset]| true,
    vec![
        vec![StereoCoset::Lit(1), StereoCoset::Lit(0), StereoCoset::Lit(0)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(1), StereoCoset::Lit(0)],
    ],
)]
#[case::unselected_open(
    vec![StereoCoset::Undetermined, StereoCoset::Undetermined],
    BTreeSet::from([1]), vec![0, 1, 0], |_: &[StereoCoset]| true,
    vec![
        vec![StereoCoset::Undetermined, StereoCoset::Lit(0)],
        vec![StereoCoset::Undetermined, StereoCoset::Lit(1)],
    ],
)]
#[case::unselected_symbolic(
    vec![StereoCoset::term(StereoTerm::var("retained")), StereoCoset::Undetermined],
    BTreeSet::from([1]), vec![1], |_: &[StereoCoset]| true,
    vec![
        vec![StereoCoset::term(StereoTerm::var("retained")), StereoCoset::Lit(0)],
        vec![StereoCoset::term(StereoTerm::var("retained")), StereoCoset::Lit(1)],
    ],
)]
#[case::asserted_value(
    vec![StereoCoset::Undetermined], BTreeSet::from([0]), vec![0],
    |values: &[StereoCoset]| values[0] == StereoCoset::Lit(1),
    vec![vec![StereoCoset::Lit(1)]],
)]
#[case::joint_admissibility(
    vec![StereoCoset::Undetermined, StereoCoset::Undetermined],
    BTreeSet::from([0, 1]), vec![0, 1], |values: &[StereoCoset]| values[0] != values[1],
    vec![
        vec![StereoCoset::Lit(0), StereoCoset::Lit(1)],
        vec![StereoCoset::Lit(1), StereoCoset::Lit(0)],
    ],
)]
#[case::inconsistent_fixed(
    vec![StereoCoset::Lit(0)], BTreeSet::from([0]), vec![0],
    |values: &[StereoCoset]| values[0] == StereoCoset::Lit(1), vec![],
)]
#[case::empty_domain(
    vec![StereoCoset::lit_set([])], BTreeSet::from([0]), vec![0], |_: &[StereoCoset]| true, vec![],
)]
#[case::empty_selection(
    vec![StereoCoset::Undetermined], BTreeSet::new(), vec![], |_: &[StereoCoset]| true,
    vec![vec![StereoCoset::Undetermined]],
)]
#[case::empty_input(vec![], BTreeSet::new(), vec![], |_: &[StereoCoset]| true, vec![vec![]])]
#[case::empty_rejected(vec![], BTreeSet::new(), vec![], |_: &[StereoCoset]| false, vec![])]
fn test_absolute_assignments(
    #[case] source: Vec<StereoCoset>,
    #[case] selected: BTreeSet<usize>,
    #[case] order: Vec<usize>,
    #[case] admissible: fn(&[StereoCoset]) -> bool,
    #[case] expected: Vec<Vec<StereoCoset>>,
) {
    assert_eq!(
        absolute_assignments(&source, &selected, &order, admissible),
        Ok(expected)
    );
}

#[rstest]
fn test_absolute_assignments_binary_space(#[values(0usize, 1, 2, 3, 4, 5, 6, 7, 8)] sites: usize) {
    let expected = if sites == 0 {
        vec![vec![]]
    } else {
        (0..(1usize << sites))
            .map(|value| {
                format!("{value:0width$b}", width = sites)
                    .bytes()
                    .map(|digit| StereoCoset::Lit(u32::from(digit - b'0')))
                    .collect()
            })
            .collect()
    };
    assert_eq!(
        absolute_assignments(
            &vec![StereoCoset::Undetermined; sites],
            &(0..sites).collect(),
            &(0..sites).collect::<Vec<_>>(),
            |_: &[StereoCoset]| true,
        ),
        Ok(expected),
    );
}

#[rstest]
#[case::unavailable(vec![StereoCoset::Undetermined], BTreeSet::from([1]), vec![1],
    ReferenceError::UnavailableSite(1))]
#[case::missing(vec![StereoCoset::Undetermined], BTreeSet::from([0]), vec![],
    ReferenceError::InvalidOrder { missing: vec![0], duplicates: vec![] })]
#[case::duplicate(vec![StereoCoset::Undetermined], BTreeSet::from([0]), vec![0, 0],
    ReferenceError::InvalidOrder { missing: vec![], duplicates: vec![0] })]
#[case::symbolic(vec![StereoCoset::term(StereoTerm::var("configuration"))],
    BTreeSet::from([0]), vec![0], ReferenceError::UnsupportedExpression(0))]
#[case::bounded(vec![StereoCoset::Undetermined; 9], (0..9).collect(), (0..9).collect(),
    ReferenceError::TooManySites(9))]
fn test_absolute_assignments_error(
    #[case] source: Vec<StereoCoset>,
    #[case] selected: BTreeSet<usize>,
    #[case] order: Vec<usize>,
    #[case] expected: ReferenceError,
) {
    assert_eq!(
        absolute_assignments(&source, &selected, &order, |_: &[StereoCoset]| true),
        Err(expected)
    );
}
