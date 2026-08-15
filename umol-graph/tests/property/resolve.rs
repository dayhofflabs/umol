//! Search-independence of aromatic assignment selection.
//!
//! Semantic property: for carriers whose components stay within the
//! assignment bound, `AromaticityResolver::select` returns the outcome of a
//! selection that enumerates every assignment of every component
//! exhaustively — the pruned search never removes a valid assignment.
//!
//! Operational domain: `strategies::select_scenario` — one- and two-ring
//! Hückel skeletons, all ring atoms in the carrier with literal
//! contributions, no stored systems, both failure policies and both
//! tie-breaks on each axis.
//!
//! Validation method: comparison with the definition-level selection below,
//! which enumerates assignments flat with no pruning. Perception
//! (`find_systems`) and the value-key comparison (`compare_by_key`) are
//! shared with production — the property targets the search and selection
//! stages, not the perception rule.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet};

use proptest::prelude::*;
use umol_graph::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityInconsistency, AromaticityPerception,
};
use umol_graph::ops::model::{AromaticityTieBreak, ValenceTieBreak};
use umol_graph::ops::resolve::{AromaticityFailurePolicy, AromaticityResolver, ResolveState};
use umol_graph::ops::valence::compare::compare_by_key;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AtomForm, AtomId, NumForm, RingConfig, RingModel,
    RingSetKind,
};
use umol_utils::solution::Solution;

use crate::strategies::{select_scenario, SelectScenario};

fn contribution(form: &AtomForm) -> Option<u8> {
    match form.constraints.aromatic_valence() {
        Some(AromaticValenceForm::Aromatic(NumForm::Lit(valence))) => u8::try_from(*valence).ok(),
        _ => None,
    }
}

fn compare_restrictions(
    a: &[(AtomId, AtomForm)],
    b: &[(AtomId, AtomForm)],
    tie_break: ValenceTieBreak,
) -> Ordering {
    let key = tie_break.key();
    let b_forms: BTreeMap<AtomId, &AtomForm> = b.iter().map(|(atom, form)| (*atom, form)).collect();
    for (atom, a_form) in a {
        let Some(b_form) = b_forms.get(atom) else {
            continue;
        };
        let ordering = compare_by_key(key, a_form, b_form);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    Ordering::Equal
}

/// Definition-level selection: flat assignment enumeration per component with
/// no pruning and no assignment bound, followed by the documented validity
/// and selection stages. Scoped to the generated domain: no stored systems,
/// no atom-level assertions, every ring atom in the carrier.
fn exhaustive_select(
    scenario: &SelectScenario,
) -> Solution<ResolveState, AromaticityContradiction> {
    let SelectScenario {
        molecule,
        completions,
        model,
        config,
        tie_break,
    } = scenario;
    let perception = AromaticityPerception::new(model);
    let mut completions = completions.clone();
    let mut systems: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
    let mut tie_breaks: Vec<AtomId> = Vec::new();

    let flexible: Vec<(AtomId, Vec<Option<u8>>)> = completions
        .iter()
        .filter_map(|(atom, disjuncts)| {
            let contributions: Vec<Option<u8>> = disjuncts.iter().map(contribution).collect();
            (contributions.iter().collect::<BTreeSet<_>>().len() > 1)
                .then_some((atom, contributions))
        })
        .collect();

    let rings = molecule
        .rings(
            RingModel {
                kind: RingSetKind::Relevant,
                max_ring_size: 22,
            },
            RingConfig::default(),
        )
        .into_ring_set();
    let capable = |atom: AtomId| {
        completions
            .get(atom)
            .is_some_and(|disjuncts| disjuncts.iter().any(|form| contribution(form).is_some()))
    };
    let mut components: Vec<BTreeSet<AtomId>> = Vec::new();
    for ring in rings.iter() {
        if !ring.atoms().iter().all(|&atom| capable(atom)) {
            continue;
        }
        let ring: BTreeSet<AtomId> = ring.atoms().iter().copied().collect();
        let (connected, disjoint): (Vec<_>, Vec<_>) = components
            .into_iter()
            .partition(|component| !component.is_disjoint(&ring));
        let mut merged = ring;
        for component in connected {
            merged.extend(component);
        }
        components = disjoint;
        components.push(merged);
    }
    components.sort();

    type Assignment = (
        Vec<(AtomId, AtomForm)>,
        Vec<(Vec<AtomId>, AromaticSystemForm)>,
    );
    let mut accepted: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
    let mut tie_break_uses: BTreeSet<AtomId> = BTreeSet::new();
    let mut claimed: BTreeSet<AtomId> = BTreeSet::new();
    for component in &components {
        let component_flexible: Vec<(AtomId, Vec<Option<u8>>)> = flexible
            .iter()
            .filter(|(atom, _)| component.contains(atom))
            .cloned()
            .collect();
        let aromatic_only: Vec<AtomId> = component
            .iter()
            .copied()
            .filter(|&atom| {
                completions.get(atom).is_some_and(|disjuncts| {
                    disjuncts.iter().all(|form| {
                        matches!(
                            form.constraints.aromatic_valence(),
                            Some(AromaticValenceForm::Aromatic(_))
                        )
                    })
                })
            })
            .collect();

        let mut assignments: Vec<Assignment> = Vec::new();
        let mut assignment_indices = vec![0usize; component_flexible.len()];
        loop {
            let choice: BTreeMap<AtomId, usize> = component_flexible
                .iter()
                .zip(&assignment_indices)
                .map(|(&(atom, _), &index)| (atom, index))
                .collect();
            let outcome =
                perception
                    .find_systems(molecule, AromaticityConfig::default(), |atom| match choice
                        .get(&atom)
                    {
                        Some(&index) => {
                            contribution(&completions.get(atom).expect("flexible atom")[index])
                        }
                        None => completions
                            .get(atom)
                            .and_then(|disjuncts| contribution(&disjuncts[0])),
                    })
                    .expect("Hückel perception is infallible");
            let found = match outcome {
                Solution::Determined(found) => found,
                Solution::Underdetermined(_) => {
                    return Solution::Underdetermined(ResolveState {
                        completions,
                        systems,
                        tie_breaks,
                    });
                }
                Solution::Contradictory(contradiction) => {
                    return Solution::Contradictory(contradiction);
                }
            };
            let mut partition: Vec<(Vec<AtomId>, AromaticSystemForm)> = found
                .into_iter()
                .filter(|(atoms, _)| atoms.iter().all(|atom| component.contains(atom)))
                .collect();
            partition.sort_by(|(a, _), (b, _)| a.cmp(b));
            let restriction: Vec<(AtomId, AtomForm)> = choice
                .iter()
                .filter(|(atom, _)| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
                .map(|(&atom, &index)| {
                    (
                        atom,
                        completions.get(atom).expect("flexible atom")[index].clone(),
                    )
                })
                .collect();
            let assignment = (restriction, partition);
            if !assignments.contains(&assignment) {
                assignments.push(assignment);
            }

            let mut position = component_flexible.len();
            loop {
                if position == 0 {
                    break;
                }
                position -= 1;
                assignment_indices[position] += 1;
                if assignment_indices[position] < component_flexible[position].1.len() {
                    break;
                }
                assignment_indices[position] = 0;
            }
            if component_flexible.is_empty() || assignment_indices.iter().all(|&index| index == 0) {
                break;
            }
        }

        let mut valid = assignments;
        if config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
            valid.retain(|(_, partition)| {
                aromatic_only
                    .iter()
                    .all(|atom| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
            });
        }
        if valid.is_empty() {
            continue;
        }
        for (_, partition) in &valid {
            for (atoms, _) in partition {
                claimed.extend(atoms.iter().copied());
            }
        }

        let mut structural_decided = false;
        if model.tie_break == AromaticityTieBreak::MaxAtomCount {
            let structure = |partition: &Vec<(Vec<AtomId>, AromaticSystemForm)>| {
                (
                    Reverse(
                        partition
                            .iter()
                            .map(|(atoms, _)| atoms.len())
                            .sum::<usize>(),
                    ),
                    partition
                        .iter()
                        .map(|(atoms, _)| atoms.clone())
                        .collect::<Vec<_>>(),
                )
            };
            let best = valid
                .iter()
                .map(|(_, partition)| structure(partition))
                .min()
                .expect("non-empty survivors");
            structural_decided = valid
                .iter()
                .any(|(_, partition)| structure(partition) != best);
            valid.retain(|(_, partition)| structure(partition) == best);
        }

        let mut restrictions: Vec<&Vec<(AtomId, AtomForm)>> = Vec::new();
        for (restriction, _) in &valid {
            if !restrictions.contains(&restriction) {
                restrictions.push(restriction);
            }
        }
        let (winner_restriction, by_key) = if restrictions.len() == 1 {
            (restrictions[0], false)
        } else {
            if tie_break.key().is_empty() {
                continue;
            }
            let domain: BTreeSet<AtomId> = restrictions[0].iter().map(|(atom, _)| *atom).collect();
            if !restrictions.iter().all(|restriction| {
                restriction
                    .iter()
                    .map(|(atom, _)| *atom)
                    .collect::<BTreeSet<_>>()
                    == domain
            }) {
                continue;
            }
            let best = restrictions
                .iter()
                .copied()
                .max_by(|a, b| compare_restrictions(a, b, *tie_break))
                .expect("non-empty restrictions");
            let unique = restrictions
                .iter()
                .filter(|restriction| compare_restrictions(restriction, best, *tie_break).is_eq())
                .count()
                == 1;
            if !unique {
                continue;
            }
            (best, true)
        };
        let (_, partition) = valid
            .iter()
            .filter(|(restriction, _)| restriction == winner_restriction)
            .min_by(|(_, a), (_, b)| {
                a.iter()
                    .map(|(atoms, _)| atoms)
                    .cmp(b.iter().map(|(atoms, _)| atoms))
            })
            .expect("winner restriction present");
        accepted.extend(partition.iter().cloned());
        if structural_decided {
            for (atoms, _) in partition {
                tie_break_uses.extend(atoms.iter().copied());
            }
        }
        let narrowed: Vec<(AtomId, AtomForm)> = winner_restriction.clone();
        for (atom, form) in &narrowed {
            if by_key && completions.get(*atom).is_some_and(|entry| entry.len() > 1) {
                tie_break_uses.insert(*atom);
            }
            completions.insert(*atom, smallvec::smallvec![form.clone()]);
        }
    }

    if config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
        for (atom, disjuncts) in completions.iter() {
            if claimed.contains(&atom) {
                continue;
            }
            if disjuncts.iter().all(|form| {
                matches!(
                    form.constraints.aromatic_valence(),
                    Some(AromaticValenceForm::Aromatic(_))
                )
            }) {
                return Solution::Contradictory(
                    AromaticityInconsistency::AromaticValenceFailure { atom }.into(),
                );
            }
        }
    }

    systems.extend(accepted);
    tie_breaks.extend(tie_break_uses);
    tie_breaks.sort_unstable();
    tie_breaks.dedup();
    Solution::Determined(ResolveState {
        completions,
        systems,
        tie_breaks,
    })
}

proptest! {
    #[test]
    fn test_aromaticity_resolver_select_search_independence(
        scenario in select_scenario()
    ) {
        let resolver = AromaticityResolver::with_config(&scenario.model, scenario.config);
        let actual = resolver
            .select(
                &scenario.molecule,
                ResolveState {
                    completions: scenario.completions.clone(),
                    ..ResolveState::default()
                },
                scenario.tie_break,
            )
            .expect("Hückel perception is infallible");
        prop_assert_eq!(actual, exhaustive_select(&scenario));
    }
}
