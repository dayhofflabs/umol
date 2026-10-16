//! Reaction-span properties over generated reactions, explicitly reindexed molecule pairs, and
//! independently generated span entries. The parent module retains the bridge and serialization
//! properties; its children exercise dense union-frame remapping and aggregate canonicalization.

use std::fmt::Debug;
use std::iter;

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_core::{Correspondence, EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomDelta, AtomFieldChange, AtomForm, AtomId, BondId, ConstraintSpan,
    DativeBondId, Delta, EntitySpan, Molecule, MoleculeCorrespondence, MoleculeEntries,
    MoleculeRemapping, MulticenterBondId, NoncovalentBondId, Normalize, Reaction, ReactionSpan,
    ReactionSpanEntries, Reframe, StereoAtomId, StereoBondId, StereoLigand,
};

use crate::strategies::{
    constraint_strategy, element_strategy, materializable_reaction_strategy,
    molecule_entries_strategy, molecule_entries_structurally_unambiguous_strategy,
    overlay_reaction_strategy, reaction_strategy, with_normalized_constraints, ConstraintCounts,
};

#[path = "span/canonicalize.rs"]
mod canonicalize;
#[path = "span/reframe.rs"]
mod reframe;
#[path = "span/remap.rs"]
mod remap;

#[derive(Clone, Copy)]
enum SpanPresence {
    Both,
    Left,
    Right,
}

fn intersection_presence(
    presences: impl IntoIterator<Item = SpanPresence>,
) -> Option<SpanPresence> {
    let mut left = true;
    let mut right = true;
    for presence in presences {
        left &= !matches!(presence, SpanPresence::Right);
        right &= !matches!(presence, SpanPresence::Left);
    }
    match (left, right) {
        (true, true) => Some(SpanPresence::Both),
        (true, false) => Some(SpanPresence::Left),
        (false, true) => Some(SpanPresence::Right),
        (false, false) => None,
    }
}

fn entity_span<T>(value: T, presence: SpanPresence) -> EntitySpan<T> {
    match presence {
        SpanPresence::Both => EntitySpan::Unchanged(value),
        SpanPresence::Left => EntitySpan::Removed(value),
        SpanPresence::Right => EntitySpan::Added(value),
    }
}

fn lhs_anchored<T>(entries: impl IntoIterator<Item = (T, SpanPresence)>) -> Vec<(T, SpanPresence)> {
    let mut entries: Vec<_> = entries.into_iter().collect();
    entries.sort_by_key(|(_, presence)| matches!(presence, SpanPresence::Right));
    entries
}

#[derive(Debug)]
struct ReactionSides {
    lhs: Molecule,
    rhs: Molecule,
    atom_correspondence: Correspondence<AtomId>,
    projected_rhs_atoms: Correspondence<AtomId>,
}

/// A molecule pair related by a crossing, non-total atom correspondence. Three isolated atoms are
/// appended before reversing the rhs atom frame: the final lhs atom and first rhs atom remain
/// unmatched, while every generated bond and overlay remains in the matched substructure.
fn crossing_reaction_sides_strategy() -> impl Strategy<Value = ReactionSides> {
    molecule_entries_structurally_unambiguous_strategy().prop_map(|mut lhs_entries| {
        lhs_entries.atoms.extend([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ]);
        let atom_count = lhs_entries.atoms.len();
        let reverse_atom = |id: AtomId| AtomId((atom_count - 1 - id.0 as usize) as u32);
        let lhs = Molecule::from_entries(lhs_entries.clone());

        let MoleculeEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        } = lhs_entries;
        debug_assert!(constraints.is_empty());
        let rhs = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.into_iter().rev().collect(),
            bonds: bonds
                .into_iter()
                .map(|(first, second, attributes)| {
                    (reverse_atom(first), reverse_atom(second), attributes)
                })
                .collect(),
            dative: dative
                .into_iter()
                .map(|(donors, acceptor, attributes)| {
                    (
                        donors.into_iter().map(reverse_atom).collect(),
                        reverse_atom(acceptor),
                        attributes,
                    )
                })
                .collect(),
            aromatic: aromatic
                .into_iter()
                .map(|(atoms, attributes)| {
                    (atoms.into_iter().map(reverse_atom).collect(), attributes)
                })
                .collect(),
            multicenter: multicenter
                .into_iter()
                .map(|(atoms, attributes)| {
                    (atoms.into_iter().map(reverse_atom).collect(), attributes)
                })
                .collect(),
            noncovalent: noncovalent
                .into_iter()
                .map(|(atoms, attributes)| (atoms.map(reverse_atom), attributes))
                .collect(),
            stereo_atoms: stereo_atoms
                .into_iter()
                .map(|(site, ligands, attributes)| {
                    (
                        reverse_atom(site),
                        ligands
                            .into_iter()
                            .map(|ligand| {
                                StereoLigand::new(reverse_atom(ligand.atom_id), ligand.kind)
                            })
                            .collect(),
                        attributes,
                    )
                })
                .collect(),
            stereo_bonds: stereo_bonds
                .into_iter()
                .map(|(site, ligands, attributes)| {
                    (
                        site,
                        ligands
                            .into_iter()
                            .map(|ligand| {
                                StereoLigand::new(reverse_atom(ligand.atom_id), ligand.kind)
                            })
                            .collect(),
                        attributes,
                    )
                })
                .collect(),
            ..Default::default()
        });

        let atom_correspondence = Correspondence::new(
            (0..atom_count - 1)
                .map(|index| (AtomId(index as u32), reverse_atom(AtomId(index as u32))))
                .collect(),
            atom_count,
            atom_count,
        )
        .expect("generated pairs form a crossing partial bijection");
        let projected_rhs_atoms = Correspondence::from_images(
            &(0..atom_count).rev().map(AtomId::from).collect::<Vec<_>>(),
            atom_count,
        );
        ReactionSides {
            lhs,
            rhs,
            atom_correspondence,
            projected_rhs_atoms,
        }
    })
}

/// Structurally valid, lhs-anchored span entries derived from a generated union molecule. Every
/// entity present on the lhs precedes rhs-only entities of its kind, matching the normalized union
/// frame produced by superimposition.
fn reaction_span_entries_strategy() -> impl Strategy<Value = ReactionSpanEntries> {
    molecule_entries_strategy()
        .prop_flat_map(|entries| {
            let atom_count = entries.atoms.len();
            (
                Just(entries),
                0usize..=atom_count,
                prop::collection::vec(any::<bool>(), atom_count),
            )
        })
        .prop_map(|(entries, lhs_atom_count, shared_atoms)| {
            let atom_presence: Vec<_> = shared_atoms
                .into_iter()
                .enumerate()
                .map(|(index, shared)| {
                    if index >= lhs_atom_count {
                        SpanPresence::Right
                    } else if shared {
                        SpanPresence::Both
                    } else {
                        SpanPresence::Left
                    }
                })
                .collect();
            let presence_of_atom = |id: AtomId| atom_presence[id.0 as usize];

            let bond_presence: Vec<_> = entries
                .bonds
                .iter()
                .map(|(first, second, _)| {
                    intersection_presence([presence_of_atom(*first), presence_of_atom(*second)])
                })
                .collect();
            let bonds = lhs_anchored(entries.bonds.into_iter().enumerate().filter_map(
                |(old_id, (first, second, attributes))| {
                    bond_presence[old_id]
                        .map(|presence| ((old_id, first, second, attributes), presence))
                },
            ));
            let mut bond_ids = vec![None; bond_presence.len()];
            for (new_id, ((old_id, ..), _)) in bonds.iter().enumerate() {
                bond_ids[*old_id] = Some(BondId(new_id as u32));
            }

            let dative = lhs_anchored(entries.dative.into_iter().filter_map(
                |(donors, acceptor, attributes)| {
                    let presence = intersection_presence(
                        donors
                            .iter()
                            .copied()
                            .chain(iter::once(acceptor))
                            .map(presence_of_atom),
                    )?;
                    Some(((donors, acceptor, attributes), presence))
                },
            ));
            let aromatic = lhs_anchored(entries.aromatic.into_iter().filter_map(
                |(atoms, attributes)| {
                    let presence =
                        intersection_presence(atoms.iter().copied().map(presence_of_atom))?;
                    Some(((atoms, attributes), presence))
                },
            ));
            let multicenter = lhs_anchored(entries.multicenter.into_iter().filter_map(
                |(atoms, attributes)| {
                    let presence =
                        intersection_presence(atoms.iter().copied().map(presence_of_atom))?;
                    Some(((atoms, attributes), presence))
                },
            ));
            let noncovalent = lhs_anchored(entries.noncovalent.into_iter().filter_map(
                |(atoms, attributes)| {
                    let presence = intersection_presence(atoms.into_iter().map(presence_of_atom))?;
                    Some(((atoms, attributes), presence))
                },
            ));
            let stereo_atoms = lhs_anchored(entries.stereo_atoms.into_iter().filter_map(
                |(site, ligands, attributes)| {
                    let presence = intersection_presence(
                        iter::once(site)
                            .chain(ligands.iter().map(|ligand| ligand.atom_id))
                            .map(presence_of_atom),
                    )?;
                    Some(((site, ligands, attributes), presence))
                },
            ));
            let stereo_bonds = lhs_anchored(entries.stereo_bonds.into_iter().filter_map(
                |(site, ligands, attributes)| {
                    let site_presence = bond_presence.get(site.0 as usize).copied().flatten()?;
                    let site = bond_ids.get(site.0 as usize).copied().flatten()?;
                    let presence = intersection_presence(
                        iter::once(site_presence).chain(
                            ligands
                                .iter()
                                .map(|ligand| presence_of_atom(ligand.atom_id)),
                        ),
                    )?;
                    Some(((site, ligands, attributes), presence))
                },
            ));

            ReactionSpanEntries {
                atoms: entries
                    .atoms
                    .into_iter()
                    .zip(atom_presence)
                    .map(|(attributes, presence)| entity_span(attributes, presence))
                    .collect(),
                bonds: bonds
                    .into_iter()
                    .map(|((_, first, second, attributes), presence)| {
                        (first, second, entity_span(attributes, presence))
                    })
                    .collect(),
                dative: dative
                    .into_iter()
                    .map(|((donors, acceptor, attributes), presence)| {
                        (donors, acceptor, entity_span(attributes, presence))
                    })
                    .collect(),
                aromatic: aromatic
                    .into_iter()
                    .map(|((atoms, attributes), presence)| {
                        (atoms, entity_span(attributes, presence))
                    })
                    .collect(),
                multicenter: multicenter
                    .into_iter()
                    .map(|((atoms, attributes), presence)| {
                        (atoms, entity_span(attributes, presence))
                    })
                    .collect(),
                noncovalent: noncovalent
                    .into_iter()
                    .map(|((atoms, attributes), presence)| {
                        (atoms, entity_span(attributes, presence))
                    })
                    .collect(),
                stereo_atoms: stereo_atoms
                    .into_iter()
                    .map(|((site, ligands, attributes), presence)| {
                        (site, ligands, entity_span(attributes, presence))
                    })
                    .collect(),
                stereo_bonds: stereo_bonds
                    .into_iter()
                    .map(|((site, ligands, attributes), presence)| {
                        (site, ligands, entity_span(attributes, presence))
                    })
                    .collect(),
                constraints: Vec::new(),
            }
        })
}

/// Constraint spans on about one in ten spans: one or two trees over the span's entity counts, each
/// Unchanged, Added, or Removed, ordered by normal form as materialization orders them. A choice
/// that fails span integrity keeps the span without constraints.
fn sparsely_constrained_reaction_span(
    entries: ReactionSpanEntries,
) -> impl Strategy<Value = ReactionSpan> {
    let span = ReactionSpan::try_from_entries(entries.clone())
        .expect("generated span entries satisfy representation integrity");
    let counts = ConstraintCounts::from_reaction_span(&span);
    prop::option::weighted(
        0.1,
        prop::collection::vec((constraint_strategy(counts), 0u8..3), 1..=2),
    )
    .prop_map(move |constraints| {
        let Some(constraints) = constraints else {
            return span.clone();
        };
        let mut constrained = entries.clone();
        constrained.constraints = constraints
            .into_iter()
            .map(|(constraint, tag)| match tag {
                0 => ConstraintSpan::Unchanged(constraint),
                1 => ConstraintSpan::Added(constraint),
                _ => ConstraintSpan::Removed(constraint),
            })
            .collect();
        constrained.constraints.sort_by_cached_key(|span| {
            let constraint = match span {
                ConstraintSpan::Unchanged(constraint)
                | ConstraintSpan::Added(constraint)
                | ConstraintSpan::Removed(constraint) => constraint,
            };
            constraint
                .clone()
                .normalize()
                .unwrap_or_else(|_| constraint.clone())
        });
        ReactionSpan::try_from_entries(constrained).unwrap_or_else(|_| span.clone())
    })
}

fn reaction_span_strategy() -> impl Strategy<Value = ReactionSpan> {
    prop_oneof![
        2 => reaction_span_entries_strategy().prop_flat_map(sparsely_constrained_reaction_span),
        3 => materializable_reaction_strategy().prop_map(|reaction| {
            reaction
                .to_reaction_span()
                .expect("generated reaction materializes a span")
        }),
    ]
}

fn union_count<Id>(correspondence: &Correspondence<Id>) -> usize
where
    Id: Copy + Ord + From<usize>,
{
    correspondence.left_count() + correspondence.right_count() - correspondence.matched_pair_count()
}

fn dense_permutation<Id>(count: usize, seed: u64, entity_kind: u32) -> Remapping<Id>
where
    Id: Copy + Debug + Into<usize> + From<usize>,
{
    if count == 0 {
        return Remapping::default();
    }
    let mut state = seed ^ (u64::from(entity_kind) + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    let mut images = (0..count).collect::<Vec<_>>();
    for index in (1..count).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        images.swap(index, state as usize % (index + 1));
    }
    let images = images.into_iter().map(Id::from).collect::<Vec<_>>();
    Remapping::new(images).unwrap()
}

fn dense_union_remapping(span: &ReactionSpan, seed: u64) -> MoleculeRemapping {
    let side = span.correspondence();
    MoleculeRemapping::new(
        GraphRemapping::new(
            dense_permutation::<NodeId>(union_count(side.atoms()), seed, 0),
            dense_permutation::<EdgeId>(union_count(side.bonds()), seed, 1),
        ),
        dense_permutation::<DativeBondId>(union_count(side.dative_bonds()), seed, 2),
        dense_permutation::<AromaticSystemId>(union_count(side.aromatic_systems()), seed, 3),
        dense_permutation::<MulticenterBondId>(union_count(side.multicenter_bonds()), seed, 4),
        dense_permutation::<NoncovalentBondId>(union_count(side.noncovalent_bonds()), seed, 5),
        dense_permutation::<StereoAtomId>(union_count(side.stereo_atoms()), seed, 6),
        dense_permutation::<StereoBondId>(union_count(side.stereo_bonds()), seed, 7),
    )
}

#[derive(Debug)]
pub(super) struct ReactionSpanScenario {
    pub(super) span: ReactionSpan,
    pub(super) first: MoleculeRemapping,
    pub(super) second: MoleculeRemapping,
}

pub(super) fn reaction_span_scenario_strategy() -> impl Strategy<Value = ReactionSpanScenario> {
    (reaction_span_strategy(), any::<u64>(), any::<u64>()).prop_map(
        |(span, first_seed, second_seed)| ReactionSpanScenario {
            first: dense_union_remapping(&span, first_seed),
            second: dense_union_remapping(&span, second_seed),
            span,
        },
    )
}

/// A reaction and a reaction derived from it for the canonical equality contract: its
/// renumbering, normal form, or reframing, its span rematerialized, or a no-op charge
/// modification of one atom (equivalent), or one added atom (usually not equivalent).
pub(super) fn reaction_canonical_pair_strategy() -> impl Strategy<Value = (Reaction, Reaction)> {
    (
        reaction_span_scenario_strategy(),
        0u8..6,
        element_strategy(),
        any::<prop::sample::Index>(),
    )
        .prop_map(|(scenario, transform, element, index)| {
            let reaction = scenario.span.to_reaction();
            let derived = match transform {
                0 => scenario.span.remap(&scenario.first).to_reaction(),
                1 => reaction
                    .clone()
                    .normalize()
                    .unwrap_or_else(|_| reaction.clone()),
                2 => reaction
                    .clone()
                    .reframe()
                    .unwrap_or_else(|_| reaction.clone()),
                3 => reaction
                    .to_reaction_span()
                    .map(|span| span.to_reaction())
                    .unwrap_or_else(|_| reaction.clone()),
                4 if reaction.lhs().atoms().count() > 0 => {
                    let id = AtomId(index.index(reaction.lhs().atoms().count()) as u32);
                    let charge = reaction.lhs().atom(id).attributes().charge.clone();
                    let mut deltas = reaction.deltas().clone();
                    deltas.push(Delta::Atom(AtomDelta::ModifyField {
                        id,
                        change: AtomFieldChange::Charge {
                            old: charge.clone(),
                            new: charge,
                        },
                    }));
                    Reaction::new(reaction.lhs().clone(), deltas)
                }
                _ => {
                    let added = reaction
                        .deltas()
                        .iter()
                        .filter(|delta| matches!(delta, Delta::Atom(AtomDelta::Add { .. })))
                        .count();
                    let mut deltas = reaction.deltas().clone();
                    deltas.push(Delta::Atom(AtomDelta::Add {
                        id: AtomId((reaction.lhs().atoms().count() + added) as u32),
                        attributes: AtomForm::from_element(element),
                    }));
                    Reaction::new(reaction.lhs().clone(), deltas)
                }
            };
            (reaction, derived)
        })
}

/// A reaction span and a span derived from it, through the reaction pairs.
pub(super) fn reaction_span_canonical_pair_strategy(
) -> impl Strategy<Value = (ReactionSpan, ReactionSpan)> {
    reaction_canonical_pair_strategy().prop_filter_map(
        "both reactions materialize",
        |(reaction, derived)| {
            Some((
                reaction.to_reaction_span().ok()?,
                derived.to_reaction_span().ok()?,
            ))
        },
    )
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    /// Both projections of every published span satisfy molecule integrity, and converting the span
    /// to its operational reaction satisfies reaction integrity.
    #[test]
    fn test_reaction_span_integrity_preservation(span in reaction_span_strategy()) {
        let lhs = span.lhs();
        let rhs = span.rhs();
        let reaction = span.to_reaction();

        prop_assert_eq!(lhs.clone().edit().finish(), Ok(lhs));
        prop_assert_eq!(rhs.clone().edit().finish(), Ok(rhs));
        prop_assert_eq!(
            Reaction::try_new(reaction.lhs().clone(), reaction.deltas().clone()),
            Ok(reaction),
        );
    }

    /// Cross-validate the two span constructions: the direct `superimpose` (Strategy A) reproduces
    /// the span the delta path (`to_reaction_span`) builds. Recover `(L, R, C)` from the delta-path
    /// span and reassemble; a mismatch flags a diff-completeness or frame gap between the paths.
    #[test]
    fn test_reaction_span_superimpose_matches_delta_path(reaction in reaction_strategy()) {
        let span = reaction
            .to_reaction_span()
            .expect("generated reaction materializes a span");
        let rebuilt = ReactionSpan::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
        prop_assert_eq!(rebuilt, Some(span));
    }

    /// `reverse` swaps the span's sides and reverses its correspondence. Constructing that span
    /// directly must reproduce the span obtained by reversing the reaction, including the union
    /// frame chosen for entities unmatched on only one side.
    #[test]
    fn test_reaction_reverse_swaps_sides(reaction in reaction_strategy()) {
        let span = reaction
            .to_reaction_span()
            .expect("generated reaction materializes a span");
        let reverse = reaction
            .reverse()
            .expect("materializable generated reaction reverses");
        let reverse_span = reverse
            .to_reaction_span()
            .expect("reversed generated reaction materializes a span");
        let expected = ReactionSpan::superimpose(
            &span.rhs(),
            &span.lhs(),
            &span.correspondence().reverse(),
        );
        prop_assert_eq!(Some(reverse_span), expected);
    }

    /// Cross-validate the two span constructions with overlays present: the direct `superimpose`
    /// reassembles the delta-path span across all overlay kinds, not just atoms/bonds.
    #[test]
    fn test_reaction_span_superimpose_matches_delta_path_overlay(
        reaction in overlay_reaction_strategy(),
    ) {
        let span = reaction
            .to_reaction_span()
            .expect("generated overlay reaction materializes a span");
        let rebuilt = ReactionSpan::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
        prop_assert_eq!(rebuilt, Some(span));
    }

    /// Reaction → span → reaction is an exact idempotent delta normal form that retains the lhs,
    /// with its molecule-level constraints compared as a set, and the materialized span, including
    /// all entity kinds.
    #[test]
    fn test_reaction_span_roundtrip(reaction in materializable_reaction_strategy()) {
        let span = reaction
            .to_reaction_span()
            .expect("generated reaction materializes a span");
        let normalized = span.to_reaction();
        let normalized_span = normalized
            .to_reaction_span()
            .expect("reaction derived from a constructed span materializes");
        let normalized_twice = normalized_span.to_reaction();

        prop_assert_eq!(
            with_normalized_constraints(normalized.lhs()),
            with_normalized_constraints(reaction.lhs()),
        );
        prop_assert_eq!(normalized_span, span);
        prop_assert_eq!(normalized_twice, normalized);
    }

    /// `from_sides` retains the lhs frame exactly. Its materialized rhs is the supplied rhs
    /// reindexed into that frame; a total induced correspondence establishes framed equivalence.
    #[test]
    fn test_reaction_from_sides_partial(sides in crossing_reaction_sides_strategy()) {
        let reaction = Reaction::from_sides(
            sides.lhs.clone(),
            sides.rhs.clone(),
            sides.atom_correspondence,
        ).expect("generated incidence is unique under the atom correspondence");
        let span = reaction.to_reaction_span().map_err(|error| {
            TestCaseError::fail(format!("reaction did not materialize: {error}"))
        })?;

        prop_assert_eq!(span.lhs(), sides.lhs);
        let projected_rhs = span.rhs();
        let projected_to_source = MoleculeCorrespondence::induce(
            &projected_rhs,
            &sides.rhs,
            sides.projected_rhs_atoms,
        ).expect("reaction-frame normalization preserves unique entity incidence");
        prop_assert!(projected_to_source.is_total());
        let remapping = MoleculeRemapping::new(
            GraphRemapping::new(
            Remapping::new(projected_to_source.atoms().matched_pairs().iter().map(|&(_, right)| NodeId::from(right)).collect()).unwrap(),
            Remapping::new(projected_to_source.bonds().matched_pairs().iter().map(|&(_, right)| EdgeId::from(right)).collect()).unwrap(),
        ),
            Remapping::new(projected_to_source.dative_bonds().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
            Remapping::new(projected_to_source.aromatic_systems().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
            Remapping::new(projected_to_source.multicenter_bonds().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
            Remapping::new(projected_to_source.noncovalent_bonds().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
            Remapping::new(projected_to_source.stereo_atoms().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
            Remapping::new(projected_to_source.stereo_bonds().matched_pairs().iter().map(|&(_, right)| right).collect()).unwrap(),
        );
        prop_assert!(projected_rhs.framed_eq_under(&sides.rhs, &remapping));
    }

    /// Independently generated, structurally valid span entries converge through direct
    /// construction, DSL render/parse, and superimposition of the two projected sides.
    #[test]
    fn test_reaction_span_from_entries_roundtrip(
        entries in reaction_span_entries_strategy(),
    ) {
        let direct = ReactionSpan::try_from_entries(entries)
            .map_err(|error| TestCaseError::fail(format!("generated entries were invalid: {error}")))?;
        let parsed = direct
            .to_string()
            .parse::<ReactionSpan>()
            .map_err(|error| TestCaseError::fail(format!("rendered span did not parse: {error}")))?;
        let superimposed =
            ReactionSpan::superimpose(&direct.lhs(), &direct.rhs(), &direct.correspondence());

        prop_assert_eq!(parsed, direct.clone());
        prop_assert_eq!(superimposed, Some(direct));
    }

    /// Span → reaction → span is idempotent. The generated domain has no redundant constraint-span
    /// occurrences, so its first reaction-derived normal form is structurally equal to the input.
    #[test]
    fn test_reaction_span_to_reaction_roundtrip(
        entries in reaction_span_entries_strategy(),
    ) {
        let span = ReactionSpan::try_from_entries(entries)
            .map_err(|error| TestCaseError::fail(format!("generated entries were invalid: {error}")))?;
        let normalized = span
            .to_reaction()
            .to_reaction_span()
            .expect("reaction derived from a constructed span materializes");
        let normalized_twice = normalized
            .to_reaction()
            .to_reaction_span()
            .expect("reaction-derived span normal form materializes");

        prop_assert_eq!(&normalized, &span);
        prop_assert_eq!(normalized_twice, normalized);
    }
}
