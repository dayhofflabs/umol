//! Reaction-span bridge properties over three domains: generated reactions cross-validate the
//! delta and superimposition paths; explicitly reindexed molecule pairs exercise lhs-anchored
//! normalization under crossing partial correspondences; generated span entries cross-validate
//! direct construction, DSL parsing, and superimposition.

use std::iter::once;

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondId, EntitySpan, Molecule, MoleculeCorrespondence, MoleculeEntries,
    Reaction, ReactionSpanAst, ReactionSpanEntries, StereoLigand,
};

use crate::strategies::{
    comprehensive_reaction_strategy, molecule_entries_strategy,
    molecule_entries_structurally_unambiguous_strategy, overlay_reaction_strategy,
    reaction_strategy,
};

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
                .map(|(first, second, ast)| (reverse_atom(first), reverse_atom(second), ast))
                .collect(),
            dative: dative
                .into_iter()
                .map(|(donors, acceptor, ast)| {
                    (
                        donors.into_iter().map(reverse_atom).collect(),
                        reverse_atom(acceptor),
                        ast,
                    )
                })
                .collect(),
            aromatic: aromatic
                .into_iter()
                .map(|(atoms, ast)| (atoms.into_iter().map(reverse_atom).collect(), ast))
                .collect(),
            multicenter: multicenter
                .into_iter()
                .map(|(atoms, ast)| (atoms.into_iter().map(reverse_atom).collect(), ast))
                .collect(),
            noncovalent: noncovalent
                .into_iter()
                .map(|(first, second, ast)| (reverse_atom(first), reverse_atom(second), ast))
                .collect(),
            stereo_atoms: stereo_atoms
                .into_iter()
                .map(|(site, ligands, ast)| {
                    (
                        reverse_atom(site),
                        ligands
                            .into_iter()
                            .map(|ligand| {
                                StereoLigand::new(reverse_atom(ligand.atom_id), ligand.kind)
                            })
                            .collect(),
                        ast,
                    )
                })
                .collect(),
            stereo_bonds: stereo_bonds
                .into_iter()
                .map(|(site, ligands, ast)| {
                    (
                        site,
                        ligands
                            .into_iter()
                            .map(|ligand| {
                                StereoLigand::new(reverse_atom(ligand.atom_id), ligand.kind)
                            })
                            .collect(),
                        ast,
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
                |(old_id, (first, second, ast))| {
                    bond_presence[old_id].map(|presence| ((old_id, first, second, ast), presence))
                },
            ));
            let mut bond_ids = vec![None; bond_presence.len()];
            for (new_id, ((old_id, ..), _)) in bonds.iter().enumerate() {
                bond_ids[*old_id] = Some(BondId(new_id as u32));
            }

            let dative = lhs_anchored(entries.dative.into_iter().filter_map(
                |(donors, acceptor, ast)| {
                    let presence = intersection_presence(
                        donors
                            .iter()
                            .copied()
                            .chain(once(acceptor))
                            .map(presence_of_atom),
                    )?;
                    Some(((donors, acceptor, ast), presence))
                },
            ));
            let aromatic = lhs_anchored(entries.aromatic.into_iter().filter_map(|(atoms, ast)| {
                let presence = intersection_presence(atoms.iter().copied().map(presence_of_atom))?;
                Some(((atoms, ast), presence))
            }));
            let multicenter =
                lhs_anchored(entries.multicenter.into_iter().filter_map(|(atoms, ast)| {
                    let presence =
                        intersection_presence(atoms.iter().copied().map(presence_of_atom))?;
                    Some(((atoms, ast), presence))
                }));
            let noncovalent = lhs_anchored(entries.noncovalent.into_iter().filter_map(
                |(first, second, ast)| {
                    let presence =
                        intersection_presence([presence_of_atom(first), presence_of_atom(second)])?;
                    Some(((first, second, ast), presence))
                },
            ));
            let stereo_atoms = lhs_anchored(entries.stereo_atoms.into_iter().filter_map(
                |(site, ligands, ast)| {
                    let presence = intersection_presence(
                        once(site)
                            .chain(ligands.iter().map(|ligand| ligand.atom_id))
                            .map(presence_of_atom),
                    )?;
                    Some(((site, ligands, ast), presence))
                },
            ));
            let stereo_bonds = lhs_anchored(entries.stereo_bonds.into_iter().filter_map(
                |(site, ligands, ast)| {
                    let site_presence = bond_presence.get(site.0 as usize).copied().flatten()?;
                    let site = bond_ids.get(site.0 as usize).copied().flatten()?;
                    let presence = intersection_presence(
                        once(site_presence).chain(
                            ligands
                                .iter()
                                .map(|ligand| presence_of_atom(ligand.atom_id)),
                        ),
                    )?;
                    Some(((site, ligands, ast), presence))
                },
            ));

            ReactionSpanEntries {
                atoms: entries
                    .atoms
                    .into_iter()
                    .zip(atom_presence)
                    .map(|(ast, presence)| entity_span(ast, presence))
                    .collect(),
                bonds: bonds
                    .into_iter()
                    .map(|((_, first, second, ast), presence)| {
                        (first, second, entity_span(ast, presence))
                    })
                    .collect(),
                dative: dative
                    .into_iter()
                    .map(|((donors, acceptor, ast), presence)| {
                        (donors, acceptor, entity_span(ast, presence))
                    })
                    .collect(),
                aromatic: aromatic
                    .into_iter()
                    .map(|((atoms, ast), presence)| (atoms, entity_span(ast, presence)))
                    .collect(),
                multicenter: multicenter
                    .into_iter()
                    .map(|((atoms, ast), presence)| (atoms, entity_span(ast, presence)))
                    .collect(),
                noncovalent: noncovalent
                    .into_iter()
                    .map(|((first, second, ast), presence)| {
                        (first, second, entity_span(ast, presence))
                    })
                    .collect(),
                stereo_atoms: stereo_atoms
                    .into_iter()
                    .map(|((site, ligands, ast), presence)| {
                        (site, ligands, entity_span(ast, presence))
                    })
                    .collect(),
                stereo_bonds: stereo_bonds
                    .into_iter()
                    .map(|((site, ligands, ast), presence)| {
                        (site, ligands, entity_span(ast, presence))
                    })
                    .collect(),
                constraints: Vec::new(),
            }
        })
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    /// Cross-validate the two span constructions: the direct `superimpose` (Strategy A) reproduces
    /// the span the delta path (`to_reaction_span`) builds. Recover `(L, R, C)` from the delta-path
    /// span and reassemble; a mismatch flags a diff-completeness or frame gap between the paths.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, Some(span));
        }
    }

    /// `reverse` swaps the span's sides and reverses its correspondence. Constructing that span
    /// directly must reproduce the span obtained by reversing the reaction, including the union
    /// frame chosen for entities unmatched on only one side.
    #[test]
    fn test_reaction_ast_reverse_swaps_sides(reaction in reaction_strategy()) {
        if let (Ok(span), Ok(reverse)) = (reaction.to_reaction_span(), reaction.reverse()) {
            if let Ok(reverse_span) = reverse.to_reaction_span() {
                let expected = ReactionSpanAst::superimpose(
                    &span.rhs(),
                    &span.lhs(),
                    &span.correspondence().reverse(),
                );
                prop_assert_eq!(Some(reverse_span), expected);
            }
        }
    }

    /// Cross-validate the two span constructions with overlays present: the direct `superimpose`
    /// reassembles the delta-path span across all overlay families, not just atoms/bonds.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path_overlay(
        reaction in overlay_reaction_strategy(),
    ) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, Some(span));
        }
    }

    /// Reaction → span → reaction may normalize relative deltas into absolute updates. The
    /// resulting reaction nevertheless materializes the same span, including all overlay families.
    #[test]
    fn test_reaction_ast_span_roundtrip(reaction in comprehensive_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            if let Ok(rebuilt) = span.to_reaction().to_reaction_span() {
                prop_assert_eq!(rebuilt, span);
            }
        }
    }

    /// `from_sides` retains the lhs frame exactly. Its materialized rhs is the supplied rhs
    /// reindexed into that frame; a total induced correspondence establishes framed equivalence.
    #[test]
    fn test_reaction_ast_from_sides_partial(sides in crossing_reaction_sides_strategy()) {
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
        prop_assert!(projected_rhs.equiv_under(&sides.rhs, &projected_to_source));
    }

    /// Independently generated, structurally valid span entries converge through direct
    /// construction, DSL render/parse, and superimposition of the two projected sides.
    #[test]
    fn test_reaction_span_ast_from_entries_roundtrip(
        entries in reaction_span_entries_strategy(),
    ) {
        let direct = ReactionSpanAst::try_from_entries(entries).map_err(|error| {
            TestCaseError::fail(format!("generated entries were invalid: {error}"))
        })?;
        let parsed = direct.to_string().parse::<ReactionSpanAst>().map_err(|error| {
            TestCaseError::fail(format!("rendered span did not parse: {error}"))
        })?;
        let superimposed = ReactionSpanAst::superimpose(
            &direct.lhs(),
            &direct.rhs(),
            &direct.correspondence(),
        );

        prop_assert_eq!(parsed, direct.clone());
        prop_assert_eq!(superimposed, Some(direct));
    }
}
