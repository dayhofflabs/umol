//! Colors and remapping checks of a reaction span.

use umol_graph_core::{EdgeId, NodeId};
use umol_perm::Permutation;

use super::super::electrons::ElectronCountsForm;
use super::super::entity::{DescriptionLevel, Entity};
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition, BondId};
use super::super::incidence::{Incidence, IncidenceGraph};
use super::super::ligand::StereoLigand;
use super::super::span::EntitySpan;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::StereoKind;
use super::super::symmetry::{rank_colors, InitialColors};
use super::super::traits::{AsColor, FrameTransport, Normalize};
use super::ReactionSpan;

/// Whether the remapping `images` carries every stereo entity of `source` onto a stereo entity of
/// `target` with the same site and the same form.
///
/// `images` holds the image of each node of `incidence_graph`, the incidence graph of `source`,
/// and `target` has the entity counts of `source`. A generator of `source` is such a remapping of
/// `source` onto itself.
pub(crate) fn remapping_preserves_stereo(
    source: &ReactionSpan,
    target: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    images: &[NodeId],
) -> bool {
    let entity_image =
        |entity| incidence_graph.entity(images[incidence_graph.node_of(entity).index()]);
    let frame_action = |ligands: &[StereoLigand], image_ligands: &[StereoLigand]| {
        let mapped = ligands
            .iter()
            .map(|ligand| {
                let Entity::Atom(atom_id) = entity_image(Entity::Atom(ligand.atom_id)) else {
                    return None;
                };
                Some(StereoLigand::new(atom_id, ligand.kind))
            })
            .collect::<Option<Vec<_>>>()?;
        Permutation::between(&mapped, image_ligands)
    };
    fn form_preserved<Form>(source: &Form, target: &Form, action: Permutation) -> bool
    where
        Form: Clone + FrameTransport<Action = Permutation> + Normalize,
    {
        source
            .clone()
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(target))
    }

    for id in source.stereo_atoms().ids() {
        let Entity::StereoAtom(image_id) = entity_image(Entity::StereoAtom(id)) else {
            return false;
        };
        if entity_image(Entity::Atom(source.stereo_atoms().site(id)))
            != Entity::Atom(target.stereo_atoms().site(image_id))
        {
            return false;
        }
        let Some(action) = frame_action(
            source.stereo_atoms().ligands(id),
            target.stereo_atoms().ligands(image_id),
        ) else {
            return false;
        };
        if !form_preserved(
            source.stereo_atoms().attributes(id),
            target.stereo_atoms().attributes(image_id),
            action,
        ) {
            return false;
        }
    }

    for id in source.stereo_bonds().ids() {
        let Entity::StereoBond(image_id) = entity_image(Entity::StereoBond(id)) else {
            return false;
        };
        if entity_image(Entity::Bond(source.stereo_bonds().site(id)))
            != Entity::Bond(target.stereo_bonds().site(image_id))
        {
            return false;
        }
        let Some(action) = frame_action(
            source.stereo_bonds().ligands(id),
            target.stereo_bonds().ligands(image_id),
        ) else {
            return false;
        };
        if !form_preserved(
            source.stereo_bonds().attributes(id),
            target.stereo_bonds().attributes(image_id),
            action,
        ) {
            return false;
        }
    }

    true
}

pub(crate) fn retain_stereo_preserving_generators(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    if span.stereo_atoms().count() == 0 && span.stereo_bonds().count() == 0 {
        return;
    }
    generators
        .retain(|generator| remapping_preserves_stereo(span, span, incidence_graph, generator));
}

/// The color of the incidence at `edge`: its role and, for an atom of an aromatic system or a
/// multicenter bond, the span of the electron counts of that atom on each side.
pub(crate) fn incidence_color<'a>(
    span: &ReactionSpan,
    incidence_graph: &'a IncidenceGraph,
    edge: EdgeId,
) -> (&'a Incidence, Option<EntitySpan<Option<i64>>>) {
    let incidence = incidence_graph.incidence(edge);
    let electron_count = match incidence {
        Incidence::AromaticAtom | Incidence::MulticenterAtom => {
            let [first, second] = incidence_graph.graph().edge_endpoints(edge);
            let (owner, atom) = match (
                incidence_graph.entity(first),
                incidence_graph.entity(second),
            ) {
                (Entity::Atom(atom), owner) | (owner, Entity::Atom(atom)) => (owner, atom),
                _ => {
                    unreachable!("an electron contribution connects an overlay to an atom")
                }
            };
            let sides = match owner {
                Entity::AromaticSystem(id) => {
                    let atom_position = AtomPosition::from(
                        span.aromatic_systems()
                            .atoms(id)
                            .position(|id| id == atom)
                            .expect("incidence atom belongs to its overlay"),
                    );
                    let attributes = span.aromatic_systems().attributes(id);
                    [attributes.lhs(), attributes.rhs()]
                        .map(|form| form.map(|form| form.electrons.at(atom_position)))
                }
                Entity::MulticenterBond(id) => {
                    let atom_position = AtomPosition::from(
                        span.multicenter_bonds()
                            .atoms(id)
                            .position(|id| id == atom)
                            .expect("incidence atom belongs to its overlay"),
                    );
                    let attributes = span.multicenter_bonds().attributes(id);
                    [attributes.lhs(), attributes.rhs()]
                        .map(|form| form.map(|form| form.electrons.at(atom_position)))
                }
                _ => unreachable!(
                    "electron contributions belong to aromatic or multicenter entities"
                ),
            };
            Some(match sides {
                [Some(lhs), Some(rhs)] if lhs == rhs => EntitySpan::Unchanged(lhs),
                [Some(lhs), Some(rhs)] => EntitySpan::Modified { lhs, rhs },
                [Some(lhs), None] => EntitySpan::Removed(lhs),
                [None, Some(rhs)] => EntitySpan::Added(rhs),
                [None, None] => unreachable!("an entity span is present on a side"),
            })
        }
        _ => None,
    };
    (incidence, electron_count)
}

/// The integers of the colors of the supplied span at `level`: one per entity, by incidence node,
/// and one per incidence, by incidence edge.
///
/// The color of an entity is the color of its entity span. The colors of each entity kind are
/// ranked in kind order, and the incidences continue after the entities. The color of an incidence
/// is its role and, for an atom of an aromatic system or a multicenter bond, the span of the
/// electron counts of that atom on each side. At Topology the incidences are not ranked, since the
/// topology search input does not read them.
/// `include_constraints` selects the frame-independent entity constraints in the colors and the
/// frame-relative stereo constraints in compatibility checks. Configuration compatibility is
/// checked on each present side in both modes. No form is normalized.
pub(crate) fn initial_colors(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    include_constraints: bool,
) -> Result<InitialColors, Contradiction> {
    let mut entities = vec![0; incidence_graph.graph().node_count()];
    let position = |entity: Entity| incidence_graph.node_of(entity).index();

    let atoms = span
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            (
                if include_constraints {
                    atom.as_color()
                } else {
                    atom.as_inherent_color()
                },
                position(Entity::Atom(AtomId(index as u32))),
            )
        })
        .collect();
    let mut next = rank_colors(atoms, 0, &mut entities);
    let bonds = span
        .bonds()
        .iter()
        .enumerate()
        .map(|(index, bond)| {
            (
                if include_constraints {
                    bond.as_color()
                } else {
                    bond.as_inherent_color()
                },
                position(Entity::Bond(BondId(index as u32))),
            )
        })
        .collect();
    next = rank_colors(bonds, next, &mut entities);
    if level == DescriptionLevel::Topology {
        return Ok(InitialColors {
            entities,
            incidences: Vec::new(),
        });
    }

    let dative = span
        .dative_bonds()
        .ids()
        .map(|id| {
            (
                if include_constraints {
                    span.dative_bonds().attributes(id).as_color()
                } else {
                    span.dative_bonds().attributes(id).as_inherent_color()
                },
                position(Entity::DativeBond(id)),
            )
        })
        .collect();
    next = rank_colors(dative, next, &mut entities);
    let mut aromatic = Vec::with_capacity(span.aromatic_systems().count());
    for id in span.aromatic_systems().ids() {
        let attributes = span.aromatic_systems().attributes(id);
        let atom_count = span.aromatic_systems().atoms(id).len();
        if [attributes.lhs(), attributes.rhs()]
            .into_iter()
            .flatten()
            .any(|form| matches!(&form.electrons, ElectronCountsForm::Lit(counts) if counts.len() != atom_count))
        {
            return Err(Contradiction);
        }
        aromatic.push((
            if include_constraints {
                attributes.as_color()
            } else {
                attributes.as_inherent_color()
            },
            position(Entity::AromaticSystem(id)),
        ));
    }
    next = rank_colors(aromatic, next, &mut entities);
    let mut multicenter = Vec::with_capacity(span.multicenter_bonds().count());
    for id in span.multicenter_bonds().ids() {
        let attributes = span.multicenter_bonds().attributes(id);
        let atom_count = span.multicenter_bonds().atoms(id).len();
        if [attributes.lhs(), attributes.rhs()]
            .into_iter()
            .flatten()
            .any(|form| matches!(&form.electrons, ElectronCountsForm::Lit(counts) if counts.len() != atom_count))
        {
            return Err(Contradiction);
        }
        multicenter.push((
            if include_constraints {
                attributes.as_color()
            } else {
                attributes.as_inherent_color()
            },
            position(Entity::MulticenterBond(id)),
        ));
    }
    next = rank_colors(multicenter, next, &mut entities);
    let noncovalent = span
        .noncovalent_bonds()
        .ids()
        .map(|id| {
            (
                if include_constraints {
                    span.noncovalent_bonds().attributes(id).as_color()
                } else {
                    span.noncovalent_bonds().attributes(id).as_inherent_color()
                },
                position(Entity::NoncovalentBond(id)),
            )
        })
        .collect();
    next = rank_colors(noncovalent, next, &mut entities);

    if level >= DescriptionLevel::Structure {
        let mut stereo_atoms = Vec::with_capacity(span.stereo_atoms().count());
        for id in span.stereo_atoms().ids() {
            let attributes = span.stereo_atoms().attributes(id);
            for form in [attributes.lhs(), attributes.rhs()].into_iter().flatten() {
                if let Some(kind) = form.configuration.kind() {
                    check_stereo_atom_kind(Entity::StereoAtom(id), kind)
                        .map_err(|_| Contradiction)?;
                }
            }
            let action = Permutation::identity(span.stereo_atoms().ligands(id).len());
            if include_constraints {
                attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            } else {
                for form in [attributes.lhs(), attributes.rhs()].into_iter().flatten() {
                    form.configuration
                        .clone()
                        .reframe_by(&action)
                        .ok_or(Contradiction)?;
                }
            }
            stereo_atoms.push((
                if include_constraints {
                    attributes.as_color()
                } else {
                    attributes.as_inherent_color()
                },
                position(Entity::StereoAtom(id)),
            ));
        }
        next = rank_colors(stereo_atoms, next, &mut entities);
        let mut stereo_bonds = Vec::with_capacity(span.stereo_bonds().count());
        for id in span.stereo_bonds().ids() {
            let attributes = span.stereo_bonds().attributes(id);
            for form in [attributes.lhs(), attributes.rhs()].into_iter().flatten() {
                if let Some(kind) = form.configuration.kind() {
                    check_stereo_bond_kind(Entity::StereoBond(id), kind)
                        .map_err(|_| Contradiction)?;
                }
            }
            let action = Permutation::identity(span.stereo_bonds().ligands(id).len());
            if include_constraints {
                attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            } else {
                if !StereoKind::CisTrans.class_key().space().allows(action) {
                    return Err(Contradiction);
                }
                for form in [attributes.lhs(), attributes.rhs()].into_iter().flatten() {
                    form.configuration
                        .clone()
                        .reframe_by(&action)
                        .ok_or(Contradiction)?;
                }
            }
            stereo_bonds.push((
                if include_constraints {
                    attributes.as_color()
                } else {
                    attributes.as_inherent_color()
                },
                position(Entity::StereoBond(id)),
            ));
        }
        next = rank_colors(stereo_bonds, next, &mut entities);
    }

    let incidence_colors = incidence_graph
        .incidences()
        .map(|(edge, _)| (incidence_color(span, incidence_graph, edge), edge.index()))
        .collect();
    let mut incidences = vec![0; incidence_graph.graph().edge_count()];
    rank_colors(incidence_colors, next, &mut incidences);

    Ok(InitialColors {
        entities,
        incidences,
    })
}

#[cfg(test)]
pub(crate) mod tests;
