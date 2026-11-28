//! Colors and remapping checks of a reaction span.

use std::collections::BTreeSet;

use umol_graph_core::{AutomorphismAlgorithm, EdgeId, NodeId};
use umol_perm::Permutation;

use super::super::edit::collect_constraint_entities;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::{DescriptionLevel, Entity};
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition, BondId, StereoAtomId, StereoBondId};
use super::super::incidence::{
    Incidence, IncidenceEdge, IncidenceEdgeId, IncidenceGraph, IncidenceLevel, IncidenceNodeId,
};
use super::super::ligand::StereoLigand;
use super::super::span::{ConstraintSpan, EntitySpan};
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::StereoKind;
use super::super::symmetry::{
    extract_subgroup_generators, rank_colors, AutomorphismAdapter, InitialColors,
};
use super::super::traits::{AsColor, FrameTransport, Normalize, Reframe};
use super::canonicalize::remapping_from_images;
use super::ReactionSpan;

#[derive(Clone, Copy)]
enum SpanStereoRecordId {
    Atom(StereoAtomId),
    Bond(StereoBondId),
}

fn form_preserved<Form>(source: &Form, target: &Form, action: Permutation) -> bool
where
    Form: Clone + FrameTransport<Action = Permutation> + Normalize,
{
    source
        .clone()
        .reframe_by(&action)
        .is_some_and(|mapped| mapped.normalized_eq(target))
}

fn stereo_record_preserved(
    source: &ReactionSpan,
    target: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    image_of: impl Fn(usize) -> u32,
    record: SpanStereoRecordId,
) -> bool {
    let entity_image = |entity: Entity| {
        let node = incidence_graph
            .node_of(entity)
            .expect("stereo remapping entity is in incidence graph");
        incidence_graph.entity(IncidenceNodeId(image_of(node.index())))
    };
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

    match record {
        SpanStereoRecordId::Atom(id) => {
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
            form_preserved(
                source.stereo_atoms().attributes(id),
                target.stereo_atoms().attributes(image_id),
                action,
            )
        }
        SpanStereoRecordId::Bond(id) => {
            let Entity::StereoBond(image_id) = entity_image(Entity::StereoBond(id)) else {
                return false;
            };
            let [s0, s1] = source
                .graph()
                .edge_endpoints(EdgeId::from(source.stereo_bonds().site(id)))
                .map(AtomId::from);
            let (Entity::Atom(m0), Entity::Atom(m1)) = (
                entity_image(Entity::Atom(s0)),
                entity_image(Entity::Atom(s1)),
            ) else {
                return false;
            };
            let [t0, t1] = target
                .graph()
                .edge_endpoints(EdgeId::from(target.stereo_bonds().site(image_id)))
                .map(AtomId::from);
            if !((m0 == t0 && m1 == t1) || (m0 == t1 && m1 == t0)) {
                return false;
            }
            let Some(action) = frame_action(
                source.stereo_bonds().ligands(id),
                target.stereo_bonds().ligands(image_id),
            ) else {
                return false;
            };
            form_preserved(
                source.stereo_bonds().attributes(id),
                target.stereo_bonds().attributes(image_id),
                action,
            )
        }
    }
}

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
    let image_of = |index: usize| images[index].0;
    for id in source.stereo_atoms().ids() {
        if !stereo_record_preserved(
            source,
            target,
            incidence_graph,
            image_of,
            SpanStereoRecordId::Atom(id),
        ) {
            return false;
        }
    }
    for id in source.stereo_bonds().ids() {
        if !stereo_record_preserved(
            source,
            target,
            incidence_graph,
            image_of,
            SpanStereoRecordId::Bond(id),
        ) {
            return false;
        }
    }
    true
}

impl ReactionSpan {
    /// Generating set for the proper (delta-, overlay-, stereo-, and constraint-preserving)
    /// automorphism subgroup of this span over `IncidenceLevel::Full`, paired with the incidence
    /// graph on which the generators act.
    ///
    /// Returns `None` when initial coloring rejects malformed span attributes or when a stereo or
    /// constraint check encounters a contradiction; callers treat `None` as the trivial subgroup
    /// `{1}`.
    pub(crate) fn proper_incidence_generators(
        &self,
        algorithm: AutomorphismAlgorithm,
    ) -> Option<(IncidenceGraph, Vec<Vec<NodeId>>)> {
        let incidence = self.incidence_graph(IncidenceLevel::Full);
        let node_count = incidence.node_count();
        if node_count <= 1 {
            return Some((incidence, Vec::new()));
        }

        let colors = initial_colors(self, &incidence, DescriptionLevel::Structure, true).ok()?;
        let adapter = AutomorphismAdapter::new(&incidence, &colors);
        let group = adapter.graph().automorphisms(
            |node| adapter.node_color(node),
            |edge| adapter.edge_color(edge),
            algorithm,
        );
        let generators = group.generators();
        if generators.is_empty() {
            return Some((incidence, Vec::new()));
        }

        let has_stereo = self.stereo_atoms().count() > 0 || self.stereo_bonds().count() > 0;
        let has_constraints = !self.constraints().is_empty();
        if !has_stereo && !has_constraints {
            return Some((incidence, generators.to_vec()));
        }

        if has_stereo {
            let stereo_normalizable = self.stereo_atoms().ids().all(|id| {
                self.stereo_atoms()
                    .attributes(id)
                    .clone()
                    .normalize()
                    .is_ok()
            }) && self.stereo_bonds().ids().all(|id| {
                self.stereo_bonds()
                    .attributes(id)
                    .clone()
                    .normalize()
                    .is_ok()
            });
            if !stereo_normalizable {
                return None;
            }
        }
        if has_constraints && self.clone().reframe().is_err() {
            return None;
        }

        let preserves_constraints = |node_images: &[NodeId]| -> bool {
            let remapping = remapping_from_images(self, self, &incidence, node_images);
            self.framed_eq_under(self, &remapping)
        };
        if generators.iter().all(|generator| {
            (!has_stereo || remapping_preserves_stereo(self, self, &incidence, generator))
                && (!has_constraints || preserves_constraints(generator))
        }) {
            return Some((incidence, generators.to_vec()));
        }

        let mut base = Vec::new();
        let mut base_pos = vec![usize::MAX; node_count];
        let mut add_base = |node: IncidenceNodeId| -> usize {
            let idx = node.index();
            if base_pos[idx] == usize::MAX {
                base_pos[idx] = base.len();
                base.push(idx);
            }
            base_pos[idx]
        };

        let mut tagged_records =
            Vec::with_capacity(self.stereo_atoms().count() + self.stereo_bonds().count());
        for id in self.stereo_atoms().ids() {
            let mut max_pos = add_base(
                incidence
                    .node_of(Entity::StereoAtom(id))
                    .expect("stereo atom is in full incidence graph"),
            );
            max_pos = max_pos.max(add_base(
                incidence
                    .node_of(Entity::Atom(self.stereo_atoms().site(id)))
                    .expect("stereo atom site is in incidence graph"),
            ));
            for ligand in self.stereo_atoms().ligands(id) {
                max_pos = max_pos.max(add_base(
                    incidence
                        .node_of(Entity::Atom(ligand.atom_id))
                        .expect("stereo atom ligand is in incidence graph"),
                ));
            }
            tagged_records.push((max_pos, SpanStereoRecordId::Atom(id)));
        }
        for id in self.stereo_bonds().ids() {
            let mut max_pos = add_base(
                incidence
                    .node_of(Entity::StereoBond(id))
                    .expect("stereo bond is in full incidence graph"),
            );
            let [first, second] = self
                .graph()
                .edge_endpoints(EdgeId::from(self.stereo_bonds().site(id)))
                .map(AtomId::from);
            max_pos = max_pos.max(add_base(
                incidence
                    .node_of(Entity::Atom(first))
                    .expect("stereo bond endpoint is in incidence graph"),
            ));
            max_pos = max_pos.max(add_base(
                incidence
                    .node_of(Entity::Atom(second))
                    .expect("stereo bond endpoint is in incidence graph"),
            ));
            for ligand in self.stereo_bonds().ligands(id) {
                max_pos = max_pos.max(add_base(
                    incidence
                        .node_of(Entity::Atom(ligand.atom_id))
                        .expect("stereo bond ligand is in incidence graph"),
                ));
            }
            tagged_records.push((max_pos, SpanStereoRecordId::Bond(id)));
        }

        if has_constraints {
            let mut referenced = BTreeSet::new();
            for span in self.constraints() {
                match span {
                    ConstraintSpan::Unchanged(constraint)
                    | ConstraintSpan::Added(constraint)
                    | ConstraintSpan::Removed(constraint) => {
                        collect_constraint_entities(constraint, &mut referenced);
                    }
                }
            }
            for entity in referenced {
                match entity {
                    Entity::Bond(bond) => {
                        let [first, second] = self
                            .graph()
                            .edge_endpoints(EdgeId::from(bond))
                            .map(AtomId::from);
                        add_base(
                            incidence
                                .node_of(Entity::Atom(first))
                                .expect("constraint bond endpoint is in incidence graph"),
                        );
                        add_base(
                            incidence
                                .node_of(Entity::Atom(second))
                                .expect("constraint bond endpoint is in incidence graph"),
                        );
                    }
                    other => {
                        add_base(
                            incidence
                                .node_of(other)
                                .expect("constraint entity is in full incidence graph"),
                        );
                    }
                }
            }
        }

        if base.is_empty() {
            return None;
        }

        let mut records_at_depth = vec![Vec::new(); base.len()];
        for (depth, record) in tagged_records {
            records_at_depth[depth].push(record);
        }
        let final_depth = base.len() - 1;
        let mut node_buf = if has_constraints {
            vec![NodeId(0); node_count]
        } else {
            Vec::new()
        };

        let proper = extract_subgroup_generators(node_count, generators, &base, |depth, images| {
            let stereo_ok = records_at_depth[depth].iter().all(|&record| {
                stereo_record_preserved(self, self, &incidence, |idx| images[idx], record)
            });
            if !stereo_ok {
                return false;
            }
            if has_constraints && depth == final_depth {
                for (slot, &img) in node_buf.iter_mut().zip(images) {
                    *slot = NodeId(img);
                }
                preserves_constraints(&node_buf)
            } else {
                true
            }
        });

        debug_assert!(proper.iter().all(|generator| {
            (!has_stereo || remapping_preserves_stereo(self, self, &incidence, generator))
                && (!has_constraints || preserves_constraints(generator))
        }));
        Some((incidence, proper))
    }
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

/// The color of the overlay incidence at `edge`: its role and, for an atom of an aromatic system
/// or a multicenter bond, the span of the electron counts of that atom on each side.
pub(crate) fn incidence_color(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    edge: IncidenceEdgeId,
) -> (Incidence, Option<EntitySpan<Option<i64>>>) {
    let IncidenceEdge::Incidence(incidence) = incidence_graph.edge(edge) else {
        unreachable!("incidence_color is called only on overlay incidence edges");
    };
    let electron_count = match incidence {
        Incidence::AromaticAtom | Incidence::MulticenterAtom => {
            let (owner, atom) = incidence_graph
                .incidence_endpoints(edge)
                .expect("overlay incidence connects an overlay to an atom");
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

/// The integers of the colors of the supplied span at `level`: one per incidence node (atoms +
/// included overlays), and one per incidence edge (localized bonds followed by overlay
/// incidences).
///
/// The colors of each entity kind are ranked in kind order, followed by localized bonds and then
/// overlay incidences in disjoint integer ranges. The color of an overlay incidence is its role
/// and, for an atom of an aromatic system or a multicenter bond, the span of the electron counts
/// of that atom on each side.
/// `include_constraints` selects the frame-independent entity constraints in the colors and the
/// frame-relative stereo constraints in compatibility checks. Configuration compatibility is
/// checked on each present side in both modes. No form is normalized.
#[inline]
pub(crate) fn initial_colors(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    include_constraints: bool,
) -> Result<InitialColors, Contradiction> {
    let mut entities = vec![0; incidence_graph.node_count()];
    let mut edges = vec![0; incidence_graph.edge_count()];
    let position = |entity: Entity| {
        incidence_graph
            .node_of(entity)
            .expect("entity belongs to incidence graph at this level")
            .index()
    };

    let atoms = span
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| (atom, position(Entity::Atom(AtomId(index as u32)))));
    let mut next = if include_constraints {
        rank_colors(
            atoms
                .map(|(form, position)| (form.as_color(), position))
                .collect(),
            0,
            &mut entities,
        )
    } else {
        rank_colors(
            atoms
                .map(|(form, position)| (form.as_field_color(), position))
                .collect(),
            0,
            &mut entities,
        )
    };

    if level >= DescriptionLevel::Constitution {
        let dative = span.dative_bonds().ids().map(|id| {
            (
                span.dative_bonds().attributes(id),
                position(Entity::DativeBond(id)),
            )
        });
        next = if include_constraints {
            rank_colors(
                dative
                    .map(|(form, position)| (form.as_color(), position))
                    .collect(),
                next,
                &mut entities,
            )
        } else {
            rank_colors(
                dative
                    .map(|(form, position)| (form.as_field_color(), position))
                    .collect(),
                next,
                &mut entities,
            )
        };

        let aromatic = span.aromatic_systems().ids().map(|id| {
            let attributes = span.aromatic_systems().attributes(id);
            let atom_count = span.aromatic_systems().atoms(id).len();
            if [attributes.lhs(), attributes.rhs()].into_iter().flatten().any(|form| {
                matches!(&form.electrons, ElectronCountsForm::Lit(counts) if counts.len() != atom_count)
            }) {
                return Err(Contradiction);
            }
            Ok((attributes, position(Entity::AromaticSystem(id))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(span.aromatic_systems().count());
            for entry in aromatic {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(span.aromatic_systems().count());
            for entry in aromatic {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };

        let multicenter = span.multicenter_bonds().ids().map(|id| {
            let attributes = span.multicenter_bonds().attributes(id);
            let atom_count = span.multicenter_bonds().atoms(id).len();
            if [attributes.lhs(), attributes.rhs()].into_iter().flatten().any(|form| {
                matches!(&form.electrons, ElectronCountsForm::Lit(counts) if counts.len() != atom_count)
            }) {
                return Err(Contradiction);
            }
            Ok((attributes, position(Entity::MulticenterBond(id))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(span.multicenter_bonds().count());
            for entry in multicenter {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(span.multicenter_bonds().count());
            for entry in multicenter {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };

        let noncovalent = span.noncovalent_bonds().ids().map(|id| {
            (
                span.noncovalent_bonds().attributes(id),
                position(Entity::NoncovalentBond(id)),
            )
        });
        next = if include_constraints {
            rank_colors(
                noncovalent
                    .map(|(form, position)| (form.as_color(), position))
                    .collect(),
                next,
                &mut entities,
            )
        } else {
            rank_colors(
                noncovalent
                    .map(|(form, position)| (form.as_field_color(), position))
                    .collect(),
                next,
                &mut entities,
            )
        };
    }

    if level >= DescriptionLevel::Structure {
        let stereo_atoms = span.stereo_atoms().ids().map(|id| {
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
            Ok((attributes, position(Entity::StereoAtom(id))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(span.stereo_atoms().count());
            for entry in stereo_atoms {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(span.stereo_atoms().count());
            for entry in stereo_atoms {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };

        let stereo_bonds = span.stereo_bonds().ids().map(|id| {
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
            Ok((attributes, position(Entity::StereoBond(id))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(span.stereo_bonds().count());
            for entry in stereo_bonds {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(span.stereo_bonds().count());
            for entry in stereo_bonds {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };
    }

    let bonds = span.bonds().iter().enumerate().map(|(index, bond)| {
        (
            bond,
            incidence_graph
                .bond_edge(BondId(index as u32))
                .expect("bond belongs to incidence graph")
                .index(),
        )
    });
    next = if include_constraints {
        rank_colors(
            bonds
                .map(|(form, position)| (form.as_color(), position))
                .collect(),
            next,
            &mut edges,
        )
    } else {
        rank_colors(
            bonds
                .map(|(form, position)| (form.as_field_color(), position))
                .collect(),
            next,
            &mut edges,
        )
    };

    if incidence_graph.incidence_count() > 0 {
        let incidence_colors = incidence_graph
            .incidences()
            .map(|(edge, _)| (incidence_color(span, incidence_graph, edge), edge.index()))
            .collect();
        rank_colors(incidence_colors, next, &mut edges);
    }

    Ok(InitialColors { entities, edges })
}

#[cfg(test)]
mod tests;
