//! Reaction-span and reaction canonicalization.

use umol_graph_core::{EdgeId, Graph, NodeId, SubdivisionNodeSource};
use umol_perm::Permutation;

use super::super::aromatic::{aromatic_system_representative_action, AromaticSystemForm};
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition, BondId};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::StereoLigand;
use super::super::multicenter::{multicenter_bond_representative_action, MulticenterBondForm};
use super::super::num::NumForm;
use super::super::reaction::Reaction;
use super::super::reaction_span::ReactionSpan;
use super::super::remap::MoleculeRemapping;
use super::super::span::EntitySpan;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{
    stereo_atom_representative_action, stereo_bond_representative_action, StereoAtomForm,
    StereoBondForm, StereoConfigurationForm,
};
use super::super::traits::{AsColor, FrameTransport, Normalize, Reframe};
use super::{
    aromatic_system_constraint_form_key, atom_constraint_form_key, atom_inherent_fields,
    bond_constraint_form_key, bond_inherent_fields, canonical_search,
    canonical_search_from_partition, canonical_search_with_generator_filter,
    constitution_partition_descriptors, constraint_key_span, constraint_span_key,
    dative_bond_constraint_form_key, electron_counts_form_key, entity_span_key, hash_value,
    incidence_key, index_key, is_empty_sequence, molecule_remapping,
    multicenter_bond_constraint_form_key, noncovalent_bond_constraint_form_key,
    noncovalent_bond_kind_form_key, normalized_entity_span_key, num_form_key, option,
    positioned_product, product, push_span_block, rank_colors, rank_initial_colors,
    reaction_canonicalize_level, reaction_span_canonicalize_level, reaction_span_counts,
    relation_span_key, sequence, stereo_atom_constraint_form_key, stereo_bond_constraint_form_key,
    stereo_configuration_form_key, stereo_kind_key, stereo_ligand_key, unpaired_electrons_form_key,
    variant, AutomorphismAdapter, AutomorphismColor, BranchOrdering, CanonicalCandidate,
    CanonicalComparisonKey, CanonicalKeyValue, CanonicalSearchOptions, CanonicalSearchResult,
    Canonicalize, CanonicalizeContext, ConstraintBlockPosition, DeferredKey, DescriptionLevel,
    EntityBlockPosition, InitialColorKey, InitialColors, OrderedPartition, PositionedKey,
    ReactionCanonicalizeError, ReactionSpanCanonicalizeError, TopologySearchInput,
};

/// The key of a leaf: the span in the entity order of the leaf, with each entity span as the
/// integer of its color.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct ColorKey {
    /// The integer of each atom.
    atoms: Vec<u32>,
    /// The two endpoint images and the integer of each bond.
    bonds: Vec<[u32; 3]>,
    /// The sorted donor images, the acceptor image, and the integer of each dative bond.
    dative: Vec<(Vec<u32>, u32, u32)>,
    /// The sorted atom images, the electron counts of each side in that order, and the integer of
    /// each aromatic system.
    aromatic: Vec<(Vec<u32>, EntitySpan<ElectronCountsForm>, u32)>,
    /// The sorted atom images, the electron counts of each side in that order, and the integer of
    /// each multicenter bond.
    multicenter: Vec<(Vec<u32>, EntitySpan<ElectronCountsForm>, u32)>,
    /// The two endpoint images and the integer of each noncovalent bond.
    noncovalent: Vec<[u32; 3]>,
    /// The site image, the ligands in the frame of the leaf, and the forms of each side in that
    /// frame of each stereo atom.
    stereo_atoms: Vec<(u32, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)>,
    /// The image of the site bond, the ligands in the frame of the leaf, and the forms of each side
    /// in that frame of each stereo bond.
    stereo_bonds: Vec<(u32, Vec<StereoLigand>, EntitySpan<StereoBondForm>)>,
}

fn entity_color_key(span: &ReactionSpan, entity: Entity) -> Result<InitialColorKey, Contradiction> {
    let (position, value) = match entity {
        Entity::Atom(id) => (
            EntityBlockPosition::ATOM,
            normalized_entity_span_key(&span.atoms()[id.index()], |attributes| {
                Ok(CanonicalKeyValue::Product(atom_inherent_fields(
                    attributes,
                )?))
            })?,
        ),
        Entity::Bond(id) => (
            EntityBlockPosition::BOND,
            normalized_entity_span_key(&span.bonds()[id.index()], |attributes| {
                Ok(CanonicalKeyValue::Product(bond_inherent_fields(
                    attributes,
                )?))
            })?,
        ),
        Entity::DativeBond(id) => (
            EntityBlockPosition::DATIVE_BOND,
            normalized_entity_span_key(span.dative_bonds().attributes(id), |attributes| {
                Ok(positioned_product([(
                    2,
                    num_form_key(&attributes.order.clone().normalize()?),
                )]))
            })?,
        ),
        Entity::AromaticSystem(id) => (
            EntityBlockPosition::AROMATIC_SYSTEM,
            normalized_entity_span_key(span.aromatic_systems().attributes(id), |attributes| {
                if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                    if counts.len() != span.aromatic_systems().atoms(id).len())
                {
                    return Err(Contradiction);
                }
                Ok(positioned_product([
                    (2, num_form_key(&attributes.charge.clone().normalize()?)),
                    (
                        3,
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ),
                ]))
            })?,
        ),
        Entity::MulticenterBond(id) => (
            EntityBlockPosition::MULTICENTER_BOND,
            normalized_entity_span_key(span.multicenter_bonds().attributes(id), |attributes| {
                if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                    if counts.len() != span.multicenter_bonds().atoms(id).len())
                {
                    return Err(Contradiction);
                }
                Ok(positioned_product([
                    (2, num_form_key(&attributes.charge.clone().normalize()?)),
                    (
                        3,
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ),
                ]))
            })?,
        ),
        Entity::NoncovalentBond(id) => (
            EntityBlockPosition::NONCOVALENT_BOND,
            normalized_entity_span_key(span.noncovalent_bonds().attributes(id), |attributes| {
                Ok(positioned_product([(
                    1,
                    noncovalent_bond_kind_form_key(&attributes.kind.clone().normalize()?),
                )]))
            })?,
        ),
        Entity::StereoAtom(id) => (
            EntityBlockPosition::STEREO_ATOM,
            normalized_entity_span_key(span.stereo_atoms().attributes(id), |attributes| {
                if let Some(kind) = attributes.configuration.kind() {
                    check_stereo_atom_kind(Entity::StereoAtom(id), kind)
                        .map_err(|_| Contradiction)?;
                }
                attributes
                    .clone()
                    .reframe_by(&Permutation::identity(
                        span.stereo_atoms().ligands(id).len(),
                    ))
                    .ok_or(Contradiction)?;
                Ok(positioned_product([(
                    2,
                    option(
                        attributes
                            .configuration
                            .clone()
                            .normalize()?
                            .kind()
                            .map(stereo_kind_key),
                    ),
                )]))
            })?,
        ),
        Entity::StereoBond(id) => (
            EntityBlockPosition::STEREO_BOND,
            normalized_entity_span_key(span.stereo_bonds().attributes(id), |attributes| {
                if let Some(kind) = attributes.configuration.kind() {
                    check_stereo_bond_kind(Entity::StereoBond(id), kind)
                        .map_err(|_| Contradiction)?;
                }
                attributes
                    .clone()
                    .reframe_by(&Permutation::identity(
                        span.stereo_bonds().ligands(id).len(),
                    ))
                    .ok_or(Contradiction)?;
                Ok(positioned_product([(
                    2,
                    option(
                        attributes
                            .configuration
                            .clone()
                            .normalize()?
                            .kind()
                            .map(stereo_kind_key),
                    ),
                )]))
            })?,
        ),
    };
    Ok(InitialColorKey::Entity { position, value })
}

fn remapping_from_order(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> MoleculeRemapping {
    let mut images = reaction_span_counts(span).map(|count| (0..count).collect::<Vec<_>>());
    let mut next = [0; 8];
    for &node in order {
        let (entity_kind, kind_id) = match incidence_graph.entity(node) {
            Entity::Atom(id) => (0, id.index()),
            Entity::Bond(id) => (1, id.index()),
            Entity::DativeBond(id) => (2, id.index()),
            Entity::AromaticSystem(id) => (3, id.index()),
            Entity::MulticenterBond(id) => (4, id.index()),
            Entity::NoncovalentBond(id) => (5, id.index()),
            Entity::StereoAtom(id) => (6, id.index()),
            Entity::StereoBond(id) => (7, id.index()),
        };
        images[entity_kind][kind_id] = next[entity_kind];
        next[entity_kind] += 1;
    }
    molecule_remapping(&images)
}

/// The remapping that sends every entity of `incidence_graph` to its image under `generator`.
fn remapping_from_generator(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    generator: &[NodeId],
) -> MoleculeRemapping {
    let mut images = reaction_span_counts(span).map(|count| (0..count).collect::<Vec<_>>());
    for node in incidence_graph.graph().node_ids() {
        let entity = incidence_graph.entity(node);
        images[entity.kind() as usize][entity.kind_id()] =
            incidence_graph.entity(generator[node.index()]).kind_id();
    }
    molecule_remapping(&images)
}

fn generator_preserves_stereo(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    generator: &[NodeId],
) -> bool {
    let entity_image =
        |entity| incidence_graph.entity(generator[incidence_graph.node_of(entity).index()]);
    let frame_action = |source: &[StereoLigand], target: &[StereoLigand]| {
        let mapped = source
            .iter()
            .map(|ligand| {
                let Entity::Atom(atom_id) = entity_image(Entity::Atom(ligand.atom_id)) else {
                    return None;
                };
                Some(StereoLigand::new(atom_id, ligand.kind))
            })
            .collect::<Option<Vec<_>>>()?;
        Permutation::between(&mapped, target)
    };
    let configuration_preserved = |source: Option<&StereoConfigurationForm>,
                                   target: Option<&StereoConfigurationForm>,
                                   action: Permutation| {
        match (source, target) {
            (Some(source), Some(target)) => source
                .clone()
                .reframe_by(&action)
                .is_some_and(|mapped| mapped.normalized_eq(target)),
            (None, None) => true,
            (Some(_), None) | (None, Some(_)) => false,
        }
    };

    for source_id in span.stereo_atoms().ids() {
        let Entity::StereoAtom(target_id) = entity_image(Entity::StereoAtom(source_id)) else {
            return false;
        };
        if entity_image(Entity::Atom(span.stereo_atoms().site(source_id)))
            != Entity::Atom(span.stereo_atoms().site(target_id))
        {
            return false;
        }
        let Some(action) = frame_action(
            span.stereo_atoms().ligands(source_id),
            span.stereo_atoms().ligands(target_id),
        ) else {
            return false;
        };
        let source = span.stereo_atoms().attributes(source_id);
        let target = span.stereo_atoms().attributes(target_id);
        if !configuration_preserved(
            source.lhs().map(|attributes| &attributes.configuration),
            target.lhs().map(|attributes| &attributes.configuration),
            action,
        ) || !configuration_preserved(
            source.rhs().map(|attributes| &attributes.configuration),
            target.rhs().map(|attributes| &attributes.configuration),
            action,
        ) {
            return false;
        }
    }

    for source_id in span.stereo_bonds().ids() {
        let Entity::StereoBond(target_id) = entity_image(Entity::StereoBond(source_id)) else {
            return false;
        };
        if entity_image(Entity::Bond(span.stereo_bonds().site(source_id)))
            != Entity::Bond(span.stereo_bonds().site(target_id))
        {
            return false;
        }
        let Some(action) = frame_action(
            span.stereo_bonds().ligands(source_id),
            span.stereo_bonds().ligands(target_id),
        ) else {
            return false;
        };
        let source = span.stereo_bonds().attributes(source_id);
        let target = span.stereo_bonds().attributes(target_id);
        if !configuration_preserved(
            source.lhs().map(|attributes| &attributes.configuration),
            target.lhs().map(|attributes| &attributes.configuration),
            action,
        ) || !configuration_preserved(
            source.rhs().map(|attributes| &attributes.configuration),
            target.rhs().map(|attributes| &attributes.configuration),
            action,
        ) {
            return false;
        }
    }

    true
}

fn retain_stereo_preserving_generators(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    if span.stereo_atoms().count() == 0 && span.stereo_bonds().count() == 0 {
        return;
    }
    generators.retain(|generator| generator_preserves_stereo(span, incidence_graph, generator));
}

/// Keep the generators that carry `span` onto itself: the span remapped by the generator equals the
/// span under framed equality.
fn retain_preserving_generators(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    generators.retain(|generator| {
        span.clone()
            .try_remap(&remapping_from_generator(span, incidence_graph, generator))
            .is_some_and(|remapped| remapped.framed_eq(span))
    });
}

fn initial_color_keys(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialColorKey>, Vec<InitialColorKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| entity_color_key(span, incidence_graph.entity(node)))
        .collect::<Result<Vec<_>, _>>()?;
    let incidence_keys = incidence_graph
        .incidences()
        .map(|(edge, incidence)| {
            let value = match incidence {
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
                    match owner {
                        Entity::AromaticSystem(id) => {
                            let position = span
                                .aromatic_systems()
                                .atoms(id)
                                .position(|id| id == atom)
                                .expect("incidence atom belongs to its overlay");
                            variant(
                                3,
                                [normalized_entity_span_key(
                                    span.aromatic_systems().attributes(id),
                                    |attributes| {
                                        let contribution = match &attributes.electrons {
                                            ElectronCountsForm::Undetermined => {
                                                NumForm::Undetermined
                                            }
                                            ElectronCountsForm::Lit(counts) => {
                                                NumForm::Lit(counts[position])
                                            }
                                        };
                                        Ok(num_form_key(&contribution))
                                    },
                                )?],
                            )
                        }
                        Entity::MulticenterBond(id) => {
                            let position = span
                                .multicenter_bonds()
                                .atoms(id)
                                .position(|id| id == atom)
                                .expect("incidence atom belongs to its overlay");
                            variant(
                                4,
                                [normalized_entity_span_key(
                                    span.multicenter_bonds().attributes(id),
                                    |attributes| {
                                        let contribution = match &attributes.electrons {
                                            ElectronCountsForm::Undetermined => {
                                                NumForm::Undetermined
                                            }
                                            ElectronCountsForm::Lit(counts) => {
                                                NumForm::Lit(counts[position])
                                            }
                                        };
                                        Ok(num_form_key(&contribution))
                                    },
                                )?],
                            )
                        }
                        _ => unreachable!(
                            "electron contributions belong to aromatic or multicenter entities"
                        ),
                    }
                }
                incidence => incidence_key(incidence),
            };
            Ok(InitialColorKey::Incidence(value))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    Ok((entity_keys, incidence_keys))
}

/// The integers of the colors of a normalized span at `level`: one per entity, by incidence node,
/// and one per incidence, by incidence edge.
///
/// The color of an entity is the color of its entity span. The colors of each entity kind are
/// ranked in kind order, and the incidences continue after the entities. The color of an incidence
/// is its role and, for an atom of an aromatic system or a multicenter bond, the span of the
/// electron counts of that atom on each side. At Topology the incidences are not ranked, since the
/// topology search input does not read them.
fn initial_colors(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
) -> Result<InitialColors, Contradiction> {
    let mut entities = vec![0; incidence_graph.graph().node_count()];
    let position = |entity: Entity| incidence_graph.node_of(entity).index();

    let atoms = span
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            (
                atom.as_color(),
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
                bond.as_color(),
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
                span.dative_bonds().attributes(id).as_color(),
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
        aromatic.push((attributes.as_color(), position(Entity::AromaticSystem(id))));
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
        multicenter.push((attributes.as_color(), position(Entity::MulticenterBond(id))));
    }
    next = rank_colors(multicenter, next, &mut entities);
    let noncovalent = span
        .noncovalent_bonds()
        .ids()
        .map(|id| {
            (
                span.noncovalent_bonds().attributes(id).as_color(),
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
            attributes
                .clone()
                .reframe_by(&Permutation::identity(
                    span.stereo_atoms().ligands(id).len(),
                ))
                .ok_or(Contradiction)?;
            stereo_atoms.push((attributes.as_color(), position(Entity::StereoAtom(id))));
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
            attributes
                .clone()
                .reframe_by(&Permutation::identity(
                    span.stereo_bonds().ligands(id).len(),
                ))
                .ok_or(Contradiction)?;
            stereo_bonds.push((attributes.as_color(), position(Entity::StereoBond(id))));
        }
        next = rank_colors(stereo_bonds, next, &mut entities);
    }

    let incidence_colors = incidence_graph
        .incidences()
        .map(|(edge, incidence)| {
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
            ((incidence, electron_count), edge.index())
        })
        .collect();
    let mut incidences = vec![0; incidence_graph.graph().edge_count()];
    rank_colors(incidence_colors, next, &mut incidences);

    Ok(InitialColors {
        entities,
        incidences,
    })
}

fn comparison_key(
    span: &ReactionSpan,
    level: DescriptionLevel,
) -> Result<CanonicalComparisonKey, Contradiction> {
    let entries = span.entries();
    let atoms = entries
        .atoms
        .iter()
        .map(|span| {
            normalized_entity_span_key(span, |attributes| {
                Ok(CanonicalKeyValue::Product(atom_inherent_fields(
                    attributes,
                )?))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = entries
        .bonds
        .iter()
        .map(|(first, second, span)| {
            let endpoints = product([
                index_key(first.index().min(second.index())),
                index_key(first.index().max(second.index())),
            ]);
            normalized_entity_span_key(span, |attributes| {
                Ok(product([
                    endpoints.clone(),
                    CanonicalKeyValue::Product(bond_inherent_fields(attributes)?),
                ]))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut entity_blocks = Vec::new();
    push_span_block(&mut entity_blocks, EntityBlockPosition::ATOM, atoms);
    push_span_block(&mut entity_blocks, EntityBlockPosition::BOND, bonds);

    if matches!(
        level,
        DescriptionLevel::Constitution | DescriptionLevel::Structure | DescriptionLevel::Full
    ) {
        let dative = entries
            .dative
            .iter()
            .map(|(donors, acceptor, span)| {
                let mut donors = donors.iter().map(|id| id.index()).collect::<Vec<_>>();
                donors.sort_unstable();
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        sequence(donors.iter().copied().map(index_key)),
                        index_key(acceptor.index()),
                        num_form_key(&attributes.order.clone().normalize()?),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let aromatic = entries
            .aromatic
            .iter()
            .map(|(atoms, span)| {
                relation_span_key(atoms, span, |attributes| {
                    Ok(product([
                        electron_counts_form_key(&attributes.electrons.clone().normalize()?),
                        num_form_key(&attributes.charge.clone().normalize()?),
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let multicenter = entries
            .multicenter
            .iter()
            .map(|(atoms, span)| {
                relation_span_key(atoms, span, |attributes| {
                    Ok(product([
                        electron_counts_form_key(&attributes.electrons.clone().normalize()?),
                        num_form_key(&attributes.charge.clone().normalize()?),
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let noncovalent = entries
            .noncovalent
            .iter()
            .map(|([first, second], span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        product([
                            index_key(first.index().min(second.index())),
                            index_key(first.index().max(second.index())),
                        ]),
                        noncovalent_bond_kind_form_key(&attributes.kind.clone().normalize()?),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        push_span_block(&mut entity_blocks, EntityBlockPosition::DATIVE_BOND, dative);
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::AROMATIC_SYSTEM,
            aromatic,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::MULTICENTER_BOND,
            multicenter,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::NONCOVALENT_BOND,
            noncovalent,
        );
    }

    if matches!(level, DescriptionLevel::Structure | DescriptionLevel::Full) {
        let stereo_atoms = entries
            .stereo_atoms
            .iter()
            .map(|(site, ligands, span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        index_key(site.index()),
                        sequence(
                            ligands
                                .iter()
                                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                        ),
                        stereo_configuration_form_key(
                            &attributes.configuration.clone().normalize()?,
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let stereo_bonds = entries
            .stereo_bonds
            .iter()
            .map(|(site, ligands, span)| {
                normalized_entity_span_key(span, |attributes| {
                    Ok(product([
                        index_key(site.index()),
                        sequence(
                            ligands
                                .iter()
                                .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                        ),
                        stereo_configuration_form_key(
                            &attributes.configuration.clone().normalize()?,
                        ),
                    ]))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::STEREO_ATOM,
            stereo_atoms,
        );
        push_span_block(
            &mut entity_blocks,
            EntityBlockPosition::STEREO_BOND,
            stereo_bonds,
        );
    }

    let constraints = if level == DescriptionLevel::Full {
        let mut blocks = Vec::new();
        macro_rules! inline_span_block {
            ($position:expr, $entries:expr, $constraints:expr, $key:expr) => {{
                let rows = $entries
                    .iter()
                    .enumerate()
                    .filter_map(|(index, entry)| {
                        let span = $constraints(entry);
                        let value = entity_span_key(&span, Clone::clone);
                        let empty = span.lhs().is_none_or(is_empty_sequence)
                            && span.rhs().is_none_or(is_empty_sequence);
                        (!empty).then(|| product([index_key(index), value]))
                    })
                    .collect::<Vec<_>>();
                if !rows.is_empty() {
                    blocks.push(PositionedKey {
                        position: $position,
                        value: sequence(rows),
                    });
                }
            }};
        }
        inline_span_block!(
            ConstraintBlockPosition::ATOM,
            entries.atoms,
            |span: &EntitySpan<AtomForm>| constraint_key_span(span, |form| form
                .constraints
                .iter()
                .map(atom_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::BOND,
            entries.bonds,
            |entry: &(AtomId, AtomId, EntitySpan<BondForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(bond_constraint_form_key)
            ),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::DATIVE_BOND,
            entries.dative,
            |entry: &(
                Vec<AtomId>,
                AtomId,
                EntitySpan<super::super::dative::DativeBondForm>
            )| constraint_key_span(&entry.2, |form| form
                .constraints
                .iter()
                .map(dative_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::AROMATIC_SYSTEM,
            entries.aromatic,
            |entry: &(
                Vec<AtomId>,
                EntitySpan<super::super::aromatic::AromaticSystemForm>
            )| {
                constraint_key_span(&entry.1, |form| {
                    form.constraints
                        .iter()
                        .map(aromatic_system_constraint_form_key)
                })
            },
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::MULTICENTER_BOND,
            entries.multicenter,
            |entry: &(
                Vec<AtomId>,
                EntitySpan<super::super::multicenter::MulticenterBondForm>
            )| constraint_key_span(&entry.1, |form| form
                .constraints
                .iter()
                .map(multicenter_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::NONCOVALENT_BOND,
            entries.noncovalent,
            |entry: &(
                [AtomId; 2],
                EntitySpan<super::super::noncovalent::NoncovalentBondForm>
            )| constraint_key_span(&entry.1, |form| form
                .constraints
                .iter()
                .map(noncovalent_bond_constraint_form_key)),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::STEREO_ATOM,
            entries.stereo_atoms,
            |entry: &(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(stereo_atom_constraint_form_key)
            ),
            Clone::clone
        );
        inline_span_block!(
            ConstraintBlockPosition::STEREO_BOND,
            entries.stereo_bonds,
            |entry: &(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)| constraint_key_span(
                &entry.2,
                |form| form.constraints.iter().map(stereo_bond_constraint_form_key)
            ),
            Clone::clone
        );

        let mut values = entries
            .constraints
            .into_iter()
            .map(|span| Ok(constraint_span_key(&span.normalize()?)))
            .collect::<Result<Vec<_>, Contradiction>>()?;
        values.sort();
        values.dedup();
        if !values.is_empty() {
            blocks.push(PositionedKey {
                position: ConstraintBlockPosition::MOLECULE,
                value: sequence(values),
            });
        }
        blocks
    } else {
        Vec::new()
    };
    Ok(CanonicalComparisonKey {
        entity_blocks,
        constraints,
    })
}

/// Construct the topology search input from the dense normalized colors of topology entities.
///
/// The vertices are the atoms of the union graph. A removed and an added bond on one atom pair
/// differ in color, so at most one of them is a direct edge.
fn topology_search_input(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    entity_colors: &[u32],
) -> TopologySearchInput {
    let bond_color = |index: usize| {
        let node = incidence_graph.node_of(Entity::Bond(BondId(index as u32)));
        entity_colors[node.index()]
    };
    let mut bond_color_counts = vec![0_usize; entity_colors.len()];
    for index in 0..span.bonds().len() {
        bond_color_counts[bond_color(index) as usize] += 1;
    }
    let direct_bond_color = bond_color_counts
        .into_iter()
        .enumerate()
        .filter(|(_, count)| *count > 0)
        .max_by(|(lhs_color, lhs_count), (rhs_color, rhs_count)| {
            lhs_count
                .cmp(rhs_count)
                .then_with(|| rhs_color.cmp(lhs_color))
        })
        .map(|(color, _)| color as u32);

    let atom_count = span.atoms().len();
    let mut colors = entity_colors[..atom_count]
        .iter()
        .copied()
        .map(AutomorphismColor::Entity)
        .collect::<Vec<_>>();
    let mut node_sources = (0..atom_count)
        .map(|index| SubdivisionNodeSource::Node(NodeId(index as u32)))
        .collect::<Vec<_>>();
    let mut edges = Vec::<[u32; 2]>::with_capacity(span.bonds().len() * 2);

    for (index, edge) in span.graph().edge_ids().enumerate() {
        let color = bond_color(index);
        let [first, second] = span.graph().edge_endpoints(edge).map(|atom| atom.0);
        if direct_bond_color == Some(color) {
            edges.push([first, second]);
        } else {
            let subdivision = NodeId(node_sources.len() as u32);
            node_sources.push(SubdivisionNodeSource::Node(subdivision));
            colors.push(AutomorphismColor::Entity(color));
            edges.push([first, subdivision.0]);
            edges.push([subdivision.0, second]);
        }
    }

    let graph = Graph::new(node_sources.len(), &edges);
    debug_assert!(graph.is_simple());
    let initial_partition = OrderedPartition::from_color_ranks(colors.iter().map(|color| {
        let AutomorphismColor::Entity(color) = *color else {
            unreachable!("the topology search input uses entity colors")
        };
        color
    }));

    TopologySearchInput {
        adapter: AutomorphismAdapter {
            graph,
            colors,
            node_sources,
            source_node_count: atom_count,
        },
        initial_partition,
    }
}

/// The Topology key of the leaf that orders the entities of `span` as `entity_order`.
///
/// It equals the comparison key of the remapped span at Topology and is read without remapping:
/// Topology has no participant frame.
fn topology_key(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    entity_order: &[NodeId],
) -> CanonicalComparisonKey {
    let mut atom_images = vec![0_usize; span.atoms().len()];
    let mut atoms = Vec::with_capacity(atom_images.len());
    let mut bonds = Vec::with_capacity(span.bonds().len());
    for &node in entity_order {
        match incidence_graph.entity(node) {
            Entity::Atom(id) => {
                atom_images[id.index()] = atoms.len();
                atoms.push(
                    normalized_entity_span_key(&span.atoms()[id.index()], |attributes| {
                        Ok(CanonicalKeyValue::Product(atom_inherent_fields(
                            attributes,
                        )?))
                    })
                    .expect("initial colors established span normalization"),
                );
            }
            Entity::Bond(id) => {
                let [first, second] = span
                    .graph()
                    .edge_endpoints(EdgeId::from(id))
                    .map(|atom| atom_images[atom.index()]);
                let endpoints =
                    product([index_key(first.min(second)), index_key(first.max(second))]);
                bonds.push(
                    normalized_entity_span_key(&span.bonds()[id.index()], |attributes| {
                        Ok(product([
                            endpoints.clone(),
                            CanonicalKeyValue::Product(bond_inherent_fields(attributes)?),
                        ]))
                    })
                    .expect("initial colors established span normalization"),
                );
            }
            _ => unreachable!("the Topology order holds atoms and bonds"),
        }
    }
    let mut entity_blocks = Vec::with_capacity(2);
    push_span_block(&mut entity_blocks, EntityBlockPosition::ATOM, atoms);
    push_span_block(&mut entity_blocks, EntityBlockPosition::BOND, bonds);
    CanonicalComparisonKey {
        entity_blocks,
        constraints: Vec::new(),
    }
}

/// The order of every entity of `level` for an atom order.
///
/// The atoms stand as the search orders them. The entities of each further kind follow with the
/// lhs-present ones first, each group by its participants. Entities of one kind and group on the
/// same participants keep their search order.
fn entity_order(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    order: &[NodeId],
) -> Vec<NodeId> {
    let node_count = incidence_graph.graph().node_count();
    let mut positions = vec![0_usize; node_count];
    let mut atom_images = vec![0_usize; span.atoms().len()];
    let mut entity_order = Vec::with_capacity(node_count);
    for (position, &node) in order.iter().enumerate() {
        positions[node.index()] = position;
        if let Entity::Atom(id) = incidence_graph.entity(node) {
            atom_images[id.index()] = entity_order.len();
            entity_order.push(node);
        }
    }
    let atom_image = |atom: AtomId| atom_images[atom.index()];
    let search_position = |entity: Entity| {
        let node = incidence_graph.node_of(entity);
        (positions[node.index()], node)
    };

    let mut bond_images = vec![0_usize; span.bonds().len()];
    let mut bonds = span
        .graph()
        .edge_ids()
        .enumerate()
        .map(|(index, edge)| {
            let [first, second] = span
                .graph()
                .edge_endpoints(edge)
                .map(|atom| atom_images[atom.index()]);
            let (position, node) = search_position(Entity::Bond(BondId(index as u32)));
            (
                matches!(span.bonds()[index], EntitySpan::Added(_)),
                (first.min(second), first.max(second)),
                position,
                node,
                index,
            )
        })
        .collect::<Vec<_>>();
    bonds.sort_unstable();
    for (image, &(_, _, _, node, index)) in bonds.iter().enumerate() {
        bond_images[index] = image;
        entity_order.push(node);
    }

    if level >= DescriptionLevel::Constitution {
        let mut dative = span
            .dative_bonds()
            .ids()
            .map(|id| {
                let mut donors = span
                    .dative_bonds()
                    .donors(id)
                    .map(atom_image)
                    .collect::<Vec<_>>();
                donors.sort_unstable();
                let (position, node) = search_position(Entity::DativeBond(id));
                (
                    matches!(span.dative_bonds().attributes(id), EntitySpan::Added(_)),
                    (donors, atom_image(span.dative_bonds().acceptor(id))),
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        dative.sort_unstable();
        entity_order.extend(dative.into_iter().map(|(_, _, _, node)| node));

        let mut aromatic = span
            .aromatic_systems()
            .ids()
            .map(|id| {
                let mut atoms = span
                    .aromatic_systems()
                    .atoms(id)
                    .map(atom_image)
                    .collect::<Vec<_>>();
                atoms.sort_unstable();
                let (position, node) = search_position(Entity::AromaticSystem(id));
                (
                    matches!(span.aromatic_systems().attributes(id), EntitySpan::Added(_)),
                    atoms,
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        aromatic.sort_unstable();
        entity_order.extend(aromatic.into_iter().map(|(_, _, _, node)| node));

        let mut multicenter = span
            .multicenter_bonds()
            .ids()
            .map(|id| {
                let mut atoms = span
                    .multicenter_bonds()
                    .atoms(id)
                    .map(atom_image)
                    .collect::<Vec<_>>();
                atoms.sort_unstable();
                let (position, node) = search_position(Entity::MulticenterBond(id));
                (
                    matches!(
                        span.multicenter_bonds().attributes(id),
                        EntitySpan::Added(_)
                    ),
                    atoms,
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        multicenter.sort_unstable();
        entity_order.extend(multicenter.into_iter().map(|(_, _, _, node)| node));

        let mut noncovalent = span
            .noncovalent_bonds()
            .ids()
            .map(|id| {
                let [first, second] = span.noncovalent_bonds().atoms(id).map(atom_image);
                let (position, node) = search_position(Entity::NoncovalentBond(id));
                (
                    matches!(
                        span.noncovalent_bonds().attributes(id),
                        EntitySpan::Added(_)
                    ),
                    (first.min(second), first.max(second)),
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        noncovalent.sort_unstable();
        entity_order.extend(noncovalent.into_iter().map(|(_, _, _, node)| node));
    }

    if level >= DescriptionLevel::Structure {
        let mut stereo_atoms = span
            .stereo_atoms()
            .ids()
            .map(|id| {
                let (position, node) = search_position(Entity::StereoAtom(id));
                (
                    matches!(span.stereo_atoms().attributes(id), EntitySpan::Added(_)),
                    atom_image(span.stereo_atoms().site(id)),
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        stereo_atoms.sort_unstable();
        entity_order.extend(stereo_atoms.into_iter().map(|(_, _, _, node)| node));

        let mut stereo_bonds = span
            .stereo_bonds()
            .ids()
            .map(|id| {
                let (position, node) = search_position(Entity::StereoBond(id));
                (
                    matches!(span.stereo_bonds().attributes(id), EntitySpan::Added(_)),
                    bond_images[span.stereo_bonds().site(id).index()],
                    position,
                    node,
                )
            })
            .collect::<Vec<_>>();
        stereo_bonds.sort_unstable();
        entity_order.extend(stereo_bonds.into_iter().map(|(_, _, _, node)| node));
    }

    entity_order
}

/// The key of the leaf that orders the entities of `span` as `entity_order`.
///
/// A row holds the participants of its entity in their new numbers, in the representative order of
/// the entity, and the integer of its color. The row of an aromatic system or a multicenter bond
/// also holds the electron counts of each side in that order, and the row of a stereo entity holds
/// the forms of each side in the frame of the leaf in place of the integer.
fn color_key(
    span: &ReactionSpan,
    incidence_graph: &IncidenceGraph,
    colors: &InitialColors,
    entity_order: &[NodeId],
) -> ColorKey {
    let mut atom_images = vec![0_u32; span.atoms().len()];
    let mut bond_images = vec![0_u32; span.bonds().len()];
    let mut key = ColorKey {
        atoms: Vec::with_capacity(atom_images.len()),
        bonds: Vec::with_capacity(bond_images.len()),
        dative: Vec::new(),
        aromatic: Vec::new(),
        multicenter: Vec::new(),
        noncovalent: Vec::new(),
        stereo_atoms: Vec::new(),
        stereo_bonds: Vec::new(),
    };
    for &node in entity_order {
        let color = colors.entities[node.index()];
        match incidence_graph.entity(node) {
            Entity::Atom(id) => {
                atom_images[id.index()] = key.atoms.len() as u32;
                key.atoms.push(color);
            }
            Entity::Bond(id) => {
                bond_images[id.index()] = key.bonds.len() as u32;
                let [first, second] = span
                    .graph()
                    .edge_endpoints(EdgeId::from(id))
                    .map(|atom| atom_images[atom.index()]);
                key.bonds
                    .push([first.min(second), first.max(second), color]);
            }
            Entity::DativeBond(id) => {
                let mut donors = span
                    .dative_bonds()
                    .donors(id)
                    .map(|atom| atom_images[atom.index()])
                    .collect::<Vec<_>>();
                donors.sort_unstable();
                key.dative.push((
                    donors,
                    atom_images[span.dative_bonds().acceptor(id).index()],
                    color,
                ));
            }
            Entity::AromaticSystem(id) => {
                let atoms = span
                    .aromatic_systems()
                    .atoms(id)
                    .map(|atom| AtomId(atom_images[atom.index()]))
                    .collect::<Vec<_>>();
                let action = aromatic_system_representative_action(atoms.clone());
                let electrons = |form: &AromaticSystemForm| {
                    form.electrons
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established the length of the electron counts")
                };
                key.aromatic.push((
                    action
                        .act(&atoms)
                        .expect("the action has the degree of the atom list")
                        .into_iter()
                        .map(|atom| atom.0)
                        .collect(),
                    match span.aromatic_systems().attributes(id) {
                        EntitySpan::Unchanged(form) => EntitySpan::Unchanged(electrons(form)),
                        EntitySpan::Removed(form) => EntitySpan::Removed(electrons(form)),
                        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
                            lhs: electrons(lhs),
                            rhs: electrons(rhs),
                        },
                        EntitySpan::Added(form) => EntitySpan::Added(electrons(form)),
                    },
                    color,
                ));
            }
            Entity::MulticenterBond(id) => {
                let atoms = span
                    .multicenter_bonds()
                    .atoms(id)
                    .map(|atom| AtomId(atom_images[atom.index()]))
                    .collect::<Vec<_>>();
                let action = multicenter_bond_representative_action(atoms.clone());
                let electrons = |form: &MulticenterBondForm| {
                    form.electrons
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established the length of the electron counts")
                };
                key.multicenter.push((
                    action
                        .act(&atoms)
                        .expect("the action has the degree of the atom list")
                        .into_iter()
                        .map(|atom| atom.0)
                        .collect(),
                    match span.multicenter_bonds().attributes(id) {
                        EntitySpan::Unchanged(form) => EntitySpan::Unchanged(electrons(form)),
                        EntitySpan::Removed(form) => EntitySpan::Removed(electrons(form)),
                        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
                            lhs: electrons(lhs),
                            rhs: electrons(rhs),
                        },
                        EntitySpan::Added(form) => EntitySpan::Added(electrons(form)),
                    },
                    color,
                ));
            }
            Entity::NoncovalentBond(id) => {
                let [first, second] = span
                    .noncovalent_bonds()
                    .atoms(id)
                    .map(|atom| atom_images[atom.index()]);
                key.noncovalent
                    .push([first.min(second), first.max(second), color]);
            }
            Entity::StereoAtom(id) => {
                let ligands = span
                    .stereo_atoms()
                    .ligands(id)
                    .iter()
                    .map(|ligand| {
                        StereoLigand::new(AtomId(atom_images[ligand.atom_id.index()]), ligand.kind)
                    })
                    .collect::<Vec<_>>();
                let action = stereo_atom_representative_action(&ligands)
                    .expect("stereo atom ligand frame has a supported degree");
                key.stereo_atoms.push((
                    atom_images[span.stereo_atoms().site(id).index()],
                    action.act(&ligands),
                    span.stereo_atoms()
                        .attributes(id)
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established compatible stereo frame actions")
                        .normalize()
                        .expect("the span is normalized"),
                ));
            }
            Entity::StereoBond(id) => {
                let ligands = span
                    .stereo_bonds()
                    .ligands(id)
                    .iter()
                    .map(|ligand| {
                        StereoLigand::new(AtomId(atom_images[ligand.atom_id.index()]), ligand.kind)
                    })
                    .collect::<Vec<_>>();
                let action = stereo_bond_representative_action(&ligands)
                    .expect("stereo bond ligand frame has a supported degree");
                key.stereo_bonds.push((
                    bond_images[span.stereo_bonds().site(id).index()],
                    action.act(&ligands),
                    span.stereo_bonds()
                        .attributes(id)
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established compatible stereo frame actions")
                        .normalize()
                        .expect("the span is normalized"),
                ));
            }
        }
    }
    key
}

/// Search the canonical entity order of `span` at `level`.
///
/// With `defer_key`, the key of the selected leaf is present only when the search compared it with
/// another leaf.
fn canonical_candidate(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
    defer_key: bool,
) -> Result<CanonicalSearchResult<Option<CanonicalComparisonKey>>, Contradiction> {
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = initial_color_keys(span, &incidence_graph)?;
    let remapped_key = |order: &[NodeId]| {
        let remapping = remapping_from_order(span, &incidence_graph, order);
        let remapped = span.clone().remap(&remapping);
        let action = remapped.representative_action();
        let reframed = remapped
            .reframe_by(&action)
            .expect("integrity established compatible reaction-span frame actions");
        comparison_key(&reframed, level).expect("initial colors established span normalization")
    };
    let build_key = |order: &[NodeId]| {
        if level == DescriptionLevel::Topology {
            let key = topology_key(span, &incidence_graph, order);
            debug_assert!(key == remapped_key(order));
            key
        } else {
            remapped_key(order)
        }
    };
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(span, &incidence_graph, level, order);
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let selected = match level {
        DescriptionLevel::Topology => {
            let entity_colors = rank_initial_colors(&entity_keys, &[]).entities;
            let TopologySearchInput {
                adapter,
                initial_partition,
            } = topology_search_input(span, &incidence_graph, &entity_colors);
            canonical_search_from_partition(
                &adapter,
                initial_partition,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
            )
        }
        DescriptionLevel::Constitution | DescriptionLevel::Full => {
            let colors = rank_initial_colors(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            let descriptors = constitution_partition_descriptors(&adapter, &entity_keys);
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
            )
        }
        DescriptionLevel::Structure => {
            let colors = rank_initial_colors(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            let descriptors = constitution_partition_descriptors(&adapter, &entity_keys);
            let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
                retain_stereo_preserving_generators(span, &incidence_graph, generators);
            };
            canonical_search_with_generator_filter(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
                &filter_generators,
            )
        }
    };
    let CanonicalCandidate { key, entity_order } = selected.candidate;
    let key = if defer_key {
        key.into_built()
    } else {
        Some(key.into_key())
    };
    Ok(CanonicalSearchResult {
        candidate: CanonicalCandidate { key, entity_order },
        #[cfg(test)]
        stats: selected.stats,
    })
}

impl ReactionSpan {
    /// The canonical form at `level` and its source-to-canonical remapping.
    fn tracked_canonicalize_at(
        self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), ReactionSpanCanonicalizeError> {
        let automorphism_pruning = level != DescriptionLevel::Full;
        let options = CanonicalSearchOptions {
            automorphism_pruning,
            branch_order: if automorphism_pruning {
                BranchOrdering::BackendCanonical
            } else {
                BranchOrdering::NodeId
            },
        };
        let normalized = self.normalize()?;
        let candidate = canonical_candidate(&normalized, level, context, options, true)?.candidate;
        let incidence_level = match level {
            DescriptionLevel::Topology => IncidenceLevel::Topology,
            DescriptionLevel::Constitution => IncidenceLevel::Constitution,
            DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
        };
        let incidence_graph = normalized.incidence_graph(incidence_level);
        let remapping =
            remapping_from_order(&normalized, &incidence_graph, &candidate.entity_order);
        let canonical = normalized.remap(&remapping).reframe()?;
        Ok((canonical, remapping))
    }

    /// The comparison key of the canonical form at `level`: the key canonical equality compares.
    fn canonical_key_at(
        &self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<CanonicalComparisonKey, ReactionSpanCanonicalizeError> {
        let automorphism_pruning = level != DescriptionLevel::Full;
        let options = CanonicalSearchOptions {
            automorphism_pruning,
            branch_order: if automorphism_pruning {
                BranchOrdering::BackendCanonical
            } else {
                BranchOrdering::NodeId
            },
        };
        if level == DescriptionLevel::Full {
            let normalized = self.clone().normalize()?;
            return Ok(
                canonical_candidate(&normalized, level, context, options, false)?
                    .candidate
                    .key
                    .expect("the key was requested"),
            );
        }
        Ok(canonical_candidate(self, level, context, options, false)?
            .candidate
            .key
            .expect("the key was requested"))
    }
}

impl Canonicalize for ReactionSpan {
    type Error = ReactionSpanCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = reaction_span_canonicalize_level(&self);
        self.tracked_canonicalize_at(level, context)
            .map(|(canonical, _)| canonical)
    }

    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error> {
        let level = reaction_span_canonicalize_level(&self);
        self.tracked_canonicalize_at(level, context)
    }

    fn canonical_hash(&self, context: &CanonicalizeContext) -> Result<u64, Self::Error> {
        let level = reaction_span_canonicalize_level(self);
        self.canonical_key_at(level, context)
            .map(|key| hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level =
            reaction_span_canonicalize_level(self).max(reaction_span_canonicalize_level(other));
        match (
            self.canonical_key_at(level, context),
            other.canonical_key_at(level, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
                Err(ReactionSpanCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

impl Reaction {
    /// The canonical form at `level`, read back from the canonical reaction span, and the
    /// remapping between the two span union frames.
    fn tracked_canonicalize_at(
        self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), ReactionCanonicalizeError> {
        let (canonical, remapping) = self
            .to_reaction_span()?
            .tracked_canonicalize_at(level, context)?;
        Ok((canonical.to_reaction(), remapping))
    }

    /// The comparison key of the canonical reaction span at `level`.
    fn canonical_key_at(
        &self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<CanonicalComparisonKey, ReactionCanonicalizeError> {
        Ok(self.to_reaction_span()?.canonical_key_at(level, context)?)
    }
}

impl Canonicalize for Reaction {
    type Error = ReactionCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = reaction_canonicalize_level(&self);
        self.tracked_canonicalize_at(level, context)
            .map(|(canonical, _)| canonical)
    }

    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error> {
        let level = reaction_canonicalize_level(&self);
        self.tracked_canonicalize_at(level, context)
    }

    fn canonical_hash(&self, context: &CanonicalizeContext) -> Result<u64, Self::Error> {
        let level = reaction_canonicalize_level(self);
        self.canonical_key_at(level, context)
            .map(|key| hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level = reaction_canonicalize_level(self).max(reaction_canonicalize_level(other));
        match (
            self.canonical_key_at(level, context),
            other.canonical_key_at(level, context),
        ) {
            (Ok(left), Ok(right)) => left == right,
            (
                Err(ReactionCanonicalizeError::Contradiction(_)),
                Err(ReactionCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests;
