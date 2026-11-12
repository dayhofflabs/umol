//! Reaction-span and reaction canonicalization.

use umol_graph_core::{EdgeId, Graph, NodeId, SubdivisionNodeSource};
use umol_perm::Permutation;

use super::super::aromatic::{aromatic_system_representative_action, AromaticSystemForm};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition, BondId};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::StereoLigand;
use super::super::multicenter::{multicenter_bond_representative_action, MulticenterBondForm};
use super::super::reaction::Reaction;
use super::super::reaction_span::ReactionSpan;
use super::super::remap::MoleculeRemapping;
use super::super::span::{ConstraintSpan, EntitySpan};
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{
    stereo_atom_representative_action, stereo_bond_representative_action, StereoAtomForm,
    StereoBondForm,
};
use super::super::traits::{AsColor, FrameTransport, Normalize, Reframe};
use super::{
    canonical_search, canonical_search_from_partition, canonical_search_with_generator_filter,
    hash_value, molecule_remapping, rank_colors, reaction_canonicalize_level,
    reaction_span_canonicalize_level, reaction_span_counts, AutomorphismAdapter, AutomorphismColor,
    CanonicalCandidate, CanonicalSearchResult, Canonicalize, CanonicalizeContext, DeferredKey,
    DescriptionLevel, InitialColors, OrderedPartition, ReactionCanonicalizeError,
    ReactionSpanCanonicalizeError, TopologySearchInput,
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
    /// The constraint spans in the entity order and the frames of the leaf.
    constraints: Vec<ConstraintSpan>,
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
    fn form_preserved<Form>(source: &Form, target: &Form, action: Permutation) -> bool
    where
        Form: Clone + FrameTransport<Action = Permutation> + Normalize,
    {
        source
            .clone()
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(target))
    }

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
        if !form_preserved(
            span.stereo_atoms().attributes(source_id),
            span.stereo_atoms().attributes(target_id),
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
        if !form_preserved(
            span.stereo_bonds().attributes(source_id),
            span.stereo_bonds().attributes(target_id),
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
/// the forms of each side in the frame of the leaf in place of the integer. The constraint spans
/// follow, read from the span remapped and reframed in the order of the leaf.
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
        constraints: Vec::new(),
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
    if !span.constraints().is_empty() {
        key.constraints = span
            .clone()
            .remap(&remapping_from_order(span, incidence_graph, entity_order))
            .reframe()
            .expect("the span is normalized")
            .constraints()
            .to_vec();
    }
    key
}

/// Search the canonical entity order of a normalized `span` at `level`.
///
/// The key of the selected leaf is present only when the search compared it with another leaf.
fn canonical_candidate(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
) -> Result<CanonicalSearchResult<Option<ColorKey>>, Contradiction> {
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let colors = initial_colors(span, &incidence_graph, level)?;
    let build_key = |order: &[NodeId]| color_key(span, &incidence_graph, &colors, order);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(span, &incidence_graph, level, order);
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let selected = match level {
        DescriptionLevel::Topology => {
            let TopologySearchInput {
                adapter,
                initial_partition,
            } = topology_search_input(span, &incidence_graph, &colors.entities);
            canonical_search_from_partition(
                &adapter,
                initial_partition,
                context.automorphism_algorithm,
                &leaf_candidate,
            )
        }
        DescriptionLevel::Constitution => {
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            canonical_search(
                &adapter,
                &adapter.colors,
                context.automorphism_algorithm,
                &leaf_candidate,
            )
        }
        DescriptionLevel::Structure => {
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
                retain_stereo_preserving_generators(span, &incidence_graph, generators);
            };
            canonical_search_with_generator_filter(
                &adapter,
                &adapter.colors,
                context.automorphism_algorithm,
                &leaf_candidate,
                &filter_generators,
            )
        }
        DescriptionLevel::Full => {
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
                retain_preserving_generators(span, &incidence_graph, generators);
            };
            canonical_search_with_generator_filter(
                &adapter,
                &adapter.colors,
                context.automorphism_algorithm,
                &leaf_candidate,
                &filter_generators,
            )
        }
    };
    let CanonicalCandidate { key, entity_order } = selected.candidate;
    Ok(CanonicalSearchResult {
        candidate: CanonicalCandidate {
            key: key.into_built(),
            entity_order,
        },
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
        let normalized = self.normalize()?;
        let candidate = canonical_candidate(&normalized, level, context)?.candidate;
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
        self.clone()
            .canonicalize(context)
            .map(|canonical| hash_value(&canonical))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        match (
            self.clone().canonicalize(context),
            other.clone().canonicalize(context),
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
        self.clone()
            .canonicalize(context)
            .map(|canonical| hash_value(&canonical))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        match (
            self.clone().canonicalize(context),
            other.clone().canonicalize(context),
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
