//! Molecule canonicalization.

use umol_graph_core::{Graph, NodeId, SubdivisionNodeSource};
use umol_perm::Permutation;

use super::super::aromatic::aromatic_system_representative_action;
use super::super::constraint::Constraints;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::{Entity, EntityKind};
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::StereoLigand;
use super::super::molecule::Molecule;
use super::super::multicenter::multicenter_bond_representative_action;
use super::super::remap::MoleculeRemapping;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{
    stereo_atom_representative_action, stereo_bond_representative_action, StereoAtomForm,
    StereoBondForm,
};
use super::super::traits::{AsColor, FrameTransport, Normalize, Reframe};
use super::{
    canonical_search, canonical_search_from_partition, canonical_search_with_generator_filter,
    hash_value, molecule_canonicalize_level, molecule_counts, molecule_remapping, rank_colors,
    AutomorphismAdapter, AutomorphismColor, BranchOrdering, CanonicalCandidate,
    CanonicalSearchOptions, Canonicalize, CanonicalizeContext, DeferredKey, DescriptionLevel,
    InitialColors, MoleculeCanonicalizeError, OrderedPartition, TopologySearchInput,
};

/// The key of a leaf: the molecule in the entity order of the leaf, with each form as the integer
/// of its color.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct ColorKey {
    /// The integer of each atom.
    atoms: Vec<u32>,
    /// The two endpoint images and the integer of each bond.
    bonds: Vec<[u32; 3]>,
    /// The sorted donor images, the acceptor image, and the integer of each dative bond.
    dative: Vec<(Vec<u32>, u32, u32)>,
    /// The sorted atom images, the electron counts in that order, and the integer of each aromatic
    /// system.
    aromatic: Vec<(Vec<u32>, ElectronCountsForm, u32)>,
    /// The sorted atom images, the electron counts in that order, and the integer of each
    /// multicenter bond.
    multicenter: Vec<(Vec<u32>, ElectronCountsForm, u32)>,
    /// The two endpoint images and the integer of each noncovalent bond.
    noncovalent: Vec<[u32; 3]>,
    /// The site image, the ligands in the frame of the leaf, and the form in that frame of each
    /// stereo atom.
    stereo_atoms: Vec<(u32, Vec<StereoLigand>, StereoAtomForm)>,
    /// The image of the site bond, the ligands in the frame of the leaf, and the form in that frame
    /// of each stereo bond.
    stereo_bonds: Vec<(u32, Vec<StereoLigand>, StereoBondForm)>,
    /// The molecule-level constraints in the entity order and the frames of the leaf.
    constraints: Constraints,
}

fn generator_preserves_stereo(
    molecule: &Molecule,
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

    for source in molecule.stereo_atoms().iter() {
        let Entity::StereoAtom(target_id) = entity_image(Entity::StereoAtom(source.id())) else {
            return false;
        };
        let target = molecule.stereo_atom(target_id);
        if entity_image(Entity::Atom(source.site_id())) != Entity::Atom(target.site_id()) {
            return false;
        }
        let Some(action) = frame_action(source.ligand_ids(), target.ligand_ids()) else {
            return false;
        };
        if !form_preserved(source.attributes(), target.attributes(), action) {
            return false;
        }
    }

    for source in molecule.stereo_bonds().iter() {
        let Entity::StereoBond(target_id) = entity_image(Entity::StereoBond(source.id())) else {
            return false;
        };
        let target = molecule.stereo_bond(target_id);
        if entity_image(Entity::Bond(source.site_id())) != Entity::Bond(target.site_id()) {
            return false;
        }
        let Some(action) = frame_action(source.ligand_ids(), target.ligand_ids()) else {
            return false;
        };
        if !form_preserved(source.attributes(), target.attributes(), action) {
            return false;
        }
    }

    true
}

fn retain_stereo_preserving_generators(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    if molecule.stereo_atoms().count() == 0 && molecule.stereo_bonds().count() == 0 {
        return;
    }
    generators.retain(|generator| generator_preserves_stereo(molecule, incidence_graph, generator));
}

/// Keep the generators that carry `molecule` onto itself: the molecule remapped by the generator
/// equals the molecule under framed equality.
fn retain_preserving_generators(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    generators.retain(|generator| {
        molecule.framed_eq_under(
            molecule,
            &remapping_from_generator(molecule, incidence_graph, generator),
        )
    });
}

/// The integers of the colors of a normalized molecule at `level`: one per entity, by incidence
/// node, and one per incidence, by incidence edge.
///
/// The colors of each entity kind are ranked in kind order, and the incidences continue after the
/// entities. The color of an incidence is its role and, for an atom of an aromatic system or a
/// multicenter bond, the electron count of that atom. At Topology the incidences are not ranked,
/// since the topology search input does not read them.
fn initial_colors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
) -> Result<InitialColors, Contradiction> {
    let mut entities = vec![0; incidence_graph.graph().node_count()];
    let position = |entity: Entity| incidence_graph.node_of(entity).index();

    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| {
            (
                atom.attributes().as_color(),
                position(Entity::Atom(atom.id())),
            )
        })
        .collect();
    let mut next = rank_colors(atoms, 0, &mut entities);
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            (
                bond.attributes().as_color(),
                position(Entity::Bond(bond.id())),
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

    let dative = molecule
        .dative_bonds()
        .iter()
        .map(|bond| {
            (
                bond.attributes().as_color(),
                position(Entity::DativeBond(bond.id())),
            )
        })
        .collect();
    next = rank_colors(dative, next, &mut entities);
    let mut aromatic = Vec::with_capacity(molecule.aromatic_systems().count());
    for system in molecule.aromatic_systems().iter() {
        let attributes = system.attributes();
        if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
            if counts.len() != system.atom_ids().len())
        {
            return Err(Contradiction);
        }
        aromatic.push((
            attributes.as_color(),
            position(Entity::AromaticSystem(system.id())),
        ));
    }
    next = rank_colors(aromatic, next, &mut entities);
    let mut multicenter = Vec::with_capacity(molecule.multicenter_bonds().count());
    for bond in molecule.multicenter_bonds().iter() {
        let attributes = bond.attributes();
        if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
            if counts.len() != bond.atom_ids().len())
        {
            return Err(Contradiction);
        }
        multicenter.push((
            attributes.as_color(),
            position(Entity::MulticenterBond(bond.id())),
        ));
    }
    next = rank_colors(multicenter, next, &mut entities);
    let noncovalent = molecule
        .noncovalent_bonds()
        .iter()
        .map(|bond| {
            (
                bond.attributes().as_color(),
                position(Entity::NoncovalentBond(bond.id())),
            )
        })
        .collect();
    next = rank_colors(noncovalent, next, &mut entities);

    if level >= DescriptionLevel::Structure {
        let mut stereo_atoms = Vec::with_capacity(molecule.stereo_atoms().count());
        for stereo in molecule.stereo_atoms().iter() {
            let attributes = stereo.attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_atom_kind(Entity::StereoAtom(stereo.id()), kind)
                    .map_err(|_| Contradiction)?;
            }
            attributes
                .clone()
                .reframe_by(&Permutation::identity(stereo.ligand_count()))
                .ok_or(Contradiction)?;
            stereo_atoms.push((
                attributes.as_color(),
                position(Entity::StereoAtom(stereo.id())),
            ));
        }
        next = rank_colors(stereo_atoms, next, &mut entities);
        let mut stereo_bonds = Vec::with_capacity(molecule.stereo_bonds().count());
        for stereo in molecule.stereo_bonds().iter() {
            let attributes = stereo.attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_bond_kind(Entity::StereoBond(stereo.id()), kind)
                    .map_err(|_| Contradiction)?;
            }
            attributes
                .clone()
                .reframe_by(&Permutation::identity(stereo.ligand_count()))
                .ok_or(Contradiction)?;
            stereo_bonds.push((
                attributes.as_color(),
                position(Entity::StereoBond(stereo.id())),
            ));
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
                    let (atom_position, electrons) = match owner {
                        Entity::AromaticSystem(id) => (
                            molecule
                                .aromatic_system(id)
                                .atom_ids()
                                .position(|id| id == atom),
                            &molecule.aromatic_system(id).attributes().electrons,
                        ),
                        Entity::MulticenterBond(id) => (
                            molecule
                                .multicenter_bond(id)
                                .atom_ids()
                                .position(|id| id == atom),
                            &molecule.multicenter_bond(id).attributes().electrons,
                        ),
                        _ => unreachable!(
                            "electron contributions belong to aromatic or multicenter entities"
                        ),
                    };
                    let atom_position =
                        atom_position.expect("incidence atom belongs to its overlay");
                    electrons.at(AtomPosition::from(atom_position))
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
fn topology_search_input(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    entity_colors: &[u32],
) -> TopologySearchInput {
    let mut bond_color_counts = vec![0_usize; entity_colors.len()];
    for bond in molecule.bonds().iter() {
        let node = incidence_graph.node_of(Entity::Bond(bond.id()));
        bond_color_counts[entity_colors[node.index()] as usize] += 1;
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

    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let mut colors = entity_colors[..atom_count]
        .iter()
        .copied()
        .map(AutomorphismColor::Entity)
        .collect::<Vec<_>>();
    let mut node_sources = (0..atom_count)
        .map(|index| SubdivisionNodeSource::Node(NodeId(index as u32)))
        .collect::<Vec<_>>();
    let mut edges = Vec::<[u32; 2]>::with_capacity(molecule.bonds().count() * 2);

    for bond in molecule.bonds().iter() {
        let node = incidence_graph.node_of(Entity::Bond(bond.id()));
        let bond_color = entity_colors[node.index()];
        let [first, second] = bond.atom_ids().map(|atom| atom.index() as u32);
        if direct_bond_color == Some(bond_color) {
            edges.push([first, second]);
        } else {
            let subdivision = NodeId(node_sources.len() as u32);
            node_sources.push(SubdivisionNodeSource::Node(subdivision));
            colors.push(AutomorphismColor::Entity(bond_color));
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

fn remapping_from_order(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> MoleculeRemapping {
    let mut images = molecule_counts(molecule).map(|count| (0..count).collect::<Vec<_>>());
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
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generator: &[NodeId],
) -> MoleculeRemapping {
    let mut images = molecule_counts(molecule).map(|count| (0..count).collect::<Vec<_>>());
    for node in incidence_graph.graph().node_ids() {
        let entity = incidence_graph.entity(node);
        images[entity.kind() as usize][entity.kind_id()] =
            incidence_graph.entity(generator[node.index()]).kind_id();
    }
    molecule_remapping(&images)
}

/// The order of every entity of `level` for an atom order.
///
/// The atoms stand as the search orders them. The entities of each further kind follow by their
/// participants. Entities of one kind on the same participants keep their search order.
fn entity_order(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    order: &[NodeId],
) -> Vec<NodeId> {
    let node_count = incidence_graph.graph().node_count();
    let mut positions = vec![0_usize; node_count];
    let mut atom_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Atom)];
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

    let mut bond_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Bond)];
    let mut bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let [first, second] = bond.atom_ids().map(atom_image);
            let (position, node) = search_position(Entity::Bond(bond.id()));
            (
                (first.min(second), first.max(second)),
                position,
                node,
                bond.id(),
            )
        })
        .collect::<Vec<_>>();
    bonds.sort_unstable();
    for (image, &(_, _, node, id)) in bonds.iter().enumerate() {
        bond_images[id.index()] = image;
        entity_order.push(node);
    }

    if level >= DescriptionLevel::Constitution {
        let mut dative = molecule
            .dative_bonds()
            .iter()
            .map(|bond| {
                let mut donors = bond.donor_ids().map(atom_image).collect::<Vec<_>>();
                donors.sort_unstable();
                let (position, node) = search_position(Entity::DativeBond(bond.id()));
                ((donors, atom_image(bond.acceptor_id())), position, node)
            })
            .collect::<Vec<_>>();
        dative.sort_unstable();
        entity_order.extend(dative.into_iter().map(|(_, _, node)| node));

        let mut aromatic = molecule
            .aromatic_systems()
            .iter()
            .map(|system| {
                let mut atoms = system.atom_ids().map(atom_image).collect::<Vec<_>>();
                atoms.sort_unstable();
                let (position, node) = search_position(Entity::AromaticSystem(system.id()));
                (atoms, position, node)
            })
            .collect::<Vec<_>>();
        aromatic.sort_unstable();
        entity_order.extend(aromatic.into_iter().map(|(_, _, node)| node));

        let mut multicenter = molecule
            .multicenter_bonds()
            .iter()
            .map(|bond| {
                let mut atoms = bond.atom_ids().map(atom_image).collect::<Vec<_>>();
                atoms.sort_unstable();
                let (position, node) = search_position(Entity::MulticenterBond(bond.id()));
                (atoms, position, node)
            })
            .collect::<Vec<_>>();
        multicenter.sort_unstable();
        entity_order.extend(multicenter.into_iter().map(|(_, _, node)| node));

        let mut noncovalent = molecule
            .noncovalent_bonds()
            .iter()
            .map(|bond| {
                let [first, second] = bond.atom_ids().map(atom_image);
                let (position, node) = search_position(Entity::NoncovalentBond(bond.id()));
                ((first.min(second), first.max(second)), position, node)
            })
            .collect::<Vec<_>>();
        noncovalent.sort_unstable();
        entity_order.extend(noncovalent.into_iter().map(|(_, _, node)| node));
    }

    if level >= DescriptionLevel::Structure {
        let mut stereo_atoms = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| {
                let (position, node) = search_position(Entity::StereoAtom(stereo.id()));
                (atom_image(stereo.site_id()), position, node)
            })
            .collect::<Vec<_>>();
        stereo_atoms.sort_unstable();
        entity_order.extend(stereo_atoms.into_iter().map(|(_, _, node)| node));

        let mut stereo_bonds = molecule
            .stereo_bonds()
            .iter()
            .map(|stereo| {
                let (position, node) = search_position(Entity::StereoBond(stereo.id()));
                (bond_images[stereo.site_id().index()], position, node)
            })
            .collect::<Vec<_>>();
        stereo_bonds.sort_unstable();
        entity_order.extend(stereo_bonds.into_iter().map(|(_, _, node)| node));
    }

    entity_order
}

/// The key of the leaf that orders the entities of `molecule` as `entity_order`.
///
/// A row holds the participants of its entity in their new numbers, in the representative order of
/// the entity, and the integer of its color. The row of an aromatic system or a multicenter bond
/// also holds its electron counts in that order, and the row of a stereo entity holds its form in
/// the frame of the leaf in place of the integer. The molecule-level constraints follow, read
/// from the molecule remapped and reframed in the order of the leaf.
fn color_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    colors: &InitialColors,
    entity_order: &[NodeId],
) -> ColorKey {
    let mut atom_images = vec![0_u32; incidence_graph.entity_count(EntityKind::Atom)];
    let mut bond_images = vec![0_u32; incidence_graph.entity_count(EntityKind::Bond)];
    let mut key = ColorKey {
        atoms: Vec::with_capacity(atom_images.len()),
        bonds: Vec::with_capacity(bond_images.len()),
        dative: Vec::new(),
        aromatic: Vec::new(),
        multicenter: Vec::new(),
        noncovalent: Vec::new(),
        stereo_atoms: Vec::new(),
        stereo_bonds: Vec::new(),
        constraints: Constraints::new(),
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
                let [first, second] = molecule
                    .bond(id)
                    .atom_ids()
                    .map(|atom| atom_images[atom.index()]);
                key.bonds
                    .push([first.min(second), first.max(second), color]);
            }
            Entity::DativeBond(id) => {
                let bond = molecule.dative_bond(id);
                let mut donors = bond
                    .donor_ids()
                    .map(|atom| atom_images[atom.index()])
                    .collect::<Vec<_>>();
                donors.sort_unstable();
                key.dative
                    .push((donors, atom_images[bond.acceptor_id().index()], color));
            }
            Entity::AromaticSystem(id) => {
                let system = molecule.aromatic_system(id);
                let atoms = system
                    .atom_ids()
                    .map(|atom| AtomId(atom_images[atom.index()]))
                    .collect::<Vec<_>>();
                let action = aromatic_system_representative_action(atoms.clone());
                key.aromatic.push((
                    action
                        .act(&atoms)
                        .expect("the action has the degree of the atom list")
                        .into_iter()
                        .map(|atom| atom.0)
                        .collect(),
                    system
                        .attributes()
                        .electrons
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established the length of the electron counts"),
                    color,
                ));
            }
            Entity::MulticenterBond(id) => {
                let bond = molecule.multicenter_bond(id);
                let atoms = bond
                    .atom_ids()
                    .map(|atom| AtomId(atom_images[atom.index()]))
                    .collect::<Vec<_>>();
                let action = multicenter_bond_representative_action(atoms.clone());
                key.multicenter.push((
                    action
                        .act(&atoms)
                        .expect("the action has the degree of the atom list")
                        .into_iter()
                        .map(|atom| atom.0)
                        .collect(),
                    bond.attributes()
                        .electrons
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established the length of the electron counts"),
                    color,
                ));
            }
            Entity::NoncovalentBond(id) => {
                let [first, second] = molecule
                    .noncovalent_bond(id)
                    .atom_ids()
                    .map(|atom| atom_images[atom.index()]);
                key.noncovalent
                    .push([first.min(second), first.max(second), color]);
            }
            Entity::StereoAtom(id) => {
                let stereo = molecule.stereo_atom(id);
                let ligands = stereo
                    .ligand_ids()
                    .iter()
                    .map(|ligand| {
                        StereoLigand::new(AtomId(atom_images[ligand.atom_id.index()]), ligand.kind)
                    })
                    .collect::<Vec<_>>();
                let action = stereo_atom_representative_action(&ligands)
                    .expect("stereo atom ligand frame has a supported degree");
                key.stereo_atoms.push((
                    atom_images[stereo.site_id().index()],
                    action.act(&ligands),
                    stereo
                        .attributes()
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established compatible stereo frame actions")
                        .normalize()
                        .expect("the molecule is normalized"),
                ));
            }
            Entity::StereoBond(id) => {
                let stereo = molecule.stereo_bond(id);
                let ligands = stereo
                    .ligand_ids()
                    .iter()
                    .map(|ligand| {
                        StereoLigand::new(AtomId(atom_images[ligand.atom_id.index()]), ligand.kind)
                    })
                    .collect::<Vec<_>>();
                let action = stereo_bond_representative_action(&ligands)
                    .expect("stereo bond ligand frame has a supported degree");
                key.stereo_bonds.push((
                    bond_images[stereo.site_id().index()],
                    action.act(&ligands),
                    stereo
                        .attributes()
                        .clone()
                        .reframe_by(&action)
                        .expect("the colors established compatible stereo frame actions")
                        .normalize()
                        .expect("the molecule is normalized"),
                ));
            }
        }
    }
    if molecule.has_constraints() {
        key.constraints = molecule
            .clone()
            .remap(&remapping_from_order(
                molecule,
                incidence_graph,
                entity_order,
            ))
            .reframe()
            .expect("the molecule is normalized")
            .constraints()
            .clone();
    }
    key
}

fn canonicalize_topology(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let molecule = molecule.normalize()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Topology)?;
    let TopologySearchInput {
        adapter,
        initial_partition,
    } = topology_search_input(&molecule, &incidence_graph, &colors.entities);
    let build_key = |order: &[NodeId]| color_key(&molecule, &incidence_graph, &colors, order);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Topology,
            order,
        );
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let entity_order = canonical_search_from_partition(
        &adapter,
        initial_partition,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
    )
    .candidate
    .entity_order;
    let remapping = remapping_from_order(&molecule, &incidence_graph, &entity_order);

    Ok((molecule.remap(&remapping), remapping))
}

fn canonicalize_constitution(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let molecule = molecule.normalize()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Constitution)?;
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let build_key = |order: &[NodeId]| color_key(&molecule, &incidence_graph, &colors, order);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Constitution,
            order,
        );
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let entity_order = canonical_search(
        &adapter,
        &adapter.colors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
    )
    .candidate
    .entity_order;
    let remapping = remapping_from_order(&molecule, &incidence_graph, &entity_order);

    let canonical = molecule.remap(&remapping).reframe()?;
    Ok((canonical, remapping))
}

fn canonicalize_structure(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let molecule = molecule.normalize()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure)?;
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let build_key = |order: &[NodeId]| color_key(&molecule, &incidence_graph, &colors, order);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        );
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let entity_order = canonical_search_with_generator_filter(
        &adapter,
        &adapter.colors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &filter_generators,
    )
    .candidate
    .entity_order;
    let remapping = remapping_from_order(&molecule, &incidence_graph, &entity_order);

    let canonical = molecule.remap(&remapping).reframe()?;
    Ok((canonical, remapping))
}

fn canonicalize_full(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let molecule = molecule.normalize()?;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Full)?;
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let build_key = |order: &[NodeId]| color_key(&molecule, &incidence_graph, &colors, order);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(&molecule, &incidence_graph, DescriptionLevel::Full, order);
        CanonicalCandidate {
            key: DeferredKey::new(&entity_order, &build_key),
            entity_order,
        }
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let entity_order = canonical_search_with_generator_filter(
        &adapter,
        &adapter.colors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &filter_generators,
    )
    .candidate
    .entity_order;
    let remapping = remapping_from_order(&molecule, &incidence_graph, &entity_order);
    let canonical = molecule.remap(&remapping).reframe()?;
    Ok((canonical, remapping))
}

impl Molecule {
    /// The canonical form at `level` and its source-to-canonical remapping.
    fn tracked_canonicalize_at(
        self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), MoleculeCanonicalizeError> {
        let options = CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        };
        match level {
            DescriptionLevel::Topology => canonicalize_topology(self, context, options),
            DescriptionLevel::Constitution => canonicalize_constitution(self, context, options),
            DescriptionLevel::Structure => canonicalize_structure(self, context, options),
            DescriptionLevel::Full => canonicalize_full(self, context, options),
        }
    }
}

impl Canonicalize for Molecule {
    type Error = MoleculeCanonicalizeError;

    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error> {
        let level = molecule_canonicalize_level(&self);
        self.tracked_canonicalize_at(level, context)
            .map(|(canonical, _)| canonical)
    }

    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error> {
        let level = molecule_canonicalize_level(&self);
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
                Err(MoleculeCanonicalizeError::Contradiction(_)),
                Err(MoleculeCanonicalizeError::Contradiction(_)),
            ) => true,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests;
