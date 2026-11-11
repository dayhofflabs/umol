//! Molecule canonicalization.

use umol_graph_core::{Graph, NodeId, SubdivisionNodeSource};
use umol_perm::Permutation;

use super::super::electrons::ElectronCountsForm;
use super::super::entity::{Entity, EntityKind};
use super::super::error::Contradiction;
use super::super::id::{AtomId, BondId};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::remap::MoleculeRemapping;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::StereoConfigurationForm;
use super::super::traits::{FrameTransport, Normalize, Reframe};
use super::{
    aromatic_system_constraint_form_key, atom_constraint_form_key, atom_inherent_fields,
    bond_constraint_form_key, bond_inherent_fields, canonical_kinded_stereo_frame,
    canonical_search, canonical_search_from_partition, canonical_search_with_generator_filter,
    constitution_partition_descriptors, constraint_key, dative_bond_constraint_form_key,
    electron_counts_form_key, field, hash_value, incidence_key, index_key,
    molecule_canonicalize_level, molecule_counts, molecule_remapping,
    multicenter_bond_constraint_form_key, noncovalent_bond_constraint_form_key,
    noncovalent_bond_kind_form_key, num_form_key, option, positioned_product, product,
    rank_initial_colors, sequence, stereo_atom_constraint_form_key,
    stereo_bond_constraint_form_key, stereo_configuration_form_key, stereo_frame_permutations,
    stereo_kind_key, stereo_ligand_key, unpaired_electrons_form_key, variant, AutomorphismAdapter,
    AutomorphismColor, BranchOrdering, CanonicalCandidate, CanonicalComparisonKey,
    CanonicalKeyValue, CanonicalSearchOptions, Canonicalize, CanonicalizeContext,
    ConstraintBlockKey, ConstraintBlockPosition, DescriptionLevel, EntityBlockPosition, FieldKey,
    InitialColorKey, MoleculeCanonicalizeError, OrderedPartition, PositionedKey,
    TopologySearchInput,
};

/// Construct the ordered constraint blocks used by complete canonicalization.
fn constraint_blocks(molecule: &Molecule) -> Vec<ConstraintBlockKey> {
    let mut blocks = Vec::new();
    macro_rules! inline_block {
        ($position:expr, $entities:expr, $key:expr) => {{
            let rows = $entities
                .iter()
                .flat_map(|entity| {
                    entity
                        .attributes()
                        .constraints
                        .iter()
                        .map(move |constraint| {
                            product([index_key(entity.id().index()), $key(constraint)])
                        })
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

    let rows = molecule
        .atoms()
        .iter()
        .flat_map(|atom| {
            atom.attributes().constraints.iter().map(move |constraint| {
                product([
                    index_key(atom.id().index()),
                    atom_constraint_form_key(constraint),
                ])
            })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::ATOM,
            value: sequence(rows),
        });
    }
    let rows = molecule
        .bonds()
        .iter()
        .flat_map(|bond| {
            bond.attributes().constraints.iter().map(move |constraint| {
                product([
                    index_key(bond.id().index()),
                    bond_constraint_form_key(constraint),
                ])
            })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::BOND,
            value: sequence(rows),
        });
    }
    let rows = molecule
        .dative_bonds()
        .iter()
        .flat_map(|bond| {
            bond.attributes().constraints.iter().map(move |constraint| {
                product([
                    index_key(bond.id().index()),
                    dative_bond_constraint_form_key(constraint),
                ])
            })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::DATIVE_BOND,
            value: sequence(rows),
        });
    }
    let rows = molecule
        .aromatic_systems()
        .iter()
        .flat_map(|system| {
            system
                .attributes()
                .constraints
                .iter()
                .map(move |constraint| {
                    product([
                        index_key(system.id().index()),
                        aromatic_system_constraint_form_key(constraint),
                    ])
                })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::AROMATIC_SYSTEM,
            value: sequence(rows),
        });
    }
    let rows = molecule
        .multicenter_bonds()
        .iter()
        .flat_map(|bond| {
            bond.attributes().constraints.iter().map(move |constraint| {
                product([
                    index_key(bond.id().index()),
                    multicenter_bond_constraint_form_key(constraint),
                ])
            })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::MULTICENTER_BOND,
            value: sequence(rows),
        });
    }
    let rows = molecule
        .noncovalent_bonds()
        .iter()
        .flat_map(|bond| {
            bond.attributes().constraints.iter().map(move |constraint| {
                product([
                    index_key(bond.id().index()),
                    noncovalent_bond_constraint_form_key(constraint),
                ])
            })
        })
        .collect::<Vec<_>>();
    if !rows.is_empty() {
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::NONCOVALENT_BOND,
            value: sequence(rows),
        });
    }
    inline_block!(
        ConstraintBlockPosition::STEREO_ATOM,
        molecule.stereo_atoms(),
        stereo_atom_constraint_form_key
    );
    inline_block!(
        ConstraintBlockPosition::STEREO_BOND,
        molecule.stereo_bonds(),
        stereo_bond_constraint_form_key
    );

    if !molecule.constraints().is_empty() {
        let mut constraints = molecule
            .constraints()
            .iter()
            .map(constraint_key)
            .collect::<Vec<_>>();
        constraints.sort();
        constraints.dedup();
        blocks.push(PositionedKey {
            position: ConstraintBlockPosition::MOLECULE,
            value: sequence(constraints),
        });
    }
    blocks
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
    let configuration_preserved = |source: &StereoConfigurationForm,
                                   target: &StereoConfigurationForm,
                                   action: Permutation| {
        source
            .clone()
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(target))
    };

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
        if !configuration_preserved(
            &source.attributes().configuration,
            &target.attributes().configuration,
            action,
        ) {
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
        if !configuration_preserved(
            &source.attributes().configuration,
            &target.attributes().configuration,
            action,
        ) {
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

/// Construct the normalized entity and incidence keys used for initial colors.
fn initial_color_keys(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<(Vec<InitialColorKey>, Vec<InitialColorKey>), Contradiction> {
    let entity_keys = incidence_graph
        .graph()
        .node_ids()
        .map(|node| entity_color_key(molecule, incidence_graph.entity(node)))
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
                    let (tag, position, electrons) = match owner {
                        Entity::AromaticSystem(id) => (
                            3,
                            molecule
                                .aromatic_system(id)
                                .atom_ids()
                                .position(|id| id == atom),
                            &molecule.aromatic_system(id).attributes().electrons,
                        ),
                        Entity::MulticenterBond(id) => (
                            4,
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
                    let position = position.expect("incidence atom belongs to its overlay");
                    let contribution = match electrons {
                        ElectronCountsForm::Undetermined => NumForm::Undetermined,
                        ElectronCountsForm::Lit(counts) => NumForm::Lit(counts[position]),
                    };
                    variant(tag, [num_form_key(&contribution)])
                }
                incidence => incidence_key(incidence),
            };
            InitialColorKey::Incidence(value)
        })
        .collect();

    Ok((entity_keys, incidence_keys))
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
            incidence_nodes: Vec::new(),
            source_node_count: atom_count,
        },
        initial_partition,
    }
}

fn constitution_entity_classes(molecule: &Molecule) -> Result<Vec<u32>, Contradiction> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let partition = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
    let classes = partition.cell_indices(adapter.graph().node_count());
    Ok(classes[..adapter.source_node_count].to_vec())
}

fn stereo_refinement_descriptor(
    site_class: u32,
    ligand_classes: &[(u32, StereoLigandKind)],
    configuration: &StereoConfigurationForm,
) -> Result<CanonicalKeyValue, Contradiction> {
    let descriptor = |ligands: Vec<(u32, StereoLigandKind)>,
                      configuration: StereoConfigurationForm| {
        product([
            CanonicalKeyValue::Unsigned(site_class.into()),
            sequence(
                ligands
                    .into_iter()
                    .map(|(class, kind)| stereo_ligand_key(class, kind)),
            ),
            stereo_configuration_form_key(&configuration),
        ])
    };

    match configuration {
        StereoConfigurationForm::Undetermined => {
            let mut ligands = ligand_classes.to_vec();
            ligands.sort_unstable();
            Ok(descriptor(ligands, StereoConfigurationForm::Undetermined))
        }
        StereoConfigurationForm::Kinded(kind, _) => {
            if ligand_classes.len() != kind.degree() {
                return Err(Contradiction);
            }
            stereo_frame_permutations(*kind)
                .map(|permutation| {
                    configuration
                        .apply(permutation)
                        .ok_or(Contradiction)?
                        .normalize()
                        .map(|configuration| {
                            descriptor(permutation.act(ligand_classes), configuration)
                        })
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .min()
                .ok_or(Contradiction)
        }
    }
}

fn structure_partition_descriptors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    entity_classes: &[u32],
) -> Result<Vec<InitialColorKey>, Contradiction> {
    let entity_class = |entity: Entity| entity_classes[incidence_graph.node_of(entity).index()];

    adapter
        .node_sources
        .iter()
        .map(|source| match *source {
            SubdivisionNodeSource::Node(node) => {
                let entity = incidence_graph.entity(node);
                let InitialColorKey::Entity { position, .. } = entity_keys[node.index()] else {
                    unreachable!("entity key corresponds to an entity node")
                };
                let value = match entity {
                    Entity::StereoAtom(id) => {
                        let stereo = molecule
                            .stereo_atoms()
                            .get(id)
                            .expect("structure incidence stereo atom is in range");
                        let ligand_classes = stereo
                            .ligands()
                            .map(|ligand| {
                                (entity_class(Entity::Atom(ligand.atom_id())), ligand.kind())
                            })
                            .collect::<Vec<_>>();
                        stereo_refinement_descriptor(
                            entity_class(Entity::Atom(stereo.site_id())),
                            &ligand_classes,
                            &stereo.attributes().configuration,
                        )?
                    }
                    Entity::StereoBond(id) => {
                        let stereo = molecule
                            .stereo_bonds()
                            .get(id)
                            .expect("structure incidence stereo bond is in range");
                        let ligand_classes = stereo
                            .ligands()
                            .map(|ligand| {
                                (entity_class(Entity::Atom(ligand.atom_id())), ligand.kind())
                            })
                            .collect::<Vec<_>>();
                        stereo_refinement_descriptor(
                            entity_class(Entity::Bond(stereo.site_id())),
                            &ligand_classes,
                            &stereo.attributes().configuration,
                        )?
                    }
                    _ => product([CanonicalKeyValue::Unsigned(entity_class(entity).into())]),
                };
                Ok(InitialColorKey::Entity { position, value })
            }
            SubdivisionNodeSource::Edge(edge) => {
                let value = match incidence_graph.incidence(edge) {
                    Incidence::DativeDonor | Incidence::DativeAcceptor => variant(1, []),
                    incidence => incidence_key(incidence),
                };
                Ok(InitialColorKey::Incidence(value))
            }
        })
        .collect()
}

fn para_stereo_partition_descriptors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    partition: &OrderedPartition,
) -> Result<Vec<(u32, Option<CanonicalKeyValue>)>, Contradiction> {
    let classes = partition.cell_indices(adapter.graph().node_count());
    let entity_class = |entity: Entity| classes[incidence_graph.node_of(entity).index()];
    let ligand_classes = |carrier: Entity| {
        let carrier_node = incidence_graph.node_of(carrier);
        let mut edges = incidence_graph
            .graph()
            .neighbors(carrier_node)
            .iter()
            .map(|neighbor| neighbor.edge)
            .collect::<Vec<_>>();
        edges.sort_unstable();
        edges
            .into_iter()
            .filter_map(|edge| {
                let Incidence::StereoLigand(kind) = incidence_graph.incidence(edge) else {
                    return None;
                };
                let occurrence = adapter
                    .node_of(SubdivisionNodeSource::Edge(edge))
                    .expect("stereo ligand incidence has an occurrence node");
                Some((classes[occurrence.index()], *kind))
            })
            .collect::<Vec<_>>()
    };

    adapter
        .node_sources
        .iter()
        .enumerate()
        .map(|(index, source)| {
            let descriptor = match *source {
                SubdivisionNodeSource::Node(node) => match incidence_graph.entity(node) {
                    Entity::StereoAtom(id) => {
                        let stereo = molecule
                            .stereo_atoms()
                            .get(id)
                            .expect("structure incidence stereo atom is in range");
                        Some(stereo_refinement_descriptor(
                            entity_class(Entity::Atom(stereo.site_id())),
                            &ligand_classes(Entity::StereoAtom(id)),
                            &stereo.attributes().configuration,
                        )?)
                    }
                    Entity::StereoBond(id) => {
                        let stereo = molecule
                            .stereo_bonds()
                            .get(id)
                            .expect("structure incidence stereo bond is in range");
                        Some(stereo_refinement_descriptor(
                            entity_class(Entity::Bond(stereo.site_id())),
                            &ligand_classes(Entity::StereoBond(id)),
                            &stereo.attributes().configuration,
                        )?)
                    }
                    _ => None,
                },
                SubdivisionNodeSource::Edge(_) => None,
            };
            Ok((classes[index], descriptor))
        })
        .collect()
}

/// Refine the structure partition, including para-stereo rounds when requested.
fn structure_partition(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    adapter: &AutomorphismAdapter,
    entity_keys: &[InitialColorKey],
    para_stereo: bool,
) -> Result<(OrderedPartition, usize), Contradiction> {
    let constitution_classes = constitution_entity_classes(molecule)?;
    let descriptors = structure_partition_descriptors(
        molecule,
        incidence_graph,
        adapter,
        entity_keys,
        &constitution_classes,
    )?;
    let mut partition = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
    let mut rounds = 1;
    let has_unresolved_stereo = |partition: &OrderedPartition| {
        partition.cells.iter().any(|cell| {
            cell.len() > 1
                && cell.iter().any(|&node| {
                    let SubdivisionNodeSource::Node(node) = adapter.node_source(node) else {
                        return false;
                    };
                    matches!(
                        incidence_graph.entity(node),
                        Entity::StereoAtom(_) | Entity::StereoBond(_)
                    )
                })
        })
    };

    if !para_stereo || !has_unresolved_stereo(&partition) {
        return Ok((partition, rounds));
    }

    loop {
        let cell_count = partition.cells.len();
        let descriptors =
            para_stereo_partition_descriptors(molecule, incidence_graph, adapter, &partition)?;
        let next = OrderedPartition::from_descriptors(&descriptors).refine(adapter.graph());
        rounds += 1;
        if next == partition {
            return Ok((partition, rounds));
        }
        debug_assert!(next.cells.len() > cell_count);
        partition = next;
        if !has_unresolved_stereo(&partition) {
            return Ok((partition, rounds));
        }
    }
}

fn entity_color_key(molecule: &Molecule, entity: Entity) -> Result<InitialColorKey, Contradiction> {
    let (position, value) = match entity {
        Entity::Atom(id) => {
            let attributes = molecule.atom(id).attributes();
            (
                EntityBlockPosition::ATOM,
                CanonicalKeyValue::Product(atom_inherent_fields(attributes)?),
            )
        }
        Entity::Bond(id) => {
            let attributes = molecule.bond(id).attributes();
            (
                EntityBlockPosition::BOND,
                CanonicalKeyValue::Product(bond_inherent_fields(attributes)?),
            )
        }
        Entity::DativeBond(id) => {
            let attributes = molecule
                .dative_bonds()
                .get(id)
                .expect("incidence dative bond is in range")
                .attributes();
            (
                EntityBlockPosition::DATIVE_BOND,
                positioned_product([(2, num_form_key(&attributes.order.clone().normalize()?))]),
            )
        }
        Entity::AromaticSystem(id) => {
            let attributes = molecule
                .aromatic_systems()
                .get(id)
                .expect("incidence aromatic system is in range")
                .attributes();
            if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                if counts.len() != molecule.aromatic_system(id).atom_ids().len())
            {
                return Err(Contradiction);
            }
            (
                EntityBlockPosition::AROMATIC_SYSTEM,
                positioned_product([
                    (2, num_form_key(&attributes.charge.clone().normalize()?)),
                    (
                        3,
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ),
                ]),
            )
        }
        Entity::MulticenterBond(id) => {
            let attributes = molecule
                .multicenter_bonds()
                .get(id)
                .expect("incidence multicenter bond is in range")
                .attributes();
            if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                if counts.len() != molecule.multicenter_bond(id).atom_ids().len())
            {
                return Err(Contradiction);
            }
            (
                EntityBlockPosition::MULTICENTER_BOND,
                positioned_product([
                    (2, num_form_key(&attributes.charge.clone().normalize()?)),
                    (
                        3,
                        unpaired_electrons_form_key(
                            &attributes.unpaired_electrons.clone().normalize()?,
                        ),
                    ),
                ]),
            )
        }
        Entity::NoncovalentBond(id) => {
            let attributes = molecule
                .noncovalent_bonds()
                .get(id)
                .expect("incidence noncovalent bond is in range")
                .attributes();
            (
                EntityBlockPosition::NONCOVALENT_BOND,
                positioned_product([(
                    1,
                    noncovalent_bond_kind_form_key(&attributes.kind.clone().normalize()?),
                )]),
            )
        }
        Entity::StereoAtom(id) => {
            let attributes = molecule
                .stereo_atoms()
                .get(id)
                .expect("incidence stereo atom is in range")
                .attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_atom_kind(Entity::StereoAtom(id), kind).map_err(|_| Contradiction)?;
            }
            attributes
                .clone()
                .reframe_by(&Permutation::identity(
                    molecule.stereo_atom(id).ligand_count(),
                ))
                .ok_or(Contradiction)?;
            (
                EntityBlockPosition::STEREO_ATOM,
                positioned_product([(
                    2,
                    option(
                        attributes
                            .configuration
                            .clone()
                            .normalize()?
                            .kind()
                            .map(stereo_kind_key),
                    ),
                )]),
            )
        }
        Entity::StereoBond(id) => {
            let attributes = molecule
                .stereo_bonds()
                .get(id)
                .expect("incidence stereo bond is in range")
                .attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_bond_kind(Entity::StereoBond(id), kind).map_err(|_| Contradiction)?;
            }
            attributes
                .clone()
                .reframe_by(&Permutation::identity(
                    molecule.stereo_bond(id).ligand_count(),
                ))
                .ok_or(Contradiction)?;
            (
                EntityBlockPosition::STEREO_BOND,
                positioned_product([(
                    2,
                    option(
                        attributes
                            .configuration
                            .clone()
                            .normalize()?
                            .kind()
                            .map(stereo_kind_key),
                    ),
                )]),
            )
        }
    };

    Ok(InitialColorKey::Entity { position, value })
}

fn topology_bond_key_rows(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    atom_images: &[usize],
) -> Result<Vec<(CanonicalKeyValue, BondId, NodeId)>, Contradiction> {
    let bond_count = incidence_graph.entity_count(EntityKind::Bond);
    let mut bonds = molecule
        .bonds()
        .iter()
        .map(|id| {
            let bond = id;
            let mut fields = Vec::with_capacity(4);
            let [first, second] = bond.atom_ids().map(|atom| atom_images[atom.index()] as u64);
            fields.push(field(
                0,
                product([
                    CanonicalKeyValue::Unsigned(first.min(second)),
                    CanonicalKeyValue::Unsigned(first.max(second)),
                ]),
            ));
            fields.extend(bond_inherent_fields(bond.attributes())?);
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id(),
                incidence_graph.node_of(Entity::Bond(bond.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    debug_assert_eq!(bonds.len(), bond_count);
    bonds.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    Ok(bonds)
}

fn topology_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let mut atom_images = vec![0_usize; atom_count];
    let mut atom_order = Vec::with_capacity(atom_count);

    for &node in order {
        if let Entity::Atom(id) = incidence_graph.entity(node) {
            atom_images[id.index()] = atom_order.len();
            atom_order.push(id);
        }
    }

    debug_assert_eq!(atom_order.len(), atom_count);

    let atoms = atom_order
        .iter()
        .copied()
        .map(|id| {
            atom_inherent_fields(molecule.atom(id).attributes()).map(CanonicalKeyValue::Product)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = topology_bond_key_rows(molecule, incidence_graph, &atom_images)?;

    let mut entity_order = atom_order
        .into_iter()
        .map(|id| incidence_graph.node_of(Entity::Atom(id)))
        .collect::<Vec<_>>();
    entity_order.extend(bonds.iter().map(|(_, _, node)| *node));
    let bonds = bonds.into_iter().map(|(key, _, _)| key).collect::<Vec<_>>();

    let mut entity_blocks = Vec::with_capacity(2);
    if !atoms.is_empty() {
        entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::ATOM,
            value: CanonicalKeyValue::Sequence(atoms),
        });
    }
    if !bonds.is_empty() {
        entity_blocks.push(PositionedKey {
            position: EntityBlockPosition::BOND,
            value: CanonicalKeyValue::Sequence(bonds),
        });
    }

    Ok(CanonicalCandidate {
        key: CanonicalComparisonKey {
            entity_blocks,
            constraints: Vec::new(),
        },
        entity_order,
    })
}

fn constitution_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    fn electron_occurrence_fields(
        atom_ids: impl IntoIterator<Item = AtomId>,
        electrons: &ElectronCountsForm,
        atom_images: &[usize],
    ) -> [FieldKey; 2] {
        match electrons {
            ElectronCountsForm::Undetermined => {
                let mut participants = atom_ids
                    .into_iter()
                    .map(|atom| atom_images[atom.index()] as u64)
                    .collect::<Vec<_>>();
                participants.sort_unstable();
                [
                    field(
                        0,
                        sequence(participants.into_iter().map(CanonicalKeyValue::Unsigned)),
                    ),
                    field(1, electron_counts_form_key(electrons)),
                ]
            }
            ElectronCountsForm::Lit(electrons) => {
                let mut occurrences = atom_ids
                    .into_iter()
                    .zip(electrons)
                    .map(|(atom, &electrons)| (atom_images[atom.index()] as u64, electrons))
                    .collect::<Vec<_>>();
                occurrences.sort_unstable_by_key(|&(atom, _)| atom);
                let (participants, electrons): (Vec<_>, Vec<_>) = occurrences.into_iter().unzip();
                [
                    field(
                        0,
                        sequence(participants.into_iter().map(CanonicalKeyValue::Unsigned)),
                    ),
                    field(
                        1,
                        electron_counts_form_key(&ElectronCountsForm::Lit(electrons)),
                    ),
                ]
            }
        }
    }

    let mut candidate = topology_candidate(molecule, incidence_graph, order)?;
    let atom_count = incidence_graph.entity_count(EntityKind::Atom);
    let mut atom_images = vec![0_usize; atom_count];
    for (image, &node) in candidate.entity_order[..atom_count].iter().enumerate() {
        let Entity::Atom(id) = incidence_graph.entity(node) else {
            unreachable!("topology candidate begins with every atom")
        };
        atom_images[id.index()] = image;
    }

    let mut dative = molecule
        .dative_bonds()
        .iter()
        .map(|bond| {
            let mut donors = bond
                .donor_ids()
                .map(|atom| atom_images[atom.index()] as u64)
                .collect::<Vec<_>>();
            donors.sort_unstable();
            let fields = vec![
                field(
                    0,
                    sequence(donors.into_iter().map(CanonicalKeyValue::Unsigned)),
                ),
                field(
                    1,
                    CanonicalKeyValue::Unsigned(atom_images[bond.acceptor_id().index()] as u64),
                ),
                field(2, num_form_key(&bond.order().clone().normalize()?)),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id(),
                incidence_graph.node_of(Entity::DativeBond(bond.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    dative.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let dative = dative
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut aromatic = molecule
        .aromatic_systems()
        .iter()
        .map(|system| {
            let electrons = system.electrons().clone().normalize()?;
            let mut fields =
                electron_occurrence_fields(system.atom_ids(), &electrons, &atom_images)
                    .into_iter()
                    .collect::<Vec<_>>();
            fields.extend([
                field(2, num_form_key(&system.charge().clone().normalize()?)),
                field(
                    3,
                    unpaired_electrons_form_key(&system.unpaired_electrons().clone().normalize()?),
                ),
            ]);
            Ok((
                CanonicalKeyValue::Product(fields),
                system.id(),
                incidence_graph.node_of(Entity::AromaticSystem(system.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    aromatic.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let aromatic = aromatic
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut multicenter = molecule
        .multicenter_bonds()
        .iter()
        .map(|bond| {
            let electrons = bond.electrons().clone().normalize()?;
            let mut fields = electron_occurrence_fields(bond.atom_ids(), &electrons, &atom_images)
                .into_iter()
                .collect::<Vec<_>>();
            fields.extend([
                field(2, num_form_key(&bond.charge().clone().normalize()?)),
                field(
                    3,
                    unpaired_electrons_form_key(&bond.unpaired_electrons().clone().normalize()?),
                ),
            ]);
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id(),
                incidence_graph.node_of(Entity::MulticenterBond(bond.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    multicenter.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let multicenter = multicenter
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut noncovalent = molecule
        .noncovalent_bonds()
        .iter()
        .map(|bond| {
            let [first, second] = bond.atom_ids().map(|atom| atom_images[atom.index()] as u64);
            let fields = vec![
                field(
                    0,
                    product([
                        CanonicalKeyValue::Unsigned(first.min(second)),
                        CanonicalKeyValue::Unsigned(first.max(second)),
                    ]),
                ),
                field(
                    1,
                    noncovalent_bond_kind_form_key(&bond.kind().clone().normalize()?),
                ),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                bond.id(),
                incidence_graph.node_of(Entity::NoncovalentBond(bond.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    noncovalent.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let noncovalent = noncovalent
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    for (position, rows) in [
        (EntityBlockPosition::DATIVE_BOND, &dative),
        (EntityBlockPosition::AROMATIC_SYSTEM, &aromatic),
        (EntityBlockPosition::MULTICENTER_BOND, &multicenter),
        (EntityBlockPosition::NONCOVALENT_BOND, &noncovalent),
    ] {
        if !rows.is_empty() {
            candidate.key.entity_blocks.push(PositionedKey {
                position,
                value: sequence(rows.iter().map(|(key, _)| key.clone())),
            });
            candidate
                .entity_order
                .extend(rows.iter().map(|(_, node)| *node));
        }
    }

    Ok(candidate)
}

fn structure_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalCandidate<CanonicalComparisonKey>, Contradiction> {
    let mut candidate = constitution_candidate(molecule, incidence_graph, order)?;
    let mut atom_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Atom)];
    let mut bond_images = vec![0_usize; incidence_graph.entity_count(EntityKind::Bond)];
    for (image, &node) in candidate.entity_order.iter().enumerate() {
        match incidence_graph.entity(node) {
            Entity::Atom(id) => atom_images[id.index()] = image,
            Entity::Bond(id) => {
                bond_images[id.index()] = image - atom_images.len();
            }
            _ => {}
        }
    }

    let remap_ligands = |ligands: &[StereoLigand]| {
        ligands
            .iter()
            .map(|ligand| {
                StereoLigand::new(
                    AtomId(atom_images[ligand.atom_id.index()] as u32),
                    ligand.kind,
                )
            })
            .collect::<Vec<_>>()
    };
    let canonical_frame =
        |ligands: &[StereoLigand],
         configuration: &StereoConfigurationForm|
         -> Result<(Vec<StereoLigand>, StereoConfigurationForm), Contradiction> {
            let mut ligands = remap_ligands(ligands);
            match configuration {
                StereoConfigurationForm::Undetermined => {
                    ligands.sort_unstable();
                    Ok((ligands, StereoConfigurationForm::Undetermined))
                }
                StereoConfigurationForm::Kinded(..) => {
                    let frame = canonical_kinded_stereo_frame(&ligands, configuration)?
                        .expect("kinded stereo frames have an allowed permutation");
                    Ok((frame.ligands, frame.configuration))
                }
            }
        };

    let mut stereo_atoms = molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| {
            let (ligands, configuration) =
                canonical_frame(stereo.ligand_ids(), &stereo.attributes().configuration)?;
            let fields = vec![
                field(
                    0,
                    CanonicalKeyValue::Unsigned(atom_images[stereo.site_id().index()] as u64),
                ),
                field(
                    1,
                    sequence(
                        ligands
                            .iter()
                            .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                    ),
                ),
                field(2, stereo_configuration_form_key(&configuration)),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                stereo.id(),
                incidence_graph.node_of(Entity::StereoAtom(stereo.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    stereo_atoms.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let stereo_atoms = stereo_atoms
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    let mut stereo_bonds = molecule
        .stereo_bonds()
        .iter()
        .map(|stereo| {
            let (ligands, configuration) =
                canonical_frame(stereo.ligand_ids(), &stereo.attributes().configuration)?;
            let fields = vec![
                field(
                    0,
                    CanonicalKeyValue::Unsigned(bond_images[stereo.site_id().index()] as u64),
                ),
                field(
                    1,
                    sequence(
                        ligands
                            .iter()
                            .map(|ligand| stereo_ligand_key(ligand.atom_id.0, ligand.kind)),
                    ),
                ),
                field(2, stereo_configuration_form_key(&configuration)),
            ];
            Ok((
                CanonicalKeyValue::Product(fields),
                stereo.id(),
                incidence_graph.node_of(Entity::StereoBond(stereo.id())),
            ))
        })
        .collect::<Result<Vec<_>, Contradiction>>()?;
    stereo_bonds.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0).then_with(|| lhs.1.cmp(&rhs.1)));
    let stereo_bonds = stereo_bonds
        .into_iter()
        .map(|(key, _, node)| (key, node))
        .collect::<Vec<_>>();

    for (position, rows) in [
        (EntityBlockPosition::STEREO_ATOM, &stereo_atoms),
        (EntityBlockPosition::STEREO_BOND, &stereo_bonds),
    ] {
        if !rows.is_empty() {
            candidate.key.entity_blocks.push(PositionedKey {
                position,
                value: sequence(rows.iter().map(|(key, _)| key.clone())),
            });
            candidate
                .entity_order
                .extend(rows.iter().map(|(_, node)| *node));
        }
    }

    Ok(candidate)
}

fn complete_candidate(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<(CanonicalCandidate<CanonicalComparisonKey>, Molecule), Contradiction> {
    let mut candidate = structure_candidate(molecule, incidence_graph, order)?;
    let remapping = remapping_from_order(molecule, incidence_graph, &candidate.entity_order);
    let complete = molecule.clone().remap(&remapping).reframe()?;
    candidate.key.constraints = constraint_blocks(&complete);
    Ok((candidate, complete))
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

fn canonicalize_topology(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, _) = initial_color_keys(&molecule, &incidence_graph)?;
    let entity_colors = rank_initial_colors(&entity_keys, &[]).entities;
    let TopologySearchInput {
        adapter,
        initial_partition,
    } = topology_search_input(&molecule, &incidence_graph, &entity_colors);
    let leaf_candidate = |order: &[NodeId]| {
        topology_candidate(&molecule, &incidence_graph, order)
            .expect("initial colors established topology normalization")
    };
    let selected = canonical_search_from_partition(
        &adapter,
        initial_partition,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
    );
    let remapping = remapping_from_order(
        &molecule,
        &incidence_graph,
        &selected.candidate.entity_order,
    );

    let canonical = molecule.remap(&remapping).reframe()?;
    Ok((canonical, remapping))
}

fn canonicalize_constitution(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let (entity_keys, incidence_keys) = initial_color_keys(&molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let leaf_candidate = |order: &[NodeId]| {
        constitution_candidate(&molecule, &incidence_graph, order)
            .expect("initial colors established constitution normalization")
    };
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
    );
    let remapping = remapping_from_order(
        &molecule,
        &incidence_graph,
        &selected.candidate.entity_order,
    );

    let canonical = molecule.remap(&remapping).reframe()?;
    Ok((canonical, remapping))
}

fn canonicalize_structure(
    molecule: Molecule,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<(Molecule, MoleculeRemapping), MoleculeCanonicalizeError> {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_color_keys(&molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        &molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        context.para_stereo,
    )?;
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        structure_candidate(&molecule, &incidence_graph, order)
            .expect("structure descriptors established entity normalization")
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let selected = canonical_search_with_generator_filter(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        options,
        &leaf_candidate,
        &filter_generators,
    );
    let remapping = remapping_from_order(
        &molecule,
        &incidence_graph,
        &selected.candidate.entity_order,
    );

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
    let (entity_keys, incidence_keys) = initial_color_keys(&molecule, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        &molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        context.para_stereo,
    )?;
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        complete_candidate(&molecule, &incidence_graph, order)
            .expect("structure descriptors established complete normalization")
            .0
    };
    let selected = canonical_search(
        &adapter,
        &descriptors,
        context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            ..options
        },
        &leaf_candidate,
    );
    let remapping = remapping_from_order(
        &molecule,
        &incidence_graph,
        &selected.candidate.entity_order,
    );
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
            automorphism_pruning: level != DescriptionLevel::Full,
            branch_order: if level == DescriptionLevel::Full {
                BranchOrdering::NodeId
            } else {
                BranchOrdering::BackendCanonical
            },
        };
        match level {
            DescriptionLevel::Topology => canonicalize_topology(self, context, options),
            DescriptionLevel::Constitution => canonicalize_constitution(self, context, options),
            DescriptionLevel::Structure => canonicalize_structure(self, context, options),
            DescriptionLevel::Full => canonicalize_full(self, context, options),
        }
    }

    /// The comparison key of the canonical form at `level`: the key canonical equality compares.
    fn canonical_key_at(
        &self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<CanonicalComparisonKey, MoleculeCanonicalizeError> {
        if level == DescriptionLevel::Full {
            let molecule = &self.clone().normalize()?;
            let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
            let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
            let colors = rank_initial_colors(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
            let (partition, _) = structure_partition(
                molecule,
                &incidence_graph,
                &adapter,
                &entity_keys,
                context.para_stereo,
            )?;
            let descriptors = partition.cell_indices(adapter.graph().node_count());
            let leaf_candidate = |order: &[NodeId]| {
                complete_candidate(molecule, &incidence_graph, order)
                    .expect("normalized structure descriptors established complete normalization")
                    .0
            };
            return Ok(canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                CanonicalSearchOptions {
                    automorphism_pruning: false,
                    branch_order: BranchOrdering::NodeId,
                },
                &leaf_candidate,
            )
            .candidate
            .key);
        }
        let molecule = self;
        let incidence_level = match level {
            DescriptionLevel::Topology => IncidenceLevel::Topology,
            DescriptionLevel::Constitution => IncidenceLevel::Constitution,
            DescriptionLevel::Structure => IncidenceLevel::Full,
            DescriptionLevel::Full => unreachable!("handled before incidence construction"),
        };
        let incidence_graph = molecule.incidence_graph(incidence_level);
        let (entity_keys, incidence_keys) = initial_color_keys(molecule, &incidence_graph)?;
        let automorphism_pruning = level != DescriptionLevel::Structure;
        let options = CanonicalSearchOptions {
            automorphism_pruning,
            branch_order: if automorphism_pruning {
                BranchOrdering::BackendCanonical
            } else {
                BranchOrdering::NodeId
            },
        };

        if level == DescriptionLevel::Topology {
            let entity_colors = rank_initial_colors(&entity_keys, &[]).entities;
            let TopologySearchInput {
                adapter,
                initial_partition,
            } = topology_search_input(molecule, &incidence_graph, &entity_colors);
            let leaf_candidate = |order: &[NodeId]| {
                topology_candidate(molecule, &incidence_graph, order)
                    .expect("initial colors established topology normalization")
            };
            return Ok(canonical_search_from_partition(
                &adapter,
                initial_partition,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
            )
            .candidate
            .key);
        }

        let colors = rank_initial_colors(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
        let key = match level {
            DescriptionLevel::Constitution => {
                let descriptors =
                    constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
                let leaf_candidate = |order: &[NodeId]| {
                    constitution_candidate(molecule, &incidence_graph, order)
                        .expect("initial colors established constitution normalization")
                };
                canonical_search(
                    &adapter,
                    &descriptors,
                    context.automorphism_algorithm,
                    options,
                    &leaf_candidate,
                )
                .candidate
                .key
            }
            DescriptionLevel::Structure => {
                let (partition, _) = structure_partition(
                    molecule,
                    &incidence_graph,
                    &adapter,
                    &entity_keys,
                    context.para_stereo,
                )?;
                let descriptors = partition.cell_indices(adapter.graph().node_count());
                let leaf_candidate = |order: &[NodeId]| {
                    structure_candidate(molecule, &incidence_graph, order)
                        .expect("structure descriptors established entity normalization")
                };
                canonical_search(
                    &adapter,
                    &descriptors,
                    context.automorphism_algorithm,
                    options,
                    &leaf_candidate,
                )
                .candidate
                .key
            }
            DescriptionLevel::Topology | DescriptionLevel::Full => {
                unreachable!("handled before selected-layer search")
            }
        };
        Ok(key)
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
        let level = molecule_canonicalize_level(self);
        self.canonical_key_at(level, context)
            .map(|key| hash_value(&key))
    }

    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool {
        if self == other {
            return true;
        }
        let level = molecule_canonicalize_level(self).max(molecule_canonicalize_level(other));
        match (
            self.canonical_key_at(level, context),
            other.canonical_key_at(level, context),
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
