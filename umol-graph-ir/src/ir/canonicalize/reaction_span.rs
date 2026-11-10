//! Reaction-span and reaction canonicalization.

use umol_graph_core::NodeId;
use umol_perm::Permutation;

use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::delta::EntitySpan;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::error::Contradiction;
use super::super::id::{AtomId, BondId};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::num::NumForm;
use super::super::reaction::Reaction;
use super::super::reaction_span::ReactionSpan;
use super::super::remap::MoleculeRemapping;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{StereoAtomForm, StereoBondForm, StereoConfigurationForm};
use super::super::traits::{FrameTransport, Normalize, Reframe};
use super::{
    aromatic_system_constraint_form_key, atom_constraint_form_key, atom_inherent_fields,
    bond_constraint_form_key, bond_inherent_fields, canonical_search,
    canonical_search_from_partition, constitution_partition_descriptors, constraint_key_span,
    constraint_span_key, dative_bond_constraint_form_key, electron_counts_form_key,
    entity_span_key, hash_value, incidence_key, index_key, is_empty_sequence, molecule_remapping,
    multicenter_bond_constraint_form_key, noncovalent_bond_constraint_form_key,
    noncovalent_bond_kind_form_key, normalized_entity_span_key, num_form_key, option,
    positioned_product, product, push_span_block, rank_initial_colors, reaction_canonicalize_level,
    reaction_span_canonicalize_level, reaction_span_counts, relation_span_key, sequence,
    stereo_atom_constraint_form_key, stereo_bond_constraint_form_key,
    stereo_configuration_form_key, stereo_frame_permutations, stereo_kind_key, stereo_ligand_key,
    topology_initial_partition, unpaired_electrons_form_key, variant, AutomorphismAdapter,
    BranchOrdering, CanonicalCandidate, CanonicalComparisonKey, CanonicalKeyValue,
    CanonicalSearchOptions, CanonicalSearchResult, Canonicalize, CanonicalizeContext,
    ConstraintBlockPosition, DescriptionLevel, EntityBlockPosition, InitialColorKey, PositionedKey,
    ReactionCanonicalizeError, ReactionSpanCanonicalizeError,
};

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
                    option(attributes.configuration.kind().map(stereo_kind_key)),
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
                    option(attributes.configuration.kind().map(stereo_kind_key)),
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

/// The stereo descriptor of a span entity, minimized over the frame permutations of its kind.
///
/// One permutation acts on the ligand classes and on the configuration of each present side.
#[cfg_attr(not(test), expect(dead_code))]
fn stereo_refinement_descriptor<T>(
    site_class: u32,
    ligand_classes: &[(u32, StereoLigandKind)],
    span: &EntitySpan<T>,
    configuration: impl Fn(&T) -> &StereoConfigurationForm,
) -> Result<CanonicalKeyValue, Contradiction> {
    let descriptor = |ligands: &[(u32, StereoLigandKind)],
                      configuration: StereoConfigurationForm| {
        product([
            CanonicalKeyValue::Unsigned(site_class.into()),
            sequence(
                ligands
                    .iter()
                    .map(|&(class, kind)| stereo_ligand_key(class, kind)),
            ),
            stereo_configuration_form_key(&configuration),
        ])
    };

    let kind = [span.lhs(), span.rhs()]
        .into_iter()
        .flatten()
        .find_map(|attributes| configuration(attributes).kind());
    match kind {
        None => {
            let mut ligands = ligand_classes.to_vec();
            ligands.sort_unstable();
            normalized_entity_span_key(span, |_| {
                Ok(descriptor(&ligands, StereoConfigurationForm::Undetermined))
            })
        }
        Some(kind) => {
            if ligand_classes.len() != kind.degree() {
                return Err(Contradiction);
            }
            stereo_frame_permutations(kind)
                .map(|permutation| {
                    let ligands = permutation.act(ligand_classes);
                    normalized_entity_span_key(span, |attributes| {
                        configuration(attributes)
                            .apply(permutation)
                            .ok_or(Contradiction)?
                            .normalize()
                            .map(|configuration| descriptor(&ligands, configuration))
                    })
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .min()
                .ok_or(Contradiction)
        }
    }
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

fn canonical_candidate(
    span: &ReactionSpan,
    level: DescriptionLevel,
    context: &CanonicalizeContext,
    options: CanonicalSearchOptions,
) -> Result<CanonicalSearchResult<CanonicalComparisonKey>, Contradiction> {
    let incidence_level = match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    };
    let incidence_graph = span.incidence_graph(incidence_level);
    let (entity_keys, incidence_keys) = initial_color_keys(span, &incidence_graph)?;
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let leaf_candidate = |order: &[NodeId]| {
        let remapping = remapping_from_order(span, &incidence_graph, order);
        let remapped = span.clone().remap(&remapping);
        let action = remapped.representative_action();
        let reframed = remapped
            .reframe_by(&action)
            .expect("integrity established compatible reaction-span frame actions");
        CanonicalCandidate {
            key: comparison_key(&reframed, level)
                .expect("initial colors established span normalization"),
            entity_order: order.to_vec(),
        }
    };
    let selected = match level {
        DescriptionLevel::Topology => canonical_search_from_partition(
            &adapter,
            topology_initial_partition(&adapter, &colors),
            context.automorphism_algorithm,
            options,
            &leaf_candidate,
        ),
        DescriptionLevel::Constitution | DescriptionLevel::Structure | DescriptionLevel::Full => {
            let descriptors =
                constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
            canonical_search(
                &adapter,
                &descriptors,
                context.automorphism_algorithm,
                options,
                &leaf_candidate,
            )
        }
    };
    Ok(selected)
}

impl ReactionSpan {
    /// The canonical form at `level` and its source-to-canonical remapping.
    fn tracked_canonicalize_at(
        self,
        level: DescriptionLevel,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), ReactionSpanCanonicalizeError> {
        let automorphism_pruning = matches!(
            level,
            DescriptionLevel::Topology | DescriptionLevel::Constitution
        );
        let options = CanonicalSearchOptions {
            automorphism_pruning,
            branch_order: if automorphism_pruning {
                BranchOrdering::BackendCanonical
            } else {
                BranchOrdering::NodeId
            },
        };
        let normalized = self.normalize()?;
        let candidate = canonical_candidate(&normalized, level, context, options)?.candidate;
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
        let automorphism_pruning = matches!(
            level,
            DescriptionLevel::Topology | DescriptionLevel::Constitution
        );
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
            return Ok(canonical_candidate(&normalized, level, context, options)?
                .candidate
                .key);
        }
        Ok(canonical_candidate(self, level, context, options)?
            .candidate
            .key)
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
