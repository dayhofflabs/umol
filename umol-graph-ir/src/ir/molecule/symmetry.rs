//! Molecule graph-symmetry data, queries, colors and remapping checks.

use umol_graph_core::{AutomorphismAlgorithm, EdgeId, NodeId};
use umol_perm::{Orientation, Permutation};

use super::super::electrons::ElectronCountsForm;
use super::super::entity::{DescriptionLevel, Entity};
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::StereoLigand;
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{
    coset_apply_permutation, StereoConfigurationForm, StereoCoset, StereoKind,
};
use super::super::symmetry::{
    rank_colors, remapping_exists, remapping_preserves_incidences, AutomorphismAdapter,
    InitialColors, OrderedPartition,
};
use super::super::traits::{AsColor, AsLit, FrameTransport, Normalize};
use super::super::view::GraphSymmetryView;
use super::Molecule;

/// Auxiliary search data computed for a molecule's graph-symmetry readings.
///
/// The value holds incidence, colors, a partition and config, without a molecule.
/// Queries read a separately supplied molecule. The data must describe that molecule for
/// meaningful results; after changing the molecule, compute new data. Pairing mismatched
/// data and a molecule does not panic but has no specified semantic result.
///
/// # Semantic properties
///
/// Cloning preserves the computed data and its contextual requirements. Computation and
/// queries preserve the supplied molecule's exact representation without normalization.
#[derive(Clone, Debug)]
pub struct GraphSymmetry {
    incidence: IncidenceGraph,
    adapter: AutomorphismAdapter,
    partition: OrderedPartition,
    cell_indices: Vec<u32>,
    config: GraphSymmetryConfig,
}

impl GraphSymmetry {
    pub(crate) fn cell(&self, entity: Entity) -> Option<u32> {
        (entity.kind_id() < self.incidence.entity_count(entity.kind()))
            .then(|| self.cell_indices[self.incidence.node_of(entity).index()])
    }
}

impl Molecule {
    /// Compute graph-symmetry data and return a view borrowing this molecule.
    ///
    /// Forms are read as supplied. Encoding-incompatible attributes receive distinct entity
    /// colors, giving singleton atom orbits. This fallback does not diagnose the attributes.
    pub fn graph_symmetry(&self, config: &GraphSymmetryConfig) -> GraphSymmetryView<'_> {
        let incidence = self.incidence_graph(IncidenceLevel::Constitution);
        let stereo_compatible = self.stereo_atoms().iter().all(|stereo| {
            let attributes = stereo.attributes();
            if attributes.configuration.kind().is_some_and(|kind| {
                check_stereo_atom_kind(Entity::StereoAtom(stereo.id()), kind).is_err()
            }) {
                return false;
            }
            let action = Permutation::identity(stereo.ligand_count());
            if config.include_constraints {
                attributes.clone().reframe_by(&action).is_some()
            } else {
                attributes
                    .configuration
                    .clone()
                    .reframe_by(&action)
                    .is_some()
            }
        }) && self.stereo_bonds().iter().all(|stereo| {
            let attributes = stereo.attributes();
            if attributes.configuration.kind().is_some_and(|kind| {
                check_stereo_bond_kind(Entity::StereoBond(stereo.id()), kind).is_err()
            }) {
                return false;
            }
            let action = Permutation::identity(stereo.ligand_count());
            if config.include_constraints {
                attributes.clone().reframe_by(&action).is_some()
            } else {
                StereoKind::CisTrans.class_key().space().allows(action)
                    && attributes
                        .configuration
                        .clone()
                        .reframe_by(&action)
                        .is_some()
            }
        });
        let colors = initial_colors(
            self,
            &incidence,
            DescriptionLevel::Constitution,
            config.include_constraints,
        )
        .and_then(|colors| {
            if stereo_compatible {
                Ok(colors)
            } else {
                Err(Contradiction)
            }
        })
        .unwrap_or_else(|_| InitialColors {
            entities: (0..incidence.graph().node_count() as u32).collect(),
            incidences: (incidence.graph().node_count() as u32
                ..(incidence.graph().node_count() + incidence.graph().edge_count()) as u32)
                .collect(),
        });
        let adapter = AutomorphismAdapter::new(&incidence, &colors);
        let partition = OrderedPartition::from_descriptors(&adapter.colors).refine(adapter.graph());
        let cell_indices = partition.cell_indices(adapter.graph().node_count());
        GraphSymmetryView::new(
            self,
            GraphSymmetry {
                incidence,
                adapter,
                partition,
                cell_indices,
                config: *config,
            },
        )
    }

    /// Whether no remapping carries every supplied configuration to its mirror image.
    ///
    /// The auxiliary data must describe this molecule. Mismatched context has no specified
    /// semantic result and does not panic. A molecule without stereo records is achiral.
    pub fn is_chiral(&self, symmetry: &GraphSymmetry) -> bool {
        !remapping_exists_at(
            self,
            symmetry,
            DescriptionLevel::Structure,
            Orientation::Improper,
            &[],
            &[],
        )
    }

    /// Whether a configuration-preserving remapping sends atom `a` to atom `b`.
    ///
    /// The data must describe this molecule. Mismatched context has no specified result
    /// and does not panic. An unavailable atom id returns false.
    pub fn same_proper_orbit(&self, symmetry: &GraphSymmetry, a: AtomId, b: AtomId) -> bool {
        if self.atoms().get(a).is_none() || self.atoms().get(b).is_none() {
            return false;
        }
        if a == b {
            return true;
        }
        remapping_exists_at(
            self,
            symmetry,
            DescriptionLevel::Structure,
            Orientation::Proper,
            &[a],
            &[b],
        )
    }

    /// Whether a preserving or mirror remapping sends atom `a` to atom `b`.
    ///
    /// The data must describe this molecule. Mismatched context has no specified result
    /// and does not panic. An unavailable atom id returns false.
    pub fn same_star_orbit(&self, symmetry: &GraphSymmetry, a: AtomId, b: AtomId) -> bool {
        self.same_proper_orbit(symmetry, a, b)
            || (self.atoms().get(a).is_some()
                && self.atoms().get(b).is_some()
                && remapping_exists_at(
                    self,
                    symmetry,
                    DescriptionLevel::Structure,
                    Orientation::Improper,
                    &[a],
                    &[b],
                ))
    }

    /// All atoms in `a`'s proper orbit, in atom-id order.
    ///
    /// The data must describe this molecule. Mismatched context has no specified result
    /// and does not panic. An unavailable atom id gives an empty orbit.
    pub fn proper_orbit_of(&self, symmetry: &GraphSymmetry, a: AtomId) -> Vec<AtomId> {
        self.atoms()
            .ids()
            .filter(|&b| self.same_proper_orbit(symmetry, a, b))
            .collect()
    }

    /// All atoms in `a`'s proper-or-mirror orbit, in atom-id order.
    ///
    /// The data must describe this molecule. Mismatched context has no specified result
    /// and does not panic. An unavailable atom id gives an empty orbit.
    pub fn star_orbit_of(&self, symmetry: &GraphSymmetry, a: AtomId) -> Vec<AtomId> {
        self.atoms()
            .ids()
            .filter(|&b| self.same_star_orbit(symmetry, a, b))
            .collect()
    }
}

pub(crate) fn remapping_exists_at(
    molecule: &Molecule,
    symmetry: &GraphSymmetry,
    level: DescriptionLevel,
    orientation: Orientation,
    source: &[AtomId],
    target: &[AtomId],
) -> bool {
    molecule_remapping_exists(molecule, symmetry, level, orientation, source, target, None)
}

pub(crate) fn remapping_exists_for(
    molecule: &Molecule,
    symmetry: &GraphSymmetry,
    entity: Entity,
    coset: u32,
    fixed: &[AtomId],
) -> bool {
    molecule_remapping_exists(
        molecule,
        symmetry,
        DescriptionLevel::Structure,
        Orientation::Proper,
        fixed,
        fixed,
        Some((entity, coset)),
    )
}

fn molecule_remapping_exists(
    molecule: &Molecule,
    symmetry: &GraphSymmetry,
    level: DescriptionLevel,
    orientation: Orientation,
    source: &[AtomId],
    target: &[AtomId],
    changed: Option<(Entity, u32)>,
) -> bool {
    if source.len() != target.len() {
        return false;
    }
    for (&a, &b) in source.iter().zip(target) {
        if molecule.atoms().get(a).is_none()
            || molecule.atoms().get(b).is_none()
            || symmetry.cell(Entity::Atom(a)).is_none()
            || symmetry.cell(Entity::Atom(a)) != symmetry.cell(Entity::Atom(b))
        {
            return false;
        }
    }
    let partition = |atoms: &[AtomId]| {
        atoms
            .iter()
            .fold(symmetry.partition.clone(), |partition, &atom| {
                let node = symmetry.incidence.node_of(Entity::Atom(atom));
                let cell = partition
                    .cells
                    .iter()
                    .position(|cell| cell.contains(&node))
                    .unwrap();
                partition.individualize(cell, node)
            })
    };
    let matches_stereo = |images: &[NodeId],
                          orientation: Orientation,
                          changed: Option<(Entity, u32)>| {
        level == DescriptionLevel::Constitution
            || if orientation == Orientation::Improper && changed.is_none() {
                remapping_mirrors_stereo(
                    molecule,
                    &symmetry.incidence,
                    images,
                    symmetry.config.include_constraints,
                )
            } else if let (Orientation::Proper, Some((entity, coset))) = (orientation, changed) {
                remapping_changes_coset(
                    molecule,
                    &symmetry.incidence,
                    images,
                    entity,
                    coset,
                    symmetry.config.include_constraints,
                )
            } else {
                remapping_matches_stereo(
                    molecule,
                    &symmetry.incidence,
                    images,
                    symmetry.config.include_constraints,
                    orientation,
                    changed,
                )
            }
    };
    remapping_exists(
        &symmetry.adapter,
        partition(source),
        &symmetry.adapter,
        partition(target),
        symmetry.config.automorphism_algorithm,
        &|images| {
            remapping_preserves_incidences(&symmetry.adapter, &symmetry.adapter, images)
                && matches_stereo(images, orientation, changed)
        },
        true,
        &|generators: &mut Vec<Vec<NodeId>>| {
            if let Some((entity, coset)) = changed {
                let record = match entity {
                    Entity::StereoAtom(id) => molecule.stereo_atoms().get(id).map(|stereo| {
                        (
                            stereo.attributes().configuration.kind(),
                            stereo.ligand_ids(),
                            Entity::Atom(stereo.site_id()),
                        )
                    }),
                    Entity::StereoBond(id) => molecule.stereo_bonds().get(id).map(|stereo| {
                        (
                            stereo.attributes().configuration.kind(),
                            stereo.ligand_ids(),
                            Entity::Bond(stereo.site_id()),
                        )
                    }),
                    _ => None,
                };
                let Some((Some(kind), ligands, site)) = record else {
                    generators.clear();
                    return;
                };
                generators.retain(|images| {
                    if entity_image(&symmetry.incidence, images, site) != Some(site)
                        || !matches_stereo(images, Orientation::Proper, None)
                    {
                        return false;
                    }
                    let mapped = ligands
                        .iter()
                        .map(|ligand| {
                            let Entity::Atom(atom) = entity_image(
                                &symmetry.incidence,
                                images,
                                Entity::Atom(ligand.atom_id),
                            )?
                            else {
                                return None;
                            };
                            Some(StereoLigand::new(atom, ligand.kind))
                        })
                        .collect::<Option<Vec<_>>>();
                    let Some(mapped) = mapped else {
                        return false;
                    };
                    let Some(frame) = Permutation::between(&mapped, ligands) else {
                        return false;
                    };
                    coset_apply_permutation(&StereoCoset::Lit(coset), frame, kind)
                        == Some(StereoCoset::Lit(coset))
                });
            } else {
                generators.retain(|images| matches_stereo(images, Orientation::Proper, None));
            }
        },
    )
}

fn entity_image(incidence: &IncidenceGraph, images: &[NodeId], entity: Entity) -> Option<Entity> {
    if entity.kind_id() >= incidence.entity_count(entity.kind()) {
        return None;
    }
    let image = *images.get(incidence.node_of(entity).index())?;
    (image.index() < incidence.graph().node_count()).then(|| incidence.entity(image))
}

fn remapping_matches_stereo(
    molecule: &Molecule,
    incidence: &IncidenceGraph,
    images: &[NodeId],
    include_constraints: bool,
    orientation: Orientation,
    changed: Option<(Entity, u32)>,
) -> bool {
    let image = |entity: Entity| {
        if incidence.entity_count(entity.kind()) != 0 {
            return entity_image(incidence, images, entity);
        }
        match entity {
            Entity::StereoAtom(id) => {
                let source = molecule.stereo_atoms().get(id)?;
                let Entity::Atom(site) =
                    entity_image(incidence, images, Entity::Atom(source.site_id()))?
                else {
                    return None;
                };
                molecule
                    .stereo_atoms()
                    .iter()
                    .find(|stereo| stereo.site_id() == site)
                    .map(|stereo| Entity::StereoAtom(stereo.id()))
            }
            Entity::StereoBond(id) => {
                let source = molecule.stereo_bonds().get(id)?;
                let Entity::Bond(site) =
                    entity_image(incidence, images, Entity::Bond(source.site_id()))?
                else {
                    return None;
                };
                molecule
                    .stereo_bonds()
                    .iter()
                    .find(|stereo| stereo.site_id() == site)
                    .map(|stereo| Entity::StereoBond(stereo.id()))
            }
            _ => entity_image(incidence, images, entity),
        }
    };
    let record = |entity| match entity {
        Entity::StereoAtom(id) => {
            let stereo = molecule.stereo_atoms().get(id)?;
            Some((
                &stereo.attributes().configuration,
                stereo.ligand_ids(),
                Entity::Atom(stereo.site_id()),
            ))
        }
        Entity::StereoBond(id) => {
            let stereo = molecule.stereo_bonds().get(id)?;
            Some((
                &stereo.attributes().configuration,
                stereo.ligand_ids(),
                Entity::Bond(stereo.site_id()),
            ))
        }
        _ => None,
    };
    for entity in molecule
        .stereo_atoms()
        .ids()
        .map(Entity::StereoAtom)
        .chain(molecule.stereo_bonds().ids().map(Entity::StereoBond))
    {
        let Some(target_entity) = image(entity) else {
            return false;
        };
        if entity.kind() != target_entity.kind() {
            return false;
        }
        let Some((configuration, ligands, site)) = record(entity) else {
            return false;
        };
        let Some((target_configuration, target_ligands, target_site)) = record(target_entity)
        else {
            return false;
        };
        if image(site) != Some(target_site) {
            return false;
        }
        let mapped = ligands
            .iter()
            .map(|ligand| {
                let Entity::Atom(atom) = image(Entity::Atom(ligand.atom_id))? else {
                    return None;
                };
                Some(StereoLigand::new(atom, ligand.kind))
            })
            .collect::<Option<Vec<_>>>();
        let Some(mapped) = mapped else {
            return false;
        };
        let Some(frame) = Permutation::between(&mapped, target_ligands) else {
            return false;
        };
        let mut configuration = configuration.clone();
        let mut target_configuration = target_configuration.clone();
        if let Some((focus, coset)) = changed {
            if entity == focus {
                let StereoConfigurationForm::Kinded(kind, own) = &configuration else {
                    return false;
                };
                configuration = StereoConfigurationForm::kinded(*kind, own.as_lit().unwrap_or(0));
            }
            if target_entity == focus {
                let Some(kind) = target_configuration.kind() else {
                    return false;
                };
                target_configuration = StereoConfigurationForm::kinded(kind, coset);
            }
        }
        let transported = match configuration {
            StereoConfigurationForm::Undetermined => Some(StereoConfigurationForm::Undetermined),
            StereoConfigurationForm::Kinded(kind, coset) => {
                if kind.degree() != frame.degree() {
                    return false;
                }
                let action = if orientation == Orientation::Improper {
                    kind.mirror_permutation().compose(frame)
                } else {
                    frame
                };
                coset_apply_permutation(&coset, action, kind)
                    .map(|coset| StereoConfigurationForm::kinded(kind, coset))
            }
        };
        if transported.as_ref() != Some(&target_configuration) {
            return false;
        }
        if include_constraints {
            let preserved = match (entity, target_entity) {
                (Entity::StereoAtom(source), Entity::StereoAtom(target)) => molecule
                    .stereo_atom(source)
                    .attributes()
                    .constraints
                    .clone()
                    .reframe_by(&frame)
                    .is_some_and(|mapped| {
                        mapped == molecule.stereo_atom(target).attributes().constraints
                    }),
                (Entity::StereoBond(source), Entity::StereoBond(target)) => molecule
                    .stereo_bond(source)
                    .attributes()
                    .constraints
                    .clone()
                    .reframe_by(&frame)
                    .is_some_and(|mapped| {
                        mapped == molecule.stereo_bond(target).attributes().constraints
                    }),
                _ => false,
            };
            if !preserved {
                return false;
            }
        }
    }
    true
}

fn remapping_mirrors_stereo(
    molecule: &Molecule,
    incidence: &IncidenceGraph,
    images: &[NodeId],
    include_constraints: bool,
) -> bool {
    remapping_matches_stereo(
        molecule,
        incidence,
        images,
        include_constraints,
        Orientation::Improper,
        None,
    )
}

fn remapping_changes_coset(
    molecule: &Molecule,
    incidence: &IncidenceGraph,
    images: &[NodeId],
    entity: Entity,
    coset: u32,
    include_constraints: bool,
) -> bool {
    remapping_matches_stereo(
        molecule,
        incidence,
        images,
        include_constraints,
        Orientation::Proper,
        Some((entity, coset)),
    )
}

/// Whether the remapping `images` carries every stereo entity of `source` onto a stereo entity of
/// `target` with the same site and the same form.
///
/// `images` holds the image of each node of `incidence_graph`, the incidence graph of `source`,
/// and `target` has the entity counts of `source`. A generator of `source` is such a remapping of
/// `source` onto itself.
pub(crate) fn remapping_preserves_stereo(
    source: &Molecule,
    target: &Molecule,
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

    for stereo in source.stereo_atoms().iter() {
        let Entity::StereoAtom(image_id) = entity_image(Entity::StereoAtom(stereo.id())) else {
            return false;
        };
        let image = target.stereo_atom(image_id);
        if entity_image(Entity::Atom(stereo.site_id())) != Entity::Atom(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        if !form_preserved(stereo.attributes(), image.attributes(), action) {
            return false;
        }
    }

    for stereo in source.stereo_bonds().iter() {
        let Entity::StereoBond(image_id) = entity_image(Entity::StereoBond(stereo.id())) else {
            return false;
        };
        let image = target.stereo_bond(image_id);
        if entity_image(Entity::Bond(stereo.site_id())) != Entity::Bond(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        if !form_preserved(stereo.attributes(), image.attributes(), action) {
            return false;
        }
    }

    true
}

pub(crate) fn retain_stereo_preserving_generators(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    generators: &mut Vec<Vec<NodeId>>,
) {
    if molecule.stereo_atoms().count() == 0 && molecule.stereo_bonds().count() == 0 {
        return;
    }
    generators.retain(|generator| {
        remapping_preserves_stereo(molecule, molecule, incidence_graph, generator)
    });
}

/// The color of the incidence at `edge`: its role and, for an atom of an aromatic system or a
/// multicenter bond, the electron count of that atom.
pub(crate) fn incidence_color<'a>(
    molecule: &Molecule,
    incidence_graph: &'a IncidenceGraph,
    edge: EdgeId,
) -> (&'a Incidence, Option<i64>) {
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
            let atom_position = atom_position.expect("incidence atom belongs to its overlay");
            electrons.at(AtomPosition::from(atom_position))
        }
        _ => None,
    };
    (incidence, electron_count)
}

/// The integers of the colors of the supplied molecule at `level`: one per entity, by incidence
/// node, and one per incidence, by incidence edge.
///
/// The colors of each entity kind are ranked in kind order, and the incidences continue after the
/// entities. The color of an incidence is its role and, for an atom of an aromatic system or a
/// multicenter bond, the electron count of that atom. At Topology the incidences are not ranked,
/// since the topology search input does not read them.
/// `include_constraints` selects the frame-independent entity constraints in the colors and the
/// frame-relative stereo constraints in compatibility checks. Configuration compatibility is
/// checked in both modes. No form is normalized.
#[inline]
pub(crate) fn initial_colors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    include_constraints: bool,
) -> Result<InitialColors, Contradiction> {
    let mut entities = vec![0; incidence_graph.graph().node_count()];
    let position = |entity: Entity| incidence_graph.node_of(entity).index();

    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| (atom.attributes(), position(Entity::Atom(atom.id()))));
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

    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| (bond.attributes(), position(Entity::Bond(bond.id()))));
    next = if include_constraints {
        rank_colors(
            bonds
                .map(|(form, position)| (form.as_color(), position))
                .collect(),
            next,
            &mut entities,
        )
    } else {
        rank_colors(
            bonds
                .map(|(form, position)| (form.as_field_color(), position))
                .collect(),
            next,
            &mut entities,
        )
    };

    if level == DescriptionLevel::Topology {
        return Ok(InitialColors {
            entities,
            incidences: Vec::new(),
        });
    }

    let dative = molecule
        .dative_bonds()
        .iter()
        .map(|bond| (bond.attributes(), position(Entity::DativeBond(bond.id()))));
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

    let aromatic = molecule.aromatic_systems().iter().map(|system| {
        let attributes = system.attributes();
        if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                if counts.len() != system.atom_ids().len())
        {
            return Err(Contradiction);
        }
        Ok((attributes, position(Entity::AromaticSystem(system.id()))))
    });
    next = if include_constraints {
        let mut colors = Vec::with_capacity(molecule.aromatic_systems().count());
        for entry in aromatic {
            let (form, position) = entry?;
            colors.push((form.as_color(), position));
        }
        rank_colors(colors, next, &mut entities)
    } else {
        let mut colors = Vec::with_capacity(molecule.aromatic_systems().count());
        for entry in aromatic {
            let (form, position) = entry?;
            colors.push((form.as_field_color(), position));
        }
        rank_colors(colors, next, &mut entities)
    };

    let multicenter = molecule.multicenter_bonds().iter().map(|bond| {
        let attributes = bond.attributes();
        if matches!(&attributes.electrons, ElectronCountsForm::Lit(counts)
                if counts.len() != bond.atom_ids().len())
        {
            return Err(Contradiction);
        }
        Ok((attributes, position(Entity::MulticenterBond(bond.id()))))
    });
    next = if include_constraints {
        let mut colors = Vec::with_capacity(molecule.multicenter_bonds().count());
        for entry in multicenter {
            let (form, position) = entry?;
            colors.push((form.as_color(), position));
        }
        rank_colors(colors, next, &mut entities)
    } else {
        let mut colors = Vec::with_capacity(molecule.multicenter_bonds().count());
        for entry in multicenter {
            let (form, position) = entry?;
            colors.push((form.as_field_color(), position));
        }
        rank_colors(colors, next, &mut entities)
    };

    let noncovalent = molecule.noncovalent_bonds().iter().map(|bond| {
        (
            bond.attributes(),
            position(Entity::NoncovalentBond(bond.id())),
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

    if level >= DescriptionLevel::Structure {
        let stereo_atoms = molecule.stereo_atoms().iter().map(|stereo| {
            let attributes = stereo.attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_atom_kind(Entity::StereoAtom(stereo.id()), kind)
                    .map_err(|_| Contradiction)?;
            }
            let action = Permutation::identity(stereo.ligand_count());
            if include_constraints {
                attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            } else {
                attributes
                    .configuration
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            }
            Ok((attributes, position(Entity::StereoAtom(stereo.id()))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(molecule.stereo_atoms().count());
            for entry in stereo_atoms {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(molecule.stereo_atoms().count());
            for entry in stereo_atoms {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };

        let stereo_bonds = molecule.stereo_bonds().iter().map(|stereo| {
            let attributes = stereo.attributes();
            if let Some(kind) = attributes.configuration.kind() {
                check_stereo_bond_kind(Entity::StereoBond(stereo.id()), kind)
                    .map_err(|_| Contradiction)?;
            }
            let action = Permutation::identity(stereo.ligand_count());
            if include_constraints {
                attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            } else {
                if !StereoKind::CisTrans.class_key().space().allows(action) {
                    return Err(Contradiction);
                }
                attributes
                    .configuration
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
            }
            Ok((attributes, position(Entity::StereoBond(stereo.id()))))
        });
        next = if include_constraints {
            let mut colors = Vec::with_capacity(molecule.stereo_bonds().count());
            for entry in stereo_bonds {
                let (form, position) = entry?;
                colors.push((form.as_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        } else {
            let mut colors = Vec::with_capacity(molecule.stereo_bonds().count());
            for entry in stereo_bonds {
                let (form, position) = entry?;
                colors.push((form.as_field_color(), position));
            }
            rank_colors(colors, next, &mut entities)
        };
    }

    let incidence_colors = incidence_graph
        .incidences()
        .map(|(edge, _)| {
            (
                incidence_color(molecule, incidence_graph, edge),
                edge.index(),
            )
        })
        .collect();
    let mut incidences = vec![0; incidence_graph.graph().edge_count()];
    rank_colors(incidence_colors, next, &mut incidences);

    Ok(InitialColors {
        entities,
        incidences,
    })
}

/// Configuration for molecule graph symmetry and derived stereo readings.
///
/// Entity constraints are excluded by default. Molecule-level constraints are always excluded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GraphSymmetryConfig {
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub max_assignments: u32,
    pub include_constraints: bool,
}

impl Default for GraphSymmetryConfig {
    fn default() -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            max_assignments: 4096,
            include_constraints: false,
        }
    }
}

#[cfg(test)]
mod tests;
