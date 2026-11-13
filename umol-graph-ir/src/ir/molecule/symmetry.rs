//! Graph-automorphism symmetry of a molecule, graded into proper vs improper.

use std::hash::{DefaultHasher, Hash, Hasher};

use umol_graph_core::{AutomorphismAlgorithm, AutomorphismOutput, EdgeId, NodeId};
use umol_perm::{Orientation, OrientedPermutation, OrientedPermutationGroup, Permutation};

use super::super::coloring::MoleculeColoring;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::{DescriptionLevel, Entity, EntityKind};
use super::super::error::Contradiction;
use super::super::id::{AtomId, AtomPosition, StereoAtomId, StereoBondId, StereoLigandPosition};
use super::super::incidence::{Incidence, IncidenceGraph, IncidenceLevel};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::stereo::integrity::{check_stereo_atom_kind, check_stereo_bond_kind};
use super::super::stereo::{
    coset_apply_permutation, StereoConfigurationForm, StereoCoset, StereoKind, Stereogenicity,
    Topicity,
};
use super::super::symmetry::{rank_colors, InitialColors};
use super::super::traits::{AsColor, AsLit, FrameTransport, Normalize};
use super::Molecule;
#[cfg(test)]
use super::MoleculeEntries;

/// Configuration for [`Molecule::graph_symmetry`].
pub struct GraphSymmetryConfig<C: MoleculeColoring> {
    pub coloring: C,
    pub iterate_to_fixpoint: bool,
    pub max_iterations: u32,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

/// The converged, orientation-graded graph-automorphism symmetry of a molecule
/// under a coloring. Owned and self-contained: molecule-level orbit queries read
/// the stored partitions; per-carrier stabilizer re-runs use the stored colors.
pub struct GraphSymmetry {
    incidence: IncidenceGraph,
    colors: Vec<u64>, // converged node colors; for site_stabilizer re-runs
    proper_orbits: Vec<NodeId>,
    star_orbits: Vec<NodeId>,
    chiral: bool,
    automorphism_algorithm: AutomorphismAlgorithm,
}

impl Molecule {
    /// Build the molecule's graded graph-automorphism symmetry under `cfg.coloring`.
    pub fn graph_symmetry<C: MoleculeColoring>(
        &self,
        cfg: &GraphSymmetryConfig<C>,
    ) -> GraphSymmetry {
        let incidence = self.incidence_graph(IncidenceLevel::Full);
        let node_count = incidence.graph().node_count();
        let base: Vec<u64> = (0..node_count)
            .map(|i| cfg.coloring.color(self, incidence.entity(NodeId(i as u32))))
            .collect();

        // Fixpoint: fold each stereo node's partition-dependent observable coset on
        // top of the static base color, re-refine until the orbit partition stabilizes.
        let mut auto = incidence
            .graph()
            .automorphisms(|node| base[node.index()], cfg.automorphism_algorithm);
        let mut orbits = orbit_vec(&auto, node_count);
        let mut iterations: u32 = 0;
        loop {
            iterations += 1;
            let colors: Vec<u64> = (0..node_count)
                .map(|i| self.folded_color(base[i], incidence.entity(NodeId(i as u32)), &orbits))
                .collect();
            let next = incidence
                .graph()
                .automorphisms(|node| colors[node.index()], cfg.automorphism_algorithm);
            let next_orbits = orbit_vec(&next, node_count);
            let stable = next_orbits == orbits;
            auto = next;
            orbits = next_orbits;
            if stable || !cfg.iterate_to_fixpoint || iterations >= cfg.max_iterations {
                break;
            }
        }

        let colors: Vec<u64> = (0..node_count)
            .map(|i| self.folded_color(base[i], incidence.entity(NodeId(i as u32)), &orbits))
            .collect();

        // Grade the converged generators by their action on stereocenter cosets.
        let mut proper: Vec<Vec<NodeId>> = Vec::new();
        let mut improper: Vec<Vec<NodeId>> = Vec::new();
        for generator in auto.generators() {
            match self.grade_generator(
                &incidence,
                generator,
                self.stereo_atoms()
                    .ids()
                    .map(Entity::StereoAtom)
                    .chain(self.stereo_bonds().ids().map(Entity::StereoBond)),
                None,
            ) {
                Some(Orientation::Proper) => proper.push(generator.clone()),
                Some(Orientation::Improper) => improper.push(generator.clone()),
                None => {} // mixed / invalid → discard
            }
        }

        let all: Vec<Vec<NodeId>> = proper.iter().chain(improper.iter()).cloned().collect();
        let proper_orbits = union_find(node_count, &proper);
        let star_orbits = union_find(node_count, &all);
        let chiral = self.has_oriented_center(&incidence) && improper.is_empty();

        GraphSymmetry {
            incidence,
            colors,
            proper_orbits,
            star_orbits,
            chiral,
            automorphism_algorithm: cfg.automorphism_algorithm,
        }
    }

    /// `base` combined with the stereo node's observable coset, or `base` unchanged.
    fn folded_color(&self, base: u64, entity: Entity, orbits: &[NodeId]) -> u64 {
        match self.observable_descriptor(entity, orbits) {
            Some(observable) => {
                let mut hasher = DefaultHasher::new();
                base.hash(&mut hasher);
                observable.hash(&mut hasher);
                hasher.finish()
            }
            None => base,
        }
    }

    /// The observable coset of a stereo node under the current orbit partition:
    /// the stored coset reduced by ligand-position swaps that preserve the ligands'
    /// current classes. `None` for non-stereo nodes or undetermined and
    /// nonliteral configurations, or a ligand count that differs from the kind's degree.
    fn observable_descriptor(&self, entity: Entity, orbits: &[NodeId]) -> Option<u32> {
        let (kind, coset, ligands) = self.stereo_center(entity)?;
        if ligands.len() != kind.degree() {
            return None;
        }
        let &StereoCoset::Lit(raw) = coset else {
            return None;
        };
        let coset_space = kind.class_key().space();
        let classes: Vec<(StereoLigandKind, NodeId)> = ligands
            .iter()
            .map(|ligand| (ligand.kind, orbits[ligand.atom_id.index()]))
            .collect();
        let mut generators = Vec::new();
        for i in 0..classes.len() {
            for j in (i + 1)..classes.len() {
                if classes[i] == classes[j] {
                    generators.push(
                        Permutation::from_cycles(coset_space.degree(), &[vec![i, j]])
                            .expect("class-pair transposition is in range and disjoint"),
                    );
                }
            }
        }
        coset_space.observable_coset(raw, &generators)
    }

    /// Combine a generator's orientation over the supplied stereo carriers.
    /// `None` ⇒ discard (mixed across centers, or inconsistent).
    fn grade_generator(
        &self,
        incidence: &IncidenceGraph,
        generator: &[NodeId],
        carriers: impl Iterator<Item = Entity>,
        mut accumulated: Option<Orientation>,
    ) -> Option<Orientation> {
        for entity in carriers {
            let contribution = match self.grade_center(incidence, entity, generator) {
                Err(()) => return None,
                Ok(None) => continue,
                Ok(Some(orientation)) => orientation,
            };
            match accumulated {
                None => accumulated = Some(contribution),
                Some(previous) if previous == contribution => {}
                Some(_) => return None,
            }
        }
        Some(accumulated.unwrap_or(Orientation::Proper))
    }

    /// One stereocenter's orientation contribution to a generator. `Ok(None)` =
    /// no constraint (non-stereogenic / undetermined); `Err` = generator invalid.
    fn grade_center(
        &self,
        incidence: &IncidenceGraph,
        entity: Entity,
        generator: &[NodeId],
    ) -> Result<Option<Orientation>, ()> {
        let Some((kind, coset, ligands)) = self.stereo_center(entity) else {
            return Ok(None);
        };
        let &StereoCoset::Lit(source) = coset else {
            return Ok(None);
        };
        let coset_space = kind.class_key().space();
        let target_entity = incidence.entity(generator[incidence.node_of(entity).index()]);
        let Some((_, target_coset, target_ligands)) = self.stereo_center(target_entity) else {
            return Err(());
        };
        if !matches!(target_coset, StereoCoset::Lit(_)) {
            return Err(());
        }

        // The generator carries each ligand's bearing atom; the kind is preserved.
        let image: Vec<StereoLigand> = ligands
            .iter()
            .map(|ligand| {
                StereoLigand::new(AtomId::from(generator[ligand.atom_id.index()]), ligand.kind)
            })
            .collect();
        let target_in_image = reexpress(kind, target_coset, &target_ligands, &image).ok_or(())?;

        if source == target_in_image {
            Ok(Some(Orientation::Proper))
        } else if coset_space.enantiomer(target_in_image) == Some(source) {
            Ok(Some(Orientation::Improper))
        } else {
            Err(())
        }
    }

    /// Whether any stereocenter establishes a definite handedness (the precondition
    /// for the molecule to be provably chiral).
    fn has_oriented_center(&self, incidence: &IncidenceGraph) -> bool {
        [EntityKind::StereoAtom, EntityKind::StereoBond]
            .into_iter()
            .flat_map(|kind| (0..incidence.entity_count(kind)).map(move |i| kind.with_id(i as u32)))
            .any(|entity| {
                matches!(
                    self.stereo_center(entity),
                    Some((_, StereoCoset::Lit(_), _))
                )
            })
    }

    /// The kind, stored coset, and ordered ligands of a stereo entity; `None` for
    /// non-stereo entities, out-of-range ids, or an undetermined kind.
    fn stereo_center(
        &self,
        entity: Entity,
    ) -> Option<(StereoKind, &StereoCoset, Vec<StereoLigand>)> {
        match entity {
            Entity::StereoAtom(id) => {
                let view = self.stereo_atoms().get(id)?;
                let StereoConfigurationForm::Kinded(kind, coset) = &view.attributes().configuration
                else {
                    return None;
                };
                let ligands = view
                    .ligands()
                    .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                    .collect();
                Some((*kind, coset, ligands))
            }
            Entity::StereoBond(id) => {
                let view = self.stereo_bonds().get(id)?;
                let StereoConfigurationForm::Kinded(kind, coset) = &view.attributes().configuration
                else {
                    return None;
                };
                let ligands = view
                    .ligands()
                    .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                    .collect();
                Some((*kind, coset, ligands))
            }
            _ => None,
        }
    }
}

impl GraphSymmetry {
    pub fn is_chiral(&self) -> bool {
        self.chiral
    }

    pub fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.proper_orbits[a.index()] == self.proper_orbits[b.index()]
    }

    pub fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.star_orbits[a.index()] == self.star_orbits[b.index()]
    }

    pub fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId> {
        self.orbit_members(&self.proper_orbits, a)
    }

    pub fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId> {
        self.orbit_members(&self.star_orbits, a)
    }

    fn orbit_members(&self, orbits: &[NodeId], a: AtomId) -> Vec<AtomId> {
        let rep = orbits[a.index()];
        (0..self.incidence.entity_count(EntityKind::Atom))
            .filter(|&i| orbits[i] == rep)
            .map(|i| AtomId(i as u32))
            .collect()
    }

    pub(crate) fn incidence(&self) -> &IncidenceGraph {
        &self.incidence
    }

    /// Generators of the stabilizer of a node: one automorphism run with the
    /// site uniquely colored. Ungraded (graph-core node space); the per-carrier
    /// projection grades them.
    pub(crate) fn site_stabilizer(&self, site: NodeId) -> Vec<Vec<NodeId>> {
        self.incidence
            .graph()
            .automorphisms(
                |node| (node == site, self.colors[node.index()]),
                self.automorphism_algorithm,
            )
            .generators()
            .to_vec()
    }
}

/// The local oriented ligand-position symmetry of one stereo carrier, plus its
/// kind and stored coset — the compact input to stereo assertions.
pub struct StereoSymmetry {
    group: OrientedPermutationGroup,
    kind: StereoKind,
    coset: StereoCoset,
}

impl Molecule {
    /// Project the molecule symmetry onto a stereo atom's ligand positions.
    ///
    /// # Panics
    ///
    /// Panics if the id is unavailable or the stereo kind is undetermined.
    pub fn stereo_atom_symmetry(&self, gs: &GraphSymmetry, id: StereoAtomId) -> StereoSymmetry {
        let site = self
            .stereo_atoms()
            .get(id)
            .expect("stereo atom id in range")
            .site_id();
        self.project_stereo(gs, Entity::StereoAtom(id), Entity::Atom(site))
    }

    /// Project the molecule symmetry onto a stereo bond's ligand positions.
    ///
    /// # Panics
    ///
    /// Panics if the id is unavailable or the stereo kind is undetermined.
    pub fn stereo_bond_symmetry(&self, gs: &GraphSymmetry, id: StereoBondId) -> StereoSymmetry {
        let site = self
            .stereo_bonds()
            .get(id)
            .expect("stereo bond id in range")
            .site_id();
        self.project_stereo(gs, Entity::StereoBond(id), Entity::Bond(site))
    }

    fn project_stereo(&self, gs: &GraphSymmetry, carrier: Entity, site: Entity) -> StereoSymmetry {
        let (kind, coset, ligands) = {
            let (k, c, l) = self
                .stereo_center(carrier)
                .expect("carrier is a stereo element");
            (k, c.clone(), l)
        };
        let degree = ligands.len();
        let site_node = gs.incidence().node_of(site);

        let mut oriented: Vec<OrientedPermutation> = Vec::new();
        // The carrier stabilizer, projected onto ligand positions and graded.
        for generator in gs.site_stabilizer(site_node) {
            let Some(permutation) = project_onto_ligands(&ligands, &generator) else {
                continue;
            };
            let local_orientation = match kind {
                StereoKind::Tetrahedral => {
                    let orientation = if permutation.sign() == 1 {
                        Orientation::Proper
                    } else {
                        Orientation::Improper
                    };
                    Some(orientation)
                }
                StereoKind::CisTrans => {
                    if !kind.class_key().space().allows(permutation) {
                        continue;
                    }
                    None
                }
                _ => None,
            };
            let binary = matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans);
            let carriers = self
                .stereo_atoms()
                .ids()
                .map(Entity::StereoAtom)
                .chain(self.stereo_bonds().ids().map(Entity::StereoBond))
                .filter(|&entity| !binary || entity != carrier);
            let Some(orientation) =
                self.grade_generator(gs.incidence(), &generator, carriers, local_orientation)
            else {
                continue;
            };
            oriented.push(oriented_permutation(orientation, permutation));
        }
        StereoSymmetry {
            group: OrientedPermutationGroup::generate(degree, &oriented),
            kind,
            coset,
        }
    }
}

impl StereoSymmetry {
    /// The local oriented ligand-position symmetry group.
    pub fn group(&self) -> &OrientedPermutationGroup {
        &self.group
    }

    /// The carrier's stereo kind.
    pub fn kind(&self) -> StereoKind {
        self.kind
    }

    /// The carrier's stored coset.
    pub fn coset(&self) -> &StereoCoset {
        &self.coset
    }

    /// Whether a literal coset has a singleton orbit under local ligand symmetry.
    ///
    /// For tetrahedral and cis/trans kinds, an undetermined coset or a nonempty
    /// finite coset set is classified from the action on both cosets.
    /// This does not choose or modify the stored coset.
    /// Other nonliteral configurations, empty or invalid finite sets, invalid
    /// literal indices, and incompatible ligand actions return false.
    ///
    /// # Semantic properties
    ///
    /// For a valid tetrahedral or cis/trans ligand action, the result agrees for
    /// an undetermined coset, any valid nonempty finite coset set, and either
    /// literal coset in the same ligand frame.
    pub fn is_stereogenic(&self) -> bool {
        let coset = self.coset.as_lit();
        let binary_domain = matches!(self.kind, StereoKind::Tetrahedral | StereoKind::CisTrans)
            && match &self.coset {
                StereoCoset::Undetermined => true,
                StereoCoset::LitSet(domain) => {
                    !domain.is_empty() && domain.iter().all(|&index| index < 2)
                }
                _ => false,
            };
        if coset.is_none() && !binary_domain {
            return false;
        }
        let perms: Vec<Permutation> = self
            .group
            .elements()
            .iter()
            .map(|op| op.permutation())
            .collect();
        let Some(classes) = self.kind.class_key().space().orbit_reps(&perms) else {
            return false;
        };
        match coset {
            Some(coset) => {
                let Some(&class) = classes.get(coset as usize) else {
                    return false;
                };
                classes.iter().filter(|&&c| c == class).count() == 1
            }
            None => classes.len() == 2 && classes[0] != classes[1],
        }
    }

    /// Stereogenic, prostereogenic (an enantiotopic ligand pair), or nonstereogenic.
    ///
    /// Graph-derived actions fix virtual ligands, so prostereogenicity does not arise
    /// from treating virtual hydrogens as distinguishable positions.
    pub fn stereogenicity(&self) -> Stereogenicity {
        if self.is_stereogenic() {
            Stereogenicity::Stereogenic
        } else if self.has_enantiotopic_pair() {
            Stereogenicity::Prostereogenic
        } else {
            Stereogenicity::Nonstereogenic
        }
    }

    /// The topicity relation between two ligand positions.
    pub fn topicity(&self, a: StereoLigandPosition, b: StereoLigandPosition) -> Topicity {
        if self.group.proper_orbit_of(a.index()).contains(&b.index()) {
            Topicity::Homotopic
        } else if self.group.star_orbit_of(a.index()).contains(&b.index()) {
            Topicity::Enantiotopic
        } else {
            Topicity::Diastereotopic
        }
    }

    fn has_enantiotopic_pair(&self) -> bool {
        let degree = self.kind.degree();
        (0..degree).any(|a| {
            (a + 1..degree).any(|b| {
                self.topicity(
                    StereoLigandPosition(a as u32),
                    StereoLigandPosition(b as u32),
                ) == Topicity::Enantiotopic
            })
        })
    }
}

fn oriented_permutation(orientation: Orientation, permutation: Permutation) -> OrientedPermutation {
    match orientation {
        Orientation::Proper => OrientedPermutation::proper(permutation),
        Orientation::Improper => OrientedPermutation::improper(permutation),
    }
}

/// A stabilizer generator's action on the carrier's ligand positions: atom ligands follow the
/// generator, while distinct virtual ligands sit on the fixed site and stay put. `None` if an atom
/// ligand's image is not a ligand atom of the carrier.
fn project_onto_ligands(ligands: &[StereoLigand], generator: &[NodeId]) -> Option<Permutation> {
    let degree = ligands.len();
    let mut image = vec![0usize; degree];
    for (i, ligand) in ligands.iter().enumerate() {
        let target = match ligand.kind {
            StereoLigandKind::Atom => {
                let atom = generator[ligand.atom_id.index()].index();
                (0..degree).find(|&j| {
                    ligands[j].kind == StereoLigandKind::Atom && ligands[j].atom_id.index() == atom
                })?
            }
            _ => i,
        };
        image[i] = target;
    }
    Permutation::try_from(image.as_slice()).ok()
}

fn orbit_vec(auto: &AutomorphismOutput, node_count: usize) -> Vec<NodeId> {
    (0..node_count)
        .map(|i| auto.orbit_of(NodeId(i as u32)))
        .collect()
}

/// Orbit representative (the component's minimum node) per node, under the group
/// generated by `generators` (each an image map over `0..node_count`).
fn union_find(node_count: usize, generators: &[Vec<NodeId>]) -> Vec<NodeId> {
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut parent: Vec<usize> = (0..node_count).collect();
    for generator in generators {
        for (i, image) in generator.iter().enumerate() {
            let a = find(&mut parent, i);
            let b = find(&mut parent, image.index());
            if a != b {
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    (0..node_count)
        .map(|i| NodeId(find(&mut parent, i) as u32))
        .collect()
}

/// `coset` (in `stored` ligand order) re-expressed in `requested` order, as a coset
/// index. `None` if the orders aren't a relabeling of one identical ligand set.
fn reexpress(
    kind: StereoKind,
    coset: &StereoCoset,
    stored: &[StereoLigand],
    requested: &[StereoLigand],
) -> Option<u32> {
    let permutation = Permutation::between(stored, requested)?;
    match coset_apply_permutation(coset, permutation, kind) {
        Some(StereoCoset::Lit(index)) => Some(index),
        _ => None,
    }
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

/// Whether the remapping `images` of `molecule` onto itself carries every stereo entity onto a
/// stereo entity with the same site and the mirror image of its form.
///
/// `images` holds the image of each node of `incidence_graph`, the incidence graph of `molecule`.
/// A form without a mirror image, as with a coset outside the range of its kind, is not carried.
#[cfg(test)]
fn remapping_mirrors_stereo(
    molecule: &Molecule,
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

    for stereo in molecule.stereo_atoms().iter() {
        let Entity::StereoAtom(image_id) = entity_image(Entity::StereoAtom(stereo.id())) else {
            return false;
        };
        let image = molecule.stereo_atom(image_id);
        if entity_image(Entity::Atom(stereo.site_id())) != Entity::Atom(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        let Some(mirrored) = stereo.attributes().mirror() else {
            return false;
        };
        if !mirrored
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(image.attributes()))
        {
            return false;
        }
    }

    for stereo in molecule.stereo_bonds().iter() {
        let Entity::StereoBond(image_id) = entity_image(Entity::StereoBond(stereo.id())) else {
            return false;
        };
        let image = molecule.stereo_bond(image_id);
        if entity_image(Entity::Bond(stereo.site_id())) != Entity::Bond(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        let Some(mirrored) = stereo.attributes().mirror() else {
            return false;
        };
        if !mirrored
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(image.attributes()))
        {
            return false;
        }
    }

    true
}

/// Whether the remapping `images` of `molecule` onto itself preserves every stereo entity except
/// `entity`, whose coset it changes to `coset`.
///
/// Every stereo entity is carried onto a stereo entity with the same site and the same form,
/// where `entity` has its own coset as a source and `coset` as an image. Its own coset is the
/// stored one, or the first coset of its kind when the stored one is not a literal. `images`
/// holds the image of each node of `incidence_graph`, the incidence graph of `molecule`. A stereo
/// entity without a kind has no coset to change.
#[cfg(test)]
fn remapping_changes_coset(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    images: &[NodeId],
    entity: Entity,
    coset: u32,
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

    for stereo in molecule.stereo_atoms().iter() {
        let source = Entity::StereoAtom(stereo.id());
        let target = entity_image(source);
        let Entity::StereoAtom(image_id) = target else {
            return false;
        };
        let image = molecule.stereo_atom(image_id);
        if entity_image(Entity::Atom(stereo.site_id())) != Entity::Atom(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        let mut form = stereo.attributes().clone();
        if source == entity {
            let Some(kind) = form.configuration.kind() else {
                return false;
            };
            let own = form
                .configuration
                .coset()
                .and_then(AsLit::as_lit)
                .unwrap_or(0);
            form.configuration = StereoConfigurationForm::kinded(kind, own);
        }
        let mut image_form = image.attributes().clone();
        if target == entity {
            let Some(kind) = image_form.configuration.kind() else {
                return false;
            };
            image_form.configuration = StereoConfigurationForm::kinded(kind, coset);
        }
        if !form
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(&image_form))
        {
            return false;
        }
    }

    for stereo in molecule.stereo_bonds().iter() {
        let source = Entity::StereoBond(stereo.id());
        let target = entity_image(source);
        let Entity::StereoBond(image_id) = target else {
            return false;
        };
        let image = molecule.stereo_bond(image_id);
        if entity_image(Entity::Bond(stereo.site_id())) != Entity::Bond(image.site_id()) {
            return false;
        }
        let Some(action) = frame_action(stereo.ligand_ids(), image.ligand_ids()) else {
            return false;
        };
        let mut form = stereo.attributes().clone();
        if source == entity {
            let Some(kind) = form.configuration.kind() else {
                return false;
            };
            let own = form
                .configuration
                .coset()
                .and_then(AsLit::as_lit)
                .unwrap_or(0);
            form.configuration = StereoConfigurationForm::kinded(kind, own);
        }
        let mut image_form = image.attributes().clone();
        if target == entity {
            let Some(kind) = image_form.configuration.kind() else {
                return false;
            };
            image_form.configuration = StereoConfigurationForm::kinded(kind, coset);
        }
        if !form
            .reframe_by(&action)
            .is_some_and(|mapped| mapped.normalized_eq(&image_form))
        {
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

/// The integers of the colors of a normalized molecule at `level`: one per entity, by incidence
/// node, and one per incidence, by incidence edge.
///
/// The colors of each entity kind are ranked in kind order, and the incidences continue after the
/// entities. The color of an incidence is its role and, for an atom of an aromatic system or a
/// multicenter bond, the electron count of that atom. At Topology the incidences are not ranked,
/// since the topology search input does not read them.
pub(crate) fn initial_colors(
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

#[cfg(test)]
pub(crate) mod tests;
