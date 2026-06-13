//! Graph-automorphism symmetry of a molecule, graded into proper vs improper.

use std::collections::HashSet;
use std::hash::{DefaultHasher, Hash, Hasher};

use umol_graph_core::{Automorphism, AutomorphismAlgorithm, NodeId};
use umol_perm::{space, Orientation, OrientedPermutation, OrientedPermutationGroup, Permutation};

use super::coloring::MoleculeColoring;
use super::entity::{Entity, EntityKind};
use super::ids::{AtomId, StereoAtomId, StereoBondId, StereoLigandId};
use super::incidence::{IncidenceGraph, IncidenceNodeSelection};
use super::ligand::{StereoLigand, StereoLigandKind};
use super::molecule::MoleculeAst;
use super::stereo::{Stereogenicity, StereoCosetAst, StereoKind, Topicity};
use super::traits::AsLit;

/// Configuration for [`MoleculeAst::graph_symmetry`].
pub struct GraphSymmetryConfig<C: MoleculeColoring> {
    pub coloring: C,
    pub iterate_to_fixpoint: bool,
    pub max_iterations: usize,
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
}

impl MoleculeAst {
    /// Build the molecule's graded graph-automorphism symmetry under `cfg.coloring`.
    pub fn graph_symmetry<C: MoleculeColoring>(
        &self,
        cfg: &GraphSymmetryConfig<C>,
    ) -> GraphSymmetry {
        let incidence = self.incidence_graph(IncidenceNodeSelection::full());
        let node_count = incidence.graph().node_count();
        let base: Vec<u64> = (0..node_count)
            .map(|i| cfg.coloring.color(self, incidence.entity(NodeId(i as u32))))
            .collect();

        // Fixpoint: fold each stereo node's partition-dependent observable coset on
        // top of the static base color, re-refine until the orbit partition stabilizes.
        let mut auto = incidence
            .graph()
            .automorphisms(|node| base[node.index()], AutomorphismAlgorithm::Nauty);
        let mut orbits = orbit_vec(&auto, node_count);
        let mut iterations = 0;
        loop {
            iterations += 1;
            let colors: Vec<u64> = (0..node_count)
                .map(|i| self.folded_color(base[i], incidence.entity(NodeId(i as u32)), &orbits))
                .collect();
            let next = incidence
                .graph()
                .automorphisms(|node| colors[node.index()], AutomorphismAlgorithm::Nauty);
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
            match self.grade_generator(&incidence, generator) {
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
    /// current classes. `None` for non-stereo / undetermined / malformed nodes.
    fn observable_descriptor(&self, entity: Entity, orbits: &[NodeId]) -> Option<u32> {
        let (kind, coset, ligands) = self.stereo_center(entity)?;
        let &StereoCosetAst::Lit(raw) = coset else {
            return None;
        };
        let coset_space = space(kind.class_key());
        if ligands.len() != coset_space.degree() {
            return None;
        }
        let classes: Vec<(StereoLigandKind, NodeId)> = ligands
            .iter()
            .map(|ligand| (ligand.kind, orbits[ligand.atom_id.index()]))
            .collect();
        let mut generators = Vec::new();
        for i in 0..classes.len() {
            for j in (i + 1)..classes.len() {
                if classes[i] == classes[j] {
                    generators.push(Permutation::from_cycles(coset_space.degree(), &[vec![i, j]]));
                }
            }
        }
        Some(coset_space.observable_coset(raw, &generators))
    }

    /// A generator's orientation = its uniform action on every stereocenter's coset.
    /// `None` ⇒ discard (mixed across centers, or inconsistent).
    fn grade_generator(&self, incidence: &IncidenceGraph, generator: &[NodeId]) -> Option<Orientation> {
        let mut net: Option<Orientation> = None;
        for kind in [EntityKind::StereoAtom, EntityKind::StereoBond] {
            for index in 0..incidence.entity_count(kind) {
                let contribution = match self.grade_center(incidence, kind.with_id(index as u32), generator)
                {
                    Err(()) => return None,
                    Ok(None) => continue,
                    Ok(Some(orientation)) => orientation,
                };
                match net {
                    None => net = Some(contribution),
                    Some(previous) if previous == contribution => {}
                    Some(_) => return None,
                }
            }
        }
        Some(net.unwrap_or(Orientation::Proper))
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
        let &StereoCosetAst::Lit(source) = coset else {
            return Ok(None);
        };
        let coset_space = space(kind.class_key());
        if ligands.len() != coset_space.degree() || !all_distinct(&ligands) {
            return Ok(None);
        }

        let target_entity = incidence.entity(generator[incidence.node_of(entity).index()]);
        let Some((_, target_coset, target_ligands)) = self.stereo_center(target_entity) else {
            return Err(());
        };
        if !matches!(target_coset, StereoCosetAst::Lit(_)) {
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
        } else if source == coset_space.enantiomer(target_in_image) {
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
            .any(|entity| match self.stereo_center(entity) {
                Some((kind, StereoCosetAst::Lit(_), ligands)) => {
                    ligands.len() == space(kind.class_key()).degree() && all_distinct(&ligands)
                }
                _ => false,
            })
    }

    /// The kind, stored coset, and ordered ligands of a stereo entity; `None` for
    /// non-stereo entities or out-of-range ids.
    fn stereo_center(&self, entity: Entity) -> Option<(StereoKind, &StereoCosetAst, Vec<StereoLigand>)> {
        match entity {
            Entity::StereoAtom(id) => {
                let view = self.stereo_atoms().get(id)?;
                let ligands = view
                    .ligands()
                    .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                    .collect();
                Some((view.kind(), view.coset(), ligands))
            }
            Entity::StereoBond(id) => {
                let view = self.stereo_bonds().get(id)?;
                let ligands = view
                    .ligands()
                    .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                    .collect();
                Some((view.kind(), view.coset(), ligands))
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

    /// Generators of the stabilizer of a node: one nauty run with the site uniquely
    /// colored. Ungraded (graph-core node space); the per-carrier projection grades them.
    pub(crate) fn site_stabilizer(&self, site: NodeId) -> Vec<Vec<NodeId>> {
        self.incidence
            .graph()
            .automorphisms(
                |node| (node == site, self.colors[node.index()]),
                AutomorphismAlgorithm::Nauty,
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
    coset: StereoCosetAst,
}

impl MoleculeAst {
    /// Project the molecule symmetry onto a stereo atom's ligand positions.
    pub fn stereo_atom_symmetry(&self, gs: &GraphSymmetry, id: StereoAtomId) -> StereoSymmetry {
        let site = self.stereo_atoms().get(id).expect("stereo atom id in range").site_id();
        self.project_stereo(gs, Entity::StereoAtom(id), Entity::Atom(site))
    }

    /// Project the molecule symmetry onto a stereo bond's ligand positions.
    pub fn stereo_bond_symmetry(&self, gs: &GraphSymmetry, id: StereoBondId) -> StereoSymmetry {
        let site = self.stereo_bonds().get(id).expect("stereo bond id in range").site_id();
        self.project_stereo(gs, Entity::StereoBond(id), Entity::Bond(site))
    }

    fn project_stereo(&self, gs: &GraphSymmetry, carrier: Entity, site: Entity) -> StereoSymmetry {
        let (kind, coset, ligands) = {
            let (k, c, l) = self.stereo_center(carrier).expect("carrier is a stereo element");
            (k, c.clone(), l)
        };
        let degree = ligands.len();
        let site_node = gs.incidence().node_of(site);

        let mut oriented: Vec<OrientedPermutation> = Vec::new();
        // The carrier stabilizer, projected onto ligand positions and graded.
        for generator in gs.site_stabilizer(site_node) {
            let Some(orientation) = self.grade_generator(gs.incidence(), &generator) else {
                continue;
            };
            let Some(perm) = project_onto_ligands(&ligands, &generator) else {
                continue;
            };
            oriented.push(oriented_permutation(orientation, perm));
        }
        // Same-kind virtual ligands are interchangeable; grade each swap locally.
        for swap in virtual_block_swaps(&ligands) {
            if let Some(orientation) = grade_local(kind, &coset, swap) {
                oriented.push(oriented_permutation(orientation, swap));
            }
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
    pub fn coset(&self) -> &StereoCosetAst {
        &self.coset
    }

    /// Whether the stored arrangement is a genuine stereocenter: its coset is not
    /// identified with any other by the local symmetry (proper *or* improper). A
    /// proper ligand symmetry (homotopic ligands) or an improper one (enantiotopic
    /// ligands) both collapse the class and make the center non-stereogenic.
    pub fn is_stereogenic(&self) -> bool {
        let Some(coset) = self.coset.as_lit() else {
            return false;
        };
        let perms: Vec<Permutation> = self.group.elements().iter().map(|op| op.perm()).collect();
        let classes = space(self.kind.class_key()).merge_under(&perms);
        let class = classes[coset as usize];
        classes.iter().filter(|&&c| c == class).count() == 1
    }

    /// The carrier's stereogenicity classification: a genuine stereocenter, a
    /// prochiral center (some enantiotopic ligand pair), or symmetric.
    pub fn stereogenicity(&self) -> Stereogenicity {
        if self.is_stereogenic() {
            Stereogenicity::Stereogenic
        } else if self.has_enantiotopic_pair() {
            Stereogenicity::Prochiral
        } else {
            Stereogenicity::Symmetric
        }
    }

    /// The topicity relation between two ligand positions.
    pub fn topicity(&self, a: StereoLigandId, b: StereoLigandId) -> Topicity {
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
                self.topicity(StereoLigandId(a as u8), StereoLigandId(b as u8))
                    == Topicity::Enantiotopic
            })
        })
    }
}

fn oriented_permutation(orientation: Orientation, perm: Permutation) -> OrientedPermutation {
    match orientation {
        Orientation::Proper => OrientedPermutation::proper(perm),
        Orientation::Improper => OrientedPermutation::improper(perm),
    }
}

/// A stabilizer generator's action on the carrier's ligand positions: atom ligands
/// follow the generator; virtual ligands sit on the (fixed) site, so they stay put
/// (their interchange is supplied separately). `None` if an atom ligand's image
/// isn't a ligand atom of the carrier.
fn project_onto_ligands(ligands: &[StereoLigand], generator: &[NodeId]) -> Option<Permutation> {
    let degree = ligands.len();
    let mut image = vec![0u8; degree];
    let mut used = vec![false; degree];
    for (i, ligand) in ligands.iter().enumerate() {
        let target = match ligand.kind {
            StereoLigandKind::Atom => {
                let atom = generator[ligand.atom_id.index()].index();
                (0..degree).find(|&j| {
                    !used[j]
                        && ligands[j].kind == StereoLigandKind::Atom
                        && ligands[j].atom_id.index() == atom
                })?
            }
            _ => i,
        };
        if used[target] {
            return None;
        }
        image[i] = target as u8;
        used[target] = true;
    }
    Some(Permutation::from_image(degree, &image))
}

/// Adjacent transpositions within each same-kind virtual-ligand block (generating
/// the symmetric group on that block).
fn virtual_block_swaps(ligands: &[StereoLigand]) -> Vec<Permutation> {
    let degree = ligands.len();
    let mut swaps = Vec::new();
    for kind in [StereoLigandKind::ImplicitHydrogen, StereoLigandKind::LonePair] {
        let positions: Vec<usize> = (0..degree).filter(|&i| ligands[i].kind == kind).collect();
        for pair in positions.windows(2) {
            swaps.push(Permutation::from_cycles(degree, &[vec![pair[0], pair[1]]]));
        }
    }
    swaps
}

/// Orientation of a local ligand-position permutation: does it preserve the coset
/// (proper) or send it to its enantiomer (improper)?
fn grade_local(kind: StereoKind, coset: &StereoCosetAst, perm: Permutation) -> Option<Orientation> {
    let index = coset.as_lit()?;
    let coset_space = space(kind.class_key());
    let StereoCosetAst::Lit(transported) =
        StereoCosetAst::Lit(index).apply_permutation(kind, perm)
    else {
        return None;
    };
    if transported == index {
        Some(Orientation::Proper)
    } else if transported == coset_space.enantiomer(index) {
        Some(Orientation::Improper)
    } else {
        None
    }
}

fn orbit_vec(auto: &Automorphism, node_count: usize) -> Vec<NodeId> {
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
    coset: &StereoCosetAst,
    stored: &[StereoLigand],
    requested: &[StereoLigand],
) -> Option<u32> {
    if stored.len() != requested.len() || !all_distinct(stored) {
        return None;
    }
    let stored_set: HashSet<StereoLigand> = stored.iter().copied().collect();
    let requested_set: HashSet<StereoLigand> = requested.iter().copied().collect();
    if stored_set != requested_set {
        return None;
    }
    match coset.apply_permutation(kind, Permutation::between(stored, requested)) {
        StereoCosetAst::Lit(index) => Some(index),
        _ => None,
    }
}

fn all_distinct(ligands: &[StereoLigand]) -> bool {
    ligands.iter().copied().collect::<HashSet<_>>().len() == ligands.len()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::coloring::ConstitutionColoring;
    use crate::ast::constraint::Constraints;
    use crate::ast::ids::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandId};
    use crate::ast::stereo::{StereoAtomAst, StereoBondAst, StereoCosetAst, StereoKind};

    fn config() -> GraphSymmetryConfig<ConstitutionColoring> {
        GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: true,
            max_iterations: 16,
        }
    }

    #[fixture]
    fn benzene_ring() -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); 6];
        let bonds = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    // A tetrahedral center on atom 0 with the four given peripheral elements.
    fn tetrahedral(peripherals: [Element; 4]) -> MoleculeAst {
        let mut atoms = vec![AtomAst::from_element(Element::C)];
        atoms.extend(peripherals.into_iter().map(AtomAst::from_element));
        let bonds = (1..=4)
            .map(|i| (AtomId(0), AtomId(i), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_parts(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
            )],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_molecule_ast_graph_symmetry_constitutional(benzene_ring: MoleculeAst) {
        let symmetry = benzene_ring.graph_symmetry(&config());
        assert!(symmetry.same_proper_orbit(AtomId(0), AtomId(3)));
        assert_eq!(
            symmetry.proper_orbit_of(AtomId(0)),
            vec![
                AtomId(0),
                AtomId(1),
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5),
            ],
        );
        assert!(!symmetry.is_chiral());
    }

    #[rstest]
    fn test_molecule_ast_graph_symmetry_chiral() {
        // Four distinct ligands ⇒ a genuine stereocenter, no compensating symmetry.
        let mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
        let symmetry = mol.graph_symmetry(&config());
        assert!(symmetry.is_chiral());
    }

    #[rstest]
    fn test_molecule_ast_graph_symmetry_enantiotopic() {
        // Two identical Cl ligands: enantiotopic — same star orbit, distinct proper.
        let mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
        let symmetry = mol.graph_symmetry(&config());
        assert!(!symmetry.is_chiral());
        assert!(!symmetry.same_proper_orbit(AtomId(1), AtomId(2)));
        assert!(symmetry.same_star_orbit(AtomId(1), AtomId(2)));
    }

    #[rstest]
    fn test_graph_symmetry_site_stabilizer(benzene_ring: MoleculeAst) {
        let symmetry = benzene_ring.graph_symmetry(&config());
        // Every stabilizer generator fixes the distinguished site.
        let generators = symmetry.site_stabilizer(NodeId(0));
        assert!(generators.iter().all(|g| g[0] == NodeId(0)));
    }

    #[rstest]
    fn test_graph_symmetry_is_chiral_no_stereo(benzene_ring: MoleculeAst) {
        // No stereocenters ⇒ trivially achiral despite no improper generators.
        assert!(!benzene_ring.graph_symmetry(&config()).is_chiral());
    }

    #[rstest]
    fn test_molecule_ast_stereo_atom_symmetry_stereogenic() {
        let mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
        let gs = mol.graph_symmetry(&config());
        let stereo = mol.stereo_atom_symmetry(&gs, StereoAtomId(0));
        assert!(stereo.is_stereogenic());
        // Four distinct ligands ⇒ no local symmetry ⇒ pairwise diastereotopic.
        assert_eq!(
            stereo.topicity(StereoLigandId(0), StereoLigandId(1)),
            Topicity::Diastereotopic
        );
    }

    #[rstest]
    fn test_molecule_ast_stereo_atom_symmetry_prochiral() {
        let mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
        let gs = mol.graph_symmetry(&config());
        let stereo = mol.stereo_atom_symmetry(&gs, StereoAtomId(0));
        // Two identical Cl ⇒ prochiral, not a genuine stereocenter.
        assert!(!stereo.is_stereogenic());
        assert_eq!(
            stereo.topicity(StereoLigandId(0), StereoLigandId(1)),
            Topicity::Enantiotopic
        ); // the two Cl
        assert_eq!(
            stereo.topicity(StereoLigandId(0), StereoLigandId(2)),
            Topicity::Diastereotopic
        ); // Cl vs F
    }

    #[rstest]
    fn test_molecule_ast_stereo_bond_symmetry() {
        // C0=C1 with four distinct substituents (F,Cl on C0; Br,I on C1): E/Z stereogenic.
        let mol = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(2)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
                (AtomId(0), AtomId(3), BondAst::from_order(1)),
                (AtomId(1), AtomId(4), BondAst::from_order(1)),
                (AtomId(1), AtomId(5), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(0)),
            )],
            Constraints::default(),
        );
        let gs = mol.graph_symmetry(&config());
        let stereo = mol.stereo_bond_symmetry(&gs, StereoBondId(0));
        assert!(stereo.is_stereogenic());
        assert_eq!(
            stereo.topicity(StereoLigandId(0), StereoLigandId(1)),
            Topicity::Diastereotopic
        ); // F vs Cl on C0
    }
}
