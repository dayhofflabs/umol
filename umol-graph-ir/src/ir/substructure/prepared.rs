//! Prepared pattern and host carriers for repeated substructure matching.
//!
//! [`PreparedSubstructureHost`] snapshots a host [`Molecule`] alongside its core target bitsets,
//! element/bond/overlay summary masks, and transparent lazy caches for incidence graphs, relevant
//! rings, and proper host symmetry. [`PreparedSubstructurePattern`] snapshots a pattern
//! [`Molecule`], precompiles the query plan, deduplicates atom/bond predicate classes, and builds
//! proper pattern symmetry constraints when configured.

use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::ControlFlow;
use std::sync::{Arc, OnceLock};

use umol_chem::element::MAX_ATOMIC_NUMBER;
use umol_graph_core::{
    AutomorphismAlgorithm, Correspondence, EdgeId, NodeId, PreparedSubgraphQuery,
    PreparedSubgraphTarget, QuerySymmetryConstraints, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismPredicate, SubgraphIsomorphismScratch, TargetSymmetry,
};

use super::super::atom::{AtomForm, ElementForm};
use super::super::bond::BondForm;
use super::super::constraint::{AtomConstraintForm, BondConstraintForm};
use super::super::correspondence::MoleculeCorrespondence;
use super::super::electrons::ElectronCountsForm;
use super::super::entity::{Entity, EntityKind};
use super::super::id::{AtomId, BondId};
use super::super::incidence::{
    Incidence, IncidenceEdge, IncidenceEdgeId, IncidenceGraph, IncidenceLevel, IncidenceNodeId,
};
use super::super::molecule::{GraphSymmetryConfig, Molecule};
use super::super::num::NumForm;
use super::super::ring::{RingConfig, RingModel, RingSet, RingSetKind};
use super::super::traits::Lattice;
use super::{
    atom_fields_match, bond_fields_match, SubstructureMatchAlgorithm, SubstructureMatchConfig,
};

const OVERLAY_DATIVE: u8 = 1 << 0;
const OVERLAY_AROMATIC: u8 = 1 << 1;
const OVERLAY_MULTICENTER: u8 = 1 << 2;
const OVERLAY_NONCOVALENT: u8 = 1 << 3;
const OVERLAY_STEREO_ATOM: u8 = 1 << 4;
const OVERLAY_STEREO_BOND: u8 = 1 << 5;

/// Bitmask covering atomic numbers `1..=MAX_ATOMIC_NUMBER`.
const ALL_ELEMENT_BITS: u128 = ((1u128 << (MAX_ATOMIC_NUMBER as u32 + 1)) - 1) ^ 1u128;

/// Maximum literal bond order accelerated by compact 8-bit order masks (`0..=7`).
const MAX_FAST_BOND_ORDER: i64 = 7;

struct SubstructureScratch {
    subiso: SubgraphIsomorphismScratch,
    atoms: Vec<AtomId>,
}

impl SubstructureScratch {
    fn new() -> Self {
        Self {
            subiso: SubgraphIsomorphismScratch::new(),
            atoms: Vec::new(),
        }
    }
}

thread_local! {
    static SUBSTRUCTURE_SCRATCH: RefCell<SubstructureScratch> =
        RefCell::new(SubstructureScratch::new());
}

fn with_substructure_scratch<R>(f: impl FnOnce(&mut SubstructureScratch) -> R) -> R {
    SUBSTRUCTURE_SCRATCH.with(|cell| {
        if let Ok(mut guard) = cell.try_borrow_mut() {
            f(&mut guard)
        } else {
            let mut fallback = SubstructureScratch::new();
            f(&mut fallback)
        }
    })
}

#[inline]
fn bit_rank_u128(mask: u128, bit: u8) -> Option<usize> {
    let bit_mask = 1u128 << bit;
    if (mask & bit_mask) == 0 {
        return None;
    }
    Some((mask & (bit_mask - 1)).count_ones() as usize)
}

#[inline]
fn bit_rank_u8(mask: u8, bit: u8) -> Option<usize> {
    let bit_mask = 1u8 << bit;
    if (mask & bit_mask) == 0 {
        return None;
    }
    Some((mask & (bit_mask - 1)).count_ones() as usize)
}

#[inline]
fn range_mask_u64(start: usize, end: usize) -> u64 {
    let len = end.saturating_sub(start);
    if len == 0 || start >= 64 {
        return 0;
    }
    let clamped = len.min(64 - start);
    let bits = if clamped == 64 {
        u64::MAX
    } else {
        (1u64 << clamped) - 1
    };
    bits << start
}

#[inline]
fn range_mask_u128(start: usize, end: usize) -> u128 {
    let len = end.saturating_sub(start);
    if len == 0 || start >= 128 {
        return 0;
    }
    let clamped = len.min(128 - start);
    let bits = if clamped == 128 {
        u128::MAX
    } else {
        (1u128 << clamped) - 1
    };
    bits << start
}

#[inline]
fn range_mask_word(start: usize, end: usize, word_idx: usize) -> u64 {
    let word_start = word_idx * 64;
    let word_end = word_start + 64;
    let lo = start.max(word_start);
    let hi = end.min(word_end);
    if lo >= hi {
        return 0;
    }
    range_mask_u64(lo - word_start, hi - word_start)
}

fn overlay_summary(molecule: &Molecule) -> (u8, [u32; 6]) {
    let counts = [
        molecule.dative_bonds().count() as u32,
        molecule.aromatic_systems().count() as u32,
        molecule.multicenter_bonds().count() as u32,
        molecule.noncovalent_bonds().count() as u32,
        molecule.stereo_atoms().count() as u32,
        molecule.stereo_bonds().count() as u32,
    ];
    let mut present = 0u8;
    if counts[0] > 0 {
        present |= OVERLAY_DATIVE;
    }
    if counts[1] > 0 {
        present |= OVERLAY_AROMATIC;
    }
    if counts[2] > 0 {
        present |= OVERLAY_MULTICENTER;
    }
    if counts[3] > 0 {
        present |= OVERLAY_NONCOVALENT;
    }
    if counts[4] > 0 {
        present |= OVERLAY_STEREO_ATOM;
    }
    if counts[5] > 0 {
        present |= OVERLAY_STEREO_BOND;
    }
    (present, counts)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum IncidenceElectron {
    NotApplicable,
    Undetermined,
    Lit(i64),
    Malformed,
}

fn precompute_incidence_electrons(
    molecule: &Molecule,
    incidence: &IncidenceGraph,
) -> Vec<IncidenceElectron> {
    incidence
        .incidences()
        .map(|(edge, role)| {
            if !matches!(role, Incidence::AromaticAtom | Incidence::MulticenterAtom) {
                return IncidenceElectron::NotApplicable;
            }
            let (owner, atom) = incidence
                .incidence_endpoints(edge)
                .expect("overlay incidence connects an overlay to an atom");
            match owner {
                Entity::AromaticSystem(id) => {
                    let system = molecule.aromatic_system(id);
                    match system.electrons() {
                        ElectronCountsForm::Undetermined => IncidenceElectron::Undetermined,
                        ElectronCountsForm::Lit(counts)
                            if counts.len() == system.atom_ids().len() =>
                        {
                            system
                                .atom_ids()
                                .position(|member| member == atom)
                                .map_or(IncidenceElectron::Malformed, |position| {
                                    IncidenceElectron::Lit(counts[position])
                                })
                        }
                        ElectronCountsForm::Lit(_) => IncidenceElectron::Malformed,
                    }
                }
                Entity::MulticenterBond(id) => {
                    let bond = molecule.multicenter_bond(id);
                    match bond.electrons() {
                        ElectronCountsForm::Undetermined => IncidenceElectron::Undetermined,
                        ElectronCountsForm::Lit(counts)
                            if counts.len() == bond.atom_ids().len() =>
                        {
                            bond.atom_ids()
                                .position(|member| member == atom)
                                .map_or(IncidenceElectron::Malformed, |position| {
                                    IncidenceElectron::Lit(counts[position])
                                })
                        }
                        ElectronCountsForm::Lit(_) => IncidenceElectron::Malformed,
                    }
                }
                _ => {
                    unreachable!(
                        "electron contributions belong to aromatic or multicenter entities"
                    )
                }
            }
        })
        .collect()
}

#[inline]
fn incidence_electrons_match(pattern: IncidenceElectron, host: IncidenceElectron) -> bool {
    match (pattern, host) {
        (IncidenceElectron::NotApplicable, IncidenceElectron::NotApplicable) => true,
        (
            IncidenceElectron::Undetermined,
            IncidenceElectron::Undetermined | IncidenceElectron::Lit(_),
        ) => true,
        (IncidenceElectron::Lit(p), IncidenceElectron::Lit(h)) => p == h,
        _ => false,
    }
}

#[derive(Clone, Debug)]
enum HostElementBitsets {
    U64 {
        by_element: Vec<u64>,
        non_lit: u64,
    },
    U128 {
        by_element: Vec<u128>,
        non_lit: u128,
    },
    Wide {
        words: usize,
        by_element: Vec<u64>,
        non_lit: Vec<u64>,
    },
}

#[derive(Clone, Debug)]
enum HostBondOrderBitsets {
    None,
    U64 { order_count: usize, adj: Vec<u64> },
    U128 { order_count: usize, adj: Vec<u128> },
}

#[derive(Debug)]
struct PreparedIncidenceHost {
    incidence: IncidenceGraph,
    target: PreparedSubgraphTarget,
    /// Electron contributions indexed by `incidence.overlay_incidence_index(edge)`.
    incidence_electrons: Vec<IncidenceElectron>,
}

impl PreparedIncidenceHost {
    fn new(molecule: &Molecule) -> Self {
        let incidence = molecule.incidence_graph(IncidenceLevel::Constitution);
        let target = PreparedSubgraphTarget::new(incidence.raw_graph());
        let incidence_electrons = precompute_incidence_electrons(molecule, &incidence);
        Self {
            incidence,
            target,
            incidence_electrons,
        }
    }

    #[inline]
    fn kind_range(&self, kind: EntityKind) -> (usize, usize) {
        self.incidence.node_range(kind)
    }
}

#[derive(Debug, Default)]
struct HostSymmetryCache {
    atom_symmetry: OnceLock<TargetSymmetry>,
    incidence_symmetry: OnceLock<TargetSymmetry>,
}

impl HostSymmetryCache {
    fn atom_symmetry(
        &self,
        molecule: &Molecule,
        algorithm: AutomorphismAlgorithm,
    ) -> &TargetSymmetry {
        self.atom_symmetry.get_or_init(|| {
            let generators = molecule.proper_atom_generators(algorithm, true);
            TargetSymmetry::try_from_generators(molecule.raw_graph(), &generators)
                .expect("proper molecule symmetry generators are valid graph automorphisms")
        })
    }

    fn incidence_symmetry(
        &self,
        molecule: &Molecule,
        algorithm: AutomorphismAlgorithm,
        incidence: &PreparedIncidenceHost,
    ) -> &TargetSymmetry {
        self.incidence_symmetry.get_or_init(|| {
            let symmetry = molecule
                .graph_symmetry(&GraphSymmetryConfig {
                    automorphism_algorithm: algorithm,
                    include_constraints: true,
                    ..GraphSymmetryConfig::default()
                })
                .into_graph_symmetry();
            let generators = molecule.proper_incidence_generators(&symmetry);
            TargetSymmetry::try_from_generators(incidence.incidence.raw_graph(), &generators)
                .expect("proper incidence symmetry generators are valid graph automorphisms")
        })
    }

    #[cfg(test)]
    fn is_populated(&self) -> bool {
        self.atom_symmetry.get().is_some() || self.incidence_symmetry.get().is_some()
    }
}

#[derive(Debug, Default)]
struct HostCaches {
    incidence: OnceLock<PreparedIncidenceHost>,
    rings_vismara: OnceLock<RingSet>,
    symmetry_trail_ir: HostSymmetryCache,
    symmetry_nauty: HostSymmetryCache,
}

/// Immutable prepared host molecule for repeated substructure matching.
///
/// Holds an `Arc`-shared snapshot of the host [`Molecule`] alongside its core target bitsets,
/// element/bond/overlay summary masks, and transparent thread-safe lazy caches for its
/// constitution incidence graph, relevant rings (through size 22), and proper stereo-preserving
/// host symmetry per [`AutomorphismAlgorithm`].
///
/// # Semantic properties
///
/// - Matching any [`PreparedSubstructurePattern`] against `PreparedSubstructureHost::new(host)`
///   produces the exact same match set, existence, and count as matching directly against `host`.
/// - Cloning shares the immutable snapshot and lazy caches in `O(1)`; mutating the source
///   [`Molecule`] after construction does not affect the prepared host.
/// - Lazy caches are observably pure and safe to share concurrently across threads.
#[derive(Clone, Debug)]
pub struct PreparedSubstructureHost {
    molecule: Molecule,
    target: PreparedSubgraphTarget,
    present_elements: u128,
    all_atoms_lit_element: bool,
    all_atoms_fields_satisfiable: bool,
    element_masks: HostElementBitsets,
    present_bond_orders: u8,
    all_bonds_fast_lit_order: bool,
    all_bonds_fields_satisfiable: bool,
    bond_orders: Vec<u8>,
    bond_order_adj: HostBondOrderBitsets,
    present_overlays: u8,
    overlay_counts: [u32; 6],
    caches: Arc<HostCaches>,
}

impl PreparedSubstructureHost {
    /// Prepare `host` for repeated substructure matching.
    pub fn new(host: &Molecule) -> Self {
        let molecule = host.clone();
        let target = PreparedSubgraphTarget::new(molecule.raw_graph());
        let atom_count = molecule.atoms().count();
        let bond_count = molecule.bonds().count();

        let default_atom = AtomForm::default();
        let mut present_elements = 0u128;
        let mut all_atoms_lit_element = true;
        let mut all_atoms_fields_satisfiable = true;
        for atom in molecule.atoms().iter() {
            let form = atom.attributes();
            if let ElementForm::Lit(element) = form.element {
                present_elements |= 1u128 << element.atomic_number();
            } else {
                all_atoms_lit_element = false;
            }
            if all_atoms_fields_satisfiable && !atom_fields_match(&default_atom, form) {
                all_atoms_fields_satisfiable = false;
            }
        }

        let elem_count = present_elements.count_ones() as usize;
        let element_masks = if atom_count <= 64 {
            let mut by_element = vec![0u64; elem_count];
            let mut non_lit = 0u64;
            for atom in molecule.atoms().iter() {
                let idx = atom.id().index();
                if let ElementForm::Lit(element) = atom.attributes().element {
                    let rank = bit_rank_u128(present_elements, element.atomic_number())
                        .expect("present element bit is set");
                    by_element[rank] |= 1u64 << idx;
                } else {
                    non_lit |= 1u64 << idx;
                }
            }
            HostElementBitsets::U64 {
                by_element,
                non_lit,
            }
        } else if atom_count <= 128 {
            let mut by_element = vec![0u128; elem_count];
            let mut non_lit = 0u128;
            for atom in molecule.atoms().iter() {
                let idx = atom.id().index();
                if let ElementForm::Lit(element) = atom.attributes().element {
                    let rank = bit_rank_u128(present_elements, element.atomic_number())
                        .expect("present element bit is set");
                    by_element[rank] |= 1u128 << idx;
                } else {
                    non_lit |= 1u128 << idx;
                }
            }
            HostElementBitsets::U128 {
                by_element,
                non_lit,
            }
        } else {
            let words = atom_count.div_ceil(64);
            let mut by_element = vec![0u64; elem_count * words];
            let mut non_lit = if all_atoms_lit_element {
                Vec::new()
            } else {
                vec![0u64; words]
            };
            for atom in molecule.atoms().iter() {
                let idx = atom.id().index();
                let word = idx / 64;
                let bit = 1u64 << (idx % 64);
                if let ElementForm::Lit(element) = atom.attributes().element {
                    let rank = bit_rank_u128(present_elements, element.atomic_number())
                        .expect("present element bit is set");
                    by_element[rank * words + word] |= bit;
                } else {
                    non_lit[word] |= bit;
                }
            }
            HostElementBitsets::Wide {
                words,
                by_element,
                non_lit,
            }
        };

        let default_bond = BondForm::default();
        let mut present_bond_orders = 0u8;
        let mut all_bonds_fast_lit_order = true;
        let mut all_bonds_fields_satisfiable = true;
        for bond in molecule.bonds().iter() {
            let form = bond.attributes();
            match form.order {
                NumForm::Lit(order) if (0..=MAX_FAST_BOND_ORDER).contains(&order) => {
                    present_bond_orders |= 1u8 << (order as u8);
                }
                _ => {
                    all_bonds_fast_lit_order = false;
                }
            }
            if all_bonds_fields_satisfiable && !bond_fields_match(&default_bond, form) {
                all_bonds_fields_satisfiable = false;
            }
        }

        let bond_orders = if all_bonds_fast_lit_order && bond_count > 0 {
            molecule
                .bonds()
                .iter()
                .map(|bond| match bond.attributes().order {
                    NumForm::Lit(order) => order as u8,
                    _ => unreachable!("validated all bonds have fast literal orders"),
                })
                .collect()
        } else {
            Vec::new()
        };

        let order_count = present_bond_orders.count_ones() as usize;
        let bond_order_adj = if all_bonds_fast_lit_order
            && all_bonds_fields_satisfiable
            && order_count > 0
            && atom_count <= 64
        {
            let mut adj = vec![0u64; atom_count * order_count];
            for bond in molecule.bonds().iter() {
                let [first, second] = bond.atom_ids();
                let order = bond_orders[bond.id().index()];
                let rank =
                    bit_rank_u8(present_bond_orders, order).expect("present bond order bit is set");
                adj[first.index() * order_count + rank] |= 1u64 << second.index();
                adj[second.index() * order_count + rank] |= 1u64 << first.index();
            }
            HostBondOrderBitsets::U64 { order_count, adj }
        } else if all_bonds_fast_lit_order
            && all_bonds_fields_satisfiable
            && order_count > 0
            && atom_count <= 128
        {
            let mut adj = vec![0u128; atom_count * order_count];
            for bond in molecule.bonds().iter() {
                let [first, second] = bond.atom_ids();
                let order = bond_orders[bond.id().index()];
                let rank =
                    bit_rank_u8(present_bond_orders, order).expect("present bond order bit is set");
                adj[first.index() * order_count + rank] |= 1u128 << second.index();
                adj[second.index() * order_count + rank] |= 1u128 << first.index();
            }
            HostBondOrderBitsets::U128 { order_count, adj }
        } else {
            HostBondOrderBitsets::None
        };

        let (present_overlays, overlay_counts) = overlay_summary(&molecule);

        Self {
            molecule,
            target,
            present_elements,
            all_atoms_lit_element,
            all_atoms_fields_satisfiable,
            element_masks,
            present_bond_orders,
            all_bonds_fast_lit_order,
            all_bonds_fields_satisfiable,
            bond_orders,
            bond_order_adj,
            present_overlays,
            overlay_counts,
            caches: Arc::new(HostCaches::default()),
        }
    }

    /// Snapshotted host [`Molecule`].
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    fn incidence(&self) -> &PreparedIncidenceHost {
        self.caches
            .incidence
            .get_or_init(|| PreparedIncidenceHost::new(&self.molecule))
    }

    fn rings(&self, algorithm: RelevantCycleEnumerationAlgorithm) -> &RingSet {
        match algorithm {
            RelevantCycleEnumerationAlgorithm::Vismara => {
                self.caches.rings_vismara.get_or_init(|| {
                    self.molecule
                        .rings(
                            RingModel {
                                kind: RingSetKind::Relevant,
                                max_ring_size: 22,
                            },
                            RingConfig {
                                relevant_cycle_algorithm: algorithm,
                                ..RingConfig::default()
                            },
                        )
                        .into_ring_set()
                })
            }
        }
    }

    fn symmetry_cache(&self, algorithm: AutomorphismAlgorithm) -> &HostSymmetryCache {
        match algorithm {
            AutomorphismAlgorithm::TrailIr => &self.caches.symmetry_trail_ir,
            AutomorphismAlgorithm::Nauty => &self.caches.symmetry_nauty,
        }
    }

    fn atom_symmetry(&self, algorithm: AutomorphismAlgorithm) -> Option<&TargetSymmetry> {
        let sym = self
            .symmetry_cache(algorithm)
            .atom_symmetry(&self.molecule, algorithm);
        (!sym.is_trivial()).then_some(sym)
    }

    fn incidence_symmetry(
        &self,
        algorithm: AutomorphismAlgorithm,
        incidence: &PreparedIncidenceHost,
    ) -> Option<&TargetSymmetry> {
        let sym =
            self.symmetry_cache(algorithm)
                .incidence_symmetry(&self.molecule, algorithm, incidence);
        (!sym.is_trivial()).then_some(sym)
    }

    #[inline]
    fn element_mask_u64(&self, atomic_number: u8) -> u64 {
        let Some(rank) = bit_rank_u128(self.present_elements, atomic_number) else {
            return 0;
        };
        match &self.element_masks {
            HostElementBitsets::U64 { by_element, .. } => by_element[rank],
            HostElementBitsets::U128 { .. } | HostElementBitsets::Wide { .. } => {
                unreachable!("u64 mask queried only when atom_count <= 64")
            }
        }
    }

    #[inline]
    fn element_mask_u128(&self, atomic_number: u8) -> u128 {
        let Some(rank) = bit_rank_u128(self.present_elements, atomic_number) else {
            return 0;
        };
        match &self.element_masks {
            HostElementBitsets::U64 { by_element, .. } => u128::from(by_element[rank]),
            HostElementBitsets::U128 { by_element, .. } => by_element[rank],
            HostElementBitsets::Wide { .. } => {
                unreachable!("u128 mask queried only when atom_count <= 128")
            }
        }
    }

    #[inline]
    fn element_set_mask_u64(&self, allowed: u128) -> u64 {
        let mut active = allowed & self.present_elements;
        let mut mask = 0u64;
        while active != 0 {
            let z = active.trailing_zeros() as u8;
            active &= active - 1;
            mask |= self.element_mask_u64(z);
        }
        mask
    }

    #[inline]
    fn element_set_mask_u128(&self, allowed: u128) -> u128 {
        let mut active = allowed & self.present_elements;
        let mut mask = 0u128;
        while active != 0 {
            let z = active.trailing_zeros() as u8;
            active &= active - 1;
            mask |= self.element_mask_u128(z);
        }
        mask
    }

    #[inline]
    fn non_lit_u64(&self) -> u64 {
        match &self.element_masks {
            HostElementBitsets::U64 { non_lit, .. } => *non_lit,
            HostElementBitsets::U128 { .. } | HostElementBitsets::Wide { .. } => {
                unreachable!("u64 non_lit queried only when atom_count <= 64")
            }
        }
    }

    #[inline]
    fn non_lit_u128(&self) -> u128 {
        match &self.element_masks {
            HostElementBitsets::U64 { non_lit, .. } => u128::from(*non_lit),
            HostElementBitsets::U128 { non_lit, .. } => *non_lit,
            HostElementBitsets::Wide { .. } => {
                unreachable!("u128 non_lit queried only when atom_count <= 128")
            }
        }
    }

    #[inline]
    fn element_word(&self, atomic_number: u8, word_idx: usize) -> u64 {
        let Some(rank) = bit_rank_u128(self.present_elements, atomic_number) else {
            return 0;
        };
        match &self.element_masks {
            HostElementBitsets::U64 { by_element, .. } => {
                if word_idx == 0 {
                    by_element[rank]
                } else {
                    0
                }
            }
            HostElementBitsets::U128 { by_element, .. } => match word_idx {
                0 => by_element[rank] as u64,
                1 => (by_element[rank] >> 64) as u64,
                _ => 0,
            },
            HostElementBitsets::Wide {
                words, by_element, ..
            } => {
                if word_idx < *words {
                    by_element[rank * *words + word_idx]
                } else {
                    0
                }
            }
        }
    }

    #[inline]
    fn element_set_word(&self, allowed: u128, word_idx: usize) -> u64 {
        let mut active = allowed & self.present_elements;
        let mut mask = 0u64;
        while active != 0 {
            let z = active.trailing_zeros() as u8;
            active &= active - 1;
            mask |= self.element_word(z, word_idx);
        }
        mask
    }

    #[inline]
    fn non_lit_word(&self, word_idx: usize) -> u64 {
        match &self.element_masks {
            HostElementBitsets::U64 { non_lit, .. } => {
                if word_idx == 0 {
                    *non_lit
                } else {
                    0
                }
            }
            HostElementBitsets::U128 { non_lit, .. } => match word_idx {
                0 => *non_lit as u64,
                1 => (*non_lit >> 64) as u64,
                _ => 0,
            },
            HostElementBitsets::Wide { non_lit, .. } => {
                if word_idx < non_lit.len() {
                    non_lit[word_idx]
                } else {
                    0
                }
            }
        }
    }

    #[inline]
    fn bond_order_neighbors_u64(&self, orders_mask: u8, target_node: NodeId) -> Option<u64> {
        if target_node.index() >= self.molecule.atoms().count() {
            return None;
        }
        let HostBondOrderBitsets::U64 { order_count, adj } = &self.bond_order_adj else {
            return None;
        };
        let mut active = orders_mask & self.present_bond_orders;
        if active == 0 {
            return Some(0);
        }
        let base = target_node.index() * *order_count;
        let mut mask = 0u64;
        while active != 0 {
            let order = active.trailing_zeros() as u8;
            active &= active - 1;
            let rank = bit_rank_u8(self.present_bond_orders, order)
                .expect("present bond order bit is set");
            mask |= adj[base + rank];
        }
        Some(mask)
    }

    #[inline]
    fn bond_order_neighbors_u128(&self, orders_mask: u8, target_node: NodeId) -> Option<u128> {
        if target_node.index() >= self.molecule.atoms().count() {
            return None;
        }
        match &self.bond_order_adj {
            HostBondOrderBitsets::None => None,
            HostBondOrderBitsets::U64 { order_count, adj } => {
                let mut active = orders_mask & self.present_bond_orders;
                if active == 0 {
                    return Some(0);
                }
                let base = target_node.index() * *order_count;
                let mut mask = 0u128;
                while active != 0 {
                    let order = active.trailing_zeros() as u8;
                    active &= active - 1;
                    let rank = bit_rank_u8(self.present_bond_orders, order)
                        .expect("present bond order bit is set");
                    mask |= u128::from(adj[base + rank]);
                }
                Some(mask)
            }
            HostBondOrderBitsets::U128 { order_count, adj } => {
                let mut active = orders_mask & self.present_bond_orders;
                if active == 0 {
                    return Some(0);
                }
                let base = target_node.index() * *order_count;
                let mut mask = 0u128;
                while active != 0 {
                    let order = active.trailing_zeros() as u8;
                    active &= active - 1;
                    let rank = bit_rank_u8(self.present_bond_orders, order)
                        .expect("present bond order bit is set");
                    mask |= adj[base + rank];
                }
                Some(mask)
            }
        }
    }

    #[cfg(test)]
    pub(crate) fn has_cached_incidence(&self) -> bool {
        self.caches.incidence.get().is_some()
    }

    #[cfg(test)]
    pub(crate) fn has_cached_rings(&self, algorithm: RelevantCycleEnumerationAlgorithm) -> bool {
        match algorithm {
            RelevantCycleEnumerationAlgorithm::Vismara => self.caches.rings_vismara.get().is_some(),
        }
    }

    #[cfg(test)]
    pub(crate) fn has_cached_symmetry(&self, algorithm: AutomorphismAlgorithm) -> bool {
        self.symmetry_cache(algorithm).is_populated()
    }
}

#[derive(Clone, Copy, Debug)]
enum AtomClassFastKind {
    Any,
    LitElementOnly(u8),
    ElementSetOnly(u128),
    Filtered { element_mask: Option<u128> },
}

#[derive(Clone, Debug)]
struct PatternAtomClass {
    representative: AtomId,
    fast_kind: AtomClassFastKind,
    constraint_slot: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
enum BondFastKind {
    Undetermined,
    OrderMask(u8),
    General,
}

#[derive(Clone, Debug)]
struct PatternBondClass {
    representative: BondId,
    fast_kind: BondFastKind,
    constraint_slot: Option<usize>,
}

#[derive(Clone, Debug)]
enum PreparedPatternEngine {
    GraphAndOverlays {
        query: PreparedSubgraphQuery,
    },
    Incidence {
        incidence: IncidenceGraph,
        query: PreparedSubgraphQuery,
        incidence_electrons: Vec<IncidenceElectron>,
        node_classes: Vec<u32>,
        node_class_count: usize,
    },
}

/// Immutable prepared pattern molecule for repeated substructure matching.
///
/// Validates molecule-scope constraints once at construction time, snapshots the pattern
/// [`Molecule`], deduplicates atom and bond predicate classes, precomputes element/bond/overlay
/// screening masks, compiles the selected strategy's [`PreparedSubgraphQuery`], and computes
/// proper stereo-, overlay-, and constraint-preserving [`QuerySymmetryConstraints`] when
/// `config.reduce_pattern_symmetry` is enabled.
///
/// # Semantic properties
///
/// - Matching `PreparedSubstructurePattern::new(pattern, config)` against `host` (or
///   `PreparedSubstructureHost::new(host)`) produces the exact same match set, existence, and
///   count as matching `pattern` directly with `config`.
/// - Cloning preserves the prepared snapshot and configuration; mutating the source [`Molecule`]
///   after construction does not affect the prepared pattern.
#[derive(Clone, Debug)]
pub struct PreparedSubstructurePattern {
    pattern: Molecule,
    config: SubstructureMatchConfig,
    needs_rings: bool,
    required_elements: u128,
    required_bond_orders: u8,
    present_overlays: u8,
    overlay_counts: [u32; 6],
    atom_classes: Vec<PatternAtomClass>,
    atom_class_of: Vec<u32>,
    constrained_atom_class_count: usize,
    bond_classes: Vec<PatternBondClass>,
    bond_class_of: Vec<u32>,
    constrained_bond_class_count: usize,
    all_bonds_undetermined: bool,
    engine: PreparedPatternEngine,
}

impl PreparedSubstructurePattern {
    /// Prepare `pattern` for substructure matching under `config`.
    ///
    /// # Panics
    ///
    /// Panics if `pattern` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn new(pattern: &Molecule, config: SubstructureMatchConfig) -> Self {
        Self::with_pattern_symmetry(pattern, config, |pattern, incidence| match incidence {
            None => pattern.proper_atom_generators(config.automorphism_algorithm, true),
            Some(_) => {
                let graph_symmetry = pattern
                    .graph_symmetry(&GraphSymmetryConfig {
                        automorphism_algorithm: config.automorphism_algorithm,
                        include_constraints: true,
                        ..GraphSymmetryConfig::default()
                    })
                    .into_graph_symmetry();
                pattern.proper_incidence_generators(&graph_symmetry)
            }
        })
    }

    /// Prepare `pattern` for substructure matching under `config`, supplying custom proper
    /// pattern-symmetry generators when `config.reduce_pattern_symmetry` is enabled.
    ///
    /// `pattern_symmetry_generators` is invoked only when `config.reduce_pattern_symmetry` is
    /// `true`, receiving `None` for `GraphAndOverlays` (generators over `pattern.raw_graph()`) or
    /// `Some(&incidence)` for `Incidence` (generators over `incidence.raw_graph()`).
    ///
    /// # Panics
    ///
    /// Panics if `pattern` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub(crate) fn with_pattern_symmetry(
        pattern: &Molecule,
        config: SubstructureMatchConfig,
        pattern_symmetry_generators: impl FnOnce(&Molecule, Option<&IncidenceGraph>) -> Vec<Vec<NodeId>>,
    ) -> Self {
        assert!(
            pattern.constraints().is_empty(),
            "substructure matching does not evaluate molecule constraints"
        );
        let pattern = pattern.clone();
        let needs_rings = pattern_needs_rings(&pattern);
        let (present_overlays, overlay_counts) = overlay_summary(&pattern);

        // Deduplicate atom predicates into exact equivalence classes.
        let mut atom_class_map: HashMap<&AtomForm, u32> = HashMap::new();
        let mut atom_classes = Vec::new();
        let mut atom_class_of = Vec::with_capacity(pattern.atoms().count());
        let mut constrained_atom_class_count = 0usize;
        let mut required_elements = 0u128;

        for atom in pattern.atoms().iter() {
            let form = atom.attributes();
            if let ElementForm::Lit(element) = form.element {
                required_elements |= 1u128 << element.atomic_number();
            }
            let class_id = if let Some(&existing) = atom_class_map.get(form) {
                existing
            } else {
                let id = atom_classes.len() as u32;
                let constraint_slot = (!form.constraints.is_empty()).then(|| {
                    let slot = constrained_atom_class_count;
                    constrained_atom_class_count += 1;
                    slot
                });
                let fast_kind = classify_atom_form(form);
                atom_classes.push(PatternAtomClass {
                    representative: atom.id(),
                    fast_kind,
                    constraint_slot,
                });
                atom_class_map.insert(form, id);
                id
            };
            atom_class_of.push(class_id);
        }

        // Deduplicate bond predicates into exact equivalence classes.
        let mut bond_class_map: HashMap<&BondForm, u32> = HashMap::new();
        let mut bond_classes = Vec::new();
        let mut bond_class_of = Vec::with_capacity(pattern.bonds().count());
        let mut constrained_bond_class_count = 0usize;
        let mut required_bond_orders = 0u8;
        let mut all_bonds_undetermined = true;

        for bond in pattern.bonds().iter() {
            let form = bond.attributes();
            if let NumForm::Lit(order) = form.order {
                if (0..=MAX_FAST_BOND_ORDER).contains(&order) {
                    required_bond_orders |= 1u8 << (order as u8);
                }
            }
            let class_id = if let Some(&existing) = bond_class_map.get(form) {
                existing
            } else {
                let id = bond_classes.len() as u32;
                let constraint_slot = (!form.constraints.is_empty()).then(|| {
                    let slot = constrained_bond_class_count;
                    constrained_bond_class_count += 1;
                    slot
                });
                let fast_kind = classify_bond_form(form);
                bond_classes.push(PatternBondClass {
                    representative: bond.id(),
                    fast_kind,
                    constraint_slot,
                });
                bond_class_map.insert(form, id);
                id
            };
            if !matches!(
                bond_classes[class_id as usize].fast_kind,
                BondFastKind::Undetermined
            ) {
                all_bonds_undetermined = false;
            }
            bond_class_of.push(class_id);
        }

        let engine = match config.match_algorithm {
            SubstructureMatchAlgorithm::GraphAndOverlays => {
                let mut query = PreparedSubgraphQuery::new(
                    pattern.raw_graph(),
                    config.subgraph_isomorphism_algorithm,
                );
                if config.reduce_pattern_symmetry {
                    let generators = pattern_symmetry_generators(&pattern, None);
                    if !generators.is_empty() {
                        let symmetry = QuerySymmetryConstraints::try_from_generators(
                            pattern.raw_graph(),
                            &generators,
                        )
                        .expect("proper pattern symmetry generators are valid graph automorphisms");
                        query = query.with_symmetry(symmetry);
                    }
                }
                PreparedPatternEngine::GraphAndOverlays { query }
            }
            SubstructureMatchAlgorithm::Incidence => {
                let incidence = pattern.incidence_graph(IncidenceLevel::Constitution);
                let incidence_electrons = precompute_incidence_electrons(&pattern, &incidence);
                let mut query = PreparedSubgraphQuery::new(
                    incidence.raw_graph(),
                    config.subgraph_isomorphism_algorithm,
                );
                if config.reduce_pattern_symmetry {
                    let generators = pattern_symmetry_generators(&pattern, Some(&incidence));
                    if !generators.is_empty() {
                        let symmetry = QuerySymmetryConstraints::try_from_generators(
                            incidence.raw_graph(),
                            &generators,
                        )
                        .expect(
                            "proper pattern incidence symmetry generators are valid graph automorphisms",
                        );
                        query = query.with_symmetry(symmetry);
                    }
                }
                let (node_classes, node_class_count) =
                    build_incidence_node_classes(&incidence, &atom_class_of, atom_classes.len());
                PreparedPatternEngine::Incidence {
                    incidence,
                    query,
                    incidence_electrons,
                    node_classes,
                    node_class_count,
                }
            }
        };

        Self {
            pattern,
            config,
            needs_rings,
            required_elements,
            required_bond_orders,
            present_overlays,
            overlay_counts,
            atom_classes,
            atom_class_of,
            constrained_atom_class_count,
            bond_classes,
            bond_class_of,
            constrained_bond_class_count,
            all_bonds_undetermined,
            engine,
        }
    }

    /// Snapshotted pattern [`Molecule`].
    pub fn pattern(&self) -> &Molecule {
        &self.pattern
    }

    /// Matching configuration bound to this prepared pattern.
    pub fn config(&self) -> SubstructureMatchConfig {
        self.config
    }

    /// Visit each occurrence of `self` within `host` as an injective pattern→host
    /// [`MoleculeCorrespondence`] until traversal completes or `visitor` returns
    /// [`ControlFlow::Break`].
    pub fn visit_matches<B, F>(&self, host: &Molecule, visitor: F) -> ControlFlow<B>
    where
        F: FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    {
        if !self.passes_molecule_prescreen(host) {
            return ControlFlow::Continue(());
        }
        let prepared_host = PreparedSubstructureHost::new(host);
        self.visit_matches_in_host(&prepared_host, visitor)
    }

    /// Visit each occurrence of `self` within `host` as a borrowed slice of matched host
    /// [`AtomId`]s indexed by pattern [`AtomId`], without materializing a
    /// [`MoleculeCorrespondence`], until traversal completes or `visitor` returns
    /// [`ControlFlow::Break`].
    pub fn visit_embeddings<B, F>(&self, host: &Molecule, visitor: F) -> ControlFlow<B>
    where
        F: FnMut(&[AtomId]) -> ControlFlow<B>,
    {
        if !self.passes_molecule_prescreen(host) {
            return ControlFlow::Continue(());
        }
        let prepared_host = PreparedSubstructureHost::new(host);
        self.visit_embeddings_in_host(&prepared_host, visitor)
    }

    /// Collect all occurrences of `self` within `host`.
    pub fn matches(&self, host: &Molecule) -> Vec<MoleculeCorrespondence> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> = self.visit_matches(host, |correspondence| {
            occurrences.push(correspondence);
            ControlFlow::Continue(())
        });
        occurrences
    }

    /// Return whether `self` has at least one occurrence within `host`.
    pub fn has_match(&self, host: &Molecule) -> bool {
        if !self.passes_molecule_prescreen(host) {
            return false;
        }
        let prepared_host = PreparedSubstructureHost::new(host);
        self.has_match_in_host(&prepared_host)
    }

    /// Return the first occurrence of `self` within `host`, or `None` if no match exists.
    pub fn first_match(&self, host: &Molecule) -> Option<MoleculeCorrespondence> {
        match self.visit_matches(host, ControlFlow::Break) {
            ControlFlow::Break(correspondence) => Some(correspondence),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Count all occurrences of `self` within `host` without materializing correspondences.
    pub fn count_matches(&self, host: &Molecule) -> u64 {
        if !self.passes_molecule_prescreen(host) {
            return 0;
        }
        let prepared_host = PreparedSubstructureHost::new(host);
        self.count_matches_in_host(&prepared_host)
    }

    /// Visit each occurrence of `self` within a [`PreparedSubstructureHost`] as an injective
    /// pattern→host [`MoleculeCorrespondence`] until traversal completes or `visitor` returns
    /// [`ControlFlow::Break`].
    pub fn visit_matches_in_host<B, F>(
        &self,
        host: &PreparedSubstructureHost,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    {
        if !self.passes_host_prescreen(host) {
            return ControlFlow::Continue(());
        }
        let pattern_atom_count = self.pattern.atoms().count();
        let host_atom_count = host.molecule.atoms().count();
        with_substructure_scratch(|scratch| {
            let SubstructureScratch { subiso, atoms } = scratch;
            self.run_visit_raw(host, subiso, |embedding| {
                fill_atom_embedding(atoms, embedding, pattern_atom_count);
                let atom_corr = Correspondence::from_images(atoms, host_atom_count);
                match self.pattern.verify_overlays(&host.molecule, atom_corr) {
                    Some(correspondence) => visitor(correspondence),
                    None => ControlFlow::Continue(()),
                }
            })
        })
    }

    /// Visit each occurrence of `self` within a [`PreparedSubstructureHost`] as a borrowed slice
    /// of matched host [`AtomId`]s indexed by pattern [`AtomId`], without materializing a
    /// [`MoleculeCorrespondence`], until traversal completes or `visitor` returns
    /// [`ControlFlow::Break`].
    pub fn visit_embeddings_in_host<B, F>(
        &self,
        host: &PreparedSubstructureHost,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[AtomId]) -> ControlFlow<B>,
    {
        if !self.passes_host_prescreen(host) {
            return ControlFlow::Continue(());
        }
        let pattern_atom_count = self.pattern.atoms().count();
        let host_atom_count = host.molecule.atoms().count();
        let has_overlays = self.present_overlays != 0;
        with_substructure_scratch(|scratch| {
            let SubstructureScratch { subiso, atoms } = scratch;
            self.run_visit_raw(host, subiso, |embedding| {
                fill_atom_embedding(atoms, embedding, pattern_atom_count);
                if has_overlays {
                    let atom_corr = Correspondence::from_images(atoms, host_atom_count);
                    if self
                        .pattern
                        .verify_overlays(&host.molecule, atom_corr)
                        .is_none()
                    {
                        return ControlFlow::Continue(());
                    }
                }
                visitor(atoms)
            })
        })
    }

    /// Collect all occurrences of `self` within a [`PreparedSubstructureHost`].
    pub fn matches_in_host(&self, host: &PreparedSubstructureHost) -> Vec<MoleculeCorrespondence> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> = self.visit_matches_in_host(host, |correspondence| {
            occurrences.push(correspondence);
            ControlFlow::Continue(())
        });
        occurrences
    }

    /// Return whether `self` has at least one occurrence within a [`PreparedSubstructureHost`].
    pub fn has_match_in_host(&self, host: &PreparedSubstructureHost) -> bool {
        if !self.passes_host_prescreen(host) {
            return false;
        }
        if self.present_overlays != 0 {
            return self
                .visit_embeddings_in_host(host, |_| ControlFlow::Break(()))
                .is_break();
        }
        with_substructure_scratch(|scratch| self.run_has_raw(host, &mut scratch.subiso))
    }

    /// Return the first occurrence of `self` within a [`PreparedSubstructureHost`], or `None` if
    /// no match exists.
    pub fn first_match_in_host(
        &self,
        host: &PreparedSubstructureHost,
    ) -> Option<MoleculeCorrespondence> {
        match self.visit_matches_in_host(host, ControlFlow::Break) {
            ControlFlow::Break(correspondence) => Some(correspondence),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Count all occurrences of `self` within a [`PreparedSubstructureHost`] without materializing
    /// correspondences.
    pub fn count_matches_in_host(&self, host: &PreparedSubstructureHost) -> u64 {
        if !self.passes_host_prescreen(host) {
            return 0;
        }
        if self.present_overlays != 0 {
            let mut count = 0u64;
            let _: ControlFlow<()> = self.visit_embeddings_in_host(host, |_| {
                count += 1;
                ControlFlow::Continue(())
            });
            return count;
        }
        with_substructure_scratch(|scratch| self.run_count_raw(host, &mut scratch.subiso))
    }

    fn passes_molecule_prescreen(&self, host: &Molecule) -> bool {
        if self.pattern.atoms().count() > host.atoms().count()
            || self.pattern.bonds().count() > host.bonds().count()
            || self.pattern.raw_graph().max_degree() > host.raw_graph().max_degree()
        {
            return false;
        }
        if self.present_overlays != 0 {
            let host_counts = [
                host.dative_bonds().count() as u32,
                host.aromatic_systems().count() as u32,
                host.multicenter_bonds().count() as u32,
                host.noncovalent_bonds().count() as u32,
                host.stereo_atoms().count() as u32,
                host.stereo_bonds().count() as u32,
            ];
            for (pat_count, host_count) in self.overlay_counts.iter().zip(host_counts) {
                if *pat_count > host_count {
                    return false;
                }
            }
        }
        true
    }

    fn passes_host_prescreen(&self, host: &PreparedSubstructureHost) -> bool {
        if self.pattern.atoms().count() > host.molecule.atoms().count()
            || self.pattern.bonds().count() > host.molecule.bonds().count()
            || self.pattern.raw_graph().max_degree() > host.molecule.raw_graph().max_degree()
            || (self.present_overlays & !host.present_overlays) != 0
        {
            return false;
        }
        if self.present_overlays != 0 {
            for (&pat_count, &host_count) in
                self.overlay_counts.iter().zip(host.overlay_counts.iter())
            {
                if pat_count > host_count {
                    return false;
                }
            }
        }
        if host.all_atoms_lit_element && (self.required_elements & !host.present_elements) != 0 {
            return false;
        }
        if host.all_bonds_fast_lit_order
            && (self.required_bond_orders & !host.present_bond_orders) != 0
        {
            return false;
        }
        true
    }

    fn run_visit_raw<B, F>(
        &self,
        host: &PreparedSubstructureHost,
        scratch: &mut SubgraphIsomorphismScratch,
        visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[NodeId]) -> ControlFlow<B>,
    {
        let host_rings = self
            .needs_rings
            .then(|| host.rings(self.config.relevant_cycle_algorithm));
        match &self.engine {
            PreparedPatternEngine::GraphAndOverlays { query } => {
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.atom_symmetry(self.config.automorphism_algorithm))
                    .flatten();
                let mut predicate = GraphSubstructurePredicate::new(self, host, host_rings);
                query.visit_matches(&host.target, host_sym, &mut predicate, scratch, visitor)
            }
            PreparedPatternEngine::Incidence {
                incidence,
                query,
                incidence_electrons,
                node_classes,
                node_class_count,
            } => {
                let inc_host = host.incidence();
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.incidence_symmetry(self.config.automorphism_algorithm, inc_host))
                    .flatten();
                let mut predicate = IncidenceSubstructurePredicate::new(
                    self,
                    incidence,
                    incidence_electrons,
                    node_classes,
                    *node_class_count,
                    host,
                    inc_host,
                    host_rings,
                );
                query.visit_matches(&inc_host.target, host_sym, &mut predicate, scratch, visitor)
            }
        }
    }

    fn run_has_raw(
        &self,
        host: &PreparedSubstructureHost,
        scratch: &mut SubgraphIsomorphismScratch,
    ) -> bool {
        let host_rings = self
            .needs_rings
            .then(|| host.rings(self.config.relevant_cycle_algorithm));
        match &self.engine {
            PreparedPatternEngine::GraphAndOverlays { query } => {
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.atom_symmetry(self.config.automorphism_algorithm))
                    .flatten();
                let mut predicate = GraphSubstructurePredicate::new(self, host, host_rings);
                query.has_match(&host.target, host_sym, &mut predicate, scratch)
            }
            PreparedPatternEngine::Incidence {
                incidence,
                query,
                incidence_electrons,
                node_classes,
                node_class_count,
            } => {
                let inc_host = host.incidence();
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.incidence_symmetry(self.config.automorphism_algorithm, inc_host))
                    .flatten();
                let mut predicate = IncidenceSubstructurePredicate::new(
                    self,
                    incidence,
                    incidence_electrons,
                    node_classes,
                    *node_class_count,
                    host,
                    inc_host,
                    host_rings,
                );
                query.has_match(&inc_host.target, host_sym, &mut predicate, scratch)
            }
        }
    }

    fn run_count_raw(
        &self,
        host: &PreparedSubstructureHost,
        scratch: &mut SubgraphIsomorphismScratch,
    ) -> u64 {
        let host_rings = self
            .needs_rings
            .then(|| host.rings(self.config.relevant_cycle_algorithm));
        match &self.engine {
            PreparedPatternEngine::GraphAndOverlays { query } => {
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.atom_symmetry(self.config.automorphism_algorithm))
                    .flatten();
                let mut predicate = GraphSubstructurePredicate::new(self, host, host_rings);
                query.count_matches(&host.target, host_sym, &mut predicate, scratch)
            }
            PreparedPatternEngine::Incidence {
                incidence,
                query,
                incidence_electrons,
                node_classes,
                node_class_count,
            } => {
                let inc_host = host.incidence();
                let host_sym = self
                    .config
                    .reduce_host_symmetry
                    .then(|| host.incidence_symmetry(self.config.automorphism_algorithm, inc_host))
                    .flatten();
                let mut predicate = IncidenceSubstructurePredicate::new(
                    self,
                    incidence,
                    incidence_electrons,
                    node_classes,
                    *node_class_count,
                    host,
                    inc_host,
                    host_rings,
                );
                query.count_matches(&inc_host.target, host_sym, &mut predicate, scratch)
            }
        }
    }
}

#[inline]
fn fill_atom_embedding(out: &mut Vec<AtomId>, embedding: &[NodeId], atom_count: usize) {
    out.clear();
    out.extend(embedding[..atom_count].iter().copied().map(AtomId::from));
}

fn pattern_needs_rings(pattern: &Molecule) -> bool {
    pattern
        .atoms()
        .iter()
        .flat_map(|atom| atom.attributes().constraints.iter())
        .any(|constraint| {
            matches!(
                constraint,
                AtomConstraintForm::RingDegree(_)
                    | AtomConstraintForm::RingValence(_)
                    | AtomConstraintForm::RingMembership(_)
            )
        })
        || pattern
            .bonds()
            .iter()
            .flat_map(|bond| bond.attributes().constraints.iter())
            .any(|constraint| matches!(constraint, BondConstraintForm::RingMembership(_)))
}

fn classify_atom_form(form: &AtomForm) -> AtomClassFastKind {
    let AtomForm {
        element,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired_electrons,
        constraints,
    } = form;
    let only_element_and_empty_constraints = isotope_mass.is_undetermined()
        && charge.is_undetermined()
        && implicit_hydrogens.is_undetermined()
        && lone_pairs.is_undetermined()
        && unpaired_electrons.count.is_undetermined()
        && unpaired_electrons.multiplicity.is_undetermined()
        && constraints.is_empty();

    match element {
        ElementForm::Undetermined => {
            if only_element_and_empty_constraints {
                AtomClassFastKind::Any
            } else {
                AtomClassFastKind::Filtered { element_mask: None }
            }
        }
        ElementForm::Lit(el) => {
            let z = el.atomic_number();
            if only_element_and_empty_constraints {
                AtomClassFastKind::LitElementOnly(z)
            } else {
                AtomClassFastKind::Filtered {
                    element_mask: Some(1u128 << z),
                }
            }
        }
        ElementForm::LitSet(elements) => {
            let mut mask = 0u128;
            for el in elements.iter() {
                mask |= 1u128 << el.atomic_number();
            }
            if only_element_and_empty_constraints && !elements.is_empty() {
                AtomClassFastKind::ElementSetOnly(mask)
            } else {
                AtomClassFastKind::Filtered {
                    element_mask: Some(mask),
                }
            }
        }
        ElementForm::NotSet(excluded) => {
            let mut mask = ALL_ELEMENT_BITS;
            for el in excluded.iter() {
                mask &= !(1u128 << el.atomic_number());
            }
            if only_element_and_empty_constraints && excluded.len() < MAX_ATOMIC_NUMBER as usize {
                AtomClassFastKind::ElementSetOnly(mask)
            } else {
                AtomClassFastKind::Filtered {
                    element_mask: Some(mask),
                }
            }
        }
        ElementForm::Var(_) => AtomClassFastKind::Filtered { element_mask: None },
    }
}

fn classify_bond_form(form: &BondForm) -> BondFastKind {
    let BondForm {
        order,
        charge,
        unpaired_electrons,
        constraints,
    } = form;
    if !charge.is_undetermined()
        || !unpaired_electrons.count.is_undetermined()
        || !unpaired_electrons.multiplicity.is_undetermined()
        || !constraints.is_empty()
    {
        return BondFastKind::General;
    }
    match order {
        NumForm::Undetermined => BondFastKind::Undetermined,
        NumForm::Lit(value) if (0..=MAX_FAST_BOND_ORDER).contains(value) => {
            BondFastKind::OrderMask(1u8 << (*value as u8))
        }
        NumForm::LitSet(values)
            if !values.is_empty()
                && values
                    .iter()
                    .all(|value| (0..=MAX_FAST_BOND_ORDER).contains(value)) =>
        {
            let mut mask = 0u8;
            for &value in values.iter() {
                mask |= 1u8 << (value as u8);
            }
            BondFastKind::OrderMask(mask)
        }
        _ => BondFastKind::General,
    }
}

fn build_incidence_node_classes(
    incidence: &IncidenceGraph,
    atom_class_of: &[u32],
    atom_class_count: usize,
) -> (Vec<u32>, usize) {
    let node_count = incidence.node_count();
    let mut node_classes = Vec::with_capacity(node_count);
    node_classes.extend_from_slice(atom_class_of);
    let mut next_class = atom_class_count;

    let overlay_kinds = [
        EntityKind::DativeBond,
        EntityKind::AromaticSystem,
        EntityKind::MulticenterBond,
        EntityKind::NoncovalentBond,
    ];
    for kind in overlay_kinds {
        let count = incidence.node_count_of(kind);
        if count > 0 {
            let class_id = next_class as u32;
            next_class += 1;
            node_classes.resize(node_classes.len() + count, class_id);
        }
    }
    debug_assert_eq!(node_classes.len(), node_count);
    (node_classes, next_class)
}

/// Shared lazy constraint-satisfaction caches and atom/bond matching helpers used by both
/// `GraphAndOverlays` and `Incidence` predicates.
struct EntityMatcher<'a> {
    pattern: &'a PreparedSubstructurePattern,
    host: &'a PreparedSubstructureHost,
    host_rings: Option<&'a RingSet>,
    /// `0 = uncached`, `1 = false`, `2 = true`.
    atom_constraint_cache: Vec<u8>,
    bond_constraint_cache: Vec<u8>,
}

impl<'a> EntityMatcher<'a> {
    fn new(
        pattern: &'a PreparedSubstructurePattern,
        host: &'a PreparedSubstructureHost,
        host_rings: Option<&'a RingSet>,
    ) -> Self {
        Self {
            pattern,
            host,
            host_rings,
            atom_constraint_cache: Vec::new(),
            bond_constraint_cache: Vec::new(),
        }
    }

    #[inline]
    fn match_atom_class(&mut self, class_id: usize, host_atom: AtomId) -> bool {
        let cls = &self.pattern.atom_classes[class_id];
        let pat_form = self.pattern.pattern.atom(cls.representative).attributes();
        let host_form = self.host.molecule.atom(host_atom).attributes();
        if !atom_fields_match(pat_form, host_form) {
            return false;
        }
        let Some(slot) = cls.constraint_slot else {
            return true;
        };
        let host_atom_count = self.host.molecule.atoms().count();
        if self.atom_constraint_cache.is_empty() {
            self.atom_constraint_cache =
                vec![0u8; self.pattern.constrained_atom_class_count * host_atom_count];
        }
        let idx = slot * host_atom_count + host_atom.index();
        match self.atom_constraint_cache[idx] {
            1 => false,
            2 => true,
            _ => {
                let mut reading = self.host.molecule.atom(host_atom).constraints();
                if let Some(rings) = self.host_rings {
                    reading = reading.with_rings(rings);
                }
                let ok = reading.satisfies(&pat_form.constraints);
                self.atom_constraint_cache[idx] = if ok { 2 } else { 1 };
                ok
            }
        }
    }

    #[inline]
    fn match_bond_class(&mut self, class_id: usize, host_bond: BondId) -> bool {
        let cls = &self.pattern.bond_classes[class_id];
        if self.host.all_bonds_fields_satisfiable {
            match cls.fast_kind {
                BondFastKind::Undetermined => return true,
                BondFastKind::OrderMask(mask) if self.host.all_bonds_fast_lit_order => {
                    return (mask & (1u8 << self.host.bond_orders[host_bond.index()])) != 0;
                }
                _ => {}
            }
        }
        let pat_form = self.pattern.pattern.bond(cls.representative).attributes();
        let host_form = self.host.molecule.bond(host_bond).attributes();
        if !bond_fields_match(pat_form, host_form) {
            return false;
        }
        let Some(slot) = cls.constraint_slot else {
            return true;
        };
        let host_bond_count = self.host.molecule.bonds().count();
        if self.bond_constraint_cache.is_empty() {
            self.bond_constraint_cache =
                vec![0u8; self.pattern.constrained_bond_class_count * host_bond_count];
        }
        let idx = slot * host_bond_count + host_bond.index();
        match self.bond_constraint_cache[idx] {
            1 => false,
            2 => true,
            _ => {
                let mut reading = self.host.molecule.bond(host_bond).constraints();
                if let Some(rings) = self.host_rings {
                    reading = reading.with_rings(rings);
                }
                let ok = reading.satisfies(&pat_form.constraints);
                self.bond_constraint_cache[idx] = if ok { 2 } else { 1 };
                ok
            }
        }
    }

    #[inline]
    fn atom_candidates_u64(&mut self, class_id: usize, eligible: u64) -> u64 {
        if eligible == 0 {
            return 0;
        }
        let fast_kind = self.pattern.atom_classes[class_id].fast_kind;
        match fast_kind {
            AtomClassFastKind::Any => {
                if self.host.all_atoms_fields_satisfiable {
                    return eligible;
                }
                self.probe_atom_u64(class_id, eligible)
            }
            AtomClassFastKind::LitElementOnly(z) => {
                let lit_mask = self.host.element_mask_u64(z);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched | self.probe_atom_u64(class_id, eligible & self.host.non_lit_u64())
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_u64());
                    self.probe_atom_u64(class_id, narrowed)
                }
            }
            AtomClassFastKind::ElementSetOnly(allowed) => {
                let lit_mask = self.host.element_set_mask_u64(allowed);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched | self.probe_atom_u64(class_id, eligible & self.host.non_lit_u64())
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_u64());
                    self.probe_atom_u64(class_id, narrowed)
                }
            }
            AtomClassFastKind::Filtered { element_mask } => {
                let narrowed = match element_mask {
                    Some(allowed) => {
                        eligible
                            & (self.host.element_set_mask_u64(allowed) | self.host.non_lit_u64())
                    }
                    None => eligible,
                };
                self.probe_atom_u64(class_id, narrowed)
            }
        }
    }

    #[inline]
    fn atom_candidates_u128(&mut self, class_id: usize, eligible: u128) -> u128 {
        if eligible == 0 {
            return 0;
        }
        let fast_kind = self.pattern.atom_classes[class_id].fast_kind;
        match fast_kind {
            AtomClassFastKind::Any => {
                if self.host.all_atoms_fields_satisfiable {
                    return eligible;
                }
                self.probe_atom_u128(class_id, eligible)
            }
            AtomClassFastKind::LitElementOnly(z) => {
                let lit_mask = self.host.element_mask_u128(z);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched
                            | self.probe_atom_u128(class_id, eligible & self.host.non_lit_u128())
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_u128());
                    self.probe_atom_u128(class_id, narrowed)
                }
            }
            AtomClassFastKind::ElementSetOnly(allowed) => {
                let lit_mask = self.host.element_set_mask_u128(allowed);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched
                            | self.probe_atom_u128(class_id, eligible & self.host.non_lit_u128())
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_u128());
                    self.probe_atom_u128(class_id, narrowed)
                }
            }
            AtomClassFastKind::Filtered { element_mask } => {
                let narrowed = match element_mask {
                    Some(allowed) => {
                        eligible
                            & (self.host.element_set_mask_u128(allowed) | self.host.non_lit_u128())
                    }
                    None => eligible,
                };
                self.probe_atom_u128(class_id, narrowed)
            }
        }
    }

    #[inline]
    fn atom_candidates_word(&mut self, class_id: usize, word_idx: usize, eligible: u64) -> u64 {
        if eligible == 0 {
            return 0;
        }
        let fast_kind = self.pattern.atom_classes[class_id].fast_kind;
        match fast_kind {
            AtomClassFastKind::Any => {
                if self.host.all_atoms_fields_satisfiable {
                    return eligible;
                }
                self.probe_atom_word(class_id, word_idx, eligible)
            }
            AtomClassFastKind::LitElementOnly(z) => {
                let lit_mask = self.host.element_word(z, word_idx);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched
                            | self.probe_atom_word(
                                class_id,
                                word_idx,
                                eligible & self.host.non_lit_word(word_idx),
                            )
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_word(word_idx));
                    self.probe_atom_word(class_id, word_idx, narrowed)
                }
            }
            AtomClassFastKind::ElementSetOnly(allowed) => {
                let lit_mask = self.host.element_set_word(allowed, word_idx);
                if self.host.all_atoms_fields_satisfiable {
                    let matched = eligible & lit_mask;
                    if self.host.all_atoms_lit_element {
                        matched
                    } else {
                        matched
                            | self.probe_atom_word(
                                class_id,
                                word_idx,
                                eligible & self.host.non_lit_word(word_idx),
                            )
                    }
                } else {
                    let narrowed = eligible & (lit_mask | self.host.non_lit_word(word_idx));
                    self.probe_atom_word(class_id, word_idx, narrowed)
                }
            }
            AtomClassFastKind::Filtered { element_mask } => {
                let narrowed = match element_mask {
                    Some(allowed) => {
                        eligible
                            & (self.host.element_set_word(allowed, word_idx)
                                | self.host.non_lit_word(word_idx))
                    }
                    None => eligible,
                };
                self.probe_atom_word(class_id, word_idx, narrowed)
            }
        }
    }

    #[inline]
    fn probe_atom_u64(&mut self, class_id: usize, mut remaining: u64) -> u64 {
        let mut mask = 0u64;
        while remaining != 0 {
            let target = remaining.trailing_zeros();
            remaining &= remaining - 1;
            if self.match_atom_class(class_id, AtomId(target)) {
                mask |= 1u64 << target;
            }
        }
        mask
    }

    #[inline]
    fn probe_atom_u128(&mut self, class_id: usize, mut remaining: u128) -> u128 {
        let mut mask = 0u128;
        while remaining != 0 {
            let target = remaining.trailing_zeros();
            remaining &= remaining - 1;
            if self.match_atom_class(class_id, AtomId(target)) {
                mask |= 1u128 << target;
            }
        }
        mask
    }

    #[inline]
    fn probe_atom_word(&mut self, class_id: usize, word_idx: usize, mut remaining: u64) -> u64 {
        let mut mask = 0u64;
        let base = (word_idx * 64) as u32;
        while remaining != 0 {
            let bit = remaining.trailing_zeros();
            remaining &= remaining - 1;
            if self.match_atom_class(class_id, AtomId(base + bit)) {
                mask |= 1u64 << bit;
            }
        }
        mask
    }
}

struct GraphSubstructurePredicate<'a> {
    matcher: EntityMatcher<'a>,
}

impl<'a> GraphSubstructurePredicate<'a> {
    fn new(
        pattern: &'a PreparedSubstructurePattern,
        host: &'a PreparedSubstructureHost,
        host_rings: Option<&'a RingSet>,
    ) -> Self {
        Self {
            matcher: EntityMatcher::new(pattern, host, host_rings),
        }
    }
}

impl SubgraphIsomorphismPredicate for GraphSubstructurePredicate<'_> {
    #[inline]
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        let class_id = self.matcher.pattern.atom_class_of[query_node.index()] as usize;
        self.matcher
            .match_atom_class(class_id, AtomId::from(target_node))
    }

    #[inline]
    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        let class_id = self.matcher.pattern.bond_class_of[query_edge.index()] as usize;
        self.matcher
            .match_bond_class(class_id, BondId::from(target_edge))
    }

    #[inline]
    fn all_edges_undetermined(&self) -> bool {
        self.matcher.pattern.all_bonds_undetermined
            && self.matcher.host.all_bonds_fields_satisfiable
    }

    #[inline]
    fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
        let class_id = self.matcher.pattern.bond_class_of[query_edge.index()] as usize;
        self.matcher.host.all_bonds_fields_satisfiable
            && matches!(
                self.matcher.pattern.bond_classes[class_id].fast_kind,
                BondFastKind::Undetermined
            )
    }

    #[inline]
    fn query_node_class_count(&self) -> Option<usize> {
        Some(self.matcher.pattern.atom_classes.len())
    }

    #[inline]
    fn query_node_class(&self, query_node: NodeId) -> usize {
        self.matcher.pattern.atom_class_of[query_node.index()] as usize
    }

    #[inline]
    fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
        let class_id = self.matcher.pattern.atom_class_of[query_node.index()] as usize;
        self.matcher.atom_candidates_u64(class_id, eligible)
    }

    #[inline]
    fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
        let class_id = self.matcher.pattern.atom_class_of[query_node.index()] as usize;
        self.matcher.atom_candidates_u128(class_id, eligible)
    }

    fn fill_node_candidates_words(
        &mut self,
        query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let class_id = self.matcher.pattern.atom_class_of[query_node.index()] as usize;
        let len = eligible.len().min(out.len());
        for word_idx in 0..len {
            out[word_idx] =
                self.matcher
                    .atom_candidates_word(class_id, word_idx, eligible[word_idx]);
        }
        out[len..].fill(0);
    }

    #[inline]
    fn edge_neighbors_u64(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u64> {
        let class_id = self.matcher.pattern.bond_class_of[query_edge.index()] as usize;
        match self.matcher.pattern.bond_classes[class_id].fast_kind {
            BondFastKind::OrderMask(orders_mask) => self
                .matcher
                .host
                .bond_order_neighbors_u64(orders_mask, target_node),
            BondFastKind::Undetermined | BondFastKind::General => None,
        }
    }

    #[inline]
    fn edge_neighbors_u128(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u128> {
        let class_id = self.matcher.pattern.bond_class_of[query_edge.index()] as usize;
        match self.matcher.pattern.bond_classes[class_id].fast_kind {
            BondFastKind::OrderMask(orders_mask) => self
                .matcher
                .host
                .bond_order_neighbors_u128(orders_mask, target_node),
            BondFastKind::Undetermined | BondFastKind::General => None,
        }
    }
}

struct IncidenceSubstructurePredicate<'a> {
    matcher: EntityMatcher<'a>,
    pattern_incidence: &'a IncidenceGraph,
    pattern_incidence_electrons: &'a [IncidenceElectron],
    node_classes: &'a [u32],
    node_class_count: usize,
    host_incidence: &'a PreparedIncidenceHost,
}

impl<'a> IncidenceSubstructurePredicate<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        pattern: &'a PreparedSubstructurePattern,
        pattern_incidence: &'a IncidenceGraph,
        pattern_incidence_electrons: &'a [IncidenceElectron],
        node_classes: &'a [u32],
        node_class_count: usize,
        host: &'a PreparedSubstructureHost,
        host_incidence: &'a PreparedIncidenceHost,
        host_rings: Option<&'a RingSet>,
    ) -> Self {
        Self {
            matcher: EntityMatcher::new(pattern, host, host_rings),
            pattern_incidence,
            pattern_incidence_electrons,
            node_classes,
            node_class_count,
            host_incidence,
        }
    }
}

impl SubgraphIsomorphismPredicate for IncidenceSubstructurePredicate<'_> {
    #[inline]
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        match (
            self.pattern_incidence
                .entity(IncidenceNodeId::from(query_node)),
            self.host_incidence
                .incidence
                .entity(IncidenceNodeId::from(target_node)),
        ) {
            (Entity::Atom(pa), Entity::Atom(ha)) => {
                let class_id = self.matcher.pattern.atom_class_of[pa.index()] as usize;
                self.matcher.match_atom_class(class_id, ha)
            }
            (pe, he) => pe.kind() == he.kind(),
        }
    }

    #[inline]
    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        let q_edge = IncidenceEdgeId::from(query_edge);
        let t_edge = IncidenceEdgeId::from(target_edge);
        match (
            self.pattern_incidence.edge(q_edge),
            self.host_incidence.incidence.edge(t_edge),
        ) {
            (IncidenceEdge::Bond(pb), IncidenceEdge::Bond(hb)) => {
                let class_id = self.matcher.pattern.bond_class_of[pb.index()] as usize;
                self.matcher.match_bond_class(class_id, hb)
            }
            (IncidenceEdge::Incidence(pi), IncidenceEdge::Incidence(hi)) => {
                if pi != hi {
                    return false;
                }
                let q_idx = self
                    .pattern_incidence
                    .overlay_incidence_index(q_edge)
                    .expect("overlay edge has incidence index");
                let t_idx = self
                    .host_incidence
                    .incidence
                    .overlay_incidence_index(t_edge)
                    .expect("overlay edge has incidence index");
                incidence_electrons_match(
                    self.pattern_incidence_electrons[q_idx],
                    self.host_incidence.incidence_electrons[t_idx],
                )
            }
            _ => false,
        }
    }

    #[inline]
    fn all_edges_undetermined(&self) -> bool {
        self.pattern_incidence.incidence_count() == 0
            && self.matcher.pattern.all_bonds_undetermined
            && self.matcher.host.all_bonds_fields_satisfiable
    }

    #[inline]
    fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
        match self
            .pattern_incidence
            .edge(IncidenceEdgeId::from(query_edge))
        {
            IncidenceEdge::Bond(pb) => {
                let class_id = self.matcher.pattern.bond_class_of[pb.index()] as usize;
                self.matcher.host.all_bonds_fields_satisfiable
                    && matches!(
                        self.matcher.pattern.bond_classes[class_id].fast_kind,
                        BondFastKind::Undetermined
                    )
            }
            IncidenceEdge::Incidence(_) => false,
        }
    }

    #[inline]
    fn query_node_class_count(&self) -> Option<usize> {
        Some(self.node_class_count)
    }

    #[inline]
    fn query_node_class(&self, query_node: NodeId) -> usize {
        self.node_classes[query_node.index()] as usize
    }

    fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
        match self
            .pattern_incidence
            .entity(IncidenceNodeId::from(query_node))
        {
            Entity::Atom(pa) => {
                let (start, end) = self.host_incidence.kind_range(EntityKind::Atom);
                let atom_eligible = eligible & range_mask_u64(start, end);
                let class_id = self.matcher.pattern.atom_class_of[pa.index()] as usize;
                self.matcher.atom_candidates_u64(class_id, atom_eligible)
            }
            entity => {
                let (start, end) = self.host_incidence.kind_range(entity.kind());
                eligible & range_mask_u64(start, end)
            }
        }
    }

    fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
        match self
            .pattern_incidence
            .entity(IncidenceNodeId::from(query_node))
        {
            Entity::Atom(pa) => {
                let (start, end) = self.host_incidence.kind_range(EntityKind::Atom);
                let atom_eligible = eligible & range_mask_u128(start, end);
                let class_id = self.matcher.pattern.atom_class_of[pa.index()] as usize;
                self.matcher.atom_candidates_u128(class_id, atom_eligible)
            }
            entity => {
                let (start, end) = self.host_incidence.kind_range(entity.kind());
                eligible & range_mask_u128(start, end)
            }
        }
    }

    fn fill_node_candidates_words(
        &mut self,
        query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let len = eligible.len().min(out.len());
        match self
            .pattern_incidence
            .entity(IncidenceNodeId::from(query_node))
        {
            Entity::Atom(pa) => {
                let (start, end) = self.host_incidence.kind_range(EntityKind::Atom);
                let class_id = self.matcher.pattern.atom_class_of[pa.index()] as usize;
                for word_idx in 0..len {
                    let atom_eligible = eligible[word_idx] & range_mask_word(start, end, word_idx);
                    out[word_idx] =
                        self.matcher
                            .atom_candidates_word(class_id, word_idx, atom_eligible);
                }
            }
            entity => {
                let (start, end) = self.host_incidence.kind_range(entity.kind());
                for word_idx in 0..len {
                    out[word_idx] = eligible[word_idx] & range_mask_word(start, end, word_idx);
                }
            }
        }
        out[len..].fill(0);
    }

    #[inline]
    fn edge_neighbors_u64(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u64> {
        let IncidenceEdge::Bond(pb) = self
            .pattern_incidence
            .edge(IncidenceEdgeId::from(query_edge))
        else {
            return None;
        };
        let class_id = self.matcher.pattern.bond_class_of[pb.index()] as usize;
        match self.matcher.pattern.bond_classes[class_id].fast_kind {
            BondFastKind::OrderMask(orders_mask) => self
                .matcher
                .host
                .bond_order_neighbors_u64(orders_mask, target_node),
            BondFastKind::Undetermined | BondFastKind::General => None,
        }
    }

    #[inline]
    fn edge_neighbors_u128(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u128> {
        let IncidenceEdge::Bond(pb) = self
            .pattern_incidence
            .edge(IncidenceEdgeId::from(query_edge))
        else {
            return None;
        };
        let class_id = self.matcher.pattern.bond_class_of[pb.index()] as usize;
        match self.matcher.pattern.bond_classes[class_id].fast_kind {
            BondFastKind::OrderMask(orders_mask) => self
                .matcher
                .host
                .bond_order_neighbors_u128(orders_mask, target_node),
            BondFastKind::Undetermined | BondFastKind::General => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeSet, VecDeque};
    use std::ops::ControlFlow;

    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
    use umol_graph_core::{
        AutomorphismAlgorithm, NodeId, RelevantCycleEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm,
    };

    use super::{PreparedSubstructureHost, PreparedSubstructurePattern};
    use crate::ir::substructure::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
    use crate::ir::substructure::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
    use crate::ir::{
        AromaticSystemForm, AtomForm, AtomId, BondForm, ElectronCountsForm, ElementForm, Molecule,
        MoleculeCorrespondence, MoleculeEntries, NumForm,
    };
    use crate::mol_dsl;

    const SUBISO_ALGS: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];
    const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];
    const AUTO_ALGS: [AutomorphismAlgorithm; 2] =
        [AutomorphismAlgorithm::TrailIr, AutomorphismAlgorithm::Nauty];

    fn group_closure(node_count: usize, generators: &[Vec<NodeId>]) -> Vec<Vec<usize>> {
        let identity: Vec<usize> = (0..node_count).collect();
        let mut seen = BTreeSet::new();
        let mut queue = VecDeque::new();
        seen.insert(identity.clone());
        queue.push_back(identity);
        while let Some(cur) = queue.pop_front() {
            for generator in generators {
                let next: Vec<usize> = cur.iter().map(|&v| generator[v].index()).collect();
                if seen.insert(next.clone()) {
                    queue.push_back(next);
                }
            }
        }
        seen.into_iter().collect()
    }

    fn orbit_set(
        all_embeddings: &[Vec<AtomId>],
        q_group: &[Vec<usize>],
        h_group: &[Vec<usize>],
    ) -> BTreeSet<Vec<Vec<AtomId>>> {
        let mut orbits = BTreeSet::new();
        for embedding in all_embeddings {
            let mut orbit = BTreeSet::new();
            for q_perm in q_group {
                for h_perm in h_group {
                    let mapped: Vec<AtomId> = q_perm
                        .iter()
                        .map(|&q| AtomId(h_perm[embedding[q].index()] as u32))
                        .collect();
                    orbit.insert(mapped);
                }
            }
            orbits.insert(orbit.into_iter().collect());
        }
        orbits
    }

    fn embedding_vectors(matches: &[MoleculeCorrespondence]) -> Vec<Vec<AtomId>> {
        let mut out: Vec<Vec<AtomId>> = matches
            .iter()
            .map(|c| {
                c.atoms()
                    .matched_pairs()
                    .iter()
                    .map(|&(_, host)| host)
                    .collect()
            })
            .collect();
        out.sort();
        out
    }

    #[rstest]
    fn test_prepared_substructure_host_new() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PreparedSubstructureHost>();

        let mut source = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let prepared = PreparedSubstructureHost::new(&source);
        let cloned = prepared.clone();

        // Mutating the source molecule must not affect the prepared host snapshot.
        source.atom_mut(AtomId(2)).attributes_mut().element = ElementForm::Lit(Element::N);
        assert_eq!(prepared.molecule().atoms().count(), 3);
        assert_eq!(
            prepared.molecule().atom(AtomId(2)).attributes().element,
            ElementForm::Lit(Element::O)
        );
        assert_eq!(cloned.molecule(), prepared.molecule());

        // Lazy caches start empty until explicitly needed.
        assert!(!prepared.has_cached_incidence());
        assert!(!prepared.has_cached_rings(RelevantCycleEnumerationAlgorithm::Vismara));
        assert!(!prepared.has_cached_symmetry(AutomorphismAlgorithm::TrailIr));
        assert!(!prepared.has_cached_symmetry(AutomorphismAlgorithm::Nauty));

        // Plain GraphAndOverlays matching populates none of the lazy caches.
        let pat = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            SubstructureMatchConfig {
                match_algorithm: GraphAndOverlays,
                subgraph_isomorphism_algorithm: BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            },
        );
        assert_eq!(pat.count_matches_in_host(&prepared), 1);
        assert!(!prepared.has_cached_incidence());
        assert!(!prepared.has_cached_rings(RelevantCycleEnumerationAlgorithm::Vismara));
        assert!(!prepared.has_cached_symmetry(AutomorphismAlgorithm::TrailIr));

        // Ring constraint populates only the ring cache (and shares it with clones).
        let ring_pat = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C#R0"] :bonds []}"#),
            SubstructureMatchConfig {
                match_algorithm: GraphAndOverlays,
                subgraph_isomorphism_algorithm: BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            },
        );
        assert_eq!(ring_pat.count_matches_in_host(&prepared), 2);
        assert!(prepared.has_cached_rings(RelevantCycleEnumerationAlgorithm::Vismara));
        assert!(cloned.has_cached_rings(RelevantCycleEnumerationAlgorithm::Vismara));
        assert!(!prepared.has_cached_incidence());
        assert!(!prepared.has_cached_symmetry(AutomorphismAlgorithm::TrailIr));

        // Host symmetry populates only the requested automorphism algorithm's cache.
        let sym_pat = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
            SubstructureMatchConfig {
                match_algorithm: GraphAndOverlays,
                subgraph_isomorphism_algorithm: BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: true,
            },
        );
        assert_eq!(sym_pat.count_matches_in_host(&prepared), 2);
        assert!(prepared.has_cached_symmetry(AutomorphismAlgorithm::TrailIr));
        assert!(!prepared.has_cached_symmetry(AutomorphismAlgorithm::Nauty));

        // Incidence matching populates the incidence cache.
        let inc_pat = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            SubstructureMatchConfig {
                match_algorithm: Incidence,
                subgraph_isomorphism_algorithm: BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            },
        );
        assert_eq!(inc_pat.count_matches_in_host(&prepared), 1);
        assert!(prepared.has_cached_incidence());
    }

    #[rstest]
    fn test_prepared_substructure_pattern_new() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PreparedSubstructurePattern>();

        let mut pattern = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let config = SubstructureMatchConfig {
            match_algorithm: GraphAndOverlays,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        let prepared = PreparedSubstructurePattern::new(&pattern, config);
        let cloned = prepared.clone();
        pattern.atom_mut(AtomId(1)).attributes_mut().element = ElementForm::Lit(Element::N);

        let host = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        assert_eq!(prepared.config(), config);
        assert_eq!(cloned.config(), config);
        assert_eq!(
            prepared.pattern().atom(AtomId(1)).attributes().element,
            ElementForm::Lit(Element::O)
        );
        assert!(prepared.has_match(&host));
        assert!(cloned.has_match(&host));
    }

    #[rstest]
    #[case::one_constraint(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#)
    )]
    #[case::two_constraints(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}} {:charge-sum {:sum 0}}]}"#)
    )]
    #[should_panic(expected = "substructure matching does not evaluate molecule constraints")]
    fn test_prepared_substructure_pattern_new_error(
        #[case] pattern: Molecule,
        #[values(GraphAndOverlays, Incidence)] strategy: SubstructureMatchAlgorithm,
    ) {
        let config = SubstructureMatchConfig {
            match_algorithm: strategy,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        let _ = PreparedSubstructurePattern::new(&pattern, config);
    }

    #[rstest]
    #[case::empty_pattern_in_empty_host(
        Molecule::new(),
        Molecule::new(),
        vec![vec![]]
    )]
    #[case::empty_pattern_in_nonempty_host(
        mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
        Molecule::new(),
        vec![vec![]]
    )]
    #[case::concrete_topology(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::element_set_and_not_set(
        mol_dsl!(r#"{:atoms ["C" "N" "O" "S"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#),
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::new(ElementForm::lit_set([Element::N, Element::O])),
                AtomForm::new(ElementForm::not(Element::C)),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::new(NumForm::lit_set([1, 2])))],
            ..Default::default()
        }),
        vec![
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
        ]
    )]
    #[case::undetermined_pattern_and_bonds(
        mol_dsl!(r#"{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        mol_dsl!(r#"{:atoms ["*" "*"] :bonds [[0 1 "*"]]}"#),
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
        ]
    )]
    #[case::prescreen_absent_element(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "S"] :bonds [[0 1 "1"]]}"#),
        vec![]
    )]
    #[case::prescreen_absent_bond_order(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "3"]]}"#),
        vec![]
    )]
    #[case::singleton_litset_host_element_matches_lit_pattern(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::new(ElementForm::lit_set([Element::C])),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::new(NumForm::lit_set([1])))],
            ..Default::default()
        }),
        mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)]]
    )]
    #[case::unsatisfiable_host_atom_rejected_by_wildcard_pattern(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::new(ElementForm::lit_set(Vec::<Element>::new())),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        mol_dsl!(r#"{:atoms ["*"] :bonds []}"#),
        vec![vec![AtomId(1)]]
    )]
    #[case::unsatisfiable_host_bond_rejected_by_undetermined_pattern_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                BondForm::new(NumForm::lit_set(Vec::<i64>::new())),
            )],
            ..Default::default()
        }),
        mol_dsl!(r#"{:atoms ["*" "*"] :bonds [[0 1 "*"]]}"#),
        vec![]
    )]
    #[case::ring_and_atom_constraints(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#R1" "C#R0"] :bonds [[0 1 "1#R0"]]}"#),
        vec![vec![AtomId(2), AtomId(3)]]
    )]
    #[case::aromatic_and_noncovalent_overlays(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "N"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}] :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![AtomId(0), AtomId(2), AtomId(1)],
        ]
    )]
    #[case::stereo_atom_overlay(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        vec![vec![AtomId(0), AtomId(3), AtomId(2), AtomId(1)]]
    )]
    fn test_prepared_substructure_pattern_matches(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let prepared_host = PreparedSubstructureHost::new(&host);
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let prepared = PreparedSubstructurePattern::new(&pattern, config);

                let direct = pattern.substructure_matches(&host, config);
                let via_mol = prepared.matches(&host);
                let via_host = prepared.matches_in_host(&prepared_host);

                assert_eq!(
                    embedding_vectors(&direct),
                    expected,
                    "{strategy:?}/{subiso:?}"
                );
                assert_eq!(
                    embedding_vectors(&via_mol),
                    expected,
                    "{strategy:?}/{subiso:?}"
                );
                assert_eq!(
                    embedding_vectors(&via_host),
                    expected,
                    "{strategy:?}/{subiso:?}"
                );
                assert_eq!(prepared.has_match(&host), !expected.is_empty());
                assert_eq!(
                    prepared.has_match_in_host(&prepared_host),
                    !expected.is_empty()
                );
                assert_eq!(prepared.count_matches(&host), expected.len() as u64);
                assert_eq!(
                    prepared.count_matches_in_host(&prepared_host),
                    expected.len() as u64
                );

                let first_mol = prepared.first_match(&host);
                let first_host = prepared.first_match_in_host(&prepared_host);
                assert_eq!(first_mol.is_some(), !expected.is_empty());
                assert_eq!(first_host, first_mol);
                if let Some(first) = first_host {
                    let atoms: Vec<AtomId> = first
                        .atoms()
                        .matched_pairs()
                        .iter()
                        .map(|&(_, h)| h)
                        .collect();
                    assert!(expected.contains(&atoms));
                }

                let mut visited_embeddings = Vec::new();
                let flow: ControlFlow<()> =
                    prepared.visit_embeddings_in_host(&prepared_host, |emb| {
                        visited_embeddings.push(emb.to_vec());
                        ControlFlow::Continue(())
                    });
                assert_eq!(flow, ControlFlow::Continue(()));
                visited_embeddings.sort();
                assert_eq!(visited_embeddings, expected);

                if !expected.is_empty() {
                    let early =
                        prepared.visit_embeddings(&host, |emb| ControlFlow::Break(emb.to_vec()));
                    let ControlFlow::Break(first_emb) = early else {
                        panic!("expected Break on first embedding");
                    };
                    assert!(expected.contains(&first_emb));
                }
            }
        }
    }

    #[rstest]
    #[case::cyclohexane_c3_chain(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        (12, 6, 1, 1)
    )]
    #[case::tetrahedral_methane_th0_self(
        mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
        (12, 1, 1, 1)
    )]
    #[case::rr_tartaric_half_pattern(
        mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0" "C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [3 4 "1"] [3 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"} {:site 3 :ligands [4 5 0 [:h 3]] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        (2, 2, 1, 1)
    )]
    #[case::meso_tartaric_half_pattern(
        mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0" "C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [3 4 "1"] [3 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"} {:site 3 :ligands [4 5 0 [:h 3]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        (2, 2, 2, 2)
    )]
    #[case::aromatic_nonuniform_triangle_self(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        (2, 1, 1, 1)
    )]
    #[case::two_disjoint_ethanes_in_three_ethanes(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [2 3 "1"] [4 5 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [2 3 "1"]]}"#),
        (24, 3, 1, 1)
    )]
    #[case::malformed_aromatic_electron_length_trivial_fallback(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1, 1])),
            )],
            ..Default::default()
        }),
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1, 1])),
            )],
            ..Default::default()
        }),
        (0, 0, 0, 0)
    )]
    fn test_prepared_substructure_pattern_matches_symmetry(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected_counts: (u64, u64, u64, u64),
    ) {
        let prepared_host = PreparedSubstructureHost::new(&host);
        let q_gens = pattern.proper_atom_generators(AutomorphismAlgorithm::TrailIr, true);
        let h_gens = host.proper_atom_generators(AutomorphismAlgorithm::TrailIr, true);
        let q_group = group_closure(pattern.atoms().count(), &q_gens);
        let h_group = group_closure(host.atoms().count(), &h_gens);
        let trivial_q = vec![(0..pattern.atoms().count()).collect::<Vec<_>>()];
        let trivial_h = vec![(0..host.atoms().count()).collect::<Vec<_>>()];

        let all_embeddings = embedding_vectors(&pattern.substructure_matches(
            &host,
            SubstructureMatchConfig {
                match_algorithm: GraphAndOverlays,
                subgraph_isomorphism_algorithm: BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            },
        ));

        let modes = [
            (false, false, expected_counts.0, &trivial_q, &trivial_h),
            (true, false, expected_counts.1, &q_group, &trivial_h),
            (false, true, expected_counts.2, &trivial_q, &h_group),
            (true, true, expected_counts.3, &q_group, &h_group),
        ];

        for (red_q, red_h, expected_count, active_q, active_h) in modes {
            let expected_orbits = orbit_set(&all_embeddings, active_q, active_h);
            assert_eq!(expected_orbits.len() as u64, expected_count);

            for strategy in STRATEGIES {
                let mut baseline_reps: Option<Vec<Vec<AtomId>>> = None;
                for auto_alg in AUTO_ALGS {
                    for subiso in SUBISO_ALGS {
                        let config = SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm: auto_alg,
                            reduce_pattern_symmetry: red_q,
                            reduce_host_symmetry: red_h,
                        };
                        let prepared = PreparedSubstructurePattern::new(&pattern, config);
                        let reps = embedding_vectors(&prepared.matches_in_host(&prepared_host));
                        let count = prepared.count_matches_in_host(&prepared_host);

                        assert_eq!(
                            count, expected_count,
                            "{strategy:?}/{auto_alg:?}/{subiso:?} ({red_q},{red_h})"
                        );
                        assert_eq!(
                            reps.len() as u64,
                            expected_count,
                            "{strategy:?}/{auto_alg:?}/{subiso:?} ({red_q},{red_h})"
                        );
                        assert_eq!(
                            prepared.has_match_in_host(&prepared_host),
                            expected_count > 0
                        );

                        // Every emitted representative must land in a distinct expected orbit,
                        // covering all expected orbits exactly once.
                        let covered = orbit_set(&reps, active_q, active_h);
                        assert_eq!(
                            covered, expected_orbits,
                            "{strategy:?}/{auto_alg:?}/{subiso:?} ({red_q},{red_h})"
                        );

                        // Within each strategy, the canonical representative set must be identical
                        // across all 4 SubgraphIsomorphismAlgorithms and both AutomorphismAlgorithms.
                        if let Some(expected_reps) = &baseline_reps {
                            assert_eq!(
                                &reps, expected_reps,
                                "{strategy:?}/{auto_alg:?}/{subiso:?} ({red_q},{red_h})"
                            );
                        } else {
                            baseline_reps = Some(reps);
                        }
                    }
                }
            }
        }
    }

    #[rstest]
    #[case::u128_chain_70(70)]
    #[case::wide_chain_140(140)]
    fn test_prepared_substructure_pattern_matches_wide(#[case] atom_count: usize) {
        let mut atoms = vec![AtomForm::from_element(Element::C); atom_count];
        atoms[atom_count / 2] = AtomForm::from_element(Element::N);
        let bonds: Vec<(AtomId, AtomId, BondForm)> = (0..atom_count - 1)
            .map(|i| {
                let order = if i % 2 == 0 { 1 } else { 2 };
                (
                    AtomId(i as u32),
                    AtomId((i + 1) as u32),
                    BondForm::from_order(order),
                )
            })
            .collect();
        let host = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        });
        let pattern = mol_dsl!(r#"{:atoms ["C" "N" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#);
        let prepared_host = PreparedSubstructureHost::new(&host);

        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let prepared = PreparedSubstructurePattern::new(&pattern, config);
                assert_eq!(prepared.count_matches_in_host(&prepared_host), 1);
                let matches = prepared.matches_in_host(&prepared_host);
                assert_eq!(matches.len(), 1);
                let center = AtomId((atom_count / 2) as u32);
                assert_eq!(matches[0].atoms().right_of(AtomId(1)), Some(center));
            }
        }
    }

    #[rstest]
    fn test_prepared_substructure_pattern_visit_embeddings_reentrant() {
        let host = PreparedSubstructureHost::new(&mol_dsl!(
            r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ));
        let config = SubstructureMatchConfig {
            match_algorithm: GraphAndOverlays,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        let outer = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
            config,
        );
        let inner = PreparedSubstructurePattern::new(
            &mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            config,
        );

        let mut nested_counts = Vec::new();
        let _: ControlFlow<()> = outer.visit_embeddings_in_host(&host, |_| {
            nested_counts.push(inner.count_matches_in_host(&host));
            ControlFlow::Continue(())
        });
        assert_eq!(nested_counts, vec![1, 1]);
    }

    #[rstest]
    fn test_prepared_substructure_pattern_incidence_u64_atoms_u128_nodes() {
        let mut atoms = vec![AtomForm::from_element(Element::C); 63];
        atoms[0] = AtomForm::from_element(Element::N);
        let bonds: Vec<(AtomId, AtomId, BondForm)> = (0..62)
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId((i + 1) as u32),
                    BondForm::from_order(if i == 0 { 2 } else { 1 }),
                )
            })
            .collect();
        let aromatic = vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(10), AtomId(11)],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            ),
            (
                vec![AtomId(20), AtomId(21)],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            ),
        ];
        let host = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            aromatic,
            ..Default::default()
        });
        let pattern = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            ..Default::default()
        });
        let prepared_host = PreparedSubstructureHost::new(&host);

        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let prepared = PreparedSubstructurePattern::new(&pattern, config);
                assert_eq!(prepared.count_matches_in_host(&prepared_host), 1);
                let matches = prepared.matches_in_host(&prepared_host);
                assert_eq!(matches.len(), 1);
                assert_eq!(matches[0].atoms().right_of(AtomId(0)), Some(AtomId(0)));
                assert_eq!(matches[0].atoms().right_of(AtomId(1)), Some(AtomId(1)));
            }
        }
    }
}
