//! Aggregate canonicalization inputs and failures.

use std::cell::OnceCell;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};

use thiserror::Error;
use umol_graph_core::{AutomorphismAlgorithm, EdgeId, GraphRemapping, NodeId, Remapping};

use super::delta::Delta;
use super::entity::DescriptionLevel;
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
    StereoBondId,
};
use super::molecule::Molecule;
use super::reaction::Reaction;
use super::reaction_span::ReactionSpan;
use super::remap::MoleculeRemapping;
use super::symmetry::{
    canonical_search, canonical_search_from_partition, canonical_search_with_generator_filter,
    rank_colors, remapping_exists, remapping_preserves_incidences, AutomorphismAdapter,
    AutomorphismColor, CanonicalCandidate, CanonicalSearchResult, InitialColors, OrderedPartition,
};
#[cfg(test)]
use super::symmetry::{same_cells, AutomorphismAdapterOutput, CanonicalSearchStats};
use super::traits::Reframe;

mod molecule;
mod reaction_span;

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeContext {
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

fn delta_canonicalize_level(delta: &Delta) -> DescriptionLevel {
    match delta {
        Delta::Atom(_) | Delta::Bond(_) => DescriptionLevel::Topology,
        Delta::DativeBond(_)
        | Delta::AromaticSystem(_)
        | Delta::MulticenterBond(_)
        | Delta::NoncovalentBond(_) => DescriptionLevel::Constitution,
        Delta::StereoAtom(_) | Delta::StereoBond(_) => DescriptionLevel::Structure,
        Delta::Constraint(_) => DescriptionLevel::Full,
    }
}

fn molecule_canonicalize_level(molecule: &Molecule) -> DescriptionLevel {
    if molecule.has_constraints() {
        DescriptionLevel::Full
    } else if molecule.has_stereo_atoms() || molecule.has_stereo_bonds() {
        DescriptionLevel::Structure
    } else if molecule.has_dative_bonds()
        || molecule.has_aromatic_systems()
        || molecule.has_multicenter_bonds()
        || molecule.has_noncovalent_bonds()
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
}

fn reaction_canonicalize_level(reaction: &Reaction) -> DescriptionLevel {
    reaction.deltas().iter().fold(
        molecule_canonicalize_level(reaction.lhs()),
        |level, delta| level.max(delta_canonicalize_level(delta)),
    )
}

fn reaction_span_canonicalize_level(span: &ReactionSpan) -> DescriptionLevel {
    if !span.constraints().is_empty() {
        DescriptionLevel::Full
    } else if span.stereo_atoms().count() > 0 || span.stereo_bonds().count() > 0 {
        DescriptionLevel::Structure
    } else if span.dative_bonds().count() > 0
        || span.aromatic_systems().count() > 0
        || span.multicenter_bonds().count() > 0
        || span.noncovalent_bonds().count() > 0
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
}

/// Canonical entity-frame selection for complete indexed graph-IR aggregates.
///
/// Canonicalization extends [`Reframe`]: normalization is the fixed-frame prefix, reframing adds
/// participant-frame selection, and canonicalization additionally selects entity ids. It preserves
/// the represented aggregate, transports every reference and position-sensitive value, and
/// normalizes carried forms in the selected frame. Construction of the canonical form is fallible;
/// canonical equality is total and compares two values that both fail by framed equality.
///
/// # Semantic properties
///
/// For every integrity-valid aggregate and fixed context:
///
/// - successful complete canonicalization is exactly idempotent and invariant under valid dense
///   entity remapping;
/// - `canonical_eq` is reflexive, symmetric, and transitive, including on values whose
///   canonicalization fails;
/// - successful canonical hashes are invariant under valid dense entity remapping, and canonical
///   equality implies equal canonical hashes when both hash operations succeed.
///
/// For a fixed umol release and context, canonicalization is deterministic. During the 0.x series,
/// the order of keys and the resulting representatives may change between releases. Returned
/// canonical forms are ordinary IR values without version provenance and must not be used as
/// persistent cross-release identifiers.
pub trait Canonicalize: Reframe {
    type Error;

    /// Construct the complete canonical form.
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific canonicalization error when intrinsic normalization is
    /// unsatisfiable or a carried attribute cannot be transported in its stored frame.
    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error>;

    /// Construct the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping maps every entity id in the input frame to its id in the returned
    /// canonical frame. Each of its eight entity-kind components is total on both sides and
    /// therefore represents a dense bijection. For [`Reaction`], the two frames are the complete
    /// union frames of the materialized input and returned reaction spans.
    ///
    /// # Errors
    ///
    /// Returns the same aggregate-specific canonicalization error as
    /// [`Self::canonicalize`]. A [`Reaction`] also reports failure to materialize its reaction span.
    ///
    /// # Semantic properties
    ///
    /// Discarding the remapping yields exactly [`Self::canonicalize`] under the same context.
    /// Remapping the input through the remapping and then applying [`Reframe::reframe`]
    /// reconstructs the returned canonical value.
    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error>;

    /// Hash the complete canonical form.
    ///
    /// The returned value uses Rust's [`DefaultHasher`]. It may collide and is not a persistent
    /// identifier: both the hasher and umol's canonical representatives may change between
    /// releases.
    ///
    /// # Errors
    ///
    /// Returns the same error as [`Self::canonicalize`].
    ///
    /// # Semantic properties
    ///
    /// Successful hashes are invariant under valid dense entity remapping. Canonically equal
    /// values have equal canonical hashes when both hash operations succeed.
    fn canonical_hash(&self, context: &CanonicalizeContext) -> Result<u64, Self::Error>;

    /// Compare complete canonical forms.
    ///
    /// Structural identity short-circuits the operation. A value whose canonicalization fails is
    /// unequal to one whose canonicalization succeeds, and two values whose canonicalizations
    /// both fail are compared by [`Reframe::framed_eq`].
    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool;
}

fn hash_value<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

/// A leaf key that is built when it is first compared or read.
struct DeferredKey<'a, K, Build> {
    order: Vec<NodeId>,
    build: &'a Build,
    key: OnceCell<K>,
}

impl<'a, K, Build: Fn(&[NodeId]) -> K> DeferredKey<'a, K, Build> {
    fn new(order: &[NodeId], build: &'a Build) -> Self {
        Self {
            order: order.to_vec(),
            build,
            key: OnceCell::new(),
        }
    }

    fn get(&self) -> &K {
        self.key.get_or_init(|| (self.build)(&self.order))
    }

    /// The key, if a comparison built it.
    fn into_built(self) -> Option<K> {
        self.key.into_inner()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialEq for DeferredKey<'_, K, Build> {
    fn eq(&self, other: &Self) -> bool {
        self.get() == other.get()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Eq for DeferredKey<'_, K, Build> {}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialOrd for DeferredKey<'_, K, Build> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Ord for DeferredKey<'_, K, Build> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.get().cmp(other.get())
    }
}

/// Topology search input with atoms and non-modal bond subdivisions as vertices.
///
/// Bonds in the largest normalized bond-form class are represented directly as edges. Every other
/// bond remains a colored subdivision vertex. The complete typed leaf key and remapping still
/// order and map every bond.
#[derive(Clone, Debug)]
struct TopologySearchInput {
    adapter: AutomorphismAdapter,
    initial_partition: OrderedPartition,
}

/// One position per color integer in each of two lists, in the order of the integers: the first
/// position that holds the integer. `None` when the lists do not hold each integer equally often.
fn color_representatives(source: &[u32], target: &[u32]) -> Option<Vec<(usize, usize)>> {
    if source.len() != target.len() {
        return None;
    }
    let (Some(&first), Some(&last)) = (source.iter().min(), source.iter().max()) else {
        return Some(Vec::new());
    };
    let integer_count = (last - first) as usize + 1;
    let mut counts = vec![0_i64; integer_count];
    let mut source_positions = vec![usize::MAX; integer_count];
    let mut target_positions = vec![usize::MAX; integer_count];
    for (position, &color) in source.iter().enumerate() {
        let index = (color - first) as usize;
        counts[index] += 1;
        if source_positions[index] == usize::MAX {
            source_positions[index] = position;
        }
    }
    for (position, &color) in target.iter().enumerate() {
        if color < first || color > last {
            return None;
        }
        let index = (color - first) as usize;
        counts[index] -= 1;
        if target_positions[index] == usize::MAX {
            target_positions[index] = position;
        }
    }
    if counts.iter().any(|&count| count != 0) {
        return None;
    }
    Some(
        source_positions
            .into_iter()
            .zip(target_positions)
            .filter(|&(source_position, _)| source_position != usize::MAX)
            .collect(),
    )
}

fn molecule_remapping(images: &[Vec<usize>; 8]) -> MoleculeRemapping {
    fn remapping<Id: Copy + Debug + Into<usize> + From<usize>>(images: &[usize]) -> Remapping<Id> {
        Remapping::new(images.iter().copied().map(Id::from).collect())
            .expect("selected entity order is a permutation")
    }
    MoleculeRemapping::new(
        GraphRemapping::new(
            remapping::<NodeId>(&images[0]),
            remapping::<EdgeId>(&images[1]),
        ),
        remapping::<DativeBondId>(&images[2]),
        remapping::<AromaticSystemId>(&images[3]),
        remapping::<MulticenterBondId>(&images[4]),
        remapping::<NoncovalentBondId>(&images[5]),
        remapping::<StereoAtomId>(&images[6]),
        remapping::<StereoBondId>(&images[7]),
    )
}

/// Failure to construct a canonical [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizeError {
    /// Normalization or span materialization failed, or literal electron counts
    /// do not match their atom list.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

impl From<ReactionSpanCanonicalizeError> for ReactionCanonicalizeError {
    fn from(error: ReactionSpanCanonicalizeError) -> Self {
        match error {
            ReactionSpanCanonicalizeError::Contradiction(contradiction) => {
                Self::Contradiction(contradiction)
            }
        }
    }
}

#[cfg(test)]
mod tests;
