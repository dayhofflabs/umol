//! Constraint reading of an entity, reached by accessor chaining:
//! `molecule.atom(id).constraints()`, `molecule.bond(id).constraints()`.
//!
//! The bare read API is inherited from the stored container and keeps its
//! meanings — typed getters, `iter`, `is_empty` all read the asserted side,
//! which is what constraints store. The keyed core names both sides
//! explicitly: `asserted(key)` reads storage, `derived(key)` reads present
//! relations (vacuous on absence), and `derived_complete(key)` adds the
//! closure — under the caller's claim that the relation set is complete,
//! absence of a resolution-written overlay closes to its definite negative.
//! Topology keys read identically under both derived accessors; only the
//! overlay keys have an absence cell for the closure to fill. Mutation never
//! routes through this view; it belongs to the stored container.

use super::super::boolean::BooleanForm;
use super::super::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemConstraintsForm,
    AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
    BondConstraintForm, BondConstraintKey, BondConstraintsForm, DativeBondConstraintForm,
    DativeBondConstraintKey, DativeBondConstraintsForm, FluxionalityForm, LigandPermutation,
    LigandSymmetryForm, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
    MulticenterBondConstraintsForm, MulticenterValenceForm, NoncovalentBondConstraintForm,
    NoncovalentBondConstraintKey, NoncovalentBondConstraintsForm, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoAtomConstraintsForm, StereoBondConstraintForm,
    StereoBondConstraintKey, StereoBondConstraintsForm, StereoLigandPair, StereogenicityForm,
    TopicityForm, TopicityRelationForm,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::ring::RingSet;
use super::super::stereo::{CisTransStereoForm, TetrahedralStereoForm};
use super::super::traits::Lattice;
use super::aromatic::{aromatic_system_asserted_constraints, aromatic_system_derived_constraint};
use super::atom::{
    atom_asserted_complete_constraint, atom_asserted_constraints, atom_derived_constraint,
};
use super::bond::{
    bond_asserted_complete_constraint, bond_asserted_constraints, bond_derived_constraint,
};
use super::dative::{
    dative_bond_asserted_complete_constraint, dative_bond_asserted_constraints,
    dative_bond_derived_constraint,
};
use super::multicenter::{
    multicenter_bond_asserted_constraints, multicenter_bond_derived_constraint,
};
use super::noncovalent::{
    noncovalent_bond_asserted_constraints, noncovalent_bond_derived_constraint,
};
use super::stereo::{
    stereo_atom_asserted_constraints, stereo_atom_derived_constraint,
    stereo_bond_asserted_constraints, stereo_bond_derived_constraint,
};

/// Constraint reading of one atom: the asserted side under the container's
/// read API, and both sides under the keyed accessors.
#[derive(Clone, Copy, Debug)]
pub struct AtomConstraintsView<'a> {
    molecule: &'a Molecule,
    atom: AtomId,
    rings: Option<&'a RingSet>,
}

impl<'a> AtomConstraintsView<'a> {
    pub(crate) fn new(molecule: &'a Molecule, atom: AtomId) -> Self {
        Self {
            molecule,
            atom,
            rings: None,
        }
    }

    /// Attach ring context for the ring keys (`#R`, `#x`, `#y`).
    pub fn with_rings(mut self, rings: &'a RingSet) -> Self {
        self.rings = Some(rings);
        self
    }

    /// The stored side of `key`; absence is the vacuous constraint.
    pub fn asserted(&self, key: AtomConstraintKey) -> Option<&'a AtomConstraintForm> {
        atom_asserted_constraints(self.molecule, self.atom).get(key)
    }

    /// The asserted side of `key` under resolution's closed-world claim:
    /// absence of a constraint is actual absence, so the absence cell of an
    /// entity-creating overlay key closes to its definite negative. Aromatic
    /// evidence merges both dialect placements: the atom's own assertion,
    /// else an incident bond asserting `#a`, else `NotAromatic`. Topology
    /// keys return the assertion or `None` — they have no absence cell.
    /// Never reads relations; matching keeps the open-world [`Self::asserted`].
    pub fn asserted_complete(&self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        atom_asserted_complete_constraint(self.molecule, self.atom, key)
    }

    /// The derived side of `key`, obtained by projection from present
    /// relations only; vacuous on absence.
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived(&self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        atom_derived_constraint(self.molecule, self.atom, self.rings, key, false)
    }

    /// The derived side of `key` under the closure: absence of a
    /// resolution-written overlay yields its definite negative. Positive
    /// incidence and the topology keys agree with [`Self::derived`].
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived_complete(&self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        atom_derived_constraint(self.molecule, self.atom, self.rings, key, true)
    }

    /// Whether this atom's constraint reading satisfies `pattern`: every
    /// pattern entry is refined by the meet of the asserted and
    /// [`Self::derived_complete`] sides at its key — the query-against-host
    /// reading. An internally conflicted key (the sides meet to `⊥`)
    /// satisfies nothing. Evaluation is driven by the pattern's keys; an
    /// empty pattern is satisfied.
    ///
    /// # Panics
    ///
    /// A ring key in `pattern` without ring context ([`Self::with_rings`]) is
    /// a caller error.
    pub fn satisfies(&self, pattern: &AtomConstraintsForm) -> bool {
        pattern.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived_complete(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => entry.as_undetermined(),
            };
            host.satisfies(entry)
        })
    }

    /// Whether `other` is compatible with this atom's constraint reading:
    /// for every key of `other`, a meet with the asserted and
    /// [`Self::derived`] sides exists — the narrowing-admissibility reading.
    /// A key on which this atom carries nothing constrains nothing; an
    /// internally conflicted key (the sides meet to `⊥`) is compatible with
    /// nothing. Evaluation is driven by `other`'s keys; an empty `other` is
    /// compatible.
    ///
    /// # Panics
    ///
    /// A ring key in `other` without ring context ([`Self::with_rings`]) is a
    /// caller error.
    pub fn is_compatible(&self, other: &AtomConstraintsForm) -> bool {
        other.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => return true,
            };
            entry.is_compatible(&host)
        })
    }

    // The stored container's read API, inherited with its meanings intact:
    // every accessor below reads the asserted side.

    pub fn valence(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).valence()
    }

    pub fn aromatic_valence(&self) -> Option<&'a AromaticValenceForm> {
        atom_asserted_constraints(self.molecule, self.atom).aromatic_valence()
    }

    pub fn multicenter_valence(&self) -> Option<&'a MulticenterValenceForm> {
        atom_asserted_constraints(self.molecule, self.atom).multicenter_valence()
    }

    pub fn tetrahedral_stereo(&self) -> Option<&'a TetrahedralStereoForm> {
        atom_asserted_constraints(self.molecule, self.atom).tetrahedral_stereo()
    }

    pub fn degree(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).degree()
    }

    pub fn total_degree(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).total_degree()
    }

    pub fn total_valence(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).total_valence()
    }

    pub fn ring_degree(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).ring_degree()
    }

    pub fn ring_valence(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).ring_valence()
    }

    pub fn total_hydrogens(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).total_hydrogens()
    }

    pub fn donated_pairs(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).donated_pairs()
    }

    pub fn accepted_pairs(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).accepted_pairs()
    }

    pub fn ring_count(&self) -> Option<&'a NumForm> {
        atom_asserted_constraints(self.molecule, self.atom).ring_count()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = &'a AtomConstraintForm> {
        atom_asserted_constraints(self.molecule, self.atom).iter()
    }

    pub fn is_empty(&self) -> bool {
        atom_asserted_constraints(self.molecule, self.atom).is_empty()
    }

    pub fn len(&self) -> usize {
        atom_asserted_constraints(self.molecule, self.atom).len()
    }
}

/// Constraint reading of one bond: the asserted side under the container's
/// read API, and both sides under the keyed accessors.
#[derive(Clone, Copy, Debug)]
pub struct BondConstraintsView<'a> {
    molecule: &'a Molecule,
    bond: BondId,
    rings: Option<&'a RingSet>,
}

impl<'a> BondConstraintsView<'a> {
    pub(crate) fn new(molecule: &'a Molecule, bond: BondId) -> Self {
        Self {
            molecule,
            bond,
            rings: None,
        }
    }

    /// Attach ring context for the ring key (`#R`).
    pub fn with_rings(mut self, rings: &'a RingSet) -> Self {
        self.rings = Some(rings);
        self
    }

    /// The stored side of `key`; absence is the vacuous constraint.
    pub fn asserted(&self, key: BondConstraintKey) -> Option<&'a BondConstraintForm> {
        bond_asserted_constraints(self.molecule, self.bond).get(key)
    }

    /// The asserted side of `key` under resolution's closed-world claim:
    /// absence of a constraint is actual absence. Never reads relations — a
    /// bond inside a stored aromatic system without its own `#a` assertion
    /// still reads `Lit(false)` here; matching keeps the open-world
    /// [`Self::asserted`].
    pub fn asserted_complete(&self, key: BondConstraintKey) -> Option<BondConstraintForm> {
        bond_asserted_complete_constraint(self.molecule, self.bond, key)
    }

    /// The derived side of `key`, obtained by projection from present
    /// relations only; vacuous on absence.
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived(&self, key: BondConstraintKey) -> Option<BondConstraintForm> {
        bond_derived_constraint(self.molecule, self.bond, self.rings, key, false)
    }

    /// The derived side of `key` under the closure: absence of a
    /// resolution-written overlay yields its definite negative. Positive
    /// incidence and the topology keys agree with [`Self::derived`].
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived_complete(&self, key: BondConstraintKey) -> Option<BondConstraintForm> {
        bond_derived_constraint(self.molecule, self.bond, self.rings, key, true)
    }

    /// Whether this bond's constraint reading satisfies `pattern`: every
    /// pattern entry is refined by the meet of the asserted and
    /// [`Self::derived_complete`] sides at its key — the query-against-host
    /// reading. An internally conflicted key (the sides meet to `⊥`)
    /// satisfies nothing. Evaluation is driven by the pattern's keys; an
    /// empty pattern is satisfied.
    ///
    /// # Panics
    ///
    /// A ring key in `pattern` without ring context ([`Self::with_rings`]) is
    /// a caller error.
    pub fn satisfies(&self, pattern: &BondConstraintsForm) -> bool {
        pattern.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived_complete(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => entry.as_undetermined(),
            };
            host.satisfies(entry)
        })
    }

    /// Whether `other` is compatible with this bond's constraint reading:
    /// for every key of `other`, a meet with the asserted and
    /// [`Self::derived`] sides exists — the narrowing-admissibility reading.
    /// A key on which this bond carries nothing constrains nothing; an
    /// internally conflicted key (the sides meet to `⊥`) is compatible with
    /// nothing. Evaluation is driven by `other`'s keys; an empty `other` is
    /// compatible.
    ///
    /// # Panics
    ///
    /// A ring key in `other` without ring context ([`Self::with_rings`]) is a
    /// caller error.
    pub fn is_compatible(&self, other: &BondConstraintsForm) -> bool {
        other.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => return true,
            };
            entry.is_compatible(&host)
        })
    }

    // The stored container's read API, inherited with its meanings intact:
    // every accessor below reads the asserted side.

    pub fn aromatic(&self) -> BooleanForm {
        bond_asserted_constraints(self.molecule, self.bond).aromatic()
    }

    pub fn cis_trans_stereo(&self) -> Option<&'a CisTransStereoForm> {
        bond_asserted_constraints(self.molecule, self.bond).cis_trans_stereo()
    }

    pub fn ring_count(&self) -> Option<&'a NumForm> {
        bond_asserted_constraints(self.molecule, self.bond).ring_count()
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&'a NumForm> {
        bond_asserted_constraints(self.molecule, self.bond).ring_size_count(s)
    }

    pub fn iter(&self) -> impl Iterator<Item = &'a BondConstraintForm> {
        bond_asserted_constraints(self.molecule, self.bond).iter()
    }

    pub fn is_empty(&self) -> bool {
        bond_asserted_constraints(self.molecule, self.bond).is_empty()
    }

    pub fn len(&self) -> usize {
        bond_asserted_constraints(self.molecule, self.bond).len()
    }
}

/// Generates the chained constraint view of a family without ring context:
/// the keyed core and the comparisons; the container's positional read API
/// (`iter`/`is_empty`/`len`) is included, and family typed getters follow in
/// per-family impl blocks.
macro_rules! constraints_view {
    ($view:ident, $entity:literal, $id:ty, $form:ty, $key:ty, $container:ty,
     $asserted:path, $derived:path) => {
        #[doc = concat!("Constraint reading of one ", $entity, ": the asserted side under")]
        #[doc = "the container's read API, and both sides under the keyed accessors."]
        #[derive(Clone, Copy, Debug)]
        pub struct $view<'a> {
            molecule: &'a Molecule,
            id: $id,
        }

        impl<'a> $view<'a> {
            pub(crate) fn new(molecule: &'a Molecule, id: $id) -> Self {
                Self { molecule, id }
            }

            /// The stored side of `key`; absence is the vacuous constraint.
            pub fn asserted(&self, key: $key) -> Option<&'a $form> {
                $asserted(self.molecule, self.id).get(key)
            }

            /// The derived side of `key`, obtained by projection from present
            /// relations only; vacuous on absence and wherever the family
            /// defines no projection.
            pub fn derived(&self, key: $key) -> Option<$form> {
                $derived(self.molecule, self.id, key, false)
            }

            /// The derived side of `key` under the closure; agrees with
            /// [`Self::derived`] wherever no absence cell exists.
            pub fn derived_complete(&self, key: $key) -> Option<$form> {
                $derived(self.molecule, self.id, key, true)
            }

            /// Whether this entity's constraint reading satisfies `pattern`:
            /// every pattern entry is refined by the meet of the asserted and
            /// [`Self::derived_complete`] sides at its key. An internally
            /// conflicted key (the sides meet to `⊥`) satisfies nothing; an
            /// empty pattern is satisfied.
            pub fn satisfies(&self, pattern: &$container) -> bool {
                pattern.iter().all(|entry| {
                    let key = entry.key();
                    let host = match (self.asserted(key), self.derived_complete(key)) {
                        (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                            Some(host) => host,
                            None => return false,
                        },
                        (Some(asserted), None) => asserted.clone(),
                        (None, Some(derived)) => derived,
                        (None, None) => entry.as_undetermined(),
                    };
                    host.satisfies(entry)
                })
            }

            /// Whether `other` is compatible with this entity's constraint
            /// reading: for every key of `other`, a meet with the asserted
            /// and [`Self::derived`] sides exists. A key on which this entity
            /// carries nothing constrains nothing; an internally conflicted
            /// key is compatible with nothing; an empty `other` is
            /// compatible.
            pub fn is_compatible(&self, other: &$container) -> bool {
                other.iter().all(|entry| {
                    let key = entry.key();
                    let host = match (self.asserted(key), self.derived(key)) {
                        (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                            Some(host) => host,
                            None => return false,
                        },
                        (Some(asserted), None) => asserted.clone(),
                        (None, Some(derived)) => derived,
                        (None, None) => return true,
                    };
                    entry.is_compatible(&host)
                })
            }

            pub fn iter(&self) -> std::slice::Iter<'a, $form> {
                $asserted(self.molecule, self.id).iter()
            }

            pub fn is_empty(&self) -> bool {
                $asserted(self.molecule, self.id).is_empty()
            }

            pub fn len(&self) -> usize {
                $asserted(self.molecule, self.id).len()
            }
        }
    };
}

constraints_view!(
    DativeBondConstraintsView,
    "dative bond",
    DativeBondId,
    DativeBondConstraintForm,
    DativeBondConstraintKey,
    DativeBondConstraintsForm,
    dative_bond_asserted_constraints,
    dative_bond_derived_constraint
);

constraints_view!(
    AromaticSystemConstraintsView,
    "aromatic system",
    AromaticSystemId,
    AromaticSystemConstraintForm,
    AromaticSystemConstraintKey,
    AromaticSystemConstraintsForm,
    aromatic_system_asserted_constraints,
    aromatic_system_derived_constraint
);

constraints_view!(
    MulticenterBondConstraintsView,
    "multicenter bond",
    MulticenterBondId,
    MulticenterBondConstraintForm,
    MulticenterBondConstraintKey,
    MulticenterBondConstraintsForm,
    multicenter_bond_asserted_constraints,
    multicenter_bond_derived_constraint
);

constraints_view!(
    NoncovalentBondConstraintsView,
    "noncovalent bond",
    NoncovalentBondId,
    NoncovalentBondConstraintForm,
    NoncovalentBondConstraintKey,
    NoncovalentBondConstraintsForm,
    noncovalent_bond_asserted_constraints,
    noncovalent_bond_derived_constraint
);

constraints_view!(
    StereoAtomConstraintsView,
    "stereo atom",
    StereoAtomId,
    StereoAtomConstraintForm,
    StereoAtomConstraintKey,
    StereoAtomConstraintsForm,
    stereo_atom_asserted_constraints,
    stereo_atom_derived_constraint
);

constraints_view!(
    StereoBondConstraintsView,
    "stereo bond",
    StereoBondId,
    StereoBondConstraintForm,
    StereoBondConstraintKey,
    StereoBondConstraintsForm,
    stereo_bond_asserted_constraints,
    stereo_bond_derived_constraint
);

// The stored containers' typed read API, inherited per family with its
// meanings intact: every getter reads the asserted side.

impl<'a> DativeBondConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim:
    /// absence of a constraint is actual absence. Never reads relations.
    pub fn asserted_complete(
        &self,
        key: DativeBondConstraintKey,
    ) -> Option<DativeBondConstraintForm> {
        dative_bond_asserted_complete_constraint(self.molecule, self.id, key)
    }

    pub fn aromatic(&self) -> BooleanForm {
        dative_bond_asserted_constraints(self.molecule, self.id).aromatic()
    }

    pub fn ring_count(&self) -> Option<&'a NumForm> {
        dative_bond_asserted_constraints(self.molecule, self.id).ring_count()
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&'a NumForm> {
        dative_bond_asserted_constraints(self.molecule, self.id).ring_size_count(s)
    }
}

impl<'a> AromaticSystemConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim; the
    /// family has no absence cell, so the reading equals the assertion.
    pub fn asserted_complete(
        &self,
        key: AromaticSystemConstraintKey,
    ) -> Option<AromaticSystemConstraintForm> {
        self.asserted(key).cloned()
    }

    pub fn electron_count(&self) -> NumForm {
        aromatic_system_asserted_constraints(self.molecule, self.id).electron_count()
    }
}

impl<'a> MulticenterBondConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim; the
    /// family has no absence cell, so the reading equals the assertion.
    pub fn asserted_complete(
        &self,
        key: MulticenterBondConstraintKey,
    ) -> Option<MulticenterBondConstraintForm> {
        self.asserted(key).cloned()
    }

    pub fn electron_count(&self) -> NumForm {
        multicenter_bond_asserted_constraints(self.molecule, self.id).electron_count()
    }
}

impl<'a> NoncovalentBondConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim; the
    /// family has no absence cell, so the reading equals the assertion.
    pub fn asserted_complete(
        &self,
        key: NoncovalentBondConstraintKey,
    ) -> Option<NoncovalentBondConstraintForm> {
        self.asserted(key).cloned()
    }

    pub fn intramolecular(&self) -> BooleanForm {
        noncovalent_bond_asserted_constraints(self.molecule, self.id).intramolecular()
    }
}

impl<'a> StereoAtomConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim; the
    /// family has no absence cell, so the reading equals the assertion.
    pub fn asserted_complete(
        &self,
        key: StereoAtomConstraintKey,
    ) -> Option<StereoAtomConstraintForm> {
        self.asserted(key).cloned()
    }

    pub fn ligand_symmetries(&self) -> impl Iterator<Item = &'a LigandSymmetryForm> {
        stereo_atom_asserted_constraints(self.molecule, self.id).ligand_symmetries()
    }

    pub fn fluxionalities(&self) -> impl Iterator<Item = &'a FluxionalityForm> {
        stereo_atom_asserted_constraints(self.molecule, self.id).fluxionalities()
    }

    pub fn fluxionality(&self, permutation: LigandPermutation) -> FluxionalityForm {
        stereo_atom_asserted_constraints(self.molecule, self.id).fluxionality(permutation)
    }

    pub fn topicities(&self) -> impl Iterator<Item = &'a TopicityForm> {
        stereo_atom_asserted_constraints(self.molecule, self.id).topicities()
    }

    pub fn topicity(&self, pair: StereoLigandPair) -> TopicityRelationForm {
        stereo_atom_asserted_constraints(self.molecule, self.id).topicity(pair)
    }

    pub fn stereogenicity(&self) -> StereogenicityForm {
        stereo_atom_asserted_constraints(self.molecule, self.id).stereogenicity()
    }
}

impl<'a> StereoBondConstraintsView<'a> {
    /// The asserted side of `key` under resolution's closed-world claim; the
    /// family has no absence cell, so the reading equals the assertion.
    pub fn asserted_complete(
        &self,
        key: StereoBondConstraintKey,
    ) -> Option<StereoBondConstraintForm> {
        self.asserted(key).cloned()
    }

    pub fn ligand_symmetries(&self) -> impl Iterator<Item = &'a LigandSymmetryForm> {
        stereo_bond_asserted_constraints(self.molecule, self.id).ligand_symmetries()
    }

    pub fn fluxionalities(&self) -> impl Iterator<Item = &'a FluxionalityForm> {
        stereo_bond_asserted_constraints(self.molecule, self.id).fluxionalities()
    }

    pub fn fluxionality(&self, permutation: LigandPermutation) -> FluxionalityForm {
        stereo_bond_asserted_constraints(self.molecule, self.id).fluxionality(permutation)
    }

    pub fn topicities(&self) -> impl Iterator<Item = &'a TopicityForm> {
        stereo_bond_asserted_constraints(self.molecule, self.id).topicities()
    }

    pub fn topicity(&self, pair: StereoLigandPair) -> TopicityRelationForm {
        stereo_bond_asserted_constraints(self.molecule, self.id).topicity(pair)
    }

    pub fn stereogenicity(&self) -> StereogenicityForm {
        stereo_bond_asserted_constraints(self.molecule, self.id).stereogenicity()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::{
        AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemConstraintsForm,
        AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
        BondConstraintForm, BondConstraintKey, BondConstraintsForm, DativeBondConstraintForm,
        DativeBondConstraintKey, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
        MulticenterBondConstraintsForm, MulticenterValenceForm, NoncovalentBondConstraintForm,
        NoncovalentBondConstraintKey, NoncovalentBondConstraintsForm, RingScope,
        StereoAtomConstraintForm, StereoAtomConstraintKey, StereoBondConstraintForm,
        StereoBondConstraintKey, StereogenicityForm,
    };
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;
    use crate::ir::ring::{RingConfig, RingModel};
    use crate::ir::stereo::{
        CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
        StereoKind, Stereogenicity, TetrahedralStereoForm,
    };
    use crate::mol_dsl;

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(true))),
    )]
    #[case::aromatic_absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(false))),
    )]
    #[case::aromatic_in_system_unasserted(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            )],
            ..Default::default()
        }),
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(false))),
    )]
    #[case::cis_trans_absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        BondConstraintKey::CisTransStereo,
        Some(BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo)),
    )]
    #[case::ring_absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintKey::RingMembership(RingScope::All),
        None,
    )]
    fn test_bond_constraints_view_asserted_complete(
        #[case] molecule: Molecule,
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        assert_eq!(
            molecule.bond(BondId(0)).constraints().asserted_complete(key),
            expected
        );
    }

    #[rstest]
    #[case::cis_trans_coset_plural(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([0, 1])),
            )],
            ..Default::default()
        }),
        BondId(1),
        BondConstraintKey::CisTransStereo,
        Some(BondConstraintForm::cis_trans_stereo(
            CisTransStereoForm::stereo(StereoCoset::lit_set([0, 1])),
        )),
    )]
    fn test_bond_constraints_view_derived_complete_determination(
        #[case] molecule: Molecule,
        #[case] bond: BondId,
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        assert_eq!(
            molecule.bond(bond).constraints().derived_complete(key),
            expected
        );
    }

    #[rstest]
    #[case::present(
        Some(AtomConstraintForm::valence(4)),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(4))
    )]
    #[case::absent(None, AtomConstraintKey::Degree, None)]
    fn test_atom_constraints_view_asserted(
        #[case] set: Option<AtomConstraintForm>,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let mut atom = AtomForm::from_element(Element::C);
        if let Some(constraint) = set {
            atom.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().asserted(key),
            expected.as_ref()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        mol_dsl!(r#"{:atoms ["C#v4"] :bonds []}"#),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(4)),
    )]
    #[case::aromatic_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic)),
    )]
    #[case::aromatic_bond_adjacent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::aromatic(NumForm::Undetermined),
        )),
    )]
    #[case::aromatic_in_system_unasserted(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::from_electrons(vec![1]))],
            ..Default::default()
        }),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic)),
    )]
    #[case::stereo_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo)),
    )]
    #[case::multicenter_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::MulticenterValence,
        Some(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter)),
    )]
    #[case::donated_pairs_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(NumForm::Lit(0))),
    )]
    #[case::accepted_pairs_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::AcceptedPairs,
        Some(AtomConstraintForm::accepted_pairs(NumForm::Lit(0))),
    )]
    #[case::topology_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::Valence,
        None,
    )]
    fn test_atom_constraints_view_asserted_complete(
        #[case] molecule: Molecule,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().asserted_complete(key),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(2)),
    )]
    #[case::degree(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintKey::Degree,
        Some(AtomConstraintForm::degree(1)),
    )]
    #[case::total_degree(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalDegree,
        Some(AtomConstraintForm::total_degree(4)),
    )]
    #[case::total_hydrogens(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalHydrogens,
        Some(AtomConstraintForm::total_hydrogens(4)),
    )]
    #[case::total_valence(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalValence,
        Some(AtomConstraintForm::total_valence(4)),
    )]
    #[case::donated_pairs_absent(
        mol_dsl!(r#"{:atoms ["N"] :bonds []}"#),
        AtomConstraintKey::DonatedPairs,
        None,
    )]
    #[case::donated_pairs_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        }),
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(1)),
    )]
    #[case::aromatic_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::AromaticValence,
        None,
    )]
    #[case::aromatic_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::from_electrons(vec![1]))],
            ..Default::default()
        }),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::aromatic(NumForm::Lit(1)),
        )),
    )]
    #[case::aromatic_bond_marked(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        AtomConstraintKey::AromaticValence,
        None,
    )]
    #[case::multicenter_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::MulticenterValence,
        None,
    )]
    #[case::tetrahedral_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::TetrahedralStereo,
        None,
    )]
    #[case::tetrahedral_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        }),
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::stereo(StereoCoset::Lit(1)),
        )),
    )]
    #[case::tetrahedral_undetermined_configuration(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm {
                    configuration: StereoConfigurationForm::Undetermined,
                    ..StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                },
            )],
            ..Default::default()
        }),
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::Undetermined,
        )),
    )]
    fn test_atom_constraints_view_derived(
        #[case] molecule: Molecule,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().derived(key),
            expected
        );
    }

    #[rstest]
    #[case::valence_agrees(AtomConstraintKey::Valence, Some(AtomConstraintForm::valence(0)))]
    #[case::donated_pairs_closed(
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(0))
    )]
    #[case::accepted_pairs_closed(
        AtomConstraintKey::AcceptedPairs,
        Some(AtomConstraintForm::accepted_pairs(0))
    )]
    #[case::aromatic_closed(
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic))
    )]
    #[case::multicenter_closed(
        AtomConstraintKey::MulticenterValence,
        Some(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter,))
    )]
    #[case::tetrahedral_closed(
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo))
    )]
    fn test_atom_constraints_view_derived_complete(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().derived_complete(key),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_bond_order_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
            ..Default::default()
        }),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(NumForm::Undetermined)),
    )]
    #[case::degree_bond_order_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
            ..Default::default()
        }),
        AtomConstraintKey::Degree,
        Some(AtomConstraintForm::degree(1)),
    )]
    #[case::total_degree_hydrogens_open(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::TotalDegree,
        Some(AtomConstraintForm::total_degree(NumForm::Undetermined)),
    )]
    #[case::total_hydrogens_explicit_neighbor(
        mol_dsl!(r#"{:atoms ["C#h3" "H#h0"] :bonds [[0 1 "1"]]}"#),
        AtomConstraintKey::TotalHydrogens,
        Some(AtomConstraintForm::total_hydrogens(4)),
    )]
    #[case::total_hydrogens_neighbor_element_open(
        mol_dsl!(r#"{:atoms ["C#h0" "*#h0"] :bonds [[0 1 "1"]]}"#),
        AtomConstraintKey::TotalHydrogens,
        Some(AtomConstraintForm::total_hydrogens(NumForm::Undetermined)),
    )]
    #[case::total_valence_bond_order_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                {
                    let mut atom = AtomForm::from_element(Element::C);
                    atom.implicit_hydrogens = NumForm::Lit(3);
                    atom
                },
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
            ..Default::default()
        }),
        AtomConstraintKey::TotalValence,
        Some(AtomConstraintForm::total_valence(NumForm::Undetermined)),
    )]
    #[case::donated_pairs_order_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::default())],
            ..Default::default()
        }),
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(NumForm::Undetermined)),
    )]
    #[case::aromatic_electrons_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::default())],
            ..Default::default()
        }),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::aromatic(NumForm::Undetermined),
        )),
    )]
    #[case::multicenter_electrons_open(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            ..Default::default()
        }),
        AtomConstraintKey::MulticenterValence,
        Some(AtomConstraintForm::multicenter_valence(
            MulticenterValenceForm::multicenter(NumForm::Undetermined),
        )),
    )]
    #[case::tetrahedral_coset_plural(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1])),
            )],
            ..Default::default()
        }),
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::stereo(StereoCoset::lit_set([0, 1])),
        )),
    )]
    fn test_atom_constraints_view_derived_complete_determination(
        #[case] molecule: Molecule,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().derived_complete(key),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_derived_complete_bond_marked() {
        // Derived reads relations only: a bond-carried aromatic mark is an
        // assertion, merged on the `asserted_complete` side; the closure
        // reads the absent relation as its definite negative.
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#);
        assert_eq!(
            molecule
                .atom(AtomId(0))
                .constraints()
                .derived_complete(AtomConstraintKey::AromaticValence),
            Some(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::NotAromatic,
            )),
        );
    }

    #[rstest]
    #[case::ring_degree(
        AtomConstraintKey::RingDegree,
        Some(AtomConstraintForm::ring_degree(2))
    )]
    #[case::ring_valence(
        AtomConstraintKey::RingValence,
        Some(AtomConstraintForm::ring_valence(2))
    )]
    #[case::membership_all(
        AtomConstraintKey::RingMembership(RingScope::All),
        Some(AtomConstraintForm::ring_membership(RingScope::All, 1))
    )]
    #[case::membership_size_match(
        AtomConstraintKey::RingMembership(RingScope::Size(6)),
        Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1))
    )]
    #[case::membership_size_no_match(
        AtomConstraintKey::RingMembership(RingScope::Size(5)),
        Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 0))
    )]
    fn test_atom_constraints_view_derived_ring(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let rings = molecule
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        assert_eq!(
            molecule
                .atom(AtomId(0))
                .constraints()
                .with_rings(&rings)
                .derived(key),
            expected
        );
    }

    #[rstest]
    #[should_panic(expected = "ring constraint key requires ring context")]
    fn test_atom_constraints_view_derived_error() {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        molecule
            .atom(AtomId(0))
            .constraints()
            .derived(AtomConstraintKey::RingDegree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::new(),
        true,
    )]
    #[case::valence_match(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        true,
    )]
    #[case::valence_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]),
        false,
    )]
    #[case::vacuous_entry(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(NumForm::Undetermined)]),
        true,
    )]
    #[case::closure_not_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        true,
    )]
    #[case::closure_aromatic_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::Lit(1))),
        ]),
        false,
    )]
    #[case::conflicted_host_mismatching_pattern(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        false,
    )]
    #[case::conflicted_host_matching_pattern(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        false,
    )]
    fn test_atom_constraints_view_satisfies(
        #[case] molecule: Molecule,
        #[case] pattern: AtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().satisfies(&pattern),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_satisfies_ring() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let rings = molecule
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        let pattern = AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(
            RingScope::Size(6),
            1,
        )]);
        assert!(molecule
            .atom(AtomId(0))
            .constraints()
            .with_rings(&rings)
            .satisfies(&pattern));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::new(),
        true,
    )]
    #[case::shared_key_compatible(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        true,
    )]
    #[case::shared_key_conflict(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]),
        false,
    )]
    #[case::absent_overlay_skipped(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        true,
    )]
    #[case::absent_stereo_skipped(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
        ]),
        true,
    )]
    #[case::asserted_conflict(
        mol_dsl!(r#"{:atoms ["C#a1"] :bonds []}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        false,
    )]
    #[case::conflicted_host(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        false,
    )]
    fn test_atom_constraints_view_is_compatible(
        #[case] molecule: Molecule,
        #[case] other: AtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().is_compatible(&other),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_container_methods() {
        // The container read API is inherited with its meanings intact: every
        // getter reads the asserted side.
        let mut atom = AtomForm::from_element(Element::C);
        for constraint in [
            AtomConstraintForm::valence(4),
            AtomConstraintForm::donated_pairs(1),
            AtomConstraintForm::accepted_pairs(2),
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::Lit(1))),
            AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
            AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
            AtomConstraintForm::degree(3),
            AtomConstraintForm::total_degree(4),
            AtomConstraintForm::total_valence(4),
            AtomConstraintForm::ring_degree(2),
            AtomConstraintForm::ring_valence(2),
            AtomConstraintForm::total_hydrogens(1),
            AtomConstraintForm::ring_membership(RingScope::All, 1),
        ] {
            atom.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        let view = molecule.atom(AtomId(0)).constraints();
        assert_eq!(view.valence(), Some(&NumForm::Lit(4)));
        assert_eq!(view.donated_pairs(), Some(&NumForm::Lit(1)));
        assert_eq!(view.accepted_pairs(), Some(&NumForm::Lit(2)));
        assert_eq!(
            view.aromatic_valence(),
            Some(&AromaticValenceForm::Aromatic(NumForm::Lit(1)))
        );
        assert_eq!(
            view.multicenter_valence(),
            Some(&MulticenterValenceForm::NotMulticenter)
        );
        assert_eq!(
            view.tetrahedral_stereo(),
            Some(&TetrahedralStereoForm::NotStereo)
        );
        assert_eq!(view.degree(), Some(&NumForm::Lit(3)));
        assert_eq!(view.total_degree(), Some(&NumForm::Lit(4)));
        assert_eq!(view.total_valence(), Some(&NumForm::Lit(4)));
        assert_eq!(view.ring_degree(), Some(&NumForm::Lit(2)));
        assert_eq!(view.ring_valence(), Some(&NumForm::Lit(2)));
        assert_eq!(view.total_hydrogens(), Some(&NumForm::Lit(1)));
        assert_eq!(view.ring_count(), Some(&NumForm::Lit(1)));
        assert_eq!(view.len(), 13);
        assert!(!view.is_empty());
        assert_eq!(view.iter().count(), 13);
    }

    #[rstest]
    #[case::present(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(true)))
    )]
    #[case::absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintKey::Aromatic,
        None
    )]
    fn test_bond_constraints_view_asserted(
        #[case] molecule: Molecule,
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        assert_eq!(
            molecule.bond(BondId(0)).constraints().asserted(key),
            expected.as_ref()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_in_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(0), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 1, 1]),
            )],
            ..Default::default()
        }),
        BondId(0),
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(true))),
    )]
    #[case::aromatic_absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondId(0),
        BondConstraintKey::Aromatic,
        None,
    )]
    #[case::cis_trans_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        }),
        BondId(1),
        BondConstraintKey::CisTransStereo,
        Some(BondConstraintForm::cis_trans_stereo(
            CisTransStereoForm::stereo(StereoCoset::Lit(1)),
        )),
    )]
    #[case::cis_trans_undetermined_configuration(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(2)),
                (AtomId(3), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(4), BondForm::from_order(1)),
                (AtomId(3), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm {
                    configuration: StereoConfigurationForm::Undetermined,
                    ..StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))
                },
            )],
            ..Default::default()
        }),
        BondId(1),
        BondConstraintKey::CisTransStereo,
        Some(BondConstraintForm::cis_trans_stereo(
            CisTransStereoForm::Undetermined,
        )),
    )]
    #[case::cis_trans_absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondId(0),
        BondConstraintKey::CisTransStereo,
        None,
    )]
    fn test_bond_constraints_view_derived(
        #[case] molecule: Molecule,
        #[case] bond: BondId,
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        assert_eq!(molecule.bond(bond).constraints().derived(key), expected);
    }

    #[rstest]
    #[case::aromatic_closed(
        BondConstraintKey::Aromatic,
        Some(BondConstraintForm::aromatic(BooleanForm::Lit(false)))
    )]
    #[case::cis_trans_closed(
        BondConstraintKey::CisTransStereo,
        Some(BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo))
    )]
    fn test_bond_constraints_view_derived_complete(
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        assert_eq!(
            molecule.bond(BondId(0)).constraints().derived_complete(key),
            expected
        );
    }

    #[rstest]
    #[case::membership_all(
        BondConstraintKey::RingMembership(RingScope::All),
        Some(BondConstraintForm::ring_membership(RingScope::All, 1))
    )]
    #[case::membership_size_match(
        BondConstraintKey::RingMembership(RingScope::Size(6)),
        Some(BondConstraintForm::ring_membership(RingScope::Size(6), 1))
    )]
    #[case::membership_size_no_match(
        BondConstraintKey::RingMembership(RingScope::Size(5)),
        Some(BondConstraintForm::ring_membership(RingScope::Size(5), 0))
    )]
    fn test_bond_constraints_view_derived_ring(
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let rings = molecule
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        assert_eq!(
            molecule
                .bond(BondId(0))
                .constraints()
                .with_rings(&rings)
                .derived(key),
            expected
        );
    }

    #[rstest]
    #[should_panic(expected = "ring constraint key requires ring context")]
    fn test_bond_constraints_view_derived_error() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        molecule
            .bond(BondId(0))
            .constraints()
            .derived(BondConstraintKey::RingMembership(RingScope::All));
    }

    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintsForm::new(),
        true
    )]
    #[case::closure_not_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(false))]),
        true
    )]
    #[case::closure_aromatic_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(true))]),
        false
    )]
    #[case::conflicted_host(
        // Asserted `#a` on a bond outside any system conflicts with the
        // closure's negative: the sides meet to ⊥, so nothing is satisfied.
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(true))]),
        false
    )]
    #[case::in_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(0), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 1, 1]),
            )],
            ..Default::default()
        }),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(true))]),
        true
    )]
    fn test_bond_constraints_view_satisfies(
        #[case] molecule: Molecule,
        #[case] pattern: BondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.bond(BondId(0)).constraints().satisfies(&pattern),
            expected
        );
    }

    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintsForm::new(),
        true
    )]
    #[case::absent_overlay_skipped(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(false))]),
        true
    )]
    #[case::asserted_conflict(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        BondConstraintsForm::from_iter([BondConstraintForm::aromatic(BooleanForm::Lit(false))]),
        false
    )]
    fn test_bond_constraints_view_is_compatible(
        #[case] molecule: Molecule,
        #[case] other: BondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.bond(BondId(0)).constraints().is_compatible(&other),
            expected
        );
    }

    #[rstest]
    fn test_bond_constraints_view_container_methods() {
        let mut bond = BondForm::from_order(1);
        for constraint in [
            BondConstraintForm::aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo),
            BondConstraintForm::ring_membership(RingScope::All, 1),
        ] {
            bond.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), bond)],
            ..Default::default()
        });
        let view = molecule.bond(BondId(0)).constraints();
        assert_eq!(view.aromatic(), BooleanForm::Lit(true));
        assert_eq!(
            view.cis_trans_stereo(),
            Some(&CisTransStereoForm::NotStereo)
        );
        assert_eq!(view.ring_count(), Some(&NumForm::Lit(1)));
        assert_eq!(view.ring_size_count(6), None);
        assert_eq!(view.len(), 3);
        assert!(!view.is_empty());
        assert_eq!(view.iter().count(), 3);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::binary_shared_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            )],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        }),
        false,
        Some(DativeBondConstraintForm::aromatic(true)),
    )]
    #[case::binary_unshared(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        }),
        false,
        None,
    )]
    #[case::binary_unshared_complete(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        }),
        true,
        Some(DativeBondConstraintForm::aromatic(false)),
    )]
    #[case::multi_donor_complete(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            dative: vec![(
                vec![AtomId(0), AtomId(1)],
                AtomId(2),
                DativeBondForm::from_order(1),
            )],
            ..Default::default()
        }),
        true,
        None,
    )]
    fn test_dative_bond_constraints_view_derived(
        #[case] molecule: Molecule,
        #[case] complete: bool,
        #[case] expected: Option<DativeBondConstraintForm>,
    ) {
        let view = molecule.dative_bond(DativeBondId(0)).constraints();
        let derived = if complete {
            view.derived_complete(DativeBondConstraintKey::Aromatic)
        } else {
            view.derived(DativeBondConstraintKey::Aromatic)
        };
        assert_eq!(derived, expected);
    }

    #[rstest]
    #[case::present(
        Some(DativeBondConstraintForm::aromatic(true)),
        DativeBondConstraintKey::Aromatic,
        Some(DativeBondConstraintForm::aromatic(true))
    )]
    #[case::aromatic_absent(
        None,
        DativeBondConstraintKey::Aromatic,
        Some(DativeBondConstraintForm::aromatic(false))
    )]
    #[case::ring_absent(None, DativeBondConstraintKey::RingMembership(RingScope::All), None)]
    fn test_dative_bond_constraints_view_asserted_complete(
        #[case] set: Option<DativeBondConstraintForm>,
        #[case] key: DativeBondConstraintKey,
        #[case] expected: Option<DativeBondConstraintForm>,
    ) {
        let mut form = DativeBondForm::from_order(1);
        if let Some(constraint) = set {
            form.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), form)],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .constraints()
                .asserted_complete(key),
            expected
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_view_derived_ring() {
        // The ring key has no projection for dative bonds.
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .constraints()
                .derived_complete(DativeBondConstraintKey::RingMembership(RingScope::All)),
            None
        );
    }

    #[rstest]
    #[case::count_match(
        AromaticSystemConstraintsForm::from_iter([
            AromaticSystemConstraintForm::electron_count(6),
        ]),
        true
    )]
    #[case::count_mismatch(
        AromaticSystemConstraintsForm::from_iter([
            AromaticSystemConstraintForm::electron_count(5),
        ]),
        false
    )]
    fn test_aromatic_system_constraints_view_satisfies(
        #[case] pattern: AromaticSystemConstraintsForm,
        #[case] expected: bool,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 2, 2]),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .constraints()
                .satisfies(&pattern),
            expected
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_derived() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 2, 2]),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .constraints()
                .derived(AromaticSystemConstraintKey::ElectronCount),
            Some(AromaticSystemConstraintForm::electron_count(6))
        );
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsForm::from_iter([
            AromaticSystemConstraintForm::electron_count(6),
        ]),
        Some(AromaticSystemConstraintForm::electron_count(6))
    )]
    #[case::absent(AromaticSystemConstraintsForm::default(), None)]
    fn test_aromatic_system_constraints_view_asserted_complete(
        #[case] constraints: AromaticSystemConstraintsForm,
        #[case] expected: Option<AromaticSystemConstraintForm>,
    ) {
        let mut form = AromaticSystemForm::from_electrons(vec![2, 2, 2]);
        form.constraints = constraints;
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], form)],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .constraints()
                .asserted_complete(AromaticSystemConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_derived_complete_determination() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .constraints()
                .derived_complete(AromaticSystemConstraintKey::ElectronCount),
            Some(AromaticSystemConstraintForm::electron_count(
                NumForm::Undetermined,
            ))
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_derived() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::B); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 0, 1]),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .constraints()
                .derived(MulticenterBondConstraintKey::ElectronCount),
            Some(MulticenterBondConstraintForm::electron_count(2))
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsForm::from_iter([
            MulticenterBondConstraintForm::electron_count(2),
        ]),
        Some(MulticenterBondConstraintForm::electron_count(2))
    )]
    #[case::absent(MulticenterBondConstraintsForm::default(), None)]
    fn test_multicenter_bond_constraints_view_asserted_complete(
        #[case] constraints: MulticenterBondConstraintsForm,
        #[case] expected: Option<MulticenterBondConstraintForm>,
    ) {
        let mut form = MulticenterBondForm::from_electrons(vec![1, 1, 0]);
        form.constraints = constraints;
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::B); 3],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], form)],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .constraints()
                .asserted_complete(MulticenterBondConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_derived_complete_determination() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::B); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .constraints()
                .derived_complete(MulticenterBondConstraintKey::ElectronCount),
            Some(MulticenterBondConstraintForm::electron_count(
                NumForm::Undetermined,
            ))
        );
    }

    #[rstest]
    #[case::same_component(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        }),
        Some(NoncovalentBondConstraintForm::intramolecular(true)),
    )]
    #[case::cross_component(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            noncovalent: vec![(
                AtomId(0),
                AtomId(2),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        }),
        Some(NoncovalentBondConstraintForm::intramolecular(false)),
    )]
    fn test_noncovalent_bond_constraints_view_derived(
        #[case] molecule: Molecule,
        #[case] expected: Option<NoncovalentBondConstraintForm>,
    ) {
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondId(0))
                .constraints()
                .derived(NoncovalentBondConstraintKey::Intramolecular),
            expected
        );
    }

    #[rstest]
    #[case::empty(NoncovalentBondConstraintsForm::new(), true)]
    #[case::conflict(
        NoncovalentBondConstraintsForm::from_iter([
            NoncovalentBondConstraintForm::intramolecular(false),
        ]),
        false
    )]
    fn test_noncovalent_bond_constraints_view_is_compatible(
        #[case] other: NoncovalentBondConstraintsForm,
        #[case] expected: bool,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondId(0))
                .constraints()
                .is_compatible(&other),
            expected
        );
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsForm::from_iter([
            NoncovalentBondConstraintForm::intramolecular(true),
        ]),
        Some(NoncovalentBondConstraintForm::intramolecular(true))
    )]
    #[case::absent(NoncovalentBondConstraintsForm::default(), None)]
    fn test_noncovalent_bond_constraints_view_asserted_complete(
        #[case] constraints: NoncovalentBondConstraintsForm,
        #[case] expected: Option<NoncovalentBondConstraintForm>,
    ) {
        let mut form = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);
        form.constraints = constraints;
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::O),
            ],
            noncovalent: vec![(AtomId(0), AtomId(1), form)],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondId(0))
                .constraints()
                .asserted_complete(NoncovalentBondConstraintKey::Intramolecular),
            expected
        );
    }

    #[rstest]
    #[case::present(
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))),
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        )))
    )]
    #[case::absent(None, None)]
    fn test_stereo_atom_constraints_view_asserted_complete(
        #[case] set: Option<StereoAtomConstraintForm>,
        #[case] expected: Option<StereoAtomConstraintForm>,
    ) {
        let mut form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
        if let Some(constraint) = set {
            form.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                form,
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .constraints()
                .asserted_complete(StereoAtomConstraintKey::Stereogenicity),
            expected
        );
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_asserted() {
        let mut form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
        form.constraints
            .set(StereoAtomConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            ));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                form,
            )],
            ..Default::default()
        });
        let view = molecule.stereo_atom(StereoAtomId(0)).constraints();
        assert_eq!(
            view.asserted(StereoAtomConstraintKey::Stereogenicity),
            Some(&StereoAtomConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            ))
        );
        // No projection is defined for stereo constraint kinds.
        assert_eq!(view.derived(StereoAtomConstraintKey::Stereogenicity), None);
        assert_eq!(
            view.derived_complete(StereoAtomConstraintKey::Stereogenicity),
            None
        );
    }

    #[rstest]
    #[case::present(
        Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))),
        Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        )))
    )]
    #[case::absent(None, None)]
    fn test_stereo_bond_constraints_view_asserted_complete(
        #[case] set: Option<StereoBondConstraintForm>,
        #[case] expected: Option<StereoBondConstraintForm>,
    ) {
        let mut form = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1));
        if let Some(constraint) = set {
            form.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                form,
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule
                .stereo_bond(StereoBondId(0))
                .constraints()
                .asserted_complete(StereoBondConstraintKey::Stereogenicity),
            expected
        );
    }

    #[rstest]
    fn test_stereo_bond_constraints_view_derived() {
        // No projection is defined for stereo constraint kinds.
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        });
        let view = molecule.stereo_bond(StereoBondId(0)).constraints();
        assert_eq!(view.derived(StereoBondConstraintKey::Stereogenicity), None);
        assert_eq!(
            view.derived_complete(StereoBondConstraintKey::Stereogenicity),
            None
        );
        assert!(view.is_empty());
    }
}
