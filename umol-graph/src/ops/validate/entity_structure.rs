//! Tier-2 entity-structure invariants over per-relation entities, independent of a chemistry
//! model.

use std::collections::{BTreeSet, HashMap, HashSet};

use thiserror::Error;
use umol_graph_ir::ir::{AtomId, BondId, Molecule};
use umol_utils::solution::Solution;

/// Structural shape checks on per-relation entities: per-relation participant
/// well-formedness (no self-loops, no duplicate or role-conflicting
/// participants), no same-type parallel relations, aromatic-system disjointness,
/// distinct stereo sites. Cross-type parallelism (a localized and a dative bond on the same atom
/// pair) is permitted. Parallel collection shape is representation integrity and is checked by
/// `Molecule::check_integrity` before this validator runs.
#[derive(Clone, Copy, Debug, Default)]
pub struct EntityStructureInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureInvariantsContradiction {
    #[error("bond: self-loop on atom {atom:?}")]
    BondSelfLoop { atom: AtomId },
    #[error("bond: parallel bonds on atoms {atoms:?}")]
    BondsParallel { atoms: [AtomId; 2] },
    #[error("dative bond: donor {donor:?} duplicated (acceptor {acceptor:?})")]
    DativeBondDonorDuplicate { acceptor: AtomId, donor: AtomId },
    #[error("dative bond: acceptor {atom:?} is also a donor")]
    DativeBondAcceptorIsDonor { atom: AtomId },
    #[error(
        "dative bond: parallel datives to acceptor {acceptor:?} sharing donor {shared_donor:?}"
    )]
    DativeBondsParallel {
        acceptor: AtomId,
        shared_donor: AtomId,
    },
    #[error("noncovalent bond: self-loop on atom {atom:?}")]
    NoncovalentBondSelfLoop { atom: AtomId },
    #[error("noncovalent bond: parallel bonds on atoms {atoms:?}")]
    NoncovalentBondsParallel { atoms: [AtomId; 2] },
    #[error("aromatic system: participant {atom:?} duplicated")]
    AromaticSystemDuplicateParticipant { atom: AtomId },
    #[error("aromatic systems: overlap on atom {atom:?}")]
    AromaticSystemsOverlap { atom: AtomId },
    #[error("multicenter bond: participant {atom:?} duplicated")]
    MulticenterBondDuplicateParticipant { atom: AtomId },
    #[error("multicenter bonds: identical participant set {atoms:?}")]
    MulticenterBondsIdentical { atoms: Vec<AtomId> },
    #[error("stereo atom: duplicate site {atom:?}")]
    StereoAtomSitesDuplicate { atom: AtomId },
    #[error("stereo bond: duplicate site {bond:?}")]
    StereoBondSitesDuplicate { bond: BondId },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureInvariantsError {}

impl EntityStructureInvariantsValidator {
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), EntityStructureInvariantsContradiction>, EntityStructureInvariantsError>
    {
        let contradiction = bond_structure_check(molecule)
            .or_else(|| dative_structure_check(molecule))
            .or_else(|| noncovalent_structure_check(molecule))
            .or_else(|| aromatic_structure_check(molecule))
            .or_else(|| multicenter_structure_check(molecule))
            .or_else(|| stereo_structure_check(molecule));
        match contradiction {
            Some(c) => Ok(Solution::Contradictory(c)),
            None => Ok(Solution::Determined(())),
        }
    }
}

/// Localized bonds: no self-loop, no two bonds on the same unordered atom pair.
/// Each atom's neighbors are stored sorted by atom id (CSR invariant), so a
/// self-loop is a neighbor equal to the atom and a parallel bond is two adjacent
/// equal neighbors — a single linear scan over the adjacency, no auxiliary set.
fn bond_structure_check(molecule: &Molecule) -> Option<EntityStructureInvariantsContradiction> {
    if !molecule.bonds().has_conflict() {
        return None;
    }
    for atom in molecule.atoms().ids() {
        let mut prev: Option<AtomId> = None;
        for neighbor in molecule.neighbors(atom) {
            let other = neighbor.atom_id();
            if other == atom {
                return Some(EntityStructureInvariantsContradiction::BondSelfLoop { atom });
            }
            if prev == Some(other) {
                let pair = if atom <= other {
                    [atom, other]
                } else {
                    [other, atom]
                };
                return Some(EntityStructureInvariantsContradiction::BondsParallel { atoms: pair });
            }
            prev = Some(other);
        }
    }
    None
}

/// Dative bonds: donors distinct, acceptor not among donors, and for any shared
/// acceptor the donor sets are vertex-disjoint.
fn dative_structure_check(molecule: &Molecule) -> Option<EntityStructureInvariantsContradiction> {
    if !molecule.dative_bonds().has_conflict() {
        return None;
    }
    let mut donors_by_acceptor: HashMap<AtomId, HashSet<AtomId>> = HashMap::new();
    for d in molecule.dative_bonds().iter() {
        let acceptor = d.acceptor_id();
        let mut donors: HashSet<AtomId> = HashSet::new();
        for donor in d.donor_ids() {
            if donor == acceptor {
                return Some(
                    EntityStructureInvariantsContradiction::DativeBondAcceptorIsDonor {
                        atom: acceptor,
                    },
                );
            }
            if !donors.insert(donor) {
                return Some(
                    EntityStructureInvariantsContradiction::DativeBondDonorDuplicate {
                        acceptor,
                        donor,
                    },
                );
            }
        }
        let accumulated = donors_by_acceptor.entry(acceptor).or_default();
        for &donor in &donors {
            if accumulated.contains(&donor) {
                return Some(
                    EntityStructureInvariantsContradiction::DativeBondsParallel {
                        acceptor,
                        shared_donor: donor,
                    },
                );
            }
        }
        accumulated.extend(donors);
    }
    None
}

/// Noncovalent bonds: endpoints distinct, at most one interaction per unordered atom pair
/// (parallel bonds of any kind are forbidden — the uniqueness that makes structural refs
/// unambiguous).
fn noncovalent_structure_check(
    molecule: &Molecule,
) -> Option<EntityStructureInvariantsContradiction> {
    if !molecule.noncovalent_bonds().has_conflict() {
        return None;
    }
    let mut seen: HashSet<[AtomId; 2]> = HashSet::new();
    for nc in molecule.noncovalent_bonds().iter() {
        let [a, b] = nc.atom_ids();
        if a == b {
            return Some(
                EntityStructureInvariantsContradiction::NoncovalentBondSelfLoop { atom: a },
            );
        }
        let pair = if a <= b { [a, b] } else { [b, a] };
        if !seen.insert(pair) {
            return Some(
                EntityStructureInvariantsContradiction::NoncovalentBondsParallel { atoms: pair },
            );
        }
    }
    None
}

/// Aromatic systems: participants distinct within a system and systems pairwise vertex-disjoint.
/// The disjointness conflict is the per-entity `has_conflict` primitive; the detailed contradiction
/// locates the offending atom.
fn aromatic_structure_check(molecule: &Molecule) -> Option<EntityStructureInvariantsContradiction> {
    if molecule.aromatic_systems().has_conflict() {
        let mut global: HashSet<AtomId> = HashSet::new();
        for view in molecule.aromatic_systems().iter() {
            let mut local: HashSet<AtomId> = HashSet::new();
            for atom in view.atom_ids() {
                if !local.insert(atom) {
                    return Some(
                        EntityStructureInvariantsContradiction::AromaticSystemDuplicateParticipant { atom },
                    );
                }
                if global.contains(&atom) {
                    return Some(
                        EntityStructureInvariantsContradiction::AromaticSystemsOverlap { atom },
                    );
                }
            }
            global.extend(local);
        }
    }
    None
}

/// Multicenter bonds: participants distinct within a bond, and no two bonds with an identical
/// participant set (partial overlap allowed). The duplicate/identical conflict is the per-entity
/// `has_conflict` primitive; the detailed contradiction locates the offender.
fn multicenter_structure_check(
    molecule: &Molecule,
) -> Option<EntityStructureInvariantsContradiction> {
    if molecule.multicenter_bonds().has_conflict() {
        let mut seen_sets: HashSet<BTreeSet<AtomId>> = HashSet::new();
        for view in molecule.multicenter_bonds().iter() {
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            let mut set: BTreeSet<AtomId> = BTreeSet::new();
            for &atom in &atoms {
                if !set.insert(atom) {
                    return Some(
                        EntityStructureInvariantsContradiction::MulticenterBondDuplicateParticipant { atom },
                    );
                }
            }
            if !seen_sets.insert(set) {
                return Some(
                    EntityStructureInvariantsContradiction::MulticenterBondsIdentical { atoms },
                );
            }
        }
    }
    None
}

/// Stereo overlays: stereo-atom sites pairwise distinct, stereo-bond sites pairwise distinct. The
/// conflict predicate is the per-entity `has_conflict` primitive (also consulted by `apply_at` /
/// `meet_pushout`); the detailed contradiction locates the offending site.
fn stereo_structure_check(molecule: &Molecule) -> Option<EntityStructureInvariantsContradiction> {
    if molecule.stereo_atoms().has_conflict() {
        let mut sites: HashSet<AtomId> = HashSet::new();
        let atom = molecule
            .stereo_atoms()
            .iter()
            .map(|sp| sp.site_id())
            .find(|&atom| !sites.insert(atom))
            .expect("has_conflict reports a repeated stereo-atom site");
        return Some(EntityStructureInvariantsContradiction::StereoAtomSitesDuplicate { atom });
    }
    if molecule.stereo_bonds().has_conflict() {
        let mut sites: HashSet<BondId> = HashSet::new();
        let bond = molecule
            .stereo_bonds()
            .iter()
            .map(|sp| sp.site_id())
            .find(|&bond| !sites.insert(bond))
            .expect("has_conflict reports a repeated stereo-bond site");
        return Some(EntityStructureInvariantsContradiction::StereoBondSitesDuplicate { bond });
    }
    None
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AtomId, BondId, Molecule};
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::aromatic_system(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "*"}]}"#))]
    #[case::cross_type_parallel(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#))]
    #[case::dative_shared_acceptor_disjoint_donors(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"} {:donors [2] :acceptor 0 :attrs "1"}]}"#))]
    #[case::dative_shared_donors_distinct_acceptors(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donors [2] :acceptor 0 :attrs "1"} {:donors [2] :acceptor 1 :attrs "1"}]}"#))]
    #[case::multicenter_partial_overlap(mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "*"} {:atoms [1 2 3] :attrs "*"}]}"#))]
    fn test_entity_structure_validator_validate(#[case] molecule: Molecule) {
        assert_eq!(
            EntityStructureInvariantsValidator
                .validate(&molecule)
                .unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    #[case::bond_self_loop(
        mol_dsl!(r#"{:atoms ["C"] :bonds [[0 0 "1"]]}"#),
        EntityStructureInvariantsContradiction::BondSelfLoop { atom: AtomId(0) }
    )]
    #[case::bonds_parallel(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"] [0 1 "1"]]}"#),
        EntityStructureInvariantsContradiction::BondsParallel { atoms: [AtomId(0), AtomId(1)] }
    )]
    #[case::dative_acceptor_is_donor(
        mol_dsl!(r#"{:atoms ["C"] :bonds [] :dative-bonds [{:donors [0] :acceptor 0 :attrs "1"}]}"#),
        EntityStructureInvariantsContradiction::DativeBondAcceptorIsDonor { atom: AtomId(0) }
    )]
    #[case::dative_parallel(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"} {:donors [1] :acceptor 0 :attrs "1"}]}"#),
        EntityStructureInvariantsContradiction::DativeBondsParallel { acceptor: AtomId(0), shared_donor: AtomId(1) }
    )]
    #[case::noncovalent_self_loop(
        mol_dsl!(r#"{:atoms ["C"] :bonds [] :noncovalent-bonds [{:atoms [0 0] :attrs "Hbd"}]}"#),
        EntityStructureInvariantsContradiction::NoncovalentBondSelfLoop { atom: AtomId(0) }
    )]
    #[case::noncovalent_parallel(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"} {:atoms [0 1] :attrs "Hbd"}]}"#),
        EntityStructureInvariantsContradiction::NoncovalentBondsParallel { atoms: [AtomId(0), AtomId(1)] }
    )]
    #[case::noncovalent_parallel_distinct_kinds(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"} {:atoms [0 1] :attrs "Vdw"}]}"#),
        EntityStructureInvariantsContradiction::NoncovalentBondsParallel { atoms: [AtomId(0), AtomId(1)] }
    )]
    #[case::aromatic_duplicate_participant(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 1] :attrs "*"}]}"#),
        EntityStructureInvariantsContradiction::AromaticSystemDuplicateParticipant { atom: AtomId(1) }
    )]
    #[case::aromatic_overlap(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "*"} {:atoms [1 2] :attrs "*"}]}"#),
        EntityStructureInvariantsContradiction::AromaticSystemsOverlap { atom: AtomId(1) }
    )]
    #[case::multicenter_duplicate_participant(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 1] :attrs "*"}]}"#),
        EntityStructureInvariantsContradiction::MulticenterBondDuplicateParticipant { atom: AtomId(1) }
    )]
    #[case::multicenter_identical(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "*"} {:atoms [0 1 2] :attrs "*"}]}"#),
        EntityStructureInvariantsContradiction::MulticenterBondsIdentical { atoms: vec![AtomId(0), AtomId(1), AtomId(2)] }
    )]
    #[case::stereo_atom_sites_duplicate(
        mol_dsl!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#),
        EntityStructureInvariantsContradiction::StereoAtomSitesDuplicate { atom: AtomId(0) }
    )]
    #[case::stereo_bond_sites_duplicate(
        mol_dsl!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"]] :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"} {:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]}"#),
        EntityStructureInvariantsContradiction::StereoBondSitesDuplicate { bond: BondId(0) }
    )]
    #[case::dative_donor_duplicate(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :dative-bonds [{:donors [1 1] :acceptor 0 :attrs "1"}]}"#),
        EntityStructureInvariantsContradiction::DativeBondDonorDuplicate { acceptor: AtomId(0), donor: AtomId(1) }
    )]
    fn test_entity_structure_validator_validate_error(
        #[case] molecule: Molecule,
        #[case] expected: EntityStructureInvariantsContradiction,
    ) {
        assert_eq!(
            EntityStructureInvariantsValidator
                .validate(&molecule)
                .unwrap(),
            Solution::Contradictory(expected)
        );
    }
}
