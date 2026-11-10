//! Stereo sites, selection, and structural perception.
//!
//! [`StereoPerceiver`] realizes `#T`/`#C` evidence, recognizes missing
//! tetrahedral and cis/trans entities, and classifies existing arrangements.
//! Resolver inconsistency policy and edit construction remain separate.

mod recognize;

use std::collections::BTreeSet;

pub use recognize::StereoPerceptionError;
use smallvec::SmallVec;
use thiserror::Error;
use umol_graph_ir::ir::{
    AsLit, AtomConstraintKey, AtomId, BondConstraintKey, BondId, CisTransStereoForm, Lattice,
    Molecule, Normalize, NumForm, RingConfig, RingModel, RingSet, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    Stereogenicity, TetrahedralStereoForm,
};
use umol_perm::Permutation;
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;

/// An atom or bond topology site for stereo processing.
///
/// The id refers to an input molecule's atom or bond table, independently of
/// whether a stereo entity exists there. It is not bound to a molecule.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoSite {
    Atom(AtomId),
    Bond(BondId),
}

/// Site, kind, and stereogenicity filters for stereo processing.
///
/// Each None is unrestricted; an empty set selects nothing. The filters
/// intersect. Default leaves all three unrestricted. Sets do not specify
/// processing order.
///
/// This is an open configuration value. Construction does not check site ids
/// against a molecule. The executing operation checks them before indexing or
/// mutation and obtains the kind and stereogenicity from the molecular context.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoSelection {
    pub sites: Option<BTreeSet<StereoSite>>,
    pub kinds: Option<BTreeSet<StereoKind>>,
    pub stereogenicities: Option<BTreeSet<Stereogenicity>>,
}

impl StereoSelection {
    /// Whether the supplied site, kind, and stereogenicity satisfy every filter.
    ///
    /// A missing kind or stereogenicity matches only an unrestricted filter for
    /// that field. This predicate does not look up ids or calculate stereo.
    ///
    /// # Semantic properties
    ///
    /// None accepts every value in its field, including a missing value. An
    /// empty set rejects every value. A match requires all three field filters
    /// to accept their supplied values, without changing the selection.
    pub fn matches(
        &self,
        site: StereoSite,
        kind: Option<StereoKind>,
        stereogenicity: Option<Stereogenicity>,
    ) -> bool {
        self.sites
            .as_ref()
            .is_none_or(|sites| sites.contains(&site))
            && self
                .kinds
                .as_ref()
                .is_none_or(|kinds| kind.is_some_and(|kind| kinds.contains(&kind)))
            && self
                .stereogenicities
                .as_ref()
                .is_none_or(|classes| stereogenicity.is_some_and(|class| classes.contains(&class)))
    }
}

/// Derived stereo entities and evidence inconsistencies.
///
/// Rows describe requested additions and supported stored frames, including
/// unresolved frames. Stored rows retain their ligand order and attributes.
/// Unique completions identify generic virtual kinds; unresolved rows retain
/// Virtual without an inconsistency. The derive
/// outcome reports determinacy explicitly, including unresolved typed frames.
/// Recognition returns only additions. This is an open carrier; independently
/// assembled rows are not checked against a molecule at construction.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation {
    /// Requested or retained tetrahedral stereo atoms.
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>,
    /// Requested or retained cis-trans stereo bonds.
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>,
    /// Constraint failures, entity failures, and independently valid mismatches.
    pub inconsistencies: Vec<StereoInconsistency>,
}

/// Policy-free classification of stereo constraint and entity inconsistencies.
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoInconsistency {
    #[error("tetrahedral stereo constraint at atom {atom:?} cannot be realized")]
    TetrahedralStereoFailure { atom: AtomId },
    #[error("stereo atom {stereo_atom:?} cannot be realized")]
    StereoAtomFailure { stereo_atom: StereoAtomId },
    #[error(
        "tetrahedral stereo constraint at atom {atom:?} disagrees with stereo atom {stereo_atom:?}"
    )]
    TetrahedralStereoMismatch {
        atom: AtomId,
        stereo_atom: StereoAtomId,
    },
    #[error("cis-trans stereo constraint at bond {bond:?} cannot be realized")]
    CisTransStereoFailure { bond: BondId },
    #[error("stereo bond {stereo_bond:?} cannot be realized")]
    StereoBondFailure { stereo_bond: StereoBondId },
    #[error(
        "cis-trans stereo constraint at bond {bond:?} disagrees with stereo bond {stereo_bond:?}"
    )]
    CisTransStereoMismatch {
        bond: BondId,
        stereo_bond: StereoBondId,
    },
}

/// A stereo inconsistency rejected by an operation.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoContradiction {
    #[error("stereo inconsistency: {0}")]
    Inconsistency(#[from] StereoInconsistency),
}

/// Stereo evidence, recognition, and classification under a selected chemistry model.
///
/// derive realizes #T/#C evidence and checks stored arrangements; recognize
/// additionally identifies missing selected tetrahedral/cis-trans entities.
/// classify evaluates existing arrangements without choosing their cosets.
/// All queries preserve the molecule. Resolver policy and transformation remain
/// separate from these read-only results. A stored arrangement may be symmetric
/// or prochiral, independently of its configuration being determined.
#[derive(Clone, Debug)]
pub struct StereoPerceiver {
    model: StereoModel,
}

impl StereoPerceiver {
    /// Capture a snapshot of the supplied stereo model.
    pub fn new(model: &StereoModel) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Derive requested stereo relations and identify stored virtual ligand kinds.
    /// Positive evidence supplies a kind with an undetermined coset; negative evidence
    /// requests no entity.
    /// Stored rows preserve their configuration, ligand order, and entity constraints.
    /// Ambiguous identification or insufficient evidence yields Underdetermined,
    /// retaining unresolved rows and established findings. Impossible completions
    /// produce entity failures. An undetermined coset does not imply ambiguity.
    /// Without a stored entity, insufficient frame evidence contributes no row.
    /// Contradictory evidence between realizable sources is a mismatch. Independent
    /// constraint and entity realizability failures retain their separate diagnostics.
    /// Missing stereo bonds on double bonds follow the model's minimum ring size.
    /// Existing entities are checked regardless of ring size. Ring detection uses ring_config.
    ///
    /// # Semantic properties
    ///
    /// The molecule is unchanged. Stored rows differ only at Virtual ligand kinds
    /// uniquely identified from compatible frames; fixed ligand kinds remain exact.
    pub fn derive(
        &self,
        molecule: &Molecule,
        ring_config: RingConfig,
    ) -> Solution<StereoDerivation, StereoContradiction> {
        let mut stereo_atoms = Vec::new();
        let mut stereo_bonds = Vec::new();
        let mut inconsistencies = BTreeSet::new();
        let mut underdetermined = false;
        let mut rings = None;

        for view in molecule.atoms().iter() {
            let atom = view.id();
            let relation = molecule.stereo_atoms().at(atom);
            let constraints = view.constraints();
            let asserted = constraints.tetrahedral_stereo();
            if relation.is_none() && !matches!(asserted, TetrahedralStereoForm::Stereo) {
                continue;
            }
            let candidate = if matches!(asserted, TetrahedralStereoForm::Stereo) {
                self.derive_stereo_atom(molecule, atom, &StereoCoset::Undetermined)
            } else {
                None
            };
            let completed = relation.and_then(|relation| {
                if relation.attributes().configuration.kind() != Some(StereoKind::Tetrahedral) {
                    return None;
                }
                let (ligands, unresolved) = match self.identify_virtual_ligands(
                    molecule,
                    StereoSite::Atom(atom),
                    relation.ligand_ids(),
                ) {
                    Solution::Determined(ligands) => (ligands, false),
                    Solution::Underdetermined(ligands) => (ligands, true),
                    Solution::Contradictory(()) => return None,
                };
                Some((ligands, unresolved))
            });
            if candidate.is_none()
                && completed.is_none()
                && matches!(asserted, TetrahedralStereoForm::Stereo)
            {
                if self.has_compatible_stereo_frame(molecule, StereoSite::Atom(atom)) {
                    underdetermined = true;
                } else {
                    inconsistencies.insert(StereoInconsistency::TetrahedralStereoFailure { atom });
                }
            }
            if let Some(relation) = relation {
                if constraints
                    .evidenced(AtomConstraintKey::TetrahedralStereo)
                    .is_none()
                {
                    inconsistencies.insert(StereoInconsistency::TetrahedralStereoMismatch {
                        atom,
                        stereo_atom: relation.id(),
                    });
                    continue;
                }
                let Some((ligands, unresolved)) = completed else {
                    inconsistencies.insert(StereoInconsistency::StereoAtomFailure {
                        stereo_atom: relation.id(),
                    });
                    continue;
                };
                underdetermined |= unresolved;
                stereo_atoms.push((atom, ligands, relation.attributes().clone()));
                continue;
            }
            if let Some((ligands, stereo)) = candidate {
                stereo_atoms.push((atom, ligands, stereo));
            }
        }

        for view in molecule.bonds().iter() {
            let bond = view.id();
            let relation = molecule.stereo_bonds().at(bond);
            let constraints = view.constraints();
            let asserted = constraints.cis_trans_stereo();
            if relation.is_none() && !matches!(asserted, CisTransStereoForm::Stereo) {
                continue;
            }
            if relation.is_none()
                && !self.stereo_bond_ring_eligible(molecule, bond, ring_config, &mut rings)
            {
                continue;
            }
            let candidate = if matches!(asserted, CisTransStereoForm::Stereo) {
                self.derive_stereo_bond(molecule, bond, &StereoCoset::Undetermined)
            } else {
                None
            };
            let completed = relation.and_then(|relation| {
                if relation.attributes().configuration.kind() != Some(StereoKind::CisTrans) {
                    return None;
                }
                let (ligands, unresolved) = match self.identify_virtual_ligands(
                    molecule,
                    StereoSite::Bond(bond),
                    relation.ligand_ids(),
                ) {
                    Solution::Determined(ligands) => (ligands, false),
                    Solution::Underdetermined(ligands) => (ligands, true),
                    Solution::Contradictory(()) => return None,
                };
                Some((ligands, unresolved))
            });
            if candidate.is_none()
                && completed.is_none()
                && matches!(asserted, CisTransStereoForm::Stereo)
            {
                if self.has_compatible_stereo_frame(molecule, StereoSite::Bond(bond)) {
                    underdetermined = true;
                } else {
                    inconsistencies.insert(StereoInconsistency::CisTransStereoFailure { bond });
                }
            }
            if let Some(relation) = relation {
                if constraints
                    .evidenced(BondConstraintKey::CisTransStereo)
                    .is_none()
                {
                    inconsistencies.insert(StereoInconsistency::CisTransStereoMismatch {
                        bond,
                        stereo_bond: relation.id(),
                    });
                    continue;
                }
                let Some((ligands, unresolved)) = completed else {
                    inconsistencies.insert(StereoInconsistency::StereoBondFailure {
                        stereo_bond: relation.id(),
                    });
                    continue;
                };
                underdetermined |= unresolved;
                stereo_bonds.push((bond, ligands, relation.attributes().clone()));
                continue;
            }
            if let Some((ligands, stereo)) = candidate {
                stereo_bonds.push((bond, ligands, stereo));
            }
        }

        let derivation = StereoDerivation {
            stereo_atoms,
            stereo_bonds,
            inconsistencies: inconsistencies.into_iter().collect(),
        };
        if underdetermined {
            Solution::Underdetermined(derivation)
        } else {
            Solution::Determined(derivation)
        }
    }

    fn has_compatible_stereo_frame(&self, molecule: &Molecule, site: StereoSite) -> bool {
        let mut sites = SmallVec::<[(AtomId, Option<AtomId>); 2]>::new();
        let degree = match site {
            StereoSite::Atom(atom) => {
                sites.push((atom, None));
                4
            }
            StereoSite::Bond(bond) => {
                let [first, second] = molecule.bond(bond).atom_ids();
                sites.extend([(first, Some(second)), (second, Some(first))]);
                2
            }
        };
        let mut ligands = Vec::new();
        for (atom, partner) in sites {
            let start = ligands.len();
            for neighbor in molecule.atom(atom).neighbors() {
                if Some(neighbor.atom_id()) == partner {
                    continue;
                }
                let ligand = StereoLigand::new(neighbor.atom_id(), StereoLigandKind::Atom);
                if ligands.contains(&ligand) {
                    return false;
                }
                ligands.push(ligand);
            }
            let actual = ligands.len() - start;
            if actual > degree {
                return false;
            }
            match degree - actual {
                0 => {}
                1 => ligands.push(StereoLigand::new(atom, StereoLigandKind::Virtual)),
                2 => {
                    ligands.push(StereoLigand::new(atom, StereoLigandKind::ImplicitHydrogen));
                    ligands.push(StereoLigand::new(atom, StereoLigandKind::LonePair));
                }
                _ => return false,
            }
        }
        !matches!(
            self.identify_virtual_ligands(molecule, site, &ligands),
            Solution::Contradictory(())
        )
    }

    fn stereo_bond_ring_eligible(
        &self,
        molecule: &Molecule,
        bond: BondId,
        config: RingConfig,
        rings: &mut Option<RingSet>,
    ) -> bool {
        let stereo_bond_min_ring_size = self.model.stereo_bond_min_ring_size;
        if stereo_bond_min_ring_size == 0 || molecule.bond(bond).order().as_lit() != Some(2) {
            return true;
        }
        let rings = rings.get_or_insert_with(|| {
            let mut model = RingModel::default();
            model.max_ring_size = model
                .max_ring_size
                .max(stereo_bond_min_ring_size as usize - 1);
            molecule.rings(model, config).into_ring_set()
        });
        rings
            .bond_smallest_ring_size(bond)
            .is_none_or(|size| size >= stereo_bond_min_ring_size as usize)
    }

    /// Derive a tetrahedral stereo atom in the canonical neighbor frame.
    /// All actual neighbors precede virtual hydrogen, then lone pair.
    /// Returns None for an excluded site, a hydrogen or lone-pair count other than
    /// literal zero or one, or a ligand list that does not have four entries.
    pub fn derive_stereo_atom(
        &self,
        molecule: &Molecule,
        atom: AtomId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoAtomForm)> {
        let view = molecule.atom(atom);
        if view.is_in_aromatic_system() {
            return None;
        }

        let kind = StereoKind::Tetrahedral;
        let model = self.model.kind_model(kind)?;
        if !model.scope.contains(view.element().as_lit()?) {
            return None;
        }

        let mut ligands: Vec<StereoLigand> = view
            .neighbors()
            .map(|neighbor| StereoLigand::new(neighbor.atom_id(), StereoLigandKind::Atom))
            .collect();
        for (kind, count) in [
            (
                StereoLigandKind::ImplicitHydrogen,
                view.implicit_hydrogens().as_lit()?,
            ),
            (StereoLigandKind::LonePair, view.lone_pairs().as_lit()?),
        ] {
            match count {
                0 => {}
                1 => ligands.push(StereoLigand::new(atom, kind)),
                _ => return None,
            }
        }
        if ligands.len() != kind.degree() {
            return None;
        }

        Some((ligands, StereoAtomForm::new(kind, coset.clone())))
    }

    /// Derive a cis-trans stereo bond in the canonical endpoint and side frames.
    /// Each side lists every actual neighbor except the opposite endpoint, then
    /// virtual hydrogen and lone pair. Returns None for an excluded kind or element,
    /// a hydrogen or lone-pair count other than literal zero or one, or an
    /// endpoint ligand list that does not have two distinct entries. The endpoint
    /// frames must not share a ligand.
    pub fn derive_stereo_bond(
        &self,
        molecule: &Molecule,
        bond: BondId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoBondForm)> {
        let view = molecule.bond(bond);
        let kind = StereoKind::CisTrans;
        let model = self.model.kind_model(kind)?;
        let [first, second] = view.atom_ids();
        if !model
            .scope
            .contains(molecule.atom(first).element().as_lit()?)
            || !model
                .scope
                .contains(molecule.atom(second).element().as_lit()?)
        {
            return None;
        }

        let first_side = self.bond_side_ligands(molecule, first, second)?;
        let second_side = self.bond_side_ligands(molecule, second, first)?;
        if first_side.iter().any(|ligand| second_side.contains(ligand)) {
            return None;
        }
        let ligands = vec![first_side[0], first_side[1], second_side[0], second_side[1]];
        Some((ligands, StereoBondForm::new(kind, coset.clone())))
    }

    fn identify_virtual_ligands(
        &self,
        molecule: &Molecule,
        site: StereoSite,
        stored: &[StereoLigand],
    ) -> Solution<Vec<StereoLigand>, ()> {
        let mut sites = SmallVec::<[(AtomId, Option<AtomId>); 2]>::new();
        let kind = match site {
            StereoSite::Atom(atom) => {
                if molecule.atom(atom).is_in_aromatic_system() {
                    return Solution::Contradictory(());
                }
                sites.push((atom, None));
                StereoKind::Tetrahedral
            }
            StereoSite::Bond(bond) => {
                let [first, second] = molecule.bond(bond).atom_ids();
                sites.extend([(first, Some(second)), (second, Some(first))]);
                StereoKind::CisTrans
            }
        };
        let Some(model) = self.model.kind_model(kind) else {
            return Solution::Contradictory(());
        };
        let mut unresolved_scope = false;
        let mut frames = SmallVec::<[(AtomId, Vec<StereoLigand>); 2]>::new();
        for (atom, partner) in sites {
            let view = molecule.atom(atom);
            match view.element().as_lit() {
                Some(element) if !model.scope.contains(element) => {
                    return Solution::Contradictory(());
                }
                Some(_) => {}
                None => unresolved_scope = true,
            }
            let actual = view
                .neighbors()
                .map(|neighbor| neighbor.atom_id())
                .filter(|neighbor| Some(*neighbor) != partner)
                .map(|neighbor| StereoLigand::new(neighbor, StereoLigandKind::Atom))
                .collect();
            frames.push((atom, actual));
        }
        let positions: SmallVec<[usize; 2]> = stored
            .iter()
            .enumerate()
            .filter_map(|(position, ligand)| {
                (ligand.kind == StereoLigandKind::Virtual).then_some(position)
            })
            .collect();
        let mut completion = None;
        let mut underdetermined = false;
        'assignments: for assignment in 0..(1 << positions.len()) {
            let mut ligands = stored.to_vec();
            for (bit, &position) in positions.iter().enumerate() {
                ligands[position].kind = if assignment & (1 << bit) == 0 {
                    StereoLigandKind::ImplicitHydrogen
                } else {
                    StereoLigandKind::LonePair
                };
            }
            let mut expected = Vec::with_capacity(kind.degree());
            let mut unresolved_count = false;
            for (atom, actual) in &frames {
                let side_start = expected.len();
                expected.extend_from_slice(actual);
                let view = molecule.atom(*atom);
                for (ligand_kind, count) in [
                    (
                        StereoLigandKind::ImplicitHydrogen,
                        view.implicit_hydrogens(),
                    ),
                    (StereoLigandKind::LonePair, view.lone_pairs()),
                ] {
                    let required = ligands
                        .iter()
                        .filter(|ligand| ligand.atom_id == *atom && ligand.kind == ligand_kind)
                        .count();
                    if required > 1 {
                        continue 'assignments;
                    }
                    let Ok(count) = count.clone().normalize() else {
                        continue 'assignments;
                    };
                    if matches!(count, NumForm::ArithExpr(_) | NumForm::PredExpr(_)) {
                        unresolved_count = true;
                    } else if !count.is_compatible(&NumForm::Lit(required as i64)) {
                        continue 'assignments;
                    }
                    if required == 1 {
                        expected.push(StereoLigand::new(*atom, ligand_kind));
                    }
                }
                if kind == StereoKind::CisTrans && expected.len() - side_start != 2 {
                    continue 'assignments;
                }
            }
            if expected.len() != kind.degree() {
                continue;
            }
            let Some(action) = Permutation::between(&ligands, &expected) else {
                continue;
            };
            if !kind.class_key().space().allows(action) {
                continue;
            }
            if unresolved_count || unresolved_scope {
                underdetermined = true;
            } else if completion.is_some() {
                return Solution::Underdetermined(stored.to_vec());
            } else {
                completion = Some(ligands);
            }
        }
        if underdetermined {
            Solution::Underdetermined(stored.to_vec())
        } else if let Some(ligands) = completion {
            Solution::Determined(ligands)
        } else {
            Solution::Contradictory(())
        }
    }

    fn bond_side_ligands(
        &self,
        molecule: &Molecule,
        atom: AtomId,
        partner: AtomId,
    ) -> Option<[StereoLigand; 2]> {
        let view = molecule.atom(atom);
        let mut ligands: SmallVec<[StereoLigand; 2]> = view
            .neighbors()
            .map(|neighbor| neighbor.atom_id())
            .filter(|&neighbor| neighbor != partner)
            .map(|neighbor| StereoLigand::new(neighbor, StereoLigandKind::Atom))
            .collect();
        for (kind, count) in [
            (
                StereoLigandKind::ImplicitHydrogen,
                view.implicit_hydrogens().as_lit()?,
            ),
            (StereoLigandKind::LonePair, view.lone_pairs().as_lit()?),
        ] {
            match count {
                0 => {}
                1 => ligands.push(StereoLigand::new(atom, kind)),
                _ => return None,
            }
        }
        ligands.into_inner().ok()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_ir::dsl::StereoAtomDsl;
    use umol_graph_ir::ir::{
        ArithExpr, AtomConstraintForm, AtomForm, BondConstraintForm, BondForm, ElementForm,
        MoleculeEntries, NumForm, StereoAtomId, StereoBondId,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::{ElementScope, StereoKindModel};

    #[rstest]
    #[case::unrestricted(
        StereoSelection::default(),
        vec![StereoSite::Bond(BondId(0)), StereoSite::Atom(AtomId(2)),
            StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(4)),
            StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(3))],
    )]
    #[case::empty_sites(
        StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() },
        vec![],
    )]
    #[case::empty_kinds(
        StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() },
        vec![],
    )]
    #[case::empty_stereogenicities(
        StereoSelection { stereogenicities: Some(BTreeSet::new()), ..Default::default() },
        vec![],
    )]
    #[case::atom(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::bond(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])),
            ..Default::default() },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::sites(
        StereoSelection { sites: Some(BTreeSet::from([
            StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(0)),
        ])), ..Default::default() },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Atom(AtomId(1))],
    )]
    #[case::tetrahedral(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1))],
    )]
    #[case::cis_trans(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
            ..Default::default() },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::square_planar(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::SquarePlanar])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(3))],
    )]
    #[case::stereogenic(
        StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            ..Default::default() },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Atom(AtomId(0))],
    )]
    #[case::prochiral(
        StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(4)), StereoSite::Atom(AtomId(1))],
    )]
    #[case::symmetric(
        StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Atom(AtomId(3))],
    )]
    #[case::site_kind_intersection(
        StereoSelection {
            sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(0))])),
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
            stereogenicities: None,
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    #[case::intersection(
        StereoSelection {
            sites: Some(BTreeSet::from([
                StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
                StereoSite::Atom(AtomId(3)), StereoSite::Bond(BondId(0)), StereoSite::Atom(AtomId(4)),
            ])),
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    fn test_stereo_selection_matches(
        #[case] selection: StereoSelection,
        #[case] expected: Vec<StereoSite>,
    ) {
        let candidates = [
            (
                StereoSite::Bond(BondId(0)),
                StereoKind::CisTrans,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoKind::Tetrahedral,
                Stereogenicity::Symmetric,
            ),
            (
                StereoSite::Atom(AtomId(0)),
                StereoKind::Tetrahedral,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(4)),
                StereoKind::Axial,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoKind::Tetrahedral,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoKind::SquarePlanar,
                Stereogenicity::Symmetric,
            ),
        ];
        let selected: Vec<_> = candidates
            .into_iter()
            .filter(|&(site, kind, stereogenicity)| {
                selection.matches(site, Some(kind), Some(stereogenicity))
            })
            .map(|(site, _, _)| site)
            .collect();
        assert_eq!(selected, expected);
    }

    #[rstest]
    #[case::unknown_values(
        StereoSelection::default(),
        StereoSite::Atom(AtomId(0)),
        None,
        None,
        true
    )]
    #[case::unknown_kind(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() },
        StereoSite::Atom(AtomId(0)), None, Some(Stereogenicity::Stereogenic), false,
    )]
    #[case::unknown_stereogenicity(
        StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            ..Default::default() },
        StereoSite::Atom(AtomId(0)), Some(StereoKind::Tetrahedral), None, false,
    )]
    #[case::known_kind(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() },
        StereoSite::Atom(AtomId(0)), Some(StereoKind::Tetrahedral), None, true,
    )]
    #[case::known_stereogenicity(
        StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
            ..Default::default() },
        StereoSite::Atom(AtomId(0)), None, Some(Stereogenicity::Prochiral), true,
    )]
    #[case::empty_kinds(
        StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() },
        StereoSite::Atom(AtomId(0)), None, None, false,
    )]
    #[case::empty_stereogenicities(
        StereoSelection { stereogenicities: Some(BTreeSet::new()), ..Default::default() },
        StereoSite::Atom(AtomId(0)), None, None, false,
    )]
    #[case::unbound_atom(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(u32::MAX))])),
            ..Default::default() },
        StereoSite::Atom(AtomId(u32::MAX)), None, None, true,
    )]
    #[case::unbound_bond(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(u32::MAX))])),
            ..Default::default() },
        StereoSite::Bond(BondId(u32::MAX)), None, None, true,
    )]
    #[case::excluded_atom(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(u32::MAX))])),
            ..Default::default() },
        StereoSite::Atom(AtomId(0)), None, None, false,
    )]
    #[case::excluded_bond(
        StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(u32::MAX))])),
            ..Default::default() },
        StereoSite::Bond(BondId(0)), None, None, false,
    )]
    fn test_stereo_selection_matches_values(
        #[case] selection: StereoSelection,
        #[case] site: StereoSite,
        #[case] kind: Option<StereoKind>,
        #[case] stereogenicity: Option<Stereogenicity>,
        #[case] expected: bool,
    ) {
        assert_eq!(selection.matches(site, kind, stereogenicity), expected);
    }

    #[rstest]
    #[case::materializations(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T+" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C+"] [6 7 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )],
            stereo_bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )],

            inconsistencies: vec![],
        },
    )]
    #[case::existing_elements(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2"] [6 7 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :attrs "Th1"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :attrs "Ct1"
            }]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],

            inconsistencies: vec![],
        },
    )]
    #[case::unrealizable_sites(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h0#T+" "C#h3"
                    "C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [3 4 "1"] [4 5 "2#C+"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                StereoInconsistency::CisTransStereoFailure { bond: BondId(3) },
            ],
        },
    )]
    #[case::relation_mismatches(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C!"] [6 7 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :attrs "Th0"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :attrs "Ct0"
            }]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoMismatch {
                    atom: AtomId(1),
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::CisTransStereoMismatch {
                    bond: BondId(4),
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    #[case::entity_failures(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"
                    "C#h3" "C#h0" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Sp1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![
                StereoInconsistency::StereoAtomFailure {
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::StereoBondFailure {
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    #[case::virtual_ligand_failures(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2"] [4 6 "1"] [5 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th1"}]
            :stereo-bonds [{:site 3 :ligands [6 [:lp 4] 7 [:h 5]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            inconsistencies: vec![
                StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) },
                StereoInconsistency::StereoBondFailure { stereo_bond: StereoBondId(0) },
            ],
            ..StereoDerivation::default()
        },
    )]
    #[case::negative_evidence(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#T!" "F" "Cl" "Br" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2#C!"] [4 6 "1"] [5 7 "1"]]
        }"#),
        StereoDerivation::default(),
    )]
    #[case::not_stereo_mismatches(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C!"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoMismatch {
                    atom: AtomId(1),
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::CisTransStereoMismatch {
                    bond: BondId(4),
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    #[case::cyclohexene_asserted(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1" "C#h1" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![],
        },
    )]
    #[case::cyclohexene_asserted_undetermined(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1" "C#h1" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],

            inconsistencies: vec![],
        },
    )]
    fn test_stereo_perceiver_derive(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            Solution::Determined(expected)
        );
    }

    #[rstest]
    #[case::virtual_atom_hydrogen(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:virt 0] 3 2 1] :attrs "Th1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(AtomId(0), vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            ], "Th1#g/#o(0,1)/#p(0,1)!#f(0,1)".parse::<StereoAtomDsl>().unwrap().0)],
            ..StereoDerivation::default()
        },
    )]
    #[case::virtual_atom_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 [:virt 0] 1 3] :attrs "Th*"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(AtomId(0), vec![
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ], StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))],
            ..StereoDerivation::default()
        },
    )]
    #[case::virtual_bond_mixed_kinds(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "N#n" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:virt 1] 3 2 [:virt 0]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_bonds: vec![(BondId(0), vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ], StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)))],
            ..StereoDerivation::default()
        },
    )]
    #[case::tetrahedral_refinement(
        mol_dsl_concrete!(r#"{
            :atoms ["C#T+" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        }"#),
        StereoDerivation {
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
            ..StereoDerivation::default()
        },
    )]
    #[case::cis_trans_refinement(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..StereoDerivation::default()
        },
    )]
    fn test_stereo_perceiver_derive_evidence(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            Solution::Determined(expected),
        );
    }

    #[rstest]
    #[case::hydrogen(
        NumForm::Lit(1),
        NumForm::Lit(0),
        Solution::Determined(StereoLigandKind::ImplicitHydrogen)
    )]
    #[case::lone_pair(
        NumForm::Lit(0),
        NumForm::Lit(1),
        Solution::Determined(StereoLigandKind::LonePair)
    )]
    #[case::open_counts(NumForm::Undetermined, NumForm::Undetermined, Solution::Underdetermined(()))]
    #[case::finite_counts(NumForm::lit_set([0, 1]), NumForm::lit_set([0, 1]), Solution::Underdetermined(()))]
    #[case::excluded_lone_pair(
        NumForm::Undetermined,
        NumForm::Lit(0),
        Solution::Determined(StereoLigandKind::ImplicitHydrogen)
    )]
    #[case::excluded_hydrogen(
        NumForm::Lit(0),
        NumForm::Undetermined,
        Solution::Determined(StereoLigandKind::LonePair)
    )]
    #[case::finite_hydrogen(NumForm::lit_set([0, 1]), NumForm::Lit(0), Solution::Determined(StereoLigandKind::ImplicitHydrogen))]
    #[case::ranged_hydrogen(
        NumForm::range_from(1),
        NumForm::Lit(0),
        Solution::Determined(StereoLigandKind::ImplicitHydrogen)
    )]
    #[case::ground_expression(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(0), ArithExpr::Lit(1)])), NumForm::Lit(0), Solution::Determined(StereoLigandKind::ImplicitHydrogen))]
    #[case::symbolic_count(NumForm::var("hydrogens"), NumForm::Lit(0), Solution::Underdetermined(()))]
    #[case::absent(NumForm::Lit(0), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::both_present(NumForm::Lit(1), NumForm::Lit(1), Solution::Contradictory(()))]
    #[case::two_hydrogens(NumForm::Lit(2), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::excluded_counts(NumForm::range_from(2), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::bottom_count(NumForm::lit_set([]), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::symbolic_with_impossible_count(NumForm::var("hydrogens"), NumForm::Lit(2), Solution::Contradictory(()))]
    fn test_stereo_perceiver_derive_virtual_ligands(
        #[case] hydrogens: NumForm,
        #[case] lone_pairs: NumForm,
        #[case] outcome: Solution<StereoLigandKind, (), ()>,
        #[values(StereoCoset::Lit(1), StereoCoset::Undetermined)] coset: StereoCoset,
    ) {
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(hydrogens)
                    .with_lone_pairs(lone_pairs),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
            ],
            bonds: (1..=3)
                .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, coset),
            )],
            ..Default::default()
        };
        let molecule = Molecule::from_entries(entries.clone());
        let original = molecule.clone();
        let mut stereo_atoms = entries.stereo_atoms;
        let expected = match outcome {
            Solution::Determined(kind) => {
                stereo_atoms[0].1[0].kind = kind;
                Solution::Determined(StereoDerivation {
                    stereo_atoms,
                    ..Default::default()
                })
            }
            Solution::Underdetermined(()) => Solution::Underdetermined(StereoDerivation {
                stereo_atoms,
                ..Default::default()
            }),
            Solution::Contradictory(()) => Solution::Determined(StereoDerivation {
                inconsistencies: vec![StereoInconsistency::StereoAtomFailure {
                    stereo_atom: StereoAtomId(0),
                }],
                ..Default::default()
            }),
        };
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            expected
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::fixed_hydrogen(StereoLigandKind::ImplicitHydrogen, StereoLigandKind::LonePair)]
    #[case::fixed_lone_pair(StereoLigandKind::LonePair, StereoLigandKind::ImplicitHydrogen)]
    fn test_stereo_perceiver_derive_virtual_ligands_fixed(
        #[case] fixed: StereoLigandKind,
        #[case] identified: StereoLigandKind,
    ) {
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::P),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                    StereoLigand::new(AtomId(0), fixed),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        };
        let molecule = Molecule::from_entries(entries.clone());
        let mut stereo_atoms = entries.stereo_atoms;
        stereo_atoms[0].1[1].kind = identified;
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            Solution::Determined(StereoDerivation {
                stereo_atoms,
                ..Default::default()
            }),
        );
    }

    #[rstest]
    #[case::mixed(NumForm::Lit(1), NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(1), Solution::Determined((StereoLigandKind::ImplicitHydrogen, StereoLigandKind::LonePair)))]
    #[case::partial(NumForm::Lit(1), NumForm::Lit(0), NumForm::Undetermined, NumForm::Undetermined, Solution::Underdetermined(()))]
    #[case::finite(NumForm::lit_set([0, 1]), NumForm::Lit(0), NumForm::Lit(0), NumForm::lit_set([0, 1]), Solution::Determined((StereoLigandKind::ImplicitHydrogen, StereoLigandKind::LonePair)))]
    #[case::impossible_endpoint(NumForm::Undetermined, NumForm::Undetermined, NumForm::Lit(0), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::symbolic_with_impossible_endpoint(NumForm::var("hydrogens"), NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(0), Solution::Contradictory(()))]
    fn test_stereo_perceiver_derive_virtual_bond_ligands(
        #[case] first_hydrogens: NumForm,
        #[case] first_lone_pairs: NumForm,
        #[case] second_hydrogens: NumForm,
        #[case] second_lone_pairs: NumForm,
        #[case] outcome: Solution<(StereoLigandKind, StereoLigandKind), (), ()>,
        #[values(StereoCoset::Lit(1), StereoCoset::Undetermined)] coset: StereoCoset,
    ) {
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(first_hydrogens)
                    .with_lone_pairs(first_lone_pairs),
                AtomForm::from_element(Element::N)
                    .with_implicit_hydrogens(second_hydrogens)
                    .with_lone_pairs(second_lone_pairs),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Virtual),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                ],
                StereoBondForm::new(StereoKind::CisTrans, coset),
            )],
            ..Default::default()
        };
        let molecule = Molecule::from_entries(entries.clone());
        let mut stereo_bonds = entries.stereo_bonds;
        let expected = match outcome {
            Solution::Determined((first, second)) => {
                stereo_bonds[0].1[0].kind = second;
                stereo_bonds[0].1[3].kind = first;
                Solution::Determined(StereoDerivation {
                    stereo_bonds,
                    ..Default::default()
                })
            }
            Solution::Underdetermined(()) => Solution::Underdetermined(StereoDerivation {
                stereo_bonds,
                ..Default::default()
            }),
            Solution::Contradictory(()) => Solution::Determined(StereoDerivation {
                inconsistencies: vec![StereoInconsistency::StereoBondFailure {
                    stereo_bond: StereoBondId(0),
                }],
                ..Default::default()
            }),
        };
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            expected
        );
    }

    #[rstest]
    #[case::symbolic_count(ElementForm::Lit(Element::C), NumForm::var("hydrogens"))]
    #[case::unknown_scope(ElementForm::Undetermined, NumForm::Lit(0))]
    fn test_stereo_perceiver_derive_partial(
        #[case] element: ElementForm,
        #[case] hydrogens: NumForm,
        #[values(false, true)] negative_assertion: bool,
    ) {
        let mut entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(hydrogens)
                    .with_lone_pairs(0),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: (1..=4)
                .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        };
        entries.atoms[0].element = element;
        if negative_assertion {
            entries.atoms[0]
                .constraints
                .set(AtomConstraintForm::tetrahedral_stereo(
                    TetrahedralStereoForm::NotStereo,
                ));
        }
        let molecule = Molecule::from_entries(entries.clone());
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            if negative_assertion {
                Solution::Determined(StereoDerivation {
                    inconsistencies: vec![StereoInconsistency::TetrahedralStereoMismatch {
                        atom: AtomId(0),
                        stereo_atom: StereoAtomId(0),
                    }],
                    ..Default::default()
                })
            } else {
                Solution::Underdetermined(StereoDerivation {
                    stereo_atoms: entries.stereo_atoms,
                    ..Default::default()
                })
            },
        );
    }

    #[rstest]
    #[case::open_counts(NumForm::Undetermined, NumForm::Undetermined, true)]
    #[case::open_hydrogen(NumForm::Undetermined, NumForm::Lit(0), true)]
    #[case::open_lone_pair(NumForm::Lit(0), NumForm::Undetermined, true)]
    #[case::finite_counts(NumForm::lit_set([0, 1]), NumForm::lit_set([0, 1]), true)]
    #[case::ranged_hydrogen(NumForm::range_from(1), NumForm::Lit(0), true)]
    #[case::symbolic_count(NumForm::var("hydrogens"), NumForm::Lit(0), true)]
    #[case::absent(NumForm::Lit(0), NumForm::Lit(0), false)]
    #[case::both_present(NumForm::Lit(1), NumForm::Lit(1), false)]
    #[case::two_hydrogens(NumForm::Lit(2), NumForm::Undetermined, false)]
    #[case::two_lone_pairs(NumForm::Undetermined, NumForm::Lit(2), false)]
    #[case::excluded_counts(NumForm::range_from(2), NumForm::Undetermined, false)]
    #[case::bottom_count(NumForm::lit_set([]), NumForm::Undetermined, false)]
    #[case::symbolic_with_impossible_count(NumForm::var("hydrogens"), NumForm::Lit(2), false)]
    fn test_stereo_perceiver_derive_requested_frames(
        #[case] hydrogens: NumForm,
        #[case] lone_pairs: NumForm,
        #[case] compatible: bool,
        #[values(false, true)] stereo_bond: bool,
    ) {
        let center = AtomForm::from_element(Element::C)
            .with_implicit_hydrogens(hydrogens)
            .with_lone_pairs(lone_pairs);
        let mut entries = if stereo_bond {
            MoleculeEntries {
                atoms: vec![
                    center,
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(1)
                        .with_lone_pairs(0),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(2)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(1), AtomId(3), BondForm::from_order(1)),
                ],
                ..Default::default()
            }
        } else {
            MoleculeEntries {
                atoms: vec![
                    center,
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                ],
                bonds: (1..=3)
                    .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                    .collect(),
                ..Default::default()
            }
        };
        let failure = if stereo_bond {
            entries.bonds[0]
                .2
                .constraints
                .set(BondConstraintForm::cis_trans_stereo(true));
            StereoInconsistency::CisTransStereoFailure { bond: BondId(0) }
        } else {
            entries.atoms[0]
                .constraints
                .set(AtomConstraintForm::tetrahedral_stereo(true));
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) }
        };
        let molecule = Molecule::from_entries(entries);
        let original = molecule.clone();
        let expected = if compatible {
            Solution::Underdetermined(StereoDerivation::default())
        } else {
            Solution::Determined(StereoDerivation {
                inconsistencies: vec![failure],
                ..Default::default()
            })
        };
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            expected,
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::too_many_atom_ligands(
        mol_dsl!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "I" "N"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]}"#),
        StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) }
    )]
    #[case::too_few_atom_ligands(
        mol_dsl!(r#"{:atoms ["C#T+" "F"] :bonds [[0 1 "1"]]}"#),
        StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) }
    )]
    #[case::too_many_bond_ligands(
        mol_dsl!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"]]}"#),
        StereoInconsistency::CisTransStereoFailure { bond: BondId(0) }
    )]
    #[case::shared_bond_ligand(
        mol_dsl!(r#"{:atoms ["C" "C" "F"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [1 2 "1"]]}"#),
        StereoInconsistency::CisTransStereoFailure { bond: BondId(0) }
    )]
    fn test_stereo_perceiver_derive_requested_topology(
        #[case] molecule: Molecule,
        #[case] failure: StereoInconsistency,
    ) {
        let model = StereoModel {
            stereo_bond_min_ring_size: 0,
            ..Default::default()
        };
        let original = molecule.clone();
        assert_eq!(
            StereoPerceiver::new(&model).derive(&molecule, RingConfig::default()),
            Solution::Determined(StereoDerivation {
                inconsistencies: vec![failure],
                ..Default::default()
            }),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_stereo_perceiver_derive_requested_evidence() {
        let mut entries = MoleculeEntries {
            atoms: [
                (NumForm::Undetermined, NumForm::Undetermined),
                (NumForm::Lit(1), NumForm::Lit(0)),
                (NumForm::Lit(0), NumForm::Lit(0)),
            ]
            .into_iter()
            .flat_map(|(hydrogens, lone_pairs)| {
                [
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(hydrogens)
                        .with_lone_pairs(lone_pairs),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                ]
            })
            .collect(),
            bonds: (0..3)
                .flat_map(|group| {
                    (1..=3).map(move |offset| {
                        (
                            AtomId(4 * group),
                            AtomId(4 * group + offset),
                            BondForm::from_order(1),
                        )
                    })
                })
                .collect(),
            ..Default::default()
        };
        for atom in [0, 4, 8] {
            entries.atoms[atom]
                .constraints
                .set(AtomConstraintForm::tetrahedral_stereo(true));
        }
        let molecule = Molecule::from_entries(entries);
        let original = molecule.clone();
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            Solution::Underdetermined(StereoDerivation {
                stereo_atoms: vec![(
                    AtomId(4),
                    vec![
                        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                )],
                inconsistencies: vec![StereoInconsistency::TetrahedralStereoFailure {
                    atom: AtomId(8),
                }],
                ..Default::default()
            }),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_stereo_perceiver_derive_partial_evidence() {
        let mut entries = MoleculeEntries {
            atoms: [
                (NumForm::Undetermined, NumForm::Undetermined),
                (NumForm::Lit(1), NumForm::Lit(0)),
                (NumForm::Lit(0), NumForm::Lit(0)),
                (NumForm::Lit(1), NumForm::Lit(0)),
            ]
            .into_iter()
            .flat_map(|(hydrogens, lone_pairs)| {
                [
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(hydrogens)
                        .with_lone_pairs(lone_pairs),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                ]
            })
            .collect(),
            bonds: (0..4)
                .flat_map(|group| {
                    (1..=3).map(move |offset| {
                        (
                            AtomId(4 * group),
                            AtomId(4 * group + offset),
                            BondForm::from_order(1),
                        )
                    })
                })
                .collect(),
            stereo_atoms: (0..4)
                .map(|group| {
                    (
                        AtomId(4 * group),
                        vec![
                            StereoLigand::new(AtomId(4 * group), StereoLigandKind::Virtual),
                            StereoLigand::new(AtomId(4 * group + 3), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(4 * group + 2), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(4 * group + 1), StereoLigandKind::Atom),
                        ],
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                    )
                })
                .collect(),
            ..Default::default()
        };
        entries.atoms[12]
            .constraints
            .set(AtomConstraintForm::tetrahedral_stereo(
                TetrahedralStereoForm::NotStereo,
            ));
        let molecule = Molecule::from_entries(entries.clone());
        let mut stereo_atoms = entries.stereo_atoms[..2].to_vec();
        stereo_atoms[1].1[0].kind = StereoLigandKind::ImplicitHydrogen;
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive(&molecule, RingConfig::default()),
            Solution::Underdetermined(StereoDerivation {
                stereo_atoms,
                inconsistencies: vec![
                    StereoInconsistency::StereoAtomFailure {
                        stereo_atom: StereoAtomId(2)
                    },
                    StereoInconsistency::TetrahedralStereoMismatch {
                        atom: AtomId(12),
                        stereo_atom: StereoAtomId(3)
                    },
                ],
                ..Default::default()
            }),
        );
    }

    #[rstest]
    #[case::tetrahedral(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T+" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        Some((
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )),
    )]
    #[case::disabled_kind(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::Tetrahedral as usize] = None;
            model
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T+" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::element_scope(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::Tetrahedral as usize] = Some(StereoKindModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                fluxionality: false,
            });
            model
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h#T+" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h0#T+" "N#h2"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        None,
    )]
    #[case::aromatic_exclusion(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C#T+" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]
        }"#),
        None,
    )]
    fn test_stereo_perceiver_derive_stereo_atom(
        #[case] model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoAtomForm)>,
    ) {
        assert_eq!(
            StereoPerceiver::new(&model).derive_stereo_atom(
                &molecule,
                AtomId(1),
                &StereoCoset::Lit(1),
            ),
            expected
        );
    }

    #[rstest]
    #[case::actual(
        4, 0, 0,
        Some(vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ]),
    )]
    #[case::hydrogen(
        3, 1, 0,
        Some(vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        ]),
    )]
    #[case::lone_pair(
        3, 0, 1,
        Some(vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        ]),
    )]
    #[case::hydrogen_and_lone_pair(
        2, 1, 1,
        Some(vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        ]),
    )]
    #[case::extra_neighbor(5, 0, 0, None)]
    #[case::missing_neighbor(3, 0, 0, None)]
    #[case::extra_hydrogen(4, 1, 0, None)]
    #[case::extra_lone_pair(4, 0, 1, None)]
    #[case::extra_virtual_ligand(3, 1, 1, None)]
    fn test_stereo_perceiver_derive_stereo_atom_ligands(
        #[case] neighbors: u32,
        #[case] hydrogens: i64,
        #[case] lone_pairs: i64,
        #[case] expected: Option<Vec<StereoLigand>>,
        #[values(StereoCoset::Undetermined, StereoCoset::Lit(1))] coset: StereoCoset,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::P)
                    .with_implicit_hydrogens(hydrogens)
                    .with_lone_pairs(lone_pairs),
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: (1..=neighbors)
                .rev()
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive_stereo_atom(
                &molecule,
                AtomId(0),
                &coset
            ),
            expected.map(|ligands| (ligands, StereoAtomForm::new(StereoKind::Tetrahedral, coset))),
        );
    }

    #[rstest]
    #[case::hydrogens(true)]
    #[case::lone_pairs(false)]
    fn test_stereo_perceiver_derive_stereo_atom_count(
        #[case] hydrogen_count: bool,
        #[values(3u32, 4)] neighbors: u32,
        #[values(NumForm::Lit(-1), NumForm::Lit(2), NumForm::Undetermined, NumForm::RangeFrom(0))]
        count: NumForm,
    ) {
        let mut site = AtomForm::from_element(Element::P)
            .with_implicit_hydrogens(0)
            .with_lone_pairs(0);
        if hydrogen_count {
            site.implicit_hydrogens = count;
        } else {
            site.lone_pairs = count;
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                site,
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: (1..=neighbors)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive_stereo_atom(
                &molecule,
                AtomId(0),
                &StereoCoset::Undetermined
            ),
            None,
        );
    }

    #[rstest]
    #[case::cis_trans(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]
        }"#),
        Some((
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )),
    )]
    #[case::disabled_kind(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::CisTrans as usize] = None;
            model
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::element_scope(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::CisTrans as usize] = Some(StereoKindModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                fluxionality: false,
            });
            model
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "N#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h0#n0" "C#h" "C#h3"]
            :bonds [[0 1 "2#C+"] [1 2 "1"]]
        }"#),
        None,
    )]
    fn test_stereo_perceiver_derive_stereo_bond(
        #[case] model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoBondForm)>,
    ) {
        assert_eq!(
            StereoPerceiver::new(&model).derive_stereo_bond(
                &molecule,
                BondId(1),
                &StereoCoset::Lit(1),
            ),
            expected
        );
    }

    #[rstest]
    #[case::actual(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "H" "F" "Cl" "Br"]
            :bonds [[0 1 "2"] [0 3 "1"] [1 5 "1"] [0 2 "1"] [1 4 "1"]]
        }"#),
        Some(vec![
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
        ]),
    )]
    #[case::lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "N#n" "H" "F"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        }"#),
        Some(vec![
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        ]),
    )]
    #[case::hydrogen_and_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h#n" "N#h#n"]
            :bonds [[0 1 "2"]]
        }"#),
        Some(vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        ]),
    )]
    #[case::extra_neighbor(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "H" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"]]
        }"#),
        None,
    )]
    #[case::shared_neighbor(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "F"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 2 "1"]]
        }"#),
        None,
    )]
    #[case::extra_virtual_ligand(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h#n" "N#h#n" "H" "F"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    fn test_stereo_perceiver_derive_stereo_bond_ligands(
        #[case] molecule: Molecule,
        #[case] expected: Option<Vec<StereoLigand>>,
        #[values(StereoCoset::Undetermined, StereoCoset::Lit(1))] coset: StereoCoset,
    ) {
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive_stereo_bond(
                &molecule,
                BondId(0),
                &coset
            ),
            expected.map(|ligands| (ligands, StereoBondForm::new(StereoKind::CisTrans, coset))),
        );
    }

    #[rstest]
    #[case::hydrogens(true)]
    #[case::lone_pairs(false)]
    fn test_stereo_perceiver_derive_stereo_bond_count(
        #[case] hydrogen_count: bool,
        #[values(AtomId(0), AtomId(1))] endpoint: AtomId,
        #[values(1u32, 2)] neighbors: u32,
        #[values(NumForm::Lit(-1), NumForm::Lit(2), NumForm::Undetermined, NumForm::RangeFrom(0))]
        count: NumForm,
    ) {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        }"#
        );
        if neighbors == 2 {
            let mut editor = molecule.edit();
            let ligand = editor.add_atom(AtomForm::from_element(Element::Br).into_concrete());
            editor.add_bond(endpoint, ligand, BondForm::from_order(1).into_concrete());
            molecule = editor.finish().unwrap();
        }
        let mut atom = molecule.atom_mut(endpoint);
        if hydrogen_count {
            atom.attributes_mut().implicit_hydrogens = count;
        } else {
            atom.attributes_mut().lone_pairs = count;
        }
        assert_eq!(
            StereoPerceiver::new(&StereoModel::default()).derive_stereo_bond(
                &molecule,
                BondId(0),
                &StereoCoset::Undetermined
            ),
            None,
        );
    }
}
