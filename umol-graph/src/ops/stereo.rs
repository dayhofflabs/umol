//! Stereo sites, selection, and structural perception.
//!
//! [`StereoPerception`] realizes `#T`/`#C` evidence, recognizes missing
//! tetrahedral and cis/trans entities, and classifies existing arrangements.
//! Resolver inconsistency policy and edit construction remain separate.

mod recognize;

use std::collections::BTreeSet;

pub use recognize::StereoPerceptionError;
use smallvec::SmallVec;
use thiserror::Error;
use umol_graph_ir::ir::{
    AsLit, AtomConstraintForm, AtomConstraintKey, AtomId, BondConstraintForm, BondConstraintKey,
    BondId, CisTransStereoForm, Molecule, RingConfig, RingModel, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    Stereogenicity, TetrahedralStereoForm,
};

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
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation {
    /// Realizable tetrahedral stereo atoms.
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>,
    /// Realizable cis-trans stereo bonds.
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>,
    /// Sites of cis-trans assertions on double bonds whose smallest ring is below the model's
    /// `stereo_bond_minimum_ring_size`: no stereo bond is realized, and nothing is inconsistent.
    pub skipped_stereo_bonds: Vec<BondId>,
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

/// Stereo evidence, recognition, and classification under a selected chemistry model.
#[derive(Clone, Debug)]
pub struct StereoPerception {
    model: StereoModel,
}

impl StereoPerception {
    pub fn new(model: &StereoModel) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Sites of cis-trans assertions the model does not realize: double bonds whose smallest
    /// ring is below `stereo_bond_minimum_ring_size`. Empty when the threshold is zero or no bond
    /// carries an assertion, without enumerating rings.
    pub fn skipped_stereo_bonds(&self, molecule: &Molecule) -> BTreeSet<BondId> {
        let minimum_ring_size = self.model.stereo_bond_minimum_ring_size as usize;
        let asserted: Vec<BondId> = molecule
            .bonds()
            .iter()
            .filter(|bond| {
                matches!(
                    bond.constraints().cis_trans_stereo(),
                    CisTransStereoForm::Stereo(_)
                )
            })
            .map(|bond| bond.id())
            .collect();
        if minimum_ring_size == 0 || asserted.is_empty() {
            return BTreeSet::new();
        }
        let rings = molecule.rings(RingModel::default(), RingConfig::default());
        asserted
            .into_iter()
            .filter(|&bond| {
                rings
                    .bond(bond)
                    .smallest_ring_size()
                    .is_some_and(|size| size < minimum_ring_size)
            })
            .collect()
    }

    /// Derive stereo relations from complete evidence and check existing relations.
    /// Positive evidence supplies the candidate configuration; negative evidence supplies none.
    /// Contradictory evidence between realizable sources is a mismatch. Independent constraint
    /// and entity realizability failures retain their separate diagnostics.
    pub fn derive(&self, molecule: &Molecule) -> StereoDerivation {
        let mut stereo_atoms = Vec::new();
        let mut stereo_bonds = Vec::new();
        let mut inconsistencies = BTreeSet::new();
        let skipped = self.skipped_stereo_bonds(molecule);

        for view in molecule.atoms().iter() {
            let atom = view.id();
            let relation = molecule.stereo_atoms().at(atom);
            let constraints = view.constraints();
            let asserted = constraints.tetrahedral_stereo();
            if relation.is_none() && !matches!(asserted, TetrahedralStereoForm::Stereo(_)) {
                continue;
            }
            let evidence = constraints.evidenced_complete(AtomConstraintKey::TetrahedralStereo);
            let coset = match &evidence {
                Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                    coset,
                ))) => Some(coset),
                _ => None,
            };
            let candidate = if coset.is_some() || relation.is_some() {
                self.derive_stereo_atom(molecule, atom, coset.unwrap_or(&StereoCoset::Undetermined))
            } else {
                None
            };
            if candidate.is_none()
                && (matches!(asserted, TetrahedralStereoForm::Stereo(_))
                    || evidence.is_none() && relation.is_none())
            {
                inconsistencies.insert(StereoInconsistency::TetrahedralStereoFailure { atom });
            }
            if let Some(relation) = relation {
                let realizable = relation.attributes().configuration.kind()
                    == Some(StereoKind::Tetrahedral)
                    && candidate.as_ref().is_some_and(|(ligands, _)| {
                        // Evidence already transports the configuration; tetrahedral permits
                        // every permutation of the same distinct ligands.
                        (evidence.is_some()
                            && relation.ligand_ids().len() == ligands.len()
                            && relation
                                .ligand_ids()
                                .iter()
                                .all(|ligand| ligands.contains(ligand)))
                            || relation.coset_for(ligands.iter().copied()).is_some()
                    });
                if !realizable {
                    inconsistencies.insert(StereoInconsistency::StereoAtomFailure {
                        stereo_atom: relation.id(),
                    });
                    continue;
                } else if evidence.is_none() {
                    inconsistencies.insert(StereoInconsistency::TetrahedralStereoMismatch {
                        atom,
                        stereo_atom: relation.id(),
                    });
                }
            }
            if coset.is_some() {
                if let Some((ligands, stereo)) = candidate {
                    stereo_atoms.push((atom, ligands, stereo));
                }
            }
        }

        for view in molecule.bonds().iter() {
            let bond = view.id();
            let relation = molecule.stereo_bonds().at(bond);
            let constraints = view.constraints();
            let asserted = constraints.cis_trans_stereo();
            if relation.is_none() && !matches!(asserted, CisTransStereoForm::Stereo(_)) {
                continue;
            }
            let evidence = constraints.evidenced_complete(BondConstraintKey::CisTransStereo);
            let coset = match &evidence {
                Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) => {
                    Some(coset)
                }
                _ => None,
            };
            let candidate = if coset.is_some() && !skipped.contains(&bond) || relation.is_some() {
                self.derive_stereo_bond(molecule, bond, coset.unwrap_or(&StereoCoset::Undetermined))
            } else {
                None
            };
            if candidate.is_none()
                && !skipped.contains(&bond)
                && (matches!(asserted, CisTransStereoForm::Stereo(_))
                    || evidence.is_none() && relation.is_none())
            {
                inconsistencies.insert(StereoInconsistency::CisTransStereoFailure { bond });
            }
            if let Some(relation) = relation {
                let realizable = relation.attributes().configuration.kind()
                    == Some(StereoKind::CisTrans)
                    && candidate.as_ref().is_some_and(|(ligands, _)| {
                        (evidence.is_some() && relation.ligand_ids() == ligands.as_slice())
                            || relation.coset_for(ligands.iter().copied()).is_some()
                    });
                if !realizable {
                    inconsistencies.insert(StereoInconsistency::StereoBondFailure {
                        stereo_bond: relation.id(),
                    });
                    continue;
                } else if evidence.is_none() {
                    inconsistencies.insert(StereoInconsistency::CisTransStereoMismatch {
                        bond,
                        stereo_bond: relation.id(),
                    });
                }
            }
            if coset.is_some() && !skipped.contains(&bond) {
                if let Some((ligands, stereo)) = candidate {
                    stereo_bonds.push((bond, ligands, stereo));
                }
            }
        }

        StereoDerivation {
            stereo_atoms,
            stereo_bonds,
            skipped_stereo_bonds: skipped.into_iter().collect(),
            inconsistencies: inconsistencies.into_iter().collect(),
        }
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
    use umol_graph_ir::ir::{
        AtomForm, BondForm, MoleculeEntries, NumForm, StereoAtomId, StereoBondId,
    };
    use umol_graph_ir::mol_dsl_concrete;

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
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C1"] [6 7 "1"]]
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
            skipped_stereo_bonds: vec![],
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
            skipped_stereo_bonds: vec![],
            inconsistencies: vec![],
        },
    )]
    #[case::unrealizable_sites(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h0#T1" "C#h3"
                    "C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [3 4 "1"] [4 5 "2#C1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![],
            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                StereoInconsistency::CisTransStereoFailure { bond: BondId(3) },
            ],
        },
    )]
    #[case::relation_mismatches(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C1"] [6 7 "1"]]
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
            skipped_stereo_bonds: vec![],
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
            skipped_stereo_bonds: vec![],
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
            skipped_stereo_bonds: vec![],
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
            :bonds [[0 5 "1"] [0 1 "2#C1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![BondId(1)],
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
            skipped_stereo_bonds: vec![BondId(1)],
            inconsistencies: vec![],
        },
    )]
    fn test_stereo_perception_derive(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive(&molecule),
            expected
        );
    }

    #[rstest]
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
    fn test_stereo_perception_derive_evidence(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive(&molecule),
            expected,
        );
    }

    #[rstest]
    #[case::tetrahedral(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
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
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
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
            :atoms ["C#h3" "S#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h0#T1" "N#h2"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        None,
    )]
    #[case::aromatic_exclusion(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C#T1" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]
        }"#),
        None,
    )]
    fn test_stereo_perception_derive_stereo_atom(
        #[case] model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoAtomForm)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model).derive_stereo_atom(
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
    fn test_stereo_perception_derive_stereo_atom_ligands(
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
            StereoPerception::new(&StereoModel::default()).derive_stereo_atom(
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
    fn test_stereo_perception_derive_stereo_atom_count(
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
            StereoPerception::new(&StereoModel::default()).derive_stereo_atom(
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
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
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
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
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
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h0#n0" "C#h" "C#h3"]
            :bonds [[0 1 "2#C1"] [1 2 "1"]]
        }"#),
        None,
    )]
    fn test_stereo_perception_derive_stereo_bond(
        #[case] model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoBondForm)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model).derive_stereo_bond(
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
    fn test_stereo_perception_derive_stereo_bond_ligands(
        #[case] molecule: Molecule,
        #[case] expected: Option<Vec<StereoLigand>>,
        #[values(StereoCoset::Undetermined, StereoCoset::Lit(1))] coset: StereoCoset,
    ) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive_stereo_bond(
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
    fn test_stereo_perception_derive_stereo_bond_count(
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
            StereoPerception::new(&StereoModel::default()).derive_stereo_bond(
                &molecule,
                BondId(0),
                &StereoCoset::Undetermined
            ),
            None,
        );
    }
}
