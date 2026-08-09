//! Connectivity (tier-3) validator — checks molecule connectivity with respect to [`ConnectivityModel`].

use std::collections::HashSet;
use std::iter;

use thiserror::Error;
use umol_graph_core::UnionFind;
use umol_utils::solution::Solution;

use super::super::id::{
    AromaticSystemId, AtomId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
    StereoBondId,
};
use super::super::molecule::MoleculeAst;
#[cfg(test)]
use super::super::molecule::MoleculeEntries;

/// Connectivity definitions.
/// - `allow_disconnected`: allow disconnected atom / bond graph
/// - `allow_disconnected_<family>`: allow straddling relations of that family
///   (`false` = its atoms must share one bond component).
///
/// The defaults permit a disconnected molecule and straddling dative / multicenter / noncovalent bonds
/// and molecule constraints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConnectivityModel {
    pub allow_disconnected: bool,
    pub allow_disconnected_dative: bool,
    pub allow_disconnected_aromatic: bool,
    pub allow_disconnected_multicenter: bool,
    pub allow_disconnected_noncovalent: bool,
    pub allow_disconnected_stereo_atom: bool,
    pub allow_disconnected_stereo_bond: bool,
    pub allow_disconnected_constraints: bool,
}

impl Default for ConnectivityModel {
    fn default() -> Self {
        Self {
            allow_disconnected: true,
            allow_disconnected_dative: true,
            allow_disconnected_aromatic: false,
            allow_disconnected_multicenter: true,
            allow_disconnected_noncovalent: true,
            allow_disconnected_stereo_atom: false,
            allow_disconnected_stereo_bond: false,
            allow_disconnected_constraints: true,
        }
    }
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConnectivityContradiction {
    #[error("molecule is disconnected: {components} bond components")]
    Disconnected { components: usize },
    #[error("dative bond {bond:?} straddles two bond components")]
    DisconnectedDativeBond { bond: DativeBondId },
    #[error("aromatic system {system:?} straddles two bond components")]
    DisconnectedAromaticSystem { system: AromaticSystemId },
    #[error("multicenter bond {bond:?} straddles two bond components")]
    DisconnectedMulticenterBond { bond: MulticenterBondId },
    #[error("noncovalent bond {bond:?} straddles two bond components")]
    DisconnectedNoncovalentBond { bond: NoncovalentBondId },
    #[error("stereo atom {stereo_atom:?} straddles two bond components")]
    DisconnectedStereoAtom { stereo_atom: StereoAtomId },
    #[error("stereo bond {stereo_bond:?} straddles two bond components")]
    DisconnectedStereoBond { stereo_bond: StereoBondId },
    #[error("molecule constraint {index} references atoms in two bond components")]
    DisconnectedConstraint { index: usize },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConnectivityError {}

pub struct ConnectivityValidator<'a> {
    model: &'a ConnectivityModel,
}

impl<'a> ConnectivityValidator<'a> {
    pub fn new(model: &'a ConnectivityModel) -> Self {
        Self { model }
    }

    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ConnectivityContradiction>, ConnectivityError> {
        let atom_count = ast.atoms().count();
        let mut union = UnionFind::new(atom_count);
        for bond in ast.bonds().iter() {
            let [a, b] = bond.atom_ids();
            union.union(a.index(), b.index());
        }
        let roots: Vec<usize> = (0..atom_count).map(|i| union.find(i)).collect();

        if !self.model.allow_disconnected {
            let components = roots.iter().collect::<HashSet<_>>().len();
            if components > 1 {
                return contradiction(ConnectivityContradiction::Disconnected { components });
            }
        }
        if !self.model.allow_disconnected_dative {
            for (index, dative) in ast.dative_bonds().iter().enumerate() {
                if spans(&roots, dative.atom_ids()) {
                    return contradiction(ConnectivityContradiction::DisconnectedDativeBond {
                        bond: DativeBondId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_aromatic {
            for (index, system) in ast.aromatic_systems().iter().enumerate() {
                if spans(&roots, system.atom_ids()) {
                    return contradiction(ConnectivityContradiction::DisconnectedAromaticSystem {
                        system: AromaticSystemId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_multicenter {
            for (index, bond) in ast.multicenter_bonds().iter().enumerate() {
                if spans(&roots, bond.atom_ids()) {
                    return contradiction(ConnectivityContradiction::DisconnectedMulticenterBond {
                        bond: MulticenterBondId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_noncovalent {
            for (index, bond) in ast.noncovalent_bonds().iter().enumerate() {
                if spans(&roots, bond.atom_ids()) {
                    return contradiction(ConnectivityContradiction::DisconnectedNoncovalentBond {
                        bond: NoncovalentBondId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_stereo_atom {
            for (index, stereo) in ast.stereo_atoms().iter().enumerate() {
                let atoms =
                    iter::once(stereo.site_id()).chain(stereo.ligands().map(|l| l.atom_id()));
                if spans(&roots, atoms) {
                    return contradiction(ConnectivityContradiction::DisconnectedStereoAtom {
                        stereo_atom: StereoAtomId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_stereo_bond {
            for (index, stereo) in ast.stereo_bonds().iter().enumerate() {
                let [a, b] = ast.bond(stereo.site_id()).atom_ids();
                let atoms = [a, b]
                    .into_iter()
                    .chain(stereo.ligands().map(|l| l.atom_id()));
                if spans(&roots, atoms) {
                    return contradiction(ConnectivityContradiction::DisconnectedStereoBond {
                        stereo_bond: StereoBondId(index as u32),
                    });
                }
            }
        }
        if !self.model.allow_disconnected_constraints {
            for (index, constraint) in ast.constraints().iter().enumerate() {
                if spans(&roots, ast.constraint_atoms(constraint)) {
                    return contradiction(ConnectivityContradiction::DisconnectedConstraint {
                        index,
                    });
                }
            }
        }
        Ok(Solution::Determined(()))
    }
}

/// Whether the atoms fall in more than one bond component.
fn spans(roots: &[usize], atoms: impl IntoIterator<Item = AtomId>) -> bool {
    let mut atoms = atoms.into_iter();
    match atoms.next() {
        Some(first) => {
            let component = roots[first.index()];
            atoms.any(|atom| roots[atom.index()] != component)
        }
        None => false,
    }
}

fn contradiction(
    contradiction: ConnectivityContradiction,
) -> Result<Solution<(), ConnectivityContradiction>, ConnectivityError> {
    Ok(Solution::Contradictory(contradiction))
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::Constraints;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};

    #[rstest]
    fn test_connectivity_model_default() {
        assert_eq!(
            ConnectivityModel::default(),
            ConnectivityModel {
                allow_disconnected: true,
                allow_disconnected_dative: true,
                allow_disconnected_aromatic: false,
                allow_disconnected_multicenter: true,
                allow_disconnected_noncovalent: true,
                allow_disconnected_stereo_atom: false,
                allow_disconnected_stereo_bond: false,
                allow_disconnected_constraints: true,
            }
        );
    }

    #[rstest]
    fn test_connectivity_validator_validate_disconnected_allowed() {
        let mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        });
        let model = ConnectivityModel::default();

        assert_eq!(
            ConnectivityValidator::new(&model).validate(&mol),
            Ok(Solution::Determined(()))
        );
    }

    #[rstest]
    fn test_connectivity_validator_validate_disconnected_forbidden() {
        let mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        });
        let model = ConnectivityModel {
            allow_disconnected: false,
            ..ConnectivityModel::default()
        };

        assert_eq!(
            ConnectivityValidator::new(&model).validate(&mol),
            Ok(Solution::Contradictory(
                ConnectivityContradiction::Disconnected { components: 2 }
            ))
        );
    }

    #[rstest]
    fn test_connectivity_validator_validate_aromatic_spanning() {
        // an aromatic system over atoms in the two separate bond components — disallowed by default
        let mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let model = ConnectivityModel::default();

        assert_eq!(
            ConnectivityValidator::new(&model).validate(&mol),
            Ok(Solution::Contradictory(
                ConnectivityContradiction::DisconnectedAromaticSystem {
                    system: AromaticSystemId(0)
                }
            ))
        );
    }

    #[rstest]
    fn test_connectivity_validator_validate_noncovalent_spanning_allowed() {
        // a noncovalent bond bridging the two components — permitted by default
        let mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            noncovalent: vec![(
                AtomId(0),
                AtomId(2),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let model = ConnectivityModel::default();

        assert_eq!(
            ConnectivityValidator::new(&model).validate(&mol),
            Ok(Solution::Determined(()))
        );
    }
}
