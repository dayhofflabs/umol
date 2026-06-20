//! Structural shape checks on per-relation entities: when the `electrons:
//! ElectronCountsAst` field on an aromatic system or multicenter bond is
//! `Lit`, its length must match the participants list.

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemView, ElectronCountsAst, MoleculeAst, MulticenterBondView,
};

use umol_shared::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct EntityStructureValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureContradiction {
    #[error("aromatic system: electrons.len() = {electrons_len} but atoms.len() = {atoms_len}")]
    AromaticSystemElectronsLengthMismatch {
        electrons_len: usize,
        atoms_len: usize,
    },
    #[error("multicenter bond: electrons.len() = {electrons_len} but atoms.len() = {atoms_len}")]
    MulticenterElectronsLengthMismatch {
        electrons_len: usize,
        atoms_len: usize,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureError {}

impl EntityStructureValidator {
    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), EntityStructureContradiction>, EntityStructureError> {
        let ast = ast.as_ref();
        for view in ast.aromatic_systems().iter() {
            if let Some(c) = aromatic_system_length_check(&view) {
                return Ok(Solution::Contradictory(c));
            }
        }
        for view in ast.multicenter_bonds().iter() {
            if let Some(c) = multicenter_length_check(&view) {
                return Ok(Solution::Contradictory(c));
            }
        }
        Ok(Solution::Determined(()))
    }
}

fn aromatic_system_length_check(
    view: &AromaticSystemView<'_>,
) -> Option<EntityStructureContradiction> {
    let ElectronCountsAst::Lit(counts) = &view.ast.electrons else {
        return None;
    };
    let atoms_len = view.atoms().count();
    let electrons_len = counts.len();
    if electrons_len != atoms_len {
        Some(
            EntityStructureContradiction::AromaticSystemElectronsLengthMismatch {
                electrons_len,
                atoms_len,
            },
        )
    } else {
        None
    }
}

fn multicenter_length_check(
    view: &MulticenterBondView<'_>,
) -> Option<EntityStructureContradiction> {
    let ElectronCountsAst::Lit(counts) = &view.ast.electrons else {
        return None;
    };
    let atoms_len = view.atoms().count();
    let electrons_len = counts.len();
    if electrons_len != atoms_len {
        Some(
            EntityStructureContradiction::MulticenterElectronsLengthMismatch {
                electrons_len,
                atoms_len,
            },
        )
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemAst, AtomAst, AtomId, Constraints, MoleculeAst, MulticenterBondAst,
    };
    use umol_shared::element::Element;

    use super::*;

    #[rstest]
    fn test_entity_structure_validator_aromatic_length_mismatch() {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let aromatic = vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::from_counts(vec![1, 1]),
        )];
        let ast = MoleculeAst::from_parts(
            atoms,
            vec![],
            vec![],
            aromatic,
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        let result = v.validate(ast).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(
                EntityStructureContradiction::AromaticSystemElectronsLengthMismatch {
                    electrons_len: 2,
                    atoms_len: 3,
                }
            )
        ));
    }

    #[rstest]
    fn test_entity_structure_validator_multicenter_length_mismatch() {
        let atoms = vec![
            AtomAst::from_element(Element::B),
            AtomAst::from_element(Element::B),
            AtomAst::from_element(Element::H),
        ];
        let multicenter = vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondAst::from_counts(vec![1]),
        )];
        let ast = MoleculeAst::from_parts(
            atoms,
            vec![],
            vec![],
            vec![],
            multicenter,
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        let result = v.validate(ast).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(
                EntityStructureContradiction::MulticenterElectronsLengthMismatch {
                    electrons_len: 1,
                    atoms_len: 3,
                }
            )
        ));
    }

    #[rstest]
    fn test_entity_structure_validator_empty_electrons_passes() {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let aromatic = vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default())];
        let ast = MoleculeAst::from_parts(
            atoms,
            vec![],
            vec![],
            aromatic,
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        assert!(matches!(v.validate(ast).unwrap(), Solution::Determined(())));
    }
}
