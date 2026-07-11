//! Aromatize: Kekulé form → aromatic-system form.
//!
//! [`Aromatizer`] runs aromaticity perception against a Kekulé-form input —
//! atoms with explicit single/double bonds and no aromatic hints. Per-atom π
//! contributions are derived from bond orders by [`electrons_from_kekule`]
//! rather than from the `AromaticValence` constraint that the resolver reads.
//! If the input AST already carries one or more aromatic systems, this is a
//! no-op: re-aromatizing requires kekulizing first.

use thiserror::Error;
use umol_ast::ast::{AtomView, ElementAst, MoleculeAst, ValueAst};
use umol_chem::element::Element;

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError, AromaticityPerception};
use crate::ops::model::AromaticityModel;
use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromatizerError {
    #[error("aromaticity setup: {0}")]
    Setup(#[from] AromaticityError),
    #[error("aromaticity contradiction: {0}")]
    Contradiction(#[from] AromaticityContradiction),
    #[error("aromatization input is underdetermined")]
    Underdetermined,
}

#[derive(Clone, Debug)]
pub struct Aromatizer {
    perception: AromaticityPerception,
}

impl Aromatizer {
    pub fn new(model: &AromaticityModel) -> Self {
        Self {
            perception: AromaticityPerception::new(model),
        }
    }
}

impl Transformer for Aromatizer {
    type Error = AromatizerError;

    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), AromatizerError> {
        if ast.aromatic_systems().count() > 0 {
            return Ok(());
        }
        let systems = self
            .perception
            .find_systems(ast, electrons_from_kekule)?
            .into_decisive(AromatizerError::Underdetermined)?;
        self.perception.add_systems(ast, systems);
        Ok(())
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a> {
        Box::new(self.transform(ast).ok().into_iter())
    }
}

/// Derive an atom's π contribution from a Kekulé bond-order layout.
///
/// - Exactly one incident double bond → 1 π electron (sp² atom on a single
///   π bond, e.g. benzene C, pyridine N).
/// - Zero incident double bonds, atom is an N/O/S/Se/P/As → 2 π electrons
///   (pyrrole-, furan-, thiophene-class heteroatom donating a lone pair).
/// - Zero incident double bonds, atom is C with charge `+1` → 0 π electrons
///   (sp² carbocation, empty p_z, e.g. tropylium C⁺).
/// - Anything else (sp³ C, two or more double bonds, undetermined data) →
///   `None`, marking the atom as not aromatic-eligible.
pub fn electrons_from_kekule(view: &AtomView<'_>) -> Option<u8> {
    let ElementAst::Lit(element) = view.ast.element else {
        return None;
    };
    let double_count = view
        .neighbors()
        .filter(|n| matches!(n.bond().ast.order, ValueAst::Lit(2)))
        .count();
    match double_count {
        1 => Some(1),
        0 => match element {
            Element::N | Element::O | Element::S | Element::Se | Element::P | Element::As => {
                Some(2)
            }
            Element::C if matches!(view.ast.charge, ValueAst::Lit(1)) => Some(0),
            _ => None,
        },
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticSystemId, AtomAst, AtomId, BondAst, BondConstraintKey, MoleculeAst, MoleculeParts,
        SpinStateAst,
    };
    use umol_chem::element::Element;

    use super::*;

    fn kekule_carbon() -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.spin = SpinStateAst::closed_shell();
        atom
    }

    fn benzene_kekule() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| kekule_carbon()).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| {
                let order = if i % 2 == 0 { 2 } else { 1 };
                (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(order))
            })
            .collect();
        MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    #[rstest]
    fn test_aromatizer_kekule_benzene_adds_aromatic_system() {
        let mut ast = benzene_kekule();
        Aromatizer::new(&AromaticityModel::daylight())
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 1);
        let view = ast.aromatic_system(AromaticSystemId(0));
        let atoms: Vec<AtomId> = view.atom_ids().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| view.ast.constraints.contains(BondConstraintKey::Aromatic))
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromatizer_already_aromatic_is_noop() {
        let original = {
            let mut ast = benzene_kekule();
            Aromatizer::new(&AromaticityModel::daylight())
                .transform_into(&mut ast)
                .unwrap();
            ast
        };
        let mut second = original.clone();
        Aromatizer::new(&AromaticityModel::daylight())
            .transform_into(&mut second)
            .unwrap();
        assert_eq!(original, second);
    }

    #[rstest]
    fn test_aromatizer_transform_returns_new_ast() {
        let ast = benzene_kekule();
        let aromatized = Aromatizer::new(&AromaticityModel::daylight())
            .transform(&ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 0);
        assert_eq!(aromatized.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatizer_generate_all_yields_one() {
        let ast = benzene_kekule();
        let transformer = Aromatizer::new(&AromaticityModel::daylight());
        let results: Vec<MoleculeAst> = transformer.generate_all(&ast).collect();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].aromatic_systems().count(), 1);
    }
}
