//! Aromatize: Kekulé form → aromatic-system form via perception.
//!
//! Wraps [`AromaticityResolver`] in the [`Transformer`] interface. If the
//! input AST already carries one or more aromatic systems, this is a no-op:
//! to re-aromatize, the user kekulizes first and then aromatizes.

use umol_ast::ast::MoleculeAst;

use crate::ops::aromaticity::AromaticityResolver;
use crate::ops::config::AromaticityModel;
use crate::ops::solution::Solution;
use crate::ops::transform::{Transformer, TransformerError};

#[derive(Clone, Debug)]
pub struct Aromatize {
    model: AromaticityModel,
}

impl Aromatize {
    pub fn new(model: AromaticityModel) -> Self {
        Self { model }
    }
}

impl Transformer for Aromatize {
    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), TransformerError> {
        if ast.aromatic_systems().count() > 0 {
            return Ok(());
        }
        let resolver = AromaticityResolver::new(&self.model);
        match resolver.resolve(ast) {
            Ok(Solution::Determined(())) => Ok(()),
            Ok(Solution::Underdetermined(())) => Err(TransformerError::AromatizeUnderdetermined),
            Ok(Solution::Contradictory(c)) => Err(TransformerError::AromatizeContradiction(c)),
            Err(e) => Err(TransformerError::AromatizeSetup(e)),
        }
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a> {
        Box::new(self.transform(ast).ok().into_iter())
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticSystemIdx, AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst,
        BondConstraintKind, Constraints, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;

    fn aromatic_carbon() -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.spin = SpinStateAst::closed_shell();
        atom.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        ));
        atom
    }

    fn benzene_kekule() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| aromatic_carbon()).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| {
                let order = if i % 2 == 0 { 2 } else { 1 };
                (
                    AtomIdx(i),
                    AtomIdx((i + 1) % 6),
                    BondAst::from_order(order),
                )
            })
            .collect();
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_aromatize_kekule_benzene_adds_aromatic_system() {
        let mut ast = benzene_kekule();
        Aromatize::new(AromaticityModel::daylight())
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 1);
        let view = ast.aromatic_system(AromaticSystemIdx(0));
        let atoms: Vec<AtomIdx> = view.atoms().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| {
                view.data
                    .constraints
                    .iter()
                    .any(|c| c.kind() == BondConstraintKind::Aromatic)
            })
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromatize_already_aromatic_is_noop() {
        let original = {
            let mut ast = benzene_kekule();
            Aromatize::new(AromaticityModel::daylight())
                .transform_into(&mut ast)
                .unwrap();
            ast
        };
        let mut second = original.clone();
        Aromatize::new(AromaticityModel::daylight())
            .transform_into(&mut second)
            .unwrap();
        assert_eq!(original, second);
    }

    #[rstest]
    fn test_aromatize_transform_returns_new_ast() {
        let ast = benzene_kekule();
        let aromatized = Aromatize::new(AromaticityModel::daylight())
            .transform(&ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 0);
        assert_eq!(aromatized.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatize_generate_all_yields_one() {
        let ast = benzene_kekule();
        let transformer = Aromatize::new(AromaticityModel::daylight());
        let results: Vec<MoleculeAst> = transformer.generate_all(&ast).collect();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].aromatic_systems().count(), 1);
    }
}
