//! Ground-invariant molecule with cached derived views.

use std::sync::{Arc, OnceLock};

use crate::ast::error::GroundError;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{RingEnumerationStrategy, RingEnumerator, RingFamily, RingSet};

#[derive(Debug)]
struct MoleculeInner {
    ast: MoleculeAst,
    rings: OnceLock<RingSet>,
}

#[derive(Clone, Debug)]
pub struct Molecule(Arc<MoleculeInner>);

impl Molecule {
    pub fn new(ast: MoleculeAst) -> Result<Self, GroundError> {
        if !ast.is_ground() {
            return Err(GroundError);
        }
        Ok(Self(Arc::new(MoleculeInner {
            ast,
            rings: OnceLock::new(),
        })))
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.0.ast
    }

    pub fn rings(&self) -> &RingSet {
        self.0.rings.get_or_init(|| {
            RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
                .enumerate(&self.0.ast)
        })
    }
}

#[cfg(test)]
mod tests {
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::atom::AtomAst;

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
            valence: ValueAst::Lit(4),
            donated_pairs: ValueAst::Lit(0),
            accepted_pairs: ValueAst::Lit(0),
            aromatic_valence: AromaticValenceAst::NotAromatic,
            multicenter_valence: ValueAst::Lit(0),
        }
    }

    #[test]
    fn test_molecule_new() {
        let ast = MoleculeAst::new(
            vec![ground_atom()],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert!(Molecule::new(ast).is_ok());
    }

    #[test]
    fn test_molecule_new_error() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert!(matches!(Molecule::new(ast), Err(GroundError)));
    }

    #[test]
    fn test_molecule_rings_caches() {
        let ast = MoleculeAst::new(
            vec![ground_atom()],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let mol = Molecule::new(ast).unwrap();
        let r1 = mol.rings() as *const RingSet;
        let r2 = mol.rings() as *const RingSet;
        assert_eq!(r1, r2);
    }
}
