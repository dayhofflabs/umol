//! Ground-invariant molecule with cached derived views.

use std::ops::Index;
use std::sync::Arc;

use umol_edn::ToEdn;
use umol_graph_core::Graph;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::molecule::{AromaticSystemAst, MoleculeAst, MulticenterBondAst};
use crate::ast::rings::{RingCache, RingEnumerationStrategy, RingFamily, RingSet};
use crate::ast::views::{
    AromaticSystemViews, AtomView, AtomViews, BondView, BondViews, DativeBondViews,
    MulticenterBondViews, NeighborView, NoncovalentBondViews,
};
use crate::ast::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::dsl::molecule::{parse_molecule_dsl, MoleculeAstWrapper};

use super::error::MoleculeEdnError;
use crate::ast::error::GroundError;

#[derive(Debug)]
struct MoleculeInner {
    ast: MoleculeAst,
    rings: RingCache,
}

#[derive(Clone, Debug)]
pub struct Molecule(Arc<MoleculeInner>);

impl Molecule {
    pub fn new(ast: MoleculeAst) -> Result<Self, GroundError> {
        Self::from_parts(ast, RingCache::new())
    }

    pub(crate) fn from_parts(
        ast: MoleculeAst,
        rings: RingCache,
    ) -> Result<Self, GroundError> {
        if !ast.is_ground() {
            return Err(GroundError);
        }
        Ok(Self(Arc::new(MoleculeInner { ast, rings })))
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.0.ast
    }

    pub fn into_ast(self) -> MoleculeAst {
        match Arc::try_unwrap(self.0) {
            Ok(inner) => inner.ast,
            Err(arc) => arc.ast.clone(),
        }
    }

    pub fn rings(&self) -> Arc<RingSet> {
        self.rings_with(RingFamily::Simple, &RingEnumerationStrategy::default())
    }

    pub fn rings_with(
        &self,
        family: RingFamily,
        strategy: &RingEnumerationStrategy,
    ) -> Arc<RingSet> {
        self.0.rings.get(&self.0.ast, family, strategy)
    }

    pub fn atoms(&self) -> AtomViews<'_> {
        self.0.ast.atoms()
    }

    pub fn bonds(&self) -> BondViews<'_> {
        self.0.ast.bonds()
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        self.0.ast.dative_bonds()
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        self.0.ast.noncovalent_bonds()
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        self.0.ast.aromatic_systems()
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        self.0.ast.multicenter_bonds()
    }

    pub fn atom(&self, idx: AtomIdx) -> AtomView<'_> {
        self.0.ast.atom(idx)
    }

    pub fn bond(&self, idx: BondIdx) -> BondView<'_> {
        self.0.ast.bond(idx)
    }

    pub fn neighbors(&self, atom: AtomIdx) -> impl Iterator<Item = NeighborView<'_>> {
        self.0.ast.neighbors(atom)
    }

    pub fn graph(&self) -> &Graph {
        self.0.ast.graph()
    }

    pub fn bond_order_sum(&self, atom: AtomIdx) -> Option<u8> {
        self.0.ast.bond_order_sum(atom)
    }

    pub fn dative_bond_order_sums(&self, atom: AtomIdx) -> (u8, u8) {
        self.0.ast.dative_bond_order_sums(atom)
    }

    pub fn is_in_aromatic_system(&self, atom: AtomIdx) -> bool {
        self.0.ast.is_in_aromatic_system(atom)
    }

    pub fn to_edn_str(&self) -> String {
        MoleculeAstWrapper::from_ast(self.0.ast.clone())
            .to_edn()
            .to_string()
    }

    pub fn from_edn_str(input: &str) -> Result<Self, MoleculeEdnError> {
        let (ast, _) = parse_molecule_dsl(input)?.into_parts();
        Ok(Self::new(ast)?)
    }
}

impl PartialEq for Molecule {
    fn eq(&self, other: &Self) -> bool {
        self.0.ast == other.0.ast
    }
}

impl Eq for Molecule {}

impl Index<AtomIdx> for Molecule {
    type Output = AtomAst;
    fn index(&self, idx: AtomIdx) -> &AtomAst {
        &self.0.ast[idx]
    }
}

impl Index<BondIdx> for Molecule {
    type Output = BondAst;
    fn index(&self, idx: BondIdx) -> &BondAst {
        &self.0.ast[idx]
    }
}

impl Index<DativeBondIdx> for Molecule {
    type Output = BondAst;
    fn index(&self, idx: DativeBondIdx) -> &BondAst {
        &self.0.ast[idx]
    }
}

impl Index<NoncovalentBondIdx> for Molecule {
    type Output = BondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &BondAst {
        &self.0.ast[idx]
    }
}

impl Index<AromaticSystemIdx> for Molecule {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        &self.0.ast[idx]
    }
}

impl Index<MulticenterBondIdx> for Molecule {
    type Output = MulticenterBondAst;
    fn index(&self, idx: MulticenterBondIdx) -> &MulticenterBondAst {
        &self.0.ast[idx]
    }
}

#[cfg(test)]
mod tests {
    use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
        }
    }

    fn ground_bond(order: i64) -> BondAst {
        BondAst {
            order: ValueAst::Lit(order),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
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
        let r1 = mol.rings();
        let r2 = mol.rings();
        assert!(Arc::ptr_eq(&r1, &r2));
    }

    #[test]
    fn test_molecule_views_delegate() {
        let ast = MoleculeAst::new(
            vec![ground_atom(), ground_atom()],
            vec![(AtomIdx(0), AtomIdx(1), ground_bond(1))],
            vec![], vec![], vec![], vec![], vec![],
        );
        let mol = Molecule::new(ast).unwrap();
        assert_eq!(mol.atoms().count(), 2);
        assert_eq!(mol.bonds().count(), 1);
        assert_eq!(mol.neighbors(AtomIdx(0)).count(), 1);
        assert_eq!(mol.bond_order_sum(AtomIdx(0)), Some(1));
    }

    #[test]
    fn test_molecule_partial_eq() {
        let ast1 = MoleculeAst::new(
            vec![ground_atom()],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let ast2 = ast1.clone();
        let m1 = Molecule::new(ast1).unwrap();
        let m2 = Molecule::new(ast2).unwrap();
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_molecule_to_edn_str() {
        let ast = MoleculeAst::new(
            vec![ground_atom()],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let mol = Molecule::new(ast).unwrap();
        let text = mol.to_edn_str();
        assert!(!text.is_empty());
        assert!(text.contains("atoms"));
    }

    #[test]
    fn test_molecule_from_edn_str_parse_error() {
        assert!(matches!(
            Molecule::from_edn_str("not valid edn"),
            Err(MoleculeEdnError::Parse(_))
        ));
    }
}
