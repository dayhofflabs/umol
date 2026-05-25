//! Counts valence resolver: per-atom electron-conservation enumeration
//! ([`ValenceInvariants`]), with implicit hydrogens inferred against the
//! element's `covalence_set` and the closed-shell state — fewest unpaired
//! electrons, then most lone pairs — selected.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomConstraints, AtomId, AtomView, Lattice,
    MoleculeAst, ValueAst,
};
use umol_shared::element::Element;

use super::invariants::ValenceInvariants;
use crate::ops::valence::table::ValenceTable;

#[derive(Clone, Debug)]
pub struct CountsValence {
    pub table: ValenceTable,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no valid valence state")]
    NoValidValenceState,
}

impl CountsValence {
    pub fn new(table: ValenceTable) -> Self {
        Self { table }
    }

    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), CountsError> {
        for i in 0..ast.atoms().count() as u32 {
            self.resolve_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    fn resolve_atom(&self, ast: &mut MoleculeAst, atom_id: AtomId) -> Result<(), CountsError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{AtomAst, AtomId, BondAst, Constraints, MoleculeAst};
    use umol_ast::{mol, mol_ground};
    use umol_shared::element::Element;

    use super::*;
    use crate::valence_table;

    fn carbon_methane_with_undetermined() -> MoleculeAst {
        // Carbon with implicit_hydrogens left undetermined; valence = 0 (no
        // bonds). Counts resolver should infer 4 implicit Hs.
        mol!(r#"{:atoms ["C"] :bonds []}"#)
    }

    fn ethane() -> MoleculeAst {
        let mut a = AtomAst::from_element(Element::C);
        let mut b = AtomAst::from_element(Element::C);
        a.implicit_hydrogens = ValueAst::Undetermined;
        b.implicit_hydrogens = ValueAst::Undetermined;
        MoleculeAst::from_parts(
            vec![a, b],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_ground_passthrough() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValence::new(table);
        let mut ast = mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_methane_implicit_h() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValence::new(table);
        let mut ast = carbon_methane_with_undetermined();
        resolver.resolve(&mut ast).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert!(matches!(atom.implicit_hydrogens, ValueAst::Lit(4)));
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_ethane_implicit_h() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValence::new(table);
        let mut ast = ethane();
        resolver.resolve(&mut ast).unwrap();
        for i in 0..2 {
            let atom = ast.atom(AtomId(i)).ast;
            assert!(matches!(atom.implicit_hydrogens, ValueAst::Lit(3)));
        }
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_out_of_scope_element() {
        let table = valence_table! { C => [4] };
        let resolver = CountsValence::new(table);
        let si = AtomAst::from_element(Element::Si);
        let mut ast = MoleculeAst::from_atoms_and_bonds(vec![si.clone()], vec![]);
        resolver.resolve(&mut ast).unwrap();
        // Si is absent from the table, so it is out of scope and left untouched.
        assert_eq!(ast.atom(AtomId(0)).ast, &si);
    }
}
