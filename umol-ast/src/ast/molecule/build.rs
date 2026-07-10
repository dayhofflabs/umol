//! Construction builder — build a fresh `MoleculeAst` with a bare-verb-adds convention:
//! every method *adds/declares* (there is no lookup — that is `MoleculeEditor`). Wraps a
//! `MoleculeEditor` and lowers each call onto it.

use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::dative::DativeBondAst;
use super::super::id::{AtomId, BondId, DativeBondId};
use super::{MoleculeAst, MoleculeEditor};

/// Build a molecule from scratch. `atom` adds an atom and hands back its handle; the
/// per-family verbs (`single`/`double`/`triple`, `dative`) add bonds between handles;
/// `build` finalizes. Contrast `MoleculeEditor`, which *mutates an existing* molecule
/// (`add_*` plus `atom(id)`/`bond(id)` lookups).
pub struct MoleculeBuilder {
    editor: MoleculeEditor,
}

impl MoleculeBuilder {
    /// A fresh, empty builder.
    pub fn new() -> Self {
        Self {
            editor: MoleculeAst::new().edit(),
        }
    }

    /// Add an atom, returning its handle. Accepts an element (`C`), an `AtomAst`, or a
    /// compact atom-string (`"C#h3"`) via `Into<AtomAst>`.
    pub fn atom(&mut self, spec: impl Into<AtomAst>) -> AtomId {
        self.editor.add_atom(spec.into())
    }

    /// Add a single (order-1) bond.
    pub fn single(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondAst::from_order(1))
    }

    /// Add a double (order-2) bond.
    pub fn double(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondAst::from_order(2))
    }

    /// Add a triple (order-3) bond.
    pub fn triple(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondAst::from_order(3))
    }

    /// Add a dative bond from `donors` to `acceptor` — its own family, not a bond order.
    pub fn dative(
        &mut self,
        donors: impl IntoIterator<Item = AtomId>,
        acceptor: AtomId,
    ) -> DativeBondId {
        self.editor.add_dative_bond(
            donors.into_iter().collect(),
            acceptor,
            DativeBondAst::default(),
        )
    }

    /// Finalize into a `MoleculeAst`. Unspecified atom fields stay open for resolution.
    pub fn build(self) -> MoleculeAst {
        self.editor.build()
    }
}

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;

    #[rstest]
    fn test_molecule_builder() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.atom(Element::C);
        let c2 = b.atom("C"); // From<&str>
        let o = b.atom(Element::O);
        b.single(c1, c2);
        b.double(c2, o);
        let mol = b.build();

        assert_eq!(mol.atoms().count(), 3);
        assert_eq!(mol.bonds().count(), 2);
        assert_eq!(mol.bond(BondId(0)).ast, &BondAst::from_order(1));
        assert_eq!(mol.bond(BondId(1)).ast, &BondAst::from_order(2));
    }
}
