//! Construction builder — build a fresh `MoleculeAst` with a bare-verb-adds convention:
//! every method *adds/declares* (there is no lookup — that is `MoleculeEditor`). Wraps a
//! `MoleculeEditor` and lowers each call onto it.

use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::{AromaticValenceAst, AtomConstraintAst};
use super::super::dative::DativeBondAst;
use super::super::id::{AtomId, BondId, DativeBondId};
use super::super::value::ValueAst;
use super::{MoleculeAst, MoleculeEditor};

/// Build a molecule from scratch. `atom` adds an atom and hands back its handle; the
/// per-family verbs (`single`/`double`/`triple`, `dative`) add bonds between handles;
/// `build` finalizes. Contrast `MoleculeEditor`, which *mutates an existing* molecule
/// (`add_*` plus `atom(id)`/`bond(id)` lookups).
pub struct MoleculeBuilder {
    editor: MoleculeEditor,
    /// Whether each added atom is grounded — its unspecified fields filled with their
    /// ground defaults (neutral, singlet, …). Applied per atom via `AtomAst::into_ground`.
    ground: bool,
}

impl MoleculeBuilder {
    /// A fresh, empty builder. Atom fields left unspecified stay open for resolution.
    pub fn new() -> Self {
        Self {
            editor: MoleculeAst::new().edit(),
            ground: false,
        }
    }

    /// A fresh builder that **grounds** each atom — unspecified fields take their ground
    /// defaults (neutral, singlet, …), reusing `AtomAst::into_ground`. The in-code
    /// analogue of the `mol_dsl_ground!` path. (Partial defaults — e.g. neutral but spin
    /// open — will be the L2 `+`-spec's `charge(0)`/`spin(…)` terms.)
    pub fn ground() -> Self {
        Self {
            editor: MoleculeAst::new().edit(),
            ground: true,
        }
    }

    /// Add an atom, returning its handle. Accepts an element (`C`), an `AtomAst`, or a
    /// compact atom-string (`"C#h3"`) via `Into<AtomAst>`.
    pub fn atom(&mut self, spec: impl Into<AtomAst>) -> AtomId {
        let atom = spec.into();
        let atom = if self.ground { atom.into_ground() } else { atom };
        self.editor.add_atom(atom)
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

    /// Close a single-bonded ring through `atoms` (each consecutive pair plus the
    /// wrap-around) and flag every ring atom aromatic —
    /// `AromaticValence(Aromatic(Undetermined))`, the "aromatic, count open" marker.
    /// The `AromaticSystemAst` overlay is a resolution-time product (ring-electron
    /// counts are not known at construction), so it is not laid down here. Returns the
    /// ring-bond handles in ring order.
    pub fn aromatic_ring(&mut self, atoms: impl IntoIterator<Item = AtomId>) -> Vec<BondId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        let mut bonds = Vec::with_capacity(atoms.len());
        for (position, &atom) in atoms.iter().enumerate() {
            self.editor.atom_mut(atom).ast.constraints.set(
                AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(
                    ValueAst::Undetermined,
                )),
            );
            let next = atoms[(position + 1) % atoms.len()];
            bonds.push(self.editor.add_bond(atom, next, BondAst::from_order(1)));
        }
        bonds
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
    use crate::ast::constraint::AtomConstraintsAst;

    #[rstest]
    fn test_molecule_builder() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.atom(Element::C);
        let c2 = builder.atom("C"); // From<&str>
        let o = builder.atom(Element::O);
        builder.single(c1, c2);
        builder.double(c2, o);
        let mol = builder.build();

        assert_eq!(mol.atoms().count(), 3);
        assert_eq!(mol.bonds().count(), 2);
        assert_eq!(mol.bond(BondId(0)).ast, &BondAst::from_order(1));
        assert_eq!(mol.bond(BondId(1)).ast, &BondAst::from_order(2));
    }

    #[rstest]
    #[case::undetermined_grounds(AtomAst::from_element(Element::C), ValueAst::Lit(0))]
    #[case::preset_charge_preserved(
        AtomAst::from_element(Element::C).with_charge(2_i64),
        ValueAst::Lit(2)
    )]
    fn test_molecule_builder_ground(#[case] spec: AtomAst, #[case] expected_charge: ValueAst) {
        let mut builder = MoleculeBuilder::ground();
        let atom = builder.atom(spec);
        let mol = builder.build();

        assert_eq!(mol.atom(atom).ast.charge, expected_charge);
        // an unspecified field is grounded regardless of the preset charge
        assert_eq!(mol.atom(atom).ast.implicit_hydrogens, ValueAst::Lit(0));
    }

    #[rstest]
    fn test_molecule_builder_aromatic_ring() {
        let mut builder = MoleculeBuilder::new();
        let ring: Vec<AtomId> = (0..6).map(|_| builder.atom(Element::C)).collect();
        let bonds = builder.aromatic_ring(ring.iter().copied());
        let mol = builder.build();

        assert_eq!(
            bonds,
            vec![
                BondId(0),
                BondId(1),
                BondId(2),
                BondId(3),
                BondId(4),
                BondId(5),
            ]
        );
        assert_eq!(mol.bonds().count(), 6);
        for &bond in &bonds {
            assert_eq!(mol.bond(bond).ast, &BondAst::from_order(1));
        }
        for &atom in &ring {
            assert_eq!(
                mol.atom(atom).ast.constraints,
                AtomConstraintsAst::from(AtomConstraintAst::AromaticValence(
                    AromaticValenceAst::Aromatic(ValueAst::Undetermined)
                ))
            );
        }
    }
}
