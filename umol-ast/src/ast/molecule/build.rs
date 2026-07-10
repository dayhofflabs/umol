//! Construction builder — build a fresh `MoleculeAst` with a bare-verb-adds convention:
//! every method *adds/declares* (there is no lookup — that is `MoleculeEditor`). Wraps a
//! `MoleculeEditor` and lowers each call onto it.

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::dative::DativeBondAst;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
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

    /// Wire `atoms` into a path of single bonds (each consecutive pair). Returns the
    /// bond handles in order — empty for fewer than two atoms.
    pub fn chain(&mut self, atoms: impl IntoIterator<Item = AtomId>) -> Vec<BondId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        atoms
            .windows(2)
            .map(|pair| self.single(pair[0], pair[1]))
            .collect()
    }

    /// Wire `atoms` into a ring of single bonds — a `chain` plus the closing bond from
    /// the last atom back to the first. Returns the bond handles in order.
    pub fn ring(&mut self, atoms: impl IntoIterator<Item = AtomId>) -> Vec<BondId> {
        let atoms: Vec<AtomId> = atoms.into_iter().collect();
        let mut bonds = self.chain(atoms.iter().copied());
        if atoms.len() > 2 {
            bonds.push(self.single(atoms[atoms.len() - 1], atoms[0]));
        }
        bonds
    }

    /// Add an aromatic-system overlay over `atoms`, one π-`electrons` count per atom
    /// (`[1, 1, 1, 1, 1, 1]` for benzene). The σ-framework is separate — add those
    /// bonds with the bond verbs.
    pub fn aromatic_system(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
        electrons: impl IntoIterator<Item = i64>,
    ) -> AromaticSystemId {
        self.editor.add_aromatic_system(
            atoms.into_iter().collect(),
            AromaticSystemAst::from_electrons(electrons.into_iter().collect()),
        )
    }

    /// Add a multicenter-bond overlay over `atoms`, one `electrons` count per atom.
    pub fn multicenter(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
        electrons: impl IntoIterator<Item = i64>,
    ) -> MulticenterBondId {
        self.editor.add_multicenter_bond(
            atoms.into_iter().collect(),
            MulticenterBondAst::from_electrons(electrons.into_iter().collect()),
        )
    }

    /// Add a noncovalent-bond overlay of `kind` between `first` and `second`.
    pub fn noncovalent(
        &mut self,
        first: AtomId,
        second: AtomId,
        kind: NoncovalentBondKind,
    ) -> NoncovalentBondId {
        self.editor
            .add_noncovalent_bond([first, second], NoncovalentBondAst::from_kind(kind))
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
    use crate::ast::value::ValueAst;

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
    #[case::empty(0, vec![])]
    #[case::single_atom(1, vec![])]
    #[case::path(4, vec![(0, 1), (1, 2), (2, 3)])]
    fn test_molecule_builder_chain(
        #[case] atom_count: usize,
        #[case] expected_edges: Vec<(usize, usize)>,
    ) {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomId> = (0..atom_count).map(|_| builder.atom(Element::C)).collect();
        let bonds = builder.chain(atoms.iter().copied());
        let mol = builder.build();

        let endpoints: Vec<[AtomId; 2]> = bonds.iter().map(|&b| mol.bond(b).atom_ids()).collect();
        let expected: Vec<[AtomId; 2]> =
            expected_edges.iter().map(|&(i, j)| [atoms[i], atoms[j]]).collect();
        assert_eq!(endpoints, expected);
    }

    #[rstest]
    #[case::two_atoms_no_closure(2, vec![(0, 1)])]
    #[case::triangle(3, vec![(0, 1), (1, 2), (0, 2)])]
    fn test_molecule_builder_ring(
        #[case] atom_count: usize,
        #[case] expected_edges: Vec<(usize, usize)>,
    ) {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomId> = (0..atom_count).map(|_| builder.atom(Element::C)).collect();
        let bonds = builder.ring(atoms.iter().copied());
        let mol = builder.build();

        let endpoints: Vec<[AtomId; 2]> = bonds.iter().map(|&b| mol.bond(b).atom_ids()).collect();
        let expected: Vec<[AtomId; 2]> =
            expected_edges.iter().map(|&(i, j)| [atoms[i], atoms[j]]).collect();
        assert_eq!(endpoints, expected);
    }

    #[rstest]
    fn test_molecule_builder_aromatic_system() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::C);
        let a1 = builder.atom(Element::C);
        let system = builder.aromatic_system([a0, a1], [1, 1]);
        let mol = builder.build();

        assert_eq!(system, AromaticSystemId(0));
        assert_eq!(mol.aromatic_systems().count(), 1);
        assert_eq!(
            mol.aromatic_system(system).ast,
            &AromaticSystemAst::from_electrons(vec![1, 1])
        );
        assert_eq!(
            mol.aromatic_system(system)
                .atoms()
                .map(|view| view.id)
                .collect::<Vec<_>>(),
            vec![a0, a1]
        );
    }

    #[rstest]
    fn test_molecule_builder_multicenter() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::B);
        let a1 = builder.atom(Element::B);
        let a2 = builder.atom(Element::H);
        let bond = builder.multicenter([a0, a1, a2], [1, 1, 1]);
        let mol = builder.build();

        assert_eq!(bond, MulticenterBondId(0));
        assert_eq!(
            mol.multicenter_bond(bond).ast,
            &MulticenterBondAst::from_electrons(vec![1, 1, 1])
        );
    }

    #[rstest]
    fn test_molecule_builder_noncovalent() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::O);
        let a1 = builder.atom(Element::H);
        let bond = builder.noncovalent(a0, a1, NoncovalentBondKind::HydrogenBond);
        let mol = builder.build();

        assert_eq!(bond, NoncovalentBondId(0));
        assert_eq!(
            mol.noncovalent_bond(bond).ast,
            &NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
        );
    }
}
