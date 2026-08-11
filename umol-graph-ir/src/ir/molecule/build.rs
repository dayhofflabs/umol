//! Construction builder — build a fresh `Molecule` with a bare-verb-adds convention:
//! every method *adds/declares* (there is no lookup — that is `MoleculeEditor`). Wraps a
//! `MoleculeEditor` and lowers each call onto it.

use super::super::aromatic::AromaticSystemForm;
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::constraint::BondConstraintForm;
use super::super::dative::DativeBondForm;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::multicenter::MulticenterBondForm;
use super::super::noncovalent::NoncovalentBondForm;
use super::super::stereo::{StereoAtomForm, StereoBondForm};
use super::{Molecule, MoleculeEditor};

/// Build a molecule from scratch. `atom` adds an atom and hands back its handle; the
/// per-family verbs (`single`/`double`/`triple`, `dative`) add bonds between handles;
/// `build` finalizes. Contrast `MoleculeEditor`, which *mutates an existing* molecule
/// (`add_*` plus `atom(id)`/`bond(id)` lookups).
pub struct MoleculeBuilder {
    editor: MoleculeEditor,
    /// Whether each added atom is grounded — its unspecified fields filled with their
    /// ground defaults (neutral, singlet, …). Applied per atom via `AtomForm::into_ground`.
    ground: bool,
}

impl MoleculeBuilder {
    /// A fresh, empty builder. Atom fields left unspecified stay open for resolution.
    pub fn new() -> Self {
        Self {
            editor: Molecule::new().edit(),
            ground: false,
        }
    }

    /// A fresh builder that **grounds** each atom — unspecified fields take their ground
    /// defaults (neutral, singlet, …), reusing `AtomForm::into_ground`. The in-code
    /// analogue of the `mol_dsl_ground!` path. (Partial defaults — e.g. neutral but spin
    /// open — will be the L2 `+`-spec's `charge(0)`/`spin(…)` terms.)
    pub fn ground() -> Self {
        Self {
            editor: Molecule::new().edit(),
            ground: true,
        }
    }

    /// Add an atom, returning its handle. Accepts an element (`C`), an `AtomForm`, or a
    /// compact atom-string (`"C#h3"`) via `Into<AtomForm>`.
    pub fn atom(&mut self, spec: impl Into<AtomForm>) -> AtomId {
        let atom = spec.into();
        let atom = if self.ground {
            atom.into_ground()
        } else {
            atom
        };
        self.editor.add_atom(atom)
    }

    /// Add a single (order-1) bond.
    pub fn single(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondForm::from_order(1))
    }

    /// Add a double (order-2) bond.
    pub fn double(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondForm::from_order(2))
    }

    /// Add a triple (order-3) bond.
    pub fn triple(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(first, second, BondForm::from_order(3))
    }

    /// Add a bond carrying an explicit `BondForm` — the escape hatch for a bond with a
    /// charge, spin, order-set, or constraint that `single`/`double`/`triple` can't set.
    pub fn bond(&mut self, first: AtomId, second: AtomId, bond: impl Into<BondForm>) -> BondId {
        self.editor.add_bond(first, second, bond.into())
    }

    /// Add an aromatic bond — order 1 with the aromatic flag set (`1#a`); resolution
    /// perceives the aromatic system. Distinct from, and not exclusive with, `aromatic_system`.
    pub fn aromatic_bond(&mut self, first: AtomId, second: AtomId) -> BondId {
        self.editor.add_bond(
            first,
            second,
            BondForm::from_order(1).with_constraint(BondConstraintForm::aromatic(true)),
        )
    }

    /// Add a dative bond from `donors` to `acceptor` — its own family, not a bond order. Carries
    /// `attributes` (a `DativeBondForm` or a DSL spec string).
    pub fn dative_bond(
        &mut self,
        donors: impl IntoIterator<Item = AtomId>,
        acceptor: AtomId,
        attributes: impl Into<DativeBondForm>,
    ) -> DativeBondId {
        self.editor
            .add_dative_bond(donors.into_iter().collect(), acceptor, attributes.into())
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

    /// Add an aromatic-system overlay over `atoms`, carrying `attributes` (an `AromaticSystemForm` — e.g.
    /// `from_electrons([1, 1, 1, 1, 1, 1])` for benzene — or a DSL spec string). The σ-framework is
    /// separate — add those bonds with the bond verbs.
    pub fn aromatic_system(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
        attributes: impl Into<AromaticSystemForm>,
    ) -> AromaticSystemId {
        self.editor
            .add_aromatic_system(atoms.into_iter().collect(), attributes.into())
    }

    /// Add a multicenter-bond overlay over `atoms`, carrying `attributes` (a `MulticenterBondForm` or a DSL
    /// spec string).
    pub fn multicenter_bond(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
        attributes: impl Into<MulticenterBondForm>,
    ) -> MulticenterBondId {
        self.editor
            .add_multicenter_bond(atoms.into_iter().collect(), attributes.into())
    }

    /// Add a noncovalent-bond overlay between `first` and `second`, carrying `attributes` (a
    /// `NoncovalentBondForm` — e.g. `from_kind(...)` — or a DSL spec string).
    pub fn noncovalent_bond(
        &mut self,
        first: AtomId,
        second: AtomId,
        attributes: impl Into<NoncovalentBondForm>,
    ) -> NoncovalentBondId {
        self.editor
            .add_noncovalent_bond([first, second], attributes.into())
    }

    /// Add a stereo-atom overlay: an atom `site` with its ordered `ligands` and a configuration
    /// `attributes`. Delegates to the editor's `add_stereo_atom`.
    pub fn stereo_atom(
        &mut self,
        site: AtomId,
        ligands: impl IntoIterator<Item = StereoLigand>,
        attributes: impl Into<StereoAtomForm>,
    ) -> StereoAtomId {
        self.editor
            .add_stereo_atom(site, ligands.into_iter().collect(), attributes.into())
    }

    /// Add a stereo-bond overlay: a bond `site` with its ordered `ligands` and a configuration
    /// `attributes`. Delegates to the editor's `add_stereo_bond`.
    pub fn stereo_bond(
        &mut self,
        site: BondId,
        ligands: impl IntoIterator<Item = StereoLigand>,
        attributes: impl Into<StereoBondForm>,
    ) -> StereoBondId {
        self.editor
            .add_stereo_bond(site, ligands.into_iter().collect(), attributes.into())
    }

    /// Finalize into a `Molecule`. Unspecified atom fields stay open for resolution.
    pub fn build(self) -> Molecule {
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
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::noncovalent::NoncovalentBondKind;
    use crate::ir::num::NumForm;
    use crate::ir::stereo::{StereoCoset, StereoKind};

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
        assert_eq!(mol.bond(BondId(0)).attributes, &BondForm::from_order(1));
        assert_eq!(mol.bond(BondId(1)).attributes, &BondForm::from_order(2));
    }

    #[rstest]
    #[case::undetermined_grounds(AtomForm::from_element(Element::C), NumForm::Lit(0))]
    #[case::preset_charge_preserved(
        AtomForm::from_element(Element::C).with_charge(2_i64),
        NumForm::Lit(2)
    )]
    fn test_molecule_builder_ground(#[case] spec: AtomForm, #[case] expected_charge: NumForm) {
        let mut builder = MoleculeBuilder::ground();
        let atom = builder.atom(spec);
        let mol = builder.build();

        assert_eq!(mol.atom(atom).attributes.charge, expected_charge);
        // an unspecified field is grounded regardless of the preset charge
        assert_eq!(
            mol.atom(atom).attributes.implicit_hydrogens,
            NumForm::Lit(0)
        );
    }

    #[rstest]
    fn test_molecule_builder_bond() {
        let mut builder = MoleculeBuilder::new();
        let a = builder.atom(Element::C);
        let b = builder.atom(Element::O);
        // charge is reachable only through the explicit-BondForm escape hatch
        let bond = builder.bond(a, b, BondForm::from_order(2).with_charge(-1_i64));
        let mol = builder.build();

        assert_eq!(bond, BondId(0));
        assert_eq!(
            mol.bond(bond).attributes,
            &BondForm::from_order(2).with_charge(-1_i64)
        );
    }

    #[rstest]
    fn test_molecule_builder_aromatic_bond() {
        let mut builder = MoleculeBuilder::new();
        let a = builder.atom(Element::C);
        let b = builder.atom(Element::C);
        let bond = builder.aromatic_bond(a, b);
        let mol = builder.build();

        assert_eq!(bond, BondId(0));
        assert_eq!(
            mol.bond(bond).attributes,
            &BondForm::from_order(1).with_constraint(BondConstraintForm::aromatic(true))
        );
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
        let expected: Vec<[AtomId; 2]> = expected_edges
            .iter()
            .map(|&(i, j)| [atoms[i], atoms[j]])
            .collect();
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
        let expected: Vec<[AtomId; 2]> = expected_edges
            .iter()
            .map(|&(i, j)| [atoms[i], atoms[j]])
            .collect();
        assert_eq!(endpoints, expected);
    }

    #[rstest]
    fn test_molecule_builder_aromatic_system() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::C);
        let a1 = builder.atom(Element::C);
        let system =
            builder.aromatic_system([a0, a1], AromaticSystemForm::from_electrons(vec![1, 1]));
        let mol = builder.build();

        assert_eq!(system, AromaticSystemId(0));
        assert_eq!(mol.aromatic_systems().count(), 1);
        assert_eq!(
            mol.aromatic_system(system).attributes,
            &AromaticSystemForm::from_electrons(vec![1, 1])
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
    fn test_molecule_builder_multicenter_bond() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::B);
        let a1 = builder.atom(Element::B);
        let a2 = builder.atom(Element::H);
        let bond = builder.multicenter_bond(
            [a0, a1, a2],
            MulticenterBondForm::from_electrons(vec![1, 1, 1]),
        );
        let mol = builder.build();

        assert_eq!(bond, MulticenterBondId(0));
        assert_eq!(
            mol.multicenter_bond(bond).attributes,
            &MulticenterBondForm::from_electrons(vec![1, 1, 1])
        );
    }

    #[rstest]
    fn test_molecule_builder_noncovalent_bond() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.atom(Element::O);
        let a1 = builder.atom(Element::H);
        let bond = builder.noncovalent_bond(
            a0,
            a1,
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let mol = builder.build();

        assert_eq!(bond, NoncovalentBondId(0));
        assert_eq!(
            mol.noncovalent_bond(bond).attributes,
            &NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
        );
    }

    #[rstest]
    fn test_molecule_builder_stereo_atom() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.atom(Element::C);
        let f = builder.atom(Element::F);
        let cl = builder.atom(Element::Cl);
        let br = builder.atom(Element::Br);
        let i = builder.atom(Element::I);
        let stereo = builder.stereo_atom(
            c,
            [
                StereoLigand::new(f, StereoLigandKind::Atom),
                StereoLigand::new(cl, StereoLigandKind::Atom),
                StereoLigand::new(br, StereoLigandKind::Atom),
                StereoLigand::new(i, StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        );
        let mol = builder.build();

        assert_eq!(stereo, StereoAtomId(0));
        assert_eq!(mol.stereo_atoms().count(), 1);
        assert_eq!(mol.stereo_atom(stereo).site_id(), c);
        assert_eq!(
            mol.stereo_atom(stereo).atom_ids().collect::<Vec<_>>(),
            vec![c, f, cl, br, i]
        );
        assert_eq!(
            mol.stereo_atom(stereo).attributes,
            &StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
        );
    }

    #[rstest]
    fn test_molecule_builder_stereo_bond() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.atom(Element::C);
        let c2 = builder.atom(Element::C);
        let f = builder.atom(Element::F);
        let h = builder.atom(Element::H);
        let bond = builder.double(c1, c2);
        let stereo = builder.stereo_bond(
            bond,
            [
                StereoLigand::new(f, StereoLigandKind::Atom),
                StereoLigand::new(c1, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(h, StereoLigandKind::Atom),
                StereoLigand::new(c2, StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let mol = builder.build();

        assert_eq!(stereo, StereoBondId(0));
        assert_eq!(mol.stereo_bonds().count(), 1);
        assert_eq!(mol.stereo_bond(stereo).site_id(), bond);
        assert_eq!(
            mol.stereo_bond(stereo).attributes,
            &StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))
        );
    }
}
