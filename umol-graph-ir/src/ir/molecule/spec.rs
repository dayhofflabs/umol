//! Declarative `+`-spec construction — build a `Molecule` by summing free-function
//! *terms* into a `MoleculeSpec`, then materializing. Each term lowers onto the L1
//! `MoleculeBuilder`; every atom position is an [`AtomArg`] that either creates a fresh atom
//! (optionally named) or references one already introduced (by position or by name).

use std::collections::HashMap;
use std::ops::Add;

use umol_chem::element::Element;

use super::super::aromatic::AromaticSystemForm;
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::constraint::BondConstraintForm;
use super::super::dative::DativeBondForm;
use super::super::id::{AtomId, BondId};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::multicenter::MulticenterBondForm;
use super::super::noncovalent::NoncovalentBondForm;
use super::super::stereo::{StereoAtomForm, StereoBondForm};
use super::{Molecule, MoleculeBuilder};

/// An atom argument to a spec term: create a fresh atom (optionally named) or reference
/// one already introduced — by creation `position` or by `name`. What you write picks the
/// variant by type: `C`/`"C#h3"`/`AtomForm` → create, `(name, spec)` tuple → create-named,
/// a bare integer → by position, [`name`] → by name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomArg {
    New {
        spec: AtomForm,
        name: Option<String>,
    },
    Index(u32),
    Name(String),
}

impl From<Element> for AtomArg {
    fn from(element: Element) -> Self {
        Self::New {
            spec: AtomForm::from_element(element),
            name: None,
        }
    }
}

impl From<&str> for AtomArg {
    fn from(spec: &str) -> Self {
        Self::New {
            spec: AtomForm::from(spec),
            name: None,
        }
    }
}

impl From<AtomForm> for AtomArg {
    fn from(spec: AtomForm) -> Self {
        Self::New { spec, name: None }
    }
}

impl From<u32> for AtomArg {
    fn from(position: u32) -> Self {
        Self::Index(position)
    }
}

impl From<i32> for AtomArg {
    fn from(position: i32) -> Self {
        Self::Index(
            u32::try_from(position)
                .unwrap_or_else(|_| panic!("atom position must be non-negative, got {position}")),
        )
    }
}

impl<S: Into<String>, T: Into<AtomForm>> From<(S, T)> for AtomArg {
    fn from((name, spec): (S, T)) -> Self {
        Self::New {
            spec: spec.into(),
            name: Some(name.into()),
        }
    }
}

/// Reference an already-introduced atom by name.
pub fn name(name: impl Into<String>) -> AtomArg {
    AtomArg::Name(name.into())
}

/// A reference to an existing bond, by name. Bonds pre-exist as path edges — a spec term never creates
/// one — so this only references (no `New` / positional `Index`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondArg {
    Name(String),
}

impl From<&str> for BondArg {
    fn from(name: &str) -> Self {
        Self::Name(name.to_string())
    }
}

impl From<String> for BondArg {
    fn from(name: String) -> Self {
        Self::Name(name)
    }
}

/// A ligand argument to a stereo term: an atom (referenced like any [`AtomArg`]) or a virtual ligand
/// (implicit hydrogen / lone pair). A virtual ligand's bearing atom is filled at build — the site atom
/// for a stereo atom, or the bond atom its position selects for a stereo bond.
// The `Atom(AtomArg)` variant dwarfs the unit ones, but this is a transient construction-time enum;
// boxing it would only clutter every ligand literal.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StereoLigandArg {
    Atom(AtomArg),
    ImplicitHydrogen,
    LonePair,
}

impl From<AtomArg> for StereoLigandArg {
    fn from(atom: AtomArg) -> Self {
        Self::Atom(atom)
    }
}

/// A single term of a [`MoleculeSpec`]: an introduction, a relation, or a default.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeSpecTerm {
    Atoms(Vec<AtomArg>),
    Bond {
        name: Option<String>,
        first: AtomArg,
        second: AtomArg,
        ast: BondForm,
    },
    Chain(Vec<AtomArg>),
    Ring(Vec<AtomArg>),
    DativeBond {
        donors: Vec<AtomArg>,
        acceptor: AtomArg,
        ast: DativeBondForm,
    },
    AromaticSystem {
        atoms: Vec<AtomArg>,
        ast: AromaticSystemForm,
    },
    MulticenterBond {
        atoms: Vec<AtomArg>,
        ast: MulticenterBondForm,
    },
    NoncovalentBond {
        first: AtomArg,
        second: AtomArg,
        ast: NoncovalentBondForm,
    },
    StereoAtom {
        site: AtomArg,
        ligands: Vec<StereoLigandArg>,
        ast: StereoAtomForm,
    },
    StereoBond {
        site: BondArg,
        ligands: Vec<StereoLigandArg>,
        ast: StereoBondForm,
    },
    Ground,
}

/// Introduce one atom.
pub fn atom(spec: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Atoms(vec![spec.into()])
}

/// Introduce a homogeneous list of atoms.
pub fn atoms(specs: impl IntoIterator<Item = impl Into<AtomArg>>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Atoms(specs.into_iter().map(Into::into).collect())
}

/// A single (order-1) bond.
pub fn single(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: None,
        first: first.into(),
        second: second.into(),
        ast: BondForm::from_order(1),
    }
}

/// A double (order-2) bond.
pub fn double(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: None,
        first: first.into(),
        second: second.into(),
        ast: BondForm::from_order(2),
    }
}

/// A triple (order-3) bond.
pub fn triple(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: None,
        first: first.into(),
        second: second.into(),
        ast: BondForm::from_order(3),
    }
}

/// A bond carrying an explicit `BondForm` — the escape hatch for a charge, spin, order-set,
/// or constraint the order verbs can't set.
pub fn bond(
    first: impl Into<AtomArg>,
    second: impl Into<AtomArg>,
    ast: impl Into<BondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: None,
        first: first.into(),
        second: second.into(),
        ast: ast.into(),
    }
}

/// A bond that binds a `name` in the molecule's shared atom/bond namespace, so an overlay (a stereo
/// bond) can reference it as its site. Otherwise identical to [`bond`].
pub fn named_bond(
    name: impl Into<String>,
    first: impl Into<AtomArg>,
    second: impl Into<AtomArg>,
    ast: impl Into<BondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: Some(name.into()),
        first: first.into(),
        second: second.into(),
        ast: ast.into(),
    }
}

/// An aromatic bond — order 1 with the aromatic flag (`1#a`); resolution perceives the
/// aromatic system. Not exclusive with [`aromatic_system`].
pub fn aromatic_bond(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        name: None,
        first: first.into(),
        second: second.into(),
        ast: BondForm::from_order(1).with_constraint(BondConstraintForm::aromatic(true)),
    }
}

/// Wire the atoms into a path of single bonds.
pub fn chain(specs: impl IntoIterator<Item = impl Into<AtomArg>>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Chain(specs.into_iter().map(Into::into).collect())
}

/// Wire the atoms into a ring of single bonds (path + closing bond).
pub fn ring(specs: impl IntoIterator<Item = impl Into<AtomArg>>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Ring(specs.into_iter().map(Into::into).collect())
}

/// A dative bond from `donors` to `acceptor`, carrying `ast` (a `DativeBondForm` or a DSL spec string).
pub fn dative_bond(
    donors: impl IntoIterator<Item = impl Into<AtomArg>>,
    acceptor: impl Into<AtomArg>,
    ast: impl Into<DativeBondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::DativeBond {
        donors: donors.into_iter().map(Into::into).collect(),
        acceptor: acceptor.into(),
        ast: ast.into(),
    }
}

/// An aromatic-system overlay over `atoms`, carrying `ast` (an `AromaticSystemForm` — e.g. from
/// `from_electrons` — or a DSL spec string).
pub fn aromatic_system(
    atoms: impl IntoIterator<Item = impl Into<AtomArg>>,
    ast: impl Into<AromaticSystemForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::AromaticSystem {
        atoms: atoms.into_iter().map(Into::into).collect(),
        ast: ast.into(),
    }
}

/// A multicenter-bond overlay over `atoms`, carrying `ast` (a `MulticenterBondForm` or a DSL spec
/// string).
pub fn multicenter_bond(
    atoms: impl IntoIterator<Item = impl Into<AtomArg>>,
    ast: impl Into<MulticenterBondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::MulticenterBond {
        atoms: atoms.into_iter().map(Into::into).collect(),
        ast: ast.into(),
    }
}

/// A noncovalent-bond overlay between `first` and `second`, carrying `ast` (a `NoncovalentBondForm` —
/// e.g. from `from_kind` — or a DSL spec string).
pub fn noncovalent_bond(
    first: impl Into<AtomArg>,
    second: impl Into<AtomArg>,
    ast: impl Into<NoncovalentBondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::NoncovalentBond {
        first: first.into(),
        second: second.into(),
        ast: ast.into(),
    }
}

/// A stereo-atom overlay on `site`, with its ordered `ligands` and configuration `ast`.
pub fn stereo_atom(
    site: impl Into<AtomArg>,
    ligands: impl IntoIterator<Item = impl Into<StereoLigandArg>>,
    ast: impl Into<StereoAtomForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::StereoAtom {
        site: site.into(),
        ligands: ligands.into_iter().map(Into::into).collect(),
        ast: ast.into(),
    }
}

/// A stereo-bond overlay on the named bond `site`, with its ordered `ligands` and configuration `ast`.
pub fn stereo_bond(
    site: impl Into<BondArg>,
    ligands: impl IntoIterator<Item = impl Into<StereoLigandArg>>,
    ast: impl Into<StereoBondForm>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::StereoBond {
        site: site.into(),
        ligands: ligands.into_iter().map(Into::into).collect(),
        ast: ast.into(),
    }
}

/// Fill every unspecified atom field with its ground default at `build`.
pub fn ground() -> MoleculeSpecTerm {
    MoleculeSpecTerm::Ground
}

/// A declarative molecule specification — a sum of [`MoleculeSpecTerm`]s composed with
/// `+`, materialized by [`build`](MoleculeSpec::build). Terms lower onto `MoleculeBuilder`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeSpec {
    terms: Vec<MoleculeSpecTerm>,
}

impl MoleculeSpec {
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    /// Materialize the spec. Atoms are created in `+`-order (fixing their positions),
    /// relations are wired against the created atoms, and — if any `ground` term is
    /// present — every unspecified atom field is filled with its ground default.
    pub fn build(self) -> Molecule {
        let grounded = self
            .terms
            .iter()
            .any(|term| matches!(term, MoleculeSpecTerm::Ground));
        let mut cx = BuildContext::new(grounded);

        for term in self.terms {
            match term {
                MoleculeSpecTerm::Ground => {}
                MoleculeSpecTerm::Atoms(specs) => {
                    for spec in specs {
                        cx.resolve(spec);
                    }
                }
                MoleculeSpecTerm::Bond {
                    name,
                    first,
                    second,
                    ast,
                } => {
                    let first = cx.resolve(first);
                    let second = cx.resolve(second);
                    let id = cx.builder.bond(first, second, ast);
                    if let Some(name) = name {
                        cx.bond_names.insert(name, (id, first, second));
                    }
                }
                MoleculeSpecTerm::Chain(specs) => {
                    let ids = cx.resolve_all(specs);
                    cx.builder.chain(ids);
                }
                MoleculeSpecTerm::Ring(specs) => {
                    let ids = cx.resolve_all(specs);
                    cx.builder.ring(ids);
                }
                MoleculeSpecTerm::DativeBond {
                    donors,
                    acceptor,
                    ast,
                } => {
                    let donors = cx.resolve_all(donors);
                    let acceptor = cx.resolve(acceptor);
                    cx.builder.dative_bond(donors, acceptor, ast);
                }
                MoleculeSpecTerm::AromaticSystem { atoms, ast } => {
                    let ids = cx.resolve_all(atoms);
                    cx.builder.aromatic_system(ids, ast);
                }
                MoleculeSpecTerm::MulticenterBond { atoms, ast } => {
                    let ids = cx.resolve_all(atoms);
                    cx.builder.multicenter_bond(ids, ast);
                }
                MoleculeSpecTerm::NoncovalentBond { first, second, ast } => {
                    let first = cx.resolve(first);
                    let second = cx.resolve(second);
                    cx.builder.noncovalent_bond(first, second, ast);
                }
                MoleculeSpecTerm::StereoAtom { site, ligands, ast } => {
                    let site = cx.resolve(site);
                    let ligands = cx.resolve_stereo_ligands(ligands, site, site);
                    cx.builder.stereo_atom(site, ligands, ast);
                }
                MoleculeSpecTerm::StereoBond { site, ligands, ast } => {
                    let (bond, first, second) = cx.resolve_bond(site);
                    let ligands = cx.resolve_stereo_ligands(ligands, first, second);
                    cx.builder.stereo_bond(bond, ligands, ast);
                }
            }
        }
        cx.builder.build()
    }
}

impl Add<MoleculeSpecTerm> for MoleculeSpec {
    type Output = MoleculeSpec;
    fn add(mut self, term: MoleculeSpecTerm) -> Self {
        self.terms.push(term);
        self
    }
}

impl Add<MoleculeSpecTerm> for MoleculeSpecTerm {
    type Output = MoleculeSpec;
    fn add(self, other: MoleculeSpecTerm) -> MoleculeSpec {
        MoleculeSpec {
            terms: vec![self, other],
        }
    }
}

/// Materialization state: the target builder plus the position/name resolution maps. `bond_names`
/// records each named bond's id and (written-order) endpoints for a later stereo-bond site.
struct BuildContext {
    builder: MoleculeBuilder,
    positions: Vec<AtomId>,
    names: HashMap<String, AtomId>,
    bond_names: HashMap<String, (BondId, AtomId, AtomId)>,
}

impl BuildContext {
    fn new(grounded: bool) -> Self {
        Self {
            builder: if grounded {
                MoleculeBuilder::ground()
            } else {
                MoleculeBuilder::new()
            },
            positions: Vec::new(),
            names: HashMap::new(),
            bond_names: HashMap::new(),
        }
    }

    /// Resolve one atom argument: create a fresh atom (recording its position and any
    /// name) or look one up. Referencing an atom not yet introduced panics — a wrong
    /// position/name is a construction bug, like an out-of-bounds index.
    fn resolve(&mut self, arg: AtomArg) -> AtomId {
        match arg {
            AtomArg::New { spec, name } => {
                let id = self.builder.atom(spec);
                self.positions.push(id);
                if let Some(name) = name {
                    self.names.insert(name, id);
                }
                id
            }
            AtomArg::Index(position) => {
                *self.positions.get(position as usize).unwrap_or_else(|| {
                    panic!("spec references atom position {position} before it is introduced")
                })
            }
            AtomArg::Name(name) => *self
                .names
                .get(&name)
                .unwrap_or_else(|| panic!("spec references unknown atom name {name:?}")),
        }
    }

    fn resolve_all(&mut self, args: Vec<AtomArg>) -> Vec<AtomId> {
        args.into_iter().map(|arg| self.resolve(arg)).collect()
    }

    /// Resolve a bond reference to its id and (written-order) endpoints. Referencing an unknown bond
    /// name panics — a construction bug, like a bad atom ref.
    fn resolve_bond(&self, arg: BondArg) -> (BondId, AtomId, AtomId) {
        match arg {
            BondArg::Name(name) => *self
                .bond_names
                .get(&name)
                .unwrap_or_else(|| panic!("spec references unknown bond name {name:?}")),
        }
    }

    /// Lower stereo ligand args to `StereoLigand`s. An atom ligand carries its substituent; a virtual
    /// ligand bears on `first_atom` for positions 0–1 and `second_atom` for 2+ (both are the site for a
    /// stereo atom; the bond's two atoms for a stereo bond).
    fn resolve_stereo_ligands(
        &mut self,
        ligands: Vec<StereoLigandArg>,
        first_atom: AtomId,
        second_atom: AtomId,
    ) -> Vec<StereoLigand> {
        ligands
            .into_iter()
            .enumerate()
            .map(|(index, ligand)| {
                let bearing = if index < 2 { first_atom } else { second_atom };
                match ligand {
                    StereoLigandArg::Atom(arg) => {
                        StereoLigand::new(self.resolve(arg), StereoLigandKind::Atom)
                    }
                    StereoLigandArg::ImplicitHydrogen => {
                        StereoLigand::new(bearing, StereoLigandKind::ImplicitHydrogen)
                    }
                    StereoLigandArg::LonePair => {
                        StereoLigand::new(bearing, StereoLigandKind::LonePair)
                    }
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::ir::id::{
        AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
        StereoBondId,
    };
    use crate::ir::noncovalent::NoncovalentBondKind;
    use crate::ir::stereo::{StereoCoset, StereoKind};
    use crate::ir::value::NumForm;

    #[rstest]
    #[case::element(
        AtomArg::from(Element::C),
        AtomArg::New { spec: AtomForm::from_element(Element::C), name: None }
    )]
    #[case::position(AtomArg::from(5_u32), AtomArg::Index(5))]
    #[case::position_signed(AtomArg::from(3_i32), AtomArg::Index(3))]
    #[case::named_tuple(
        AtomArg::from(("carbonyl", Element::O)),
        AtomArg::New { spec: AtomForm::from_element(Element::O), name: Some("carbonyl".to_string()) }
    )]
    #[case::name_fn(name("amide"), AtomArg::Name("amide".to_string()))]
    fn test_atom_arg_from(#[case] arg: AtomArg, #[case] expected: AtomArg) {
        assert_eq!(arg, expected);
    }

    #[rstest]
    #[should_panic(expected = "atom position must be non-negative")]
    fn test_atom_arg_from_error() {
        let _ = AtomArg::from(-1_i32);
    }

    #[rstest]
    #[case::str(BondArg::from("ring"), BondArg::Name("ring".to_string()))]
    #[case::string(BondArg::from("ring".to_string()), BondArg::Name("ring".to_string()))]
    fn test_bond_arg_from(#[case] arg: BondArg, #[case] expected: BondArg) {
        assert_eq!(arg, expected);
    }

    #[rstest]
    fn test_stereo_ligand_arg_from() {
        assert_eq!(
            StereoLigandArg::from(AtomArg::from(2_u32)),
            StereoLigandArg::Atom(AtomArg::Index(2))
        );
    }

    #[rstest]
    fn test_molecule_spec_build() {
        // named create (tuple) + name reference + explicit double bond
        let spec = atoms([("c", Element::C), ("o", Element::O)]) + double(name("c"), name("o"));
        let mol = spec.build();

        assert_eq!(mol.atoms().count(), 2);
        assert_eq!(mol.bonds().count(), 1);
        assert_eq!(mol.bond(BondId(0)).ast, &BondForm::from_order(2));
        assert_eq!(mol.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    }

    #[rstest]
    fn test_molecule_spec_build_position_reference() {
        // bare integer literals resolve to positions via From<i32>
        let spec = atoms([Element::C, Element::O]) + single(0, 1);
        let mol = spec.build();

        assert_eq!(mol.bond(BondId(0)).ast, &BondForm::from_order(1));
        assert_eq!(mol.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    }

    #[rstest]
    fn test_molecule_spec_build_create_and_wire() {
        // specs in a relation term create the atoms then wire them
        let spec = MoleculeSpec::new() + single(Element::C, Element::O);
        let mol = spec.build();

        assert_eq!(mol.atoms().count(), 2);
        assert_eq!(mol.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    }

    #[rstest]
    #[case::single(single(0_u32, 1_u32), BondForm::from_order(1))]
    #[case::double(double(0_u32, 1_u32), BondForm::from_order(2))]
    #[case::triple(triple(0_u32, 1_u32), BondForm::from_order(3))]
    #[case::aromatic(
        aromatic_bond(0_u32, 1_u32),
        BondForm::from_order(1).with_constraint(BondConstraintForm::aromatic(true))
    )]
    #[case::explicit(
        bond(0_u32, 1_u32, BondForm::from_order(2).with_charge(-1_i64)),
        BondForm::from_order(2).with_charge(-1_i64)
    )]
    fn test_molecule_spec_bond_terms(
        #[case] bond_term: MoleculeSpecTerm,
        #[case] expected: BondForm,
    ) {
        let spec = atoms([Element::C, Element::C]) + bond_term;
        let mol = spec.build();

        assert_eq!(mol.bond(BondId(0)).ast, &expected);
    }

    #[rstest]
    fn test_molecule_spec_chain() {
        let spec = MoleculeSpec::new() + chain([Element::C, Element::C, Element::O]);
        let mol = spec.build();

        assert_eq!(mol.bonds().count(), 2);
        assert_eq!(mol.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
        assert_eq!(mol.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_molecule_spec_ring() {
        let spec = MoleculeSpec::new() + ring([Element::C, Element::C, Element::C]);
        let mol = spec.build();

        assert_eq!(mol.bonds().count(), 3);
        // closing bond (atom 2 → atom 0), stored normalized
        assert_eq!(mol.bond(BondId(2)).atom_ids(), [AtomId(0), AtomId(2)]);
    }

    #[rstest]
    fn test_molecule_spec_dative_bond() {
        let spec = atoms([Element::N, Element::B])
            + dative_bond([0_u32], 1_u32, DativeBondForm::default());
        let mol = spec.build();

        assert_eq!(mol.dative_bonds().count(), 1);
        assert_eq!(
            mol.dative_bond(DativeBondId(0)).ast,
            &DativeBondForm::default()
        );
    }

    #[rstest]
    fn test_molecule_spec_aromatic_system() {
        let spec = atoms([Element::C, Element::C])
            + aromatic_system(
                [0_u32, 1_u32],
                AromaticSystemForm::from_electrons(vec![1, 1]),
            );
        let mol = spec.build();

        assert_eq!(
            mol.aromatic_system(AromaticSystemId(0)).ast,
            &AromaticSystemForm::from_electrons(vec![1, 1])
        );
    }

    #[rstest]
    fn test_molecule_spec_multicenter_bond() {
        let spec = atoms([Element::B, Element::B, Element::H])
            + multicenter_bond(
                [0_u32, 1_u32, 2_u32],
                MulticenterBondForm::from_electrons(vec![1, 1, 1]),
            );
        let mol = spec.build();

        assert_eq!(
            mol.multicenter_bond(MulticenterBondId(0)).ast,
            &MulticenterBondForm::from_electrons(vec![1, 1, 1])
        );
    }

    #[rstest]
    fn test_molecule_spec_noncovalent_bond() {
        let spec = atoms([Element::O, Element::H])
            + noncovalent_bond(
                0_u32,
                1_u32,
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            );
        let mol = spec.build();

        assert_eq!(
            mol.noncovalent_bond(NoncovalentBondId(0)).ast,
            &NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
        );
    }

    #[rstest]
    fn test_molecule_spec_named_bond() {
        // a named bond carries its spec and binds a label (inert without a reference)
        let spec = atoms([Element::C, Element::C])
            + named_bond("db", 0_u32, 1_u32, BondForm::from_order(2));
        let mol = spec.build();

        assert_eq!(mol.bond(BondId(0)).ast, &BondForm::from_order(2));
    }

    #[rstest]
    fn test_molecule_spec_stereo_atom() {
        // a tetrahedral center: three atom ligands plus an implicit hydrogen borne by the site
        let spec = atoms([Element::C, Element::F, Element::Cl, Element::Br])
            + stereo_atom(
                0_u32,
                [
                    StereoLigandArg::Atom(1_u32.into()),
                    StereoLigandArg::Atom(2_u32.into()),
                    StereoLigandArg::Atom(3_u32.into()),
                    StereoLigandArg::ImplicitHydrogen,
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            );
        let mol = spec.build();

        assert_eq!(mol.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
        // the implicit-H bears on the site (atom 0), so atom incidence is site + 3 substituents
        assert_eq!(
            mol.stereo_atom(StereoAtomId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]
        );
        assert_eq!(
            mol.stereo_atom(StereoAtomId(0)).ast,
            &StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
        );
    }

    #[rstest]
    fn test_molecule_spec_stereo_bond() {
        // stereo bond referencing the named double bond by name; virtual Hs bear on the two bond
        // atoms by ligand position (0–1 → first atom, 2–3 → second)
        let spec = atoms([Element::C, Element::C, Element::F, Element::Cl])
            + named_bond("db", 0_u32, 1_u32, BondForm::from_order(2))
            + stereo_bond(
                "db",
                [
                    StereoLigandArg::Atom(2_u32.into()),
                    StereoLigandArg::ImplicitHydrogen,
                    StereoLigandArg::Atom(3_u32.into()),
                    StereoLigandArg::ImplicitHydrogen,
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            );
        let mol = spec.build();

        assert_eq!(mol.stereo_bond(StereoBondId(0)).site_id(), BondId(0));
        assert_eq!(
            mol.stereo_bond(StereoBondId(0))
                .implicit_hydrogen_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1)]
        );
        assert_eq!(
            mol.stereo_bond(StereoBondId(0)).ast,
            &StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))
        );
    }

    #[rstest]
    #[case::grounded(atom(Element::C) + ground(), NumForm::Lit(0))]
    #[case::ungrounded(MoleculeSpec::new() + atom(Element::C), NumForm::Undetermined)]
    fn test_molecule_spec_ground(#[case] spec: MoleculeSpec, #[case] expected_charge: NumForm) {
        let mol = spec.build();

        assert_eq!(mol.atom(AtomId(0)).ast.charge, expected_charge);
    }
}
