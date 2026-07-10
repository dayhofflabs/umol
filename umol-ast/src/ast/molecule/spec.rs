//! Declarative `+`-spec construction — build a `MoleculeAst` by summing free-function
//! *terms* into a `MoleculeSpec`, then materializing. Each term lowers onto the L1
//! `MoleculeBuilder`; every atom slot is an [`AtomArg`] that either creates a fresh atom
//! (optionally named) or references one already introduced (by position or by name).

use std::collections::HashMap;
use std::ops::Add;

use umol_chem::element::Element;

use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::BondConstraintAst;
use super::super::id::AtomId;
use super::super::noncovalent::NoncovalentBondKind;
use super::{MoleculeAst, MoleculeBuilder};

/// An atom argument to a spec term: create a fresh atom (optionally named) or reference
/// one already introduced — by creation `position` or by `name`. What you write picks the
/// variant by type: `C`/`"C#h3"`/`AtomAst` → create, `(name, spec)` tuple → create-named,
/// a bare integer → by position, [`name`] → by name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomArg {
    New { spec: AtomAst, name: Option<String> },
    ByPosition(u32),
    ByName(String),
}

impl From<Element> for AtomArg {
    fn from(element: Element) -> Self {
        Self::New {
            spec: AtomAst::from_element(element),
            name: None,
        }
    }
}

impl From<&str> for AtomArg {
    fn from(spec: &str) -> Self {
        Self::New {
            spec: AtomAst::from(spec),
            name: None,
        }
    }
}

impl From<AtomAst> for AtomArg {
    fn from(spec: AtomAst) -> Self {
        Self::New { spec, name: None }
    }
}

impl From<u32> for AtomArg {
    fn from(position: u32) -> Self {
        Self::ByPosition(position)
    }
}

impl From<i32> for AtomArg {
    fn from(position: i32) -> Self {
        Self::ByPosition(
            u32::try_from(position)
                .unwrap_or_else(|_| panic!("atom position must be non-negative, got {position}")),
        )
    }
}

impl<S: Into<String>, T: Into<AtomAst>> From<(S, T)> for AtomArg {
    fn from((name, spec): (S, T)) -> Self {
        Self::New {
            spec: spec.into(),
            name: Some(name.into()),
        }
    }
}

/// Reference an already-introduced atom by name.
pub fn name(name: impl Into<String>) -> AtomArg {
    AtomArg::ByName(name.into())
}

/// A single term of a [`MoleculeSpec`]: an introduction, a relation, or a default.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeSpecTerm {
    Atoms(Vec<AtomArg>),
    Bond {
        first: AtomArg,
        second: AtomArg,
        ast: BondAst,
    },
    Chain(Vec<AtomArg>),
    Ring(Vec<AtomArg>),
    DativeBond {
        donors: Vec<AtomArg>,
        acceptor: AtomArg,
    },
    AromaticSystem {
        atoms: Vec<AtomArg>,
        electrons: Vec<i64>,
    },
    MulticenterBond {
        atoms: Vec<AtomArg>,
        electrons: Vec<i64>,
    },
    NoncovalentBond {
        first: AtomArg,
        second: AtomArg,
        kind: NoncovalentBondKind,
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
        first: first.into(),
        second: second.into(),
        ast: BondAst::from_order(1),
    }
}

/// A double (order-2) bond.
pub fn double(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        first: first.into(),
        second: second.into(),
        ast: BondAst::from_order(2),
    }
}

/// A triple (order-3) bond.
pub fn triple(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        first: first.into(),
        second: second.into(),
        ast: BondAst::from_order(3),
    }
}

/// A bond carrying an explicit `BondAst` — the escape hatch for a charge, spin, order-set,
/// or constraint the order verbs can't set.
pub fn bond(
    first: impl Into<AtomArg>,
    second: impl Into<AtomArg>,
    ast: impl Into<BondAst>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        first: first.into(),
        second: second.into(),
        ast: ast.into(),
    }
}

/// An aromatic bond — order 1 with the aromatic flag (`1#a`); resolution perceives the
/// aromatic system. Not exclusive with [`aromatic_system`].
pub fn aromatic_bond(first: impl Into<AtomArg>, second: impl Into<AtomArg>) -> MoleculeSpecTerm {
    MoleculeSpecTerm::Bond {
        first: first.into(),
        second: second.into(),
        ast: BondAst::from_order(1).with_constraint(BondConstraintAst::aromatic(true)),
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

/// A dative bond from `donors` to `acceptor`.
pub fn dative_bond(
    donors: impl IntoIterator<Item = impl Into<AtomArg>>,
    acceptor: impl Into<AtomArg>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::DativeBond {
        donors: donors.into_iter().map(Into::into).collect(),
        acceptor: acceptor.into(),
    }
}

/// An aromatic-system overlay over `atoms`, one π-`electrons` count per atom.
pub fn aromatic_system(
    atoms: impl IntoIterator<Item = impl Into<AtomArg>>,
    electrons: impl IntoIterator<Item = i64>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::AromaticSystem {
        atoms: atoms.into_iter().map(Into::into).collect(),
        electrons: electrons.into_iter().collect(),
    }
}

/// A multicenter-bond overlay over `atoms`, one `electrons` count per atom.
pub fn multicenter_bond(
    atoms: impl IntoIterator<Item = impl Into<AtomArg>>,
    electrons: impl IntoIterator<Item = i64>,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::MulticenterBond {
        atoms: atoms.into_iter().map(Into::into).collect(),
        electrons: electrons.into_iter().collect(),
    }
}

/// A noncovalent-bond overlay of `kind` between `first` and `second`.
pub fn noncovalent_bond(
    first: impl Into<AtomArg>,
    second: impl Into<AtomArg>,
    kind: NoncovalentBondKind,
) -> MoleculeSpecTerm {
    MoleculeSpecTerm::NoncovalentBond {
        first: first.into(),
        second: second.into(),
        kind,
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
    pub fn build(self) -> MoleculeAst {
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
                MoleculeSpecTerm::Bond { first, second, ast } => {
                    let first = cx.resolve(first);
                    let second = cx.resolve(second);
                    cx.builder.bond(first, second, ast);
                }
                MoleculeSpecTerm::Chain(specs) => {
                    let ids = cx.resolve_all(specs);
                    cx.builder.chain(ids);
                }
                MoleculeSpecTerm::Ring(specs) => {
                    let ids = cx.resolve_all(specs);
                    cx.builder.ring(ids);
                }
                MoleculeSpecTerm::DativeBond { donors, acceptor } => {
                    let donors = cx.resolve_all(donors);
                    let acceptor = cx.resolve(acceptor);
                    cx.builder.dative_bond(donors, acceptor);
                }
                MoleculeSpecTerm::AromaticSystem { atoms, electrons } => {
                    let ids = cx.resolve_all(atoms);
                    cx.builder.aromatic_system(ids, electrons);
                }
                MoleculeSpecTerm::MulticenterBond { atoms, electrons } => {
                    let ids = cx.resolve_all(atoms);
                    cx.builder.multicenter_bond(ids, electrons);
                }
                MoleculeSpecTerm::NoncovalentBond {
                    first,
                    second,
                    kind,
                } => {
                    let first = cx.resolve(first);
                    let second = cx.resolve(second);
                    cx.builder.noncovalent_bond(first, second, kind);
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

/// Materialization state: the target builder plus the position/name resolution maps.
struct BuildContext {
    builder: MoleculeBuilder,
    positions: Vec<AtomId>,
    names: HashMap<String, AtomId>,
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
            AtomArg::ByPosition(position) => *self.positions.get(position as usize).unwrap_or_else(
                || panic!("spec references atom position {position} before it is introduced"),
            ),
            AtomArg::ByName(name) => *self
                .names
                .get(&name)
                .unwrap_or_else(|| panic!("spec references unknown atom name {name:?}")),
        }
    }

    fn resolve_all(&mut self, args: Vec<AtomArg>) -> Vec<AtomId> {
        args.into_iter().map(|arg| self.resolve(arg)).collect()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::id::{
        AromaticSystemId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    };
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::NoncovalentBondAst;
    use crate::ast::value::ValueAst;

    #[rstest]
    #[case::element(
        AtomArg::from(Element::C),
        AtomArg::New { spec: AtomAst::from_element(Element::C), name: None }
    )]
    #[case::position(AtomArg::from(5_u32), AtomArg::ByPosition(5))]
    #[case::position_signed(AtomArg::from(3_i32), AtomArg::ByPosition(3))]
    #[case::named_tuple(
        AtomArg::from(("carbonyl", Element::O)),
        AtomArg::New { spec: AtomAst::from_element(Element::O), name: Some("carbonyl".to_string()) }
    )]
    #[case::name_fn(name("amide"), AtomArg::ByName("amide".to_string()))]
    fn test_atom_arg_from(#[case] arg: AtomArg, #[case] expected: AtomArg) {
        assert_eq!(arg, expected);
    }

    #[rstest]
    #[should_panic(expected = "atom position must be non-negative")]
    fn test_atom_arg_from_error() {
        let _ = AtomArg::from(-1_i32);
    }

    #[rstest]
    fn test_molecule_spec_build() {
        // named create (tuple) + name reference + explicit double bond
        let spec = atoms([("c", Element::C), ("o", Element::O)]) + double(name("c"), name("o"));
        let mol = spec.build();

        assert_eq!(mol.atoms().count(), 2);
        assert_eq!(mol.bonds().count(), 1);
        assert_eq!(mol.bond(BondId(0)).ast, &BondAst::from_order(2));
        assert_eq!(mol.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    }

    #[rstest]
    fn test_molecule_spec_build_position_reference() {
        // bare integer literals resolve to positions via From<i32>
        let spec = atoms([Element::C, Element::O]) + single(0, 1);
        let mol = spec.build();

        assert_eq!(mol.bond(BondId(0)).ast, &BondAst::from_order(1));
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
    #[case::single(single(0_u32, 1_u32), BondAst::from_order(1))]
    #[case::double(double(0_u32, 1_u32), BondAst::from_order(2))]
    #[case::triple(triple(0_u32, 1_u32), BondAst::from_order(3))]
    #[case::aromatic(
        aromatic_bond(0_u32, 1_u32),
        BondAst::from_order(1).with_constraint(BondConstraintAst::aromatic(true))
    )]
    #[case::explicit(
        bond(0_u32, 1_u32, BondAst::from_order(2).with_charge(-1_i64)),
        BondAst::from_order(2).with_charge(-1_i64)
    )]
    fn test_molecule_spec_bond_terms(#[case] bond_term: MoleculeSpecTerm, #[case] expected: BondAst) {
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
        let spec = atoms([Element::N, Element::B]) + dative_bond([0_u32], 1_u32);
        let mol = spec.build();

        assert_eq!(mol.dative_bonds().count(), 1);
        assert_eq!(mol.dative_bond(DativeBondId(0)).ast, &DativeBondAst::default());
    }

    #[rstest]
    fn test_molecule_spec_aromatic_system() {
        let spec = atoms([Element::C, Element::C]) + aromatic_system([0_u32, 1_u32], [1_i64, 1_i64]);
        let mol = spec.build();

        assert_eq!(
            mol.aromatic_system(AromaticSystemId(0)).ast,
            &AromaticSystemAst::from_electrons(vec![1, 1])
        );
    }

    #[rstest]
    fn test_molecule_spec_multicenter_bond() {
        let spec = atoms([Element::B, Element::B, Element::H])
            + multicenter_bond([0_u32, 1_u32, 2_u32], [1_i64, 1_i64, 1_i64]);
        let mol = spec.build();

        assert_eq!(
            mol.multicenter_bond(MulticenterBondId(0)).ast,
            &MulticenterBondAst::from_electrons(vec![1, 1, 1])
        );
    }

    #[rstest]
    fn test_molecule_spec_noncovalent_bond() {
        let spec =
            atoms([Element::O, Element::H]) + noncovalent_bond(0_u32, 1_u32, NoncovalentBondKind::HydrogenBond);
        let mol = spec.build();

        assert_eq!(
            mol.noncovalent_bond(NoncovalentBondId(0)).ast,
            &NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
        );
    }

    #[rstest]
    #[case::grounded(atom(Element::C) + ground(), ValueAst::Lit(0))]
    #[case::ungrounded(MoleculeSpec::new() + atom(Element::C), ValueAst::Undetermined)]
    fn test_molecule_spec_ground(#[case] spec: MoleculeSpec, #[case] expected_charge: ValueAst) {
        let mol = spec.build();

        assert_eq!(mol.atom(AtomId(0)).ast.charge, expected_charge);
    }
}
