//! Molecule DSL.
//!
//! `MoleculeDsl` wraps a `MoleculeAst` together with the `Metadata` that records
//! the surface-form id/alias bindings (atom ids, bond ids, etc.). The EDN
//! form is a map keyed by `:atoms`, `:bonds`, `:dative`, `:aromatic`,
//! `:multicenter`, `:noncovalent`, `:atom-aliases`/`:aliases`, and
//! `:constraints`. Each entity delegates to its own entity DSL. Constraints
//! parse directly into the typed `Constraint` tree.

use std::fmt::{self, Display};
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, FromEdn, ToEdn};

use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::dative::DativeBondDsl;
use super::error::ParseError;
use super::noncovalent::NoncovalentBondDsl;
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::config::MoleculeAstConfig;
use crate::ast::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, DativeBondConstraint,
    MulticenterBondConstraint, NoncovalentBondConstraint,
};
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::spin::SpinStateAst;
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::ValueAst;

/// Surface-form metadata paired with a `MoleculeAst`. Records atom ids,
/// per-entity ids, and the atom-alias table. Never drifts: rewrapped
/// atomically through `MoleculeDsl::from_parts`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_ids: IndexMap<AtomIdx, String>,
    pub atom_aliases: BiMap<String, AtomDsl>,
    pub bond_ids: IndexMap<BondIdx, String>,
    pub dative_bond_ids: IndexMap<DativeBondIdx, String>,
    pub aromatic_system_ids: IndexMap<AromaticSystemIdx, String>,
    pub multicenter_bond_ids: IndexMap<MulticenterBondIdx, String>,
    pub noncovalent_bond_ids: IndexMap<NoncovalentBondIdx, String>,
}

/// Surface DSL for a whole molecule. Pairs `MoleculeAst` with `Metadata`;
/// fields are private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    ast: MoleculeAst,
    metadata: Metadata,
}

impl MoleculeDsl {
    pub fn from_parts(ast: MoleculeAst, metadata: Metadata) -> Self {
        Self { ast, metadata }
    }

    pub fn from_ast(ast: MoleculeAst) -> Self {
        Self {
            ast,
            metadata: Metadata::default(),
        }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (MoleculeAst, Metadata) {
        (self.ast, self.metadata)
    }
}

impl PartialEq for MoleculeDsl {
    fn eq(&self, other: &Self) -> bool {
        self.ast == other.ast && self.metadata == other.metadata
    }
}

impl Eq for MoleculeDsl {}

impl FromStr for MoleculeDsl {
    type Err = ParseError;

    fn from_str(_s: &str) -> Result<Self, Self::Err> {
        todo!("MoleculeDsl::from_str (Phase 5)")
    }
}

impl Display for MoleculeDsl {
    fn fmt(&self, _f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!("MoleculeDsl::fmt (Phase 2)")
    }
}

impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(_edn: &Edn<'de>) -> Result<Self, DeError> {
        todo!("MoleculeDsl::from_edn (Phase 3)")
    }

    fn from_edn_str(_input: &'de str) -> Result<Self, EdnError> {
        todo!("MoleculeDsl::from_edn_str (Phase 4)")
    }
}

impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        todo!("MoleculeDsl::to_edn (Phase 2)")
    }
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    type Error = ParseError;

    fn from_ast(_ast: &MoleculeAst, _cfg: &MoleculeAstConfig) -> Result<Self, ParseError> {
        todo!("MoleculeDsl::from_ast (Phase 5)")
    }
}

impl ToAst<MoleculeAst> for MoleculeDsl {
    type Error = ParseError;

    fn to_ast(&self, _cfg: &MoleculeAstConfig) -> Result<MoleculeAst, ParseError> {
        todo!("MoleculeDsl::to_ast (Phase 5)")
    }
}

// -- Private parse intermediate ---------------------------------------------
//
// Unresolved, owned-by-value tree that mirrors the EDN shape. Atom entries and
// per-bond endpoints carry `AtomRef` (index or id); constraint leaves carry
// typed per-entity `Constraint*` variants already parsed from their EDN form.
// Lowered destructively via `into_ast(self, cfg)` so that allocations move
// into the final `MoleculeAst`.

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct MoleculeInput {
    pub(crate) atoms: Vec<AtomEntryInput>,
    pub(crate) bonds: Vec<BondEntryInput>,
    pub(crate) dative_bonds: Vec<DativeBondEntryInput>,
    pub(crate) aromatic_systems: Vec<AromaticSystemEntryInput>,
    pub(crate) multicenter_bonds: Vec<MulticenterBondEntryInput>,
    pub(crate) noncovalent_bonds: Vec<NoncovalentBondEntryInput>,
    pub(crate) atom_aliases: Vec<(String, AtomDsl)>,
    pub(crate) constraints: Vec<ConstraintInput>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomEntryInput {
    Bare(AtomDsl),
    Alias(String),
    WithId {
        id: String,
        inner: Box<AtomEntryInput>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) donor: AtomRef,
    pub(crate) acceptor: AtomRef,
    pub(crate) bond: DativeBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AromaticSystemEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) system: AromaticSystemAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct MulticenterBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) bond: MulticenterBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct NoncovalentBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: NoncovalentBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomRef {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum BondRef {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum DativeBondRef {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AromaticSystemRef {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum MulticenterBondRef {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum NoncovalentBondRef {
    Index(usize),
    Id(String),
}

/// Pre-resolution mirror of `Constraint`. Entity leaves carry the unresolved
/// reference; combinators recurse; molecule-scope leaves are parsed in place.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum ConstraintInput {
    Atom(AtomRef, AtomConstraint),
    Bond(BondRef, BondConstraint),
    DativeBond(DativeBondRef, DativeBondConstraint),
    AromaticSystem(AromaticSystemRef, AromaticSystemConstraint),
    MulticenterBond(MulticenterBondRef, MulticenterBondConstraint),
    NoncovalentBond(NoncovalentBondRef, NoncovalentBondConstraint),
    Molecule(MoleculeConstraintInput),
    And(Vec<ConstraintInput>),
    Or(Vec<ConstraintInput>),
    Not(Box<ConstraintInput>),
}

/// Pre-resolution mirror of `MoleculeConstraint`. `SubPattern.pattern` holds a
/// nested `MoleculeInput` so the sub-molecule is resolved against its own
/// local id scope before the outer molecule completes lowering.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum MoleculeConstraintInput {
    ChargeSum {
        atoms: Vec<AtomRef>,
        sum: ValueAst,
    },
    SpinSum {
        atoms: Vec<AtomRef>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Vec<BondRef>,
        sum: ValueAst,
    },
    Connected(Vec<AtomRef>),
    SubPattern {
        anchor: SubPatternAnchorInput,
        pattern: Box<MoleculeInput>,
    },
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct SubPatternAnchorInput {
    pub(crate) atoms: Vec<(AtomRef, AtomRef)>,
    pub(crate) bonds: Vec<(BondRef, BondRef)>,
    pub(crate) dative_bonds: Vec<(DativeBondRef, DativeBondRef)>,
    pub(crate) aromatic_systems: Vec<(AromaticSystemRef, AromaticSystemRef)>,
    pub(crate) multicenter_bonds: Vec<(MulticenterBondRef, MulticenterBondRef)>,
    pub(crate) noncovalent_bonds: Vec<(NoncovalentBondRef, NoncovalentBondRef)>,
}

impl MoleculeInput {
    /// Destructive lowering: consumes the input, resolves refs against the
    /// built id scopes, lifts bare entity-leaf constraints onto their entity
    /// AST's inline store, and produces the final `MoleculeAst` with its
    /// `Metadata`. Called from `FromEdn::from_edn` and the streaming path.
    #[allow(dead_code)]
    pub(crate) fn into_ast(
        self,
        _cfg: &MoleculeAstConfig,
    ) -> Result<(MoleculeAst, Metadata), ParseError> {
        todo!("MoleculeInput::into_ast (Phase 3)")
    }
}

