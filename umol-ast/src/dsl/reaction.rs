//! Reaction DSL.
//!
//! `ReactionDsl` wraps a `ReactionAst` together with the `ReactionMetadata` that records
//! the surface-form bindings (the lhs molecule metadata plus created-entity id ↔ name and
//! atom-alias bindings). The EDN form is a map keyed by `:lhs` (a molecule map) and
//! `:deltas` (a vector of `:add` / `:remove` / `:modify` / `:constraint` operations). Each
//! entity delegates to its own entity DSL.

use bimap::BiBTreeMap;
use indexmap::IndexMap;

use super::atom::AtomDsl;
use super::constraint::ConstraintDsl;
use super::molecule::{MoleculeInput, MoleculeMetadata};
use super::refs::{AtomRefDsl, BondRefDsl};
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::id::{AtomId, BondId};

/// Surface-form metadata paired with a `ReactionAst`: the lhs molecule metadata plus the
/// created-entity id bindings and atom aliases introduced by the deltas. Mirrors
/// `MoleculeMetadata` for the atom/bond entities (the reaction admits the `[:C "C#h3"]`
/// alias notation for added atoms).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    atom_ids: IndexMap<AtomId, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondId, String>,
}

impl ReactionMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn lhs(&self) -> &MoleculeMetadata {
        &self.lhs
    }

    pub fn atom_id(&self, id: AtomId) -> Option<&str> {
        self.atom_ids.get(&id).map(String::as_str)
    }

    pub fn bond_id(&self, id: BondId) -> Option<&str> {
        self.bond_ids.get(&id).map(String::as_str)
    }

    /// Name of the alias bound to this atom DSL, if any.
    pub fn atom_alias_for(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn has_atom_alias(&self, name: &str) -> bool {
        self.atom_aliases.contains_left(name)
    }

    pub fn has_atom_aliases(&self) -> bool {
        !self.atom_aliases.is_empty()
    }

    pub fn atom_aliases_len(&self) -> usize {
        self.atom_aliases.len()
    }

    pub fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_atom_id(&mut self, id: AtomId, name: impl Into<String>) {
        self.atom_ids.insert(id, name.into());
    }

    pub fn set_bond_id(&mut self, id: BondId, name: impl Into<String>) {
        self.bond_ids.insert(id, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases.insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_id(mut self, id: AtomId, name: impl Into<String>) -> Self {
        self.set_atom_id(id, name);
        self
    }

    pub fn with_bond_id(mut self, id: BondId, name: impl Into<String>) -> Self {
        self.set_bond_id(id, name);
        self
    }

    pub fn with_atom_alias(mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

/// One unresolved delta parsed from a `:deltas` entry. Refs stay symbolic
/// (`AtomRefDsl`/`BondRefDsl`); a `:modify` RHS is a partial entity AST
/// (unspecified fields `Undetermined`); `id` is the created entity's optional
/// `:id` name.
#[derive(Debug)]
pub(crate) enum DeltaInput {
    AtomAdd {
        id: Option<String>,
        value: AtomAst,
    },
    AtomRemove(AtomRefDsl),
    AtomModify(AtomRefDsl, AtomAst),
    BondAdd {
        id: Option<String>,
        atoms: [AtomRefDsl; 2],
        value: BondAst,
    },
    BondRemove(BondRefDsl),
    BondModify(BondRefDsl, BondAst),
    ConstraintAdd(ConstraintDsl),
    ConstraintRemove(ConstraintDsl),
}

/// Raw parse target for a reaction: the lhs molecule input plus the unresolved
/// deltas. Resolution (`into_ast`, R7) lifts this to `(ReactionAst, ReactionMetadata)`.
#[derive(Debug)]
pub(crate) struct ReactionInput {
    lhs: MoleculeInput,
    deltas: Vec<DeltaInput>,
}
