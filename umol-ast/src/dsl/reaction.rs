//! Reaction DSL.
//!
//! `ReactionDsl` wraps a `ReactionAst` together with the `ReactionMetadata` that records
//! the surface-form bindings (the lhs molecule metadata plus created-entity id ↔ name and
//! atom-alias bindings). The EDN form is a map keyed by `:lhs` (a molecule map) and
//! `:deltas` (a vector of `:add` / `:remove` / `:modify` / `:constraint` operations). Each
//! entity delegates to its own entity DSL.

use bimap::BiBTreeMap;
use indexmap::IndexMap;

use super::atom::{lower_atom, raise_atom, AtomDsl};
use super::bond::{lower_bond, raise_bond};
use super::config::{DeltaDefaults, ReactionDefaults};
use super::constraint::ConstraintDsl;
use super::molecule::{MoleculeDsl, MoleculeInput, MoleculeMetadata};
use super::refs::{AtomRefDsl, BondRefDsl};
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::delta::{AtomDelta, BondDelta, Delta};
use crate::ast::id::{AtomId, BondId};
use crate::ast::reaction::ReactionAst;
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL for a reaction. Pairs `ReactionAst` with `ReactionMetadata`; fields are
/// private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionDsl {
    ast: ReactionAst,
    metadata: ReactionMetadata,
}

impl ReactionDsl {
    pub fn from_parts(ast: ReactionAst, metadata: ReactionMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &ReactionAst {
        &self.ast
    }

    pub fn metadata(&self) -> &ReactionMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionAst, ReactionMetadata) {
        (self.ast, self.metadata)
    }
}

impl FromAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn from_ast(ast: &ReactionAst, cfg: &Self::Ctx) -> Self {
        let lhs = MoleculeDsl::from_ast(&ast.lhs, &cfg.molecule_defaults())
            .into_parts()
            .0;
        let delta_cfg = cfg.delta_defaults();
        let mut deltas = ast.deltas.clone();
        for delta in deltas.iter_mut() {
            lower_delta(delta, &delta_cfg);
        }
        ReactionDsl {
            ast: ReactionAst { lhs, deltas },
            metadata: ReactionMetadata::default(),
        }
    }
}

impl IntoAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> ReactionAst {
        let ReactionAst { lhs, mut deltas } = self.ast;
        let lhs = MoleculeDsl::from_parts(lhs, MoleculeMetadata::default())
            .into_ast(&cfg.molecule_defaults());
        let delta_cfg = cfg.delta_defaults();
        for delta in deltas.iter_mut() {
            raise_delta(delta, &delta_cfg);
        }
        ReactionAst { lhs, deltas }
    }
}

/// Lower a delta's embedded entity AST to DSL-display form (AST → DSL).
fn lower_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            lower_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            lower_bond(ast, &cfg.bond)
        }
        _ => {}
    }
}

/// Raise a delta's embedded entity AST from DSL-display form (DSL → AST).
fn raise_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            raise_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            raise_bond(ast, &cfg.bond)
        }
        _ => {}
    }
}

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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::constraint::Constraint;
    use crate::ast::delta::{ConstraintDelta, Deltas};
    use crate::ast::edit::AtomFieldChange;
    use crate::ast::value::ValueAst;
    use crate::mol;

    #[rstest]
    #[case::sn2(ReactionAst::new(
        mol!(r##"{:atoms ["C" "Br"] :bonds [[0 1 "1"]]}"##),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::O) }),
            Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(2)],
                ast: BondAst::from_order(1),
            }),
            Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(1),
                change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
            }),
            Delta::Constraint(ConstraintDelta::Add(Constraint::And(vec![]))),
        ]),
    ))]
    fn test_reaction_dsl_from_ast_roundtrip(#[case] reaction: ReactionAst) {
        let cfg = ReactionDefaults::ground();
        let dsl = ReactionDsl::from_parts(reaction, ReactionMetadata::default());
        let lowered = ReactionDsl::from_ast(&dsl.clone().into_ast(&cfg), &cfg);
        assert_eq!(lowered.ast(), dsl.ast());
    }
}
