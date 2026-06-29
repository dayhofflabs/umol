//! Reaction span DSL: the surface form of `ReactionSpanAst`, where each entity carries its
//! complete before/after value (`EntitySpan`) rather than a delta. Entity ids, bond endpoints,
//! and constraint topology refs are resolved in `into_ast`.

use super::constraint::ConstraintDsl;
use super::molecule::MoleculeMetadata;
use super::refs::AtomRef;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::{EntitySpan, ReactionSpanAst};

/// Surface DSL for a reaction span. Pairs `ReactionSpanAst` with the `MoleculeMetadata` recording
/// its span-frame id↔name bindings; fields private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanDsl {
    ast: ReactionSpanAst,
    metadata: MoleculeMetadata,
}

impl ReactionSpanDsl {
    pub fn from_parts(ast: ReactionSpanAst, metadata: MoleculeMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &ReactionSpanAst {
        &self.ast
    }

    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionSpanAst, MoleculeMetadata) {
        (self.ast, self.metadata)
    }
}

/// One molecule-level constraint's span, with its refs still unresolved.
#[derive(Debug)]
pub(crate) enum ConstraintSpanInput {
    Unchanged(ConstraintDsl),
    Added(ConstraintDsl),
    Removed(ConstraintDsl),
}

/// A parsed reaction span before ref resolution: each atom/bond carries its complete `EntitySpan`
/// value plus an optional surface id; bond endpoints stay unresolved until `into_ast`.
#[derive(Debug)]
pub(crate) struct SpanInput {
    atoms: Vec<(Option<String>, EntitySpan<AtomAst>)>,
    bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>)>,
    constraints: Vec<ConstraintSpanInput>,
}
