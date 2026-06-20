//! Valence-invariants validator: a thin adapter over [`ValenceInvariants`],
//! surfacing its per-atom electron-conservation verdict (orbital occupancy vs
//! `Z − q` plus the electrons the atom shares into bonds and implicit H) as a
//! validator `Solution`. The counts are per-atom and are never summed across
//! the molecule — shared bonding electrons are counted on each endpoint.

use thiserror::Error;
use umol_ast::ast::{AtomAst, MoleculeAst};

use crate::ops::invariant::{ValenceInvariants, ValenceMismatch};
use umol_shared::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct ValenceInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceInvariantsError {}

impl ValenceInvariantsValidator {
    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ValenceMismatch>, ValenceInvariantsError> {
        Ok(ValenceInvariants::check(ast.as_ref()))
    }

    pub fn validate_atom(
        &self,
        atom: &AtomAst,
    ) -> Result<Solution<(), ValenceMismatch>, ValenceInvariantsError> {
        Ok(ValenceInvariants::check_atom(atom))
    }
}
