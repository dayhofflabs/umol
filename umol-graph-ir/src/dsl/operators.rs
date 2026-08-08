//! DSL parse / render for the membership and relational operators — the string-DSL
//! surface of `ast::operators` (`MemOp` / `RelOp`), shared by the value predicates and
//! the element bind.

use winnow::combinator::alt;
use winnow::Parser;

use super::error::PResult;
use crate::ast::operators::{MemOp, RelOp};

pub(crate) fn rel_op(i: &mut &str) -> PResult<RelOp> {
    alt((
        "<=".value(RelOp::Le),
        ">=".value(RelOp::Ge),
        "==".value(RelOp::Eq),
        "!=".value(RelOp::Ne),
        '<'.value(RelOp::Lt),
        '>'.value(RelOp::Gt),
    ))
    .parse_next(i)
}

pub(crate) fn rel_op_str(op: RelOp) -> &'static str {
    match op {
        RelOp::Le => "<=",
        RelOp::Ge => ">=",
        RelOp::Eq => "==",
        RelOp::Lt => "<",
        RelOp::Gt => ">",
        RelOp::Ne => "!=",
    }
}

pub(crate) fn mem_op(i: &mut &str) -> PResult<MemOp> {
    alt(("!:".value(MemOp::NotIn), "::".value(MemOp::In))).parse_next(i)
}

pub(crate) fn mem_op_str(op: MemOp) -> &'static str {
    match op {
        MemOp::In => "::",
        MemOp::NotIn => "!:",
    }
}
