//! Multicenter bond data structs for molecule AST.

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {}

impl MulticenterBondAst {
    pub fn is_ground(&self) -> bool {
        true
    }
}
