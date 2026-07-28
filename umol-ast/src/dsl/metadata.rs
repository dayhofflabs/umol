//! Shared query surface for rendering entity references.

use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// The rendering counterpart to [`crate::dsl::Namespace`].
///
/// Entity-keyword lookup renders an AST id as a keyword reference when one is
/// available, or as its positional index otherwise. Rendering never emits
/// structural references, so this surface needs neither counts nor participant
/// indexes.
pub trait Metadata {
    fn atom_keyword(&self, id: AtomId) -> Option<&str>;
    fn bond_keyword(&self, id: BondId) -> Option<&str>;
    fn dative_bond_keyword(&self, id: DativeBondId) -> Option<&str>;
    fn aromatic_system_keyword(&self, id: AromaticSystemId) -> Option<&str>;
    fn multicenter_bond_keyword(&self, id: MulticenterBondId) -> Option<&str>;
    fn noncovalent_bond_keyword(&self, id: NoncovalentBondId) -> Option<&str>;
    fn stereo_atom_keyword(&self, id: StereoAtomId) -> Option<&str>;
    fn stereo_bond_keyword(&self, id: StereoBondId) -> Option<&str>;
}
