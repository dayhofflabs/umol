//! Transformer properties grouped by representation and operation.
//!
//! Aromaticity properties cover publication, valence, and assertion preservation.
//! Hydrogen properties use independent representations, conservation laws, and
//! count/stereo rejection domains.

#[path = "transform/aromaticity.rs"]
mod aromaticity;
#[path = "transform/hydrogens.rs"]
mod hydrogens;
