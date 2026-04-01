//! Common bond data types for graph-based molecular models.

use serde::{Deserialize, Serialize};

/// Non-covalent interaction type for weak bonds (H-bonds, halogen bonds, etc.)
/// These bonds do not contribute to valence calculations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondNoncovalent {
    Hydrogen,
}
