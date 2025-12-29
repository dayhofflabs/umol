//! Registry of default bond specs for bond spec matching

use std::collections::HashMap;
use std::sync::LazyLock;

use map_macro::hash_map;

use crate::{b, BondOrder, BondSpec};

/// Registry of bond specs by bond order
pub struct BondSpecRegistry;

impl BondSpecRegistry {
    pub fn by_order(order: BondOrder) -> Vec<BondSpec> {
        BOND_SPEC_DATA.get(&order).unwrap_or(&vec![]).clone()
    }
}

// Bond specs for bond typing
static BOND_SPEC_DATA: LazyLock<HashMap<BondOrder, Vec<BondSpec>>> = LazyLock::new(|| {
    hash_map! {
        BondOrder::Zero => vec![b!(".")],
        BondOrder::Single => vec![b!("-"), b!("->"), b!("-<")],
        BondOrder::Double => vec![b!("="), b!("=>"), b!("=<")],
        BondOrder::Triple => vec![b!("#")],
        BondOrder::Quadruple => vec![b!("$")],
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_spec_registry() {
        assert_eq!(
            BondSpecRegistry::by_order(BondOrder::Single),
            vec![b!("-"), b!("->"), b!("-<")]
        );
    }
}
