//! Registry of default bond specs for GraphIR.

use std::collections::HashMap;

use map_macro::hash_map;
use std::sync::LazyLock;

use super::bond_spec::{BondOrder, BondSpec};

pub struct BondSpecRegistry;

impl BondSpecRegistry {
    pub fn by_order(order: BondOrder) -> Vec<BondSpec> {
        BOND_SPEC_DATA.get(&order).unwrap_or(&vec![]).clone()
    }
}

static BOND_SPEC_DATA: LazyLock<HashMap<BondOrder, Vec<BondSpec>>> = LazyLock::new(|| {
    fn spec(s: &str) -> BondSpec {
        s.parse::<BondSpec>().unwrap()
    }

    hash_map! {
        BondOrder::Zero => vec![spec(".")],
        BondOrder::Single => vec![spec("-"), spec("->"), spec("-<")],
        BondOrder::Double => vec![spec("="), spec("=>"), spec("=<")],
        BondOrder::Triple => vec![spec("#")],
        BondOrder::Quadruple => vec![spec("$")],
    }
});
