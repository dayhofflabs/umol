//! Registry of default atom specs for GraphIR atom typing.

use std::collections::HashMap;

use once_cell::sync::Lazy;
use umol_data::{e, Element};

use super::atom_spec::AtomSpec;

/// Registry of default atom specs to be used for atom spec matching
pub struct AtomSpecRegistry;

impl AtomSpecRegistry {
    pub fn by_element_and_charge(element: Element, charge: i32) -> Vec<AtomSpec> {
        ATOM_SPEC_DATA
            .get(&element)
            .and_then(|inner_map| inner_map.get(&charge))
            .cloned()
            .unwrap_or_else(Vec::new)
    }

    pub fn by_element(element: Element) -> Vec<AtomSpec> {
        ATOM_SPEC_DATA.get(&element).map_or_else(Vec::new, |inner| {
            inner.values().flatten().cloned().collect()
        })
    }
}

static ATOM_SPEC_DATA: Lazy<HashMap<Element, HashMap<i32, Vec<AtomSpec>>>> = Lazy::new(|| {
    fn spec(s: &str) -> AtomSpec {
        s.parse::<AtomSpec>().unwrap()
    }

    let mut data = HashMap::new();

    macro_rules! insert_specs {
        ($map:expr, $element:expr, $charge:expr, [$($spec:expr),* $(,)?] $(,)?) => {
            $map.entry($element)
                .or_insert_with(HashMap::new)
                .insert($charge, vec![$($spec),*]);
        };
    }

    insert_specs!(data, e!(H), 0, [spec("[H+0v1]"), spec("[H+0^1v0]")]);
    insert_specs!(data, e!(H), 1, [spec("[H+1v0]")]);
    insert_specs!(data, e!(H), -1, [spec("[H-1/1v0]")]);
    insert_specs!(data, e!(He), 0, [spec("[He+0v0]")]);
    insert_specs!(data, e!(Li), 0, [spec("[Li+0v1]"), spec("[Li+0^1v0]")]);
    insert_specs!(data, e!(Li), 1, [spec("[Li+1v0]")]);
    insert_specs!(data, e!(Be), 0, [spec("[Be+0v2]"), spec("[Be+0/1v0]")]);
    insert_specs!(data, e!(Be), 2, [spec("[Be+2v0]")]);
    insert_specs!(
        data,
        e!(B),
        0,
        [
            spec("[B+0v3]"),
            spec("[B+0^1v2]"),
            spec("[B+0/1v1]"),
            spec("[B+0/1^1v0]"),
        ],
    );
    insert_specs!(data, e!(B), -1, [spec("[B-1v4]")]);
    insert_specs!(
        data,
        e!(C),
        0,
        [
            spec("[C+0v4]"),
            spec("[C+0^1v3]"),
            spec("[C+0/1^2v2]"),
            spec("[C+0/1^2*1v2]"),
            spec("[C+0/1^2v0]"),
            spec("[C+0/1^2*1v0]"),
        ],
    );
    insert_specs!(data, e!(C), 1, [spec("[C+1^3v3]")]);
    insert_specs!(data, e!(C), -1, [spec("[C-1/1v3]")]);
    insert_specs!(
        data,
        e!(N),
        0,
        [
            spec("[N+0/1v3]"),
            spec("[N+0/1^1v2]"),
            spec("[N+0/2^2*3v1]"),
            spec("[N+0/2^2*1v1]"),
            spec("[N+0/1^3v0]"),
            spec("[N+0/1^3*2v0]"),
        ],
    );
    insert_specs!(data, e!(N), 1, [spec("[N+1v4]"), spec("[N+1/1v2]")]);
    insert_specs!(data, e!(N), -1, [spec("[N-1/2v2]")]);
    insert_specs!(data, e!(N), -3, [spec("[N-3/4v0]")]);
    insert_specs!(
        data,
        e!(O),
        0,
        [
            spec("[O+0/2v2]"),
            spec("[O+0/2^1v1]"),
            spec("[O+0/2^2v0]"),
            spec("[O+0/2^2*1v0]"),
        ],
    );
    insert_specs!(data, e!(O), 1, [spec("[O+1/1v3]"), spec("[O+1/2v1]")]);
    insert_specs!(data, e!(O), -1, [spec("[O-1/3v1]")]);
    insert_specs!(data, e!(O), -2, [spec("[O-2/4v0]")]);
    insert_specs!(data, e!(F), 0, [spec("[F+0/3v1]"), spec("[F+0/3^1v0]")]);
    insert_specs!(data, e!(F), -1, [spec("[F-1/4v0]")]);
    insert_specs!(data, e!(Ne), 0, [spec("[Ne+0/4v0]")]);

    data
});
