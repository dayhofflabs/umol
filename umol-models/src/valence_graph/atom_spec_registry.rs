//! Registry of default atom specs to be used for atom spec matching

use crate::{a, AtomSpec};
use once_cell::sync::Lazy;
use std::collections::HashMap;
use umol_data::{e, Element};

/// Registry of default atom specs to be used for atom spec matching
pub struct AtomSpecRegistry;

impl AtomSpecRegistry {
    /// Retrieves predefined AtomSpec instances based on element and charge.
    /// Returns an empty vector if no specs are found for the given combination.
    pub fn by_element_and_charge(element: Element, charge: i8) -> Vec<AtomSpec> {
        ATOM_SPEC_DATA
            .get(&element) // Get the inner map for the element
            .and_then(|inner_map| inner_map.get(&charge)) // Get the vec for the charge
            .cloned() // Clone the Vec<AtomSpec> if found
            .unwrap_or_else(Vec::new) // Return empty vec if element or charge not found
    }

    /// Retrieves all predefined AtomSpec instances for a given element, across all charges.
    /// Returns an empty vector if no specs are found for the given element.
    pub fn by_element(element: Element) -> Vec<AtomSpec> {
        ATOM_SPEC_DATA
            .get(&element)
            .map_or_else(
                Vec::new, // Return empty vec if element not found
                |inner_map| inner_map.values().flatten().cloned().collect() // Collect specs from all charges
            )
    }
}

// Atom specs for atom typing, nested by Element then charge
static ATOM_SPEC_DATA: Lazy<HashMap<Element, HashMap<i8, Vec<AtomSpec>>>> = Lazy::new(|| {
    let mut data = HashMap::new();

    // Helper macro for inserting into nested map
    macro_rules! insert_specs {
        ($map:expr, $element:expr, $charge:expr, [$($spec:expr),* $(,)?]) => {
            $map.entry($element)
                .or_insert_with(HashMap::new)
                .insert($charge, vec![$($spec),*]);
        };
    }

    insert_specs!(data, e!(H), 0, [a!("[H+0v1]"), a!("[H+0^1v0]")]);
    insert_specs!(data, e!(H), 1, [a!("[H+1v0]")]);
    insert_specs!(data, e!(H), -1, [a!("[H-1/1v0]")]);
    insert_specs!(data, e!(He), 0, [a!("[He+0v0]")]);
    insert_specs!(data, e!(Li), 0, [a!("[Li+0v1]"), a!("[Li+0^1v0]")]);
    insert_specs!(data, e!(Li), 1, [a!("[Li+1v0]")]);
    insert_specs!(data, e!(Be), 0, [a!("[Be+0v2]"), a!("[Be+0/1v0]")]);
    insert_specs!(data, e!(Be), 2, [a!("[Be+2v0]")]);
    insert_specs!(data, e!(B), 0, [a!("[B+0v3]"), a!("[B+0^1v2]"), a!("[B+0/1v1]"), a!("[B+0/1^1v0]")]);
    insert_specs!(data, e!(B), -1, [a!("[B-1v4]")]);
    insert_specs!(data, e!(C), 0, [a!("[C+0v4]"), a!("[C+0^1v3]"), a!("[C+0/1^2v2]"), a!("[C+0/1^2*1v2]"), a!("[C+0/1^2v0]"), a!("[C+0/1^2*1v0]")]);
    insert_specs!(data, e!(C), 1, [a!("[C+1^3v3]")]);
    insert_specs!(data, e!(C), -1, [a!("[C-1/1v3]")]);
    insert_specs!(data, e!(N), 0, [a!("[N+0/1v3]"), a!("[N+0/1^1v2]"), a!("[N+0/2^2*3v1]"), a!("[N+0/2^2*1v1]"), a!("[N+0/1^3v0]"), a!("[N+0/1^3*2v0]")]);
    insert_specs!(data, e!(N), 1, [a!("[N+1v4]"), a!("[N+1/1v2]")]);
    insert_specs!(data, e!(N), -1, [a!("[N-1/2v2]")]);
    insert_specs!(data, e!(N), -3, [a!("[N-3/4v0]")]);
    insert_specs!(data, e!(O), 0, [a!("[O+0/2v2]"), a!("[O+0/2^1v1]"), a!("[O+0/2^2v0]"), a!("[O+0/2^2*1v0]")]);
    insert_specs!(data, e!(O), 1, [a!("[O+1/1v3]"), a!("[O+1/2v1]")]);
    insert_specs!(data, e!(O), -1, [a!("[O-1/3v1]")]);
    insert_specs!(data, e!(O), -2, [a!("[O-2/4v0]")]);
    insert_specs!(data, e!(F), 0, [a!("[F+0/3v1]"), a!("[F+0/3^1v0]")]);
    insert_specs!(data, e!(F), -1, [a!("[F-1/4v0]")]);
    insert_specs!(data, e!(Ne), 0, [a!("[Ne+0/4v0]")]);

    data
});

#[cfg(test)]
mod tests {
    use super::*;
    use umol_data::e;

    #[test]
    fn test_atom_spec_registry_by_element_and_charge() {
        let atom_specs = AtomSpecRegistry::by_element_and_charge(e!(C), 0);
        assert_eq!(atom_specs.len(), 6);
        assert_eq!(atom_specs[0], a!("[C+0v4]"));

        // Test case where element exists but charge doesn't
        let atom_specs_no_charge = AtomSpecRegistry::by_element_and_charge(e!(C), 10);
        assert!(atom_specs_no_charge.is_empty());

        // Test case where element doesn't exist
        let atom_specs_no_element = AtomSpecRegistry::by_element_and_charge(e!(Og), 0);
        assert!(atom_specs_no_element.is_empty());
    }

    #[test]
    fn test_atom_spec_registry_by_element() {
        // Test Nitrogen, which has specs for charge 0, +1, -1, -3
        let nitrogen_specs = AtomSpecRegistry::by_element(e!(N));
        // Expected counts: N(0): 6, N(1): 2, N(-1): 1, N(-3): 1 => Total 10
        assert_eq!(nitrogen_specs.len(), 10);
        // Check if a spec from each charge state is present (optional sanity check)
        assert!(nitrogen_specs.contains(&a!("[N+0/1v3]"))); // N(0)
        assert!(nitrogen_specs.contains(&a!("[N+1v4]")));   // N(1)
        assert!(nitrogen_specs.contains(&a!("[N-1/2v2]"))); // N(-1)
        assert!(nitrogen_specs.contains(&a!("[N-3/4v0]"))); // N(-3)

        // Test element not in registry
        let argon_specs = AtomSpecRegistry::by_element(e!(Og));
        assert!(argon_specs.is_empty());
    }
}
