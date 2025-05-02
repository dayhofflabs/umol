//! Atom matchers

use crate::{AtomBuilder, AtomSpec, AtomSpecRegistry};
use once_cell::sync::Lazy;
use umol::Result;

/// Matchers for atom typing. Default matcher uses the `AtomSpecRegistry` but custom
/// matchers can be used.
pub struct AtomMatcher {
    matcher: Box<dyn Fn(&AtomBuilder) -> Result<Vec<AtomSpec>> + Send + Sync>,
}

impl AtomMatcher {
    pub fn new(matcher: Box<dyn Fn(&AtomBuilder) -> Result<Vec<AtomSpec>> + Send + Sync>) -> Self {
        Self { matcher }
    }

    pub fn default() -> Self {
        Self::new(Box::new(|builder| {
            let element = builder.element();

            // Get candidate AtomSpecs based on element and potentially charge.
            // Note: We collect into a Vec here because the subsequent filter needs ownership
            // or a stable reference, and the iterators returned by the registry methods
            // might depend on the Lazy static's internal state.
            let candidates = match builder.charge() {
                Some(charge) => AtomSpecRegistry::by_element_and_charge(element, charge),
                None => AtomSpecRegistry::by_element(element),
            };

            // Filter the candidates based on the remaining Option fields in the builder
            let matches = candidates
                .iter() // Iterate over the collected Vec<&AtomSpec>
                .filter(|spec| {
                    builder.charge().map_or(true, |c| c == spec.charge()) &&
                    builder.lone_pairs().map_or(true, |lp| lp == spec.lone_pairs()) &&
                    builder.donated_pairs().map_or(true, |dp| dp == spec.donated_pairs()) &&
                    builder.accepted_pairs().map_or(true, |ap| ap == spec.accepted_pairs()) &&
                    builder.unpaired_electrons().map_or(true, |u| u == spec.unpaired_electrons()) &&
                    builder.multiplicity().map_or(true, |m| m == spec.multiplicity()) &&
                    builder.implicit_hydrogens().map_or(true, |h| h == spec.implicit_hydrogens()) &&
                    builder.valence().map_or(true, |v| v == spec.valence())
                })
                .cloned() // Clone the matching &AtomSpec to AtomSpec
                .collect::<Vec<AtomSpec>>();

            Ok(matches)
        }))
    }

    pub fn with_matcher(
        mut self,
        matcher: impl Fn(&AtomBuilder) -> Result<Vec<AtomSpec>> + Send + Sync + 'static,
    ) -> Self {
        self.matcher = Box::new(matcher);
        self
    }

    pub fn find(&self, builder: &AtomBuilder) -> Result<Vec<AtomSpec>> {
        (self.matcher)(builder)
    }
}

impl Default for AtomMatcher {
    fn default() -> Self {
        Self::default()
    }
}

pub(super) static DEFAULT_ATOM_MATCHER: Lazy<AtomMatcher> = Lazy::new(AtomMatcher::default);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{a, AtomBuilder};
    use umol_data::{e, Element};

    #[test]
    fn test_atom_matcher() {
        let matcher = AtomMatcher::default();
        let mut atom_builder = AtomBuilder::new(Element::C);
        atom_builder.set_charge(0);
        let atom_types = matcher.find(&atom_builder);
        assert_eq!(atom_types.unwrap().len(), 6);
    }

    #[test]
    fn test_atom_matcher_custom() {
        let matcher = AtomMatcher::default().with_matcher(|_| Ok(vec![a!("[C]")]));
        let mut atom_builder = AtomBuilder::new(Element::C);
        atom_builder.set_charge(0);
        let atom_types = matcher.find(&atom_builder).unwrap();
        assert_eq!(atom_types.len(), 1);
        assert_eq!(atom_types[0].to_string(), "[C]");
    }

    #[test]
    fn test_default_atom_matcher() {
        let atom_builder = AtomBuilder::new(Element::C);
        let atom_types = DEFAULT_ATOM_MATCHER.find(&atom_builder).unwrap();
        assert_eq!(atom_types.len(), 8);
    }

    #[test]
    fn test_atom_matcher_default_element_only() {
        let matcher = AtomMatcher::default();
        let builder = AtomBuilder::new(e!(H));
        let matches = matcher.find(&builder).unwrap();
        // Should match specs for H atom, H(0), H(1), H(-1)
        assert_eq!(matches.len(), 4);
        assert!(matches.contains(&a!("[H+0v1]")));
        assert!(matches.contains(&a!("[H+0^1v0]")));
        assert!(matches.contains(&a!("[H+1v0]")));
        assert!(matches.contains(&a!("[H-1/1v0]")));
    }

    #[test]
    fn test_atom_matcher_default_element_charge() {
        let matcher = AtomMatcher::default();
        let mut builder = AtomBuilder::new(e!(C));
        builder.set_charge(0);
        let matches = matcher.find(&builder).unwrap();
        // Should match only specs for C(0)
        assert_eq!(matches.len(), 6);
        assert!(matches.iter().all(|spec| spec.element() == e!(C) && spec.charge() == 0));
    }

    #[test]
    fn test_atom_matcher_default_partial_match() {
        let matcher = AtomMatcher::default();
        let mut builder = AtomBuilder::new(e!(N));
        builder.set_charge(0);
        builder.set_lone_pairs(1);
        builder.set_valence(3);
        let matches = matcher.find(&builder).unwrap();
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], a!("[N+0/1v3]"));
    }

    #[test]
    fn test_atom_matcher_default_partial_no_match_invalid_prop() {
        let matcher = AtomMatcher::default();
        let mut builder = AtomBuilder::new(e!(O));
        builder.set_charge(5); // No O(5) specs defined
        let matches = matcher.find(&builder).unwrap();
        assert!(matches.is_empty());
    }

    #[test]
    fn test_atom_matcher_default_partial_no_match_inconsistent() {
        let matcher = AtomMatcher::default();
        let mut builder = AtomBuilder::new(e!(C));
        builder.set_charge(0);
        builder.set_valence(5); // No C(0) spec with valence 5 defined
        let matches = matcher.find(&builder).unwrap();
        assert!(matches.is_empty());
    }

    #[test]
    fn test_atom_matcher_default_partial_multiple_matches() {
        let matcher = AtomMatcher::default();
        let mut builder = AtomBuilder::new(e!(C));
        builder.set_charge(0);
        builder.set_valence(2);
        let matches = matcher.find(&builder).unwrap();
        // Should match triplet [C+0/1^2v2] and singlet [C+0/1^2*1v2]
        assert_eq!(matches.len(), 2);
        assert!(matches.contains(&a!("[C+0/1^2v2]")));
        assert!(matches.contains(&a!("[C+0/1^2*1v2]")));
    }

    #[test]
    fn test_atom_matcher_default_element_not_in_registry() {
        let matcher = AtomMatcher::default();
        let builder = AtomBuilder::new(e!(Og));
        let matches = matcher.find(&builder).unwrap();
        assert!(matches.is_empty());
    }

    #[test]
    fn test_default_atom_matcher_lazy_static() {
        let atom_builder = AtomBuilder::new(e!(C));
        let atom_types = DEFAULT_ATOM_MATCHER.find(&atom_builder).unwrap();
        assert_eq!(atom_types.len(), 8);
    }
}
