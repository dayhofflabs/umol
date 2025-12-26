//! Atom matcher infrastructure copied from `umol-models-valence`.

use once_cell::sync::Lazy;

use super::atom::AtomBuilder;
use super::atom_spec::AtomSpec;
use super::atom_spec_registry::AtomSpecRegistry;
use super::error::GraphError;

type Result<T> = std::result::Result<T, GraphError>;

/// Matchers for atom typing. Default matcher uses the `AtomSpecRegistry` but custom
/// matchers can be used.
pub struct AtomMatcher {
    matcher: Box<dyn Fn(&AtomBuilder) -> Result<Vec<AtomSpec>> + Send + Sync>,
}

impl AtomMatcher {
    pub fn new(matcher: Box<dyn Fn(&AtomBuilder) -> Result<Vec<AtomSpec>> + Send + Sync>) -> Self {
        Self { matcher }
    }

    /// Strict matcher that uses the `AtomSpecRegistry` to match atom builders.
    pub fn strict() -> Self {
        Self::new(Box::new(|builder| {
            let element = builder.element();

            let candidates = match builder.charge() {
                Some(charge) => AtomSpecRegistry::by_element_and_charge(element, charge),
                None => AtomSpecRegistry::by_element(element),
            };

            let matches = candidates
                .iter()
                .filter(|spec| {
                    Self::matches_charge(builder.charge(), spec.charge())
                        && Self::matches_count(builder.lone_pairs(), spec.lone_pairs())
                        && Self::matches_count(builder.donated_pairs(), spec.donated_pairs())
                        && Self::matches_count(builder.accepted_pairs(), spec.accepted_pairs())
                        && Self::matches_count(builder.unpaired_e(), spec.unpaired_e())
                        && Self::matches_count(builder.multiplicity(), spec.multiplicity())
                        && Self::matches_count(builder.implicit_h(), spec.implicit_h())
                        && Self::matches_count(builder.valence(), spec.valence())
                })
                .cloned()
                .collect::<Vec<AtomSpec>>();

            Ok(matches)
        }))
    }

    /// Lenient matcher that matches all atom builders.
    pub fn lenient() -> Self {
        Self::always()
    }

    /// Trivial matcher that matches all atom builders.
    pub fn always() -> Self {
        Self::new(Box::new(|builder| {
            let charge = builder.charge().unwrap_or(0);
            let lp = builder.lone_pairs().unwrap_or(0);
            let dp = builder.donated_pairs().unwrap_or(0);
            let ap = builder.accepted_pairs().unwrap_or(0);
            let up = builder.unpaired_e().unwrap_or(0);
            let mult = builder.multiplicity().unwrap_or(up + 1);
            let ih = builder.implicit_h().unwrap_or(0);
            let val = builder.valence().unwrap_or(0);
            Ok(vec![AtomSpec::new(
                builder.element(),
                charge,
                lp,
                dp,
                ap,
                up,
                mult,
                ih,
                val,
            )])
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

    fn matches_count(actual: Option<u32>, expected: u32) -> bool {
        actual.map_or(true, |v| v == expected)
    }

    fn matches_charge(actual: Option<i32>, expected: i32) -> bool {
        actual.map_or(true, |v| v == expected)
    }
}

impl Default for AtomMatcher {
    fn default() -> Self {
        Self::strict()
    }
}

pub static DEFAULT_ATOM_MATCHER: Lazy<AtomMatcher> = Lazy::new(AtomMatcher::default);
pub static STRICT_ATOM_MATCHER: Lazy<AtomMatcher> = Lazy::new(AtomMatcher::strict);
pub static LENIENT_ATOM_MATCHER: Lazy<AtomMatcher> = Lazy::new(AtomMatcher::lenient);
pub static ALWAYS_ATOM_MATCHER: Lazy<AtomMatcher> = Lazy::new(AtomMatcher::always);
