//! Bond matchers

#![allow(clippy::type_complexity)]

use std::sync::LazyLock;

use umol::Result;

use crate::{BondBuilder, BondDonation, BondSpec, BondSpecRegistry};

pub struct BondMatcher {
    matcher: Box<dyn Fn(&BondBuilder) -> Result<Vec<BondSpec>> + Send + Sync>,
}

impl BondMatcher {
    pub fn new(matcher: Box<dyn Fn(&BondBuilder) -> Result<Vec<BondSpec>> + Send + Sync>) -> Self {
        Self { matcher }
    }

    pub fn strict() -> Self {
        Self::new(Box::new(|builder| {
            let order = builder.order();
            let candidates = BondSpecRegistry::by_order(order);
            let matches = candidates
                .iter()
                .filter(|spec| builder.donation().is_none_or(|d| d == spec.donation()))
                .cloned()
                .collect::<Vec<BondSpec>>();

            Ok(matches)
        }))
    }

    pub fn lenient() -> Self {
        Self::always()
    }

    pub fn always() -> Self {
        Self::new(Box::new(|builder| {
            Ok(vec![BondSpec::new(
                builder.order(),
                builder.donation().unwrap_or(BondDonation::Shared),
            )])
        }))
    }

    pub fn with_matcher(
        mut self,
        matcher: impl Fn(&BondBuilder) -> Result<Vec<BondSpec>> + Send + Sync + 'static,
    ) -> Self {
        self.matcher = Box::new(matcher);
        self
    }

    pub fn find(&self, builder: &BondBuilder) -> Result<Vec<BondSpec>> {
        (self.matcher)(builder)
    }
}

impl Default for BondMatcher {
    fn default() -> Self {
        Self::always()
    }
}

pub static DEFAULT_BOND_MATCHER: LazyLock<BondMatcher> = LazyLock::new(BondMatcher::default);
pub static STRICT_BOND_MATCHER: LazyLock<BondMatcher> = LazyLock::new(BondMatcher::strict);
pub static LENIENT_BOND_MATCHER: LazyLock<BondMatcher> = LazyLock::new(BondMatcher::lenient);
pub static ALWAYS_BOND_MATCHER: LazyLock<BondMatcher> = LazyLock::new(BondMatcher::always);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{b, BondBuilder, BondOrder};

    #[test]
    fn test_bond_matcher() {
        let matcher = BondMatcher::default();
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn test_bond_matcher_custom() {
        let matcher =
            BondMatcher::default().with_matcher(|_| Ok(vec![b!("-"), b!("->"), b!("-<")]));
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 3);
    }

    #[test]
    fn test_default_bond_matcher() {
        let matcher = BondMatcher::default();
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn test_strict_bond_matcher() {
        let matcher = BondMatcher::strict();
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 3);
    }

    #[test]
    fn test_lenient_bond_matcher() {
        let matcher = BondMatcher::lenient();
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn test_always_bond_matcher() {
        let matcher = BondMatcher::always();
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = matcher.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn test_default_bond_matcher_lazy_static() {
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = DEFAULT_BOND_MATCHER.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn test_always_bond_matcher_lazy_static() {
        let bond_builder = BondBuilder::new(BondOrder::Single);
        let matches = ALWAYS_BOND_MATCHER.find(&bond_builder).unwrap();
        assert_eq!(matches.len(), 1);
    }
}
