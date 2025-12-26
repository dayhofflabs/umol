//! Bond matcher infrastructure copied from `umol-models-valence`.

use once_cell::sync::Lazy;

use super::bond::BondBuilder;
use super::bond_spec::{BondDonation, BondSpec};
use super::bond_spec_registry::BondSpecRegistry;
use super::error::GraphError;

type Result<T> = std::result::Result<T, GraphError>;

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

pub static DEFAULT_BOND_MATCHER: Lazy<BondMatcher> = Lazy::new(BondMatcher::default);
pub static STRICT_BOND_MATCHER: Lazy<BondMatcher> = Lazy::new(BondMatcher::strict);
pub static LENIENT_BOND_MATCHER: Lazy<BondMatcher> = Lazy::new(BondMatcher::lenient);
pub static ALWAYS_BOND_MATCHER: Lazy<BondMatcher> = Lazy::new(BondMatcher::always);
