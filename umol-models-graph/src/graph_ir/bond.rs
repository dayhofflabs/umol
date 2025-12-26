//! Bond type and builder copied from `umol-models-valence`.

use std::fmt;

use serde::{Deserialize, Serialize};

use super::bond_matcher::{BondMatcher, DEFAULT_BOND_MATCHER};
use super::bond_spec::{BondDonation, BondOrder, BondSpec};
use super::error::GraphError;
use crate::simple_ir::{BondDir, BondStereo, BondSymbol};
use crate::span::Span;

type Result<T> = std::result::Result<T, GraphError>;

/// Valence bond type including strict typing.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    order: BondOrder,
    donation: BondDonation,
    symbol: Option<BondSymbol>,
    direction: Option<BondDir>,
    stereo: Option<BondStereo>,
    ring: Option<u32>,
    span: Option<Span>,
}

impl Bond {
    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> BondDonation {
        self.donation
    }

    pub fn symbol(&self) -> Option<BondSymbol> {
        self.symbol
    }

    pub fn direction(&self) -> Option<BondDir> {
        self.direction
    }

    pub fn stereo(&self) -> Option<BondStereo> {
        self.stereo
    }

    pub fn ring(&self) -> Option<u32> {
        self.ring
    }

    pub fn span(&self) -> Option<Span> {
        self.span
    }

    pub fn span_bytes(&self) -> Option<(u32, u32)> {
        self.span.and_then(|s| s.bytes_range())
    }

    pub fn to_builder(self) -> BondBuilder {
        BondBuilder {
            order: self.order,
            donation: Some(self.donation),
            symbol: self.symbol,
            direction: self.direction,
            stereo: self.stereo,
            ring: self.ring,
            span: self.span,
        }
    }

    pub fn to_spec(&self) -> BondSpec {
        BondSpec::new(self.order, self.donation)
    }
}

impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_spec())
    }
}

impl From<Bond> for BondBuilder {
    fn from(bond: Bond) -> Self {
        bond.to_builder()
    }
}

/// Builder type for creating and mutating `Bond` types including strict typing.
#[derive(Debug)]
pub struct BondBuilder {
    order: BondOrder,
    donation: Option<BondDonation>,
    symbol: Option<BondSymbol>,
    direction: Option<BondDir>,
    stereo: Option<BondStereo>,
    ring: Option<u32>,
    span: Option<Span>,
}

impl BondBuilder {
    pub fn new(order: BondOrder) -> Self {
        Self {
            order,
            donation: None,
            symbol: None,
            direction: None,
            stereo: None,
            ring: None,
            span: None,
        }
    }

    pub fn from_spec(bond_spec: BondSpec) -> Self {
        Self {
            order: bond_spec.order(),
            donation: Some(bond_spec.donation()),
            symbol: None,
            direction: None,
            stereo: None,
            ring: None,
            span: None,
        }
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }

    pub fn symbol(&self) -> Option<BondSymbol> {
        self.symbol
    }

    pub fn direction(&self) -> Option<BondDir> {
        self.direction
    }

    pub fn stereo(&self) -> Option<BondStereo> {
        self.stereo
    }

    pub fn ring(&self) -> Option<u32> {
        self.ring
    }

    pub fn span(&self) -> Option<Span> {
        self.span
    }

    pub fn span_bytes(&self) -> Option<(u32, u32)> {
        self.span.and_then(|s| s.bytes_range())
    }

    pub fn set_order(&mut self, order: BondOrder) -> &mut Self {
        self.order = order;
        self
    }

    pub fn set_donation(&mut self, donation: BondDonation) -> &mut Self {
        self.donation = Some(donation);
        self
    }

    pub fn set_symbol(&mut self, symbol: BondSymbol) -> &mut Self {
        self.symbol = Some(symbol);
        self
    }

    pub fn set_direction(&mut self, direction: Option<BondDir>) -> &mut Self {
        self.direction = direction;
        self
    }

    pub fn set_stereo(&mut self, stereo: Option<BondStereo>) -> &mut Self {
        self.stereo = stereo;
        self
    }

    pub fn set_ring(&mut self, ring: Option<u32>) -> &mut Self {
        self.ring = ring;
        self
    }

    pub fn set_span(&mut self, start: Option<u32>, end: Option<u32>) -> &mut Self {
        self.span = Span::from_bytes_opt(start, end);
        self
    }

    pub fn set_span_opt(&mut self, span: Option<Span>) -> &mut Self {
        self.span = span;
        self
    }

    pub fn update_order(&mut self, f: impl FnOnce(BondOrder) -> BondOrder) -> &mut Self {
        self.order = f(self.order);
        self
    }

    pub fn update_donation(&mut self, f: impl FnOnce(BondDonation) -> BondDonation) -> &mut Self {
        self.donation = Some(f(self.donation.unwrap_or(BondDonation::Shared)));
        self
    }

    pub fn build(self) -> Result<Bond> {
        self.build_with(&DEFAULT_BOND_MATCHER)
    }

    pub fn build_with(self, matcher: &BondMatcher) -> Result<Bond> {
        let bond_specs = matcher.find(&self)?;
        if bond_specs.is_empty() {
            return Err(GraphError::InvalidBondSpec(format!("{:?}", self)));
        } else if bond_specs.len() > 1 {
            return Err(GraphError::InvalidBondSpec(format!("{:?}", self)));
        }
        let bond_spec = bond_specs.first().unwrap();
        Ok(Bond {
            order: bond_spec.order(),
            donation: bond_spec.donation(),
            symbol: self.symbol,
            direction: self.direction,
            stereo: self.stereo,
            ring: self.ring,
            span: self.span,
        })
    }
}

impl From<BondSpec> for BondBuilder {
    fn from(bond_spec: BondSpec) -> Self {
        BondBuilder::from_spec(bond_spec)
    }
}

impl From<BondOrder> for BondBuilder {
    fn from(order: BondOrder) -> Self {
        BondBuilder::new(order)
    }
}
