//! Configuration for TableIR → GraphIR resolution.

use bitflags::bitflags;
use umol_shared::element::Element;

pub use crate::solver::ValenceStrategy;

use super::config_data::AtomTypeRegistry;

bitflags! {
    /// Topology resolution options.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct TopologyResolveFlags: u32 {
        /// Allow molecules with more than one connected component.
        const DISCONNECTED_MOLECULES = 1;
    }
}

impl Default for TopologyResolveFlags {
    fn default() -> Self {
        Self::empty()
    }
}

/// Configuration controlling which validation phases run and with what strictness.
#[derive(Clone, Debug)]
pub struct ResolveConfig {
    pub topology: TopologyConfig,
    pub valence: ValenceResolveConfig,
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}

#[derive(Clone, Debug)]
pub struct TopologyConfig {
    pub flags: TopologyResolveFlags,
}

#[derive(Clone, Debug)]
pub struct ValenceResolveConfig {
    pub strategy: ValenceStrategy,
    pub no_match_policy: ValenceMatchPolicy,
    pub ambiguous_policy: ValenceMatchPolicy,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceMatchPolicy {
    Error,
    Ignore,
}

/// Elements eligible for aromaticity perception
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementScope {
    Any,
    AllowList(Vec<Element>),
}

/// Ring size and fused-ring constraints for HueckelRule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingLimits {
    pub min_ring_size: usize,
    pub max_ring_size: usize,
    pub include_fused: bool,
    pub max_fused_combination: usize,
    pub max_fused_search: usize,
}

impl Default for RingLimits {
    fn default() -> Self {
        Self {
            min_ring_size: 3,
            max_ring_size: 22,
            include_fused: true,
            max_fused_combination: 6,
            max_fused_search: 10_000,
        }
    }
}

/// Ring enumeration parameters, independent of aromaticity model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingEnumerationStrategy {
    /// If true, restrict ring enumeration to atoms with aromatic hints
    /// (builder) or aromatic types (molecule).
    pub aromatic_only: bool,
    pub max_ring_size: usize,
    pub max_rings_per_component: usize,
}

impl Default for RingEnumerationStrategy {
    fn default() -> Self {
        Self {
            aromatic_only: false,
            max_ring_size: 22,
            max_rings_per_component: 2000,
        }
    }
}

#[derive(Clone, Debug)]
pub enum AromaticityStrategy {
    HueckelRule {
        element_scope: ElementScope,
        ring_limits: RingLimits,
    },
    Hmo {
        element_scope: ElementScope,
        /// Delocalization energy per pi-electron (in units of |beta|) required
        /// for classification as aromatic. Benzene: dE/n ~ 0.33|beta|.
        stabilization_threshold: f64,
    },
    Clar,
}

/// Policy for mismatches between aromatic hints and detected aromatic systems.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityHintPolicy {
    Strict,
    Ignore,
}

impl AromaticityStrategy {
    /// Daylight (SMILES) aromaticity: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::AllowList(vec![
                Element::C,
                Element::N,
                Element::O,
                Element::S,
                Element::Se,
                Element::As,
            ]),
            ring_limits: RingLimits::default(),
        }
    }

    /// MDL (MOL/SDF) aromaticity: C and N only. Minimum ring size 6.
    pub fn mdl() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        }
    }

    /// Permissive aromaticity: any element with aromatic valence states.
    pub fn permissive() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct AromaticityResolveConfig {
    pub aromaticity_strategy: AromaticityStrategy,
    pub enumeration_strategy: RingEnumerationStrategy,
    pub hint_policy: AromaticityHintPolicy,
}

#[derive(Clone, Debug)]
pub struct StereoResolveConfig {}

impl Default for ResolveConfig {
    fn default() -> Self {
        Self {
            topology: TopologyConfig {
                flags: TopologyResolveFlags::empty(),
            },
            valence: ValenceResolveConfig {
                strategy: ValenceStrategy::AtomTyping {
                    registry: AtomTypeRegistry::default_registry().clone(),
                },
                no_match_policy: ValenceMatchPolicy::Error,
                ambiguous_policy: ValenceMatchPolicy::Error,
            },
            aromaticity: AromaticityResolveConfig {
                aromaticity_strategy: AromaticityStrategy::daylight(),
                enumeration_strategy: RingEnumerationStrategy::default(),
                hint_policy: AromaticityHintPolicy::Strict,
            },
            stereo: StereoResolveConfig {},
        }
    }
}
