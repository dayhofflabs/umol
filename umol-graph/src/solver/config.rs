//! Top-level configuration for the resolution pipeline.

use bitflags::bitflags;

use crate::ast::rings::RingEnumerationStrategy;
use crate::solver::aromaticity::{AromaticityHintPolicy, AromaticityStrategy};
use crate::solver::propagate::ValenceStrategy;
use crate::solver::valence::AtomTypeRegistry;

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
