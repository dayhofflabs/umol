//! Configuration for TableIR → GraphIR resolution.

use bitflags::bitflags;

use super::valence::AtomTypeRegistry;

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
    pub enabled: bool,
    pub flags: TopologyResolveFlags,
}

#[derive(Clone, Debug)]
pub struct ValenceResolveConfig {
    pub enabled: bool,
    pub strategy: ValenceStrategyKind,
    pub no_match_policy: ValenceMatchPolicy,
    pub ambiguous_policy: ValenceMatchPolicy,
    pub registry: AtomTypeRegistry,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceStrategyKind {
    AtomTyping,
    Counts,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceMatchPolicy {
    Error,
    Ignore,
}

#[derive(Clone, Debug)]
pub struct AromaticityResolveConfig {
    pub enabled: bool,
}

#[derive(Clone, Debug)]
pub struct StereoResolveConfig {
    pub enabled: bool,
}

impl Default for ResolveConfig {
    fn default() -> Self {
        Self {
            topology: TopologyConfig {
                enabled: true,
                flags: TopologyResolveFlags::empty(),
            },
            valence: ValenceResolveConfig {
                enabled: true,
                strategy: ValenceStrategyKind::AtomTyping,
                no_match_policy: ValenceMatchPolicy::Error,
                ambiguous_policy: ValenceMatchPolicy::Error,
                registry: AtomTypeRegistry::default_registry().clone(),
            },
            aromaticity: AromaticityResolveConfig { enabled: true },
            stereo: StereoResolveConfig { enabled: true },
        }
    }
}
