//! Configuration for TableIR → GraphIR resolution.

use bitflags::bitflags;

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
            valence: ValenceResolveConfig { enabled: true },
            aromaticity: AromaticityResolveConfig { enabled: true },
            stereo: StereoResolveConfig { enabled: true },
        }
    }
}
