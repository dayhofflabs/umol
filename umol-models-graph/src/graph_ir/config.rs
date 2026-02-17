//! Configuration for TableIR → GraphIR resolution.

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
            topology: TopologyConfig { enabled: true },
            valence: ValenceResolveConfig { enabled: true },
            aromaticity: AromaticityResolveConfig { enabled: true },
            stereo: StereoResolveConfig { enabled: true },
        }
    }
}
