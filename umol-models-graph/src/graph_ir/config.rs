//! Configuration for TableIR → GraphIR resolution.

use bitflags::bitflags;
use umol_data::Element;

use super::config_data::{AtomTypeRegistry, ValenceTable};

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
    pub strategy: ValenceStrategy,
    pub no_match_policy: ValenceMatchPolicy,
    pub ambiguous_policy: ValenceMatchPolicy,
    pub enable_implicit_hydrogens: bool,
    pub atom_type_registry: AtomTypeRegistry,
    pub valence_table: ValenceTable,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceStrategy {
    AtomTyping,
    Counts,
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

#[derive(Clone, Debug)]
pub struct AromaticityResolveConfig {
    pub enabled: bool,
    pub strategy: AromaticityStrategy,
}

impl AromaticityResolveConfig {
    /// Daylight (SMILES) aromaticity: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self {
            enabled: true,
            strategy: AromaticityStrategy::HueckelRule {
                element_scope: ElementScope::AllowList(vec![
                    Element::C,
                    Element::N,
                    Element::O,
                    Element::S,
                    Element::Se,
                    Element::As,
                ]),
                ring_limits: RingLimits::default(),
            },
        }
    }

    /// MDL (MOL/SDF) aromaticity: C and N only. Minimum ring size 6.
    pub fn mdl() -> Self {
        Self {
            enabled: true,
            strategy: AromaticityStrategy::HueckelRule {
                element_scope: ElementScope::AllowList(vec![Element::C, Element::N]),
                ring_limits: RingLimits {
                    min_ring_size: 6,
                    ..RingLimits::default()
                },
            },
        }
    }

    /// Permissive aromaticity: any element with aromatic valence states.
    pub fn permissive() -> Self {
        Self {
            enabled: true,
            strategy: AromaticityStrategy::HueckelRule {
                element_scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
        }
    }
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
                strategy: ValenceStrategy::AtomTyping,
                no_match_policy: ValenceMatchPolicy::Error,
                ambiguous_policy: ValenceMatchPolicy::Error,
                enable_implicit_hydrogens: true,
                atom_type_registry: AtomTypeRegistry::default_registry().clone(),
                valence_table: ValenceTable::default_table().clone(),
            },
            aromaticity: AromaticityResolveConfig::daylight(),
            stereo: StereoResolveConfig { enabled: true },
        }
    }
}
