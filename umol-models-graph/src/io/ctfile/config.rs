//! Parsing configuration for CTab-based formats

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing CTab-based formats
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct CtabParseFlags: u32 {
        // Parser capabilities
        const WILDCARDS = 1;                // A, Q, *, X, M
        const CHEMAXON_WILDCARDS = 1 << 1;  // AH, QH, XH, MH
        const ELECTRONS = 1 << 2;           // LP (lone pairs)
        const RGROUPS = 1 << 3;             // APO, AAL, RGP, LOG
        const SGROUPS = 1 << 4;             // STY, SST, SLB, SAL, SBL, SMT, SCN,
                                            // SDS, SPA, CRS, SDI, SBV, SDT, SDD, SCD, SED, SPL, SNC
        const QUERY_BONDS = 1 << 5;         // Bond types 5-8
        const QUERY_PROPERTIES = 1 << 6;    // RBC, SUB, UNS, LIN, ALS, APO, AAL
        const LEGACY_ATOM_LISTS = 1 << 7;   // Legacy atom lists
        // const RESERVED_8 = 1 << 8;
        // const RESERVED_9 = 1 << 9;

        // Data validation strictness
        const EXTENDED_RANGE = 1 << 10;     // Extended range of values (bond orders, H counts)
        const EXTENDED_ISOTOPES = 1 << 11;  // Extended isotopes (no catalog check)
        // const RESERVED_12 = 1 << 12;
        // const RESERVED_13 = 1 << 13;

        // Format extensions
        const NAMED_ISOTOPES = 1 << 14;         // Allow D, T as element symbols
        const PSEUDOATOMS = 1 << 15;            // General pseudoatoms (any non-reserved string)
        const ATOM_MAP_HCOUNT_FIELDS = 1 << 16; // Allow atom mapping and H count in basic parser
        const CLARK_EXTENSIONS = 1 << 17;       // ZBO, ZCH, HYD
        const EDITOR_EXTENSIONS = 1 << 18;      // ZZC (ACD/ChemSketch)
        // const RESERVED_19 = 1 << 19;
        // const RESERVED_20 = 1 << 20;
        // const RESERVED_21 = 1 << 21;

        // Input validation strictness
        const UNICODE = 1 << 22;             // Allow Unicode whitespace
        const SKIP_UNUSED_FIELDS = 1 << 23;  // Skip validation of unused fields
        const NO_V2000_END_TAGS = 1 << 24;   // V2000 tag and M  END tag may be omitted
        // const RESERVED_25 = 1 << 25;

        // Parser behavior
        const IGNORE_POSITIONS = 1 << 26;   // Ignore position data
        // const RESERVED_27 = 1 << 27;
        // const RESERVED_28 = 1 << 28;
        // const RESERVED_29 = 1 << 29;
        // const RESERVED_30 = 1 << 30;
        // const RESERVED_31 = 1 << 31;

        // Presets
        const MINIMAL = 0;

        // Maximum capabilities for basic parser
        const BASIC_MAX = Self::EXTENDED_RANGE.bits() | Self::EXTENDED_ISOTOPES.bits() |
            Self::NAMED_ISOTOPES.bits() | Self::ATOM_MAP_HCOUNT_FIELDS.bits() | Self::CLARK_EXTENSIONS.bits() |
            Self::EDITOR_EXTENSIONS.bits() | Self::UNICODE.bits() | Self::SKIP_UNUSED_FIELDS.bits() |
            Self::NO_V2000_END_TAGS.bits() | Self::IGNORE_POSITIONS.bits();

        // Maximum capabilities for extended parser (everything)
        const EXTENDED_MAX = Self::BASIC_MAX.bits() | Self::WILDCARDS.bits() | Self::CHEMAXON_WILDCARDS.bits() |
            Self::ELECTRONS.bits() | Self::PSEUDOATOMS.bits() | Self::RGROUPS.bits() | Self::SGROUPS.bits() |
            Self::QUERY_BONDS.bits() | Self::QUERY_PROPERTIES.bits() | Self::LEGACY_ATOM_LISTS.bits();

        // Strict parser: only additional capabilities of extended parser over basic parser according to spec
        // Basic strict: BASIC & STRICT = MINIMAL
        const STRICT = Self::WILDCARDS.bits() | Self::ELECTRONS.bits() |
            Self::RGROUPS.bits() | Self::SGROUPS.bits() | Self::QUERY_BONDS.bits() | Self::QUERY_PROPERTIES.bits();

        // Default for basic parser
        const BASIC = Self::NAMED_ISOTOPES.bits() | Self::ATOM_MAP_HCOUNT_FIELDS.bits() | Self::SKIP_UNUSED_FIELDS.bits();

        // Default for extended parser
        const EXTENDED = Self::BASIC.bits() | Self::STRICT.bits();

        // Lenient parser
        const LENIENT = Self::EXTENDED.bits() | Self::CHEMAXON_WILDCARDS.bits() | Self::LEGACY_ATOM_LISTS.bits() |
            Self::EXTENDED_RANGE.bits() | Self::EXTENDED_ISOTOPES.bits() | Self::NAMED_ISOTOPES.bits() |
            Self::CLARK_EXTENSIONS.bits() | Self::EDITOR_EXTENSIONS.bits() | Self::UNICODE.bits() |
            Self::SKIP_UNUSED_FIELDS.bits() | Self::NO_V2000_END_TAGS.bits();

        // Graph-only parser
        const GRAPH_ONLY = Self::BASIC.bits() | Self::EXTENDED_RANGE.bits() | Self::EXTENDED_ISOTOPES.bits() |
            Self::IGNORE_POSITIONS.bits();

        // Smiles-compatible parser
        const SMILES_COMPAT = Self::BASIC.bits() | Self::EXTENDED_ISOTOPES.bits() | Self::IGNORE_POSITIONS.bits();

        // Smarts-compatible parser
        const SMARTS_COMPAT = Self::SMILES_COMPAT.bits() | Self::WILDCARDS.bits() | Self::QUERY_PROPERTIES.bits() |
            Self::QUERY_BONDS.bits();
    }
}

impl fmt::Display for CtabParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();

        if *self == CtabParseFlags::MINIMAL {
            parts.push("MINIMAL");
        } else if *self == CtabParseFlags::BASIC {
            parts.push("BASIC");
        } else if *self == CtabParseFlags::BASIC_MAX {
            parts.push("BASIC_MAX");
        } else if *self == CtabParseFlags::STRICT {
            parts.push("STRICT");
        } else if *self == CtabParseFlags::EXTENDED {
            parts.push("EXTENDED");
        } else if *self == CtabParseFlags::EXTENDED_MAX {
            parts.push("EXTENDED_MAX");
        } else if *self == CtabParseFlags::LENIENT {
            parts.push("LENIENT");
        } else if *self == CtabParseFlags::GRAPH_ONLY {
            parts.push("GRAPH_ONLY");
        } else if *self == CtabParseFlags::SMILES_COMPAT {
            parts.push("SMILES_COMPAT");
        } else if *self == CtabParseFlags::SMARTS_COMPAT {
            parts.push("SMARTS_COMPAT");
        } else {
            // Show individual flags
            if self.contains(CtabParseFlags::WILDCARDS) {
                parts.push("WILDCARDS");
            }
            if self.contains(CtabParseFlags::CHEMAXON_WILDCARDS) {
                parts.push("CHEMAXON_WILDCARDS");
            }
            if self.contains(CtabParseFlags::ELECTRONS) {
                parts.push("ELECTRONS");
            }
            if self.contains(CtabParseFlags::RGROUPS) {
                parts.push("RGROUPS");
            }
            if self.contains(CtabParseFlags::SGROUPS) {
                parts.push("SGROUPS");
            }
            if self.contains(CtabParseFlags::QUERY_BONDS) {
                parts.push("QUERY_BONDS");
            }
            if self.contains(CtabParseFlags::QUERY_PROPERTIES) {
                parts.push("QUERY_PROPERTIES");
            }
            if self.contains(CtabParseFlags::LEGACY_ATOM_LISTS) {
                parts.push("LEGACY_ATOM_LISTS");
            }
            if self.contains(CtabParseFlags::EXTENDED_RANGE) {
                parts.push("EXTENDED_RANGE");
            }
            if self.contains(CtabParseFlags::EXTENDED_ISOTOPES) {
                parts.push("EXTENDED_ISOTOPES");
            }
            if self.contains(CtabParseFlags::NAMED_ISOTOPES) {
                parts.push("NAMED_ISOTOPES");
            }
            if self.contains(CtabParseFlags::PSEUDOATOMS) {
                parts.push("PSEUDOATOMS");
            }
            if self.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS) {
                parts.push("ATOM_MAP_HCOUNT_FIELDS");
            }
            if self.contains(CtabParseFlags::CLARK_EXTENSIONS) {
                parts.push("CLARK_EXTENSIONS");
            }
            if self.contains(CtabParseFlags::EDITOR_EXTENSIONS) {
                parts.push("EDITOR_EXTENSIONS");
            }
            if self.contains(CtabParseFlags::UNICODE) {
                parts.push("UNICODE");
            }
            if self.contains(CtabParseFlags::SKIP_UNUSED_FIELDS) {
                parts.push("SKIP_UNUSED_FIELDS");
            }
            if self.contains(CtabParseFlags::NO_V2000_END_TAGS) {
                parts.push("NO_V2000_END_TAGS");
            }
            if self.contains(CtabParseFlags::IGNORE_POSITIONS) {
                parts.push("IGNORE_POSITIONS");
            }
        }
        write!(f, "{}", parts.join(" | "))
    }
}

impl Default for CtabParseFlags {
    fn default() -> Self {
        Self::MINIMAL
    }
}

/// Configuration for CTFile (MOL/SDF) parsing/writing
#[derive(Debug, Clone, Default)]
pub struct CtfileIoConfig {
    pub parse_flags: CtabParseFlags,
}

/// Backwards compatibility alias
pub type MolIoConfig = CtfileIoConfig;

impl CtfileIoConfig {
    pub fn with_parse_flags(flags: CtabParseFlags) -> Self {
        Self { parse_flags: flags }
    }

    pub fn minimal() -> Self {
        Self::with_parse_flags(CtabParseFlags::MINIMAL)
    }

    pub fn basic_strict() -> Self {
        Self::with_parse_flags(CtabParseFlags::BASIC & CtabParseFlags::STRICT)
    }

    pub fn basic() -> Self {
        Self::with_parse_flags(CtabParseFlags::BASIC)
    }

    pub fn basic_lenient() -> Self {
        Self::with_parse_flags(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT)
    }

    pub fn basic_max() -> Self {
        Self::with_parse_flags(CtabParseFlags::BASIC_MAX)
    }

    pub fn strict() -> Self {
        Self::with_parse_flags(CtabParseFlags::STRICT)
    }

    pub fn extended() -> Self {
        Self::with_parse_flags(CtabParseFlags::EXTENDED)
    }

    pub fn extended_strict() -> Self {
        Self::with_parse_flags(CtabParseFlags::EXTENDED & CtabParseFlags::STRICT)
    }

    pub fn lenient() -> Self {
        Self::with_parse_flags(CtabParseFlags::LENIENT)
    }

    pub fn extended_lenient() -> Self {
        Self::with_parse_flags(CtabParseFlags::EXTENDED_MAX & CtabParseFlags::LENIENT)
    }

    pub fn extended_max() -> Self {
        Self::with_parse_flags(CtabParseFlags::EXTENDED_MAX)
    }
}

impl fmt::Display for CtfileIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "CtfileIoConfig({})", self.parse_flags)
    }
}
