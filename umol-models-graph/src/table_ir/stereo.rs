//! Stereochemistry metadata for TableIR.

/// Molecule-wide stereochemistry interpretation context.
///
/// This is a small, typed signal that can be populated from format-specific flags:
/// - CTFile counts chiral flag (`ccc`)
/// - CXSMILES enhanced stereo markers (`a:` and `r`)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoInterpretation {
    /// Stereochemistry is intended to be absolute (a specific stereoisomer).
    Absolute,
    /// Stereochemistry is intended to be relative (relationships without absolute assignment).
    Relative,
}

