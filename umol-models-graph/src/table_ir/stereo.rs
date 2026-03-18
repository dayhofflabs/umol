//! Stereochemistry metadata for TableIR.

/// Molecule-wide stereochemistry interpretation context.
///
/// This is a small, typed signal that can be populated from format-specific flags:
/// - CTFile counts chiral flag (`ccc`)
/// - CXSMILES enhanced stereo markers (`a:` and `r`)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoInterpretation {
    Absolute,
    Relative,
}
