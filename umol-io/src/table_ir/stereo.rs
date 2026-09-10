//! Stereochemistry records and metadata for TableIR.

/// Whether the molecule's stereo descriptors fix the absolute configuration or
/// only the relative one. Populated from format-specific flags:
/// - CTFile counts chiral flag (`ccc`)
/// - CXSMILES enhanced stereo markers (`a:` and `r`)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConfigurationScope {
    Absolute,
    Relative,
}

/// An explicit tetrahedral configuration in an ordered ligand frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoAtom {
    pub atom: u32,
    pub ligands: Vec<StereoLigand>,
    pub winding: Winding,
}

/// An actual neighbor or a virtual ligand borne by a stereo atom's site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoLigand {
    Atom(u32),
    ImplicitHydrogen,
    LonePair,
}

/// Winding of the last three tetrahedral ligands with the first toward the viewer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Winding {
    Clockwise,
    CounterClockwise,
}
