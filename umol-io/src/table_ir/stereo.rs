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

/// The frame in which a per-atom chirality descriptor is read into a 3D
/// arrangement. It governs tetrahedral atom chirality only, not other
/// stereogenic elements (e.g. E/Z bonds). It is present only when the molecule
/// contains a raw atom chirality descriptor whose source convention must be
/// retained.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChiralityFrame {
    /// First-listed neighbor points toward the viewer; remaining neighbors,
    /// in order, wind counterclockwise for the negative token (SMILES `@`).
    FirstNeighborToward,
    /// CTfile atom parity: neighbors numbered by increasing atom number with a
    /// hydrogen last, viewed with the last neighbor away; parity 2
    /// (counterclockwise) is the configuration SMILES writes as `@` over the
    /// same order. Retained as parsed and not read by the raise, which follows
    /// the specification's "ignored when read".
    LastNeighborAway,
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
