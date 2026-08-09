//! AST lowering/raising configuration: mode enums and config structs.

use umol_edn::{FromEdn, ToEdn};

/// Aggregated lowering/raising defaults for molecule DSL <-> AST interconversion. per-entity-kind
/// defaults bundle; consumed by the molecule-level `FromIr` / `IntoIr` implementations.
#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
pub struct MoleculeDefaults {
    pub atom: AtomDefaults,
    pub bond: BondDefaults,
    pub dative_bond: DativeBondDefaults,
    pub aromatic_system: AromaticSystemDefaults,
    pub multicenter_bond: MulticenterBondDefaults,
    pub noncovalent_bond: NoncovalentBondDefaults,
    pub stereo_atom: StereoAtomDefaults,
    pub stereo_bond: StereoBondDefaults,
}

impl MoleculeDefaults {
    /// Requires every configurable entity field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            atom: AtomDefaults::new(),
            bond: BondDefaults::new(),
            dative_bond: DativeBondDefaults::new(),
            aromatic_system: AromaticSystemDefaults::new(),
            multicenter_bond: MulticenterBondDefaults::new(),
            noncovalent_bond: NoncovalentBondDefaults::new(),
            stereo_atom: StereoAtomDefaults::new(),
            stereo_bond: StereoBondDefaults::new(),
        }
    }

    /// Composes `*Defaults::ground()` for each entity.
    pub fn ground() -> Self {
        Self {
            atom: AtomDefaults::ground(),
            bond: BondDefaults::ground(),
            dative_bond: DativeBondDefaults::ground(),
            aromatic_system: AromaticSystemDefaults::ground(),
            multicenter_bond: MulticenterBondDefaults::ground(),
            noncovalent_bond: NoncovalentBondDefaults::ground(),
            stereo_atom: StereoAtomDefaults::ground(),
            stereo_bond: StereoBondDefaults::ground(),
        }
    }

    /// Composes `*Defaults::zeroed()` for each entity. Atom topology-derived
    /// fields (`valence`, `donated_pairs`, `accepted_pairs`, `aromatic_valence`,
    /// `multicenter_valence`) are skipped during `into_ir` when the molecule
    /// has incident topology — see `MoleculeDsl::into_ir`.
    pub fn zeroed() -> Self {
        Self {
            atom: AtomDefaults::zeroed(),
            bond: BondDefaults::zeroed(),
            dative_bond: DativeBondDefaults::zeroed(),
            aromatic_system: AromaticSystemDefaults::zeroed(),
            multicenter_bond: MulticenterBondDefaults::zeroed(),
            noncovalent_bond: NoncovalentBondDefaults::zeroed(),
            stereo_atom: StereoAtomDefaults::zeroed(),
            stereo_bond: StereoBondDefaults::zeroed(),
        }
    }

    /// Add overrides.
    pub fn with_overrides(self, ov: MoleculeOverrides) -> Self {
        Self {
            atom: self.atom.with_overrides(ov.atom),
            bond: self.bond.with_overrides(ov.bond),
            dative_bond: self.dative_bond.with_overrides(ov.dative_bond),
            aromatic_system: self.aromatic_system.with_overrides(ov.aromatic_system),
            multicenter_bond: self.multicenter_bond.with_overrides(ov.multicenter_bond),
            noncovalent_bond: self.noncovalent_bond.with_overrides(ov.noncovalent_bond),
            stereo_atom: self.stereo_atom.with_overrides(ov.stereo_atom),
            stereo_bond: self.stereo_bond.with_overrides(ov.stereo_bond),
        }
    }
}

impl Default for MoleculeDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `MoleculeDefaults`. Each field is the
/// corresponding per-entity `*Overrides` bundle.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct MoleculeOverrides {
    #[edn(default)]
    pub atom: AtomOverrides,
    #[edn(default)]
    pub bond: BondOverrides,
    #[edn(default)]
    pub dative_bond: DativeBondOverrides,
    #[edn(default)]
    pub aromatic_system: AromaticSystemOverrides,
    #[edn(default)]
    pub multicenter_bond: MulticenterBondOverrides,
    #[edn(default)]
    pub noncovalent_bond: NoncovalentBondOverrides,
    #[edn(default)]
    pub stereo_atom: StereoAtomOverrides,
    #[edn(default)]
    pub stereo_bond: StereoBondOverrides,
}

/// Defaults for reaction DSL <-> AST interconversion.
#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
pub struct ReactionDefaults {
    pub atom: AtomDefaults,
    pub bond: BondDefaults,
    pub dative_bond: DativeBondDefaults,
    pub aromatic_system: AromaticSystemDefaults,
    pub multicenter_bond: MulticenterBondDefaults,
    pub noncovalent_bond: NoncovalentBondDefaults,
    pub stereo_atom: StereoAtomDefaults,
    pub stereo_bond: StereoBondDefaults,
}

impl ReactionDefaults {
    /// Requires every configurable entity field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            atom: AtomDefaults::new(),
            bond: BondDefaults::new(),
            dative_bond: DativeBondDefaults::new(),
            aromatic_system: AromaticSystemDefaults::new(),
            multicenter_bond: MulticenterBondDefaults::new(),
            noncovalent_bond: NoncovalentBondDefaults::new(),
            stereo_atom: StereoAtomDefaults::new(),
            stereo_bond: StereoBondDefaults::new(),
        }
    }

    /// Composes `*Defaults::ground()` for each entity.
    pub fn ground() -> Self {
        Self {
            atom: AtomDefaults::ground(),
            bond: BondDefaults::ground(),
            dative_bond: DativeBondDefaults::ground(),
            aromatic_system: AromaticSystemDefaults::ground(),
            multicenter_bond: MulticenterBondDefaults::ground(),
            noncovalent_bond: NoncovalentBondDefaults::ground(),
            stereo_atom: StereoAtomDefaults::ground(),
            stereo_bond: StereoBondDefaults::ground(),
        }
    }

    /// Add overrides.
    pub fn with_overrides(self, ov: ReactionOverrides) -> Self {
        Self {
            atom: self.atom.with_overrides(ov.atom),
            bond: self.bond.with_overrides(ov.bond),
            dative_bond: self.dative_bond.with_overrides(ov.dative_bond),
            aromatic_system: self.aromatic_system.with_overrides(ov.aromatic_system),
            multicenter_bond: self.multicenter_bond.with_overrides(ov.multicenter_bond),
            noncovalent_bond: self.noncovalent_bond.with_overrides(ov.noncovalent_bond),
            stereo_atom: self.stereo_atom.with_overrides(ov.stereo_atom),
            stereo_bond: self.stereo_bond.with_overrides(ov.stereo_bond),
        }
    }

    /// Molecule-level defaults for converting the `lhs`.
    pub fn molecule_defaults(&self) -> MoleculeDefaults {
        MoleculeDefaults {
            atom: self.atom.clone(),
            bond: self.bond.clone(),
            dative_bond: self.dative_bond.clone(),
            aromatic_system: self.aromatic_system.clone(),
            multicenter_bond: self.multicenter_bond.clone(),
            noncovalent_bond: self.noncovalent_bond.clone(),
            stereo_atom: self.stereo_atom.clone(),
            stereo_bond: self.stereo_bond.clone(),
        }
    }

    /// Defaults for converting delta entity snapshots.
    pub fn delta_defaults(&self) -> DeltaDefaults {
        DeltaDefaults {
            atom: self.atom.clone(),
            bond: self.bond.clone(),
            dative_bond: self.dative_bond.clone(),
            aromatic_system: self.aromatic_system.clone(),
            multicenter_bond: self.multicenter_bond.clone(),
            noncovalent_bond: self.noncovalent_bond.clone(),
            stereo_atom: self.stereo_atom.clone(),
            stereo_bond: self.stereo_bond.clone(),
        }
    }
}

impl Default for ReactionDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `ReactionDefaults`. Each field is the corresponding
/// per-entity `*Overrides` bundle.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct ReactionOverrides {
    #[edn(default)]
    pub atom: AtomOverrides,
    #[edn(default)]
    pub bond: BondOverrides,
    #[edn(default)]
    pub dative_bond: DativeBondOverrides,
    #[edn(default)]
    pub aromatic_system: AromaticSystemOverrides,
    #[edn(default)]
    pub multicenter_bond: MulticenterBondOverrides,
    #[edn(default)]
    pub noncovalent_bond: NoncovalentBondOverrides,
    #[edn(default)]
    pub stereo_atom: StereoAtomOverrides,
    #[edn(default)]
    pub stereo_bond: StereoBondOverrides,
}

/// Defaults for converting reaction delta entity snapshots.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeltaDefaults {
    pub atom: AtomDefaults,
    pub bond: BondDefaults,
    pub dative_bond: DativeBondDefaults,
    pub aromatic_system: AromaticSystemDefaults,
    pub multicenter_bond: MulticenterBondDefaults,
    pub noncovalent_bond: NoncovalentBondDefaults,
    pub stereo_atom: StereoAtomDefaults,
    pub stereo_bond: StereoBondDefaults,
}

/// Lowering/raising defaults for atoms: describe how `AtomForm`
/// struct fields and constraintsa are treated when converting between DSL and AST.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct AtomDefaults {
    pub isotope: IsotopeDefault,
    pub charge: NumericDefault,
    pub implicit_hydrogens: NumericDefault,
    pub lone_pairs: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub valence: NumericDefault,
    pub donated_pairs: NumericDefault,
    pub accepted_pairs: NumericDefault,
    pub multicenter_valence: MulticenterValenceDefault,
    pub aromatic_valence: AromaticValenceDefault,
    pub tetrahedral_stereo: StereoDefault,
}

impl AtomDefaults {
    /// Requires every atom field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            isotope: IsotopeDefault::Required,
            charge: NumericDefault::Required,
            implicit_hydrogens: NumericDefault::Required,
            lone_pairs: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
            valence: NumericDefault::Required,
            donated_pairs: NumericDefault::Required,
            accepted_pairs: NumericDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
            tetrahedral_stereo: StereoDefault::Required,
        }
    }

    /// Grounds struct fields, no constraints
    pub fn ground() -> Self {
        Self {
            isotope: IsotopeDefault::Natural,
            charge: NumericDefault::Zero,
            implicit_hydrogens: NumericDefault::Zero,
            lone_pairs: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumericDefault::Required,
            donated_pairs: NumericDefault::Required,
            accepted_pairs: NumericDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
            tetrahedral_stereo: StereoDefault::Required,
        }
    }

    /// Grounds struct fields and sets the constraint defaults to their zero values.
    /// Used as the omission threshold when *lowering* an AST to DSL (a constraint
    /// equal to its zeroed value is omitted, keeping output compact). Raising
    /// (DSL → AST) uses `ground()` plus an explicit per-constraint selection.
    pub fn zeroed() -> Self {
        Self {
            isotope: IsotopeDefault::Natural,
            charge: NumericDefault::Zero,
            implicit_hydrogens: NumericDefault::Zero,
            lone_pairs: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumericDefault::Zero,
            donated_pairs: NumericDefault::Zero,
            accepted_pairs: NumericDefault::Zero,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
            tetrahedral_stereo: StereoDefault::NotStereo,
        }
    }

    /// Add overrides
    pub fn with_overrides(mut self, ov: AtomOverrides) -> Self {
        if let Some(v) = ov.isotope {
            self.isotope = v;
        }
        if let Some(v) = ov.charge {
            self.charge = v;
        }
        if let Some(v) = ov.implicit_hydrogens {
            self.implicit_hydrogens = v;
        }
        if let Some(v) = ov.lone_pairs {
            self.lone_pairs = v;
        }
        if let Some(v) = ov.unpaired_electrons {
            self.unpaired_electrons = v;
        }
        if let Some(v) = ov.multiplicity {
            self.multiplicity = v;
        }
        if let Some(v) = ov.valence {
            self.valence = v;
        }
        if let Some(v) = ov.donated_pairs {
            self.donated_pairs = v;
        }
        if let Some(v) = ov.accepted_pairs {
            self.accepted_pairs = v;
        }
        if let Some(v) = ov.multicenter_valence {
            self.multicenter_valence = v;
        }
        if let Some(v) = ov.aromatic_valence {
            self.aromatic_valence = v;
        }
        if let Some(v) = ov.tetrahedral_stereo {
            self.tetrahedral_stereo = v;
        }
        self
    }
}

impl Default for AtomDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `AtomDefaults`. Fields set to `Some(..)` replace
/// the corresponding `AtomDefaults` field; `None` leaves it unchanged.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct AtomOverrides {
    pub isotope: Option<IsotopeDefault>,
    pub charge: Option<NumericDefault>,
    pub implicit_hydrogens: Option<NumericDefault>,
    pub lone_pairs: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
    pub valence: Option<NumericDefault>,
    pub donated_pairs: Option<NumericDefault>,
    pub accepted_pairs: Option<NumericDefault>,
    pub multicenter_valence: Option<MulticenterValenceDefault>,
    pub aromatic_valence: Option<AromaticValenceDefault>,
    pub tetrahedral_stereo: Option<StereoDefault>,
}

/// Lowering/raising defaults for localized bonds.
/// See `AtomDefaults` for semantics.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct BondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub cis_trans_stereo: StereoDefault,
}

impl BondDefaults {
    /// Requires every bond field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
            cis_trans_stereo: StereoDefault::Required,
        }
    }

    /// Grounds all struct fields
    pub fn ground() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            cis_trans_stereo: StereoDefault::Required,
        }
    }

    /// Like `ground()` but additionally sets `cis_trans_stereo` to `NotStereo`.
    /// The omission threshold for *lowering* a bond to DSL (compact output).
    pub fn zeroed() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            cis_trans_stereo: StereoDefault::NotStereo,
        }
    }

    /// Add overrides
    pub fn with_overrides(mut self, ov: BondOverrides) -> Self {
        if let Some(v) = ov.charge {
            self.charge = v;
        }
        if let Some(v) = ov.unpaired_electrons {
            self.unpaired_electrons = v;
        }
        if let Some(v) = ov.multiplicity {
            self.multiplicity = v;
        }
        self
    }
}

impl Default for BondDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `BondDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct BondOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

/// Lowering/raising defaults for dative bonds. Currently empty (no defaultable fields).
#[derive(Clone, Debug, Default, PartialEq, Eq, FromEdn, ToEdn)]
pub struct DativeBondDefaults {}

impl DativeBondDefaults {
    pub fn new() -> Self {
        Self {}
    }

    pub fn ground() -> Self {
        Self {}
    }

    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn with_overrides(self, _ov: DativeBondOverrides) -> Self {
        self
    }
}

/// Sparse overrides on `DativeBondDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct DativeBondOverrides {}

/// Lowering/raising defaults for aromatic systems.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct AromaticSystemDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl AromaticSystemDefaults {
    /// Requires every aromatic-system field to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
        }
    }

    /// Grounds all struct fields
    pub fn ground() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
        }
    }

    /// Equivalent to ground()
    pub fn zeroed() -> Self {
        Self::ground()
    }

    /// Add overrides
    pub fn with_overrides(mut self, ov: AromaticSystemOverrides) -> Self {
        if let Some(v) = ov.charge {
            self.charge = v;
        }
        if let Some(v) = ov.unpaired_electrons {
            self.unpaired_electrons = v;
        }
        if let Some(v) = ov.multiplicity {
            self.multiplicity = v;
        }
        self
    }
}

impl Default for AromaticSystemDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `AromaticSystemDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct AromaticSystemOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

/// Lowering/raising defaults for multicenter bonds.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct MulticenterBondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl MulticenterBondDefaults {
    /// Requires every multicenter-bond field to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
        }
    }

    /// Grounds all struct fields.
    pub fn ground() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
        }
    }

    /// Equivalent to ground()
    pub fn zeroed() -> Self {
        Self::ground()
    }

    /// Add overrides
    pub fn with_overrides(mut self, ov: MulticenterBondOverrides) -> Self {
        if let Some(v) = ov.charge {
            self.charge = v;
        }
        if let Some(v) = ov.unpaired_electrons {
            self.unpaired_electrons = v;
        }
        if let Some(v) = ov.multiplicity {
            self.multiplicity = v;
        }
        self
    }
}

impl Default for MulticenterBondDefaults {
    fn default() -> Self {
        Self::new()
    }
}

/// Sparse overrides on `MulticenterBondDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct MulticenterBondOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

/// Lowering/raising defaults for noncovalent bonds. Currently empty
/// (no defaultable fields); exists for API uniformity.
#[derive(Clone, Debug, Default, PartialEq, Eq, FromEdn, ToEdn)]
pub struct NoncovalentBondDefaults {}

impl NoncovalentBondDefaults {
    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn ground() -> Self {
        Self {}
    }

    pub fn new() -> Self {
        Self {}
    }

    pub fn with_overrides(self, _ov: NoncovalentBondOverrides) -> Self {
        self
    }
}

/// Sparse overrides on `NoncovalentBondDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct NoncovalentBondOverrides {}

/// Lowering/raising default for stereo atoms. Currently empty (no defaultable fields).
#[derive(Clone, Debug, Default, PartialEq, Eq, FromEdn, ToEdn)]
pub struct StereoAtomDefaults {}

impl StereoAtomDefaults {
    pub fn new() -> Self {
        Self {}
    }

    pub fn ground() -> Self {
        Self {}
    }

    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn with_overrides(self, _ov: StereoAtomOverrides) -> Self {
        self
    }
}

/// Sparse overrides on `StereoAtomDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct StereoAtomOverrides {}

/// Lowering/raising default for stereo bonds. Currently empty (no defaultable fields).
#[derive(Clone, Debug, Default, PartialEq, Eq, FromEdn, ToEdn)]
pub struct StereoBondDefaults {}

impl StereoBondDefaults {
    pub fn new() -> Self {
        Self {}
    }

    pub fn ground() -> Self {
        Self {}
    }

    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn with_overrides(self, _ov: StereoBondOverrides) -> Self {
        self
    }
}

/// Sparse overrides on `StereoBondDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct StereoBondOverrides {}

/// Isotope default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum IsotopeDefault {
    Natural,
    Required,
}

/// Numeric field default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum NumericDefault {
    Zero,
    Required,
}

/// Unpaired electrons default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum UnpairedElectronsDefault {
    Zero,
    Required,
    Derived,
}

/// Spin multiplicity default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum MultiplicityDefault {
    Derived,
    Required,
}

/// Aromatic valence default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum AromaticValenceDefault {
    NotAromatic,
    Required,
}

/// Multicenter valence default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum MulticenterValenceDefault {
    NotMulticenter,
    Required,
}

/// Stereo default, shared by tetrahedral chirality and cis/trans constraints.
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum StereoDefault {
    NotStereo,
    Required,
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_edn::{read_string, FromEdn, ToEdn};

    use super::*;

    #[rstest]
    #[case::verbatim(MoleculeDefaults::new())]
    #[case::ground(MoleculeDefaults::ground())]
    #[case::zeroed(MoleculeDefaults::zeroed())]
    fn test_molecule_ast_config_roundtrip(#[case] cfg: MoleculeDefaults) {
        let edn = cfg.to_edn();
        let back = MoleculeDefaults::from_edn(&edn).unwrap();
        assert_eq!(back, cfg);
    }

    #[rstest]
    #[case::atom(
        AtomDefaults::new().unpaired_electrons,
        AtomDefaults::new().multiplicity
    )]
    #[case::bond(
        BondDefaults::new().unpaired_electrons,
        BondDefaults::new().multiplicity
    )]
    #[case::aromatic_system(
        AromaticSystemDefaults::new().unpaired_electrons,
        AromaticSystemDefaults::new().multiplicity
    )]
    #[case::multicenter_bond(
        MulticenterBondDefaults::new().unpaired_electrons,
        MulticenterBondDefaults::new().multiplicity
    )]
    fn test_new_spin_defaults_required(
        #[case] unpaired_electrons: UnpairedElectronsDefault,
        #[case] multiplicity: MultiplicityDefault,
    ) {
        assert_eq!(unpaired_electrons, UnpairedElectronsDefault::Required);
        assert_eq!(multiplicity, MultiplicityDefault::Required);
    }

    #[rstest]
    fn test_molecule_defaults_with_overrides_routes_to_per_entity() {
        let cfg = MoleculeDefaults::ground().with_overrides(MoleculeOverrides {
            atom: AtomOverrides {
                charge: Some(NumericDefault::Required),
                ..AtomOverrides::default()
            },
            bond: BondOverrides {
                multiplicity: Some(MultiplicityDefault::Required),
                ..BondOverrides::default()
            },
            dative_bond: DativeBondOverrides::default(),
            aromatic_system: AromaticSystemOverrides {
                charge: Some(NumericDefault::Required),
                ..AromaticSystemOverrides::default()
            },
            multicenter_bond: MulticenterBondOverrides {
                charge: Some(NumericDefault::Required),
                ..MulticenterBondOverrides::default()
            },
            noncovalent_bond: NoncovalentBondOverrides::default(),
            stereo_atom: StereoAtomOverrides::default(),
            stereo_bond: StereoBondOverrides::default(),
        });
        assert_eq!(cfg.atom.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.atom.charge, NumericDefault::Required);
        assert_eq!(cfg.bond.charge, NumericDefault::Zero);
        assert_eq!(cfg.bond.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.aromatic_system.charge, NumericDefault::Required);
        assert_eq!(cfg.multicenter_bond.charge, NumericDefault::Required);
    }

    #[rstest]
    #[case::required(ReactionDefaults::new())]
    #[case::ground(ReactionDefaults::ground())]
    fn test_reaction_defaults_roundtrip(#[case] defaults: ReactionDefaults) {
        assert_eq!(
            ReactionDefaults::from_edn(&defaults.to_edn()).unwrap(),
            defaults
        );
    }

    #[rstest]
    #[case::required(ReactionDefaults::new(), MoleculeDefaults::new())]
    #[case::ground(ReactionDefaults::ground(), MoleculeDefaults::ground())]
    fn test_reaction_defaults_molecule_defaults(
        #[case] defaults: ReactionDefaults,
        #[case] expected: MoleculeDefaults,
    ) {
        assert_eq!(defaults.molecule_defaults(), expected);
    }

    #[rstest]
    #[case::required(
        ReactionDefaults::new(),
        DeltaDefaults {
            atom: AtomDefaults::new(),
            bond: BondDefaults::new(),
            dative_bond: DativeBondDefaults::new(),
            aromatic_system: AromaticSystemDefaults::new(),
            multicenter_bond: MulticenterBondDefaults::new(),
            noncovalent_bond: NoncovalentBondDefaults::new(),
            stereo_atom: StereoAtomDefaults::new(),
            stereo_bond: StereoBondDefaults::new(),
        }
    )]
    #[case::ground(
        ReactionDefaults::ground(),
        DeltaDefaults {
            atom: AtomDefaults::ground(),
            bond: BondDefaults::ground(),
            dative_bond: DativeBondDefaults::ground(),
            aromatic_system: AromaticSystemDefaults::ground(),
            multicenter_bond: MulticenterBondDefaults::ground(),
            noncovalent_bond: NoncovalentBondDefaults::ground(),
            stereo_atom: StereoAtomDefaults::ground(),
            stereo_bond: StereoBondDefaults::ground(),
        }
    )]
    fn test_reaction_defaults_delta_defaults(
        #[case] defaults: ReactionDefaults,
        #[case] expected: DeltaDefaults,
    ) {
        assert_eq!(defaults.delta_defaults(), expected);
    }

    #[rstest]
    #[case::explicit_derived_spin(
        "{:isotope :required :charge :required :implicit-hydrogens :required \
         :lone-pairs :required :unpaired-electrons :derived :multiplicity :derived \
         :valence :required :donated-pairs :required :accepted-pairs :required \
         :multicenter-valence :required :aromatic-valence :required \
         :tetrahedral-stereo :required}",
        NumericDefault::Required,
        NumericDefault::Required,
        MulticenterValenceDefault::Required,
        AromaticValenceDefault::Required
    )]
    #[case::zeroed_atom(
        "{:isotope :natural :charge :zero :implicit-hydrogens :zero \
         :lone-pairs :zero :unpaired-electrons :zero :multiplicity :derived \
         :valence :zero :donated-pairs :zero :accepted-pairs :zero \
         :multicenter-valence :not-multicenter :aromatic-valence :not-aromatic \
         :tetrahedral-stereo :not-stereo}",
        NumericDefault::Zero,
        NumericDefault::Zero,
        MulticenterValenceDefault::NotMulticenter,
        AromaticValenceDefault::NotAromatic
    )]
    fn test_atom_form_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumericDefault,
        #[case] expected_h: NumericDefault,
        #[case] expected_multicenter: MulticenterValenceDefault,
        #[case] expected_aromatic: AromaticValenceDefault,
    ) {
        let tree = read_string(edn).unwrap();
        let cfg = AtomDefaults::from_edn(&tree).unwrap();
        assert_eq!(cfg.charge, expected_charge);
        assert_eq!(cfg.implicit_hydrogens, expected_h);
        assert_eq!(cfg.multicenter_valence, expected_multicenter);
        assert_eq!(cfg.aromatic_valence, expected_aromatic);
    }

    #[rstest]
    fn test_atom_defaults_ground_constraints_required() {
        let g = AtomDefaults::ground();
        assert_eq!(g.valence, NumericDefault::Required);
        assert_eq!(g.donated_pairs, NumericDefault::Required);
        assert_eq!(g.accepted_pairs, NumericDefault::Required);
        assert_eq!(g.multicenter_valence, MulticenterValenceDefault::Required);
        assert_eq!(g.aromatic_valence, AromaticValenceDefault::Required);
    }

    #[rstest]
    fn test_atom_defaults_ground_struct_fields_match_zeroed() {
        let g = AtomDefaults::ground();
        let z = AtomDefaults::zeroed();
        assert_eq!(g.isotope, z.isotope);
        assert_eq!(g.charge, z.charge);
        assert_eq!(g.implicit_hydrogens, z.implicit_hydrogens);
        assert_eq!(g.lone_pairs, z.lone_pairs);
        assert_eq!(g.unpaired_electrons, z.unpaired_electrons);
        assert_eq!(g.multiplicity, z.multiplicity);
    }

    #[rstest]
    fn test_atom_defaults_with_overrides() {
        let cfg = AtomDefaults::zeroed().with_overrides(AtomOverrides {
            isotope: Some(IsotopeDefault::Required),
            charge: Some(NumericDefault::Required),
            implicit_hydrogens: Some(NumericDefault::Required),
            lone_pairs: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
            valence: Some(NumericDefault::Required),
            donated_pairs: Some(NumericDefault::Required),
            accepted_pairs: Some(NumericDefault::Required),
            multicenter_valence: Some(MulticenterValenceDefault::Required),
            aromatic_valence: Some(AromaticValenceDefault::Required),
            tetrahedral_stereo: Some(StereoDefault::Required),
        });
        assert_eq!(cfg.isotope, IsotopeDefault::Required);
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.implicit_hydrogens, NumericDefault::Required);
        assert_eq!(cfg.lone_pairs, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.valence, NumericDefault::Required);
        assert_eq!(cfg.donated_pairs, NumericDefault::Required);
        assert_eq!(cfg.accepted_pairs, NumericDefault::Required);
        assert_eq!(cfg.multicenter_valence, MulticenterValenceDefault::Required);
        assert_eq!(cfg.aromatic_valence, AromaticValenceDefault::Required);
        assert_eq!(cfg.tetrahedral_stereo, StereoDefault::Required);
    }

    #[rstest]
    fn test_atom_defaults_with_overrides_partial() {
        let cfg = AtomDefaults::zeroed().with_overrides(AtomOverrides {
            charge: Some(NumericDefault::Required),
            ..AtomOverrides::default()
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        // Untouched fields retain the zeroed() defaults.
        assert_eq!(cfg.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.implicit_hydrogens, NumericDefault::Zero);
        assert_eq!(cfg.valence, NumericDefault::Zero);
        assert_eq!(cfg.aromatic_valence, AromaticValenceDefault::NotAromatic);
    }

    #[rstest]
    fn test_bond_defaults_with_overrides() {
        let cfg = BondDefaults::zeroed().with_overrides(BondOverrides {
            charge: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }

    #[rstest]
    fn test_aromatic_system_defaults_with_overrides() {
        let cfg = AromaticSystemDefaults::zeroed().with_overrides(AromaticSystemOverrides {
            charge: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }

    #[rstest]
    fn test_multicenter_bond_defaults_with_overrides() {
        let cfg = MulticenterBondDefaults::zeroed().with_overrides(MulticenterBondOverrides {
            charge: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }
}
