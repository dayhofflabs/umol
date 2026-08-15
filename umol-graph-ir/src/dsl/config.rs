//! Graph-IR lowering/raising configuration: mode enums and config structs.

use umol_edn::{FromEdn, ToEdn};

/// Aggregated lowering/raising defaults for molecule DSL <-> IR interconversion. Per-entity-kind
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

    /// Composes `*Defaults::concrete()` for each entity.
    pub fn concrete() -> Self {
        Self {
            atom: AtomDefaults::concrete(),
            bond: BondDefaults::concrete(),
            dative_bond: DativeBondDefaults::concrete(),
            aromatic_system: AromaticSystemDefaults::concrete(),
            multicenter_bond: MulticenterBondDefaults::concrete(),
            noncovalent_bond: NoncovalentBondDefaults::concrete(),
            stereo_atom: StereoAtomDefaults::concrete(),
            stereo_bond: StereoBondDefaults::concrete(),
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

/// Defaults for reaction DSL <-> IR interconversion.
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

    /// Composes `*Defaults::concrete()` for each entity.
    pub fn concrete() -> Self {
        Self {
            atom: AtomDefaults::concrete(),
            bond: BondDefaults::concrete(),
            dative_bond: DativeBondDefaults::concrete(),
            aromatic_system: AromaticSystemDefaults::concrete(),
            multicenter_bond: MulticenterBondDefaults::concrete(),
            noncovalent_bond: NoncovalentBondDefaults::concrete(),
            stereo_atom: StereoAtomDefaults::concrete(),
            stereo_bond: StereoBondDefaults::concrete(),
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
/// struct fields and constraints are treated when converting between DSL and IR.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct AtomDefaults {
    pub isotope: IsotopeDefault,
    pub charge: NumDefault,
    pub implicit_hydrogens: NumDefault,
    pub lone_pairs: NumDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub valence: NumDefault,
    pub donated_pairs: NumDefault,
    pub accepted_pairs: NumDefault,
    pub multicenter_valence: MulticenterValenceDefault,
    pub aromatic_valence: AromaticValenceDefault,
    pub tetrahedral_stereo: StereoDefault,
}

impl AtomDefaults {
    /// Requires every atom field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            isotope: IsotopeDefault::Required,
            charge: NumDefault::Required,
            implicit_hydrogens: NumDefault::Required,
            lone_pairs: NumDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
            valence: NumDefault::Required,
            donated_pairs: NumDefault::Required,
            accepted_pairs: NumDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
            tetrahedral_stereo: StereoDefault::Required,
        }
    }

    /// Concrete struct-field defaults; constraints stay required.
    pub fn concrete() -> Self {
        Self {
            isotope: IsotopeDefault::Natural,
            charge: NumDefault::Zero,
            implicit_hydrogens: NumDefault::Zero,
            lone_pairs: NumDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumDefault::Required,
            donated_pairs: NumDefault::Required,
            accepted_pairs: NumDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
            tetrahedral_stereo: StereoDefault::Required,
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
    pub charge: Option<NumDefault>,
    pub implicit_hydrogens: Option<NumDefault>,
    pub lone_pairs: Option<NumDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
    pub valence: Option<NumDefault>,
    pub donated_pairs: Option<NumDefault>,
    pub accepted_pairs: Option<NumDefault>,
    pub multicenter_valence: Option<MulticenterValenceDefault>,
    pub aromatic_valence: Option<AromaticValenceDefault>,
    pub tetrahedral_stereo: Option<StereoDefault>,
}

/// Lowering/raising defaults for localized bonds.
/// See `AtomDefaults` for semantics.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct BondDefaults {
    pub charge: NumDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub cis_trans_stereo: StereoDefault,
}

impl BondDefaults {
    /// Requires every bond field and constraint to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
            cis_trans_stereo: StereoDefault::Required,
        }
    }

    /// Concrete struct-field defaults.
    pub fn concrete() -> Self {
        Self {
            charge: NumDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            cis_trans_stereo: StereoDefault::Required,
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
    pub charge: Option<NumDefault>,
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

    pub fn concrete() -> Self {
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
    pub charge: NumDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl AromaticSystemDefaults {
    /// Requires every aromatic-system field to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
        }
    }

    /// Concrete struct-field defaults.
    pub fn concrete() -> Self {
        Self {
            charge: NumDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
        }
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
    pub charge: Option<NumDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

/// Lowering/raising defaults for multicenter bonds.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct MulticenterBondDefaults {
    pub charge: NumDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl MulticenterBondDefaults {
    /// Requires every multicenter-bond field to be explicit.
    pub fn new() -> Self {
        Self {
            charge: NumDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
        }
    }

    /// Concrete struct-field defaults.
    pub fn concrete() -> Self {
        Self {
            charge: NumDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
        }
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
    pub charge: Option<NumDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

/// Lowering/raising defaults for noncovalent bonds. Currently empty
/// (no defaultable fields); exists for API uniformity.
#[derive(Clone, Debug, Default, PartialEq, Eq, FromEdn, ToEdn)]
pub struct NoncovalentBondDefaults {}

impl NoncovalentBondDefaults {
    pub fn concrete() -> Self {
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

    pub fn concrete() -> Self {
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

    pub fn concrete() -> Self {
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
pub enum NumDefault {
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
    #[case::concrete(MoleculeDefaults::concrete())]
    fn test_molecule_config_roundtrip(#[case] cfg: MoleculeDefaults) {
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
        let cfg = MoleculeDefaults::concrete().with_overrides(MoleculeOverrides {
            atom: AtomOverrides {
                charge: Some(NumDefault::Required),
                ..AtomOverrides::default()
            },
            bond: BondOverrides {
                multiplicity: Some(MultiplicityDefault::Required),
                ..BondOverrides::default()
            },
            dative_bond: DativeBondOverrides::default(),
            aromatic_system: AromaticSystemOverrides {
                charge: Some(NumDefault::Required),
                ..AromaticSystemOverrides::default()
            },
            multicenter_bond: MulticenterBondOverrides {
                charge: Some(NumDefault::Required),
                ..MulticenterBondOverrides::default()
            },
            noncovalent_bond: NoncovalentBondOverrides::default(),
            stereo_atom: StereoAtomOverrides::default(),
            stereo_bond: StereoBondOverrides::default(),
        });
        assert_eq!(cfg.atom.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.atom.charge, NumDefault::Required);
        assert_eq!(cfg.bond.charge, NumDefault::Zero);
        assert_eq!(cfg.bond.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.aromatic_system.charge, NumDefault::Required);
        assert_eq!(cfg.multicenter_bond.charge, NumDefault::Required);
    }

    #[rstest]
    #[case::required(ReactionDefaults::new())]
    #[case::concrete(ReactionDefaults::concrete())]
    fn test_reaction_defaults_roundtrip(#[case] defaults: ReactionDefaults) {
        assert_eq!(
            ReactionDefaults::from_edn(&defaults.to_edn()).unwrap(),
            defaults
        );
    }

    #[rstest]
    #[case::required(ReactionDefaults::new(), MoleculeDefaults::new())]
    #[case::concrete(ReactionDefaults::concrete(), MoleculeDefaults::concrete())]
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
    #[case::concrete(
        ReactionDefaults::concrete(),
        DeltaDefaults {
            atom: AtomDefaults::concrete(),
            bond: BondDefaults::concrete(),
            dative_bond: DativeBondDefaults::concrete(),
            aromatic_system: AromaticSystemDefaults::concrete(),
            multicenter_bond: MulticenterBondDefaults::concrete(),
            noncovalent_bond: NoncovalentBondDefaults::concrete(),
            stereo_atom: StereoAtomDefaults::concrete(),
            stereo_bond: StereoBondDefaults::concrete(),
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
        NumDefault::Required,
        NumDefault::Required,
        MulticenterValenceDefault::Required,
        AromaticValenceDefault::Required
    )]
    #[case::constraint_defaulted_atom(
        "{:isotope :natural :charge :zero :implicit-hydrogens :zero \
         :lone-pairs :zero :unpaired-electrons :zero :multiplicity :derived \
         :valence :zero :donated-pairs :zero :accepted-pairs :zero \
         :multicenter-valence :not-multicenter :aromatic-valence :not-aromatic \
         :tetrahedral-stereo :not-stereo}",
        NumDefault::Zero,
        NumDefault::Zero,
        MulticenterValenceDefault::NotMulticenter,
        AromaticValenceDefault::NotAromatic
    )]
    fn test_atom_form_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumDefault,
        #[case] expected_h: NumDefault,
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
    fn test_atom_defaults_concrete_constraints_required() {
        let g = AtomDefaults::concrete();
        assert_eq!(g.valence, NumDefault::Required);
        assert_eq!(g.donated_pairs, NumDefault::Required);
        assert_eq!(g.accepted_pairs, NumDefault::Required);
        assert_eq!(g.multicenter_valence, MulticenterValenceDefault::Required);
        assert_eq!(g.aromatic_valence, AromaticValenceDefault::Required);
    }

    #[rstest]
    fn test_atom_defaults_with_overrides() {
        let cfg = AtomDefaults::concrete().with_overrides(AtomOverrides {
            isotope: Some(IsotopeDefault::Required),
            charge: Some(NumDefault::Required),
            implicit_hydrogens: Some(NumDefault::Required),
            lone_pairs: Some(NumDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
            valence: Some(NumDefault::Required),
            donated_pairs: Some(NumDefault::Required),
            accepted_pairs: Some(NumDefault::Required),
            multicenter_valence: Some(MulticenterValenceDefault::Required),
            aromatic_valence: Some(AromaticValenceDefault::Required),
            tetrahedral_stereo: Some(StereoDefault::Required),
        });
        assert_eq!(cfg.isotope, IsotopeDefault::Required);
        assert_eq!(cfg.charge, NumDefault::Required);
        assert_eq!(cfg.implicit_hydrogens, NumDefault::Required);
        assert_eq!(cfg.lone_pairs, NumDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.valence, NumDefault::Required);
        assert_eq!(cfg.donated_pairs, NumDefault::Required);
        assert_eq!(cfg.accepted_pairs, NumDefault::Required);
        assert_eq!(cfg.multicenter_valence, MulticenterValenceDefault::Required);
        assert_eq!(cfg.aromatic_valence, AromaticValenceDefault::Required);
        assert_eq!(cfg.tetrahedral_stereo, StereoDefault::Required);
    }

    #[rstest]
    fn test_atom_defaults_with_overrides_partial() {
        let cfg = AtomDefaults::concrete().with_overrides(AtomOverrides {
            charge: Some(NumDefault::Required),
            ..AtomOverrides::default()
        });
        assert_eq!(cfg.charge, NumDefault::Required);
        // Untouched fields retain the concrete() defaults.
        assert_eq!(cfg.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.implicit_hydrogens, NumDefault::Zero);
        assert_eq!(cfg.valence, NumDefault::Required);
        assert_eq!(cfg.aromatic_valence, AromaticValenceDefault::Required);
    }

    #[rstest]
    fn test_bond_defaults_with_overrides() {
        let cfg = BondDefaults::concrete().with_overrides(BondOverrides {
            charge: Some(NumDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }

    #[rstest]
    fn test_aromatic_system_defaults_with_overrides() {
        let cfg = AromaticSystemDefaults::concrete().with_overrides(AromaticSystemOverrides {
            charge: Some(NumDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }

    #[rstest]
    fn test_multicenter_bond_defaults_with_overrides() {
        let cfg = MulticenterBondDefaults::concrete().with_overrides(MulticenterBondOverrides {
            charge: Some(NumDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
        });
        assert_eq!(cfg.charge, NumDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
    }
}
