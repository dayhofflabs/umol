//! AST lowering/raising configuration: mode enums and config structs.

use umol_edn::{FromEdn, ToEdn};

/// Isotope interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum IsotopeMode {
    Natural,
    Required,
}

/// Numeric field interpretation mode.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum NumericMode {
    Zero,
    Required,
}

/// Implicit hydrogen interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum ImplicitHydrogenMode {
    Zero,
    Normal,
    Required,
}

/// Unpaired electrons interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum UnpairedElectronsMode {
    Zero,
    Required,
    Derived,
}

/// Multiplicity interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum MultiplicityMode {
    Derived,
    Required,
}

/// Aromatic interpretation mode (molecule-level: applied to atoms by injecting
/// `AromaticValence` constraints).
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum AromaticValenceMode {
    NotAromatic,
    Aromatic,
    Required,
}

/// Multicenter valence interpretation mode.
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum MulticenterValenceMode {
    NotMulticenter,
    Multicenter,
    Required,
}

#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct AtomAstConfig {
    pub isotope_mode: IsotopeMode,
    pub charge_mode: NumericMode,
    pub implicit_h_mode: ImplicitHydrogenMode,
    pub lone_pairs_mode: NumericMode,
    pub unpaired_electrons_mode: UnpairedElectronsMode,
    pub multiplicity_mode: MultiplicityMode,
    pub valence_mode: NumericMode,
    pub donated_pairs_mode: NumericMode,
    pub accepted_pairs_mode: NumericMode,
    pub multicenter_valence_mode: MulticenterValenceMode,
    pub aromatic_valence_mode: AromaticValenceMode,
}

impl AtomAstConfig {
    pub fn zeroed() -> Self {
        Self {
            isotope_mode: IsotopeMode::Natural,
            charge_mode: NumericMode::Zero,
            implicit_h_mode: ImplicitHydrogenMode::Zero,
            lone_pairs_mode: NumericMode::Zero,
            unpaired_electrons_mode: UnpairedElectronsMode::Zero,
            multiplicity_mode: MultiplicityMode::Derived,
            valence_mode: NumericMode::Zero,
            donated_pairs_mode: NumericMode::Zero,
            accepted_pairs_mode: NumericMode::Zero,
            multicenter_valence_mode: MulticenterValenceMode::NotMulticenter,
            aromatic_valence_mode: AromaticValenceMode::NotAromatic,
        }
    }

    pub fn open() -> Self {
        Self {
            isotope_mode: IsotopeMode::Required,
            charge_mode: NumericMode::Required,
            implicit_h_mode: ImplicitHydrogenMode::Required,
            lone_pairs_mode: NumericMode::Required,
            unpaired_electrons_mode: UnpairedElectronsMode::Derived,
            multiplicity_mode: MultiplicityMode::Derived,
            valence_mode: NumericMode::Required,
            donated_pairs_mode: NumericMode::Required,
            accepted_pairs_mode: NumericMode::Required,
            multicenter_valence_mode: MulticenterValenceMode::Required,
            aromatic_valence_mode: AromaticValenceMode::Required,
        }
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct AtomAstConfigOverrides {
    pub isotope_mode: Option<IsotopeMode>,
    pub charge_mode: Option<NumericMode>,
    pub implicit_h_mode: Option<ImplicitHydrogenMode>,
    pub lone_pairs_mode: Option<NumericMode>,
    pub unpaired_electrons_mode: Option<UnpairedElectronsMode>,
    pub multiplicity_mode: Option<MultiplicityMode>,
    pub valence_mode: Option<NumericMode>,
    pub donated_pairs_mode: Option<NumericMode>,
    pub accepted_pairs_mode: Option<NumericMode>,
    pub multicenter_valence_mode: Option<MulticenterValenceMode>,
    pub aromatic_valence_mode: Option<AromaticValenceMode>,
}

impl AtomAstConfig {
    pub fn with_overrides(mut self, ov: AtomAstConfigOverrides) -> Self {
        if let Some(v) = ov.isotope_mode {
            self.isotope_mode = v;
        }
        if let Some(v) = ov.charge_mode {
            self.charge_mode = v;
        }
        if let Some(v) = ov.implicit_h_mode {
            self.implicit_h_mode = v;
        }
        if let Some(v) = ov.lone_pairs_mode {
            self.lone_pairs_mode = v;
        }
        if let Some(v) = ov.unpaired_electrons_mode {
            self.unpaired_electrons_mode = v;
        }
        if let Some(v) = ov.multiplicity_mode {
            self.multiplicity_mode = v;
        }
        if let Some(v) = ov.valence_mode {
            self.valence_mode = v;
        }
        if let Some(v) = ov.donated_pairs_mode {
            self.donated_pairs_mode = v;
        }
        if let Some(v) = ov.accepted_pairs_mode {
            self.accepted_pairs_mode = v;
        }
        if let Some(v) = ov.multicenter_valence_mode {
            self.multicenter_valence_mode = v;
        }
        if let Some(v) = ov.aromatic_valence_mode {
            self.aromatic_valence_mode = v;
        }
        self
    }
}

#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct BondAstConfig {
    pub charge_mode: NumericMode,
    pub unpaired_electrons_mode: UnpairedElectronsMode,
    pub multiplicity_mode: MultiplicityMode,
}

impl BondAstConfig {
    pub fn zeroed() -> Self {
        Self {
            charge_mode: NumericMode::Zero,
            unpaired_electrons_mode: UnpairedElectronsMode::Zero,
            multiplicity_mode: MultiplicityMode::Derived,
        }
    }

    pub fn open() -> Self {
        Self {
            charge_mode: NumericMode::Required,
            unpaired_electrons_mode: UnpairedElectronsMode::Derived,
            multiplicity_mode: MultiplicityMode::Derived,
        }
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct BondAstConfigOverrides {
    pub charge_mode: Option<NumericMode>,
    pub unpaired_electrons_mode: Option<UnpairedElectronsMode>,
    pub multiplicity_mode: Option<MultiplicityMode>,
}

impl BondAstConfig {
    pub fn with_overrides(mut self, ov: BondAstConfigOverrides) -> Self {
        if let Some(v) = ov.charge_mode {
            self.charge_mode = v;
        }
        if let Some(v) = ov.unpaired_electrons_mode {
            self.unpaired_electrons_mode = v;
        }
        if let Some(v) = ov.multiplicity_mode {
            self.multiplicity_mode = v;
        }
        self
    }
}

#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct AromaticSystemAstConfig {
    pub charge_mode: NumericMode,
    pub unpaired_electrons_mode: UnpairedElectronsMode,
    pub multiplicity_mode: MultiplicityMode,
    pub electrons_mode: NumericMode,
}

impl AromaticSystemAstConfig {
    pub fn zeroed() -> Self {
        Self {
            charge_mode: NumericMode::Zero,
            unpaired_electrons_mode: UnpairedElectronsMode::Zero,
            multiplicity_mode: MultiplicityMode::Derived,
            electrons_mode: NumericMode::Zero,
        }
    }

    pub fn open() -> Self {
        Self {
            charge_mode: NumericMode::Required,
            unpaired_electrons_mode: UnpairedElectronsMode::Derived,
            multiplicity_mode: MultiplicityMode::Derived,
            electrons_mode: NumericMode::Required,
        }
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct AromaticSystemAstConfigOverrides {
    pub charge_mode: Option<NumericMode>,
    pub unpaired_electrons_mode: Option<UnpairedElectronsMode>,
    pub multiplicity_mode: Option<MultiplicityMode>,
    pub electrons_mode: Option<NumericMode>,
}

impl AromaticSystemAstConfig {
    pub fn with_overrides(mut self, ov: AromaticSystemAstConfigOverrides) -> Self {
        if let Some(v) = ov.charge_mode {
            self.charge_mode = v;
        }
        if let Some(v) = ov.unpaired_electrons_mode {
            self.unpaired_electrons_mode = v;
        }
        if let Some(v) = ov.multiplicity_mode {
            self.multiplicity_mode = v;
        }
        if let Some(v) = ov.electrons_mode {
            self.electrons_mode = v;
        }
        self
    }
}

#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct MulticenterBondAstConfig {
    pub charge_mode: NumericMode,
    pub unpaired_electrons_mode: UnpairedElectronsMode,
    pub multiplicity_mode: MultiplicityMode,
    pub electrons_mode: NumericMode,
}

impl MulticenterBondAstConfig {
    pub fn zeroed() -> Self {
        Self {
            charge_mode: NumericMode::Zero,
            unpaired_electrons_mode: UnpairedElectronsMode::Zero,
            multiplicity_mode: MultiplicityMode::Derived,
            electrons_mode: NumericMode::Zero,
        }
    }

    pub fn open() -> Self {
        Self {
            charge_mode: NumericMode::Required,
            unpaired_electrons_mode: UnpairedElectronsMode::Derived,
            multiplicity_mode: MultiplicityMode::Derived,
            electrons_mode: NumericMode::Required,
        }
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct MulticenterBondAstConfigOverrides {
    pub charge_mode: Option<NumericMode>,
    pub unpaired_electrons_mode: Option<UnpairedElectronsMode>,
    pub multiplicity_mode: Option<MultiplicityMode>,
    pub electrons_mode: Option<NumericMode>,
}

impl MulticenterBondAstConfig {
    pub fn with_overrides(mut self, ov: MulticenterBondAstConfigOverrides) -> Self {
        if let Some(v) = ov.charge_mode {
            self.charge_mode = v;
        }
        if let Some(v) = ov.unpaired_electrons_mode {
            self.unpaired_electrons_mode = v;
        }
        if let Some(v) = ov.multiplicity_mode {
            self.multiplicity_mode = v;
        }
        if let Some(v) = ov.electrons_mode {
            self.electrons_mode = v;
        }
        self
    }
}

#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
pub struct DativeBondAstConfig {}

impl DativeBondAstConfig {
    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn open() -> Self {
        Self {}
    }
}

#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
pub struct NoncovalentBondAstConfig {}

impl NoncovalentBondAstConfig {
    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn open() -> Self {
        Self {}
    }
}

#[derive(Debug, Clone, FromEdn, ToEdn)]
pub struct MoleculeAstConfig {
    pub atom: AtomAstConfig,
    pub bond: BondAstConfig,
    pub aromatic_system: AromaticSystemAstConfig,
    pub multicenter_bond: MulticenterBondAstConfig,
    pub dative_bond: DativeBondAstConfig,
    pub noncovalent_bond: NoncovalentBondAstConfig,
}

impl MoleculeAstConfig {
    pub fn zeroed() -> Self {
        Self {
            atom: AtomAstConfig::zeroed(),
            bond: BondAstConfig::zeroed(),
            aromatic_system: AromaticSystemAstConfig::zeroed(),
            multicenter_bond: MulticenterBondAstConfig::zeroed(),
            dative_bond: DativeBondAstConfig::zeroed(),
            noncovalent_bond: NoncovalentBondAstConfig::zeroed(),
        }
    }

    pub fn open() -> Self {
        Self {
            atom: AtomAstConfig::open(),
            bond: BondAstConfig::open(),
            aromatic_system: AromaticSystemAstConfig::open(),
            multicenter_bond: MulticenterBondAstConfig::open(),
            dative_bond: DativeBondAstConfig::open(),
            noncovalent_bond: NoncovalentBondAstConfig::open(),
        }
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct DativeBondAstConfigOverrides {}

impl DativeBondAstConfig {
    pub fn with_overrides(self, _ov: DativeBondAstConfigOverrides) -> Self {
        self
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct NoncovalentBondAstConfigOverrides {}

impl NoncovalentBondAstConfig {
    pub fn with_overrides(self, _ov: NoncovalentBondAstConfigOverrides) -> Self {
        self
    }
}

#[derive(Clone, Debug, Default, FromEdn)]
pub struct MoleculeAstConfigOverrides {
    #[edn(default)]
    pub atom: AtomAstConfigOverrides,
    #[edn(default)]
    pub bond: BondAstConfigOverrides,
    #[edn(default)]
    pub aromatic_system: AromaticSystemAstConfigOverrides,
    #[edn(default)]
    pub multicenter_bond: MulticenterBondAstConfigOverrides,
    #[edn(default)]
    pub dative_bond: DativeBondAstConfigOverrides,
    #[edn(default)]
    pub noncovalent_bond: NoncovalentBondAstConfigOverrides,
}

impl MoleculeAstConfig {
    pub fn with_overrides(self, ov: MoleculeAstConfigOverrides) -> Self {
        Self {
            atom: self.atom.with_overrides(ov.atom),
            bond: self.bond.with_overrides(ov.bond),
            aromatic_system: self.aromatic_system.with_overrides(ov.aromatic_system),
            multicenter_bond: self.multicenter_bond.with_overrides(ov.multicenter_bond),
            dative_bond: self.dative_bond.with_overrides(ov.dative_bond),
            noncovalent_bond: self.noncovalent_bond.with_overrides(ov.noncovalent_bond),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_edn::{read_string, FromEdn, ToEdn};

    use super::*;

    #[rstest]
    #[case::zeroed(MoleculeAstConfig::zeroed())]
    #[case::open(MoleculeAstConfig::open())]
    fn test_molecule_ast_config_roundtrip(#[case] cfg: MoleculeAstConfig) {
        let edn = cfg.to_edn();
        let back = MoleculeAstConfig::from_edn(&edn).unwrap();
        assert_eq!(cfg.atom.charge_mode, back.atom.charge_mode);
        assert_eq!(cfg.atom.implicit_h_mode, back.atom.implicit_h_mode);
        assert_eq!(cfg.bond.charge_mode, back.bond.charge_mode);
        assert_eq!(cfg.aromatic_system.electrons_mode, back.aromatic_system.electrons_mode);
        assert_eq!(cfg.multicenter_bond.electrons_mode, back.multicenter_bond.electrons_mode);
    }

    #[rstest]
    #[case::zeroed_atom(
        "{:isotope-mode :natural :charge-mode :zero :implicit-h-mode :zero \
         :lone-pairs-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived \
         :valence-mode :zero :donated-pairs-mode :zero :accepted-pairs-mode :zero \
         :multicenter-valence-mode :not-multicenter :aromatic-valence-mode :not-aromatic}",
        NumericMode::Zero, ImplicitHydrogenMode::Zero,
        MulticenterValenceMode::NotMulticenter, AromaticValenceMode::NotAromatic,
    )]
    #[case::open_atom(
        "{:isotope-mode :required :charge-mode :required :implicit-h-mode :required \
         :lone-pairs-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived \
         :valence-mode :required :donated-pairs-mode :required :accepted-pairs-mode :required \
         :multicenter-valence-mode :required :aromatic-valence-mode :required}",
        NumericMode::Required, ImplicitHydrogenMode::Required,
        MulticenterValenceMode::Required, AromaticValenceMode::Required,
    )]
    fn test_atom_ast_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumericMode,
        #[case] expected_h: ImplicitHydrogenMode,
        #[case] expected_multicenter: MulticenterValenceMode,
        #[case] expected_aromatic: AromaticValenceMode,
    ) {
        let tree = read_string(edn).unwrap();
        let cfg = AtomAstConfig::from_edn(&tree).unwrap();
        assert_eq!(cfg.charge_mode, expected_charge);
        assert_eq!(cfg.implicit_h_mode, expected_h);
        assert_eq!(cfg.multicenter_valence_mode, expected_multicenter);
        assert_eq!(cfg.aromatic_valence_mode, expected_aromatic);
    }
}
