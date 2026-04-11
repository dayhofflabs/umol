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

/// Aromatic interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum AromaticValenceMode {
    NotAromatic,
    Aromatic,
    Required,
}

/// Atom AST configuration for lowering and raising.
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
    pub aromatic_valence_mode: AromaticValenceMode,
    pub multicenter_valence_mode: NumericMode,
}

impl AtomAstConfig {
    /// Ground config: absent fields → zero/natural. For `Atom` lowering/raising.
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
            aromatic_valence_mode: AromaticValenceMode::NotAromatic,
            multicenter_valence_mode: NumericMode::Zero,
        }
    }

    /// Pattern config: absent fields → Any. For `AtomPattern` lowering/raising.
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
            aromatic_valence_mode: AromaticValenceMode::Required,
            multicenter_valence_mode: NumericMode::Required,
        }
    }
}

/// Partial atom AST config: all fields optional, for merging onto a base config.
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
    pub aromatic_valence_mode: Option<AromaticValenceMode>,
    pub multicenter_valence_mode: Option<NumericMode>,
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
        if let Some(v) = ov.aromatic_valence_mode {
            self.aromatic_valence_mode = v;
        }
        if let Some(v) = ov.multicenter_valence_mode {
            self.multicenter_valence_mode = v;
        }
        self
    }
}

/// Bond AST configuration for lowering and raising.
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

/// Partial bond AST config: all fields optional, for merging onto a base config.
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

/// Partial molecule AST config: all fields optional, for merging onto a base config.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct MoleculeAstConfigOverrides {
    #[edn(default)]
    pub atom: AtomAstConfigOverrides,
    #[edn(default)]
    pub bond: BondAstConfigOverrides,
}

impl MoleculeAstConfig {
    pub fn with_overrides(self, ov: MoleculeAstConfigOverrides) -> Self {
        Self {
            atom: self.atom.with_overrides(ov.atom),
            bond: self.bond.with_overrides(ov.bond),
        }
    }
}

/// Molecule AST configuration (combines atom + bond configs).
#[derive(Debug, Clone, FromEdn, ToEdn)]
pub struct MoleculeAstConfig {
    pub atom: AtomAstConfig,
    pub bond: BondAstConfig,
}

impl MoleculeAstConfig {
    pub fn zeroed() -> Self {
        Self {
            atom: AtomAstConfig::zeroed(),
            bond: BondAstConfig::zeroed(),
        }
    }

    pub fn open() -> Self {
        Self {
            atom: AtomAstConfig::open(),
            bond: BondAstConfig::open(),
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
        assert_eq!(cfg.atom.aromatic_valence_mode, back.atom.aromatic_valence_mode);
        assert_eq!(cfg.bond.charge_mode, back.bond.charge_mode);
    }

    #[rstest]
    #[case::zeroed(concat!("{:atom {:isotope-mode :natural :charge-mode :zero :implicit-h-mode :zero :lone-pairs-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived :valence-mode :zero :donated-pairs-mode :zero ",
                           ":accepted-pairs-mode :zero :aromatic-valence-mode :not-aromatic :multicenter-valence-mode :zero} :bond {:charge-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived}}"),
        NumericMode::Zero, ImplicitHydrogenMode::Zero, AromaticValenceMode::NotAromatic)]
    #[case::open(concat!("{:atom {:isotope-mode :required :charge-mode :required :implicit-h-mode :required :lone-pairs-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived :valence-mode :required ",
                         ":donated-pairs-mode :required :accepted-pairs-mode :required :aromatic-valence-mode :required :multicenter-valence-mode :required} :bond {:charge-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived}}"),
        NumericMode::Required, ImplicitHydrogenMode::Required, AromaticValenceMode::Required)]
    fn test_molecule_ast_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumericMode,
        #[case] expected_h: ImplicitHydrogenMode,
        #[case] expected_aromatic: AromaticValenceMode,
    ) {
        let tree = read_string(edn).unwrap();
        let cfg = MoleculeAstConfig::from_edn(&tree).unwrap();
        assert_eq!(cfg.atom.charge_mode, expected_charge);
        assert_eq!(cfg.atom.implicit_h_mode, expected_h);
        assert_eq!(cfg.atom.aromatic_valence_mode, expected_aromatic);
    }
}
