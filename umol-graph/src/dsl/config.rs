//! DSL configuration: mode enums and config structs for lowering/raising.

use serde::{Deserialize, Serialize};

/// Isotope interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum IsotopeMode {
    Natural,  // absent → Natural
    Required, // absent → Any
}

/// Numeric field interpretation mode.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum NumericMode {
    Zero,     // absent → Lit(0), field optional
    Required, // absent → Any/wildcard, field required for grounding
}

/// Implicit hydrogen interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ImplicitHydrogenMode {
    Zero,     // absent → Lit(0)
    Normal,   // absent → Normal (deferred constraint)
    Required, // absent → Any
}

/// Unpaired electrons interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum UnpairedElectronsMode {
    Zero,     // absent → Lit(0)
    Required, // absent → Any
    Derived,  // absent + m present → derive from m (m-1); absent + m absent → Any
}

/// Multiplicity interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MultiplicityMode {
    Derived,  // absent → derive from unpaired electrons (u+1); absent + u absent → Any
    Required, // absent → Any
}

/// Aromatic interpretation mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum AromaticValenceMode {
    NotAromatic, // absent → AromaticExpr::NotAromatic (#a!)
    Aromatic,    // absent → AromaticExpr::Value(Wildcard) (#a*)
    Required,    // absent → Any
}

/// Atom DSL configuration for lowering and raising.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct AtomDslConfig {
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

impl AtomDslConfig {
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

/// Partial atom DSL config: all fields optional, for merging onto a base config.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, rename_all = "kebab-case")]
pub struct AtomDslConfigOverrides {
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

impl AtomDslConfig {
    pub fn with_overrides(mut self, ov: AtomDslConfigOverrides) -> Self {
        if let Some(v) = ov.isotope_mode { self.isotope_mode = v; }
        if let Some(v) = ov.charge_mode { self.charge_mode = v; }
        if let Some(v) = ov.implicit_h_mode { self.implicit_h_mode = v; }
        if let Some(v) = ov.lone_pairs_mode { self.lone_pairs_mode = v; }
        if let Some(v) = ov.unpaired_electrons_mode { self.unpaired_electrons_mode = v; }
        if let Some(v) = ov.multiplicity_mode { self.multiplicity_mode = v; }
        if let Some(v) = ov.valence_mode { self.valence_mode = v; }
        if let Some(v) = ov.donated_pairs_mode { self.donated_pairs_mode = v; }
        if let Some(v) = ov.accepted_pairs_mode { self.accepted_pairs_mode = v; }
        if let Some(v) = ov.aromatic_valence_mode { self.aromatic_valence_mode = v; }
        if let Some(v) = ov.multicenter_valence_mode { self.multicenter_valence_mode = v; }
        self
    }
}

/// Bond DSL configuration for lowering and raising.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct BondDslConfig {
    pub charge_mode: NumericMode,
    pub unpaired_electrons_mode: UnpairedElectronsMode,
    pub multiplicity_mode: MultiplicityMode,
}

impl BondDslConfig {
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

/// Partial bond DSL config: all fields optional, for merging onto a base config.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, rename_all = "kebab-case")]
pub struct BondDslConfigOverrides {
    pub charge_mode: Option<NumericMode>,
    pub unpaired_electrons_mode: Option<UnpairedElectronsMode>,
    pub multiplicity_mode: Option<MultiplicityMode>,
}

impl BondDslConfig {
    pub fn with_overrides(mut self, ov: BondDslConfigOverrides) -> Self {
        if let Some(v) = ov.charge_mode { self.charge_mode = v; }
        if let Some(v) = ov.unpaired_electrons_mode { self.unpaired_electrons_mode = v; }
        if let Some(v) = ov.multiplicity_mode { self.multiplicity_mode = v; }
        self
    }
}

/// Partial molecule DSL config: all fields optional, for merging onto a base config.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, rename_all = "kebab-case")]
pub struct MoleculeDslConfigOverrides {
    pub atom: AtomDslConfigOverrides,
    pub bond: BondDslConfigOverrides,
}

impl MoleculeDslConfig {
    pub fn with_overrides(self, ov: MoleculeDslConfigOverrides) -> Self {
        Self {
            atom: self.atom.with_overrides(ov.atom),
            bond: self.bond.with_overrides(ov.bond),
        }
    }
}

/// Molecule DSL configuration (combines atom + bond configs).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct MoleculeDslConfig {
    pub atom: AtomDslConfig,
    pub bond: BondDslConfig,
}

impl MoleculeDslConfig {
    pub fn zeroed() -> Self {
        Self {
            atom: AtomDslConfig::zeroed(),
            bond: BondDslConfig::zeroed(),
        }
    }

    pub fn open() -> Self {
        Self {
            atom: AtomDslConfig::open(),
            bond: BondDslConfig::open(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use umol_edn::{from_str as edn_from_str, to_string as edn_to_string};

    #[rstest]
    #[case::zeroed(MoleculeDslConfig::zeroed(),
        concat!("{:atom {:isotope-mode :natural :charge-mode :zero :implicit-h-mode :zero :lone-pairs-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived :valence-mode :zero :donated-pairs-mode :zero ",
                ":accepted-pairs-mode :zero :aromatic-valence-mode :not-aromatic :multicenter-valence-mode :zero} :bond {:charge-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived}}"))]
    #[case::open(MoleculeDslConfig::open(),
        concat!("{:atom {:isotope-mode :required :charge-mode :required :implicit-h-mode :required :lone-pairs-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived :valence-mode :required :donated-pairs-mode :required ",
                ":accepted-pairs-mode :required :aromatic-valence-mode :required :multicenter-valence-mode :required} :bond {:charge-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived}}"))]
    fn test_molecule_dsl_config_to_edn(#[case] cfg: MoleculeDslConfig, #[case] expected: &str) {
        assert_eq!(edn_to_string(&cfg).unwrap(), expected);
    }

    #[rstest]
    #[case::zeroed(concat!("{:atom {:isotope-mode :natural :charge-mode :zero :implicit-h-mode :zero :lone-pairs-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived :valence-mode :zero :donated-pairs-mode :zero ",
                           ":accepted-pairs-mode :zero :aromatic-valence-mode :not-aromatic :multicenter-valence-mode :zero} :bond {:charge-mode :zero :unpaired-electrons-mode :zero :multiplicity-mode :derived}}"),
        NumericMode::Zero, ImplicitHydrogenMode::Zero, AromaticValenceMode::NotAromatic)]
    #[case::open(concat!("{:atom {:isotope-mode :required :charge-mode :required :implicit-h-mode :required :lone-pairs-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived :valence-mode :required ",
                         ":donated-pairs-mode :required :accepted-pairs-mode :required :aromatic-valence-mode :required :multicenter-valence-mode :required} :bond {:charge-mode :required :unpaired-electrons-mode :derived :multiplicity-mode :derived}}"),
        NumericMode::Required, ImplicitHydrogenMode::Required, AromaticValenceMode::Required)]
    fn test_molecule_dsl_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumericMode,
        #[case] expected_h: ImplicitHydrogenMode,
        #[case] expected_aromatic: AromaticValenceMode,
    ) {
        let cfg: MoleculeDslConfig = edn_from_str(edn).unwrap();
        assert_eq!(cfg.atom.charge_mode, expected_charge);
        assert_eq!(cfg.atom.implicit_h_mode, expected_h);
        assert_eq!(cfg.atom.aromatic_valence_mode, expected_aromatic);
    }
}
