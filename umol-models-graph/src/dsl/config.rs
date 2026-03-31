//! DSL configuration: mode enums and config structs for lowering/raising.

/// Numeric field interpretation mode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NumericMode {
    Zero,     // absent → Lit(0), field optional
    Required, // absent → Any/wildcard, field required for grounding
}

/// Unpaired electrons interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum UnpairedElectronsMode {
    Zero,     // absent → Lit(0)
    Required, // absent → Any
    Derived,  // absent + m present → derive from m (m-1); absent + m absent → Any
}

/// Multiplicity interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MultiplicityMode {
    Derived,  // absent → derive from unpaired electrons (u+1); absent + u absent → Any
    Required, // absent → Any
}

/// Isotope interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsotopeMode {
    Natural,  // absent → Natural
    Required, // absent → Any
}

/// Implicit hydrogen interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ImplicitHydrogenMode {
    Zero,     // absent → Lit(0)
    Normal,   // absent → Normal (deferred constraint)
    Required, // absent → Any
}

/// Aromatic interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticValenceMode {
    NotAromatic, // absent → AromaticExpr::NotAromatic (#a!)
    Aromatic,    // absent → AromaticExpr::Value(Wildcard) (#a*)
    Required,    // absent → Any
}

/// Atom DSL configuration for lowering and raising.
#[derive(Clone, Debug)]
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

/// Bond DSL configuration for lowering and raising.
#[derive(Clone, Debug)]
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

/// Molecule DSL configuration (combines atom + bond configs).
#[derive(Debug, Clone)]
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
