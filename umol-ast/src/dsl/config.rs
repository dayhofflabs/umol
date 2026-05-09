//! AST lowering/raising configuration: mode enums and config structs.

use umol_edn::{FromEdn, ToEdn};

/// Aggregated lowering/raising defaults for molecule DSL <-> AST interconversion. per-entity-kind
/// defaults bundle; consumed by the molecule-level `FromAst` / `IntoAst` implementations.
#[derive(Debug, Clone, FromEdn, ToEdn)]
pub struct MoleculeDefaults {
    pub atom: AtomDefaults,
    pub bond: BondDefaults,
    pub aromatic_system: AromaticSystemDefaults,
    pub multicenter_bond: MulticenterBondDefaults,
    pub dative_bond: DativeBondDefaults,
    pub noncovalent_bond: NoncovalentBondDefaults,
}

impl MoleculeDefaults {
    /// No-op defaults.
    pub fn new() -> Self {
        Self {
            atom: AtomDefaults::new(),
            bond: BondDefaults::new(),
            aromatic_system: AromaticSystemDefaults::new(),
            multicenter_bond: MulticenterBondDefaults::new(),
            dative_bond: DativeBondDefaults::new(),
            noncovalent_bond: NoncovalentBondDefaults::new(),
        }
    }

    /// Composes `*Defaults::ground()` for each entity.
    pub fn ground() -> Self {
        Self {
            atom: AtomDefaults::ground(),
            bond: BondDefaults::ground(),
            aromatic_system: AromaticSystemDefaults::ground(),
            multicenter_bond: MulticenterBondDefaults::ground(),
            dative_bond: DativeBondDefaults::ground(),
            noncovalent_bond: NoncovalentBondDefaults::ground(),
        }
    }

    /// Composes `*Defaults::zeroed()` for each entity.
    pub fn zeroed() -> Self {
        Self {
            atom: AtomDefaults::zeroed(),
            bond: BondDefaults::zeroed(),
            aromatic_system: AromaticSystemDefaults::zeroed(),
            multicenter_bond: MulticenterBondDefaults::zeroed(),
            dative_bond: DativeBondDefaults::zeroed(),
            noncovalent_bond: NoncovalentBondDefaults::zeroed(),
        }
    }

    /// Add overrides.
    pub fn with_overrides(self, ov: MoleculeOverrides) -> Self {
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
    pub aromatic_system: AromaticSystemOverrides,
    #[edn(default)]
    pub multicenter_bond: MulticenterBondOverrides,
    #[edn(default)]
    pub dative_bond: DativeBondOverrides,
    #[edn(default)]
    pub noncovalent_bond: NoncovalentBondOverrides,
}

/// Lowering/raising defaults for atoms: describe how `AtomAst`
/// struct fields and constraintsa are treated when converting between DSL and AST.
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct AtomDefaults {
    pub isotope: IsotopeDefault,
    pub charge: NumericDefault,
    pub implicit_hydrogens: ImplicitHydrogensDefault,
    pub lone_pairs: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub valence: NumericDefault,
    pub donated_pairs: NumericDefault,
    pub accepted_pairs: NumericDefault,
    pub multicenter_valence: MulticenterValenceDefault,
    pub aromatic_valence: AromaticValenceDefault,
}

impl AtomDefaults {
    /// No-op defaults
    pub fn new() -> Self {
        Self {
            isotope: IsotopeDefault::Required,
            charge: NumericDefault::Required,
            implicit_hydrogens: ImplicitHydrogensDefault::Required,
            lone_pairs: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumericDefault::Required,
            donated_pairs: NumericDefault::Required,
            accepted_pairs: NumericDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
        }
    }

    /// Grounds struct fields, no constraints
    pub fn ground() -> Self {
        Self {
            isotope: IsotopeDefault::Natural,
            charge: NumericDefault::Zero,
            implicit_hydrogens: ImplicitHydrogensDefault::Zero,
            lone_pairs: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumericDefault::Required,
            donated_pairs: NumericDefault::Required,
            accepted_pairs: NumericDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
        }
    }

    /// Grounds struct fields and sets all constarints to zero values (for registry entries)
    pub fn zeroed() -> Self {
        Self {
            isotope: IsotopeDefault::Natural,
            charge: NumericDefault::Zero,
            implicit_hydrogens: ImplicitHydrogensDefault::Zero,
            lone_pairs: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            valence: NumericDefault::Zero,
            donated_pairs: NumericDefault::Zero,
            accepted_pairs: NumericDefault::Zero,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
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
    pub implicit_hydrogens: Option<ImplicitHydrogensDefault>,
    pub lone_pairs: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
    pub valence: Option<NumericDefault>,
    pub donated_pairs: Option<NumericDefault>,
    pub accepted_pairs: Option<NumericDefault>,
    pub multicenter_valence: Option<MulticenterValenceDefault>,
    pub aromatic_valence: Option<AromaticValenceDefault>,
}

/// Lowering/raising defaults for localized bonds.
/// See `AtomDefaults` for semantics.
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct BondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl BondDefaults {
    /// No-op defaults
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
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
#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
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
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct AromaticSystemDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl AromaticSystemDefaults {
    /// No-op defaults
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
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
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct MulticenterBondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl MulticenterBondDefaults {
    /// No-op defaults
    pub fn new() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
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
#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
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

/// Implicit hydrogen default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum ImplicitHydrogensDefault {
    Zero,
    Normal,
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
    Aromatic,
    Required,
}

/// Multicenter valence default
#[derive(Clone, Copy, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub enum MulticenterValenceDefault {
    NotMulticenter,
    Multicenter,
    Required,
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_edn::{read_string, FromEdn, ToEdn};

    use super::*;

    #[rstest]
    #[case::zeroed(MoleculeDefaults::zeroed())]
    #[case::verbatim(MoleculeDefaults::new())]
    #[case::ground(MoleculeDefaults::ground())]
    fn test_molecule_ast_config_roundtrip(#[case] cfg: MoleculeDefaults) {
        let edn = cfg.to_edn();
        let back = MoleculeDefaults::from_edn(&edn).unwrap();
        assert_eq!(cfg.atom.charge, back.atom.charge);
        assert_eq!(cfg.atom.implicit_hydrogens, back.atom.implicit_hydrogens);
        assert_eq!(cfg.bond.charge, back.bond.charge);
        assert_eq!(cfg.aromatic_system.charge, back.aromatic_system.charge);
        assert_eq!(cfg.multicenter_bond.charge, back.multicenter_bond.charge);
    }

    #[rstest]
    fn test_molecule_defaults_with_overrides_routes_to_per_entity() {
        let cfg = MoleculeDefaults::zeroed().with_overrides(MoleculeOverrides {
            atom: AtomOverrides {
                charge: Some(NumericDefault::Required),
                ..AtomOverrides::default()
            },
            bond: BondOverrides {
                multiplicity: Some(MultiplicityDefault::Required),
                ..BondOverrides::default()
            },
            aromatic_system: AromaticSystemOverrides {
                charge: Some(NumericDefault::Required),
                ..AromaticSystemOverrides::default()
            },
            multicenter_bond: MulticenterBondOverrides {
                charge: Some(NumericDefault::Required),
                ..MulticenterBondOverrides::default()
            },
            dative_bond: DativeBondOverrides::default(),
            noncovalent_bond: NoncovalentBondOverrides::default(),
        });
        assert_eq!(cfg.atom.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.atom.charge, NumericDefault::Required);
        assert_eq!(cfg.bond.charge, NumericDefault::Zero);
        assert_eq!(cfg.bond.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.aromatic_system.charge, NumericDefault::Required);
        assert_eq!(cfg.multicenter_bond.charge, NumericDefault::Required);
    }

    #[rstest]
    #[case::required_atom(
        "{:isotope :required :charge :required :implicit-hydrogens :required \
         :lone-pairs :required :unpaired-electrons :derived :multiplicity :derived \
         :valence :required :donated-pairs :required :accepted-pairs :required \
         :multicenter-valence :required :aromatic-valence :required}",
        NumericDefault::Required,
        ImplicitHydrogensDefault::Required,
        MulticenterValenceDefault::Required,
        AromaticValenceDefault::Required
    )]
    #[case::zeroed_atom(
        "{:isotope :natural :charge :zero :implicit-hydrogens :zero \
         :lone-pairs :zero :unpaired-electrons :zero :multiplicity :derived \
         :valence :zero :donated-pairs :zero :accepted-pairs :zero \
         :multicenter-valence :not-multicenter :aromatic-valence :not-aromatic}",
        NumericDefault::Zero,
        ImplicitHydrogensDefault::Zero,
        MulticenterValenceDefault::NotMulticenter,
        AromaticValenceDefault::NotAromatic
    )]
    fn test_atom_ast_config_from_edn(
        #[case] edn: &str,
        #[case] expected_charge: NumericDefault,
        #[case] expected_h: ImplicitHydrogensDefault,
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
            implicit_hydrogens: Some(ImplicitHydrogensDefault::Normal),
            lone_pairs: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
            valence: Some(NumericDefault::Required),
            donated_pairs: Some(NumericDefault::Required),
            accepted_pairs: Some(NumericDefault::Required),
            multicenter_valence: Some(MulticenterValenceDefault::Required),
            aromatic_valence: Some(AromaticValenceDefault::Aromatic),
        });
        assert_eq!(cfg.isotope, IsotopeDefault::Required);
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.implicit_hydrogens, ImplicitHydrogensDefault::Normal);
        assert_eq!(cfg.lone_pairs, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.valence, NumericDefault::Required);
        assert_eq!(cfg.donated_pairs, NumericDefault::Required);
        assert_eq!(cfg.accepted_pairs, NumericDefault::Required);
        assert_eq!(cfg.multicenter_valence, MulticenterValenceDefault::Required);
        assert_eq!(cfg.aromatic_valence, AromaticValenceDefault::Aromatic);
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
        assert_eq!(cfg.implicit_hydrogens, ImplicitHydrogensDefault::Zero);
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
