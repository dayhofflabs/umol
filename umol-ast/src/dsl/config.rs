//! AST lowering/raising configuration: mode enums and config structs.

use umol_edn::{FromEdn, ToEdn};

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

/// Complete set of lowering/raising defaults for atoms: how each
/// `AtomAst` field is treated when converting between DSL and AST.
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

    pub fn verbatim() -> Self {
        Self {
            isotope: IsotopeDefault::Required,
            charge: NumericDefault::Required,
            implicit_hydrogens: ImplicitHydrogensDefault::Required,
            lone_pairs: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Required,
            multiplicity: MultiplicityDefault::Required,
            valence: NumericDefault::Required,
            donated_pairs: NumericDefault::Required,
            accepted_pairs: NumericDefault::Required,
            multicenter_valence: MulticenterValenceDefault::Required,
            aromatic_valence: AromaticValenceDefault::Required,
        }
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

impl AtomDefaults {
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

/// Lowering/raising defaults for covalent bonds.
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct BondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
}

impl BondDefaults {
    pub fn zeroed() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
        }
    }

    pub fn verbatim() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
        }
    }
}

/// Sparse overrides on `BondDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct BondOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
}

impl BondDefaults {
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

/// Lowering/raising defaults for aromatic systems.
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct AromaticSystemDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub electrons: NumericDefault,
}

impl AromaticSystemDefaults {
    pub fn zeroed() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            electrons: NumericDefault::Zero,
        }
    }

    pub fn verbatim() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
            electrons: NumericDefault::Required,
        }
    }
}

/// Sparse overrides on `AromaticSystemDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct AromaticSystemOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
    pub electrons: Option<NumericDefault>,
}

impl AromaticSystemDefaults {
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
        if let Some(v) = ov.electrons {
            self.electrons = v;
        }
        self
    }
}

/// Lowering/raising defaults for multicenter bonds.
#[derive(Clone, Debug, FromEdn, ToEdn)]
pub struct MulticenterBondDefaults {
    pub charge: NumericDefault,
    pub unpaired_electrons: UnpairedElectronsDefault,
    pub multiplicity: MultiplicityDefault,
    pub electrons: NumericDefault,
}

impl MulticenterBondDefaults {
    pub fn zeroed() -> Self {
        Self {
            charge: NumericDefault::Zero,
            unpaired_electrons: UnpairedElectronsDefault::Zero,
            multiplicity: MultiplicityDefault::Derived,
            electrons: NumericDefault::Zero,
        }
    }

    pub fn verbatim() -> Self {
        Self {
            charge: NumericDefault::Required,
            unpaired_electrons: UnpairedElectronsDefault::Derived,
            multiplicity: MultiplicityDefault::Derived,
            electrons: NumericDefault::Required,
        }
    }
}

/// Sparse overrides on `MulticenterBondDefaults`.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct MulticenterBondOverrides {
    pub charge: Option<NumericDefault>,
    pub unpaired_electrons: Option<UnpairedElectronsDefault>,
    pub multiplicity: Option<MultiplicityDefault>,
    pub electrons: Option<NumericDefault>,
}

impl MulticenterBondDefaults {
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
        if let Some(v) = ov.electrons {
            self.electrons = v;
        }
        self
    }
}

/// Lowering/raising defaults for dative bonds. Currently empty
/// (no defaultable fields); exists for API uniformity.
#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
pub struct DativeBondDefaults {}

impl DativeBondDefaults {
    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn verbatim() -> Self {
        Self {}
    }
}

/// Lowering/raising defaults for noncovalent bonds. Currently empty
/// (no defaultable fields); exists for API uniformity.
#[derive(Clone, Debug, Default, FromEdn, ToEdn)]
pub struct NoncovalentBondDefaults {}

impl NoncovalentBondDefaults {
    pub fn zeroed() -> Self {
        Self {}
    }

    pub fn verbatim() -> Self {
        Self {}
    }
}

/// Aggregated lowering/raising defaults for an entire molecule. One
/// per-entity-kind defaults bundle; consumed by the molecule-level
/// `FromAst` / `IntoAst` implementations.
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

    pub fn verbatim() -> Self {
        Self {
            atom: AtomDefaults::verbatim(),
            bond: BondDefaults::verbatim(),
            aromatic_system: AromaticSystemDefaults::verbatim(),
            multicenter_bond: MulticenterBondDefaults::verbatim(),
            dative_bond: DativeBondDefaults::verbatim(),
            noncovalent_bond: NoncovalentBondDefaults::verbatim(),
        }
    }
}

/// Sparse overrides on `DativeBondDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct DativeBondOverrides {}

impl DativeBondDefaults {
    pub fn with_overrides(self, _ov: DativeBondOverrides) -> Self {
        self
    }
}

/// Sparse overrides on `NoncovalentBondDefaults`. Currently empty.
#[derive(Clone, Debug, Default, FromEdn)]
pub struct NoncovalentBondOverrides {}

impl NoncovalentBondDefaults {
    pub fn with_overrides(self, _ov: NoncovalentBondOverrides) -> Self {
        self
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

impl MoleculeDefaults {
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

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_edn::{read_string, FromEdn, ToEdn};

    use super::*;

    #[rstest]
    #[case::zeroed(MoleculeDefaults::zeroed())]
    #[case::verbatim(MoleculeDefaults::verbatim())]
    fn test_molecule_ast_config_roundtrip(#[case] cfg: MoleculeDefaults) {
        let edn = cfg.to_edn();
        let back = MoleculeDefaults::from_edn(&edn).unwrap();
        assert_eq!(cfg.atom.charge, back.atom.charge);
        assert_eq!(cfg.atom.implicit_hydrogens, back.atom.implicit_hydrogens);
        assert_eq!(cfg.bond.charge, back.bond.charge);
        assert_eq!(
            cfg.aromatic_system.electrons,
            back.aromatic_system.electrons
        );
        assert_eq!(
            cfg.multicenter_bond.electrons,
            back.multicenter_bond.electrons
        );
    }

    #[rstest]
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
    #[case::open_atom(
        "{:isotope :required :charge :required :implicit-hydrogens :required \
         :lone-pairs :required :unpaired-electrons :derived :multiplicity :derived \
         :valence :required :donated-pairs :required :accepted-pairs :required \
         :multicenter-valence :required :aromatic-valence :required}",
        NumericDefault::Required,
        ImplicitHydrogensDefault::Required,
        MulticenterValenceDefault::Required,
        AromaticValenceDefault::Required
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

    // -- with_overrides --------

    #[rstest]
    fn test_atom_defaults_with_overrides_all_fields() {
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
    fn test_atom_defaults_with_overrides_partial_preserves_unset() {
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
    fn test_bond_defaults_with_overrides_all_fields() {
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
    fn test_aromatic_system_defaults_with_overrides_all_fields() {
        let cfg = AromaticSystemDefaults::zeroed().with_overrides(AromaticSystemOverrides {
            charge: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
            electrons: Some(NumericDefault::Required),
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.electrons, NumericDefault::Required);
    }

    #[rstest]
    fn test_multicenter_bond_defaults_with_overrides_all_fields() {
        let cfg = MulticenterBondDefaults::zeroed().with_overrides(MulticenterBondOverrides {
            charge: Some(NumericDefault::Required),
            unpaired_electrons: Some(UnpairedElectronsDefault::Derived),
            multiplicity: Some(MultiplicityDefault::Required),
            electrons: Some(NumericDefault::Required),
        });
        assert_eq!(cfg.charge, NumericDefault::Required);
        assert_eq!(cfg.unpaired_electrons, UnpairedElectronsDefault::Derived);
        assert_eq!(cfg.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.electrons, NumericDefault::Required);
    }

    #[rstest]
    fn test_dative_bond_defaults_with_overrides_is_noop() {
        // `DativeBondDefaults` has no fields; `with_overrides` is an identity.
        // Call-site coverage only; no field comparison.
        let _ = DativeBondDefaults::zeroed().with_overrides(DativeBondOverrides::default());
    }

    #[rstest]
    fn test_noncovalent_bond_defaults_with_overrides_is_noop() {
        let _ =
            NoncovalentBondDefaults::zeroed().with_overrides(NoncovalentBondOverrides::default());
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
                electrons: Some(NumericDefault::Required),
                ..AromaticSystemOverrides::default()
            },
            multicenter_bond: MulticenterBondOverrides {
                charge: Some(NumericDefault::Required),
                ..MulticenterBondOverrides::default()
            },
            dative_bond: DativeBondOverrides::default(),
            noncovalent_bond: NoncovalentBondOverrides::default(),
        });
        assert_eq!(cfg.atom.charge, NumericDefault::Required);
        assert_eq!(cfg.bond.multiplicity, MultiplicityDefault::Required);
        assert_eq!(cfg.aromatic_system.electrons, NumericDefault::Required);
        assert_eq!(cfg.multicenter_bond.charge, NumericDefault::Required);
        // Untouched per-entity fields retain zeroed() values.
        assert_eq!(cfg.atom.isotope, IsotopeDefault::Natural);
        assert_eq!(cfg.bond.charge, NumericDefault::Zero);
    }
}
