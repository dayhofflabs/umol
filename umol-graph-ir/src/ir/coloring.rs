//! Round-0 graph coloring: one `u64` per graph node — an atom or a relation
//! pseudonode — saying what counts as distinguishable. Used in automorphism
//! symmetry computation; the impl is the policy of "same".

use std::hash::{DefaultHasher, Hash, Hasher};

use bitflags::bitflags;

use super::entity::Entity;
use super::molecule::Molecule;

/// Round-0 color of any graph-participating molecule entity.
pub trait MoleculeColoring {
    fn color(&self, mol: &Molecule, entity: Entity) -> u64;
}

bitflags! {
    /// The inherent fields a `ConstitutionColoring` folds into a color, across
    /// every graph-participating entity kind. Derived constitution predicates
    /// (ring, degree, valence) are excluded here, they are folded in by the
    /// automorphism.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct MoleculeColoringFeatures: u32 {
        // atom
        const ELEMENT = 1 << 0;
        const ISOTOPE = 1 << 1;
        const CHARGE = 1 << 2;
        const IMPLICIT_HYDROGENS = 1 << 3;
        const LONE_PAIRS = 1 << 4;
        const UNPAIRED_ELECTRONS = 1 << 5;
        // bond
        const BOND_ORDER = 1 << 6;
        const BOND_CHARGE = 1 << 7;
        const BOND_UNPAIRED_ELECTRONS = 1 << 8;
        // dative bond
        const DATIVE_ORDER = 1 << 9;
        // aromatic system
        const AROMATIC_ELECTRONS = 1 << 10;
        const AROMATIC_CHARGE = 1 << 11;
        const AROMATIC_UNPAIRED_ELECTRONS = 1 << 12;
        // multicenter bond
        const MULTICENTER_ELECTRONS = 1 << 13;
        const MULTICENTER_CHARGE = 1 << 14;
        const MULTICENTER_UNPAIRED_ELECTRONS = 1 << 15;
        // noncovalent bond
        const NONCOVALENT_KIND = 1 << 16;
        // stereo atom / bond
        const STEREO_KIND = 1 << 17;
    }
}

/// Constitution coloring includes only inherent fields on all entities.
#[derive(Clone, Copy, Debug)]
pub struct ConstitutionColoring {
    features: MoleculeColoringFeatures,
}

impl ConstitutionColoring {
    pub fn new(features: MoleculeColoringFeatures) -> Self {
        Self { features }
    }

    pub fn entity_only() -> Self {
        Self::new(MoleculeColoringFeatures::empty())
    }

    /// Every inherent field selected.
    pub fn full() -> Self {
        Self::new(MoleculeColoringFeatures::all())
    }
}

impl MoleculeColoring for ConstitutionColoring {
    fn color(&self, mol: &Molecule, entity: Entity) -> u64 {
        // The leading `EntityKind` tag keeps kinds in disjoint color ranges; then
        // the entity's inherent fields, each gated by `features`.
        let mut hasher = DefaultHasher::new();
        (entity.kind() as u8).hash(&mut hasher);
        match entity {
            Entity::Atom(id) => {
                let atom = mol.atom(id);
                if self.features.contains(MoleculeColoringFeatures::ELEMENT) {
                    atom.element().hash(&mut hasher);
                }
                if self.features.contains(MoleculeColoringFeatures::ISOTOPE) {
                    atom.isotope_mass().hash(&mut hasher);
                }
                if self.features.contains(MoleculeColoringFeatures::CHARGE) {
                    atom.charge().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::IMPLICIT_HYDROGENS)
                {
                    atom.implicit_hydrogens().hash(&mut hasher);
                }
                if self.features.contains(MoleculeColoringFeatures::LONE_PAIRS) {
                    atom.lone_pairs().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::UNPAIRED_ELECTRONS)
                {
                    atom.unpaired_electrons().hash(&mut hasher);
                }
            }
            Entity::Bond(id) => {
                let bond = mol.bond(id);
                if self.features.contains(MoleculeColoringFeatures::BOND_ORDER) {
                    bond.order().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::BOND_CHARGE)
                {
                    bond.charge().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::BOND_UNPAIRED_ELECTRONS)
                {
                    bond.unpaired_electrons().hash(&mut hasher);
                }
            }
            Entity::AromaticSystem(id) => {
                let system = mol
                    .aromatic_systems()
                    .get(id)
                    .expect("aromatic id in range");
                if self
                    .features
                    .contains(MoleculeColoringFeatures::AROMATIC_ELECTRONS)
                {
                    system.electron_count().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::AROMATIC_CHARGE)
                {
                    system.charge().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::AROMATIC_UNPAIRED_ELECTRONS)
                {
                    system.unpaired_electrons().hash(&mut hasher);
                }
            }
            Entity::MulticenterBond(id) => {
                let bond = mol
                    .multicenter_bonds()
                    .get(id)
                    .expect("multicenter id in range");
                if self
                    .features
                    .contains(MoleculeColoringFeatures::MULTICENTER_ELECTRONS)
                {
                    bond.electron_count().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::MULTICENTER_CHARGE)
                {
                    bond.charge().hash(&mut hasher);
                }
                if self
                    .features
                    .contains(MoleculeColoringFeatures::MULTICENTER_UNPAIRED_ELECTRONS)
                {
                    bond.unpaired_electrons().hash(&mut hasher);
                }
            }
            Entity::DativeBond(id) => {
                let bond = mol.dative_bonds().get(id).expect("dative id in range");
                if self
                    .features
                    .contains(MoleculeColoringFeatures::DATIVE_ORDER)
                {
                    bond.order().hash(&mut hasher);
                }
            }
            Entity::NoncovalentBond(id) => {
                let bond = mol
                    .noncovalent_bonds()
                    .get(id)
                    .expect("noncovalent id in range");
                if self
                    .features
                    .contains(MoleculeColoringFeatures::NONCOVALENT_KIND)
                {
                    bond.kind().hash(&mut hasher);
                }
            }
            Entity::StereoAtom(id) => {
                if self
                    .features
                    .contains(MoleculeColoringFeatures::STEREO_KIND)
                {
                    mol.stereo_atoms()
                        .get(id)
                        .expect("stereo atom id in range")
                        .kind()
                        .hash(&mut hasher);
                }
            }
            Entity::StereoBond(id) => {
                if self
                    .features
                    .contains(MoleculeColoringFeatures::STEREO_KIND)
                {
                    mol.stereo_bonds()
                        .get(id)
                        .expect("stereo bond id in range")
                        .kind()
                        .hash(&mut hasher);
                }
            }
        }
        hasher.finish()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::id::{AtomId, BondId, StereoAtomId, StereoBondId};
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::num::NumForm;
    use crate::ir::spin::UnpairedElectronsForm;
    use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoCoset, StereoKind};

    #[fixture]
    fn ethanol_fragment() -> Molecule {
        // C-C-O: the two carbons share an element; oxygen differs.
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct_element( MoleculeColoringFeatures::all(), Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
    #[case::same_element( MoleculeColoringFeatures::all(), Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), true)]
    #[case::kind_tag_disjoint( MoleculeColoringFeatures::empty(), Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
    #[case::element_not_selected( MoleculeColoringFeatures::empty(), Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), true)]
    fn test_constitution_coloring_color(
        ethanol_fragment: Molecule,
        #[case] features: MoleculeColoringFeatures,
        #[case] left: Entity,
        #[case] right: Entity,
        #[case] expected_equal: bool,
    ) {
        let coloring = ConstitutionColoring::new(features);
        let equal =
            coloring.color(&ethanol_fragment, left) == coloring.color(&ethanol_fragment, right);
        assert_eq!(equal, expected_equal);
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        // Two stereo atoms of different kinds (Tetrahedral, SquarePlanar) and a
        // stereo bond, on a C₆ chain — enough to exercise kind + tag distinction.
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
                (
                    AtomId(3),
                    vec![
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(1)),
                ),
            ],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_distinct(Entity::StereoAtom(StereoAtomId(0)), Entity::StereoAtom(StereoAtomId(1)), MoleculeColoringFeatures::STEREO_KIND, false)]
    #[case::kind_blind(Entity::StereoAtom(StereoAtomId(0)), Entity::StereoAtom(StereoAtomId(1)), MoleculeColoringFeatures::empty(), true)]
    #[case::atom_bond_tag_disjoint(Entity::StereoAtom(StereoAtomId(0)), Entity::StereoBond(StereoBondId(0)), MoleculeColoringFeatures::empty(), false)]
    fn test_constitution_coloring_color_stereo(
        stereo_molecule: Molecule,
        #[case] left: Entity,
        #[case] right: Entity,
        #[case] features: MoleculeColoringFeatures,
        #[case] expected_equal: bool,
    ) {
        let coloring = ConstitutionColoring::new(features);
        let equal = coloring.color(&stereo_molecule, left) == coloring.color(&stereo_molecule, right);
        assert_eq!(equal, expected_equal);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined( UnpairedElectronsForm::default(), UnpairedElectronsForm::default(), MoleculeColoringFeatures::UNPAIRED_ELECTRONS, true)]
    #[case::equal_triplet((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), MoleculeColoringFeatures::UNPAIRED_ELECTRONS, true)]
    #[case::unpaired_differs((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), MoleculeColoringFeatures::UNPAIRED_ELECTRONS, false)]
    #[case::multiplicity_differs((2_u8, 3_u8).into(), (2_u8, 1_u8).into(), MoleculeColoringFeatures::UNPAIRED_ELECTRONS, false)]
    #[case::partial_vs_undetermined(
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined },
        UnpairedElectronsForm::default(),
        MoleculeColoringFeatures::UNPAIRED_ELECTRONS,
        false,
    )]
    #[case::unpaired_electrons_not_selected((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), MoleculeColoringFeatures::empty(), true)]
    fn test_constitution_coloring_color_unpaired_electrons(
        #[case] unpaired_electrons_a: UnpairedElectronsForm,
        #[case] unpaired_electrons_b: UnpairedElectronsForm,
        #[case] features: MoleculeColoringFeatures,
        #[case] expected_equal: bool,
    ) {
        // Two same-element atoms differing only in unpaired electrons; the color matches iff
        // the unpaired-electron components are indistinguishable under `features`.
        let mol = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_unpaired_electrons(unpaired_electrons_a),
                AtomForm::from_element(Element::C)
                    .with_unpaired_electrons(unpaired_electrons_b),
            ],
            bonds: vec![],
            ..Default::default()
        });
        let coloring = ConstitutionColoring::new(features);
        let equal = coloring.color(&mol, Entity::Atom(AtomId(0)))
            == coloring.color(&mol, Entity::Atom(AtomId(1)));
        assert_eq!(equal, expected_equal);
    }
}
