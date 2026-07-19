//! Morgan featurizer: a bit-exact replica of RDKit's connectivity Morgan
//! fingerprint (`GetMorganFingerprint`, useFeatures=false), pinned to RDKit
//! 2026.03.x. The invariant, bond typing, iteration, and 32-bit boost hash all
//! follow RDKit; the dedup is the shared circular refinement in `umol-graph-core`.
//!
//! Public entry points reject non-ground molecules before invariant extraction.
//! `deltaMass` follows RDKit: natural isotope → 0; a labelled isotope →
//! `int(exact isotope mass − standard atomic weight)`. That integer truncation
//! is table-sensitive, so isotope bit-exactness depends on our mass tables
//! agreeing with RDKit's; natural atoms are always 0.

use umol_ast::ast::{AsLit, AtomId, BondId, IsotopeMassAst, MoleculeAst, RingSet};
use umol_chem::isotope::Isotope;
use umol_graph_core::CircularRefinementAlgorithm;

use super::feature_set::{CountedFeatureSet, FeatureSet};
use super::featurizer::FingerprintError;
use crate::hash::Morgan;

/// RDKit Morgan fingerprint of `radius` iterations (ECFP_{2·radius} equivalent).
#[derive(Clone, Copy, Debug)]
pub struct MorganFeaturizer {
    pub radius: u32,
}

impl MorganFeaturizer {
    pub fn new(radius: u32) -> Self {
        Self { radius }
    }

    /// Returns the deduplicated set of identifiers.
    pub fn featurize(&self, mol: &MoleculeAst) -> Result<FeatureSet<u64>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        Ok(FeatureSet::from_features(self.identifiers(mol)))
    }

    /// Compute per-identifier occurrences.
    pub fn featurize_counted(
        &self,
        mol: &MoleculeAst,
    ) -> Result<CountedFeatureSet<u64>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        Ok(CountedFeatureSet::from_features(self.identifiers(mol)))
    }

    /// The circular-refinement identifier multiset (one per surviving environment);
    /// dedup yields the binary set, counting yields the counted set.
    fn identifiers(&self, mol: &MoleculeAst) -> Vec<u64> {
        let rings = mol.rings().into_ring_set();
        mol.raw_graph().circular_refine(
            |node| atom_components(mol, &rings, AtomId::from(node)),
            |edge| bond_type(mol, BondId::from(edge)),
            CircularRefinementAlgorithm::Ec {
                radius: self.radius,
                scheme: Morgan,
            },
        )
    }
}

/// RDKit `getConnectivityInvariants` component vector: atomic number, total degree
/// (heavy + H), attached H count, formal charge, deltaMass, and a trailing `1` if
/// the atom is in a ring.
fn atom_components(mol: &MoleculeAst, rings: &RingSet, id: AtomId) -> Vec<u32> {
    let atom = mol.atom(id);
    let element = atom.element().as_lit().expect("ground atom");
    let heavy_degree = atom.heavy_atom_degree().as_lit().expect("ground atom");
    let hydrogens = atom.total_hydrogens().as_lit().expect("ground atom");
    let charge = atom.charge().as_lit().expect("ground atom");
    let delta_mass: i32 = match atom.isotope_mass() {
        IsotopeMassAst::Natural => 0,
        IsotopeMassAst::Lit(mass_number) => {
            let isotope = Isotope::checked_new(element, *mass_number).expect("valid isotope");
            (isotope.mass() - element.mass()) as i32
        }
        _ => unreachable!("ground atom isotope is Natural or Lit"),
    };

    let mut components = vec![
        u32::from(element.atomic_number()),
        (heavy_degree + hydrogens) as u32,
        hydrogens as u32,
        charge as i32 as u32,
        delta_mass as u32,
    ];
    if rings.contains_atom(id) {
        components.push(1);
    }
    components
}

/// RDKit `Bond::BondType` integer: aromatic bonds are `12`, otherwise the order.
fn bond_type(mol: &MoleculeAst, id: BondId) -> u32 {
    let bond = mol.bond(id);
    if bond.is_in_aromatic_system() {
        12
    } else {
        bond.order().as_lit().expect("ground bond") as u32
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::ingest::ingest_smiles;

    // RDKit 2026.03.3 `GetMorganFingerprint` sparse identifiers (sorted) per
    // (SMILES, radius). Bit-exact target for the replica.
    #[rstest]
    #[case::propane_r0("CCC", 0, &[2245384272, 2246728737])]
    #[case::propane_r1("CCC", 1, &[2068133184, 2245384272, 2246728737, 3542456614])]
    #[case::ethanol_r0("CCO", 0, &[864662311, 2245384272, 2246728737])]
    #[case::ethanol_r1("CCO", 1, &[864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386])]
    #[case::acetic_r1("CC(=O)O", 1, &[864662311, 864942730, 1510328189, 1533864325, 2205501948, 2246699815, 2246728737, 3545365497])]
    #[case::isobutane_r1("CC(C)C", 1, &[1511968919, 2245273601, 2246728737, 3537119515])]
    #[case::cyclopropane_r2("C1CC1", 2, &[171200514, 2142032900, 2968968094])]
    #[case::ethylamine_r1("CCN", 1, &[772817685, 847957139, 2245384272, 2246728737, 2592785365, 3542456614])]
    #[case::benzene_r2("c1ccccc1", 2, &[98513984, 2763854213, 3218693969])]
    #[case::pyridine_r2("c1ccncc1", 2, &[98513984, 1207774339, 1343371647, 1821698485, 2041434490, 2763854213, 3118255683, 3218693969, 3776905034])]
    #[case::pyrrole_r2("c1cc[nH]c1", 2, &[98513984, 116898731, 1482649460, 2132511834, 2293755984, 2654043257, 2753863138, 3218693969])]
    #[case::furan_r2("c1ccoc1", 2, &[98513984, 1325841767, 1832581338, 3143719699, 3189457552, 3218693969, 3872712528, 4278515623])]
    #[case::chlorobenzene_r2("Clc1ccccc1", 2, &[98513984, 951226070, 1016841875, 2246340824, 2604440622, 2763854213, 3084241488, 3217380708, 3218693969, 3452535345, 3999906991])]
    #[case::acetate_r1("CC(=O)[O-]", 1, &[864942730, 864942795, 1510323402, 1510328189, 2246699815, 2246728737, 3219326737, 3545365497])]
    #[case::ammonium_r1("C[NH3+]", 1, &[847694221, 2246728737, 2567926842])]
    #[case::naphthalene_r2("c1ccc2ccccc2c1", 2, &[98513984, 951226070, 2126281302, 2360741695, 3217380708, 3218693969, 3976623167, 3999906991])]
    fn test_morgan_featurizer_featurize(
        #[case] smiles: &str,
        #[case] radius: u32,
        #[case] expected: &[u64],
    ) {
        let mol = ingest_smiles(smiles).expect("ingest");
        let fingerprint = MorganFeaturizer::new(radius).featurize(&mol).unwrap();
        assert_eq!(fingerprint.ids(), expected);
    }

    // RDKit 2026.03.3 `GetMorganFingerprint().GetNonzeroElements()` — the same
    // identifiers with their occurrence counts. Validates the counted path.
    #[rstest]
    #[case::propane_r1("CCC", 1, &[(2068133184, 1), (2245384272, 1), (2246728737, 2), (3542456614, 2)])]
    #[case::ethanol_r1("CCO", 1, &[(864662311, 1), (1535166686, 1), (2245384272, 1), (2246728737, 1), (3542456614, 1), (4018048386, 1)])]
    #[case::benzene_r2("c1ccccc1", 2, &[(98513984, 6), (2763854213, 6), (3218693969, 6)])]
    fn test_morgan_featurizer_featurize_counted(
        #[case] smiles: &str,
        #[case] radius: u32,
        #[case] expected: &[(u64, u32)],
    ) {
        let mol = ingest_smiles(smiles).expect("ingest");
        let fingerprint = MorganFeaturizer::new(radius)
            .featurize_counted(&mol)
            .unwrap();
        assert_eq!(fingerprint.entries(), expected);
    }
}
