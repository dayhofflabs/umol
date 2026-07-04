//! A molecule's parse-time **namespace**: per entity kind, a running count, an id-keyword lookup, and
//! a participant lookup, plus the atom-alias table. Grown while parsing a molecule (or applying
//! reaction deltas); the roundtrip subset projects out as [`MoleculeMetadata`]. Reshapes the former
//! `EntityCounts` so index bounds are checked as entities are parsed rather than only at the end, and
//! structural refs (a non-atom entity named by its constituent atoms/bonds) resolve against it.
//!
//! Cost splits by kind: atoms carry no participant lookup (the base kind), bonds an O(1)
//! `(min,max) → id` endpoint map (a bond is a graph edge), the five overlays a participant index
//! over their small collections. §4.1 uniqueness (no two same-constituent entries) makes every
//! participant lookup a ≤1 hit.

// Wired into molecule parsing (S2b) and the reaction delta loop (S2d); until then the participant/name
// query surface is exercised only by the unit tests below.
#![allow(dead_code)]

use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;

use bimap::BiBTreeMap;
use indexmap::IndexMap;

use super::atom::AtomDsl;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;

/// Count + id-keyword lookup for one entity kind. Atoms — the base kind — use it directly; the
/// other kinds wrap it in a [`KeyedRegistry`] that adds the participant lookup.
#[derive(Debug)]
struct NamedRegistry<Id> {
    count: usize,
    by_name: IndexMap<String, Id>,
}

impl<Id> Default for NamedRegistry<Id> {
    fn default() -> Self {
        Self {
            count: 0,
            by_name: IndexMap::new(),
        }
    }
}

impl<Id: Copy + From<usize>> NamedRegistry<Id> {
    /// Reserve the next id (growing the count) and, if named, record its `:id` keyword.
    fn register(&mut self, name: Option<String>) -> Id {
        let id = Id::from(self.count);
        self.count += 1;
        if let Some(name) = name {
            self.by_name.insert(name, id);
        }
        id
    }

    fn by_name(&self, name: &str) -> Option<Id> {
        self.by_name.get(name).copied()
    }

    /// The named entities of this kind as `(id, name)` pairs — the inverse of `by_name`, the
    /// projection `MoleculeMetadata` needs for rendering.
    fn names(&self) -> impl Iterator<Item = (Id, &str)> {
        self.by_name.iter().map(|(name, &id)| (id, name.as_str()))
    }

    fn count(&self) -> usize {
        self.count
    }

    /// A registry whose count starts at `count` (so `register` hands out ids from there on) with an
    /// empty name map — the shape a delta namespace takes over its lhs's id space.
    fn with_count(count: usize) -> Self {
        Self {
            count,
            by_name: IndexMap::new(),
        }
    }
}

/// Count + id-keyword + participant lookup for one non-atom entity kind. `Key` is the entity's
/// canonical participant key (a normalized endpoint pair, atom set, donor-set + acceptor, or stereo
/// site); §4.1 uniqueness makes it injective, so a hit is unique.
#[derive(Debug)]
struct KeyedRegistry<Id, Key> {
    named: NamedRegistry<Id>,
    by_participants: HashMap<Key, Id>,
}

impl<Id, Key> Default for KeyedRegistry<Id, Key> {
    fn default() -> Self {
        Self {
            named: NamedRegistry::default(),
            by_participants: HashMap::new(),
        }
    }
}

impl<Id: Copy + From<usize>, Key: Eq + Hash> KeyedRegistry<Id, Key> {
    /// Reserve the next id and record its `:id` keyword and canonical participant key.
    fn register(&mut self, name: Option<String>, key: Key) -> Id {
        let id = self.named.register(name);
        self.by_participants.insert(key, id);
        id
    }

    fn by_participants(&self, key: &Key) -> Option<Id> {
        self.by_participants.get(key).copied()
    }

    fn by_name(&self, name: &str) -> Option<Id> {
        self.named.by_name(name)
    }

    fn names(&self) -> impl Iterator<Item = (Id, &str)> {
        self.named.names()
    }

    fn count(&self) -> usize {
        self.named.count()
    }

    fn with_count(count: usize) -> Self {
        Self {
            named: NamedRegistry::with_count(count),
            by_participants: HashMap::new(),
        }
    }
}

/// The unordered endpoint pair of a bond / noncovalent bond, in ascending order — the canonical key.
fn atom_pair_key(a: AtomId, b: AtomId) -> [AtomId; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

/// The canonical key of a stereo element's ligand frame: the ligand **multiset** (sorted, so frame
/// order doesn't matter but repeats — virtual ligands — do), matching `connecting_id`'s semantics.
fn ligand_multiset(ligands: &[StereoLigand]) -> Vec<StereoLigand> {
    let mut ligands = ligands.to_vec();
    ligands.sort_unstable();
    ligands
}

/// The eight per-kind registries built while parsing a molecule or applying reaction deltas.
#[derive(Debug, Default)]
pub(crate) struct MoleculeNamespace {
    atoms: NamedRegistry<AtomId>,
    bonds: KeyedRegistry<BondId, [AtomId; 2]>,
    dative_bonds: KeyedRegistry<DativeBondId, (BTreeSet<AtomId>, AtomId)>,
    aromatic_systems: KeyedRegistry<AromaticSystemId, BTreeSet<AtomId>>,
    multicenter_bonds: KeyedRegistry<MulticenterBondId, BTreeSet<AtomId>>,
    noncovalent_bonds: KeyedRegistry<NoncovalentBondId, [AtomId; 2]>,
    stereo_atoms: KeyedRegistry<StereoAtomId, (AtomId, Vec<StereoLigand>)>,
    stereo_bonds: KeyedRegistry<StereoBondId, (BondId, Vec<StereoLigand>)>,
    /// The bijective atom-alias table (name ↔ atom-spec template) — part of the atom name namespace
    /// (an `:id` may not collide with an alias name), so the namespace owns it.
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

impl MoleculeNamespace {
    /// A namespace continuing another's id space: each kind's count starts at `other`'s count (so
    /// `register_*` hands out ids following it), with empty name / participant / alias maps — it holds
    /// only the entities registered into it. Used for a reaction's delta namespace over its lhs.
    pub(crate) fn continuation(other: &MoleculeNamespace) -> Self {
        Self {
            atoms: NamedRegistry::with_count(other.atoms.count()),
            bonds: KeyedRegistry::with_count(other.bonds.count()),
            dative_bonds: KeyedRegistry::with_count(other.dative_bonds.count()),
            aromatic_systems: KeyedRegistry::with_count(other.aromatic_systems.count()),
            multicenter_bonds: KeyedRegistry::with_count(other.multicenter_bonds.count()),
            noncovalent_bonds: KeyedRegistry::with_count(other.noncovalent_bonds.count()),
            stereo_atoms: KeyedRegistry::with_count(other.stereo_atoms.count()),
            stereo_bonds: KeyedRegistry::with_count(other.stereo_bonds.count()),
            atom_aliases: BiBTreeMap::new(),
        }
    }

    pub(crate) fn register_atom(&mut self, name: Option<String>) -> AtomId {
        self.atoms.register(name)
    }

    pub(crate) fn register_bond(&mut self, name: Option<String>, a: AtomId, b: AtomId) -> BondId {
        self.bonds.register(name, atom_pair_key(a, b))
    }

    pub(crate) fn register_dative_bond(
        &mut self,
        name: Option<String>,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> DativeBondId {
        self.dative_bonds
            .register(name, (donors.iter().copied().collect(), acceptor))
    }

    pub(crate) fn register_aromatic_system(
        &mut self,
        name: Option<String>,
        atoms: &[AtomId],
    ) -> AromaticSystemId {
        self.aromatic_systems
            .register(name, atoms.iter().copied().collect())
    }

    pub(crate) fn register_multicenter_bond(
        &mut self,
        name: Option<String>,
        atoms: &[AtomId],
    ) -> MulticenterBondId {
        self.multicenter_bonds
            .register(name, atoms.iter().copied().collect())
    }

    pub(crate) fn register_noncovalent_bond(
        &mut self,
        name: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> NoncovalentBondId {
        self.noncovalent_bonds.register(name, atom_pair_key(a, b))
    }

    pub(crate) fn register_stereo_atom(
        &mut self,
        name: Option<String>,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> StereoAtomId {
        self.stereo_atoms
            .register(name, (site, ligand_multiset(ligands)))
    }

    pub(crate) fn register_stereo_bond(
        &mut self,
        name: Option<String>,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> StereoBondId {
        self.stereo_bonds
            .register(name, (site, ligand_multiset(ligands)))
    }

    pub(crate) fn atom_count(&self) -> usize {
        self.atoms.count()
    }

    pub(crate) fn bond_count(&self) -> usize {
        self.bonds.count()
    }

    pub(crate) fn dative_bond_count(&self) -> usize {
        self.dative_bonds.count()
    }

    pub(crate) fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.count()
    }

    pub(crate) fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.count()
    }

    pub(crate) fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.count()
    }

    pub(crate) fn stereo_atom_count(&self) -> usize {
        self.stereo_atoms.count()
    }

    pub(crate) fn stereo_bond_count(&self) -> usize {
        self.stereo_bonds.count()
    }

    pub(crate) fn atom_by_name(&self, name: &str) -> Option<AtomId> {
        self.atoms.by_name(name)
    }

    pub(crate) fn bond_by_name(&self, name: &str) -> Option<BondId> {
        self.bonds.by_name(name)
    }

    pub(crate) fn dative_bond_by_name(&self, name: &str) -> Option<DativeBondId> {
        self.dative_bonds.by_name(name)
    }

    pub(crate) fn aromatic_system_by_name(&self, name: &str) -> Option<AromaticSystemId> {
        self.aromatic_systems.by_name(name)
    }

    pub(crate) fn multicenter_bond_by_name(&self, name: &str) -> Option<MulticenterBondId> {
        self.multicenter_bonds.by_name(name)
    }

    pub(crate) fn noncovalent_bond_by_name(&self, name: &str) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds.by_name(name)
    }

    pub(crate) fn stereo_atom_by_name(&self, name: &str) -> Option<StereoAtomId> {
        self.stereo_atoms.by_name(name)
    }

    pub(crate) fn stereo_bond_by_name(&self, name: &str) -> Option<StereoBondId> {
        self.stereo_bonds.by_name(name)
    }

    pub(crate) fn bond_by_participants(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.bonds.by_participants(&atom_pair_key(a, b))
    }

    pub(crate) fn dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId> {
        self.dative_bonds
            .by_participants(&(donors.iter().copied().collect(), acceptor))
    }

    pub(crate) fn aromatic_system_by_participants(
        &self,
        atoms: &[AtomId],
    ) -> Option<AromaticSystemId> {
        self.aromatic_systems
            .by_participants(&atoms.iter().copied().collect())
    }

    pub(crate) fn multicenter_bond_by_participants(
        &self,
        atoms: &[AtomId],
    ) -> Option<MulticenterBondId> {
        self.multicenter_bonds
            .by_participants(&atoms.iter().copied().collect())
    }

    pub(crate) fn noncovalent_bond_by_participants(
        &self,
        a: AtomId,
        b: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds.by_participants(&atom_pair_key(a, b))
    }

    pub(crate) fn stereo_atom_by_participants(
        &self,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Option<StereoAtomId> {
        self.stereo_atoms
            .by_participants(&(site, ligand_multiset(ligands)))
    }

    pub(crate) fn stereo_bond_by_participants(
        &self,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Option<StereoBondId> {
        self.stereo_bonds
            .by_participants(&(site, ligand_multiset(ligands)))
    }

    pub(crate) fn register_atom_alias(&mut self, name: String, dsl: Box<AtomDsl>) {
        self.atom_aliases.insert(name, dsl);
    }

    pub(crate) fn atom_names(&self) -> impl Iterator<Item = (AtomId, &str)> {
        self.atoms.names()
    }

    pub(crate) fn bond_names(&self) -> impl Iterator<Item = (BondId, &str)> {
        self.bonds.names()
    }

    pub(crate) fn dative_bond_names(&self) -> impl Iterator<Item = (DativeBondId, &str)> {
        self.dative_bonds.names()
    }

    pub(crate) fn aromatic_system_names(&self) -> impl Iterator<Item = (AromaticSystemId, &str)> {
        self.aromatic_systems.names()
    }

    pub(crate) fn multicenter_bond_names(&self) -> impl Iterator<Item = (MulticenterBondId, &str)> {
        self.multicenter_bonds.names()
    }

    pub(crate) fn noncovalent_bond_names(&self) -> impl Iterator<Item = (NoncovalentBondId, &str)> {
        self.noncovalent_bonds.names()
    }

    pub(crate) fn stereo_atom_names(&self) -> impl Iterator<Item = (StereoAtomId, &str)> {
        self.stereo_atoms.names()
    }

    pub(crate) fn stereo_bond_names(&self) -> impl Iterator<Item = (StereoBondId, &str)> {
        self.stereo_bonds.names()
    }

    pub(crate) fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(name, dsl)| (name.as_str(), dsl.as_ref()))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::ast::ligand::StereoLigandKind;

    #[rstest]
    fn test_molecule_namespace_continuation() {
        let mut lhs = MoleculeNamespace::default();
        lhs.register_atom(None);
        lhs.register_atom(Some("c1".into()));
        lhs.register_bond(None, AtomId(0), AtomId(1));

        let mut delta = MoleculeNamespace::continuation(&lhs);
        // Counts continue the lhs id space; names/participants start empty.
        assert_eq!(delta.atom_count(), 2);
        assert_eq!(delta.bond_count(), 1);
        assert_eq!(delta.atom_by_name("c1"), None);
        assert_eq!(delta.bond_by_participants(AtomId(0), AtomId(1)), None);
        // `register_*` hands out global ids following the lhs.
        assert_eq!(delta.register_atom(None), AtomId(2));
        assert_eq!(delta.register_bond(None, AtomId(0), AtomId(2)), BondId(1));
        assert_eq!(delta.atom_count(), 3);
        assert_eq!(
            delta.bond_by_participants(AtomId(0), AtomId(2)),
            Some(BondId(1))
        );
    }

    #[rstest]
    fn test_molecule_namespace_register_atom() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(namespace.register_atom(None), AtomId(0));
        assert_eq!(namespace.register_atom(Some("c1".into())), AtomId(1));
        assert_eq!(namespace.register_atom(None), AtomId(2));
        assert_eq!(namespace.atom_count(), 3);
        assert_eq!(namespace.atom_by_name("c1"), Some(AtomId(1)));
        assert_eq!(namespace.atom_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace.register_bond(None, AtomId(2), AtomId(0)),
            BondId(0)
        );
        assert_eq!(
            namespace.register_bond(Some("b1".into()), AtomId(1), AtomId(3)),
            BondId(1)
        );
        assert_eq!(namespace.bond_count(), 2);
        assert_eq!(namespace.bond_by_name("b1"), Some(BondId(1)));
        assert_eq!(namespace.bond_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_dative_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace.register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0)),
            DativeBondId(0)
        );
        assert_eq!(
            namespace.register_dative_bond(Some("d1".into()), &[AtomId(4)], AtomId(3)),
            DativeBondId(1)
        );
        assert_eq!(namespace.dative_bond_count(), 2);
        assert_eq!(namespace.dative_bond_by_name("d1"), Some(DativeBondId(1)));
        assert_eq!(namespace.dative_bond_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_aromatic_system() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace.register_aromatic_system(None, &[AtomId(0), AtomId(1), AtomId(2)]),
            AromaticSystemId(0)
        );
        assert_eq!(
            namespace.register_aromatic_system(Some("a1".into()), &[AtomId(3), AtomId(4)]),
            AromaticSystemId(1)
        );
        assert_eq!(namespace.aromatic_system_count(), 2);
        assert_eq!(
            namespace.aromatic_system_by_name("a1"),
            Some(AromaticSystemId(1))
        );
        assert_eq!(namespace.aromatic_system_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_multicenter_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace
                .register_multicenter_bond(Some("m".into()), &[AtomId(0), AtomId(1), AtomId(2)]),
            MulticenterBondId(0)
        );
        assert_eq!(
            namespace.register_multicenter_bond(None, &[AtomId(3), AtomId(4)]),
            MulticenterBondId(1)
        );
        assert_eq!(namespace.multicenter_bond_count(), 2);
        assert_eq!(
            namespace.multicenter_bond_by_name("m"),
            Some(MulticenterBondId(0))
        );
        assert_eq!(namespace.multicenter_bond_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_noncovalent_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace.register_noncovalent_bond(None, AtomId(3), AtomId(1)),
            NoncovalentBondId(0)
        );
        assert_eq!(
            namespace.register_noncovalent_bond(Some("n1".into()), AtomId(0), AtomId(4)),
            NoncovalentBondId(1)
        );
        assert_eq!(namespace.noncovalent_bond_count(), 2);
        assert_eq!(
            namespace.noncovalent_bond_by_name("n1"),
            Some(NoncovalentBondId(1))
        );
        assert_eq!(namespace.noncovalent_bond_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_stereo_atom() {
        let mut namespace = MoleculeNamespace::default();
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        assert_eq!(
            namespace.register_stereo_atom(None, AtomId(4), &ligands),
            StereoAtomId(0)
        );
        assert_eq!(
            namespace.register_stereo_atom(Some("s1".into()), AtomId(5), &ligands),
            StereoAtomId(1)
        );
        assert_eq!(namespace.stereo_atom_count(), 2);
        assert_eq!(namespace.stereo_atom_by_name("s1"), Some(StereoAtomId(1)));
        assert_eq!(namespace.stereo_atom_by_name("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_stereo_bond() {
        let mut namespace = MoleculeNamespace::default();
        let ligands = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        assert_eq!(
            namespace.register_stereo_bond(None, BondId(2), &ligands),
            StereoBondId(0)
        );
        assert_eq!(
            namespace.register_stereo_bond(Some("sb1".into()), BondId(0), &ligands),
            StereoBondId(1)
        );
        assert_eq!(namespace.stereo_bond_count(), 2);
        assert_eq!(namespace.stereo_bond_by_name("sb1"), Some(StereoBondId(1)));
        assert_eq!(namespace.stereo_bond_by_name("nope"), None);
    }

    #[rstest]
    #[case::forward(AtomId(0), AtomId(2), Some(BondId(0)))]
    #[case::reversed(AtomId(2), AtomId(0), Some(BondId(0)))]
    #[case::absent(AtomId(0), AtomId(4), None)]
    fn test_molecule_namespace_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<BondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_bond(None, AtomId(2), AtomId(0));
        assert_eq!(namespace.bond_by_participants(a, b), expected);
    }

    #[rstest]
    #[case::donors_reordered(&[AtomId(2), AtomId(1)], AtomId(0), Some(DativeBondId(0)))]
    #[case::wrong_acceptor(&[AtomId(1), AtomId(2)], AtomId(3), None)]
    #[case::wrong_donors(&[AtomId(1), AtomId(3)], AtomId(0), None)]
    fn test_molecule_namespace_dative_bond_by_participants(
        #[case] donors: &[AtomId],
        #[case] acceptor: AtomId,
        #[case] expected: Option<DativeBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0));
        assert_eq!(
            namespace.dative_bond_by_participants(donors, acceptor),
            expected
        );
    }

    #[rstest]
    #[case::reordered(&[AtomId(0), AtomId(1), AtomId(2)], Some(AromaticSystemId(0)))]
    #[case::subset(&[AtomId(0), AtomId(1)], None)]
    #[case::superset(&[AtomId(0), AtomId(1), AtomId(2), AtomId(3)], None)]
    fn test_molecule_namespace_aromatic_system_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)]);
        assert_eq!(namespace.aromatic_system_by_participants(atoms), expected);
    }

    #[rstest]
    #[case::reordered(&[AtomId(2), AtomId(1), AtomId(0)], Some(MulticenterBondId(0)))]
    #[case::absent(&[AtomId(0), AtomId(1), AtomId(3)], None)]
    fn test_molecule_namespace_multicenter_bond_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<MulticenterBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_multicenter_bond(None, &[AtomId(0), AtomId(1), AtomId(2)]);
        assert_eq!(namespace.multicenter_bond_by_participants(atoms), expected);
    }

    #[rstest]
    #[case::reversed(AtomId(1), AtomId(3), Some(NoncovalentBondId(0)))]
    #[case::absent(AtomId(1), AtomId(2), None)]
    fn test_molecule_namespace_noncovalent_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<NoncovalentBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_noncovalent_bond(None, AtomId(3), AtomId(1));
        assert_eq!(namespace.noncovalent_bond_by_participants(a, b), expected);
    }

    #[rstest]
    #[case::reordered_ligands(AtomId(4), &[AtomId(2), AtomId(1)], Some(StereoAtomId(0)))]
    #[case::wrong_ligands(AtomId(4), &[AtomId(1)], None)]
    #[case::wrong_site(AtomId(0), &[AtomId(1), AtomId(2)], None)]
    fn test_molecule_namespace_stereo_atom_by_participants(
        #[case] site: AtomId,
        #[case] ligand_atoms: &[AtomId],
        #[case] expected: Option<StereoAtomId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        // All test ligands are `Atom`-kind; only the atom set varies per case.
        let registered = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        namespace.register_stereo_atom(None, AtomId(4), &registered);
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            namespace.stereo_atom_by_participants(site, &query),
            expected
        );
    }

    #[rstest]
    #[case::matching(BondId(2), &[AtomId(3)], Some(StereoBondId(0)))]
    #[case::wrong_site(BondId(0), &[AtomId(3)], None)]
    #[case::empty_ligands(BondId(2), &[], None)]
    fn test_molecule_namespace_stereo_bond_by_participants(
        #[case] site: BondId,
        #[case] ligand_atoms: &[AtomId],
        #[case] expected: Option<StereoBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        let registered = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        namespace.register_stereo_bond(None, BondId(2), &registered);
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            namespace.stereo_bond_by_participants(site, &query),
            expected
        );
    }
}
