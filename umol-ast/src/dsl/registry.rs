//! Growing per-entity registry built while parsing a molecule (or applying reaction deltas): a
//! running count, an id-keyword lookup, and a participant lookup for each entity kind. Reshapes the
//! former `EntityCounts` so index bounds are checked as entities are parsed rather than only at the
//! end, and structural refs (a non-atom entity named by its constituent atoms/bonds) resolve
//! against it.
//!
//! Cost splits by kind: atoms carry no participant lookup (the base kind), bonds an O(1)
//! `(min,max) → id` endpoint map (a bond is a graph edge), the five overlays a participant index
//! over their small collections. §4.1 uniqueness (no two same-constituent entries) makes every
//! participant lookup a ≤1 hit.

// Wired into molecule parsing (S2b) and the reaction delta loop (S2c); until then the registry API
// is exercised only by the unit tests below.
#![allow(dead_code)]

use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;

use indexmap::IndexMap;

use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

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

    fn count(&self) -> usize {
        self.count
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

    fn count(&self) -> usize {
        self.named.count()
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

/// The eight per-kind registries built while parsing a molecule or applying reaction deltas.
#[derive(Debug, Default)]
pub(crate) struct EntityRegistry {
    atoms: NamedRegistry<AtomId>,
    bonds: KeyedRegistry<BondId, [AtomId; 2]>,
    dative_bonds: KeyedRegistry<DativeBondId, (BTreeSet<AtomId>, AtomId)>,
    aromatic_systems: KeyedRegistry<AromaticSystemId, BTreeSet<AtomId>>,
    multicenter_bonds: KeyedRegistry<MulticenterBondId, BTreeSet<AtomId>>,
    noncovalent_bonds: KeyedRegistry<NoncovalentBondId, [AtomId; 2]>,
    stereo_atoms: KeyedRegistry<StereoAtomId, AtomId>,
    stereo_bonds: KeyedRegistry<StereoBondId, BondId>,
}

impl EntityRegistry {
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
    ) -> StereoAtomId {
        self.stereo_atoms.register(name, site)
    }

    pub(crate) fn register_stereo_bond(
        &mut self,
        name: Option<String>,
        site: BondId,
    ) -> StereoBondId {
        self.stereo_bonds.register(name, site)
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

    pub(crate) fn stereo_atom_by_participants(&self, site: AtomId) -> Option<StereoAtomId> {
        self.stereo_atoms.by_participants(&site)
    }

    pub(crate) fn stereo_bond_by_participants(&self, site: BondId) -> Option<StereoBondId> {
        self.stereo_bonds.by_participants(&site)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_entity_registry_register_atom() {
        let mut registry = EntityRegistry::default();
        assert_eq!(registry.register_atom(None), AtomId(0));
        assert_eq!(registry.register_atom(Some("c1".into())), AtomId(1));
        assert_eq!(registry.register_atom(None), AtomId(2));
        assert_eq!(registry.atom_count(), 3);
        assert_eq!(registry.atom_by_name("c1"), Some(AtomId(1)));
        assert_eq!(registry.atom_by_name("nope"), None);
    }

    #[rstest]
    fn test_entity_registry_register_bond() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_bond(None, AtomId(2), AtomId(0)),
            BondId(0)
        );
        assert_eq!(
            registry.register_bond(Some("b1".into()), AtomId(1), AtomId(3)),
            BondId(1)
        );
        assert_eq!(registry.bond_count(), 2);
        assert_eq!(registry.bond_by_name("b1"), Some(BondId(1)));
        // Endpoint key is order-independent.
        assert_eq!(
            registry.bond_by_participants(AtomId(0), AtomId(2)),
            Some(BondId(0))
        );
        assert_eq!(
            registry.bond_by_participants(AtomId(2), AtomId(0)),
            Some(BondId(0))
        );
        assert_eq!(registry.bond_by_participants(AtomId(0), AtomId(4)), None);
    }

    #[rstest]
    fn test_entity_registry_register_dative_bond() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0)),
            DativeBondId(0)
        );
        assert_eq!(registry.dative_bond_count(), 1);
        // Donor set is order-independent; the acceptor is distinguished.
        assert_eq!(
            registry.dative_bond_by_participants(&[AtomId(2), AtomId(1)], AtomId(0)),
            Some(DativeBondId(0))
        );
        assert_eq!(
            registry.dative_bond_by_participants(&[AtomId(1), AtomId(2)], AtomId(3)),
            None
        );
    }

    #[rstest]
    fn test_entity_registry_register_aromatic_system() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)]),
            AromaticSystemId(0)
        );
        assert_eq!(registry.aromatic_system_count(), 1);
        // Atom set is order-independent.
        assert_eq!(
            registry.aromatic_system_by_participants(&[AtomId(0), AtomId(1), AtomId(2)]),
            Some(AromaticSystemId(0))
        );
        assert_eq!(
            registry.aromatic_system_by_participants(&[AtomId(0), AtomId(1)]),
            None
        );
    }

    #[rstest]
    fn test_entity_registry_register_multicenter_bond() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry
                .register_multicenter_bond(Some("m".into()), &[AtomId(0), AtomId(1), AtomId(2)]),
            MulticenterBondId(0)
        );
        assert_eq!(registry.multicenter_bond_count(), 1);
        assert_eq!(
            registry.multicenter_bond_by_name("m"),
            Some(MulticenterBondId(0))
        );
        assert_eq!(
            registry.multicenter_bond_by_participants(&[AtomId(2), AtomId(1), AtomId(0)]),
            Some(MulticenterBondId(0))
        );
    }

    #[rstest]
    fn test_entity_registry_register_noncovalent_bond() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_noncovalent_bond(None, AtomId(3), AtomId(1)),
            NoncovalentBondId(0)
        );
        assert_eq!(registry.noncovalent_bond_count(), 1);
        assert_eq!(
            registry.noncovalent_bond_by_participants(AtomId(1), AtomId(3)),
            Some(NoncovalentBondId(0))
        );
        assert_eq!(
            registry.noncovalent_bond_by_participants(AtomId(1), AtomId(2)),
            None
        );
    }

    #[rstest]
    fn test_entity_registry_register_stereo_atom() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_stereo_atom(None, AtomId(4)),
            StereoAtomId(0)
        );
        assert_eq!(registry.stereo_atom_count(), 1);
        assert_eq!(
            registry.stereo_atom_by_participants(AtomId(4)),
            Some(StereoAtomId(0))
        );
        assert_eq!(registry.stereo_atom_by_participants(AtomId(0)), None);
    }

    #[rstest]
    fn test_entity_registry_register_stereo_bond() {
        let mut registry = EntityRegistry::default();
        assert_eq!(
            registry.register_stereo_bond(None, BondId(2)),
            StereoBondId(0)
        );
        assert_eq!(registry.stereo_bond_count(), 1);
        assert_eq!(
            registry.stereo_bond_by_participants(BondId(2)),
            Some(StereoBondId(0))
        );
        assert_eq!(registry.stereo_bond_by_participants(BondId(0)), None);
    }
}
