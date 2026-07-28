//! A molecule's parse-time context: the roundtrip [`MoleculeMetadata`], per-kind running counts,
//! and participant indexes. It grows while parsing a molecule (or applying reaction deltas).
//! Index bounds are checked as entities are registered (not only at the end), and structural refs
//! (a non-atom entity named by its constituent atoms/bonds) resolve against it.
//!
//! Cost splits by kind: atoms carry no participant lookup (the base kind), bonds an O(1)
//! `(min,max) → id` endpoint map (a bond is a graph edge), the five overlays a participant index
//! over their small collections. §4.1 uniqueness (no two same-constituent entries) makes every
//! participant lookup a ≤1 hit.

use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;
use std::marker::PhantomData;

use super::atom::AtomDsl;
use super::error::ParseError;
use super::metadata::{MetadataError, MoleculeMetadata};
use crate::ast::entity::Entity;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;
use crate::ast::molecule::MoleculeAst;

/// The metadata, counters, and participant indexes built while parsing a molecule or applying
/// reaction deltas. Atoms need only a counter; the seven non-atom kinds add participant lookup.
#[derive(Debug, Default)]
pub struct MoleculeContext {
    metadata: MoleculeMetadata,
    atoms: EntityCounter<AtomId>,
    bonds: EntityRegistry<BondId, [AtomId; 2]>,
    dative_bonds: EntityRegistry<DativeBondId, (BTreeSet<AtomId>, AtomId)>,
    aromatic_systems: EntityRegistry<AromaticSystemId, BTreeSet<AtomId>>,
    multicenter_bonds: EntityRegistry<MulticenterBondId, BTreeSet<AtomId>>,
    noncovalent_bonds: EntityRegistry<NoncovalentBondId, [AtomId; 2]>,
    stereo_atoms: EntityRegistry<StereoAtomId, (AtomId, Vec<StereoLigand>)>,
    stereo_bonds: EntityRegistry<StereoBondId, (BondId, Vec<StereoLigand>)>,
}

/// The parse-time **resolution** query surface — everything `ref::resolve` and the constraint /
/// relational resolvers read to turn a surface ref into a numerical AST id
/// (keyword / index / participants → id). Written once, generic over this
/// trait; implemented for [`MoleculeContext`] (a molecule), a reaction's
/// `ReactionContext`, and sub-pattern contexts. The inverse direction
/// (numerical AST id → keyword, for rendering) is the separate `Metadata`
/// trait.
pub trait Namespace {
    fn atom_count(&self) -> usize;
    fn bond_count(&self) -> usize;
    fn dative_bond_count(&self) -> usize;
    fn aromatic_system_count(&self) -> usize;
    fn multicenter_bond_count(&self) -> usize;
    fn noncovalent_bond_count(&self) -> usize;
    fn stereo_atom_count(&self) -> usize;
    fn stereo_bond_count(&self) -> usize;

    fn find_atom_by_keyword(&self, keyword: &str) -> Option<AtomId>;
    fn find_bond_by_keyword(&self, keyword: &str) -> Option<BondId>;
    fn find_dative_bond_by_keyword(&self, keyword: &str) -> Option<DativeBondId>;
    fn find_aromatic_system_by_keyword(&self, keyword: &str) -> Option<AromaticSystemId>;
    fn find_multicenter_bond_by_keyword(&self, keyword: &str) -> Option<MulticenterBondId>;
    fn find_noncovalent_bond_by_keyword(&self, keyword: &str) -> Option<NoncovalentBondId>;
    fn find_stereo_atom_by_keyword(&self, keyword: &str) -> Option<StereoAtomId>;
    fn find_stereo_bond_by_keyword(&self, keyword: &str) -> Option<StereoBondId>;

    fn find_bond_by_participants(&self, first: AtomId, second: AtomId) -> Option<BondId>;
    fn find_dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId>;
    fn find_aromatic_system_by_participants(&self, atoms: &[AtomId]) -> Option<AromaticSystemId>;
    fn find_multicenter_bond_by_participants(&self, atoms: &[AtomId]) -> Option<MulticenterBondId>;
    fn find_noncovalent_bond_by_participants(
        &self,
        first: AtomId,
        second: AtomId,
    ) -> Option<NoncovalentBondId>;
    fn find_stereo_atom_by_participants(
        &self,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Option<StereoAtomId>;
    fn find_stereo_bond_by_participants(
        &self,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Option<StereoBondId>;

    /// Whether `keyword` is already taken as an entity keyword or atom-alias name.
    fn contains_keyword(&self, keyword: &str) -> bool;

    /// The atom-spec template registered under alias `name`, for resolving `<alias>` atom specs.
    fn find_atom_alias(&self, name: &str) -> Option<&AtomDsl>;
}

impl MoleculeContext {
    /// A context continuing another's id space: each kind's count starts at `other`'s count (so
    /// `register_*` hands out ids following it), with empty metadata and participant indexes. It
    /// holds only the entities registered into it and is used for reaction deltas over an lhs.
    pub(crate) fn continuation(other: &MoleculeContext) -> Self {
        Self {
            metadata: MoleculeMetadata::new(),
            atoms: EntityCounter::from_count(other.atoms.count()),
            bonds: EntityRegistry::from_count(other.bonds.count()),
            dative_bonds: EntityRegistry::from_count(other.dative_bonds.count()),
            aromatic_systems: EntityRegistry::from_count(other.aromatic_systems.count()),
            multicenter_bonds: EntityRegistry::from_count(other.multicenter_bonds.count()),
            noncovalent_bonds: EntityRegistry::from_count(other.noncovalent_bonds.count()),
            stereo_atoms: EntityRegistry::from_count(other.stereo_atoms.count()),
            stereo_bonds: EntityRegistry::from_count(other.stereo_bonds.count()),
        }
    }

    /// The context of an already-resolved molecule: every entity registered anonymously (no
    /// keyword) with its participants, so a sub-pattern's index and structural refs resolve against
    /// it. The ids are anonymous, so registration cannot collide.
    pub fn from_ast(ast: &MoleculeAst) -> Self {
        let free = "anonymous entity registration never collides";
        let mut context = Self::default();
        for _ in ast.atoms().ids() {
            context.register_atom(None).expect(free);
        }
        for view in ast.bonds().iter() {
            let [a, b] = view.atom_ids();
            context.register_bond(None, a, b).expect(free);
        }
        for view in ast.dative_bonds().iter() {
            let donors: Vec<AtomId> = view.donor_ids().collect();
            context
                .register_dative_bond(None, &donors, view.acceptor_id())
                .expect(free);
        }
        for view in ast.aromatic_systems().iter() {
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            context.register_aromatic_system(None, &atoms).expect(free);
        }
        for view in ast.multicenter_bonds().iter() {
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            context.register_multicenter_bond(None, &atoms).expect(free);
        }
        for view in ast.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            context.register_noncovalent_bond(None, a, b).expect(free);
        }
        for view in ast.stereo_atoms().iter() {
            let ligands: Vec<StereoLigand> = view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect();
            context
                .register_stereo_atom(None, view.site_id(), &ligands)
                .expect(free);
        }
        for view in ast.stereo_bonds().iter() {
            let ligands: Vec<StereoLigand> = view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect();
            context
                .register_stereo_bond(None, view.site_id(), &ligands)
                .expect(free);
        }
        context
    }

    /// Metadata accumulated by this context.
    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    /// Consume the parse-time indexes and return the accumulated roundtrip metadata.
    pub fn into_metadata(self) -> MoleculeMetadata {
        self.metadata
    }

    fn set_keyword(&mut self, entity: Entity, keyword: Option<String>) -> Result<(), ParseError> {
        if let Some(keyword) = keyword {
            self.metadata
                .set_keyword(entity, keyword)
                .map_err(metadata_parse_error)?;
        }
        Ok(())
    }

    pub(crate) fn register_atom(&mut self, keyword: Option<String>) -> Result<AtomId, ParseError> {
        let id = self.atoms.next_id();
        self.set_keyword(Entity::Atom(id), keyword)?;
        Ok(self.atoms.register())
    }

    pub(crate) fn register_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<BondId, ParseError> {
        let id = self.bonds.next_id();
        self.set_keyword(Entity::Bond(id), keyword)?;
        Ok(self.bonds.register(atom_pair_key(a, b)))
    }

    pub(crate) fn register_dative_bond(
        &mut self,
        keyword: Option<String>,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Result<DativeBondId, ParseError> {
        let id = self.dative_bonds.next_id();
        self.set_keyword(Entity::DativeBond(id), keyword)?;
        Ok(self
            .dative_bonds
            .register((donors.iter().copied().collect(), acceptor)))
    }

    pub(crate) fn register_aromatic_system(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<AromaticSystemId, ParseError> {
        let id = self.aromatic_systems.next_id();
        self.set_keyword(Entity::AromaticSystem(id), keyword)?;
        Ok(self
            .aromatic_systems
            .register(atoms.iter().copied().collect()))
    }

    pub(crate) fn register_multicenter_bond(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<MulticenterBondId, ParseError> {
        let id = self.multicenter_bonds.next_id();
        self.set_keyword(Entity::MulticenterBond(id), keyword)?;
        Ok(self
            .multicenter_bonds
            .register(atoms.iter().copied().collect()))
    }

    pub(crate) fn register_noncovalent_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<NoncovalentBondId, ParseError> {
        let id = self.noncovalent_bonds.next_id();
        self.set_keyword(Entity::NoncovalentBond(id), keyword)?;
        Ok(self.noncovalent_bonds.register(atom_pair_key(a, b)))
    }

    pub(crate) fn register_stereo_atom(
        &mut self,
        keyword: Option<String>,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Result<StereoAtomId, ParseError> {
        let id = self.stereo_atoms.next_id();
        self.set_keyword(Entity::StereoAtom(id), keyword)?;
        Ok(self.stereo_atoms.register((site, ligand_multiset(ligands))))
    }

    pub(crate) fn register_stereo_bond(
        &mut self,
        keyword: Option<String>,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Result<StereoBondId, ParseError> {
        let id = self.stereo_bonds.next_id();
        self.set_keyword(Entity::StereoBond(id), keyword)?;
        Ok(self.stereo_bonds.register((site, ligand_multiset(ligands))))
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

    pub(crate) fn find_atom_by_keyword(&self, keyword: &str) -> Option<AtomId> {
        match self.metadata.entity(keyword) {
            Some(Entity::Atom(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_bond_by_keyword(&self, keyword: &str) -> Option<BondId> {
        match self.metadata.entity(keyword) {
            Some(Entity::Bond(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_dative_bond_by_keyword(&self, keyword: &str) -> Option<DativeBondId> {
        match self.metadata.entity(keyword) {
            Some(Entity::DativeBond(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_aromatic_system_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<AromaticSystemId> {
        match self.metadata.entity(keyword) {
            Some(Entity::AromaticSystem(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_multicenter_bond_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<MulticenterBondId> {
        match self.metadata.entity(keyword) {
            Some(Entity::MulticenterBond(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_noncovalent_bond_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<NoncovalentBondId> {
        match self.metadata.entity(keyword) {
            Some(Entity::NoncovalentBond(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_stereo_atom_by_keyword(&self, keyword: &str) -> Option<StereoAtomId> {
        match self.metadata.entity(keyword) {
            Some(Entity::StereoAtom(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_stereo_bond_by_keyword(&self, keyword: &str) -> Option<StereoBondId> {
        match self.metadata.entity(keyword) {
            Some(Entity::StereoBond(id)) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn find_bond_by_participants(
        &self,
        first: AtomId,
        second: AtomId,
    ) -> Option<BondId> {
        self.bonds
            .find_by_participants(&atom_pair_key(first, second))
    }

    pub(crate) fn find_dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId> {
        self.dative_bonds
            .find_by_participants(&(donors.iter().copied().collect(), acceptor))
    }

    pub(crate) fn find_aromatic_system_by_participants(
        &self,
        atoms: &[AtomId],
    ) -> Option<AromaticSystemId> {
        self.aromatic_systems
            .find_by_participants(&atoms.iter().copied().collect())
    }

    pub(crate) fn find_multicenter_bond_by_participants(
        &self,
        atoms: &[AtomId],
    ) -> Option<MulticenterBondId> {
        self.multicenter_bonds
            .find_by_participants(&atoms.iter().copied().collect())
    }

    pub(crate) fn find_noncovalent_bond_by_participants(
        &self,
        first: AtomId,
        second: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds
            .find_by_participants(&atom_pair_key(first, second))
    }

    pub(crate) fn find_stereo_atom_by_participants(
        &self,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Option<StereoAtomId> {
        self.stereo_atoms
            .find_by_participants(&(site, ligand_multiset(ligands)))
    }

    pub(crate) fn find_stereo_bond_by_participants(
        &self,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Option<StereoBondId> {
        self.stereo_bonds
            .find_by_participants(&(site, ligand_multiset(ligands)))
    }

    pub(crate) fn register_atom_alias(
        &mut self,
        name: String,
        dsl: AtomDsl,
    ) -> Result<(), ParseError> {
        self.metadata
            .add_atom_alias(name, dsl)
            .map_err(metadata_parse_error)
    }
}

impl Namespace for MoleculeContext {
    fn atom_count(&self) -> usize {
        self.atom_count()
    }
    fn bond_count(&self) -> usize {
        self.bond_count()
    }
    fn dative_bond_count(&self) -> usize {
        self.dative_bond_count()
    }
    fn aromatic_system_count(&self) -> usize {
        self.aromatic_system_count()
    }
    fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bond_count()
    }
    fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bond_count()
    }
    fn stereo_atom_count(&self) -> usize {
        self.stereo_atom_count()
    }
    fn stereo_bond_count(&self) -> usize {
        self.stereo_bond_count()
    }

    fn find_atom_by_keyword(&self, keyword: &str) -> Option<AtomId> {
        self.find_atom_by_keyword(keyword)
    }
    fn find_bond_by_keyword(&self, keyword: &str) -> Option<BondId> {
        self.find_bond_by_keyword(keyword)
    }
    fn find_dative_bond_by_keyword(&self, keyword: &str) -> Option<DativeBondId> {
        self.find_dative_bond_by_keyword(keyword)
    }
    fn find_aromatic_system_by_keyword(&self, keyword: &str) -> Option<AromaticSystemId> {
        self.find_aromatic_system_by_keyword(keyword)
    }
    fn find_multicenter_bond_by_keyword(&self, keyword: &str) -> Option<MulticenterBondId> {
        self.find_multicenter_bond_by_keyword(keyword)
    }
    fn find_noncovalent_bond_by_keyword(&self, keyword: &str) -> Option<NoncovalentBondId> {
        self.find_noncovalent_bond_by_keyword(keyword)
    }
    fn find_stereo_atom_by_keyword(&self, keyword: &str) -> Option<StereoAtomId> {
        self.find_stereo_atom_by_keyword(keyword)
    }
    fn find_stereo_bond_by_keyword(&self, keyword: &str) -> Option<StereoBondId> {
        self.find_stereo_bond_by_keyword(keyword)
    }

    fn find_bond_by_participants(&self, first: AtomId, second: AtomId) -> Option<BondId> {
        self.find_bond_by_participants(first, second)
    }
    fn find_dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId> {
        self.find_dative_bond_by_participants(donors, acceptor)
    }
    fn find_aromatic_system_by_participants(&self, atoms: &[AtomId]) -> Option<AromaticSystemId> {
        self.find_aromatic_system_by_participants(atoms)
    }
    fn find_multicenter_bond_by_participants(&self, atoms: &[AtomId]) -> Option<MulticenterBondId> {
        self.find_multicenter_bond_by_participants(atoms)
    }
    fn find_noncovalent_bond_by_participants(
        &self,
        first: AtomId,
        second: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.find_noncovalent_bond_by_participants(first, second)
    }
    fn find_stereo_atom_by_participants(
        &self,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Option<StereoAtomId> {
        self.find_stereo_atom_by_participants(site, ligands)
    }
    fn find_stereo_bond_by_participants(
        &self,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Option<StereoBondId> {
        self.find_stereo_bond_by_participants(site, ligands)
    }

    fn contains_keyword(&self, keyword: &str) -> bool {
        self.find_atom_by_keyword(keyword).is_some()
            || self.find_bond_by_keyword(keyword).is_some()
            || self.find_dative_bond_by_keyword(keyword).is_some()
            || self.find_aromatic_system_by_keyword(keyword).is_some()
            || self.find_multicenter_bond_by_keyword(keyword).is_some()
            || self.find_noncovalent_bond_by_keyword(keyword).is_some()
            || self.find_stereo_atom_by_keyword(keyword).is_some()
            || self.find_stereo_bond_by_keyword(keyword).is_some()
            || self.metadata.atom_alias(keyword).is_some()
    }

    fn find_atom_alias(&self, name: &str) -> Option<&AtomDsl> {
        self.metadata.atom_alias(name)
    }
}

/// Running count for one entity kind. Atoms use it directly.
#[derive(Debug)]
struct EntityCounter<Id> {
    count: usize,
    marker: PhantomData<Id>,
}

impl<Id> Default for EntityCounter<Id> {
    fn default() -> Self {
        Self {
            count: 0,
            marker: PhantomData,
        }
    }
}

impl<Id: Copy + From<usize>> EntityCounter<Id> {
    fn next_id(&self) -> Id {
        Id::from(self.count)
    }

    fn register(&mut self) -> Id {
        let id = self.next_id();
        self.count += 1;
        id
    }

    fn count(&self) -> usize {
        self.count
    }

    /// A counter whose next id starts at `count`.
    fn from_count(count: usize) -> Self {
        Self {
            count,
            marker: PhantomData,
        }
    }
}

/// Count + participant lookup for one non-atom entity kind. `Key` is the entity's
/// canonical participant key (a normalized endpoint pair, atom set, donor-set + acceptor, or stereo
/// site + ligand multiset); §4.1 uniqueness makes it injective, so a hit is unique.
#[derive(Debug)]
struct EntityRegistry<Id, Key> {
    count: usize,
    by_participants: HashMap<Key, Id>,
}

impl<Id, Key> Default for EntityRegistry<Id, Key> {
    fn default() -> Self {
        Self {
            count: 0,
            by_participants: HashMap::new(),
        }
    }
}

impl<Id: Copy + From<usize>, Key: Eq + Hash> EntityRegistry<Id, Key> {
    fn next_id(&self) -> Id {
        Id::from(self.count)
    }

    /// Reserve the next id and record its canonical participant key.
    fn register(&mut self, key: Key) -> Id {
        let id = self.next_id();
        self.count += 1;
        self.by_participants.insert(key, id);
        id
    }

    fn find_by_participants(&self, key: &Key) -> Option<Id> {
        self.by_participants.get(key).copied()
    }

    fn count(&self) -> usize {
        self.count
    }

    fn from_count(count: usize) -> Self {
        Self {
            count,
            by_participants: HashMap::new(),
        }
    }
}

fn metadata_parse_error(error: MetadataError) -> ParseError {
    match error {
        MetadataError::DuplicateKeyword(keyword) => ParseError::DuplicateKeyword(keyword),
        MetadataError::DuplicateAtomAlias(_) => ParseError::InvalidValue(
            "atom-aliases must be bijective: two names map to the same atom".into(),
        ),
        MetadataError::EntityOutOfRange(entity) => {
            ParseError::InvalidValue(format!("metadata entity is out of range: {entity}"))
        }
        MetadataError::EntityNotAdded(entity) => ParseError::InvalidValue(format!(
            "metadata entity is not introduced by an add delta: {entity}"
        )),
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

/// The canonical key of a stereo element's ligand frame: the ligand sorted multiset
/// (can have duplicate virtual ligands).
fn ligand_multiset(ligands: &[StereoLigand]) -> Vec<StereoLigand> {
    let mut ligands = ligands.to_vec();
    ligands.sort_unstable();
    ligands
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::ast::ligand::StereoLigandKind;

    #[rstest]
    fn test_molecule_context_continuation() {
        let mut lhs = MoleculeContext::default();
        lhs.register_atom(None).unwrap();
        lhs.register_atom(Some("c1".into())).unwrap();
        lhs.register_bond(None, AtomId(0), AtomId(1)).unwrap();

        let mut delta = MoleculeContext::continuation(&lhs);
        // Counts continue the lhs id space; keywords / participants start empty.
        assert_eq!(delta.atom_count(), 2);
        assert_eq!(delta.bond_count(), 1);
        assert_eq!(delta.find_atom_by_keyword("c1"), None);
        assert_eq!(delta.find_bond_by_participants(AtomId(0), AtomId(1)), None);
        // `register_*` hands out global ids following the lhs.
        assert_eq!(delta.register_atom(None).unwrap(), AtomId(2));
        assert_eq!(
            delta.register_bond(None, AtomId(0), AtomId(2)).unwrap(),
            BondId(1)
        );
        assert_eq!(delta.atom_count(), 3);
        assert_eq!(
            delta.find_bond_by_participants(AtomId(0), AtomId(2)),
            Some(BondId(1))
        );
    }

    #[rstest]
    fn test_molecule_context_register_atom() {
        let mut context = MoleculeContext::default();
        assert_eq!(context.register_atom(None).unwrap(), AtomId(0));
        assert_eq!(context.register_atom(Some("c1".into())).unwrap(), AtomId(1));
        assert_eq!(context.register_atom(None).unwrap(), AtomId(2));
        assert_eq!(context.atom_count(), 3);
        assert_eq!(context.find_atom_by_keyword("c1"), Some(AtomId(1)));
        assert_eq!(context.find_atom_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_atom_error() {
        let mut context = MoleculeContext::default();
        context.register_atom(Some("a".into())).unwrap();
        assert_eq!(
            context.register_atom(Some("a".into())).unwrap_err(),
            ParseError::DuplicateKeyword("a".into())
        );
        assert_eq!(
            context
                .register_bond(Some("a".into()), AtomId(0), AtomId(0))
                .unwrap_err(),
            ParseError::DuplicateKeyword("a".into())
        );
        assert_eq!(context.atom_count(), 1);
        assert_eq!(context.bond_count(), 0);
        assert_eq!(context.find_atom_by_keyword("a"), Some(AtomId(0)));
        assert_eq!(
            context.find_bond_by_participants(AtomId(0), AtomId(0)),
            None
        );
    }

    #[rstest]
    fn test_molecule_context_register_bond() {
        let mut context = MoleculeContext::default();
        assert_eq!(
            context.register_bond(None, AtomId(2), AtomId(0)).unwrap(),
            BondId(0)
        );
        assert_eq!(
            context
                .register_bond(Some("b1".into()), AtomId(1), AtomId(3))
                .unwrap(),
            BondId(1)
        );
        assert_eq!(context.bond_count(), 2);
        assert_eq!(context.find_bond_by_keyword("b1"), Some(BondId(1)));
        assert_eq!(context.find_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_dative_bond() {
        let mut context = MoleculeContext::default();
        assert_eq!(
            context
                .register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0))
                .unwrap(),
            DativeBondId(0)
        );
        assert_eq!(
            context
                .register_dative_bond(Some("d1".into()), &[AtomId(4)], AtomId(3))
                .unwrap(),
            DativeBondId(1)
        );
        assert_eq!(context.dative_bond_count(), 2);
        assert_eq!(
            context.find_dative_bond_by_keyword("d1"),
            Some(DativeBondId(1))
        );
        assert_eq!(context.find_dative_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_aromatic_system() {
        let mut context = MoleculeContext::default();
        assert_eq!(
            context
                .register_aromatic_system(None, &[AtomId(0), AtomId(1), AtomId(2)])
                .unwrap(),
            AromaticSystemId(0)
        );
        assert_eq!(
            context
                .register_aromatic_system(Some("a1".into()), &[AtomId(3), AtomId(4)])
                .unwrap(),
            AromaticSystemId(1)
        );
        assert_eq!(context.aromatic_system_count(), 2);
        assert_eq!(
            context.find_aromatic_system_by_keyword("a1"),
            Some(AromaticSystemId(1))
        );
        assert_eq!(context.find_aromatic_system_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_multicenter_bond() {
        let mut context = MoleculeContext::default();
        assert_eq!(
            context
                .register_multicenter_bond(Some("m".into()), &[AtomId(0), AtomId(1), AtomId(2)])
                .unwrap(),
            MulticenterBondId(0)
        );
        assert_eq!(
            context
                .register_multicenter_bond(None, &[AtomId(3), AtomId(4)])
                .unwrap(),
            MulticenterBondId(1)
        );
        assert_eq!(context.multicenter_bond_count(), 2);
        assert_eq!(
            context.find_multicenter_bond_by_keyword("m"),
            Some(MulticenterBondId(0))
        );
        assert_eq!(context.find_multicenter_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_noncovalent_bond() {
        let mut context = MoleculeContext::default();
        assert_eq!(
            context
                .register_noncovalent_bond(None, AtomId(3), AtomId(1))
                .unwrap(),
            NoncovalentBondId(0)
        );
        assert_eq!(
            context
                .register_noncovalent_bond(Some("n1".into()), AtomId(0), AtomId(4))
                .unwrap(),
            NoncovalentBondId(1)
        );
        assert_eq!(context.noncovalent_bond_count(), 2);
        assert_eq!(
            context.find_noncovalent_bond_by_keyword("n1"),
            Some(NoncovalentBondId(1))
        );
        assert_eq!(context.find_noncovalent_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_stereo_atom() {
        let mut context = MoleculeContext::default();
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        assert_eq!(
            context
                .register_stereo_atom(None, AtomId(4), &ligands)
                .unwrap(),
            StereoAtomId(0)
        );
        assert_eq!(
            context
                .register_stereo_atom(Some("s1".into()), AtomId(5), &ligands)
                .unwrap(),
            StereoAtomId(1)
        );
        assert_eq!(context.stereo_atom_count(), 2);
        assert_eq!(
            context.find_stereo_atom_by_keyword("s1"),
            Some(StereoAtomId(1))
        );
        assert_eq!(context.find_stereo_atom_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_stereo_bond() {
        let mut context = MoleculeContext::default();
        let ligands = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        assert_eq!(
            context
                .register_stereo_bond(None, BondId(2), &ligands)
                .unwrap(),
            StereoBondId(0)
        );
        assert_eq!(
            context
                .register_stereo_bond(Some("sb1".into()), BondId(0), &ligands)
                .unwrap(),
            StereoBondId(1)
        );
        assert_eq!(context.stereo_bond_count(), 2);
        assert_eq!(
            context.find_stereo_bond_by_keyword("sb1"),
            Some(StereoBondId(1))
        );
        assert_eq!(context.find_stereo_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_context_register_atom_alias() {
        let mut context = MoleculeContext::default();
        let dsl = "C".parse::<AtomDsl>().unwrap();
        context
            .register_atom_alias("me".into(), dsl.clone())
            .unwrap();
        assert_eq!(context.find_atom_alias("me"), Some(&dsl));
        assert!(context.contains_keyword("me"));
        assert_eq!(context.find_atom_alias("nope"), None);
    }

    #[rstest]
    #[case::name_taken_by_entity("a", "N", ParseError::DuplicateKeyword("a".into()))]
    #[case::name_taken_by_alias("me", "O", ParseError::DuplicateKeyword("me".into()))]
    #[case::duplicate_spec(
        "me2",
        "C",
        ParseError::InvalidValue(
            "atom-aliases must be bijective: two names map to the same atom".into()
        )
    )]
    fn test_molecule_context_register_atom_alias_error(
        #[case] name: &str,
        #[case] spec: &str,
        #[case] expected: ParseError,
    ) {
        let mut context = MoleculeContext::default();
        context.register_atom(Some("a".into())).unwrap();
        context
            .register_atom_alias("me".into(), "C".parse::<AtomDsl>().unwrap())
            .unwrap();
        let before = context.metadata().clone();
        assert_eq!(
            context
                .register_atom_alias(name.into(), spec.parse::<AtomDsl>().unwrap())
                .unwrap_err(),
            expected
        );
        assert_eq!(context.metadata(), &before);
    }

    #[rstest]
    #[case::atom("a", true)]
    #[case::bond("b", true)]
    #[case::dative_bond("d", true)]
    #[case::aromatic_system("ar", true)]
    #[case::multicenter_bond("m", true)]
    #[case::noncovalent_bond("nc", true)]
    #[case::stereo_atom("sa", true)]
    #[case::stereo_bond("sb", true)]
    #[case::alias("al", true)]
    #[case::absent("nope", false)]
    fn test_molecule_context_contains_keyword(#[case] keyword: &str, #[case] expected: bool) {
        let mut context = MoleculeContext::default();
        context.register_atom(Some("a".into())).unwrap();
        context
            .register_bond(Some("b".into()), AtomId(0), AtomId(0))
            .unwrap();
        context
            .register_dative_bond(Some("d".into()), &[AtomId(0)], AtomId(0))
            .unwrap();
        context
            .register_aromatic_system(Some("ar".into()), &[AtomId(0)])
            .unwrap();
        context
            .register_multicenter_bond(Some("m".into()), &[AtomId(0)])
            .unwrap();
        context
            .register_noncovalent_bond(Some("nc".into()), AtomId(0), AtomId(0))
            .unwrap();
        context
            .register_stereo_atom(Some("sa".into()), AtomId(0), &[])
            .unwrap();
        context
            .register_stereo_bond(Some("sb".into()), BondId(0), &[])
            .unwrap();
        context
            .register_atom_alias("al".into(), "C".parse::<AtomDsl>().unwrap())
            .unwrap();
        assert_eq!(context.contains_keyword(keyword), expected);
    }

    #[rstest]
    #[case::forward(AtomId(0), AtomId(2), Some(BondId(0)))]
    #[case::reversed(AtomId(2), AtomId(0), Some(BondId(0)))]
    #[case::absent(AtomId(0), AtomId(4), None)]
    fn test_molecule_context_find_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<BondId>,
    ) {
        let mut context = MoleculeContext::default();
        context.register_bond(None, AtomId(2), AtomId(0)).unwrap();
        assert_eq!(context.find_bond_by_participants(a, b), expected);
    }

    #[rstest]
    #[case::donors_reordered(&[AtomId(2), AtomId(1)], AtomId(0), Some(DativeBondId(0)))]
    #[case::wrong_acceptor(&[AtomId(1), AtomId(2)], AtomId(3), None)]
    #[case::wrong_donors(&[AtomId(1), AtomId(3)], AtomId(0), None)]
    fn test_molecule_context_find_dative_bond_by_participants(
        #[case] donors: &[AtomId],
        #[case] acceptor: AtomId,
        #[case] expected: Option<DativeBondId>,
    ) {
        let mut context = MoleculeContext::default();
        context
            .register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0))
            .unwrap();
        assert_eq!(
            context.find_dative_bond_by_participants(donors, acceptor),
            expected
        );
    }

    #[rstest]
    #[case::reordered(&[AtomId(0), AtomId(1), AtomId(2)], Some(AromaticSystemId(0)))]
    #[case::subset(&[AtomId(0), AtomId(1)], None)]
    #[case::superset(&[AtomId(0), AtomId(1), AtomId(2), AtomId(3)], None)]
    fn test_molecule_context_find_aromatic_system_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let mut context = MoleculeContext::default();
        context
            .register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)])
            .unwrap();
        assert_eq!(
            context.find_aromatic_system_by_participants(atoms),
            expected
        );
    }

    #[rstest]
    #[case::reordered(&[AtomId(2), AtomId(1), AtomId(0)], Some(MulticenterBondId(0)))]
    #[case::absent(&[AtomId(0), AtomId(1), AtomId(3)], None)]
    fn test_molecule_context_find_multicenter_bond_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<MulticenterBondId>,
    ) {
        let mut context = MoleculeContext::default();
        context
            .register_multicenter_bond(None, &[AtomId(0), AtomId(1), AtomId(2)])
            .unwrap();
        assert_eq!(
            context.find_multicenter_bond_by_participants(atoms),
            expected
        );
    }

    #[rstest]
    #[case::reversed(AtomId(1), AtomId(3), Some(NoncovalentBondId(0)))]
    #[case::absent(AtomId(1), AtomId(2), None)]
    fn test_molecule_context_find_noncovalent_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<NoncovalentBondId>,
    ) {
        let mut context = MoleculeContext::default();
        context
            .register_noncovalent_bond(None, AtomId(3), AtomId(1))
            .unwrap();
        assert_eq!(
            context.find_noncovalent_bond_by_participants(a, b),
            expected
        );
    }

    #[rstest]
    #[case::reordered_ligands(AtomId(4), &[AtomId(2), AtomId(1)], Some(StereoAtomId(0)))]
    #[case::wrong_ligands(AtomId(4), &[AtomId(1)], None)]
    #[case::wrong_site(AtomId(0), &[AtomId(1), AtomId(2)], None)]
    fn test_molecule_context_find_stereo_atom_by_participants(
        #[case] site: AtomId,
        #[case] ligand_atoms: &[AtomId],
        #[case] expected: Option<StereoAtomId>,
    ) {
        let mut context = MoleculeContext::default();
        // All test ligands are `Atom`-kind; only the atom set varies per case.
        let registered = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        context
            .register_stereo_atom(None, AtomId(4), &registered)
            .unwrap();
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            context.find_stereo_atom_by_participants(site, &query),
            expected
        );
    }

    #[rstest]
    #[case::matching(BondId(2), &[AtomId(3)], Some(StereoBondId(0)))]
    #[case::wrong_site(BondId(0), &[AtomId(3)], None)]
    #[case::empty_ligands(BondId(2), &[], None)]
    fn test_molecule_context_find_stereo_bond_by_participants(
        #[case] site: BondId,
        #[case] ligand_atoms: &[AtomId],
        #[case] expected: Option<StereoBondId>,
    ) {
        let mut context = MoleculeContext::default();
        let registered = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        context
            .register_stereo_bond(None, BondId(2), &registered)
            .unwrap();
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            context.find_stereo_bond_by_participants(site, &query),
            expected
        );
    }
}
