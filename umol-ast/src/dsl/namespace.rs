//! A molecule's parse-time **namespace**: per entity kind, a running count, an id-keyword lookup, and
//! a participant lookup, plus the atom-alias table. Grown while parsing a molecule (or applying
//! reaction deltas); the roundtrip subset projects out as [`MoleculeMetadata`]. Index bounds are
//! checked as entities are registered (not only at the end), and structural refs (a non-atom entity
//! named by its constituent atoms/bonds) resolve against it.
//!
//! Cost splits by kind: atoms carry no participant lookup (the base kind), bonds an O(1)
//! `(min,max) → id` endpoint map (a bond is a graph edge), the five overlays a participant index
//! over their small collections. §4.1 uniqueness (no two same-constituent entries) makes every
//! participant lookup a ≤1 hit.

// Wired into molecule parsing (S2b) and the reaction delta loop (S2d); until then the participant /
// keyword query surface is exercised only by the unit tests below.
#![allow(dead_code)]

use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;

use bimap::BiBTreeMap;
use indexmap::IndexMap;

use super::atom::AtomDsl;
use super::error::ParseError;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;
use crate::ast::molecule::MoleculeAst;

/// The eight per-kind registries built while parsing a molecule or applying reaction deltas. Atoms —
/// the base kind, no participants — use a [`KeywordRegistry`]; the seven non-atom kinds use an
/// [`EntityRegistry`], which adds the participant lookup.
#[derive(Debug, Default)]
pub struct MoleculeNamespace {
    atoms: KeywordRegistry<AtomId>,
    bonds: EntityRegistry<BondId, [AtomId; 2]>,
    dative_bonds: EntityRegistry<DativeBondId, (BTreeSet<AtomId>, AtomId)>,
    aromatic_systems: EntityRegistry<AromaticSystemId, BTreeSet<AtomId>>,
    multicenter_bonds: EntityRegistry<MulticenterBondId, BTreeSet<AtomId>>,
    noncovalent_bonds: EntityRegistry<NoncovalentBondId, [AtomId; 2]>,
    stereo_atoms: EntityRegistry<StereoAtomId, (AtomId, Vec<StereoLigand>)>,
    stereo_bonds: EntityRegistry<StereoBondId, (BondId, Vec<StereoLigand>)>,
    /// The bijective atom-alias table (alias name ↔ atom-spec template) — part of the atom name
    /// namespace (an `:id` may not collide with an alias name), so the namespace owns it.
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

/// The parse-time **resolution** query surface — everything `ref::resolve` and the constraint /
/// relational resolvers read to turn a surface ref into an AST id (keyword / index / participants →
/// id). Written once, generic over this trait; implemented for [`MoleculeNamespace`] (a molecule),
/// and later for a reaction's `ReactionNamespace` and a sub-pattern's namespaces. The inverse
/// direction (id → keyword, for rendering) is the separate `Metadata` trait.
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

    fn find_bond_by_participants(&self, a: AtomId, b: AtomId) -> Option<BondId>;
    fn find_dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId>;
    fn find_aromatic_system_by_participants(&self, atoms: &[AtomId]) -> Option<AromaticSystemId>;
    fn find_multicenter_bond_by_participants(&self, atoms: &[AtomId]) -> Option<MulticenterBondId>;
    fn find_noncovalent_bond_by_participants(
        &self,
        a: AtomId,
        b: AtomId,
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

    /// Whether `id` is already taken as any entity's `:id` keyword or an atom-alias name — the
    /// id-uniqueness check the delta / entry loops need across the whole id namespace.
    fn contains_id(&self, id: &str) -> bool;

    /// The atom-spec template registered under alias `name`, for resolving `<alias>` atom specs.
    fn find_atom_alias(&self, name: &str) -> Option<&AtomDsl>;
}

impl MoleculeNamespace {
    /// A namespace continuing another's id space: each kind's count starts at `other`'s count (so
    /// `register_*` hands out ids following it), with empty keyword / participant / alias maps — it
    /// holds only the entities registered into it. Used for a reaction's delta namespace over its lhs.
    pub(crate) fn continuation(other: &MoleculeNamespace) -> Self {
        Self {
            atoms: KeywordRegistry::from_count(other.atoms.count()),
            bonds: EntityRegistry::from_count(other.bonds.count()),
            dative_bonds: EntityRegistry::from_count(other.dative_bonds.count()),
            aromatic_systems: EntityRegistry::from_count(other.aromatic_systems.count()),
            multicenter_bonds: EntityRegistry::from_count(other.multicenter_bonds.count()),
            noncovalent_bonds: EntityRegistry::from_count(other.noncovalent_bonds.count()),
            stereo_atoms: EntityRegistry::from_count(other.stereo_atoms.count()),
            stereo_bonds: EntityRegistry::from_count(other.stereo_bonds.count()),
            atom_aliases: BiBTreeMap::new(),
        }
    }

    /// The namespace of an already-resolved molecule: every entity registered anonymously (no
    /// keyword) with its participants, so a sub-pattern's index and structural refs resolve against
    /// it. The ids are anonymous, so registration cannot collide.
    pub fn from_ast(ast: &MoleculeAst) -> Self {
        let free = "anonymous entity registration never collides";
        let mut ns = Self::default();
        for _ in ast.atoms().ids() {
            ns.register_atom(None).expect(free);
        }
        for view in ast.bonds().iter() {
            let [a, b] = view.atom_ids();
            ns.register_bond(None, a, b).expect(free);
        }
        for view in ast.dative_bonds().iter() {
            let donors: Vec<AtomId> = view.donor_ids().collect();
            ns.register_dative_bond(None, &donors, view.acceptor_id())
                .expect(free);
        }
        for view in ast.aromatic_systems().iter() {
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            ns.register_aromatic_system(None, &atoms).expect(free);
        }
        for view in ast.multicenter_bonds().iter() {
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            ns.register_multicenter_bond(None, &atoms).expect(free);
        }
        for view in ast.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            ns.register_noncovalent_bond(None, a, b).expect(free);
        }
        for view in ast.stereo_atoms().iter() {
            let ligands: Vec<StereoLigand> = view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect();
            ns.register_stereo_atom(None, view.site_id(), &ligands)
                .expect(free);
        }
        for view in ast.stereo_bonds().iter() {
            let ligands: Vec<StereoLigand> = view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect();
            ns.register_stereo_bond(None, view.site_id(), &ligands)
                .expect(free);
        }
        ns
    }

    /// Whether a keyword is free across the whole namespace (every entity kind + aliases) — the
    /// disjointness check every `register_*` runs before handing out an id.
    fn check_keyword_free(&self, keyword: Option<&str>) -> Result<(), ParseError> {
        match keyword {
            Some(kw) if self.contains_id(kw) => Err(ParseError::DuplicateId(kw.to_string())),
            _ => Ok(()),
        }
    }

    pub(crate) fn register_atom(&mut self, keyword: Option<String>) -> Result<AtomId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self.atoms.register(keyword))
    }

    pub(crate) fn register_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<BondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self.bonds.register(keyword, atom_pair_key(a, b)))
    }

    pub(crate) fn register_dative_bond(
        &mut self,
        keyword: Option<String>,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Result<DativeBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .dative_bonds
            .register(keyword, (donors.iter().copied().collect(), acceptor)))
    }

    pub(crate) fn register_aromatic_system(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<AromaticSystemId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .aromatic_systems
            .register(keyword, atoms.iter().copied().collect()))
    }

    pub(crate) fn register_multicenter_bond(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<MulticenterBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .multicenter_bonds
            .register(keyword, atoms.iter().copied().collect()))
    }

    pub(crate) fn register_noncovalent_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<NoncovalentBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .noncovalent_bonds
            .register(keyword, atom_pair_key(a, b)))
    }

    pub(crate) fn register_stereo_atom(
        &mut self,
        keyword: Option<String>,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Result<StereoAtomId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .stereo_atoms
            .register(keyword, (site, ligand_multiset(ligands))))
    }

    pub(crate) fn register_stereo_bond(
        &mut self,
        keyword: Option<String>,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Result<StereoBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        Ok(self
            .stereo_bonds
            .register(keyword, (site, ligand_multiset(ligands))))
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
        self.atoms.find_by_keyword(keyword)
    }

    pub(crate) fn find_bond_by_keyword(&self, keyword: &str) -> Option<BondId> {
        self.bonds.find_by_keyword(keyword)
    }

    pub(crate) fn find_dative_bond_by_keyword(&self, keyword: &str) -> Option<DativeBondId> {
        self.dative_bonds.find_by_keyword(keyword)
    }

    pub(crate) fn find_aromatic_system_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<AromaticSystemId> {
        self.aromatic_systems.find_by_keyword(keyword)
    }

    pub(crate) fn find_multicenter_bond_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<MulticenterBondId> {
        self.multicenter_bonds.find_by_keyword(keyword)
    }

    pub(crate) fn find_noncovalent_bond_by_keyword(
        &self,
        keyword: &str,
    ) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds.find_by_keyword(keyword)
    }

    pub(crate) fn find_stereo_atom_by_keyword(&self, keyword: &str) -> Option<StereoAtomId> {
        self.stereo_atoms.find_by_keyword(keyword)
    }

    pub(crate) fn find_stereo_bond_by_keyword(&self, keyword: &str) -> Option<StereoBondId> {
        self.stereo_bonds.find_by_keyword(keyword)
    }

    pub(crate) fn find_bond_by_participants(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.bonds.find_by_participants(&atom_pair_key(a, b))
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
        a: AtomId,
        b: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds
            .find_by_participants(&atom_pair_key(a, b))
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
        dsl: Box<AtomDsl>,
    ) -> Result<(), ParseError> {
        if self.contains_id(&name) {
            return Err(ParseError::DuplicateId(name));
        }
        if self.atom_aliases.contains_right(&dsl) {
            return Err(ParseError::InvalidValue(
                "atom-aliases must be bijective: two names map to the same atom".into(),
            ));
        }
        self.atom_aliases.insert(name, dsl);
        Ok(())
    }

    pub(crate) fn atom_keywords(&self) -> impl Iterator<Item = (AtomId, &str)> {
        self.atoms.keywords()
    }

    pub(crate) fn bond_keywords(&self) -> impl Iterator<Item = (BondId, &str)> {
        self.bonds.keywords()
    }

    pub(crate) fn dative_bond_keywords(&self) -> impl Iterator<Item = (DativeBondId, &str)> {
        self.dative_bonds.keywords()
    }

    pub(crate) fn aromatic_system_keywords(
        &self,
    ) -> impl Iterator<Item = (AromaticSystemId, &str)> {
        self.aromatic_systems.keywords()
    }

    pub(crate) fn multicenter_bond_keywords(
        &self,
    ) -> impl Iterator<Item = (MulticenterBondId, &str)> {
        self.multicenter_bonds.keywords()
    }

    pub(crate) fn noncovalent_bond_keywords(
        &self,
    ) -> impl Iterator<Item = (NoncovalentBondId, &str)> {
        self.noncovalent_bonds.keywords()
    }

    pub(crate) fn stereo_atom_keywords(&self) -> impl Iterator<Item = (StereoAtomId, &str)> {
        self.stereo_atoms.keywords()
    }

    pub(crate) fn stereo_bond_keywords(&self) -> impl Iterator<Item = (StereoBondId, &str)> {
        self.stereo_bonds.keywords()
    }

    pub(crate) fn atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(name, dsl)| (name.as_str(), dsl.as_ref()))
    }
}

impl Namespace for MoleculeNamespace {
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

    fn find_bond_by_participants(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.find_bond_by_participants(a, b)
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
        a: AtomId,
        b: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.find_noncovalent_bond_by_participants(a, b)
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

    fn contains_id(&self, id: &str) -> bool {
        self.find_atom_by_keyword(id).is_some()
            || self.find_bond_by_keyword(id).is_some()
            || self.find_dative_bond_by_keyword(id).is_some()
            || self.find_aromatic_system_by_keyword(id).is_some()
            || self.find_multicenter_bond_by_keyword(id).is_some()
            || self.find_noncovalent_bond_by_keyword(id).is_some()
            || self.find_stereo_atom_by_keyword(id).is_some()
            || self.find_stereo_bond_by_keyword(id).is_some()
            || self.atom_aliases.get_by_left(id).is_some()
    }

    fn find_atom_alias(&self, name: &str) -> Option<&AtomDsl> {
        self.atom_aliases.get_by_left(name).map(|dsl| dsl.as_ref())
    }
}

/// Count + id-keyword lookup for one entity kind. Atoms — the base kind, which have no participants —
/// use it directly.
#[derive(Debug)]
struct KeywordRegistry<Id> {
    count: usize,
    by_keyword: IndexMap<String, Id>,
}

impl<Id> Default for KeywordRegistry<Id> {
    fn default() -> Self {
        Self {
            count: 0,
            by_keyword: IndexMap::new(),
        }
    }
}

impl<Id: Copy + From<usize>> KeywordRegistry<Id> {
    /// Reserve the next id (growing the count) and, if the entry carries an `:id`, record its keyword.
    fn register(&mut self, keyword: Option<String>) -> Id {
        let id = Id::from(self.count);
        self.count += 1;
        if let Some(keyword) = keyword {
            self.by_keyword.insert(keyword, id);
        }
        id
    }

    fn find_by_keyword(&self, keyword: &str) -> Option<Id> {
        self.by_keyword.get(keyword).copied()
    }

    /// The keyworded entities of this kind as `(id, keyword)` pairs — the inverse of
    /// `find_by_keyword`, the projection `MoleculeMetadata` needs for rendering.
    fn keywords(&self) -> impl Iterator<Item = (Id, &str)> {
        self.by_keyword
            .iter()
            .map(|(keyword, &id)| (id, keyword.as_str()))
    }

    fn count(&self) -> usize {
        self.count
    }

    /// A registry whose count starts at `count` (so `register` hands out ids from there on) with an
    /// empty keyword map — the shape a delta namespace takes over its lhs's id space.
    fn from_count(count: usize) -> Self {
        Self {
            count,
            by_keyword: IndexMap::new(),
        }
    }
}

/// Count + id-keyword + participant lookup for one non-atom entity kind. `Key` is the entity's
/// canonical participant key (a normalized endpoint pair, atom set, donor-set + acceptor, or stereo
/// site + ligand multiset); §4.1 uniqueness makes it injective, so a hit is unique.
#[derive(Debug)]
struct EntityRegistry<Id, Key> {
    count: usize,
    by_keyword: IndexMap<String, Id>,
    by_participants: HashMap<Key, Id>,
}

impl<Id, Key> Default for EntityRegistry<Id, Key> {
    fn default() -> Self {
        Self {
            count: 0,
            by_keyword: IndexMap::new(),
            by_participants: HashMap::new(),
        }
    }
}

impl<Id: Copy + From<usize>, Key: Eq + Hash> EntityRegistry<Id, Key> {
    /// Reserve the next id, record its `:id` keyword (if any) and its canonical participant key.
    fn register(&mut self, keyword: Option<String>, key: Key) -> Id {
        let id = Id::from(self.count);
        self.count += 1;
        if let Some(keyword) = keyword {
            self.by_keyword.insert(keyword, id);
        }
        self.by_participants.insert(key, id);
        id
    }

    fn find_by_keyword(&self, keyword: &str) -> Option<Id> {
        self.by_keyword.get(keyword).copied()
    }

    fn find_by_participants(&self, key: &Key) -> Option<Id> {
        self.by_participants.get(key).copied()
    }

    fn keywords(&self) -> impl Iterator<Item = (Id, &str)> {
        self.by_keyword
            .iter()
            .map(|(keyword, &id)| (id, keyword.as_str()))
    }

    fn count(&self) -> usize {
        self.count
    }

    fn from_count(count: usize) -> Self {
        Self {
            count,
            by_keyword: IndexMap::new(),
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::ast::ligand::StereoLigandKind;

    #[rstest]
    fn test_molecule_namespace_continuation() {
        let mut lhs = MoleculeNamespace::default();
        lhs.register_atom(None).unwrap();
        lhs.register_atom(Some("c1".into())).unwrap();
        lhs.register_bond(None, AtomId(0), AtomId(1)).unwrap();

        let mut delta = MoleculeNamespace::continuation(&lhs);
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
    fn test_molecule_namespace_register_atom() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(namespace.register_atom(None).unwrap(), AtomId(0));
        assert_eq!(
            namespace.register_atom(Some("c1".into())).unwrap(),
            AtomId(1)
        );
        assert_eq!(namespace.register_atom(None).unwrap(), AtomId(2));
        assert_eq!(namespace.atom_count(), 3);
        assert_eq!(namespace.find_atom_by_keyword("c1"), Some(AtomId(1)));
        assert_eq!(namespace.find_atom_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_atom_error() {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_atom(Some("a".into())).unwrap();
        assert_eq!(
            namespace.register_atom(Some("a".into())).unwrap_err(),
            ParseError::DuplicateId("a".into())
        );
        // The disjointness check spans the whole namespace: a bond may not reuse the atom's keyword.
        assert_eq!(
            namespace
                .register_bond(Some("a".into()), AtomId(0), AtomId(0))
                .unwrap_err(),
            ParseError::DuplicateId("a".into())
        );
    }

    #[rstest]
    fn test_molecule_namespace_register_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace.register_bond(None, AtomId(2), AtomId(0)).unwrap(),
            BondId(0)
        );
        assert_eq!(
            namespace
                .register_bond(Some("b1".into()), AtomId(1), AtomId(3))
                .unwrap(),
            BondId(1)
        );
        assert_eq!(namespace.bond_count(), 2);
        assert_eq!(namespace.find_bond_by_keyword("b1"), Some(BondId(1)));
        assert_eq!(namespace.find_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_dative_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace
                .register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0))
                .unwrap(),
            DativeBondId(0)
        );
        assert_eq!(
            namespace
                .register_dative_bond(Some("d1".into()), &[AtomId(4)], AtomId(3))
                .unwrap(),
            DativeBondId(1)
        );
        assert_eq!(namespace.dative_bond_count(), 2);
        assert_eq!(
            namespace.find_dative_bond_by_keyword("d1"),
            Some(DativeBondId(1))
        );
        assert_eq!(namespace.find_dative_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_aromatic_system() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace
                .register_aromatic_system(None, &[AtomId(0), AtomId(1), AtomId(2)])
                .unwrap(),
            AromaticSystemId(0)
        );
        assert_eq!(
            namespace
                .register_aromatic_system(Some("a1".into()), &[AtomId(3), AtomId(4)])
                .unwrap(),
            AromaticSystemId(1)
        );
        assert_eq!(namespace.aromatic_system_count(), 2);
        assert_eq!(
            namespace.find_aromatic_system_by_keyword("a1"),
            Some(AromaticSystemId(1))
        );
        assert_eq!(namespace.find_aromatic_system_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_multicenter_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace
                .register_multicenter_bond(Some("m".into()), &[AtomId(0), AtomId(1), AtomId(2)])
                .unwrap(),
            MulticenterBondId(0)
        );
        assert_eq!(
            namespace
                .register_multicenter_bond(None, &[AtomId(3), AtomId(4)])
                .unwrap(),
            MulticenterBondId(1)
        );
        assert_eq!(namespace.multicenter_bond_count(), 2);
        assert_eq!(
            namespace.find_multicenter_bond_by_keyword("m"),
            Some(MulticenterBondId(0))
        );
        assert_eq!(namespace.find_multicenter_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_noncovalent_bond() {
        let mut namespace = MoleculeNamespace::default();
        assert_eq!(
            namespace
                .register_noncovalent_bond(None, AtomId(3), AtomId(1))
                .unwrap(),
            NoncovalentBondId(0)
        );
        assert_eq!(
            namespace
                .register_noncovalent_bond(Some("n1".into()), AtomId(0), AtomId(4))
                .unwrap(),
            NoncovalentBondId(1)
        );
        assert_eq!(namespace.noncovalent_bond_count(), 2);
        assert_eq!(
            namespace.find_noncovalent_bond_by_keyword("n1"),
            Some(NoncovalentBondId(1))
        );
        assert_eq!(namespace.find_noncovalent_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_stereo_atom() {
        let mut namespace = MoleculeNamespace::default();
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        assert_eq!(
            namespace
                .register_stereo_atom(None, AtomId(4), &ligands)
                .unwrap(),
            StereoAtomId(0)
        );
        assert_eq!(
            namespace
                .register_stereo_atom(Some("s1".into()), AtomId(5), &ligands)
                .unwrap(),
            StereoAtomId(1)
        );
        assert_eq!(namespace.stereo_atom_count(), 2);
        assert_eq!(
            namespace.find_stereo_atom_by_keyword("s1"),
            Some(StereoAtomId(1))
        );
        assert_eq!(namespace.find_stereo_atom_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_stereo_bond() {
        let mut namespace = MoleculeNamespace::default();
        let ligands = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        assert_eq!(
            namespace
                .register_stereo_bond(None, BondId(2), &ligands)
                .unwrap(),
            StereoBondId(0)
        );
        assert_eq!(
            namespace
                .register_stereo_bond(Some("sb1".into()), BondId(0), &ligands)
                .unwrap(),
            StereoBondId(1)
        );
        assert_eq!(namespace.stereo_bond_count(), 2);
        assert_eq!(
            namespace.find_stereo_bond_by_keyword("sb1"),
            Some(StereoBondId(1))
        );
        assert_eq!(namespace.find_stereo_bond_by_keyword("nope"), None);
    }

    #[rstest]
    fn test_molecule_namespace_register_atom_alias() {
        let mut namespace = MoleculeNamespace::default();
        let dsl = Box::new("C".parse::<AtomDsl>().unwrap());
        namespace
            .register_atom_alias("me".into(), dsl.clone())
            .unwrap();
        assert_eq!(namespace.find_atom_alias("me"), Some(dsl.as_ref()));
        assert!(namespace.contains_id("me"));
        assert_eq!(namespace.find_atom_alias("nope"), None);
    }

    #[rstest]
    #[case::name_taken_by_entity("a", "N", ParseError::DuplicateId("a".into()))]
    #[case::name_taken_by_alias("me", "O", ParseError::DuplicateId("me".into()))]
    #[case::duplicate_spec(
        "me2",
        "C",
        ParseError::InvalidValue(
            "atom-aliases must be bijective: two names map to the same atom".into()
        )
    )]
    fn test_molecule_namespace_register_atom_alias_error(
        #[case] name: &str,
        #[case] spec: &str,
        #[case] expected: ParseError,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_atom(Some("a".into())).unwrap();
        namespace
            .register_atom_alias("me".into(), Box::new("C".parse::<AtomDsl>().unwrap()))
            .unwrap();
        assert_eq!(
            namespace
                .register_atom_alias(name.into(), Box::new(spec.parse::<AtomDsl>().unwrap()))
                .unwrap_err(),
            expected
        );
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
    fn test_molecule_namespace_contains_id(#[case] id: &str, #[case] expected: bool) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_atom(Some("a".into())).unwrap();
        namespace
            .register_bond(Some("b".into()), AtomId(0), AtomId(0))
            .unwrap();
        namespace
            .register_dative_bond(Some("d".into()), &[AtomId(0)], AtomId(0))
            .unwrap();
        namespace
            .register_aromatic_system(Some("ar".into()), &[AtomId(0)])
            .unwrap();
        namespace
            .register_multicenter_bond(Some("m".into()), &[AtomId(0)])
            .unwrap();
        namespace
            .register_noncovalent_bond(Some("nc".into()), AtomId(0), AtomId(0))
            .unwrap();
        namespace
            .register_stereo_atom(Some("sa".into()), AtomId(0), &[])
            .unwrap();
        namespace
            .register_stereo_bond(Some("sb".into()), BondId(0), &[])
            .unwrap();
        namespace
            .register_atom_alias("al".into(), Box::new("C".parse::<AtomDsl>().unwrap()))
            .unwrap();
        assert_eq!(namespace.contains_id(id), expected);
    }

    #[rstest]
    #[case::forward(AtomId(0), AtomId(2), Some(BondId(0)))]
    #[case::reversed(AtomId(2), AtomId(0), Some(BondId(0)))]
    #[case::absent(AtomId(0), AtomId(4), None)]
    fn test_molecule_namespace_find_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<BondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_bond(None, AtomId(2), AtomId(0)).unwrap();
        assert_eq!(namespace.find_bond_by_participants(a, b), expected);
    }

    #[rstest]
    #[case::donors_reordered(&[AtomId(2), AtomId(1)], AtomId(0), Some(DativeBondId(0)))]
    #[case::wrong_acceptor(&[AtomId(1), AtomId(2)], AtomId(3), None)]
    #[case::wrong_donors(&[AtomId(1), AtomId(3)], AtomId(0), None)]
    fn test_molecule_namespace_find_dative_bond_by_participants(
        #[case] donors: &[AtomId],
        #[case] acceptor: AtomId,
        #[case] expected: Option<DativeBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace
            .register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0))
            .unwrap();
        assert_eq!(
            namespace.find_dative_bond_by_participants(donors, acceptor),
            expected
        );
    }

    #[rstest]
    #[case::reordered(&[AtomId(0), AtomId(1), AtomId(2)], Some(AromaticSystemId(0)))]
    #[case::subset(&[AtomId(0), AtomId(1)], None)]
    #[case::superset(&[AtomId(0), AtomId(1), AtomId(2), AtomId(3)], None)]
    fn test_molecule_namespace_find_aromatic_system_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace
            .register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)])
            .unwrap();
        assert_eq!(
            namespace.find_aromatic_system_by_participants(atoms),
            expected
        );
    }

    #[rstest]
    #[case::reordered(&[AtomId(2), AtomId(1), AtomId(0)], Some(MulticenterBondId(0)))]
    #[case::absent(&[AtomId(0), AtomId(1), AtomId(3)], None)]
    fn test_molecule_namespace_find_multicenter_bond_by_participants(
        #[case] atoms: &[AtomId],
        #[case] expected: Option<MulticenterBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace
            .register_multicenter_bond(None, &[AtomId(0), AtomId(1), AtomId(2)])
            .unwrap();
        assert_eq!(
            namespace.find_multicenter_bond_by_participants(atoms),
            expected
        );
    }

    #[rstest]
    #[case::reversed(AtomId(1), AtomId(3), Some(NoncovalentBondId(0)))]
    #[case::absent(AtomId(1), AtomId(2), None)]
    fn test_molecule_namespace_find_noncovalent_bond_by_participants(
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: Option<NoncovalentBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        namespace
            .register_noncovalent_bond(None, AtomId(3), AtomId(1))
            .unwrap();
        assert_eq!(
            namespace.find_noncovalent_bond_by_participants(a, b),
            expected
        );
    }

    #[rstest]
    #[case::reordered_ligands(AtomId(4), &[AtomId(2), AtomId(1)], Some(StereoAtomId(0)))]
    #[case::wrong_ligands(AtomId(4), &[AtomId(1)], None)]
    #[case::wrong_site(AtomId(0), &[AtomId(1), AtomId(2)], None)]
    fn test_molecule_namespace_find_stereo_atom_by_participants(
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
        namespace
            .register_stereo_atom(None, AtomId(4), &registered)
            .unwrap();
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            namespace.find_stereo_atom_by_participants(site, &query),
            expected
        );
    }

    #[rstest]
    #[case::matching(BondId(2), &[AtomId(3)], Some(StereoBondId(0)))]
    #[case::wrong_site(BondId(0), &[AtomId(3)], None)]
    #[case::empty_ligands(BondId(2), &[], None)]
    fn test_molecule_namespace_find_stereo_bond_by_participants(
        #[case] site: BondId,
        #[case] ligand_atoms: &[AtomId],
        #[case] expected: Option<StereoBondId>,
    ) {
        let mut namespace = MoleculeNamespace::default();
        let registered = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        namespace
            .register_stereo_bond(None, BondId(2), &registered)
            .unwrap();
        let query: Vec<StereoLigand> = ligand_atoms
            .iter()
            .map(|&a| StereoLigand::new(a, StereoLigandKind::Atom))
            .collect();
        assert_eq!(
            namespace.find_stereo_bond_by_participants(site, &query),
            expected
        );
    }
}
