//! Surface-form metadata for molecule and reaction DSLs.

use bimap::BiBTreeMap;
use indexmap::IndexMap;
use thiserror::Error;

use super::atom::AtomDsl;
use super::namespace::MoleculeNamespace;
use super::reaction::ReactionNamespace;
use crate::ast::correspondence::MoleculeCorrespondence;
use crate::ast::entity::Entity;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Error raised when metadata would no longer define disjoint entity keywords
/// and bijective atom aliases.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MetadataError {
    #[error("duplicate keyword: {0}")]
    DuplicateKeyword(String),
    #[error("atom DSL already has alias: {0}")]
    DuplicateAtomAlias(String),
}

/// The rendering counterpart to [`crate::dsl::Namespace`].
///
/// Entity-keyword lookup renders an AST id as a keyword reference when one is
/// available, or as its positional index otherwise. Rendering never emits
/// structural references, so this surface needs neither counts nor participant
/// indexes.
pub trait Metadata {
    fn keyword(&self, entity: Entity) -> Option<&str>;
    fn entity(&self, keyword: &str) -> Option<Entity>;
}

/// Surface-form metadata paired with a `MoleculeAst`. Records entity keywords and atom aliases.
/// `MoleculeDsl` keeps both fields private and rewraps atomically
/// through `from_parts`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeMetadata {
    keywords: BiBTreeMap<Entity, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

impl MoleculeMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether this metadata binds no entity keywords and no atom aliases — the shape an anonymous
    /// molecule (e.g. a sub-pattern) projects.
    pub fn is_empty(&self) -> bool {
        self.keywords.is_empty() && self.atom_aliases.is_empty()
    }

    pub fn keyword(&self, entity: Entity) -> Option<&str> {
        self.keywords.get_by_left(&entity).map(String::as_str)
    }

    pub fn entity(&self, keyword: &str) -> Option<Entity> {
        self.keywords.get_by_right(keyword).copied()
    }

    /// Whether `name` is already bound as an entity keyword.
    pub fn contains_keyword(&self, name: &str) -> bool {
        self.keywords.contains_right(name)
    }

    /// Name of the alias bound to this atom DSL, if any.
    pub fn atom_alias_for(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn has_atom_alias(&self, name: &str) -> bool {
        self.atom_aliases.contains_left(name)
    }

    pub fn has_atom_aliases(&self) -> bool {
        !self.atom_aliases.is_empty()
    }

    pub fn atom_aliases_len(&self) -> usize {
        self.atom_aliases.len()
    }

    pub fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_keyword(
        &mut self,
        entity: Entity,
        name: impl Into<String>,
    ) -> Result<(), MetadataError> {
        let name = name.into();

        if self.keywords.get_by_left(&entity) == Some(&name) {
            return Ok(());
        }
        if self.keywords.contains_right(name.as_str())
            || self.atom_aliases.contains_left(name.as_str())
        {
            return Err(MetadataError::DuplicateKeyword(name));
        }

        self.keywords.insert(entity, name);
        Ok(())
    }

    pub fn add_atom_alias(
        &mut self,
        name: impl Into<String>,
        atom: impl Into<AtomDsl>,
    ) -> Result<(), MetadataError> {
        let name = name.into();
        let atom = Box::new(atom.into());

        if self.keywords.contains_right(name.as_str()) {
            return Err(MetadataError::DuplicateKeyword(name));
        }
        if let Some(existing) = self.atom_aliases.get_by_left(name.as_str()) {
            return if existing == &atom {
                Ok(())
            } else {
                Err(MetadataError::DuplicateKeyword(name))
            };
        }
        if let Some(existing) = self.atom_aliases.get_by_right(atom.as_ref()) {
            return Err(MetadataError::DuplicateAtomAlias(existing.clone()));
        }

        self.atom_aliases.insert(name, atom);
        Ok(())
    }

    /// Move entity keywords from the left id space of `correspondence` to its
    /// matched right entities. Keywords on unmatched left entities are omitted;
    /// atom aliases are independent of molecule ids and remain unchanged.
    pub fn remap(self, correspondence: &MoleculeCorrespondence) -> Self {
        let Self {
            keywords,
            atom_aliases,
        } = self;
        let keywords = keywords
            .into_iter()
            .filter_map(|(entity, keyword)| {
                correspondence
                    .right_of(entity)
                    .map(|right| (right, keyword))
            })
            .collect();

        Self {
            keywords,
            atom_aliases,
        }
    }
}

impl Metadata for MoleculeMetadata {
    fn keyword(&self, entity: Entity) -> Option<&str> {
        self.keyword(entity)
    }

    fn entity(&self, keyword: &str) -> Option<Entity> {
        self.entity(keyword)
    }
}

impl From<&MoleculeNamespace> for MoleculeMetadata {
    /// Project the namespace to its roundtrip subset: the eight `id → keyword` maps (the inverse of
    /// the namespace's `find_by_keyword`) plus the atom aliases. The namespace is the source of truth;
    /// this is the derived view — parse-only data (participant indexes, counts) is dropped.
    fn from(namespace: &MoleculeNamespace) -> Self {
        let mut metadata = MoleculeMetadata::new();
        for (id, keyword) in namespace.atom_keywords() {
            metadata
                .set_keyword(Entity::Atom(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.bond_keywords() {
            metadata
                .set_keyword(Entity::Bond(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.dative_bond_keywords() {
            metadata
                .set_keyword(Entity::DativeBond(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.aromatic_system_keywords() {
            metadata
                .set_keyword(Entity::AromaticSystem(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.multicenter_bond_keywords() {
            metadata
                .set_keyword(Entity::MulticenterBond(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.noncovalent_bond_keywords() {
            metadata
                .set_keyword(Entity::NoncovalentBond(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.stereo_atom_keywords() {
            metadata
                .set_keyword(Entity::StereoAtom(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (id, keyword) in namespace.stereo_bond_keywords() {
            metadata
                .set_keyword(Entity::StereoBond(id), keyword)
                .expect("namespace keywords are disjoint");
        }
        for (name, dsl) in namespace.atom_aliases() {
            metadata
                .add_atom_alias(name, dsl.clone())
                .expect("namespace aliases are bijective and disjoint from keywords");
        }
        metadata
    }
}

/// Surface-form metadata paired with a `ReactionAst`: the lhs molecule metadata plus the
/// created-entity keyword bindings and atom aliases introduced by the deltas. Mirrors
/// `MoleculeMetadata` for the atom/bond entities (the reaction admits the `[:C "C#h3"]`
/// alias notation for added atoms).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    atom_ids: IndexMap<AtomId, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondId, String>,
    dative_bond_ids: IndexMap<DativeBondId, String>,
    aromatic_system_ids: IndexMap<AromaticSystemId, String>,
    multicenter_bond_ids: IndexMap<MulticenterBondId, String>,
    noncovalent_bond_ids: IndexMap<NoncovalentBondId, String>,
    stereo_atom_ids: IndexMap<StereoAtomId, String>,
    stereo_bond_ids: IndexMap<StereoBondId, String>,
}

impl ReactionMetadata {
    pub fn lhs(&self) -> &MoleculeMetadata {
        &self.lhs
    }

    pub fn combined_metadata(&self) -> MoleculeMetadata {
        let mut combined = self.lhs.clone();
        for (&id, name) in &self.atom_ids {
            combined
                .set_keyword(Entity::Atom(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.bond_ids {
            combined
                .set_keyword(Entity::Bond(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.dative_bond_ids {
            combined
                .set_keyword(Entity::DativeBond(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.aromatic_system_ids {
            combined
                .set_keyword(Entity::AromaticSystem(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.multicenter_bond_ids {
            combined
                .set_keyword(Entity::MulticenterBond(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.noncovalent_bond_ids {
            combined
                .set_keyword(Entity::NoncovalentBond(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.stereo_atom_ids {
            combined
                .set_keyword(Entity::StereoAtom(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        for (&id, name) in &self.stereo_bond_ids {
            combined
                .set_keyword(Entity::StereoBond(id), name)
                .expect("reaction metadata keywords are disjoint");
        }
        combined
    }

    pub fn atom_keyword(&self, id: AtomId) -> Option<&str> {
        self.atom_ids.get(&id).map(String::as_str)
    }

    pub fn bond_keyword(&self, id: BondId) -> Option<&str> {
        self.bond_ids.get(&id).map(String::as_str)
    }

    /// Name of the alias bound to this atom DSL, if any.
    pub fn atom_alias_for(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn has_atom_alias(&self, name: &str) -> bool {
        self.atom_aliases.contains_left(name)
    }

    pub fn has_atom_aliases(&self) -> bool {
        !self.atom_aliases.is_empty()
    }

    pub fn atom_aliases_len(&self) -> usize {
        self.atom_aliases.len()
    }

    pub fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_atom_keyword(&mut self, id: AtomId, name: impl Into<String>) {
        self.atom_ids.insert(id, name.into());
    }

    pub fn set_bond_keyword(&mut self, id: BondId, name: impl Into<String>) {
        self.bond_ids.insert(id, name.into());
    }

    pub fn dative_bond_keyword(&self, id: DativeBondId) -> Option<&str> {
        self.dative_bond_ids.get(&id).map(String::as_str)
    }

    pub fn set_dative_bond_keyword(&mut self, id: DativeBondId, name: impl Into<String>) {
        self.dative_bond_ids.insert(id, name.into());
    }

    pub fn aromatic_system_keyword(&self, id: AromaticSystemId) -> Option<&str> {
        self.aromatic_system_ids.get(&id).map(String::as_str)
    }

    pub fn set_aromatic_system_keyword(&mut self, id: AromaticSystemId, name: impl Into<String>) {
        self.aromatic_system_ids.insert(id, name.into());
    }

    pub fn multicenter_bond_keyword(&self, id: MulticenterBondId) -> Option<&str> {
        self.multicenter_bond_ids.get(&id).map(String::as_str)
    }

    pub fn set_multicenter_bond_keyword(&mut self, id: MulticenterBondId, name: impl Into<String>) {
        self.multicenter_bond_ids.insert(id, name.into());
    }

    pub fn noncovalent_bond_keyword(&self, id: NoncovalentBondId) -> Option<&str> {
        self.noncovalent_bond_ids.get(&id).map(String::as_str)
    }

    pub fn set_noncovalent_bond_keyword(&mut self, id: NoncovalentBondId, name: impl Into<String>) {
        self.noncovalent_bond_ids.insert(id, name.into());
    }

    pub fn stereo_atom_keyword(&self, id: StereoAtomId) -> Option<&str> {
        self.stereo_atom_ids.get(&id).map(String::as_str)
    }

    pub fn set_stereo_atom_keyword(&mut self, id: StereoAtomId, name: impl Into<String>) {
        self.stereo_atom_ids.insert(id, name.into());
    }

    pub fn stereo_bond_keyword(&self, id: StereoBondId) -> Option<&str> {
        self.stereo_bond_ids.get(&id).map(String::as_str)
    }

    pub fn set_stereo_bond_keyword(&mut self, id: StereoBondId, name: impl Into<String>) {
        self.stereo_bond_ids.insert(id, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases.insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_keyword(mut self, id: AtomId, name: impl Into<String>) -> Self {
        self.set_atom_keyword(id, name);
        self
    }

    pub fn with_bond_keyword(mut self, id: BondId, name: impl Into<String>) -> Self {
        self.set_bond_keyword(id, name);
        self
    }

    pub fn with_atom_alias(mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

impl From<MoleculeMetadata> for ReactionMetadata {
    fn from(lhs: MoleculeMetadata) -> Self {
        Self {
            lhs,
            ..Self::default()
        }
    }
}

impl From<&ReactionNamespace> for ReactionMetadata {
    /// Project the roundtrip metadata: the lhs molecule's metadata, the delta-introduced entity
    /// keywords (any delta that binds a name, not only `:add`), and the reaction's top-level aliases.
    fn from(ns: &ReactionNamespace) -> Self {
        let mut metadata = ReactionMetadata {
            lhs: MoleculeMetadata::from(ns.lhs()),
            ..Default::default()
        };
        for (id, name) in ns.deltas().atom_keywords() {
            metadata.set_atom_keyword(id, name);
        }
        for (id, name) in ns.deltas().bond_keywords() {
            metadata.set_bond_keyword(id, name);
        }
        for (id, name) in ns.deltas().dative_bond_keywords() {
            metadata.set_dative_bond_keyword(id, name);
        }
        for (id, name) in ns.deltas().aromatic_system_keywords() {
            metadata.set_aromatic_system_keyword(id, name);
        }
        for (id, name) in ns.deltas().multicenter_bond_keywords() {
            metadata.set_multicenter_bond_keyword(id, name);
        }
        for (id, name) in ns.deltas().noncovalent_bond_keywords() {
            metadata.set_noncovalent_bond_keyword(id, name);
        }
        for (id, name) in ns.deltas().stereo_atom_keywords() {
            metadata.set_stereo_atom_keyword(id, name);
        }
        for (id, name) in ns.deltas().stereo_bond_keywords() {
            metadata.set_stereo_bond_keyword(id, name);
        }
        for (name, dsl) in ns.atom_aliases() {
            metadata.add_atom_alias(name.to_string(), dsl.clone());
        }
        metadata
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{Correspondence, NodeId};

    use super::*;
    use crate::ast::atom::AtomAst;

    #[rstest]
    fn test_molecule_metadata_new() {
        let actual = MoleculeMetadata::new();

        assert_eq!(actual, MoleculeMetadata::default());
        assert!(actual.is_empty());
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(1)))]
    #[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
    #[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
    #[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
    #[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
    fn test_molecule_metadata_keyword(#[case] entity: Entity) {
        let metadata = MoleculeMetadata {
            keywords: [(entity, "key".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.keyword(entity), Some("key"));
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(1)))]
    #[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
    #[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
    #[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
    #[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
    fn test_molecule_metadata_entity(#[case] entity: Entity) {
        let metadata = MoleculeMetadata {
            keywords: [(entity, "key".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.entity("key"), Some(entity));
        assert_eq!(metadata.entity("missing"), None);
    }

    #[rstest]
    #[case::present("key", true)]
    #[case::absent("other", false)]
    fn test_molecule_metadata_contains_keyword(#[case] keyword: &str, #[case] expected: bool) {
        let metadata = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "key".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.contains_keyword(keyword), expected);
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(1)))]
    #[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
    #[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
    #[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
    #[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
    fn test_molecule_metadata_set_keyword(#[case] entity: Entity) {
        let mut actual = MoleculeMetadata::new();
        let result = actual.set_keyword(entity, "key");
        let expected = MoleculeMetadata {
            keywords: [(entity, "key".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_molecule_metadata_set_keyword_idempotent() {
        let mut actual = MoleculeMetadata::new();
        actual.set_keyword(Entity::Atom(AtomId(0)), "key").unwrap();
        let expected = actual.clone();

        let result = actual.set_keyword(Entity::Atom(AtomId(0)), "key");

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_molecule_metadata_set_keyword_rebinding() {
        let mut actual = MoleculeMetadata::new();
        actual.set_keyword(Entity::Atom(AtomId(0)), "old").unwrap();

        let result = actual.set_keyword(Entity::Atom(AtomId(0)), "new");

        assert_eq!(result, Ok(()));
        assert_eq!(actual.keyword(Entity::Atom(AtomId(0))), Some("new"));
        assert_eq!(actual.entity("old"), None);
    }

    #[rstest]
    fn test_molecule_metadata_set_keyword_error() {
        let mut actual = MoleculeMetadata::new();
        actual.set_keyword(Entity::Atom(AtomId(0)), "key").unwrap();
        let expected = actual.clone();

        let result = actual.set_keyword(Entity::Bond(BondId(0)), "key");

        assert_eq!(
            result,
            Err(MetadataError::DuplicateKeyword("key".to_string()))
        );
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_molecule_metadata_atom_alias_for() {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let metadata = MoleculeMetadata {
            keywords: BiBTreeMap::new(),
            atom_aliases: [("carbon".to_string(), Box::new(atom.clone()))]
                .into_iter()
                .collect(),
        };

        assert_eq!(metadata.atom_alias_for(&atom), Some("carbon"));
    }

    #[rstest]
    fn test_molecule_metadata_add_atom_alias() {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let mut actual = MoleculeMetadata::new();

        let result = actual.add_atom_alias("carbon", atom.clone());

        assert_eq!(result, Ok(()));
        assert_eq!(
            actual,
            MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [("carbon".to_string(), Box::new(atom))]
                    .into_iter()
                    .collect(),
            }
        );
    }

    #[rstest]
    fn test_molecule_metadata_add_atom_alias_idempotent() {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let mut actual = MoleculeMetadata::new();
        actual.add_atom_alias("carbon", atom.clone()).unwrap();
        let expected = actual.clone();

        let result = actual.add_atom_alias("carbon", atom);

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::keyword(
        MetadataError::DuplicateKeyword("carbon".to_string()),
        AtomDsl(AtomAst::from_element(Element::N)),
        "carbon"
    )]
    #[case::atom(
        MetadataError::DuplicateAtomAlias("carbon".to_string()),
        AtomDsl(AtomAst::from_element(Element::C)),
        "other"
    )]
    fn test_molecule_metadata_add_atom_alias_error(
        #[case] expected_error: MetadataError,
        #[case] atom: AtomDsl,
        #[case] alias: &str,
    ) {
        let mut actual = MoleculeMetadata::new();
        actual
            .add_atom_alias("carbon", AtomAst::from_element(Element::C))
            .unwrap();
        let expected = actual.clone();

        let result = actual.add_atom_alias(alias, atom);

        assert_eq!(result, Err(expected_error));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_molecule_metadata_remap() {
        let input = MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(0)), "atom".to_string()),
                (Entity::Bond(BondId(0)), "bond".to_string()),
                (Entity::DativeBond(DativeBondId(0)), "dative".to_string()),
                (
                    Entity::AromaticSystem(AromaticSystemId(0)),
                    "aromatic".to_string(),
                ),
                (
                    Entity::MulticenterBond(MulticenterBondId(0)),
                    "multicenter".to_string(),
                ),
                (
                    Entity::NoncovalentBond(NoncovalentBondId(0)),
                    "noncovalent".to_string(),
                ),
                (
                    Entity::StereoAtom(StereoAtomId(0)),
                    "stereo-atom".to_string(),
                ),
                (
                    Entity::StereoBond(StereoBondId(0)),
                    "stereo-bond".to_string(),
                ),
            ]
            .into_iter()
            .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&[NodeId(1), NodeId(0)], 2),
            Correspondence::from_images(&[BondId(1), BondId(0)], 2),
            Correspondence::from_images(&[DativeBondId(1), DativeBondId(0)], 2),
            Correspondence::from_images(&[AromaticSystemId(1), AromaticSystemId(0)], 2),
            Correspondence::from_images(&[MulticenterBondId(1), MulticenterBondId(0)], 2),
            Correspondence::from_images(&[NoncovalentBondId(1), NoncovalentBondId(0)], 2),
            Correspondence::from_images(&[StereoAtomId(1), StereoAtomId(0)], 2),
            Correspondence::from_images(&[StereoBondId(1), StereoBondId(0)], 2),
        );
        let expected = MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(1)), "atom".to_string()),
                (Entity::Bond(BondId(1)), "bond".to_string()),
                (Entity::DativeBond(DativeBondId(1)), "dative".to_string()),
                (
                    Entity::AromaticSystem(AromaticSystemId(1)),
                    "aromatic".to_string(),
                ),
                (
                    Entity::MulticenterBond(MulticenterBondId(1)),
                    "multicenter".to_string(),
                ),
                (
                    Entity::NoncovalentBond(NoncovalentBondId(1)),
                    "noncovalent".to_string(),
                ),
                (
                    Entity::StereoAtom(StereoAtomId(1)),
                    "stereo-atom".to_string(),
                ),
                (
                    Entity::StereoBond(StereoBondId(1)),
                    "stereo-bond".to_string(),
                ),
            ]
            .into_iter()
            .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(input.remap(&correspondence), expected);
    }

    #[rstest]
    fn test_molecule_metadata_remap_identity() {
        let input = MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(0)), "atom".to_string()),
                (Entity::Bond(BondId(0)), "bond".to_string()),
            ]
            .into_iter()
            .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&[NodeId(0)], 1),
            Correspondence::from_images(&[BondId(0)], 1),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
        );

        assert_eq!(input.clone().remap(&correspondence), input);
    }

    #[rstest]
    fn test_molecule_metadata_remap_partial() {
        let input = MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(0)), "removed-atom".to_string()),
                (Entity::Atom(AtomId(1)), "retained-atom".to_string()),
                (Entity::Bond(BondId(0)), "removed-bond".to_string()),
            ]
            .into_iter()
            .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(0))], 2, 1),
            Correspondence::new(Vec::new(), 1, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
        );
        let expected = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "retained-atom".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(input.remap(&correspondence), expected);
    }

    #[rstest]
    fn test_molecule_metadata_remap_roundtrip() {
        let input = MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(0)), "first".to_string()),
                (Entity::Atom(AtomId(1)), "second".to_string()),
            ]
            .into_iter()
            .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&[NodeId(1), NodeId(0)], 2),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
        );

        assert_eq!(
            input
                .clone()
                .remap(&correspondence)
                .remap(&correspondence.reverse()),
            input
        );
    }

    #[rstest]
    fn test_molecule_metadata_remap_composition() {
        let input = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "atom".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };
        let first = MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(1))], 1, 2),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
        );
        let second = MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(2))], 2, 3),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
            Correspondence::new(Vec::new(), 0, 0),
        );
        let expected = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(2)), "atom".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };
        let sequential = input.clone().remap(&first).remap(&second);
        let composed = input.remap(&first.compose(&second));

        assert_eq!(sequential, composed);
        assert_eq!(composed, expected);
    }

    #[rstest]
    #[case::keyword_then_alias(false)]
    #[case::alias_then_keyword(true)]
    fn test_molecule_metadata_keyword_alias_collision(#[case] alias_first: bool) {
        let atom = AtomAst::from_element(Element::C);
        let mut actual = MoleculeMetadata::new();
        let result = if alias_first {
            actual.add_atom_alias("carbon", atom).unwrap();
            actual.set_keyword(Entity::Atom(AtomId(0)), "carbon")
        } else {
            actual
                .set_keyword(Entity::Atom(AtomId(0)), "carbon")
                .unwrap();
            actual.add_atom_alias("carbon", atom)
        };

        assert_eq!(
            result,
            Err(MetadataError::DuplicateKeyword("carbon".to_string()))
        );
        assert_eq!(actual.contains_keyword("carbon"), !alias_first);
        assert_eq!(actual.has_atom_alias("carbon"), alias_first);
    }

    #[rstest]
    fn test_molecule_metadata_iter_atom_aliases() {
        let metadata = MoleculeMetadata {
            keywords: BiBTreeMap::new(),
            atom_aliases: [
                (
                    "carbon".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::C))),
                ),
                (
                    "nitrogen".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::N))),
                ),
            ]
            .into_iter()
            .collect(),
        };
        let actual: Vec<(&str, &AtomDsl)> = metadata.iter_atom_aliases().collect();

        assert_eq!(
            actual,
            vec![
                ("carbon", &AtomDsl(AtomAst::from_element(Element::C))),
                ("nitrogen", &AtomDsl(AtomAst::from_element(Element::N))),
            ]
        );
    }
}
