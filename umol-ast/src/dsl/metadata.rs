//! Surface-form metadata for molecule and reaction DSLs.

use bimap::BiBTreeMap;
use thiserror::Error;

use super::atom::AtomDsl;
use crate::ast::correspondence::MoleculeCorrespondence;
use crate::ast::entity::Entity;
#[cfg(test)]
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Error raised when metadata would violate its namespace invariants or refer
/// outside the AST with which it is paired.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MetadataError {
    #[error("duplicate keyword: {0}")]
    DuplicateKeyword(String),
    #[error("atom DSL already has alias: {0}")]
    DuplicateAtomAlias(String),
    #[error("metadata entity is out of range: {0}")]
    EntityOutOfRange(Entity),
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
/// `MoleculeDsl` keeps both fields private and validates their coherence during
/// checked construction.
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

    pub fn iter_keywords(&self) -> impl ExactSizeIterator<Item = (Entity, &str)> {
        self.keywords
            .iter()
            .map(|(entity, keyword)| (*entity, keyword.as_str()))
    }

    pub fn atom_alias(&self, name: &str) -> Option<&AtomDsl> {
        self.atom_aliases.get_by_left(name).map(Box::as_ref)
    }

    pub fn atom_alias_name(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn iter_atom_aliases(&self) -> impl ExactSizeIterator<Item = (&str, &AtomDsl)> {
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

/// Surface-form metadata paired with a `ReactionAst`: lhs molecule metadata,
/// entity keywords introduced by deltas, and reaction-scope atom aliases.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    delta_keywords: BiBTreeMap<Entity, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

struct ReactionKeywordIter<D, L> {
    delta: D,
    lhs: L,
}

impl<D, L> ReactionKeywordIter<D, L>
where
    D: ExactSizeIterator,
    L: ExactSizeIterator<Item = D::Item>,
{
    fn remaining_len(&self) -> usize {
        self.delta
            .len()
            .checked_add(self.lhs.len())
            .expect("reaction metadata keyword count exceeds usize")
    }
}

impl<D, L> Iterator for ReactionKeywordIter<D, L>
where
    D: ExactSizeIterator,
    L: ExactSizeIterator<Item = D::Item>,
{
    type Item = D::Item;

    fn next(&mut self) -> Option<Self::Item> {
        self.delta.next().or_else(|| self.lhs.next())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.remaining_len();
        (len, Some(len))
    }
}

impl<D, L> ExactSizeIterator for ReactionKeywordIter<D, L>
where
    D: ExactSizeIterator,
    L: ExactSizeIterator<Item = D::Item>,
{
    fn len(&self) -> usize {
        self.remaining_len()
    }
}

impl ReactionMetadata {
    pub fn lhs(&self) -> &MoleculeMetadata {
        &self.lhs
    }

    pub fn keyword(&self, entity: Entity) -> Option<&str> {
        self.delta_keyword(entity)
            .or_else(|| self.lhs.keyword(entity))
    }

    pub fn entity(&self, keyword: &str) -> Option<Entity> {
        self.delta_entity(keyword)
            .or_else(|| self.lhs.entity(keyword))
    }

    pub fn iter_keywords(&self) -> impl ExactSizeIterator<Item = (Entity, &str)> {
        ReactionKeywordIter {
            delta: self.iter_delta_keywords(),
            lhs: self.lhs.iter_keywords(),
        }
    }

    pub fn delta_keyword(&self, entity: Entity) -> Option<&str> {
        self.delta_keywords.get_by_left(&entity).map(String::as_str)
    }

    pub fn delta_entity(&self, keyword: &str) -> Option<Entity> {
        self.delta_keywords.get_by_right(keyword).copied()
    }

    pub fn iter_delta_keywords(&self) -> impl ExactSizeIterator<Item = (Entity, &str)> {
        self.delta_keywords
            .iter()
            .map(|(entity, keyword)| (*entity, keyword.as_str()))
    }

    pub fn atom_alias(&self, name: &str) -> Option<&AtomDsl> {
        self.atom_aliases
            .get_by_left(name)
            .map(Box::as_ref)
            .or_else(|| self.lhs.atom_alias(name))
    }

    pub fn atom_alias_name(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases
            .get_by_right(dsl)
            .map(String::as_str)
            .or_else(|| self.lhs.atom_alias_name(dsl))
    }

    pub fn iter_reaction_atom_aliases(&self) -> impl ExactSizeIterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_delta_keyword(
        &mut self,
        entity: Entity,
        name: impl Into<String>,
    ) -> Result<(), MetadataError> {
        let name = name.into();

        if self.delta_keywords.get_by_left(&entity) == Some(&name) {
            return Ok(());
        }
        if self.delta_keywords.contains_right(name.as_str())
            || self.lhs.entity(name.as_str()).is_some()
            || self.lhs.atom_alias(name.as_str()).is_some()
            || self.atom_aliases.contains_left(name.as_str())
        {
            return Err(MetadataError::DuplicateKeyword(name));
        }

        self.delta_keywords.insert(entity, name);
        Ok(())
    }

    pub fn add_atom_alias(
        &mut self,
        name: impl Into<String>,
        atom: impl Into<AtomDsl>,
    ) -> Result<(), MetadataError> {
        let name = name.into();
        let atom = Box::new(atom.into());

        if self.entity(name.as_str()).is_some() {
            return Err(MetadataError::DuplicateKeyword(name));
        }
        if let Some(existing) = self.atom_aliases.get_by_left(name.as_str()) {
            return if existing == &atom {
                Ok(())
            } else {
                Err(MetadataError::DuplicateKeyword(name))
            };
        }
        if self.lhs.atom_alias(name.as_str()).is_some() {
            return Err(MetadataError::DuplicateKeyword(name));
        }
        if let Some(existing) = self.atom_aliases.get_by_right(atom.as_ref()) {
            return Err(MetadataError::DuplicateAtomAlias(existing.clone()));
        }
        if let Some(existing) = self.lhs.atom_alias_name(atom.as_ref()) {
            return Err(MetadataError::DuplicateAtomAlias(existing.to_string()));
        }

        self.atom_aliases.insert(name, atom);
        Ok(())
    }
}

impl Metadata for ReactionMetadata {
    fn keyword(&self, entity: Entity) -> Option<&str> {
        self.keyword(entity)
    }

    fn entity(&self, keyword: &str) -> Option<Entity> {
        self.entity(keyword)
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
    #[case::empty(MoleculeMetadata::new(), Vec::new())]
    #[case::populated(
        MoleculeMetadata {
            keywords: [
                (Entity::Atom(AtomId(0)), "atom".to_string()),
                (Entity::Bond(BondId(0)), "bond".to_string()),
            ]
            .into_iter()
            .collect(),
            atom_aliases: BiBTreeMap::new(),
        },
        vec![
            (Entity::Atom(AtomId(0)), "atom"),
            (Entity::Bond(BondId(0)), "bond"),
        ]
    )]
    fn test_molecule_metadata_iter_keywords(
        #[case] metadata: MoleculeMetadata,
        #[case] expected: Vec<(Entity, &str)>,
    ) {
        let keywords = metadata.iter_keywords();

        assert_eq!(keywords.len(), expected.len());
        assert_eq!(keywords.collect::<Vec<_>>(), expected);
    }

    #[rstest]
    #[case::present("carbon", Some(AtomDsl(AtomAst::from_element(Element::C))))]
    #[case::absent("nitrogen", None)]
    fn test_molecule_metadata_atom_alias(#[case] name: &str, #[case] expected: Option<AtomDsl>) {
        let metadata = MoleculeMetadata {
            keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(metadata.atom_alias(name), expected.as_ref());
    }

    #[rstest]
    #[case::present(AtomDsl(AtomAst::from_element(Element::C)), Some("carbon"))]
    #[case::absent(AtomDsl(AtomAst::from_element(Element::N)), None)]
    fn test_molecule_metadata_atom_alias_name(
        #[case] atom: AtomDsl,
        #[case] expected: Option<&str>,
    ) {
        let metadata = MoleculeMetadata {
            keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(metadata.atom_alias_name(&atom), expected);
    }

    #[rstest]
    #[case::empty(MoleculeMetadata::new(), Vec::new())]
    #[case::populated(
        MoleculeMetadata {
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
        },
        vec![
            ("carbon", AtomDsl(AtomAst::from_element(Element::C))),
            ("nitrogen", AtomDsl(AtomAst::from_element(Element::N))),
        ]
    )]
    fn test_molecule_metadata_iter_atom_aliases(
        #[case] metadata: MoleculeMetadata,
        #[case] expected: Vec<(&str, AtomDsl)>,
    ) {
        let aliases = metadata.iter_atom_aliases();

        assert_eq!(aliases.len(), expected.len());
        assert_eq!(
            aliases
                .map(|(name, atom)| (name, atom.clone()))
                .collect::<Vec<_>>(),
            expected
        );
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
    #[case::keyword_then_alias(false)]
    #[case::alias_then_keyword(true)]
    fn test_molecule_metadata_keyword_alias_collision(#[case] alias_first: bool) {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let mut actual = MoleculeMetadata::new();
        let result = if alias_first {
            actual.add_atom_alias("carbon", atom.clone()).unwrap();
            actual.set_keyword(Entity::Atom(AtomId(0)), "carbon")
        } else {
            actual
                .set_keyword(Entity::Atom(AtomId(0)), "carbon")
                .unwrap();
            actual.add_atom_alias("carbon", atom.clone())
        };

        assert_eq!(
            result,
            Err(MetadataError::DuplicateKeyword("carbon".to_string()))
        );
        assert_eq!(
            actual.entity("carbon"),
            (!alias_first).then_some(Entity::Atom(AtomId(0)))
        );
        assert_eq!(actual.atom_alias("carbon"), alias_first.then_some(&atom));
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
    fn test_reaction_metadata_lhs() {
        let lhs = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "lhs".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };
        let metadata = ReactionMetadata {
            lhs: lhs.clone(),
            ..Default::default()
        };

        assert_eq!(metadata.lhs(), &lhs);
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)))]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(1))
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(1))
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Entity::MulticenterBond(MulticenterBondId(1))
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1))
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(1))
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        Entity::StereoBond(StereoBondId(1))
    )]
    fn test_reaction_metadata_keyword(#[case] lhs_entity: Entity, #[case] delta_entity: Entity) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(lhs_entity, "lhs".to_string())].into_iter().collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [(delta_entity, "delta".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.keyword(lhs_entity), Some("lhs"));
        assert_eq!(metadata.keyword(delta_entity), Some("delta"));
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)))]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(1))
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(1))
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Entity::MulticenterBond(MulticenterBondId(1))
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1))
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(1))
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        Entity::StereoBond(StereoBondId(1))
    )]
    fn test_reaction_metadata_entity(#[case] lhs_entity: Entity, #[case] delta_entity: Entity) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(lhs_entity, "lhs".to_string())].into_iter().collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [(delta_entity, "delta".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.entity("lhs"), Some(lhs_entity));
        assert_eq!(metadata.entity("delta"), Some(delta_entity));
        assert_eq!(metadata.entity("missing"), None);
    }

    #[rstest]
    #[case::empty(ReactionMetadata::default(), Vec::new())]
    #[case::populated(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(Entity::Atom(AtomId(0)), "lhs".to_string())]
                    .into_iter()
                    .collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [(Entity::Bond(BondId(0)), "delta".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        },
        vec![
            (Entity::Bond(BondId(0)), "delta"),
            (Entity::Atom(AtomId(0)), "lhs"),
        ]
    )]
    fn test_reaction_metadata_iter_keywords(
        #[case] metadata: ReactionMetadata,
        #[case] expected: Vec<(Entity, &str)>,
    ) {
        let mut keywords = metadata.iter_keywords();

        assert_eq!(keywords.len(), expected.len());
        assert_eq!(keywords.size_hint(), (expected.len(), Some(expected.len())));
        assert_eq!(keywords.next(), expected.first().copied());
        assert_eq!(keywords.len(), expected.len().saturating_sub(1));
        assert_eq!(
            keywords.by_ref().collect::<Vec<_>>(),
            expected.get(1..).unwrap_or_default()
        );
        assert_eq!(keywords.len(), 0);
        assert_eq!(keywords.size_hint(), (0, Some(0)));
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)))]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(1))
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(1))
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Entity::MulticenterBond(MulticenterBondId(1))
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1))
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(1))
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        Entity::StereoBond(StereoBondId(1))
    )]
    fn test_reaction_metadata_delta_keyword(
        #[case] lhs_entity: Entity,
        #[case] delta_entity: Entity,
    ) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(lhs_entity, "lhs".to_string())].into_iter().collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [(delta_entity, "delta".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.delta_keyword(delta_entity), Some("delta"));
        assert_eq!(metadata.delta_keyword(lhs_entity), None);
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)))]
    #[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)))]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        Entity::DativeBond(DativeBondId(1))
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Entity::AromaticSystem(AromaticSystemId(1))
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Entity::MulticenterBond(MulticenterBondId(1))
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1))
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        Entity::StereoAtom(StereoAtomId(1))
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        Entity::StereoBond(StereoBondId(1))
    )]
    fn test_reaction_metadata_delta_entity(
        #[case] lhs_entity: Entity,
        #[case] delta_entity: Entity,
    ) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(lhs_entity, "lhs".to_string())].into_iter().collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [(delta_entity, "delta".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(metadata.delta_entity("delta"), Some(delta_entity));
        assert_eq!(metadata.delta_entity("lhs"), None);
        assert_eq!(metadata.delta_entity("missing"), None);
    }

    #[rstest]
    #[case::empty(ReactionMetadata::default(), Vec::new())]
    #[case::populated(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(Entity::Atom(AtomId(0)), "lhs".to_string())]
                    .into_iter()
                    .collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: [
                (Entity::Atom(AtomId(1)), "delta-atom".to_string()),
                (Entity::Bond(BondId(0)), "delta-bond".to_string()),
            ]
            .into_iter()
            .collect(),
            atom_aliases: BiBTreeMap::new(),
        },
        vec![
            (Entity::Atom(AtomId(1)), "delta-atom"),
            (Entity::Bond(BondId(0)), "delta-bond"),
        ]
    )]
    fn test_reaction_metadata_iter_delta_keywords(
        #[case] metadata: ReactionMetadata,
        #[case] expected: Vec<(Entity, &str)>,
    ) {
        let mut keywords = metadata.iter_delta_keywords();

        assert_eq!(keywords.len(), expected.len());
        assert_eq!(keywords.size_hint(), (expected.len(), Some(expected.len())));
        assert_eq!(keywords.next(), expected.first().copied());
        assert_eq!(keywords.len(), expected.len().saturating_sub(1));
        assert_eq!(
            keywords.by_ref().collect::<Vec<_>>(),
            expected.get(1..).unwrap_or_default()
        );
        assert_eq!(keywords.len(), 0);
        assert_eq!(keywords.size_hint(), (0, Some(0)));
    }

    #[rstest]
    #[case::reaction("reaction", Some(AtomDsl(AtomAst::from_element(Element::C))))]
    #[case::lhs("lhs", Some(AtomDsl(AtomAst::from_element(Element::N))))]
    #[case::missing("missing", None)]
    fn test_reaction_metadata_atom_alias(#[case] name: &str, #[case] expected: Option<AtomDsl>) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [(
                    "lhs".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::N))),
                )]
                .into_iter()
                .collect(),
            },
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "reaction".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(metadata.atom_alias(name), expected.as_ref());
    }

    #[rstest]
    #[case::reaction(AtomDsl(AtomAst::from_element(Element::C)), Some("reaction"))]
    #[case::lhs(AtomDsl(AtomAst::from_element(Element::N)), Some("lhs"))]
    #[case::missing(AtomDsl(AtomAst::from_element(Element::O)), None)]
    fn test_reaction_metadata_atom_alias_name(
        #[case] atom: AtomDsl,
        #[case] expected: Option<&str>,
    ) {
        let metadata = ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [(
                    "lhs".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::N))),
                )]
                .into_iter()
                .collect(),
            },
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "reaction".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };

        assert_eq!(metadata.atom_alias_name(&atom), expected);
    }

    #[rstest]
    #[case::empty(ReactionMetadata::default(), Vec::new())]
    #[case::populated(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [(
                    "lhs".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::O))),
                )]
                .into_iter()
                .collect(),
            },
            delta_keywords: BiBTreeMap::new(),
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
        },
        vec![
            ("carbon", AtomDsl(AtomAst::from_element(Element::C))),
            ("nitrogen", AtomDsl(AtomAst::from_element(Element::N))),
        ]
    )]
    fn test_reaction_metadata_iter_reaction_atom_aliases(
        #[case] metadata: ReactionMetadata,
        #[case] expected: Vec<(&str, AtomDsl)>,
    ) {
        let mut aliases = metadata.iter_reaction_atom_aliases();

        assert_eq!(aliases.len(), expected.len());
        assert_eq!(aliases.size_hint(), (expected.len(), Some(expected.len())));
        assert_eq!(
            aliases.next().map(|(name, atom)| (name, atom.clone())),
            expected.first().cloned()
        );
        assert_eq!(aliases.len(), expected.len().saturating_sub(1));
        assert_eq!(
            aliases
                .by_ref()
                .map(|(name, atom)| (name, atom.clone()))
                .collect::<Vec<_>>(),
            expected.get(1..).unwrap_or_default()
        );
        assert_eq!(aliases.len(), 0);
        assert_eq!(aliases.size_hint(), (0, Some(0)));
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
    fn test_reaction_metadata_set_delta_keyword(#[case] entity: Entity) {
        let mut actual = ReactionMetadata::default();
        let result = actual.set_delta_keyword(entity, "key");
        let expected = ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: [(entity, "key".to_string())].into_iter().collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_reaction_metadata_set_delta_keyword_idempotent() {
        let mut actual = ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: [(Entity::Atom(AtomId(1)), "key".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };
        let expected = actual.clone();

        let result = actual.set_delta_keyword(Entity::Atom(AtomId(1)), "key");

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_reaction_metadata_set_delta_keyword_rebinding() {
        let mut actual = ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: [(Entity::Atom(AtomId(1)), "old".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };
        let expected = ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: [(Entity::Atom(AtomId(1)), "new".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        };

        let result = actual.set_delta_keyword(Entity::Atom(AtomId(1)), "new");

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::delta_keyword(ReactionMetadata {
        lhs: MoleculeMetadata::new(),
        delta_keywords: [(Entity::Atom(AtomId(1)), "used".to_string())]
            .into_iter()
            .collect(),
        atom_aliases: BiBTreeMap::new(),
    })]
    #[case::lhs_keyword(ReactionMetadata {
        lhs: MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "used".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        },
        delta_keywords: BiBTreeMap::new(),
        atom_aliases: BiBTreeMap::new(),
    })]
    #[case::lhs_alias(ReactionMetadata {
        lhs: MoleculeMetadata {
            keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "used".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        },
        delta_keywords: BiBTreeMap::new(),
        atom_aliases: BiBTreeMap::new(),
    })]
    #[case::reaction_alias(ReactionMetadata {
        lhs: MoleculeMetadata::new(),
        delta_keywords: BiBTreeMap::new(),
        atom_aliases: [(
            "used".to_string(),
            Box::new(AtomDsl(AtomAst::from_element(Element::C))),
        )]
        .into_iter()
        .collect(),
    })]
    fn test_reaction_metadata_set_delta_keyword_error(#[case] mut actual: ReactionMetadata) {
        let expected = actual.clone();

        let result = actual.set_delta_keyword(Entity::Bond(BondId(1)), "used");

        assert_eq!(
            result,
            Err(MetadataError::DuplicateKeyword("used".to_string()))
        );
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_reaction_metadata_add_atom_alias() {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let mut actual = ReactionMetadata::default();

        let result = actual.add_atom_alias("carbon", atom.clone());

        assert_eq!(result, Ok(()));
        assert_eq!(
            actual,
            ReactionMetadata {
                lhs: MoleculeMetadata::new(),
                delta_keywords: BiBTreeMap::new(),
                atom_aliases: [("carbon".to_string(), Box::new(atom))]
                    .into_iter()
                    .collect(),
            }
        );
    }

    #[rstest]
    fn test_reaction_metadata_add_atom_alias_identity() {
        let atom = AtomDsl(AtomAst::from_element(Element::C));
        let mut actual = ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: [("carbon".to_string(), Box::new(atom.clone()))]
                .into_iter()
                .collect(),
        };
        let expected = actual.clone();

        let result = actual.add_atom_alias("carbon", atom);

        assert_eq!(result, Ok(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::lhs_keyword(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: [(Entity::Atom(AtomId(0)), "used".to_string())]
                    .into_iter()
                    .collect(),
                atom_aliases: BiBTreeMap::new(),
            },
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: BiBTreeMap::new(),
        },
        "used",
        AtomDsl(AtomAst::from_element(Element::C)),
        MetadataError::DuplicateKeyword("used".to_string())
    )]
    #[case::delta_keyword(
        ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: [(Entity::Atom(AtomId(1)), "used".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: BiBTreeMap::new(),
        },
        "used",
        AtomDsl(AtomAst::from_element(Element::C)),
        MetadataError::DuplicateKeyword("used".to_string())
    )]
    #[case::lhs_alias_name(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [(
                    "used".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::C))),
                )]
                .into_iter()
                .collect(),
            },
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: BiBTreeMap::new(),
        },
        "used",
        AtomDsl(AtomAst::from_element(Element::C)),
        MetadataError::DuplicateKeyword("used".to_string())
    )]
    #[case::reaction_alias_name(
        ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "used".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        },
        "used",
        AtomDsl(AtomAst::from_element(Element::N)),
        MetadataError::DuplicateKeyword("used".to_string())
    )]
    #[case::lhs_alias_target(
        ReactionMetadata {
            lhs: MoleculeMetadata {
                keywords: BiBTreeMap::new(),
                atom_aliases: [(
                    "used".to_string(),
                    Box::new(AtomDsl(AtomAst::from_element(Element::C))),
                )]
                .into_iter()
                .collect(),
            },
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: BiBTreeMap::new(),
        },
        "other",
        AtomDsl(AtomAst::from_element(Element::C)),
        MetadataError::DuplicateAtomAlias("used".to_string())
    )]
    #[case::reaction_alias_target(
        ReactionMetadata {
            lhs: MoleculeMetadata::new(),
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: [(
                "used".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        },
        "other",
        AtomDsl(AtomAst::from_element(Element::C)),
        MetadataError::DuplicateAtomAlias("used".to_string())
    )]
    fn test_reaction_metadata_add_atom_alias_error(
        #[case] mut actual: ReactionMetadata,
        #[case] name: &str,
        #[case] atom: AtomDsl,
        #[case] expected_error: MetadataError,
    ) {
        let expected = actual.clone();

        let result = actual.add_atom_alias(name, atom);

        assert_eq!(result, Err(expected_error));
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_reaction_metadata_from_molecule_metadata() {
        let lhs = MoleculeMetadata {
            keywords: [(Entity::Atom(AtomId(0)), "lhs".to_string())]
                .into_iter()
                .collect(),
            atom_aliases: [(
                "carbon".to_string(),
                Box::new(AtomDsl(AtomAst::from_element(Element::C))),
            )]
            .into_iter()
            .collect(),
        };
        let expected = ReactionMetadata {
            lhs: lhs.clone(),
            delta_keywords: BiBTreeMap::new(),
            atom_aliases: BiBTreeMap::new(),
        };

        assert_eq!(ReactionMetadata::from(lhs), expected);
    }
}
