//! Python values for persistent molecule and reaction DSL metadata.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph_ir::dsl::{
    MoleculeMetadata as GraphIrMoleculeMetadata, ReactionMetadata as GraphIrReactionMetadata,
};
use umol_graph_ir::ir::{
    AromaticSystemId as GraphIrAromaticSystemId, AtomId as GraphIrAtomId, BondId as GraphIrBondId,
    DativeBondId as GraphIrDativeBondId, Entity as GraphIrEntity,
    MulticenterBondId as GraphIrMulticenterBondId, NoncovalentBondId as GraphIrNoncovalentBondId,
    StereoAtomId as GraphIrStereoAtomId, StereoBondId as GraphIrStereoBondId,
};

use crate::convert::{hash_rust, variant_repr};
use crate::correspondence::MoleculeCorrespondence;
use crate::error::metadata_error;

/// A typed numerical identifier for one of the eight molecule entity kinds.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Entity {
    Atom(u32),
    Bond(u32),
    DativeBond(u32),
    AromaticSystem(u32),
    MulticenterBond(u32),
    NoncovalentBond(u32),
    StereoAtom(u32),
    StereoBond(u32),
}

#[pymethods]
impl Entity {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __hash__(&self) -> u64 {
        hash_rust(self)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            Self::Atom(_) => "Atom",
            Self::Bond(_) => "Bond",
            Self::DativeBond(_) => "DativeBond",
            Self::AromaticSystem(_) => "AromaticSystem",
            Self::MulticenterBond(_) => "MulticenterBond",
            Self::NoncovalentBond(_) => "NoncovalentBond",
            Self::StereoAtom(_) => "StereoAtom",
            Self::StereoBond(_) => "StereoBond",
        };
        variant_repr(slf.bind(py).as_any(), "Entity", variant, 1)
    }
}

impl Entity {
    pub(crate) fn from_rust(entity: GraphIrEntity) -> Self {
        match entity {
            GraphIrEntity::Atom(id) => Self::Atom(id.0),
            GraphIrEntity::Bond(id) => Self::Bond(id.0),
            GraphIrEntity::DativeBond(id) => Self::DativeBond(id.0),
            GraphIrEntity::AromaticSystem(id) => Self::AromaticSystem(id.0),
            GraphIrEntity::MulticenterBond(id) => Self::MulticenterBond(id.0),
            GraphIrEntity::NoncovalentBond(id) => Self::NoncovalentBond(id.0),
            GraphIrEntity::StereoAtom(id) => Self::StereoAtom(id.0),
            GraphIrEntity::StereoBond(id) => Self::StereoBond(id.0),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrEntity {
        match self {
            Self::Atom(id) => GraphIrEntity::Atom(GraphIrAtomId(id)),
            Self::Bond(id) => GraphIrEntity::Bond(GraphIrBondId(id)),
            Self::DativeBond(id) => GraphIrEntity::DativeBond(GraphIrDativeBondId(id)),
            Self::AromaticSystem(id) => GraphIrEntity::AromaticSystem(GraphIrAromaticSystemId(id)),
            Self::MulticenterBond(id) => {
                GraphIrEntity::MulticenterBond(GraphIrMulticenterBondId(id))
            }
            Self::NoncovalentBond(id) => {
                GraphIrEntity::NoncovalentBond(GraphIrNoncovalentBondId(id))
            }
            Self::StereoAtom(id) => GraphIrEntity::StereoAtom(GraphIrStereoAtomId(id)),
            Self::StereoBond(id) => GraphIrEntity::StereoBond(GraphIrStereoBondId(id)),
        }
    }

    fn repr(self) -> String {
        match self {
            Self::Atom(id) => format!("Entity.Atom({id})"),
            Self::Bond(id) => format!("Entity.Bond({id})"),
            Self::DativeBond(id) => format!("Entity.DativeBond({id})"),
            Self::AromaticSystem(id) => format!("Entity.AromaticSystem({id})"),
            Self::MulticenterBond(id) => format!("Entity.MulticenterBond({id})"),
            Self::NoncovalentBond(id) => format!("Entity.NoncovalentBond({id})"),
            Self::StereoAtom(id) => format!("Entity.StereoAtom({id})"),
            Self::StereoBond(id) => format!("Entity.StereoBond({id})"),
        }
    }
}

/// Persistent entity-keyword and atom-alias metadata for a molecule DSL value.
///
/// The metadata may be detached from its molecule for transport, but
/// `Molecule.render_with_metadata` checks that every keyword still names an
/// entity present in the molecule.
#[pyclass(eq, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeMetadata(GraphIrMoleculeMetadata);

#[pymethods]
impl MoleculeMetadata {
    #[new]
    fn new() -> Self {
        Self(GraphIrMoleculeMetadata::new())
    }

    /// Keyword assigned to `entity`, if any.
    fn keyword(&self, entity: Entity) -> Option<String> {
        self.0.keyword(entity.to_rust()).map(str::to_owned)
    }

    /// Entity assigned to `keyword`, if any.
    fn entity(&self, keyword: &str) -> Option<Entity> {
        self.0.entity(keyword).map(Entity::from_rust)
    }

    /// Assign `keyword` to `entity` without violating metadata namespace invariants.
    fn set_keyword(&mut self, entity: Entity, keyword: String) -> PyResult<()> {
        self.0
            .set_keyword(entity.to_rust(), keyword)
            .map_err(metadata_error)
    }

    /// Move entity keywords through a molecule correspondence.
    fn remap(&self, correspondence: &MoleculeCorrespondence) -> Self {
        Self(self.0.clone().remap(correspondence.to_rust()))
    }

    fn __repr__(&self) -> String {
        self.repr()
    }
}

impl MoleculeMetadata {
    pub(crate) fn from_rust(metadata: GraphIrMoleculeMetadata) -> Self {
        Self(metadata)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrMoleculeMetadata {
        &self.0
    }

    fn repr(&self) -> String {
        format!(
            "MoleculeMetadata(keywords={}, atom_alias_count={})",
            keyword_repr(self.0.iter_keywords()),
            self.0.iter_atom_aliases().len(),
        )
    }
}

/// Persistent lhs, delta-keyword, and atom-alias metadata for a reaction DSL value.
///
/// Combined lookup searches delta metadata before lhs metadata. Rendering
/// checks that lhs keywords name lhs entities and delta keywords name entities
/// introduced by matching add deltas.
#[pyclass(eq, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionMetadata(GraphIrReactionMetadata);

#[pymethods]
impl ReactionMetadata {
    #[new]
    #[pyo3(signature = (lhs=None))]
    fn new(lhs: Option<MoleculeMetadata>) -> Self {
        lhs.map_or_else(
            || Self(GraphIrReactionMetadata::default()),
            |lhs| Self(GraphIrReactionMetadata::from(lhs.0)),
        )
    }

    /// Detached lhs metadata snapshot.
    #[getter]
    fn lhs(&self) -> MoleculeMetadata {
        MoleculeMetadata::from_rust(self.0.lhs().clone())
    }

    /// Keyword assigned to `entity` in the delta or lhs scope, if any.
    fn keyword(&self, entity: Entity) -> Option<String> {
        self.0.keyword(entity.to_rust()).map(str::to_owned)
    }

    /// Entity assigned to `keyword` in the delta or lhs scope, if any.
    fn entity(&self, keyword: &str) -> Option<Entity> {
        self.0.entity(keyword).map(Entity::from_rust)
    }

    /// Delta-scope keyword assigned to `entity`, if any.
    fn delta_keyword(&self, entity: Entity) -> Option<String> {
        self.0.delta_keyword(entity.to_rust()).map(str::to_owned)
    }

    /// Delta-scope entity assigned to `keyword`, if any.
    fn delta_entity(&self, keyword: &str) -> Option<Entity> {
        self.0.delta_entity(keyword).map(Entity::from_rust)
    }

    /// Assign a delta-scope keyword without violating reaction metadata namespace invariants.
    fn set_delta_keyword(&mut self, entity: Entity, keyword: String) -> PyResult<()> {
        self.0
            .set_delta_keyword(entity.to_rust(), keyword)
            .map_err(metadata_error)
    }

    fn __repr__(&self) -> String {
        format!(
            "ReactionMetadata(lhs={}, delta_keywords={}, reaction_atom_alias_count={})",
            MoleculeMetadata::from_rust(self.0.lhs().clone()).repr(),
            keyword_repr(self.0.iter_delta_keywords()),
            self.0.iter_reaction_atom_aliases().len(),
        )
    }
}

impl ReactionMetadata {
    pub(crate) fn from_rust(metadata: GraphIrReactionMetadata) -> Self {
        Self(metadata)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrReactionMetadata {
        &self.0
    }
}

fn keyword_repr<'a>(keywords: impl Iterator<Item = (GraphIrEntity, &'a str)>) -> String {
    let entries = keywords
        .map(|(entity, keyword)| format!("({}, {keyword:?})", Entity::from_rust(entity).repr(),))
        .collect::<Vec<_>>();
    format!("[{}]", entries.join(", "))
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyAnyMethods;
    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::Correspondence as GraphCoreCorrespondence;
    use umol_graph_ir::dsl::AtomDsl;
    use umol_graph_ir::ir::{
        AromaticSystemId, AtomForm, AtomId, BondId, DativeBondId,
        MoleculeCorrespondence as GraphIrMoleculeCorrespondence, MulticenterBondId,
        NoncovalentBondId, StereoAtomId, StereoBondId,
    };

    use super::*;
    use crate::error::MetadataError;

    #[rstest]
    #[case::atom(Entity::Atom(1), GraphIrEntity::Atom(AtomId(1)), "Entity.Atom(1)")]
    #[case::bond(Entity::Bond(2), GraphIrEntity::Bond(BondId(2)), "Entity.Bond(2)")]
    #[case::dative_bond(
        Entity::DativeBond(3),
        GraphIrEntity::DativeBond(DativeBondId(3)),
        "Entity.DativeBond(3)"
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(4),
        GraphIrEntity::AromaticSystem(AromaticSystemId(4)),
        "Entity.AromaticSystem(4)"
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(5),
        GraphIrEntity::MulticenterBond(MulticenterBondId(5)),
        "Entity.MulticenterBond(5)"
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(6),
        GraphIrEntity::NoncovalentBond(NoncovalentBondId(6)),
        "Entity.NoncovalentBond(6)"
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(7),
        GraphIrEntity::StereoAtom(StereoAtomId(7)),
        "Entity.StereoAtom(7)"
    )]
    #[case::stereo_bond(
        Entity::StereoBond(8),
        GraphIrEntity::StereoBond(StereoBondId(8)),
        "Entity.StereoBond(8)"
    )]
    fn test_entity_roundtrip(
        #[case] entity: Entity,
        #[case] expected: GraphIrEntity,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(entity.to_rust(), expected);
        assert_eq!(Entity::from_rust(expected), entity);
        assert_eq!(entity.repr(), expected_repr);
    }

    #[rstest]
    #[case::atom(Entity::Atom(0), "atom")]
    #[case::bond(Entity::Bond(0), "bond")]
    #[case::dative_bond(Entity::DativeBond(0), "dative")]
    #[case::aromatic_system(Entity::AromaticSystem(0), "aromatic")]
    #[case::multicenter_bond(Entity::MulticenterBond(0), "multicenter")]
    #[case::noncovalent_bond(Entity::NoncovalentBond(0), "noncovalent")]
    #[case::stereo_atom(Entity::StereoAtom(0), "stereo_atom")]
    #[case::stereo_bond(Entity::StereoBond(0), "stereo_bond")]
    fn test_molecule_metadata_set_keyword(#[case] entity: Entity, #[case] keyword: &str) {
        let mut metadata = MoleculeMetadata::new();

        metadata.set_keyword(entity, keyword.to_string()).unwrap();

        assert_eq!(metadata.keyword(entity), Some(keyword.to_string()));
        assert_eq!(metadata.entity(keyword), Some(entity));
        assert_eq!(metadata.to_rust().keyword(entity.to_rust()), Some(keyword));
    }

    #[rstest]
    fn test_molecule_metadata_set_keyword_error() {
        Python::attach(|py| {
            let mut metadata = MoleculeMetadata::new();
            metadata
                .set_keyword(Entity::Atom(0), "site".to_string())
                .unwrap();

            let error = metadata
                .set_keyword(Entity::Bond(0), "site".to_string())
                .unwrap_err();

            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "duplicate keyword: site"
            );
            assert_eq!(metadata.entity("site"), Some(Entity::Atom(0)));
            assert_eq!(metadata.keyword(Entity::Bond(0)), None);
        });
    }

    #[rstest]
    fn test_molecule_metadata_remap() {
        let mut metadata = GraphIrMoleculeMetadata::new();
        for (entity, keyword) in [
            (GraphIrEntity::Atom(AtomId(0)), "atom"),
            (GraphIrEntity::Bond(BondId(0)), "bond"),
            (GraphIrEntity::DativeBond(DativeBondId(0)), "dative"),
            (
                GraphIrEntity::AromaticSystem(AromaticSystemId(0)),
                "aromatic",
            ),
            (
                GraphIrEntity::MulticenterBond(MulticenterBondId(0)),
                "multicenter",
            ),
            (
                GraphIrEntity::NoncovalentBond(NoncovalentBondId(0)),
                "noncovalent",
            ),
            (GraphIrEntity::StereoAtom(StereoAtomId(0)), "stereo_atom"),
            (GraphIrEntity::StereoBond(StereoBondId(0)), "stereo_bond"),
        ] {
            metadata.set_keyword(entity, keyword).unwrap();
        }
        metadata
            .add_atom_alias("carbon", AtomDsl(AtomForm::from_element(Element::C)))
            .unwrap();
        let correspondence = MoleculeCorrespondence::from_rust(GraphIrMoleculeCorrespondence::new(
            GraphCoreCorrespondence::new(vec![(AtomId(0), AtomId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(BondId(0), BondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(DativeBondId(0), DativeBondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(AromaticSystemId(0), AromaticSystemId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(MulticenterBondId(0), MulticenterBondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(StereoAtomId(0), StereoAtomId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ));

        let remapped = MoleculeMetadata::from_rust(metadata).remap(&correspondence);

        for (entity, keyword) in [
            (Entity::Atom(1), "atom"),
            (Entity::Bond(1), "bond"),
            (Entity::DativeBond(1), "dative"),
            (Entity::AromaticSystem(1), "aromatic"),
            (Entity::MulticenterBond(1), "multicenter"),
            (Entity::NoncovalentBond(1), "noncovalent"),
            (Entity::StereoAtom(1), "stereo_atom"),
            (Entity::StereoBond(1), "stereo_bond"),
        ] {
            assert_eq!(remapped.keyword(entity), Some(keyword.to_string()));
        }
        assert_eq!(
            remapped
                .to_rust()
                .atom_alias("carbon")
                .map(|alias| alias.0.clone()),
            Some(AtomForm::from_element(Element::C))
        );
    }

    #[rstest]
    #[case::empty(
        MoleculeMetadata::new(),
        "MoleculeMetadata(keywords=[], atom_alias_count=0)"
    )]
    #[case::populated(
        {
            let mut metadata = MoleculeMetadata::new();
            metadata.set_keyword(Entity::Atom(0), "carbon".to_string()).unwrap();
            metadata
        },
        r#"MoleculeMetadata(keywords=[(Entity.Atom(0), "carbon")], atom_alias_count=0)"#,
    )]
    fn test_molecule_metadata_repr(#[case] metadata: MoleculeMetadata, #[case] expected: &str) {
        assert_eq!(metadata.__repr__(), expected);
    }

    #[rstest]
    fn test_reaction_metadata_new() {
        let mut lhs = MoleculeMetadata::new();
        lhs.set_keyword(Entity::Atom(0), "lhs".to_string()).unwrap();

        let metadata = ReactionMetadata::new(Some(lhs.clone()));
        let mut detached_lhs = metadata.lhs();
        detached_lhs
            .set_keyword(Entity::Bond(0), "detached".to_string())
            .unwrap();

        assert_eq!(metadata.lhs(), lhs);
        assert_eq!(metadata.keyword(Entity::Atom(0)), Some("lhs".to_string()));
        assert_eq!(metadata.entity("lhs"), Some(Entity::Atom(0)));
        assert_eq!(metadata.delta_keyword(Entity::Atom(0)), None);
        assert_eq!(metadata.delta_entity("lhs"), None);
        assert_eq!(metadata.entity("detached"), None);
    }

    #[rstest]
    #[case::atom(Entity::Atom(0), "atom")]
    #[case::bond(Entity::Bond(0), "bond")]
    #[case::dative_bond(Entity::DativeBond(0), "dative")]
    #[case::aromatic_system(Entity::AromaticSystem(0), "aromatic")]
    #[case::multicenter_bond(Entity::MulticenterBond(0), "multicenter")]
    #[case::noncovalent_bond(Entity::NoncovalentBond(0), "noncovalent")]
    #[case::stereo_atom(Entity::StereoAtom(0), "stereo_atom")]
    #[case::stereo_bond(Entity::StereoBond(0), "stereo_bond")]
    fn test_reaction_metadata_set_delta_keyword(#[case] entity: Entity, #[case] keyword: &str) {
        let mut metadata = ReactionMetadata::new(None);

        metadata
            .set_delta_keyword(entity, keyword.to_string())
            .unwrap();

        assert_eq!(metadata.keyword(entity), Some(keyword.to_string()));
        assert_eq!(metadata.entity(keyword), Some(entity));
        assert_eq!(metadata.delta_keyword(entity), Some(keyword.to_string()));
        assert_eq!(metadata.delta_entity(keyword), Some(entity));
        assert_eq!(
            metadata.to_rust().delta_keyword(entity.to_rust()),
            Some(keyword)
        );
    }

    #[rstest]
    fn test_reaction_metadata_set_delta_keyword_error() {
        Python::attach(|py| {
            let mut lhs = MoleculeMetadata::new();
            lhs.set_keyword(Entity::Atom(0), "site".to_string())
                .unwrap();
            let mut metadata = ReactionMetadata::new(Some(lhs));

            let error = metadata
                .set_delta_keyword(Entity::Bond(0), "site".to_string())
                .unwrap_err();

            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "duplicate keyword: site"
            );
            assert_eq!(metadata.entity("site"), Some(Entity::Atom(0)));
            assert_eq!(metadata.delta_entity("site"), None);
            assert_eq!(metadata.delta_keyword(Entity::Bond(0)), None);
        });
    }

    #[rstest]
    fn test_reaction_metadata_alias_preservation() {
        let mut lhs = GraphIrMoleculeMetadata::new();
        lhs.add_atom_alias("carbon", AtomDsl(AtomForm::from_element(Element::C)))
            .unwrap();
        let mut metadata = GraphIrReactionMetadata::from(lhs);
        metadata
            .add_atom_alias("nitrogen", AtomDsl(AtomForm::from_element(Element::N)))
            .unwrap();

        let metadata = ReactionMetadata::from_rust(metadata);

        assert_eq!(
            metadata
                .to_rust()
                .atom_alias("carbon")
                .map(|alias| alias.0.clone()),
            Some(AtomForm::from_element(Element::C))
        );
        assert_eq!(
            metadata
                .to_rust()
                .atom_alias("nitrogen")
                .map(|alias| alias.0.clone()),
            Some(AtomForm::from_element(Element::N))
        );
    }

    #[rstest]
    #[case::empty(
        ReactionMetadata::new(None),
        concat!(
            "ReactionMetadata(lhs=MoleculeMetadata(keywords=[], atom_alias_count=0), ",
            "delta_keywords=[], reaction_atom_alias_count=0)"
        ),
    )]
    #[case::scoped(
        {
            let mut lhs = MoleculeMetadata::new();
            lhs.set_keyword(Entity::Atom(0), "lhs".to_string()).unwrap();
            let mut metadata = ReactionMetadata::new(Some(lhs));
            metadata
                .set_delta_keyword(Entity::Bond(1), "delta".to_string())
                .unwrap();
            metadata
        },
        concat!(
            "ReactionMetadata(lhs=MoleculeMetadata(keywords=[(Entity.Atom(0), \"lhs\")], ",
            "atom_alias_count=0), delta_keywords=[(Entity.Bond(1), \"delta\")], ",
            "reaction_atom_alias_count=0)"
        ),
    )]
    fn test_reaction_metadata_repr(#[case] metadata: ReactionMetadata, #[case] expected: &str) {
        assert_eq!(metadata.__repr__(), expected);
    }
}
