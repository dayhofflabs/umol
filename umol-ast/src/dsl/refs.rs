//! Surface-level entity references. Each ref is a positional index (`Edn::Int`),
//! a symbolic id keyword (`Edn::Keyword`), or — for non-atom entities — a
//! *structural* form (`Edn::Map`) naming the entity by its constituent atoms /
//! bonds. `resolve` turns a ref into an AST id against the parse-time
//! `MoleculeNamespace` (count for index bounds, `by_name` for id keywords,
//! `by_participants` for the structural form); `from_ast` renders an id back to a
//! ref against the `MoleculeMetadata` roundtrip projection.

// Only atom- and bond-ref resolution (the molecule entry loops) consumes the namespace path so far;
// the other five refs' `resolve` and their structural resolvers are wired when constraint /
// relational / reaction resolution migrates off `into_ast(metadata)`.
#![allow(dead_code)]

use umol_edn::{DeError, Edn, EdnError, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::edn_utils::{atoms_pair, atoms_vec, eof_err, parse_vec, required_key};
use super::error::ParseError;
use super::molecule::MoleculeMetadata;
use super::namespace::MoleculeNamespace;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::{StereoLigand, StereoLigandKind};

macro_rules! define_ref {
    ($name:ident, $id:ident, $accessor:ident, $kind:literal, $reader:ident,
        $count:ident, $by_name:ident
        $(, structural = $payload:ty, $parse_structural:ident, $resolve_structural:ident)?) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
            $( Structural($payload), )?
        }

        impl $name {
            /// Build a ref from an AST index, preferring an id from `metadata`
            /// if one is recorded for this index.
            pub fn from_ast(id: $id, metadata: &MoleculeMetadata) -> Self {
                if let Some(name) = metadata.$accessor(id) {
                    Self::Id(name.to_string())
                } else {
                    Self::Index(id.index())
                }
            }

            /// Resolve this ref to an AST index against `metadata`. Fails on
            /// unknown id or out-of-range numeric index.
            pub fn into_ast(
                self,
                count: usize,
                metadata: &MoleculeMetadata,
            ) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => {
                        for i in 0..count {
                            let id = $id::from(i);
                            if metadata.$accessor(id) == Some(name.as_str()) {
                                return Ok(id);
                            }
                        }
                        Err(ParseError::InvalidRef {
                            kind: $kind,
                            value: name,
                        })
                    }
                    $( Self::Structural(_) => {
                        let _phantom: fn($payload) = |_| {};
                        Err(ParseError::InvalidRef {
                            kind: $kind,
                            value: "structural".to_string(),
                        })
                    } )?
                }
            }

            /// Resolve this ref to an AST id against the parse-time `namespace`
            /// (the source of truth: count for index bounds, `by_name` for id
            /// keywords, `by_participants` for the structural form).
            pub(crate) fn resolve(self, namespace: &MoleculeNamespace) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < namespace.$count() {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => {
                        namespace.$by_name(&name).ok_or(ParseError::InvalidRef {
                            kind: $kind,
                            value: name,
                        })
                    }
                    $( Self::Structural(participants) => {
                        $resolve_structural(participants, namespace)
                    } )?
                }
            }
        }

        impl<'de> FromEdn<'de> for $name {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                match edn {
                    Edn::Int(n) => {
                        let i = usize::try_from(*n).map_err(|_| DeError::OutOfRange {
                            value: n.to_string(),
                            target: "usize",
                            path: Vec::new(),
                        })?;
                        Ok(Self::Index(i))
                    }
                    Edn::Keyword(k) => Ok(Self::Id(k.name().to_string())),
                    $( Edn::Map(m) => {
                        if m.get_keyword("type").is_some() || m.get_keyword("id").is_some() {
                            return Err(DeError::Custom(format!(
                                "{} structural ref must not carry :type or :id",
                                $kind
                            )));
                        }
                        Ok(Self::Structural($parse_structural(m)?))
                    } )?
                    other => {
                        #[allow(unused_mut)]
                        let mut expected = concat!($kind, " ref (int or keyword)");
                        $(
                            let _: fn($payload) = |_| {};
                            expected = concat!($kind, " ref (int, keyword, or structural map)");
                        )?
                        Err(DeError::TypeMismatch {
                            expected,
                            got: other.kind(),
                            path: Vec::new(),
                        })
                    }
                }
            }
        }

        impl ToEdn for $name {
            fn to_edn(&self) -> Edn<'static> {
                match self {
                    Self::Index(i) => Edn::Int(*i as i64),
                    Self::Id(name) => Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone())),
                    $( Self::Structural(_) => {
                        let _phantom: fn($payload) = |_| {};
                        unreachable!(concat!(
                            $kind,
                            " structural refs are input-only and never rendered"
                        ))
                    } )?
                }
            }
        }

        pub(super) fn $reader(de: &mut EdnStreamDeserializer<'_>) -> Result<$name, EdnError> {
            match de.peek_byte()?.ok_or_else(eof_err)? {
                b':' => Ok($name::Id(de.read_keyword_name()?.into_owned())),
                _ => {
                    let n = de.read_i64()?;
                    let i = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                        value: n.to_string(),
                        target: "usize",
                        path: Vec::new(),
                    })?;
                    Ok($name::Index(i))
                }
            }
        }
    };
}

define_ref!(
    AtomRef,
    AtomId,
    atom_id,
    "atom",
    read_atom_ref,
    atom_count,
    atom_by_name
);
define_ref!(
    BondRef,
    BondId,
    bond_id,
    "bond",
    read_bond_ref,
    bond_count,
    bond_by_name,
    structural = [AtomRef; 2],
    parse_bond_structural,
    resolve_bond_structural
);
define_ref!(
    DativeBondRef,
    DativeBondId,
    dative_bond_id,
    "dative-bond",
    read_dative_bond_ref,
    dative_bond_count,
    dative_bond_by_name,
    structural = DativeBondParticipants,
    parse_dative_structural,
    resolve_dative_structural
);
define_ref!(
    AromaticSystemRef, AromaticSystemId, aromatic_system_id, "aromatic-system",
    read_aromatic_system_ref, aromatic_system_count, aromatic_system_by_name,
    structural = Vec<AtomRef>, parse_aromatic_structural, resolve_aromatic_structural
);
define_ref!(
    MulticenterBondRef, MulticenterBondId, multicenter_bond_id, "multicenter-bond",
    read_multicenter_bond_ref, multicenter_bond_count, multicenter_bond_by_name,
    structural = Vec<AtomRef>, parse_multicenter_structural, resolve_multicenter_structural
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondId,
    noncovalent_bond_id,
    "noncovalent-bond",
    read_noncovalent_bond_ref,
    noncovalent_bond_count,
    noncovalent_bond_by_name,
    structural = [AtomRef; 2],
    parse_noncovalent_structural,
    resolve_noncovalent_structural
);
define_ref!(
    StereoAtomRef,
    StereoAtomId,
    stereo_atom_id,
    "stereo-atom",
    read_stereo_atom_ref,
    stereo_atom_count,
    stereo_atom_by_name,
    structural = StereoAtomParticipants,
    parse_stereo_atom_structural,
    resolve_stereo_atom_structural
);
define_ref!(
    StereoBondRef,
    StereoBondId,
    stereo_bond_id,
    "stereo-bond",
    read_stereo_bond_ref,
    stereo_bond_count,
    stereo_bond_by_name,
    structural = StereoBondParticipants,
    parse_stereo_bond_structural,
    resolve_stereo_bond_structural
);

/// The constituent atoms of a dative bond named structurally (`{:donors [..] :acceptor a}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DativeBondParticipants {
    pub(crate) donors: Vec<AtomRef>,
    pub(crate) acceptor: AtomRef,
}

/// The site + ligand frame of a stereo atom named structurally (`{:site a :ligands [..]}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoAtomParticipants {
    pub(crate) site: AtomRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
}

/// The site bond + ligand frame of a stereo bond named structurally
/// (`{:site bond-ref :ligands [..]}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoBondParticipants {
    pub(crate) site: BondRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
}

fn parse_bond_structural(m: &EdnMap<'_>) -> Result<[AtomRef; 2], DeError> {
    atoms_pair(m, "bond structural ref")
}

fn parse_noncovalent_structural(m: &EdnMap<'_>) -> Result<[AtomRef; 2], DeError> {
    atoms_pair(m, "noncovalent-bond structural ref")
}

fn parse_aromatic_structural(m: &EdnMap<'_>) -> Result<Vec<AtomRef>, DeError> {
    atoms_vec(m, "aromatic-system structural ref")
}

fn parse_multicenter_structural(m: &EdnMap<'_>) -> Result<Vec<AtomRef>, DeError> {
    atoms_vec(m, "multicenter-bond structural ref")
}

fn parse_dative_structural(m: &EdnMap<'_>) -> Result<DativeBondParticipants, DeError> {
    Ok(DativeBondParticipants {
        // `AtomRef::from_edn` can't be passed directly here — HRTB inference fails
        // ("FnMut not general enough"); the closure is required despite clippy.
        #[allow(clippy::redundant_closure)]
        donors: parse_vec(
            required_key(m, "donors", "dative-bond structural ref")?,
            ":donors",
            |e| AtomRef::from_edn(e),
        )?,
        acceptor: AtomRef::from_edn(required_key(m, "acceptor", "dative-bond structural ref")?)?,
    })
}

fn parse_stereo_atom_structural(m: &EdnMap<'_>) -> Result<StereoAtomParticipants, DeError> {
    Ok(StereoAtomParticipants {
        site: AtomRef::from_edn(required_key(m, "site", "stereo-atom structural ref")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-atom structural ref")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
    })
}

fn parse_stereo_bond_structural(m: &EdnMap<'_>) -> Result<StereoBondParticipants, DeError> {
    Ok(StereoBondParticipants {
        site: BondRef::from_edn(required_key(m, "site", "stereo-bond structural ref")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-bond structural ref")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
    })
}

/// Resolve a vector of atom refs against the namespace, preserving order.
fn resolve_atom_refs(
    refs: Vec<AtomRef>,
    namespace: &MoleculeNamespace,
) -> Result<Vec<AtomId>, ParseError> {
    refs.into_iter().map(|r| r.resolve(namespace)).collect()
}

/// Resolve a stereo ligand frame (each ligand keeps its kind) against the namespace.
fn resolve_ligands(
    ligands: Vec<StereoLigandRef>,
    namespace: &MoleculeNamespace,
) -> Result<Vec<StereoLigand>, ParseError> {
    ligands
        .into_iter()
        .map(|l| Ok(StereoLigand::new(l.atom.resolve(namespace)?, l.kind)))
        .collect()
}

fn format_atom_ids(ids: &[AtomId]) -> String {
    let joined = ids
        .iter()
        .map(|id| id.index().to_string())
        .collect::<Vec<_>>()
        .join(" ");
    format!("[{joined}]")
}

fn resolve_bond_structural(
    atoms: [AtomRef; 2],
    namespace: &MoleculeNamespace,
) -> Result<BondId, ParseError> {
    let [a, b] = atoms;
    let a = a.resolve(namespace)?;
    let b = b.resolve(namespace)?;
    namespace
        .bond_by_participants(a, b)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "bond",
            value: format_atom_ids(&[a, b]),
        })
}

fn resolve_noncovalent_structural(
    atoms: [AtomRef; 2],
    namespace: &MoleculeNamespace,
) -> Result<NoncovalentBondId, ParseError> {
    let [a, b] = atoms;
    let a = a.resolve(namespace)?;
    let b = b.resolve(namespace)?;
    namespace
        .noncovalent_bond_by_participants(a, b)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "noncovalent-bond",
            value: format_atom_ids(&[a, b]),
        })
}

fn resolve_aromatic_structural(
    atoms: Vec<AtomRef>,
    namespace: &MoleculeNamespace,
) -> Result<AromaticSystemId, ParseError> {
    let atoms = resolve_atom_refs(atoms, namespace)?;
    namespace
        .aromatic_system_by_participants(&atoms)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "aromatic-system",
            value: format_atom_ids(&atoms),
        })
}

fn resolve_multicenter_structural(
    atoms: Vec<AtomRef>,
    namespace: &MoleculeNamespace,
) -> Result<MulticenterBondId, ParseError> {
    let atoms = resolve_atom_refs(atoms, namespace)?;
    namespace
        .multicenter_bond_by_participants(&atoms)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "multicenter-bond",
            value: format_atom_ids(&atoms),
        })
}

fn resolve_dative_structural(
    participants: DativeBondParticipants,
    namespace: &MoleculeNamespace,
) -> Result<DativeBondId, ParseError> {
    let donors = resolve_atom_refs(participants.donors, namespace)?;
    let acceptor = participants.acceptor.resolve(namespace)?;
    namespace
        .dative_bond_by_participants(&donors, acceptor)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "dative-bond",
            value: format_atom_ids(&donors),
        })
}

fn resolve_stereo_atom_structural(
    participants: StereoAtomParticipants,
    namespace: &MoleculeNamespace,
) -> Result<StereoAtomId, ParseError> {
    let site = participants.site.resolve(namespace)?;
    let ligands = resolve_ligands(participants.ligands, namespace)?;
    namespace
        .stereo_atom_by_participants(site, &ligands)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "stereo-atom",
            value: format!("site {}", site.index()),
        })
}

fn resolve_stereo_bond_structural(
    participants: StereoBondParticipants,
    namespace: &MoleculeNamespace,
) -> Result<StereoBondId, ParseError> {
    let site = participants.site.resolve(namespace)?;
    let ligands = resolve_ligands(participants.ligands, namespace)?;
    namespace
        .stereo_bond_by_participants(site, &ligands)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "stereo-bond",
            value: format!("site {}", site.index()),
        })
}

/// One ligand of a stereo element: an atom ref tagged with its kind
/// (`Atom` for a plain `<atom-ref>`, `ImplicitHydrogen` for `[:h <ref>]`,
/// `LonePair` for `[:lp <ref>]`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct StereoLigandRef {
    pub(crate) kind: StereoLigandKind,
    pub(crate) atom: AtomRef,
}

fn stereo_ligand_kind(tag: &str) -> Result<StereoLigandKind, DeError> {
    match tag {
        "h" => Ok(StereoLigandKind::ImplicitHydrogen),
        "lp" => Ok(StereoLigandKind::LonePair),
        other => Err(DeError::Custom(format!(
            "unknown stereo ligand tag :{other}"
        ))),
    }
}

pub(super) fn parse_stereo_ligand(edn: &Edn<'_>) -> Result<StereoLigandRef, DeError> {
    match edn {
        Edn::Vector(v) if v.len() == 2 => {
            let Edn::Keyword(tag) = &v[0] else {
                return Err(DeError::TypeMismatch {
                    expected: "ligand tag keyword",
                    got: v[0].kind(),
                    path: vec!["stereo-ligand".into()],
                });
            };
            Ok(StereoLigandRef {
                kind: stereo_ligand_kind(tag.name())?,
                atom: AtomRef::from_edn(&v[1])?,
            })
        }
        _ => Ok(StereoLigandRef {
            kind: StereoLigandKind::Atom,
            atom: AtomRef::from_edn(edn)?,
        }),
    }
}

pub(super) fn read_stereo_ligand(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoLigandRef, EdnError> {
    if de.peek_byte()?.ok_or_else(eof_err)? == b'[' {
        de.consume_byte(b'[')?;
        let kind = stereo_ligand_kind(de.read_keyword_name()?.as_ref())?;
        let atom = read_atom_ref(de)?;
        de.consume_byte(b']')?;
        Ok(StereoLigandRef { kind, atom })
    } else {
        Ok(StereoLigandRef {
            kind: StereoLigandKind::Atom,
            atom: read_atom_ref(de)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnKeyword};

    use super::*;

    #[fixture]
    fn meta_with_atom_id() -> MoleculeMetadata {
        MoleculeMetadata::new().with_atom_id(AtomId(2), "c1")
    }

    #[rstest]
    #[case::int(Edn::Int(3), AtomRef::Index(3))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("c1".into())), AtomRef::Id("c1".into()))]
    fn test_atom_ref_from_edn(#[case] input: Edn<'static>, #[case] expected: AtomRef) {
        assert_eq!(AtomRef::from_edn(&input).unwrap(), expected);
    }

    #[rstest]
    fn test_atom_ref_from_edn_rejects_other_kinds() {
        let err = AtomRef::from_edn(&Edn::Str("x".into())).unwrap_err();
        assert!(matches!(
            err,
            DeError::TypeMismatch {
                expected: "atom ref (int or keyword)",
                ..
            }
        ));
    }

    #[rstest]
    #[case::index(AtomRef::Index(5), Edn::Int(5))]
    #[case::id(AtomRef::Id("c1".into()), Edn::Keyword(EdnKeyword::owned("c1".into())))]
    fn test_atom_ref_to_edn(#[case] input: AtomRef, #[case] expected: Edn<'static>) {
        assert_eq!(input.to_edn(), expected);
    }

    #[rstest]
    #[case::int("3", AtomRef::Index(3))]
    #[case::keyword(":c1", AtomRef::Id("c1".into()))]
    fn test_atom_ref_roundtrip_edn_string(#[case] input: &str, #[case] expected: AtomRef) {
        let tree = read_string(input).unwrap();
        let parsed = AtomRef::from_edn(&tree).unwrap();
        assert_eq!(parsed, expected);
        let rendered = parsed.to_edn();
        let reparsed = AtomRef::from_edn(&rendered).unwrap();
        assert_eq!(reparsed, expected);
    }

    #[rstest]
    fn test_atom_ref_from_ast_uses_id_when_present(meta_with_atom_id: MoleculeMetadata) {
        let r = AtomRef::from_ast(AtomId(2), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Id("c1".into()));
    }

    #[rstest]
    fn test_atom_ref_from_ast_falls_back_to_index_without_id(meta_with_atom_id: MoleculeMetadata) {
        let r = AtomRef::from_ast(AtomId(4), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Index(4));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_id(meta_with_atom_id: MoleculeMetadata) {
        let id = AtomRef::Id("c1".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(id, AtomId(2));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_index(meta_with_atom_id: MoleculeMetadata) {
        let id = AtomRef::Index(3).into_ast(5, &meta_with_atom_id).unwrap();
        assert_eq!(id, AtomId(3));
    }

    #[rstest]
    fn test_atom_ref_into_ast_out_of_range_index(meta_with_atom_id: MoleculeMetadata) {
        let err = AtomRef::Index(9)
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "9".into(),
            }
        );
    }

    #[rstest]
    fn test_atom_ref_into_ast_unknown_id(meta_with_atom_id: MoleculeMetadata) {
        let err = AtomRef::Id("nope".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "nope".into(),
            }
        );
    }

    #[rstest]
    fn test_bond_ref_from_edn_structural() {
        let tree = read_string("{:atoms [0 1]}").unwrap();
        assert_eq!(
            BondRef::from_edn(&tree).unwrap(),
            BondRef::Structural([AtomRef::Index(0), AtomRef::Index(1)])
        );
    }

    #[rstest]
    #[case::type_key("{:atoms [0 1] :type \"c-c\"}")]
    #[case::id_key("{:atoms [0 1] :id :b1}")]
    fn test_bond_ref_from_edn_structural_error(#[case] input: &str) {
        let tree = read_string(input).unwrap();
        let DeError::Custom(msg) = BondRef::from_edn(&tree).unwrap_err() else {
            panic!("expected a Custom error");
        };
        assert_eq!(msg, "bond structural ref must not carry :type or :id");
    }

    #[rstest]
    fn test_bond_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        ns.register_atom(None);
        ns.register_atom(None);
        ns.register_bond(None, AtomId(0), AtomId(1));
        // Endpoint order is immaterial.
        let r = BondRef::Structural([AtomRef::Index(1), AtomRef::Index(0)]);
        assert_eq!(r.resolve(&ns).unwrap(), BondId(0));
    }

    #[rstest]
    #[case::no_matching_bond([AtomRef::Index(0), AtomRef::Index(1)], "bond", "[0 1]")]
    #[case::unknown_atom([AtomRef::Index(0), AtomRef::Index(5)], "atom", "5")]
    fn test_bond_ref_resolve_structural_error(
        #[case] atoms: [AtomRef; 2],
        #[case] kind: &'static str,
        #[case] value: &str,
    ) {
        let mut ns = MoleculeNamespace::default();
        ns.register_atom(None);
        ns.register_atom(None);
        assert_eq!(
            BondRef::Structural(atoms).resolve(&ns).unwrap_err(),
            ParseError::InvalidRef {
                kind,
                value: value.into(),
            }
        );
    }

    #[rstest]
    fn test_dative_bond_ref_from_edn_structural() {
        let tree = read_string("{:donors [1 2] :acceptor 0}").unwrap();
        assert_eq!(
            DativeBondRef::from_edn(&tree).unwrap(),
            DativeBondRef::Structural(DativeBondParticipants {
                donors: vec![AtomRef::Index(1), AtomRef::Index(2)],
                acceptor: AtomRef::Index(0),
            })
        );
    }

    #[rstest]
    fn test_dative_bond_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..3 {
            ns.register_atom(None);
        }
        ns.register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0));
        let r = DativeBondRef::Structural(DativeBondParticipants {
            donors: vec![AtomRef::Index(2), AtomRef::Index(1)],
            acceptor: AtomRef::Index(0),
        });
        assert_eq!(r.resolve(&ns).unwrap(), DativeBondId(0));
    }

    #[rstest]
    fn test_aromatic_system_ref_from_edn_structural() {
        let tree = read_string("{:atoms [0 1 :c2]}").unwrap();
        assert_eq!(
            AromaticSystemRef::from_edn(&tree).unwrap(),
            AromaticSystemRef::Structural(vec![
                AtomRef::Index(0),
                AtomRef::Index(1),
                AtomRef::Id("c2".into()),
            ])
        );
    }

    #[rstest]
    fn test_aromatic_system_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..3 {
            ns.register_atom(None);
        }
        ns.register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)]);
        // Atom order is immaterial.
        let r = AromaticSystemRef::Structural(vec![
            AtomRef::Index(0),
            AtomRef::Index(1),
            AtomRef::Index(2),
        ]);
        assert_eq!(r.resolve(&ns).unwrap(), AromaticSystemId(0));
    }

    #[rstest]
    fn test_multicenter_bond_ref_from_edn_structural() {
        let tree = read_string("{:atoms [0 1 2]}").unwrap();
        assert_eq!(
            MulticenterBondRef::from_edn(&tree).unwrap(),
            MulticenterBondRef::Structural(vec![
                AtomRef::Index(0),
                AtomRef::Index(1),
                AtomRef::Index(2),
            ])
        );
    }

    #[rstest]
    fn test_multicenter_bond_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..3 {
            ns.register_atom(None);
        }
        ns.register_multicenter_bond(None, &[AtomId(0), AtomId(1), AtomId(2)]);
        let r = MulticenterBondRef::Structural(vec![
            AtomRef::Index(2),
            AtomRef::Index(1),
            AtomRef::Index(0),
        ]);
        assert_eq!(r.resolve(&ns).unwrap(), MulticenterBondId(0));
    }

    #[rstest]
    fn test_noncovalent_bond_ref_from_edn_structural() {
        let tree = read_string("{:atoms [3 1]}").unwrap();
        assert_eq!(
            NoncovalentBondRef::from_edn(&tree).unwrap(),
            NoncovalentBondRef::Structural([AtomRef::Index(3), AtomRef::Index(1)])
        );
    }

    #[rstest]
    fn test_noncovalent_bond_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..4 {
            ns.register_atom(None);
        }
        ns.register_noncovalent_bond(None, AtomId(3), AtomId(1));
        let r = NoncovalentBondRef::Structural([AtomRef::Index(1), AtomRef::Index(3)]);
        assert_eq!(r.resolve(&ns).unwrap(), NoncovalentBondId(0));
    }

    #[rstest]
    fn test_stereo_atom_ref_from_edn_structural() {
        let tree = read_string("{:site 0 :ligands [1 2 [:h 3]]}").unwrap();
        assert_eq!(
            StereoAtomRef::from_edn(&tree).unwrap(),
            StereoAtomRef::Structural(StereoAtomParticipants {
                site: AtomRef::Index(0),
                ligands: vec![
                    StereoLigandRef {
                        kind: StereoLigandKind::Atom,
                        atom: AtomRef::Index(1)
                    },
                    StereoLigandRef {
                        kind: StereoLigandKind::Atom,
                        atom: AtomRef::Index(2)
                    },
                    StereoLigandRef {
                        kind: StereoLigandKind::ImplicitHydrogen,
                        atom: AtomRef::Index(3)
                    },
                ],
            })
        );
    }

    #[rstest]
    fn test_stereo_atom_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..5 {
            ns.register_atom(None);
        }
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        ns.register_stereo_atom(None, AtomId(4), &ligands);
        // Ligand frame order is immaterial.
        let r = StereoAtomRef::Structural(StereoAtomParticipants {
            site: AtomRef::Index(4),
            ligands: vec![
                StereoLigandRef {
                    kind: StereoLigandKind::Atom,
                    atom: AtomRef::Index(2),
                },
                StereoLigandRef {
                    kind: StereoLigandKind::Atom,
                    atom: AtomRef::Index(1),
                },
            ],
        });
        assert_eq!(r.resolve(&ns).unwrap(), StereoAtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_ref_resolve_structural_error() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..5 {
            ns.register_atom(None);
        }
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        ns.register_stereo_atom(None, AtomId(4), &ligands);
        // Same site, wrong ligand set.
        let r = StereoAtomRef::Structural(StereoAtomParticipants {
            site: AtomRef::Index(4),
            ligands: vec![StereoLigandRef {
                kind: StereoLigandKind::Atom,
                atom: AtomRef::Index(1),
            }],
        });
        assert_eq!(
            r.resolve(&ns).unwrap_err(),
            ParseError::InvalidRef {
                kind: "stereo-atom",
                value: "site 4".into(),
            }
        );
    }

    #[rstest]
    fn test_stereo_bond_ref_from_edn_structural() {
        let tree = read_string("{:site 0 :ligands [1 [:lp 2]]}").unwrap();
        assert_eq!(
            StereoBondRef::from_edn(&tree).unwrap(),
            StereoBondRef::Structural(StereoBondParticipants {
                site: BondRef::Index(0),
                ligands: vec![
                    StereoLigandRef {
                        kind: StereoLigandKind::Atom,
                        atom: AtomRef::Index(1)
                    },
                    StereoLigandRef {
                        kind: StereoLigandKind::LonePair,
                        atom: AtomRef::Index(2)
                    },
                ],
            })
        );
    }

    #[rstest]
    fn test_stereo_bond_ref_resolve_structural() {
        let mut ns = MoleculeNamespace::default();
        for _ in 0..4 {
            ns.register_atom(None);
        }
        ns.register_bond(None, AtomId(0), AtomId(1));
        let ligands = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        ns.register_stereo_bond(None, BondId(0), &ligands);
        let r = StereoBondRef::Structural(StereoBondParticipants {
            site: BondRef::Index(0),
            ligands: vec![StereoLigandRef {
                kind: StereoLigandKind::Atom,
                atom: AtomRef::Index(3),
            }],
        });
        assert_eq!(r.resolve(&ns).unwrap(), StereoBondId(0));
    }
}
