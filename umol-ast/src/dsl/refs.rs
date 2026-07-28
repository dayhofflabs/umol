//! Surface-level entity references. Each ref is a positional index (`Edn::Int`),
//! a keyword (`Edn::Keyword`), or — for non-atom entities — a
//! *structural* form (`Edn::Map`) naming the entity by its constituent atoms /
//! bonds. `resolve` turns a ref into an AST id against any parse-time `Namespace`
//! (count for index bounds, `find_by_keyword` for keyword references,
//! `find_by_participants` for the structural form); `denote` renders an id back
//! to a ref against any `Metadata` view.

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::edn_utils::{
    atoms_pair, atoms_vec, eof_err, parse_vec, read_map, read_vec, required_key, two_atom_refs,
};
use super::error::ParseError;
use super::metadata::Metadata;
use super::namespace::Namespace;
use crate::ast::entity::Entity;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::{StereoLigand, StereoLigandKind};

macro_rules! define_ref {
    ($name:ident, $id:ident, $variant:ident, $kind:literal, $reader:ident,
        $count:ident, $find_by_keyword:ident
        $(, structural = $payload:ty, $parse_structural:ident, $read_structural:ident, $resolve_structural:ident)?) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Keyword(String),
            $( Structural($payload), )?
        }

        impl $name {
            /// Render a numerical AST id back to a ref: its keyword from `metadata` if one is
            /// recorded for this id, else the bare index. This is the `id → ref` inverse of
            /// `resolve` over the entity-keyword bijection.
            pub fn denote<M: Metadata>(id: $id, metadata: &M) -> Self {
                if let Some(name) = metadata.keyword(Entity::$variant(id)) {
                    Self::Keyword(name.to_string())
                } else {
                    Self::Index(id.index())
                }
            }

            /// Resolve this ref to an AST id against the parse-time `namespace`
            /// (the source of truth: count for index bounds, keyword lookup, and
            /// participant lookup for the structural form).
            pub fn resolve<N: Namespace>(self, namespace: &N) -> Result<$id, ParseError> {
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
                    Self::Keyword(keyword) => {
                        namespace
                            .$find_by_keyword(&keyword)
                            .ok_or(ParseError::InvalidRef {
                                kind: $kind,
                                value: keyword,
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
                    Edn::Keyword(k) => Ok(Self::Keyword(k.name().to_string())),
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
                    Self::Keyword(name) => Edn::Keyword(EdnKeyword::owned(name.clone())),
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
                b':' => Ok($name::Keyword(de.read_keyword_name()?.into_owned())),
                $( b'{' => Ok($name::Structural($read_structural(de)?)), )?
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
    Atom,
    "atom",
    read_atom_ref,
    atom_count,
    find_atom_by_keyword
);
define_ref!(
    BondRef,
    BondId,
    Bond,
    "bond",
    read_bond_ref,
    bond_count,
    find_bond_by_keyword,
    structural = [AtomRef; 2],
    parse_bond_structural,
    read_bond_structural,
    resolve_bond_structural
);
define_ref!(
    DativeBondRef,
    DativeBondId,
    DativeBond,
    "dative-bond",
    read_dative_bond_ref,
    dative_bond_count,
    find_dative_bond_by_keyword,
    structural = DativeBondParticipants,
    parse_dative_structural,
    read_dative_structural,
    resolve_dative_structural
);
define_ref!(
    AromaticSystemRef, AromaticSystemId, AromaticSystem, "aromatic-system",
    read_aromatic_system_ref, aromatic_system_count, find_aromatic_system_by_keyword,
    structural = Vec<AtomRef>, parse_aromatic_structural, read_aromatic_structural, resolve_aromatic_structural
);
define_ref!(
    MulticenterBondRef, MulticenterBondId, MulticenterBond, "multicenter-bond",
    read_multicenter_bond_ref, multicenter_bond_count, find_multicenter_bond_by_keyword,
    structural = Vec<AtomRef>, parse_multicenter_structural, read_multicenter_structural, resolve_multicenter_structural
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondId,
    NoncovalentBond,
    "noncovalent-bond",
    read_noncovalent_bond_ref,
    noncovalent_bond_count,
    find_noncovalent_bond_by_keyword,
    structural = [AtomRef; 2],
    parse_noncovalent_structural,
    read_noncovalent_structural,
    resolve_noncovalent_structural
);
define_ref!(
    StereoAtomRef,
    StereoAtomId,
    StereoAtom,
    "stereo-atom",
    read_stereo_atom_ref,
    stereo_atom_count,
    find_stereo_atom_by_keyword,
    structural = StereoAtomParticipants,
    parse_stereo_atom_structural,
    read_stereo_atom_structural,
    resolve_stereo_atom_structural
);
define_ref!(
    StereoBondRef,
    StereoBondId,
    StereoBond,
    "stereo-bond",
    read_stereo_bond_ref,
    stereo_bond_count,
    find_stereo_bond_by_keyword,
    structural = StereoBondParticipants,
    parse_stereo_bond_structural,
    read_stereo_bond_structural,
    resolve_stereo_bond_structural
);

/// The constituent atoms of a dative bond named structurally (`{:donors [..] :acceptor a}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DativeBondParticipants {
    pub donors: Vec<AtomRef>,
    pub acceptor: AtomRef,
}

/// The site + ligand frame of a stereo atom named structurally (`{:site a :ligands [..]}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoAtomParticipants {
    pub site: AtomRef,
    pub ligands: Vec<StereoLigandRef>,
}

/// The site bond + ligand frame of a stereo bond named structurally
/// (`{:site bond-ref :ligands [..]}`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoBondParticipants {
    pub site: BondRef,
    pub ligands: Vec<StereoLigandRef>,
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

// The streaming counterparts of the `parse_*_structural` tree readers, over the deserializer cursor
// (they must not delegate to the tree path). Each reads the structural map's keys in order, rejecting
// an entity map's `:type` / `:id` (the same guard as `from_edn`).

fn reject_structural_key(key: &str, context: &'static str) -> EdnError {
    if key == "type" || key == "id" {
        DeError::Custom(format!("{context} must not carry :type or :id")).into()
    } else {
        DeError::Custom(format!("unexpected key :{key} in {context}")).into()
    }
}

fn missing_structural_key(key: &'static str, context: &'static str) -> EdnError {
    DeError::MissingField {
        key: key.to_string(),
        path: vec![context.into()],
    }
    .into()
}

/// The `:atoms` key of a structural map as a vector of atom refs (streaming).
fn read_structural_atoms(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<Vec<AtomRef>, EdnError> {
    let mut atoms: Option<Vec<AtomRef>> = None;
    read_map(de, |de, key| match key {
        "atoms" => {
            atoms = Some(read_vec(de, read_atom_ref)?);
            Ok(())
        }
        other => Err(reject_structural_key(other, context)),
    })?;
    atoms.ok_or_else(|| missing_structural_key("atoms", context))
}

fn read_bond_structural(de: &mut EdnStreamDeserializer<'_>) -> Result<[AtomRef; 2], EdnError> {
    let atoms = read_structural_atoms(de, "bond structural ref")?;
    two_atom_refs(atoms, "bond structural ref").map_err(Into::into)
}

fn read_noncovalent_structural(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<[AtomRef; 2], EdnError> {
    let atoms = read_structural_atoms(de, "noncovalent-bond structural ref")?;
    two_atom_refs(atoms, "noncovalent-bond structural ref").map_err(Into::into)
}

fn read_aromatic_structural(de: &mut EdnStreamDeserializer<'_>) -> Result<Vec<AtomRef>, EdnError> {
    read_structural_atoms(de, "aromatic-system structural ref")
}

fn read_multicenter_structural(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<Vec<AtomRef>, EdnError> {
    read_structural_atoms(de, "multicenter-bond structural ref")
}

fn read_dative_structural(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DativeBondParticipants, EdnError> {
    let context = "dative-bond structural ref";
    let mut donors: Option<Vec<AtomRef>> = None;
    let mut acceptor: Option<AtomRef> = None;
    read_map(de, |de, key| match key {
        "donors" => {
            donors = Some(read_vec(de, read_atom_ref)?);
            Ok(())
        }
        "acceptor" => {
            acceptor = Some(read_atom_ref(de)?);
            Ok(())
        }
        other => Err(reject_structural_key(other, context)),
    })?;
    Ok(DativeBondParticipants {
        donors: donors.ok_or_else(|| missing_structural_key("donors", context))?,
        acceptor: acceptor.ok_or_else(|| missing_structural_key("acceptor", context))?,
    })
}

fn read_stereo_atom_structural(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoAtomParticipants, EdnError> {
    let context = "stereo-atom structural ref";
    let mut site: Option<AtomRef> = None;
    let mut ligands: Option<Vec<StereoLigandRef>> = None;
    read_map(de, |de, key| match key {
        "site" => {
            site = Some(read_atom_ref(de)?);
            Ok(())
        }
        "ligands" => {
            ligands = Some(read_vec(de, read_stereo_ligand)?);
            Ok(())
        }
        other => Err(reject_structural_key(other, context)),
    })?;
    Ok(StereoAtomParticipants {
        site: site.ok_or_else(|| missing_structural_key("site", context))?,
        ligands: ligands.ok_or_else(|| missing_structural_key("ligands", context))?,
    })
}

fn read_stereo_bond_structural(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoBondParticipants, EdnError> {
    let context = "stereo-bond structural ref";
    let mut site: Option<BondRef> = None;
    let mut ligands: Option<Vec<StereoLigandRef>> = None;
    read_map(de, |de, key| match key {
        "site" => {
            site = Some(read_bond_ref(de)?);
            Ok(())
        }
        "ligands" => {
            ligands = Some(read_vec(de, read_stereo_ligand)?);
            Ok(())
        }
        other => Err(reject_structural_key(other, context)),
    })?;
    Ok(StereoBondParticipants {
        site: site.ok_or_else(|| missing_structural_key("site", context))?,
        ligands: ligands.ok_or_else(|| missing_structural_key("ligands", context))?,
    })
}

/// Resolve a vector of atom refs against the namespace, preserving order.
fn resolve_atom_refs<N: Namespace>(
    refs: Vec<AtomRef>,
    namespace: &N,
) -> Result<Vec<AtomId>, ParseError> {
    refs.into_iter().map(|r| r.resolve(namespace)).collect()
}

/// Resolve a stereo ligand frame (each ligand keeps its kind) against the namespace.
fn resolve_ligands<N: Namespace>(
    ligands: Vec<StereoLigandRef>,
    namespace: &N,
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

fn resolve_bond_structural<N: Namespace>(
    atoms: [AtomRef; 2],
    namespace: &N,
) -> Result<BondId, ParseError> {
    let [a, b] = atoms;
    let a = a.resolve(namespace)?;
    let b = b.resolve(namespace)?;
    namespace
        .find_bond_by_participants(a, b)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "bond",
            value: format_atom_ids(&[a, b]),
        })
}

fn resolve_noncovalent_structural<N: Namespace>(
    atoms: [AtomRef; 2],
    namespace: &N,
) -> Result<NoncovalentBondId, ParseError> {
    let [a, b] = atoms;
    let a = a.resolve(namespace)?;
    let b = b.resolve(namespace)?;
    namespace
        .find_noncovalent_bond_by_participants(a, b)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "noncovalent-bond",
            value: format_atom_ids(&[a, b]),
        })
}

fn resolve_aromatic_structural<N: Namespace>(
    atoms: Vec<AtomRef>,
    namespace: &N,
) -> Result<AromaticSystemId, ParseError> {
    let atoms = resolve_atom_refs(atoms, namespace)?;
    namespace
        .find_aromatic_system_by_participants(&atoms)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "aromatic-system",
            value: format_atom_ids(&atoms),
        })
}

fn resolve_multicenter_structural<N: Namespace>(
    atoms: Vec<AtomRef>,
    namespace: &N,
) -> Result<MulticenterBondId, ParseError> {
    let atoms = resolve_atom_refs(atoms, namespace)?;
    namespace
        .find_multicenter_bond_by_participants(&atoms)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "multicenter-bond",
            value: format_atom_ids(&atoms),
        })
}

fn resolve_dative_structural<N: Namespace>(
    participants: DativeBondParticipants,
    namespace: &N,
) -> Result<DativeBondId, ParseError> {
    let donors = resolve_atom_refs(participants.donors, namespace)?;
    let acceptor = participants.acceptor.resolve(namespace)?;
    namespace
        .find_dative_bond_by_participants(&donors, acceptor)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "dative-bond",
            value: format_atom_ids(&donors),
        })
}

fn resolve_stereo_atom_structural<N: Namespace>(
    participants: StereoAtomParticipants,
    namespace: &N,
) -> Result<StereoAtomId, ParseError> {
    let site = participants.site.resolve(namespace)?;
    let ligands = resolve_ligands(participants.ligands, namespace)?;
    namespace
        .find_stereo_atom_by_participants(site, &ligands)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "stereo-atom",
            value: format!("site {}", site.index()),
        })
}

fn resolve_stereo_bond_structural<N: Namespace>(
    participants: StereoBondParticipants,
    namespace: &N,
) -> Result<StereoBondId, ParseError> {
    let site = participants.site.resolve(namespace)?;
    let ligands = resolve_ligands(participants.ligands, namespace)?;
    namespace
        .find_stereo_bond_by_participants(site, &ligands)
        .ok_or_else(|| ParseError::InvalidRef {
            kind: "stereo-bond",
            value: format!("site {}", site.index()),
        })
}

/// One ligand of a stereo element: an atom ref tagged with its kind
/// (`Atom` for a plain `<atom-ref>`, `ImplicitHydrogen` for `[:h <ref>]`,
/// `LonePair` for `[:lp <ref>]`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoLigandRef {
    pub kind: StereoLigandKind,
    pub atom: AtomRef,
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

    use super::super::metadata::MoleculeMetadata;
    use super::super::namespace::MoleculeContext;
    use super::*;

    #[fixture]
    fn meta_with_atom_keyword() -> MoleculeMetadata {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(2)), "c1").unwrap();
        metadata
    }

    #[fixture]
    fn namespace_with_atom_keyword() -> MoleculeContext {
        let mut context = MoleculeContext::default();
        for i in 0..5 {
            context
                .register_atom((i == 2).then(|| "c1".to_string()))
                .unwrap();
        }
        context
    }

    #[rstest]
    #[case::int(Edn::Int(3), AtomRef::Index(3))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("c1".into())), AtomRef::Keyword("c1".into()))]
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
    #[case::keyword(
        AtomRef::Keyword("c1".into()),
        Edn::Keyword(EdnKeyword::owned("c1".into()))
    )]
    fn test_atom_ref_to_edn(#[case] input: AtomRef, #[case] expected: Edn<'static>) {
        assert_eq!(input.to_edn(), expected);
    }

    #[rstest]
    #[case::int("3", AtomRef::Index(3))]
    #[case::keyword(":c1", AtomRef::Keyword("c1".into()))]
    fn test_atom_ref_roundtrip_edn_string(#[case] input: &str, #[case] expected: AtomRef) {
        let tree = read_string(input).unwrap();
        let parsed = AtomRef::from_edn(&tree).unwrap();
        assert_eq!(parsed, expected);
        let rendered = parsed.to_edn();
        let reparsed = AtomRef::from_edn(&rendered).unwrap();
        assert_eq!(reparsed, expected);
    }

    #[rstest]
    #[case::keyword_present(AtomId(2), AtomRef::Keyword("c1".into()))]
    #[case::no_keyword(AtomId(4), AtomRef::Index(4))]
    fn test_atom_ref_denote(
        meta_with_atom_keyword: MoleculeMetadata,
        #[case] id: AtomId,
        #[case] expected: AtomRef,
    ) {
        assert_eq!(AtomRef::denote(id, &meta_with_atom_keyword), expected);
    }

    #[rstest]
    #[case::keyword(AtomRef::Keyword("c1".into()), AtomId(2))]
    #[case::index(AtomRef::Index(3), AtomId(3))]
    fn test_atom_ref_resolve(
        namespace_with_atom_keyword: MoleculeContext,
        #[case] r: AtomRef,
        #[case] expected: AtomId,
    ) {
        assert_eq!(r.resolve(&namespace_with_atom_keyword).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_index(AtomRef::Index(9), "9")]
    #[case::unknown_keyword(AtomRef::Keyword("nope".into()), "nope")]
    fn test_atom_ref_resolve_error(
        namespace_with_atom_keyword: MoleculeContext,
        #[case] r: AtomRef,
        #[case] value: &str,
    ) {
        assert_eq!(
            r.resolve(&namespace_with_atom_keyword).unwrap_err(),
            ParseError::InvalidRef {
                kind: "atom",
                value: value.into(),
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
        let mut context = MoleculeContext::default();
        context.register_atom(None).unwrap();
        context.register_atom(None).unwrap();
        context.register_bond(None, AtomId(0), AtomId(1)).unwrap();
        // Endpoint order is immaterial.
        let r = BondRef::Structural([AtomRef::Index(1), AtomRef::Index(0)]);
        assert_eq!(r.resolve(&context).unwrap(), BondId(0));
    }

    #[rstest]
    #[case::no_matching_bond([AtomRef::Index(0), AtomRef::Index(1)], "bond", "[0 1]")]
    #[case::unknown_atom([AtomRef::Index(0), AtomRef::Index(5)], "atom", "5")]
    fn test_bond_ref_resolve_structural_error(
        #[case] atoms: [AtomRef; 2],
        #[case] kind: &'static str,
        #[case] value: &str,
    ) {
        let mut context = MoleculeContext::default();
        context.register_atom(None).unwrap();
        context.register_atom(None).unwrap();
        assert_eq!(
            BondRef::Structural(atoms).resolve(&context).unwrap_err(),
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
        let mut context = MoleculeContext::default();
        for _ in 0..3 {
            context.register_atom(None).unwrap();
        }
        context
            .register_dative_bond(None, &[AtomId(1), AtomId(2)], AtomId(0))
            .unwrap();
        let r = DativeBondRef::Structural(DativeBondParticipants {
            donors: vec![AtomRef::Index(2), AtomRef::Index(1)],
            acceptor: AtomRef::Index(0),
        });
        assert_eq!(r.resolve(&context).unwrap(), DativeBondId(0));
    }

    #[rstest]
    fn test_aromatic_system_ref_from_edn_structural() {
        let tree = read_string("{:atoms [0 1 :c2]}").unwrap();
        assert_eq!(
            AromaticSystemRef::from_edn(&tree).unwrap(),
            AromaticSystemRef::Structural(vec![
                AtomRef::Index(0),
                AtomRef::Index(1),
                AtomRef::Keyword("c2".into()),
            ])
        );
    }

    #[rstest]
    fn test_aromatic_system_ref_resolve_structural() {
        let mut context = MoleculeContext::default();
        for _ in 0..3 {
            context.register_atom(None).unwrap();
        }
        context
            .register_aromatic_system(None, &[AtomId(2), AtomId(0), AtomId(1)])
            .unwrap();
        // Atom order is immaterial.
        let r = AromaticSystemRef::Structural(vec![
            AtomRef::Index(0),
            AtomRef::Index(1),
            AtomRef::Index(2),
        ]);
        assert_eq!(r.resolve(&context).unwrap(), AromaticSystemId(0));
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
        let mut context = MoleculeContext::default();
        for _ in 0..3 {
            context.register_atom(None).unwrap();
        }
        context
            .register_multicenter_bond(None, &[AtomId(0), AtomId(1), AtomId(2)])
            .unwrap();
        let r = MulticenterBondRef::Structural(vec![
            AtomRef::Index(2),
            AtomRef::Index(1),
            AtomRef::Index(0),
        ]);
        assert_eq!(r.resolve(&context).unwrap(), MulticenterBondId(0));
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
        let mut context = MoleculeContext::default();
        for _ in 0..4 {
            context.register_atom(None).unwrap();
        }
        context
            .register_noncovalent_bond(None, AtomId(3), AtomId(1))
            .unwrap();
        let r = NoncovalentBondRef::Structural([AtomRef::Index(1), AtomRef::Index(3)]);
        assert_eq!(r.resolve(&context).unwrap(), NoncovalentBondId(0));
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
        let mut context = MoleculeContext::default();
        for _ in 0..5 {
            context.register_atom(None).unwrap();
        }
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        context
            .register_stereo_atom(None, AtomId(4), &ligands)
            .unwrap();
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
        assert_eq!(r.resolve(&context).unwrap(), StereoAtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_ref_resolve_structural_error() {
        let mut context = MoleculeContext::default();
        for _ in 0..5 {
            context.register_atom(None).unwrap();
        }
        let ligands = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ];
        context
            .register_stereo_atom(None, AtomId(4), &ligands)
            .unwrap();
        // Same site, wrong ligand set.
        let r = StereoAtomRef::Structural(StereoAtomParticipants {
            site: AtomRef::Index(4),
            ligands: vec![StereoLigandRef {
                kind: StereoLigandKind::Atom,
                atom: AtomRef::Index(1),
            }],
        });
        assert_eq!(
            r.resolve(&context).unwrap_err(),
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
        let mut context = MoleculeContext::default();
        for _ in 0..4 {
            context.register_atom(None).unwrap();
        }
        context.register_bond(None, AtomId(0), AtomId(1)).unwrap();
        let ligands = [StereoLigand::new(AtomId(3), StereoLigandKind::Atom)];
        context
            .register_stereo_bond(None, BondId(0), &ligands)
            .unwrap();
        let r = StereoBondRef::Structural(StereoBondParticipants {
            site: BondRef::Index(0),
            ligands: vec![StereoLigandRef {
                kind: StereoLigandKind::Atom,
                atom: AtomRef::Index(3),
            }],
        });
        assert_eq!(r.resolve(&context).unwrap(), StereoBondId(0));
    }
}
