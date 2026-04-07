//! Molecule map DSL: parser and AST

use std::borrow::Cow;
use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::{self, SerializeSeq, Serializer};
use serde::{Deserialize, Serialize};
use umol_data::SpinState;
use umol_edn::EdnKeyword;

use super::ast::DslAst;
use super::atom::{parse_atom_dsl, AtomAst};
use super::bond::BondAst;
use super::config::MoleculeDslConfig;
use super::error::ParseError;

/// Atom reference in bond endpoints.
///
/// `Index(n)` references an atom by positional index (serializes as integer).
/// `Tag(s)` references a named atom by tag (serializes as EDN keyword via `EdnKeyword`).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomRef {
    Index(usize),
    Tag(String),
}

impl AtomRef {
    pub fn name(&self) -> Cow<'_, str> {
        match self {
            AtomRef::Index(n) => Cow::Owned(n.to_string()),
            AtomRef::Tag(s) => Cow::Borrowed(s),
        }
    }
}

impl fmt::Display for AtomRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AtomRef::Index(n) => write!(f, "{n}"),
            AtomRef::Tag(s) => f.write_str(s),
        }
    }
}

impl From<usize> for AtomRef {
    fn from(n: usize) -> Self {
        AtomRef::Index(n)
    }
}

impl From<String> for AtomRef {
    fn from(s: String) -> Self {
        AtomRef::Tag(s)
    }
}

impl From<&str> for AtomRef {
    fn from(s: &str) -> Self {
        AtomRef::Tag(s.to_string())
    }
}

impl Serialize for AtomRef {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            AtomRef::Index(n) => serializer.serialize_u64(*n as u64),
            AtomRef::Tag(s) => EdnKeyword::new(s).serialize(serializer),
        }
    }
}

impl<'de> Deserialize<'de> for AtomRef {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(AtomRefVisitor)
    }
}

struct AtomRefVisitor;

impl<'de> Visitor<'de> for AtomRefVisitor {
    type Value = AtomRef;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an atom reference (integer or keyword/string tag)")
    }

    fn visit_i64<E: de::Error>(self, v: i64) -> Result<AtomRef, E> {
        let n = usize::try_from(v).map_err(|_| de::Error::custom("negative atom index"))?;
        Ok(AtomRef::Index(n))
    }

    fn visit_u64<E: de::Error>(self, v: u64) -> Result<AtomRef, E> {
        Ok(AtomRef::Index(v as usize))
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<AtomRef, E> {
        Ok(AtomRef::Tag(v.to_string()))
    }

    fn visit_string<E: de::Error>(self, v: String) -> Result<AtomRef, E> {
        Ok(AtomRef::Tag(v))
    }
}

/// Atom collection: a positional list of atoms with optional tags and aliases.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Atoms {
    pub entries: Vec<AtomAst>,
    pub tags: IndexMap<String, usize>,
    pub aliases: BiMap<String, AtomAst>,
}

impl Default for Atoms {
    fn default() -> Self {
        Self {
            entries: vec![],
            tags: IndexMap::new(),
            aliases: BiMap::new(),
        }
    }
}

impl Atoms {
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn indexed(entries: Vec<AtomAst>) -> Self {
        Self {
            entries,
            tags: IndexMap::new(),
            aliases: BiMap::new(),
        }
    }

    pub fn named(entries: IndexMap<String, AtomAst>) -> Self {
        let tags: IndexMap<String, usize> = entries
            .keys()
            .enumerate()
            .map(|(i, k)| (k.clone(), i))
            .collect();
        let atom_vec = entries.into_values().collect();
        Self {
            entries: atom_vec,
            tags,
            aliases: BiMap::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalizedBond {
    pub id: Option<String>,
    pub a: AtomRef,
    pub b: AtomRef,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DativeBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub donor: AtomRef,
    pub acceptor: AtomRef,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AromaticSystem {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub atoms: Vec<AtomRef>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MulticenterBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub atoms: Vec<AtomRef>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NoncovalentBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub a: AtomRef,
    pub b: AtomRef,
    pub bond: BondAst,
}

/// Parsed molecule map AST
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct MoleculeAst {
    pub atoms: Atoms,
    pub bonds: Vec<LocalizedBond>,
    pub dative_bonds: Vec<DativeBond>,
    pub aromatic_systems: Vec<AromaticSystem>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub noncovalent_bonds: Vec<NoncovalentBond>,
    pub charge: Option<i64>,
    pub spin: Option<SpinState>,
}

impl DslAst for MoleculeAst {
    type Config = MoleculeDslConfig;
}

/// DSL formatting metadata extracted during AST→Builder conversion.
/// Travels alongside MoleculeBuilder, NOT inside Molecule (which is a semantic object).
#[derive(Clone, Debug, Default)]
pub struct Metadata {
    pub atom_tags: IndexMap<usize, String>,
    pub atom_aliases: BiMap<String, AtomAst>,
}

impl Metadata {
    pub fn extract(atoms: &Atoms) -> Self {
        let atom_tags: IndexMap<usize, String> = atoms
            .tags
            .iter()
            .map(|(name, &pos)| (pos, name.clone()))
            .collect();
        Self {
            atom_tags,
            atom_aliases: atoms.aliases.clone(),
        }
    }

    pub fn apply(&self, atoms_entries: Vec<AtomAst>) -> Atoms {
        let tags: IndexMap<String, usize> = self
            .atom_tags
            .iter()
            .map(|(&pos, name)| (name.clone(), pos))
            .collect();
        Atoms {
            entries: atoms_entries,
            tags,
            aliases: self.atom_aliases.clone(),
        }
    }
}

// Serde: Atoms
impl Serialize for Atoms {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let tag_for_pos: IndexMap<usize, &str> = self
            .tags
            .iter()
            .map(|(name, &pos)| (pos, name.as_str()))
            .collect();

        let mut s = serializer.serialize_seq(Some(self.entries.len()))?;
        for (i, atom) in self.entries.iter().enumerate() {
            let alias_name = self.aliases.get_by_right(atom);
            if let Some(tag) = tag_for_pos.get(&i) {
                // Tagged entry: [:tag <def-or-alias>]
                if let Some(alias) = alias_name {
                    s.serialize_element(&(EdnKeyword::new(*tag), EdnKeyword::new(alias)))?;
                } else {
                    s.serialize_element(&(EdnKeyword::new(*tag), atom))?;
                }
            } else if let Some(alias) = alias_name {
                // Aliased entry (no tag): emit alias name as keyword
                s.serialize_element(&EdnKeyword::new(alias))?;
            } else {
                // Plain entry
                s.serialize_element(atom)?;
            }
        }
        s.end()
    }
}

impl<'de> Deserialize<'de> for Atoms {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(AtomsVisitor)
    }
}

struct AtomsVisitor;

impl<'de> Visitor<'de> for AtomsVisitor {
    type Value = Atoms;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a vector of atoms")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Atoms, A::Error> {
        let mut named = IndexMap::new();
        while let Some((label, atom)) = map.next_entry::<String, AtomAst>()? {
            named.insert(label, atom);
        }
        Ok(Atoms::named(named))
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Atoms, A::Error> {
        let mut atoms = Vec::new();
        while let Some(atom) = seq.next_element::<AtomAst>()? {
            atoms.push(atom);
        }
        Ok(Atoms::indexed(atoms))
    }

    fn visit_unit<E: de::Error>(self) -> Result<Atoms, E> {
        Ok(Atoms::default())
    }
}

// Serde: LocalizedBond
impl Serialize for LocalizedBond {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if self.id.is_none() {
            // Shorthand: [:a :b :bond-spec]
            let mut t = serializer.serialize_tuple(3)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.a)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.b)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.bond)?;
            ser::SerializeTuple::end(t)
        } else {
            // Map form: {:id :b1 :a :0 :b :1 :bond :single}
            let mut m = serializer.serialize_struct("LocalizedBond", 4)?;
            if let Some(id) = &self.id {
                ser::SerializeStruct::serialize_field(&mut m, "id", id)?;
            }
            ser::SerializeStruct::serialize_field(&mut m, "a", &self.a)?;
            ser::SerializeStruct::serialize_field(&mut m, "b", &self.b)?;
            ser::SerializeStruct::serialize_field(&mut m, "bond", &self.bond)?;
            ser::SerializeStruct::end(m)
        }
    }
}

impl<'de> Deserialize<'de> for LocalizedBond {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(LocalizedBondVisitor)
    }
}

struct LocalizedBondVisitor;

impl<'de> Visitor<'de> for LocalizedBondVisitor {
    type Value = LocalizedBond;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a bond vector [:a :b :spec] or map {:a :A :b :B :bond :spec}")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<LocalizedBond, A::Error> {
        let a: AtomRef = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &"3-element vector"))?;
        let b: AtomRef = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &"3-element vector"))?;
        let bond: BondAst = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(2, &"3-element vector"))?;
        Ok(LocalizedBond {
            id: None,
            a,
            b,
            bond,
        })
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<LocalizedBond, A::Error> {
        let mut id = None;
        let mut a = None;
        let mut b = None;
        let mut bond = None;
        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "id" => id = Some(map.next_value::<String>()?),
                "a" => a = Some(map.next_value::<AtomRef>()?),
                "b" => b = Some(map.next_value::<AtomRef>()?),
                "bond" => bond = Some(map.next_value::<BondAst>()?),
                _ => {
                    let _ = map.next_value::<de::IgnoredAny>()?;
                }
            }
        }
        Ok(LocalizedBond {
            id,
            a: a.ok_or_else(|| de::Error::missing_field("a"))?,
            b: b.ok_or_else(|| de::Error::missing_field("b"))?,
            bond: bond.ok_or_else(|| de::Error::missing_field("bond"))?,
        })
    }
}

// Serde: MoleculeAst
impl Serialize for MoleculeAst {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let has_aliases = !self.atoms.aliases.is_empty();
        let mut count = 2;
        if !self.dative_bonds.is_empty() {
            count += 1;
        }
        if !self.aromatic_systems.is_empty() {
            count += 1;
        }
        if !self.multicenter_bonds.is_empty() {
            count += 1;
        }
        if !self.noncovalent_bonds.is_empty() {
            count += 1;
        }
        if self.charge.is_some() {
            count += 1;
        }
        if self.spin.is_some() {
            count += 1;
        }
        if has_aliases {
            count += 1;
        }

        let mut m = serializer.serialize_struct("MoleculeAst", count)?;
        ser::SerializeStruct::serialize_field(&mut m, "atoms", &self.atoms)?;
        ser::SerializeStruct::serialize_field(&mut m, "bonds", &self.bonds)?;
        if !self.dative_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "dative", &self.dative_bonds)?;
        }
        if !self.aromatic_systems.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "aromatic", &self.aromatic_systems)?;
        }
        if !self.multicenter_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "multicenter", &self.multicenter_bonds)?;
        }
        if !self.noncovalent_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "noncovalent", &self.noncovalent_bonds)?;
        }
        if let Some(charge) = self.charge {
            ser::SerializeStruct::serialize_field(&mut m, "charge", &charge)?;
        }
        if let Some(spin) = &self.spin {
            ser::SerializeStruct::serialize_field(&mut m, "spin", &spin.to_string())?;
        }
        if has_aliases {
            let alias_vec: Vec<AliasEntry<'_>> = self
                .atoms
                .aliases
                .iter()
                .flat_map(|(name, atom)| {
                    [
                        AliasEntry::Name(EdnKeyword::new(name)),
                        AliasEntry::Def(atom),
                    ]
                })
                .collect();
            ser::SerializeStruct::serialize_field(&mut m, "atom-aliases", &alias_vec)?;
        }
        ser::SerializeStruct::end(m)
    }
}

enum AliasEntry<'a> {
    Name(EdnKeyword),
    Def(&'a AtomAst),
}

impl Serialize for AliasEntry<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            AliasEntry::Name(kw) => kw.serialize(serializer),
            AliasEntry::Def(atom) => atom.serialize(serializer),
        }
    }
}

impl<'de> Deserialize<'de> for MoleculeAst {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(MoleculeAstVisitor)
    }
}

struct MoleculeAstVisitor;

impl<'de> Visitor<'de> for MoleculeAstVisitor {
    type Value = MoleculeAst;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a molecule map with :atoms and :bonds")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<MoleculeAst, A::Error> {
        let mut atoms = None;
        let mut bonds = None;
        let mut dative_bonds = None;
        let mut aromatic_systems = None;
        let mut multicenter_bonds = None;
        let mut noncovalent_bonds = None;
        let mut charge = None;
        let mut spin = None;
        // Aliases are resolved inline during atom deserialization (not stored).
        // We collect them as raw strings here, then resolve after all fields are read.
        let mut aliases: Option<Vec<(String, String)>> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "atoms" => atoms = Some(map.next_value::<Atoms>()?),
                "bonds" => bonds = Some(map.next_value::<Vec<LocalizedBond>>()?),
                "dative" => dative_bonds = Some(map.next_value::<Vec<DativeBond>>()?),
                "aromatic" => aromatic_systems = Some(map.next_value::<Vec<AromaticSystem>>()?),
                "multicenter" => {
                    multicenter_bonds = Some(map.next_value::<Vec<MulticenterBond>>()?)
                }
                "noncovalent" => {
                    noncovalent_bonds = Some(map.next_value::<Vec<NoncovalentBond>>()?)
                }
                "charge" => charge = map.next_value::<Option<i64>>()?,
                "spin" => {
                    let s: Option<String> = map.next_value()?;
                    spin = match s {
                        Some(ref s) => Some(
                            SpinState::from_str(s).map_err(|e| de::Error::custom(e.to_string()))?,
                        ),
                        None => None,
                    };
                }
                "aliases" => aliases = Some(map.next_value::<Vec<(String, String)>>()?),
                _ => {
                    let _ = map.next_value::<de::IgnoredAny>()?;
                }
            }
        }

        let atoms = atoms.ok_or_else(|| de::Error::missing_field("atoms"))?;
        let bonds = bonds.ok_or_else(|| de::Error::missing_field("bonds"))?;
        let _ = aliases;

        Ok(MoleculeAst {
            atoms,
            bonds,
            dative_bonds: dative_bonds.unwrap_or_default(),
            aromatic_systems: aromatic_systems.unwrap_or_default(),
            multicenter_bonds: multicenter_bonds.unwrap_or_default(),
            noncovalent_bonds: noncovalent_bonds.unwrap_or_default(),
            charge,
            spin,
        })
    }
}

fn validate(ast: &MoleculeAst) -> Result<(), ParseError> {
    let mut atom_labels: HashSet<String> = (0..ast.atoms.entries.len())
        .map(|i| i.to_string())
        .collect();
    for tag in ast.atoms.tags.keys() {
        atom_labels.insert(tag.clone());
    }

    let mut seen_ids: HashSet<&str> = HashSet::new();
    for tag in ast.atoms.tags.keys() {
        seen_ids.insert(tag.as_str());
    }
    for id in ast
        .bonds
        .iter()
        .filter_map(|e| e.id.as_deref())
        .chain(ast.dative_bonds.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.aromatic_systems.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.multicenter_bonds.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.noncovalent_bonds.iter().filter_map(|e| e.id.as_deref()))
    {
        if !seen_ids.insert(id) {
            return Err(ParseError::DuplicateId(id.to_string()));
        }
    }

    let check = |atom_ref: &AtomRef| -> Result<(), ParseError> {
        let label = atom_ref.name();
        if atom_labels.contains(label.as_ref()) {
            Ok(())
        } else {
            Err(ParseError::InvalidAtomIndex(label.into_owned()))
        }
    };

    for e in &ast.bonds {
        check(&e.a)?;
        check(&e.b)?;
    }
    for e in &ast.dative_bonds {
        check(&e.donor)?;
        check(&e.acceptor)?;
    }
    for e in &ast.aromatic_systems {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.multicenter_bonds {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.noncovalent_bonds {
        check(&e.a)?;
        check(&e.b)?;
    }

    Ok(())
}

// MoleculeInput: streaming read adapter

/// Raw atom entry collected during streaming deserialization.
#[derive(Clone, Debug)]
enum AtomEntryInput {
    Str(String),
    Tagged(String, Box<AtomEntryInput>),
}

impl<'de> Deserialize<'de> for AtomEntryInput {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(AtomEntryInputVisitor)
    }
}

struct AtomEntryInputVisitor;

impl<'de> Visitor<'de> for AtomEntryInputVisitor {
    type Value = AtomEntryInput;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an atom definition string, alias keyword, or [:tag def] vector")
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<AtomEntryInput, E> {
        Ok(AtomEntryInput::Str(v.to_string()))
    }

    fn visit_string<E: de::Error>(self, v: String) -> Result<AtomEntryInput, E> {
        Ok(AtomEntryInput::Str(v))
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<AtomEntryInput, A::Error> {
        let tag: String = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &"2-element [tag def] vector"))?;
        let entry: AtomEntryInput = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &"2-element [tag def] vector"))?;
        Ok(AtomEntryInput::Tagged(tag, Box::new(entry)))
    }
}

/// Intermediate struct for streaming deserialization of the molecule DSL.
/// Collects raw strings and defers alias resolution to `into_ast()`.
struct MoleculeInput {
    atoms: Vec<AtomEntryInput>,
    bonds: Vec<LocalizedBond>,
    dative_bonds: Vec<DativeBond>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
    noncovalent_bonds: Vec<NoncovalentBond>,
    atom_aliases: Vec<String>,
    charge: Option<i64>,
    spin: Option<SpinState>,
}

impl<'de> Deserialize<'de> for MoleculeInput {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(MoleculeInputVisitor)
    }
}

struct MoleculeInputVisitor;

impl<'de> Visitor<'de> for MoleculeInputVisitor {
    type Value = MoleculeInput;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a molecule map with :atoms and :bonds")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<MoleculeInput, A::Error> {
        let mut atoms = None;
        let mut bonds = None;
        let mut dative_bonds = None;
        let mut aromatic_systems = None;
        let mut multicenter_bonds = None;
        let mut noncovalent_bonds = None;
        let mut charge = None;
        let mut spin = None;
        let mut atom_aliases = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "atoms" => atoms = Some(map.next_value::<Vec<AtomEntryInput>>()?),
                "bonds" => bonds = Some(map.next_value::<Vec<LocalizedBond>>()?),
                "dative" => dative_bonds = Some(map.next_value::<Vec<DativeBond>>()?),
                "aromatic" => aromatic_systems = Some(map.next_value::<Vec<AromaticSystem>>()?),
                "multicenter" => {
                    multicenter_bonds = Some(map.next_value::<Vec<MulticenterBond>>()?)
                }
                "noncovalent" => {
                    noncovalent_bonds = Some(map.next_value::<Vec<NoncovalentBond>>()?)
                }
                "charge" => charge = map.next_value::<Option<i64>>()?,
                "spin" => {
                    let s: Option<String> = map.next_value()?;
                    spin = match s {
                        Some(ref s) => Some(
                            SpinState::from_str(s).map_err(|e| de::Error::custom(e.to_string()))?,
                        ),
                        None => None,
                    };
                }
                "atom-aliases" | "aliases" => atom_aliases = Some(map.next_value::<Vec<String>>()?),
                _ => {
                    let _ = map.next_value::<de::IgnoredAny>()?;
                }
            }
        }

        let atoms = atoms.ok_or_else(|| de::Error::missing_field("atoms"))?;
        let bonds = bonds.ok_or_else(|| de::Error::missing_field("bonds"))?;

        Ok(MoleculeInput {
            atoms,
            bonds,
            dative_bonds: dative_bonds.unwrap_or_default(),
            aromatic_systems: aromatic_systems.unwrap_or_default(),
            multicenter_bonds: multicenter_bonds.unwrap_or_default(),
            noncovalent_bonds: noncovalent_bonds.unwrap_or_default(),
            atom_aliases: atom_aliases.unwrap_or_default(),
            charge,
            spin,
        })
    }
}

impl MoleculeInput {
    fn into_ast(self) -> Result<MoleculeAst, ParseError> {
        // 1. Build alias table from flat vector [name, def, name, def, ...]
        let raw_aliases = &self.atom_aliases;
        if raw_aliases.len() % 2 != 0 {
            return Err(ParseError::WrongFieldType {
                field: "atom-aliases".to_string(),
                expected: "flat vector of keyword/atom-spec pairs (even length)".to_string(),
            });
        }
        let mut alias_table: IndexMap<String, String> = IndexMap::new();
        for pair in raw_aliases.chunks(2) {
            let name = &pair[0];
            let def = &pair[1];
            if alias_table.contains_key(name) {
                return Err(ParseError::DuplicateId(name.clone()));
            }
            alias_table.insert(name.clone(), def.clone());
        }

        // 2. Process atom entries: resolve aliases, extract tags
        let mut entries = Vec::new();
        let mut tags: IndexMap<String, usize> = IndexMap::new();
        let mut aliases_bimap: bimap::BiMap<String, AtomAst> = bimap::BiMap::new();

        for entry in self.atoms {
            let pos = entries.len();
            let (tag, atom_str) = Self::resolve_entry(entry, &alias_table)?;

            let atom_ast = parse_atom_dsl(&atom_str)?;

            if let Some(tag_name) = tag {
                if tags.contains_key(&tag_name) {
                    return Err(ParseError::DuplicateId(tag_name));
                }
                tags.insert(tag_name, pos);
            }

            entries.push(atom_ast);
        }

        // 3. Build BiMap from alias definitions
        for (name, def) in &alias_table {
            let atom_ast = parse_atom_dsl(def)?;
            aliases_bimap.insert(name.clone(), atom_ast);
        }

        // 4. Validate tag/alias name disjointness
        for tag_name in tags.keys() {
            if alias_table.contains_key(tag_name) {
                return Err(ParseError::DuplicateId(tag_name.clone()));
            }
        }

        let ast = MoleculeAst {
            atoms: Atoms {
                entries,
                tags,
                aliases: aliases_bimap,
            },
            bonds: self.bonds,
            dative_bonds: self.dative_bonds,
            aromatic_systems: self.aromatic_systems,
            multicenter_bonds: self.multicenter_bonds,
            noncovalent_bonds: self.noncovalent_bonds,
            charge: self.charge,
            spin: self.spin,
        };

        validate(&ast)?;
        Ok(ast)
    }

    fn resolve_entry(
        entry: AtomEntryInput,
        alias_table: &IndexMap<String, String>,
    ) -> Result<(Option<String>, String), ParseError> {
        match entry {
            AtomEntryInput::Str(s) => {
                if let Some(def) = alias_table.get(&s) {
                    Ok((None, def.clone()))
                } else {
                    Ok((None, s))
                }
            }
            AtomEntryInput::Tagged(tag, inner) => {
                let (_, atom_str) = Self::resolve_entry(*inner, alias_table)?;
                Ok((Some(tag), atom_str))
            }
        }
    }
}

// FromStr / Display / parse_molecule_dsl

/// Parse a molecule AST from an EDN string using streaming deserialization.
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeAst, ParseError> {
    let mol_input: MoleculeInput =
        umol_edn::from_str(input).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    mol_input.into_ast()
}

impl FromStr for MoleculeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_molecule_dsl(s)
    }
}

impl fmt::Display for MoleculeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = umol_edn::to_string(self).map_err(|_| fmt::Error)?;
        f.write_str(&s)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{e, Element};

    use super::super::predicates::{ElementExpr, HydrogenExpr};
    use super::super::value::ValueAst;
    use super::super::ast::{FromAst, ToAst};
    use super::*;
    use crate::graph_ir::molecule_builder::MoleculeBuilder;

    #[rstest]
    #[case::empty(r#"{:atoms [] :bonds []}"#, MoleculeAst::default())]
    #[case::atom(r#"{:atoms ["C"] :bonds []}"#, MoleculeAst { atoms: Atoms::indexed(vec![AtomAst::from_element(e!(C))]), ..Default::default() })]
    #[case::atom_tagged(
        r#"{:atoms [[:C "C"]] :bonds []}"#,
        MoleculeAst {
            atoms: Atoms { entries: vec![AtomAst::from_element(e!(C))], tags: IndexMap::from([("C".to_string(), 0)]), aliases: bimap::BiMap::new() },
            ..Default::default()
        }
    )]
    #[case::bond(r#"{:atoms ["N" "N"] :bonds [[0 1 :triple]]}"#, MoleculeAst { atoms: Atoms::indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![LocalizedBond { id: None, a: AtomRef::Index(0), b: AtomRef::Index(1), bond: BondAst::from_order(3) }], ..Default::default() })]
    #[case::bond_with_tags(
        r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#,
        MoleculeAst {
            atoms: Atoms { entries: vec![AtomAst::from_element(e!(C)), AtomAst::from_element(e!(O))], tags: IndexMap::from([("C".to_string(), 0), ("O".to_string(), 1)]), aliases: bimap::BiMap::new() },
            bonds: vec![LocalizedBond { id: None, a: AtomRef::Tag("C".into()), b: AtomRef::Tag("O".into()), bond: BondAst::from_order(1) }],
            ..Default::default()
        }
    )]
    #[case::bond_id(r#"{:atoms ["H" "F"] :bonds [{:id :b1 :a 0 :b 1 :bond :single}]}"#,
        MoleculeAst { atoms: Atoms::indexed(vec![AtomAst::from_element(e!(H)), AtomAst::from_element(e!(F))]),
        bonds: vec![LocalizedBond { id: Some("b1".to_string()), a: AtomRef::Index(0), b: AtomRef::Index(1), bond: BondAst::from_order(1) }], ..Default::default() })]
    #[case::charge(
        r#"{:atoms [[:F "F#c-"]] :bonds [] :charge -1}"#,
        MoleculeAst {
            atoms: Atoms { entries: vec![AtomAst { element: ElementExpr::Lit(Element::F), isotope_mass: None, implicit_hydrogens: None, charge: Some(ValueAst::Lit(-1)), lone_pairs: None, unpaired_electrons: None,
                multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None }],
                tags: IndexMap::from([("F".to_string(), 0)]), aliases: bimap::BiMap::new() },
            charge: Some(-1), ..Default::default()
        }
    )]
    #[case::alias_indexed(r#"{:atoms [:ch] :bonds [] :aliases [:ch "C #h1"]}"#,
        MoleculeAst { atoms: Atoms { entries: vec![AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None }],
        tags: IndexMap::new(), aliases: bimap::BiMap::from_iter([("ch".to_string(), AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None })]) },
        ..Default::default() })]
    #[case::alias_reused(r#"{:atoms [:n :n] :bonds [[0 1 :single]] :aliases [:n "N"]}"#,
        MoleculeAst { atoms: Atoms { entries: vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))],
        tags: IndexMap::new(), aliases: bimap::BiMap::from_iter([("n".to_string(), AtomAst::from_element(e!(N)))]) },
        bonds: vec![LocalizedBond { id: None, a: AtomRef::Index(0), b: AtomRef::Index(1), bond: BondAst::from_order(1) }], ..Default::default() })]
    #[case::alias_tagged(
        r#"{:atoms [[:C :ch]] :bonds [] :aliases [:ch "C #h1"]}"#,
        MoleculeAst { atoms: Atoms { entries: vec![AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None }],
        tags: IndexMap::from([("C".to_string(), 0)]), aliases: bimap::BiMap::from_iter([("ch".to_string(), AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None })]) },
        ..Default::default() }
    )]
    fn test_parse_molecule_dsl(#[case] input: &str, #[case] expected: MoleculeAst) {
        let result = parse_molecule_dsl(input);
        assert!(
            result.is_ok(),
            "{input:?} should succeed, got {:?}",
            result.unwrap_err()
        );
        let ast = result.unwrap();
        assert_eq!(ast, expected);
    }

    #[test]
    fn test_parse_molecule_dsl_dative() {
        let ast = parse_molecule_dsl(
            r#"{:atoms [[:B "B #h3"] [:N "N #h3"]]
                :bonds []
                :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#,
        )
        .unwrap();
        assert_eq!(ast.dative_bonds.len(), 1);
        assert_eq!(ast.dative_bonds[0].donor, AtomRef::Tag("N".into()));
        assert_eq!(ast.dative_bonds[0].acceptor, AtomRef::Tag("B".into()));
    }

    #[test]
    fn test_parse_molecule_dsl_aromatic() {
        let ast = parse_molecule_dsl(
            r#"{:atoms [[:C1 :ch] [:C2 :ch] [:C3 :ch] [:C4 :ch] [:C5 :ch] [:C6 :ch]]
                :bonds [[:C1 :C2 :single] [:C2 :C3 :single] [:C3 :C4 :single] [:C4 :C5 :single] [:C5 :C6 :single] [:C6 :C1 :single]]
                :aromatic [{:id :ar1 :atoms [:C1 :C2 :C3 :C4 :C5 :C6]}]
                :aliases [:ch "C #h1 #v2 #a1"]}"#,
        )
        .unwrap();
        assert_eq!(ast.aromatic_systems.len(), 1);
        assert_eq!(ast.aromatic_systems[0].id, Some("ar1".to_string()));
        let refs: Vec<AtomRef> = vec![
            AtomRef::Tag("C1".into()),
            AtomRef::Tag("C2".into()),
            AtomRef::Tag("C3".into()),
            AtomRef::Tag("C4".into()),
            AtomRef::Tag("C5".into()),
            AtomRef::Tag("C6".into()),
        ];
        assert_eq!(ast.aromatic_systems[0].atoms, refs);
    }

    #[rstest]
    #[case::empty(MoleculeAst::default())]
    #[case::indexed_atoms(MoleculeAst {
        atoms: Atoms::indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![LocalizedBond { id: None, a: AtomRef::Index(0), b: AtomRef::Index(1), bond: BondAst::from_order(3) }],
        ..Default::default()
    })]
    fn test_molecule_ast_json_roundtrip(#[case] ast: MoleculeAst) {
        let json = serde_json::to_string(&ast).expect("serialize to JSON");
        let back: MoleculeAst =
            serde_json::from_str(&json).expect(&format!("deserialize from JSON: {json}"));
        assert_eq!(ast, back);
    }

    #[rstest]
    #[case::non_map("3")]
    #[case::missing_atoms(r#"{:bonds []}"#)]
    #[case::missing_bonds(r#"{:atoms ["C"]}"#)]
    #[case::unknown_endpoint(r#"{:atoms [[:C "C"]] :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#)]
    #[case::bad_atom_string(r##"{:atoms ["#h3"] :bonds []}"##)]
    #[case::trailing_content(r#"{:atoms ["C"] :bonds []} :extra :junk"#)]
    #[case::duplicate_tag(r#"{:atoms [[:C "C"] [:C "O"]] :bonds []}"#)]
    #[case::duplicate_alias(r#"{:aliases [:ch "C #h1" :ch "C #h2"] :atoms [] :bonds []}"#)]
    fn test_parse_molecule_dsl_invalid(#[case] input: &str) {
        assert!(
            parse_molecule_dsl(input).is_err(),
            "{input:?} should fail but succeeded"
        );
    }

    #[rstest]
    #[case::plain_vector(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::tagged_vector(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
    #[case::alias_vector(r#"{:atoms [:ch :ch "O"] :bonds [[0 2 :single]] :aliases [:ch "C #h1"]}"#)]
    #[case::mixed(
        r#"{:atoms [:ch [:ring-O :ch] "O"] :bonds [[0 :ring-O :single]] :aliases [:ch "C #h1"]}"#
    )]
    fn test_edn_roundtrip(#[case] input: &str) {
        let ast1 = parse_molecule_dsl(input).unwrap();
        let edn = ast1.to_string();
        let ast2 = parse_molecule_dsl(&edn).unwrap();
        assert_eq!(ast1.atoms.entries, ast2.atoms.entries);
        assert_eq!(ast1.bonds, ast2.bonds);
    }

    #[test]
    fn test_builder_metadata_roundtrip() {
        let ast =
            parse_molecule_dsl(r#"{:atoms [[:C "C #h3"] [:O "O #h1"]] :bonds [[:C :O :single]]}"#)
                .unwrap();

        let cfg = MoleculeDslConfig::zeroed();
        let (builder, meta) = MoleculeBuilder::from_ast_with_metadata(ast.clone(), &cfg).unwrap();
        let ast2 = builder.to_ast_with_metadata(&meta, &cfg);

        assert_eq!(ast.atoms.entries.len(), ast2.atoms.entries.len());
        assert!(ast2.atoms.tags.contains_key("C"));
        assert!(ast2.atoms.tags.contains_key("O"));
        assert_eq!(ast2.bonds[0].a, AtomRef::Tag("C".into()));
        assert_eq!(ast2.bonds[0].b, AtomRef::Tag("O".into()));
    }

    #[test]
    fn test_builder_no_metadata_produces_indexed() {
        let ast =
            parse_molecule_dsl(r#"{:atoms [[:C "C #h3"] [:O "O #h1"]] :bonds [[:C :O :single]]}"#)
                .unwrap();

        let cfg = MoleculeDslConfig::zeroed();
        let builder = MoleculeBuilder::from_ast(ast, &cfg).unwrap();
        let ast2 = builder.to_ast(&cfg);

        assert_eq!(ast2.bonds[0].a, AtomRef::Index(0));
        assert_eq!(ast2.bonds[0].b, AtomRef::Index(1));
        assert!(ast2.atoms.tags.is_empty());
    }

    #[test]
    fn test_alias_tag_disjointness_error() {
        let result = parse_molecule_dsl(r#"{:atoms [[:ch "C"]] :bonds [] :aliases [:ch "C #h1"]}"#);
        assert!(result.is_err(), "alias and tag with same name should fail");
    }
}
