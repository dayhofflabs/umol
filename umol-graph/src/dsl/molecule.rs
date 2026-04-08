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
use umol_edn::{Edn, EdnError, EdnKeyword, EdnMapHelper, EdnStreamDeserializer, FromEdn, ToEdn};

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
    pub id: Option<EdnKeyword>,
    pub a: AtomRef,
    pub b: AtomRef,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, FromEdn, ToEdn)]
pub struct DativeBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<EdnKeyword>,
    pub donor: AtomRef,
    pub acceptor: AtomRef,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, FromEdn, ToEdn)]
pub struct AromaticSystem {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<EdnKeyword>,
    pub atoms: Vec<AtomRef>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, FromEdn, ToEdn)]
pub struct MulticenterBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<EdnKeyword>,
    pub atoms: Vec<AtomRef>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, FromEdn, ToEdn)]
pub struct NoncovalentBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<EdnKeyword>,
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
                "id" => id = Some(map.next_value::<EdnKeyword>()?),
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

/// Parse a molecule AST from an EDN string.
///
/// Uses the native single-pass parser-fusion path (no `Edn` tree, no
/// serde data model). See `discussion/77-umol-edn-usability-review-2026-04-07.md`
/// "Phase 2.5 measurement" for the architectural decision.
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeAst, ParseError> {
    let mut de = EdnStreamDeserializer::new(input);
    let mol_input =
        read_molecule_input(&mut de).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    de.expect_eof().map_err(|e| ParseError::EdnParse(e.to_string()))?;
    mol_input.into_ast()
}

impl FromStr for MoleculeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_molecule_dsl(s)
    }
}

// ---------------------------------------------------------------------------
// Alternative parsing paths (kept for benchmarking and regression detection).
// ---------------------------------------------------------------------------

/// Parse via serde streaming. Retained for benchmark comparison only.
pub fn parse_molecule_dsl_serde(input: &str) -> Result<MoleculeAst, ParseError> {
    let mol_input: MoleculeInput =
        umol_edn::from_str(input).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    mol_input.into_ast()
}

/// Parse via native `FromEdn` walking an intermediate `Edn` tree.
/// Retained for benchmark comparison only.
pub fn parse_molecule_dsl_tree(input: &str) -> Result<MoleculeAst, ParseError> {
    let tree =
        umol_edn::read_string(input).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    let mol_input = MoleculeInput::from_edn(&tree)
        .map_err(|e| ParseError::EdnParse(e.to_string()))?;
    mol_input.into_ast()
}

fn read_molecule_input(de: &mut EdnStreamDeserializer<'_>) -> Result<MoleculeInput, EdnError> {
    de.consume_byte(b'{')?;

    let mut atoms: Option<Vec<AtomEntryInput>> = None;
    let mut bonds: Option<Vec<LocalizedBond>> = None;
    let mut dative_bonds: Vec<DativeBond> = Vec::new();
    let mut aromatic_systems: Vec<AromaticSystem> = Vec::new();
    let mut multicenter_bonds: Vec<MulticenterBond> = Vec::new();
    let mut noncovalent_bonds: Vec<NoncovalentBond> = Vec::new();
    let mut charge: Option<i64> = None;
    let mut spin: Option<SpinState> = None;
    let mut atom_aliases: Vec<String> = Vec::new();

    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "atoms" => atoms = Some(read_seq(de, read_atom_entry)?),
            "bonds" => bonds = Some(read_seq(de, read_localized_bond)?),
            "dative" => dative_bonds = read_seq(de, read_dative_bond)?,
            "aromatic" => aromatic_systems = read_seq(de, read_aromatic_system)?,
            "multicenter" => multicenter_bonds = read_seq(de, read_multicenter_bond)?,
            "noncovalent" => noncovalent_bonds = read_seq(de, read_noncovalent_bond)?,
            "charge" => charge = Some(de.read_i64()?),
            "spin" => {
                let s = de.read_string_or_keyword()?;
                spin = Some(
                    SpinState::from_str(s.as_ref())
                        .map_err(|e| EdnError::Custom(e.to_string()))?,
                );
            }
            "atom-aliases" | "aliases" => {
                atom_aliases = read_seq(de, |d| {
                    d.read_string_or_keyword().map(|c| c.into_owned())
                })?;
            }
            _ => de.read_skip_value()?,
        }
    }

    let atoms = atoms.ok_or_else(|| EdnError::Custom("missing field: atoms".to_string()))?;
    let bonds = bonds.ok_or_else(|| EdnError::Custom("missing field: bonds".to_string()))?;

    Ok(MoleculeInput {
        atoms,
        bonds,
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        atom_aliases,
        charge,
        spin,
    })
}

fn read_seq<T, F>(
    de: &mut EdnStreamDeserializer<'_>,
    mut element: F,
) -> Result<Vec<T>, EdnError>
where
    F: FnMut(&mut EdnStreamDeserializer<'_>) -> Result<T, EdnError>,
{
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            return Ok(out);
        }
        out.push(element(de)?);
    }
}

fn read_atom_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomEntryInput, EdnError> {
    match de.peek_byte()? {
        Some(b'"') | Some(b':') => {
            let s = de.read_string_or_keyword()?;
            Ok(AtomEntryInput::Str(s.into_owned()))
        }
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let tag = de.read_string_or_keyword()?.into_owned();
            let inner = read_atom_entry(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryInput::Tagged(tag, Box::new(inner)))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_atom_ref(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomRef, EdnError> {
    match de.peek_byte()? {
        Some(b':') => Ok(AtomRef::Tag(de.read_keyword_name()?.into_owned())),
        Some(b'"') => Ok(AtomRef::Tag(de.read_string()?.into_owned())),
        Some(b) if b.is_ascii_digit() || b == b'-' || b == b'+' => {
            let n = de.read_i64()?;
            let idx = usize::try_from(n).map_err(|_| EdnError::OutOfRange {
                value: n.to_string(),
                target: "AtomRef::Index",
                path: Vec::new(),
            })?;
            Ok(AtomRef::Index(idx))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_bond_spec(de: &mut EdnStreamDeserializer<'_>) -> Result<BondAst, EdnError> {
    let s = de.read_string_or_keyword()?;
    let aliases = super::bond::builtin_bond_aliases();
    if let Some(ast) = aliases.get_by_left(s.as_ref()) {
        return Ok(ast.clone());
    }
    BondAst::from_str(s.as_ref()).map_err(|e| EdnError::Custom(e.to_string()))
}

fn read_localized_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<LocalizedBond, EdnError> {
    match de.peek_byte()? {
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_spec(de)?;
            de.consume_byte(b']')?;
            Ok(LocalizedBond { id: None, a, b, bond })
        }
        Some(b'{') => {
            let (id, a, b, bond, _, _) =
                read_endpoint_bond_map(de, EndpointBondKind::Localized)?;
            Ok(LocalizedBond { id, a, b, bond })
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_dative_bond(de: &mut EdnStreamDeserializer<'_>) -> Result<DativeBond, EdnError> {
    let (id, donor, acceptor, bond, _, _) =
        read_endpoint_bond_map(de, EndpointBondKind::Dative)?;
    Ok(DativeBond {
        id,
        donor,
        acceptor,
        bond,
    })
}

fn read_noncovalent_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBond, EdnError> {
    let (id, a, b, bond, _, _) =
        read_endpoint_bond_map(de, EndpointBondKind::Noncovalent)?;
    Ok(NoncovalentBond { id, a, b, bond })
}

fn read_aromatic_system(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystem, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(AromaticSystem { id, atoms })
}

fn read_multicenter_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBond, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(MulticenterBond { id, atoms })
}

#[derive(Copy, Clone)]
enum EndpointBondKind {
    Localized,
    Dative,
    Noncovalent,
}

impl EndpointBondKind {
    fn first_key(self) -> &'static str {
        match self {
            EndpointBondKind::Dative => "donor",
            _ => "a",
        }
    }
    fn second_key(self) -> &'static str {
        match self {
            EndpointBondKind::Dative => "acceptor",
            _ => "b",
        }
    }
}

/// Read a map with `:id? :a :b :bond` (or donor/acceptor for dative).
fn read_endpoint_bond_map(
    de: &mut EdnStreamDeserializer<'_>,
    kind: EndpointBondKind,
) -> Result<(Option<EdnKeyword>, AtomRef, AtomRef, BondAst, (), ()), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<EdnKeyword> = None;
    let mut a: Option<AtomRef> = None;
    let mut b: Option<AtomRef> = None;
    let mut bond: Option<BondAst> = None;
    let first_key = kind.first_key();
    let second_key = kind.second_key();

    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        let key_ref = key.as_ref();
        if key_ref == "id" {
            id = Some(EdnKeyword::new(de.read_keyword_name()?.into_owned()));
        } else if key_ref == first_key {
            a = Some(read_atom_ref(de)?);
        } else if key_ref == second_key {
            b = Some(read_atom_ref(de)?);
        } else if key_ref == "bond" {
            bond = Some(read_bond_spec(de)?);
        } else {
            de.read_skip_value()?;
        }
    }

    let a = a.ok_or_else(|| EdnError::Custom(format!("missing field: {}", first_key)))?;
    let b = b.ok_or_else(|| EdnError::Custom(format!("missing field: {}", second_key)))?;
    let bond = bond.ok_or_else(|| EdnError::Custom("missing field: bond".to_string()))?;
    Ok((id, a, b, bond, (), ()))
}

/// Read a map with `:id? :atoms` (aromatic / multicenter).
fn read_atoms_bond_map(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<(Option<EdnKeyword>, Vec<AtomRef>), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<EdnKeyword> = None;
    let mut atoms: Option<Vec<AtomRef>> = None;
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "id" => id = Some(EdnKeyword::new(de.read_keyword_name()?.into_owned())),
            "atoms" => atoms = Some(read_seq(de, read_atom_ref)?),
            _ => de.read_skip_value()?,
        }
    }
    let atoms = atoms.ok_or_else(|| EdnError::Custom("missing field: atoms".to_string()))?;
    Ok((id, atoms))
}

fn unexpected(offset: usize, b: Option<u8>) -> EdnError {
    match b {
        Some(b) => EdnError::UnexpectedToken {
            offset,
            found: b as char,
        },
        None => EdnError::UnexpectedEof { offset },
    }
}

impl<'de> FromEdn<'de> for AtomRef {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Int(n) => {
                let idx = usize::try_from(*n).map_err(|_| EdnError::OutOfRange {
                    value: n.to_string(),
                    target: "AtomRef::Index",
                    path: Vec::new(),
                })?;
                Ok(AtomRef::Index(idx))
            }
            Edn::Keyword(k) => Ok(AtomRef::Tag(k.as_str().to_string())),
            Edn::Str(s) => Ok(AtomRef::Tag(s.to_string())),
            other => Err(EdnError::TypeMismatch {
                expected: "int or keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AtomRef {
    fn to_edn(&self) -> Edn<'_> {
        match self {
            AtomRef::Index(n) => Edn::Int(*n as i64),
            AtomRef::Tag(s) => Edn::keyword(s),
        }
    }
}

impl<'de> FromEdn<'de> for AtomEntryInput {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Str(s) => Ok(AtomEntryInput::Str(s.to_string())),
            Edn::Keyword(k) => Ok(AtomEntryInput::Str(k.as_str().to_string())),
            Edn::Vector(v) | Edn::List(v) => {
                if v.len() != 2 {
                    return Err(EdnError::Custom(format!(
                        "atom entry vector must have length 2, got {}",
                        v.len()
                    )));
                }
                let tag = match &v[0] {
                    Edn::Keyword(k) => k.as_str().to_string(),
                    Edn::Str(s) => s.to_string(),
                    other => {
                        return Err(EdnError::TypeMismatch {
                            expected: "keyword (atom tag)",
                            got: other.kind(),
                            path: Vec::new(),
                        });
                    }
                };
                let inner = AtomEntryInput::from_edn(&v[1])?;
                Ok(AtomEntryInput::Tagged(tag, Box::new(inner)))
            }
            other => Err(EdnError::TypeMismatch {
                expected: "string, keyword, or [tag def] vector",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'de> FromEdn<'de> for LocalizedBond {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Vector(v) | Edn::List(v) => {
                if v.len() != 3 {
                    return Err(EdnError::Custom(format!(
                        "bond shorthand vector must have length 3, got {}",
                        v.len()
                    )));
                }
                Ok(LocalizedBond {
                    id: None,
                    a: AtomRef::from_edn(&v[0])?,
                    b: AtomRef::from_edn(&v[1])?,
                    bond: BondAst::from_edn(&v[2])?,
                })
            }
            Edn::Map(m) => {
                let id = read_optional_id(m)?;
                let mut h = EdnMapHelper::new(m);
                let a: AtomRef = h.required("a")?;
                let b: AtomRef = h.required("b")?;
                let bond: BondAst = h.required("bond")?;
                Ok(LocalizedBond { id, a, b, bond })
            }
            other => Err(EdnError::TypeMismatch {
                expected: "bond vector or map",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

/// Read a flat alias vector entry as `String`. Accepts both keyword and
/// string forms — keywords for alias names, strings for atom-spec defs.
fn alias_entry_to_string(edn: &Edn<'_>) -> Result<String, EdnError> {
    match edn {
        Edn::Keyword(k) => Ok(k.as_str().to_string()),
        Edn::Str(s) => Ok(s.to_string()),
        other => Err(EdnError::TypeMismatch {
            expected: "keyword or string in alias vector",
            got: other.kind(),
            path: Vec::new(),
        }),
    }
}

/// Read an `id` field, which must be a keyword.
fn read_optional_id(map: &umol_edn::EdnMap<'_>) -> Result<Option<EdnKeyword>, EdnError> {
    match map.get_ref(umol_edn::EdnKeyRef::keyword("id")) {
        Some(edn) => Ok(Some(EdnKeyword::from_edn(edn)?)),
        None => Ok(None),
    }
}

impl<'de> FromEdn<'de> for MoleculeInput {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let m = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "molecule map",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        let mut h = EdnMapHelper::new(m);
        let atoms: Vec<AtomEntryInput> = h.required("atoms")?;
        let bonds: Vec<LocalizedBond> = h.required("bonds")?;
        let dative_bonds: Vec<DativeBond> = h.optional("dative")?.unwrap_or_default();
        let aromatic_systems: Vec<AromaticSystem> =
            h.optional("aromatic")?.unwrap_or_default();
        let multicenter_bonds: Vec<MulticenterBond> =
            h.optional("multicenter")?.unwrap_or_default();
        let noncovalent_bonds: Vec<NoncovalentBond> =
            h.optional("noncovalent")?.unwrap_or_default();
        let charge: Option<i64> = h.optional("charge")?;
        let spin_str: Option<String> = h.optional("spin")?;
        let spin = match spin_str {
            Some(s) => Some(SpinState::from_str(&s).map_err(|e| EdnError::Custom(e.to_string()))?),
            None => None,
        };

        // atom-aliases: a flat vector alternating keyword/string entries.
        // Read as Vec<Edn> and project each element to a String via the
        // helper. Falls back to the legacy `aliases` key for compatibility
        // with the existing serde path's tolerance.
        let alias_edn: Option<Vec<Edn<'_>>> = h
            .optional("atom-aliases")?
            .or(h.optional("aliases")?);
        let atom_aliases = match alias_edn {
            Some(v) => v
                .iter()
                .map(alias_entry_to_string)
                .collect::<Result<Vec<_>, _>>()?,
            None => Vec::new(),
        };

        Ok(MoleculeInput {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            atom_aliases,
            charge,
            spin,
        })
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
        bonds: vec![LocalizedBond { id: Some(EdnKeyword::new("b1")), a: AtomRef::Index(0), b: AtomRef::Index(1), bond: BondAst::from_order(1) }], ..Default::default() })]
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
        assert_eq!(ast.aromatic_systems[0].id, Some(EdnKeyword::new("ar1")));
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

    /// All three parsing paths (canonical fused, serde streaming, native
    /// tree) must produce identical ASTs on every accepted input.
    #[rstest]
    #[case::empty(r#"{:atoms [] :bonds []}"#)]
    #[case::plain_atoms(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::tagged_atoms(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
    #[case::bond_id(r#"{:atoms ["H" "F"] :bonds [{:id :b1 :a 0 :b 1 :bond :single}]}"#)]
    #[case::charge(r#"{:atoms [[:F "F#c-"]] :bonds [] :charge -1}"#)]
    #[case::aliases(r#"{:atoms [:ch :ch] :bonds [[0 1 :single]] :aliases [:ch "C #h1"]}"#)]
    #[case::dative(
        r#"{:atoms [[:B "B #h3"] [:N "N #h3"]] :bonds [] :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#
    )]
    fn test_parse_molecule_dsl_path_parity(#[case] input: &str) {
        let canonical = parse_molecule_dsl(input).unwrap();
        let serde_ast = parse_molecule_dsl_serde(input).unwrap();
        let tree_ast = parse_molecule_dsl_tree(input).unwrap();
        assert_eq!(canonical, serde_ast);
        assert_eq!(canonical, tree_ast);
    }

    /// `#[derive(ToEdn)]` output must round-trip through `FromEdn` for the
    /// derived bond types — proves the derived encoder and decoder agree.
    #[rstest]
    #[case::dative_with_id(DativeBond {
        id: Some(EdnKeyword::new("d1")),
        donor: AtomRef::Tag("N".into()),
        acceptor: AtomRef::Tag("B".into()),
        bond: BondAst::from_order(1),
    })]
    #[case::dative_no_id(DativeBond {
        id: None,
        donor: AtomRef::Index(0),
        acceptor: AtomRef::Index(1),
        bond: BondAst::from_order(1),
    })]
    fn test_dative_bond_to_edn_roundtrip(#[case] bond: DativeBond) {
        let edn = umol_edn::ToEdn::to_edn(&bond);
        let back = DativeBond::from_edn(&edn).unwrap();
        assert_eq!(bond, back);
    }

    #[rstest]
    #[case::aromatic(AromaticSystem {
        id: Some(EdnKeyword::new("ar1")),
        atoms: vec![AtomRef::Tag("C1".into()), AtomRef::Tag("C2".into())],
    })]
    #[case::aromatic_empty(AromaticSystem { id: None, atoms: vec![] })]
    fn test_aromatic_system_to_edn_roundtrip(#[case] sys: AromaticSystem) {
        let edn = umol_edn::ToEdn::to_edn(&sys);
        let back = AromaticSystem::from_edn(&edn).unwrap();
        assert_eq!(sys, back);
    }

    #[rstest]
    #[case::noncov(NoncovalentBond {
        id: None,
        a: AtomRef::Index(0),
        b: AtomRef::Index(3),
        bond: BondAst::from_order(1),
    })]
    fn test_noncovalent_bond_to_edn_roundtrip(#[case] bond: NoncovalentBond) {
        let edn = umol_edn::ToEdn::to_edn(&bond);
        let back = NoncovalentBond::from_edn(&edn).unwrap();
        assert_eq!(bond, back);
    }
}
