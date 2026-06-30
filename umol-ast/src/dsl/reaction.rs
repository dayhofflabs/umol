//! Reaction DSL.
//!
//! `ReactionDsl` wraps a `ReactionAst` together with the `ReactionMetadata` that records
//! the surface-form bindings (the lhs molecule metadata plus created-entity id ↔ name and
//! atom-alias bindings). The EDN form is a map keyed by `:lhs` (a molecule map) and
//! `:deltas` (a vector of `:add` / `:remove` / `:modify` / `:constraint` operations). Each
//! entity delegates to its own entity DSL.

use std::fmt::{self, Display};
use std::str::FromStr;

use bimap::BiBTreeMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::atom::{lower_atom, raise_atom, AtomDsl, PartialAtomDsl};
use super::bond::{lower_bond, raise_bond, BondDsl, PartialBondDsl};
use super::config::{DeltaDefaults, ReactionDefaults};
use super::constraint::{read_constraint_dsl, ConstraintDsl, EntityCounts};
use super::edn_utils::{
    consume_single_key_map_close, missing, parse_single_key_map, parse_vec,
    read_single_key_map_header, read_vec, single_key_map,
};
use super::error::ParseError;
use super::molecule::{
    parse_atom_aliases, parse_atom_entry, parse_bond_entry, parse_molecule_input,
    read_atom_aliases, read_atom_entry, read_bond_entry, read_molecule_input, render_molecule_edn,
    resolve_atom_spec, AtomEntryInput, BondEntryInput, MoleculeDsl, MoleculeInput, MoleculeMetadata,
};
use super::refs::{read_atom_ref, read_bond_ref, AtomRef, BondRef};
use crate::ast::atom::{AtomAst, ElementAst};
use crate::ast::bond::BondAst;
use crate::ast::delta::{AtomDelta, BondDelta, ConstraintDelta, Delta, Deltas};
use crate::ast::edit::{AtomFieldChange, BondFieldChange};
use crate::ast::id::{AtomId, BondId};
use crate::ast::reaction::ReactionAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::EntityPatch;

/// Surface DSL for a reaction. Pairs `ReactionAst` with `ReactionMetadata`; fields are
/// private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionDsl {
    ast: ReactionAst,
    metadata: ReactionMetadata,
}

impl ReactionDsl {
    pub fn from_parts(ast: ReactionAst, metadata: ReactionMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &ReactionAst {
        &self.ast
    }

    pub fn metadata(&self) -> &ReactionMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionAst, ReactionMetadata) {
        (self.ast, self.metadata)
    }
}

impl FromAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn from_ast(ast: &ReactionAst, cfg: &Self::Ctx) -> Self {
        let lhs = MoleculeDsl::from_ast(&ast.lhs, &cfg.molecule_defaults())
            .into_parts()
            .0;
        let delta_cfg = cfg.delta_defaults();
        let mut deltas = ast.deltas.clone();
        for delta in deltas.iter_mut() {
            lower_delta(delta, &delta_cfg);
        }
        ReactionDsl {
            ast: ReactionAst { lhs, deltas },
            metadata: ReactionMetadata::default(),
        }
    }
}

impl IntoAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> ReactionAst {
        let ReactionAst { lhs, mut deltas } = self.ast;
        let lhs = MoleculeDsl::from_parts(lhs, MoleculeMetadata::default())
            .into_ast(&cfg.molecule_defaults());
        let delta_cfg = cfg.delta_defaults();
        for delta in deltas.iter_mut() {
            raise_delta(delta, &delta_cfg);
        }
        ReactionAst { lhs, deltas }
    }
}

/// Lower a delta's embedded entity AST to DSL-display form (AST → DSL).
fn lower_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            lower_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            lower_bond(ast, &cfg.bond)
        }
        _ => {}
    }
}

/// Raise a delta's embedded entity AST from DSL-display form (DSL → AST).
fn raise_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            raise_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            raise_bond(ast, &cfg.bond)
        }
        _ => {}
    }
}

/// Surface-form metadata paired with a `ReactionAst`: the lhs molecule metadata plus the
/// created-entity id bindings and atom aliases introduced by the deltas. Mirrors
/// `MoleculeMetadata` for the atom/bond entities (the reaction admits the `[:C "C#h3"]`
/// alias notation for added atoms).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    atom_ids: IndexMap<AtomId, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondId, String>,
}

impl ReactionMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn lhs(&self) -> &MoleculeMetadata {
        &self.lhs
    }

    pub fn combined_metadata(&self) -> MoleculeMetadata {
        let mut combined = self.lhs.clone();
        for (&id, name) in &self.atom_ids {
            combined.set_atom_id(id, name);
        }
        for (&id, name) in &self.bond_ids {
            combined.set_bond_id(id, name);
        }
        combined
    }

    pub fn atom_id(&self, id: AtomId) -> Option<&str> {
        self.atom_ids.get(&id).map(String::as_str)
    }

    pub fn bond_id(&self, id: BondId) -> Option<&str> {
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

    pub fn set_atom_id(&mut self, id: AtomId, name: impl Into<String>) {
        self.atom_ids.insert(id, name.into());
    }

    pub fn set_bond_id(&mut self, id: BondId, name: impl Into<String>) {
        self.bond_ids.insert(id, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases.insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_id(mut self, id: AtomId, name: impl Into<String>) -> Self {
        self.set_atom_id(id, name);
        self
    }

    pub fn with_bond_id(mut self, id: BondId, name: impl Into<String>) -> Self {
        self.set_bond_id(id, name);
        self
    }

    pub fn with_atom_alias(mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

/// One unresolved delta parsed from a `:deltas` entry. Refs stay symbolic
/// (`AtomRef`/`BondRef`) and an `:add` carries the molecule entry verbatim
/// (bare atom / alias / created-id resolved in R7); a `:modify` RHS is a
/// partial entity AST (unspecified fields `Undetermined`).
#[derive(Debug, PartialEq)]
pub(crate) enum DeltaInput {
    AtomAdd(AtomEntryInput),
    AtomRemove(AtomRef),
    AtomModify(AtomRef, AtomAst),
    BondAdd(BondEntryInput),
    BondRemove(BondRef),
    BondModify(BondRef, BondAst),
    ConstraintAdd(ConstraintDsl),
    ConstraintRemove(ConstraintDsl),
}

/// Raw parse target for a reaction: the lhs molecule input plus the unresolved
/// deltas. Resolution (`into_ast`, R7) lifts this to `(ReactionAst, ReactionMetadata)`.
#[derive(Debug, PartialEq)]
pub(crate) struct ReactionInput {
    lhs: MoleculeInput,
    atom_aliases: Vec<(String, Box<AtomDsl>)>,
    deltas: Vec<DeltaInput>,
}

impl ReactionInput {
    pub(crate) fn into_ast(self) -> Result<(ReactionAst, ReactionMetadata), ParseError> {
        let ReactionInput {
            lhs,
            atom_aliases,
            deltas,
        } = self;
        let (lhs, lhs_meta) = lhs.into_ast()?;

        // Alias table for `:add` = lhs aliases ∪ reaction aliases (bijective; collisions error).
        let mut aliases: IndexMap<String, Box<AtomDsl>> = lhs_meta
            .iter_atom_aliases()
            .map(|(name, dsl)| (name.to_string(), Box::new(dsl.clone())))
            .collect();
        let mut metadata = ReactionMetadata {
            lhs: lhs_meta,
            ..Default::default()
        };
        for (name, dsl) in atom_aliases {
            if aliases.contains_key(&name) {
                return Err(ParseError::DuplicateId(name));
            }
            if aliases.values().any(|existing| existing == &dsl) {
                return Err(ParseError::InvalidValue(
                    "atom-aliases must be bijective: two names map to the same atom".into(),
                ));
            }
            metadata.add_atom_alias(name.clone(), (*dsl).clone());
            aliases.insert(name, dsl);
        }

        // Resolution namespace: lhs entity ids (all kinds) and running counts, both seeded from the
        // lhs molecule and grown in delta order as entities are defined. Every ref — entity and
        // constraint — resolves against this pair; `metadata` separately records created-entity ids
        // for roundtrip. No forward refs: only adds processed earlier are visible.
        let mut counts = EntityCounts::from_ast(&lhs);
        let mut namespace = metadata.lhs().clone();
        let lhs_atom_count = counts.atom_count;
        let lhs_bond_count = counts.bond_count;

        let mut resolved = Deltas::new();
        for delta in deltas {
            match delta {
                DeltaInput::AtomAdd(entry) => {
                    let id = counts.allocate_atom();
                    if let Some(name) = entry.id {
                        if namespace.contains_id(&name) || aliases.contains_key(&name) {
                            return Err(ParseError::DuplicateId(name));
                        }
                        namespace.set_atom_id(id, name.clone());
                        metadata.set_atom_id(id, name);
                    }
                    let ast = resolve_atom_spec(entry.spec, &aliases)?;
                    resolved.push(Delta::Atom(AtomDelta::Add { id, ast }));
                }
                DeltaInput::AtomRemove(r) => {
                    let id = r.into_ast(counts.atom_count, &namespace)?;
                    if id.index() >= lhs_atom_count {
                        return Err(ParseError::InvalidValue(format!(
                            "cannot remove atom :{} added in the same reaction",
                            id.index()
                        )));
                    }
                    resolved.push(Delta::Atom(AtomDelta::Remove {
                        id,
                        ast: lhs[id].clone(),
                    }));
                }
                DeltaInput::AtomModify(r, rhs) => {
                    let id = r.into_ast(counts.atom_count, &namespace)?;
                    if id.index() >= lhs_atom_count {
                        return Err(ParseError::InvalidValue(format!(
                            "cannot modify atom :{} added in the same reaction",
                            id.index()
                        )));
                    }
                    let new = lhs[id].update(&rhs);
                    for d in AtomDelta::diff(id, &lhs[id], &new) {
                        resolved.push(Delta::Atom(d));
                    }
                }
                DeltaInput::BondAdd(entry) => {
                    let id = counts.allocate_bond();
                    if let Some(name) = entry.id {
                        if namespace.contains_id(&name) || aliases.contains_key(&name) {
                            return Err(ParseError::DuplicateId(name));
                        }
                        namespace.set_bond_id(id, name.clone());
                        metadata.set_bond_id(id, name);
                    }
                    let a = entry.first.into_ast(counts.atom_count, &namespace)?;
                    let b = entry.second.into_ast(counts.atom_count, &namespace)?;
                    resolved.push(Delta::Bond(BondDelta::Add {
                        id,
                        atoms: [a, b],
                        ast: entry.bond.0,
                    }));
                }
                DeltaInput::BondRemove(r) => {
                    let id = r.into_ast(counts.bond_count, &namespace)?;
                    if id.index() >= lhs_bond_count {
                        return Err(ParseError::InvalidValue(format!(
                            "cannot remove bond :{} added in the same reaction",
                            id.index()
                        )));
                    }
                    resolved.push(Delta::Bond(BondDelta::Remove {
                        id,
                        atoms: lhs.bond(id).atom_ids(),
                        ast: lhs[id].clone(),
                    }));
                }
                DeltaInput::BondModify(r, rhs) => {
                    let id = r.into_ast(counts.bond_count, &namespace)?;
                    if id.index() >= lhs_bond_count {
                        return Err(ParseError::InvalidValue(format!(
                            "cannot modify bond :{} added in the same reaction",
                            id.index()
                        )));
                    }
                    let new = lhs[id].update(&rhs);
                    for d in BondDelta::diff(id, &lhs[id], &new) {
                        resolved.push(Delta::Bond(d));
                    }
                }
                DeltaInput::ConstraintAdd(dsl) => {
                    let c = dsl.into_ast(&counts, &namespace)?;
                    resolved.push(Delta::Constraint(ConstraintDelta::Add(c)));
                }
                DeltaInput::ConstraintRemove(dsl) => {
                    let c = dsl.into_ast(&counts, &namespace)?;
                    resolved.push(Delta::Constraint(ConstraintDelta::Remove(c)));
                }
            }
        }
        Ok((
            ReactionAst {
                lhs,
                deltas: resolved,
            },
            metadata,
        ))
    }
}

impl<'de> FromEdn<'de> for ReactionDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let input = parse_reaction_input(edn)?;
        let (ast, metadata) = input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionDsl::from_parts(ast, metadata))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let ri = read_reaction_input(&mut de)?;
        de.expect_eof()?;
        let (ast, metadata) = ri.into_ast().map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionDsl::from_parts(ast, metadata))
    }
}

impl FromStr for ReactionDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ReactionDsl::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl ToEdn for ReactionDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_reaction_edn(&self.ast, &self.metadata)
    }
}

impl Display for ReactionDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for ReactionAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        ReactionDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        ReactionDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for ReactionAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for ReactionAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `ReactionAst`. The lhs and refs emit canonical positional
/// form (no id keywords, no aliases) since the AST carries no metadata. For id/alias-bearing
/// surface output, wrap in [`ReactionDsl`] with the appropriate [`ReactionMetadata`].
impl ToEdn for ReactionAst {
    fn to_edn(&self) -> Edn<'static> {
        render_reaction_edn(self, &ReactionMetadata::default())
    }
}

fn read_reaction_input(de: &mut EdnStreamDeserializer<'_>) -> Result<ReactionInput, EdnError> {
    de.consume_byte(b'{')?;
    let mut lhs = None;
    let mut atom_aliases = Vec::new();
    let mut deltas = Vec::new();
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?.into_owned();
        match key.as_str() {
            "lhs" => lhs = Some(read_molecule_input(de)?),
            "atom-aliases" => atom_aliases = read_atom_aliases(de)?,
            "deltas" => deltas = read_vec(de, read_delta_input)?,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["reaction".into()],
                }
                .into())
            }
        }
    }
    Ok(ReactionInput {
        lhs: lhs.ok_or_else(|| missing("lhs", "reaction"))?,
        atom_aliases,
        deltas,
    })
}

fn parse_reaction_input(edn: &Edn<'_>) -> Result<ReactionInput, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "reaction map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut lhs = None;
    let mut atom_aliases = Vec::new();
    let mut deltas = Vec::new();
    for (k, v) in m.iter() {
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["reaction".into()],
            });
        };
        match key.name() {
            "lhs" => lhs = Some(parse_molecule_input(v)?),
            "atom-aliases" => atom_aliases = parse_atom_aliases(v)?,
            "deltas" => deltas = parse_vec(v, ":deltas", parse_delta_input)?,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["reaction".into()],
                })
            }
        }
    }
    Ok(ReactionInput {
        lhs: lhs.ok_or(DeError::MissingField {
            key: "lhs".to_string(),
            path: vec!["reaction".into()],
        })?,
        atom_aliases,
        deltas,
    })
}

fn read_delta_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let entity = read_single_key_map_header(de)?;
    let input = match entity.as_str() {
        "atom" => read_delta_atom_input(de)?,
        "bond" => read_delta_bond_input(de)?,
        "constraint" => read_delta_constraint_input(de)?,
        e => return Err(DeError::Custom(format!("unknown reaction delta :{e}")).into()),
    };
    consume_single_key_map_close(de, "delta")?;
    Ok(input)
}

fn read_delta_atom_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::AtomAdd(read_atom_entry(de)?),
        "remove" => DeltaInput::AtomRemove(read_atom_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_atom_ref(de)?;
            let s = de.read_string()?;
            let dsl: PartialAtomDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("partial-atom", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("atom :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::AtomModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown atom delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "atom delta")?;
    Ok(input)
}

fn read_delta_bond_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::BondAdd(read_bond_entry(de)?),
        "remove" => DeltaInput::BondRemove(read_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: PartialBondDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("partial-bond", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("bond :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::BondModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "bond delta")?;
    Ok(input)
}

fn read_delta_constraint_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::ConstraintAdd(read_constraint_dsl(de)?),
        "remove" => DeltaInput::ConstraintRemove(read_constraint_dsl(de)?),
        o => return Err(DeError::Custom(format!("unknown constraint delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "constraint delta")?;
    Ok(input)
}

fn parse_delta_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (entity, body) = parse_single_key_map(edn, "delta")?;
    match entity {
        "atom" => parse_delta_atom_input(body),
        "bond" => parse_delta_bond_input(body),
        "constraint" => parse_delta_constraint_input(body),
        e => Err(DeError::Custom(format!("unknown reaction delta :{e}"))),
    }
}

fn parse_delta_atom_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "atom delta")?;
    match op {
        "add" => Ok(DeltaInput::AtomAdd(parse_atom_entry(payload)?)),
        "remove" => Ok(DeltaInput::AtomRemove(AtomRef::from_edn(payload)?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "atom :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["atom delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "atom :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::AtomModify(
                AtomRef::from_edn(&v[0])?,
                PartialAtomDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!("unknown atom delta op :{o}"))),
    }
}

fn parse_delta_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "bond delta")?;
    match op {
        "add" => Ok(DeltaInput::BondAdd(parse_bond_entry(payload)?)),
        "remove" => Ok(DeltaInput::BondRemove(BondRef::from_edn(payload)?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::BondModify(
                BondRef::from_edn(&v[0])?,
                PartialBondDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!("unknown bond delta op :{o}"))),
    }
}

fn parse_delta_constraint_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "constraint delta")?;
    match op {
        "add" => Ok(DeltaInput::ConstraintAdd(ConstraintDsl::from_edn(payload)?)),
        "remove" => Ok(DeltaInput::ConstraintRemove(ConstraintDsl::from_edn(
            payload,
        )?)),
        o => Err(DeError::Custom(format!("unknown constraint delta op :{o}"))),
    }
}

/// Render a reaction to its EDN map: `:lhs` (the molecule via the molecule renderer), `:deltas`, and
/// `:atom-aliases` (reaction-level, when present). Aliases render last, as in the molecule surface.
fn render_reaction_edn(ast: &ReactionAst, meta: &ReactionMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(3);
    map.insert(
        Edn::keyword("lhs"),
        render_molecule_edn(&ast.lhs, meta.lhs()),
    );
    map.insert(
        Edn::keyword("deltas"),
        Edn::Vector(render_deltas(&ast.deltas, meta).into()),
    );
    if meta.has_atom_aliases() {
        let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases_len() * 2);
        for (name, dsl) in meta.iter_atom_aliases() {
            pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
            pairs.push(dsl.to_edn());
        }
        map.insert(Edn::keyword("atom-aliases"), Edn::Vector(pairs.into()));
    }
    Edn::Map(map)
}

/// Render `deltas` to the `:deltas` vector entries, resolving refs against the reaction `meta`.
/// Consecutive `ModifyField` / `ModifyConstraint` for one atom coalesce into a single
/// `{:atom {:modify [<ref> <partial>]}}`.
fn render_deltas(deltas: &Deltas, meta: &ReactionMetadata) -> Vec<Edn<'static>> {
    let deltas = deltas.as_slice();
    let mut out = Vec::new();
    // Built once on the first constraint delta (constraint refs resolve against lhs ∪ created).
    let mut combined_metadata: Option<MoleculeMetadata> = None;
    let mut i = 0;
    while i < deltas.len() {
        match &deltas[i] {
            Delta::Atom(AtomDelta::Add { id, ast }) => {
                out.push(single_key_map(
                    "atom",
                    single_key_map("add", render_atom_entry(*id, ast, meta)),
                ));
                i += 1;
            }
            Delta::Atom(AtomDelta::Remove { id, .. }) => {
                out.push(single_key_map(
                    "atom",
                    single_key_map("remove", render_atom_ref(*id, meta)),
                ));
                i += 1;
            }
            Delta::Atom(
                AtomDelta::ModifyField { id, .. } | AtomDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut partial = AtomAst::new(ElementAst::Undetermined);
                while let Some(Delta::Atom(delta)) = deltas.get(i) {
                    match delta {
                        AtomDelta::ModifyField { id: j, change } if *j == id => match change {
                            AtomFieldChange::Element { new, .. } => partial.element = new.clone(),
                            AtomFieldChange::IsotopeMass { new, .. } => {
                                partial.isotope_mass = new.clone()
                            }
                            AtomFieldChange::Charge { new, .. } => partial.charge = new.clone(),
                            AtomFieldChange::ImplicitHydrogens { new, .. } => {
                                partial.implicit_hydrogens = new.clone()
                            }
                            AtomFieldChange::LonePairs { new, .. } => {
                                partial.lone_pairs = new.clone()
                            }
                            AtomFieldChange::Spin { new, .. } => partial.spin = new.clone(),
                        },
                        AtomDelta::ModifyConstraint { id: j, old, new } if *j == id => match new {
                            Some(c) => {
                                partial.constraints.add(c.clone());
                            }
                            None => {
                                if let Some(old) = old {
                                    partial.constraints.add(old.as_undetermined());
                                }
                            }
                        },
                        _ => break,
                    }
                    i += 1;
                }
                let payload = Edn::Vector(
                    vec![render_atom_ref(id, meta), PartialAtomDsl(partial).to_edn()].into(),
                );
                out.push(single_key_map("atom", single_key_map("modify", payload)));
            }
            Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                out.push(single_key_map(
                    "bond",
                    single_key_map("add", render_bond_entry(*id, *atoms, ast, meta)),
                ));
                i += 1;
            }
            Delta::Bond(BondDelta::Remove { id, .. }) => {
                out.push(single_key_map(
                    "bond",
                    single_key_map("remove", render_bond_ref(*id, meta)),
                ));
                i += 1;
            }
            Delta::Bond(
                BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut partial = BondAst::default();
                while let Some(Delta::Bond(delta)) = deltas.get(i) {
                    match delta {
                        BondDelta::ModifyField { id: j, change } if *j == id => match change {
                            BondFieldChange::Order { new, .. } => partial.order = new.clone(),
                            BondFieldChange::Charge { new, .. } => partial.charge = new.clone(),
                            BondFieldChange::Spin { new, .. } => partial.spin = new.clone(),
                        },
                        BondDelta::ModifyConstraint { id: j, old, new } if *j == id => match new {
                            Some(c) => {
                                partial.constraints.add(c.clone());
                            }
                            None => {
                                if let Some(old) = old {
                                    partial.constraints.add(old.as_undetermined());
                                }
                            }
                        },
                        _ => break,
                    }
                    i += 1;
                }
                let payload = Edn::Vector(
                    vec![render_bond_ref(id, meta), PartialBondDsl(partial).to_edn()].into(),
                );
                out.push(single_key_map("bond", single_key_map("modify", payload)));
            }
            Delta::Constraint(delta) => {
                let (op, constraint) = match delta {
                    ConstraintDelta::Add(c) => ("add", c),
                    ConstraintDelta::Remove(c) => ("remove", c),
                };
                let combined = combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let dsl = ConstraintDsl::from_ast(constraint, combined)
                    .expect("ConstraintDsl::from_ast is infallible for a well-formed AST");
                out.push(single_key_map(
                    "constraint",
                    single_key_map(op, dsl.to_edn()),
                ));
                i += 1;
            }
            // TODO(I4): overlay-delta DSL rendering. Bridged so the lib links and the non-DSL
            // tests run; no current test renders an overlay delta.
            Delta::DativeBond(_)
            | Delta::AromaticSystem(_)
            | Delta::MulticenterBond(_)
            | Delta::NoncovalentBond(_) => todo!("I4: overlay delta DSL rendering"),
        }
    }
    out
}

/// A delta ref (`:remove` / `:modify`) names an existing lhs entity — resolved against the lhs frame.
fn render_atom_ref(id: AtomId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta.lhs().atom_id(id) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A created atom (`:add`) — its id and aliases live in the reaction frame (the alias namespace is
/// the lhs ∪ reaction union). Renders `<atom-dsl>` or `[<id> <atom-dsl>]`.
fn render_atom_entry(id: AtomId, atom: &AtomAst, meta: &ReactionMetadata) -> Edn<'static> {
    let dsl = AtomDsl::from_ref(atom);
    let spec = match meta
        .atom_alias_for(dsl)
        .or_else(|| meta.lhs().atom_alias_for(dsl))
    {
        Some(alias) => Edn::Keyword(EdnKeyword::owned(alias.to_string())),
        None => dsl.to_edn(),
    };
    match meta.atom_id(id) {
        Some(name) => {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(name.to_string())), spec].into())
        }
        None => spec,
    }
}

/// An atom named as a bond endpoint — resolved against the union namespace (lhs ∪ created), since a
/// bond may attach to a same-reaction atom. Unlike a delta target ref, which is lhs-only.
fn render_atom_endpoint(id: AtomId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta.atom_id(id).or_else(|| meta.lhs().atom_id(id)) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A bond delta target (`:remove` / `:modify`) names an existing lhs bond — resolved lhs-frame only.
fn render_bond_ref(id: BondId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta.lhs().bond_id(id) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A created bond (`:add`). Renders `[<a> <b> <bond-dsl>]`, or `{:id <id> :atoms [<a> <b>] :type
/// <bond-dsl>}` when the bond carries an id. Endpoints resolve against the union namespace.
fn render_bond_entry(
    id: BondId,
    atoms: [AtomId; 2],
    ast: &BondAst,
    meta: &ReactionMetadata,
) -> Edn<'static> {
    let bond_edn = BondDsl::from_ref(ast).to_edn();
    let first = render_atom_endpoint(atoms[0], meta);
    let second = render_atom_endpoint(atoms[1], meta);
    match meta.bond_id(id) {
        Some(name) => {
            let mut m = EdnMap::with_capacity(3);
            m.insert(
                Edn::keyword("id"),
                Edn::Keyword(EdnKeyword::owned(name.to_string())),
            );
            m.insert(
                Edn::keyword("atoms"),
                Edn::Vector(vec![first, second].into()),
            );
            m.insert(Edn::keyword("type"), bond_edn);
            Edn::Map(m)
        }
        None => Edn::Vector(vec![first, second, bond_edn].into()),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::atom::ElementAst;
    use crate::ast::boolean::BooleanAst;
    use crate::ast::constraint::{AtomConstraint, BondConstraint, Constraint, MoleculeConstraint};
    use crate::ast::delta::{ConstraintDelta, Deltas};
    use crate::ast::edit::{AtomFieldChange, BondFieldChange};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::value::ValueAst;
    use crate::dsl::bond::BondDsl;
    use crate::dsl::constraint::MoleculeConstraintDsl;
    use crate::dsl::molecule::AtomSpecInput;
    use crate::mol;

    #[rstest]
    #[case::sn2(ReactionAst::new(
        mol!(r##"{:atoms ["C" "Br"] :bonds [[0 1 "1"]]}"##),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::O) }),
            Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(2)],
                ast: BondAst::from_order(1),
            }),
            Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(1),
                change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
            }),
            Delta::Constraint(ConstraintDelta::Add(Constraint::And(vec![]))),
        ]),
    ))]
    fn test_reaction_dsl_from_ast_roundtrip(#[case] reaction: ReactionAst) {
        let cfg = ReactionDefaults::ground();
        let dsl = ReactionDsl::from_parts(reaction, ReactionMetadata::default());
        let lowered = ReactionDsl::from_ast(&dsl.clone().into_ast(&cfg), &cfg);
        assert_eq!(lowered.ast(), dsl.ast());
    }

    #[rstest]
    #[case::add_bare(
        r##"{:atom {:add "C#h3"}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_id(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: Some("nu".into()),
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::O));
                a.implicit_hydrogens = ValueAst::Lit(1);
                a
            }))),
        })
    )]
    #[case::add_alias(
        "{:atom {:add :foo}}",
        DeltaInput::AtomAdd(AtomEntryInput {
            id: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_id("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Id("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(AtomRef::Id("br".into()), {
            let mut a = AtomAst::new(ElementAst::Undetermined);
            a.charge = ValueAst::Lit(-1);
            a
        })
    )]
    fn test_parse_delta_input_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_bare(
        r##"{:atom {:add "C#h3"}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_id(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: Some("nu".into()),
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::O));
                a.implicit_hydrogens = ValueAst::Lit(1);
                a
            }))),
        })
    )]
    #[case::add_alias(
        "{:atom {:add :foo}}",
        DeltaInput::AtomAdd(AtomEntryInput {
            id: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_id("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Id("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(AtomRef::Id("br".into()), {
            let mut a = AtomAst::new(ElementAst::Undetermined);
            a.charge = ValueAst::Lit(-1);
            a
        })
    )]
    fn test_read_delta_input_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_vec(
        r##"{:bond {:add [0 1 "1"]}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: None,
            first: AtomRef::Index(0),
            second: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_id(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: Some("b1".into()),
            first: AtomRef::Id("c".into()),
            second: AtomRef::Id("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_id("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Id("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(BondRef::Id("b1".into()), BondAst::from_order(2))
    )]
    fn test_parse_delta_input_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_vec(
        r##"{:bond {:add [0 1 "1"]}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: None,
            first: AtomRef::Index(0),
            second: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_id(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: Some("b1".into()),
            first: AtomRef::Id("c".into()),
            second: AtomRef::Id("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_id("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Id("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(BondRef::Id("b1".into()), BondAst::from_order(2))
    )]
    fn test_read_delta_input_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add(
        "{:constraint {:add {:connected {}}}}",
        DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    #[case::remove(
        "{:constraint {:remove {:connected {}}}}",
        DeltaInput::ConstraintRemove(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    fn test_parse_delta_input_constraint(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add(
        "{:constraint {:add {:connected {}}}}",
        DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    #[case::remove(
        "{:constraint {:remove {:connected {}}}}",
        DeltaInput::ConstraintRemove(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    fn test_read_delta_input_constraint(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_parse_reaction_input() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:remove 0}} {:constraint {:add {:connected {}}}}]}"##;
        let expected = ReactionInput {
            lhs: MoleculeInput {
                atoms: vec![AtomEntryInput {
                    id: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
                }],
                ..Default::default()
            },
            atom_aliases: Vec::new(),
            deltas: vec![
                DeltaInput::AtomAdd(AtomEntryInput {
                    id: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))),
                }),
                DeltaInput::BondRemove(BondRef::Index(0)),
                DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(
                    MoleculeConstraintDsl::Connected { atoms: None },
                )),
            ],
        };
        assert_eq!(
            parse_reaction_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_read_reaction_input() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:remove 0}} {:constraint {:add {:connected {}}}}]}"##;
        let expected = ReactionInput {
            lhs: MoleculeInput {
                atoms: vec![AtomEntryInput {
                    id: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
                }],
                ..Default::default()
            },
            atom_aliases: Vec::new(),
            deltas: vec![
                DeltaInput::AtomAdd(AtomEntryInput {
                    id: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))),
                }),
                DeltaInput::BondRemove(BondRef::Index(0)),
                DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(
                    MoleculeConstraintDsl::Connected { atoms: None },
                )),
            ],
        };
        assert_eq!(
            read_reaction_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast() {
        let input = r##"{:lhs {:atoms ["C"]} :atom-aliases [:me "C#h3"] :deltas [{:atom {:add [:nu :me]}} {:atom {:add "O"}}]}"##;
        let (ast, meta) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: {
                        let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                        a.implicit_hydrogens = ValueAst::Lit(3);
                        a
                    },
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::O),
                }),
            ]),
        );
        assert_eq!(meta.atom_id(AtomId(1)), Some("nu"));
        assert_eq!(meta.atom_id(AtomId(2)), None);
        assert!(meta.has_atom_alias("me"));
    }

    #[rstest]
    fn test_reaction_input_into_ast_alias_union() {
        // `:lo` is an lhs alias, `:hi` a reaction alias; both `:add` resolve (union),
        // but each set stays in its own metadata slot for independent round-trip.
        let input = r##"{:lhs {:atoms ["C"] :atom-aliases [:lo "N"]} :atom-aliases [:hi "C#h3"] :deltas [{:atom {:add :lo}} {:atom {:add :hi}}]}"##;
        let (ast, meta) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: {
                        let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                        a.implicit_hydrogens = ValueAst::Lit(3);
                        a
                    },
                }),
            ]),
        );
        assert_eq!(meta.lhs().atom_aliases_len(), 1);
        assert!(meta.lhs().has_atom_alias("lo"));
        assert_eq!(meta.atom_aliases_len(), 1);
        assert!(meta.has_atom_alias("hi"));
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_remove() {
        let input = r##"{:lhs {:atoms [[:br "Br"] "C"]} :deltas [{:atom {:remove :br}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::Br),
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_remove_error() {
        // Adding then removing the same id is prohibited — recover-from-lhs cannot reach an added atom.
        let input =
            r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:x "O"]}} {:atom {:remove :x}}]}"##;
        let err = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap_err();
        assert!(matches!(err, ParseError::InvalidValue(_)));
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_modify() {
        // lhs :br is Br#c0; modify charge to -1 → one ModifyField (old recovered from lhs).
        let input = r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(0),
                    new: ValueAst::Lit(-1),
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_add() {
        // Bond-add attaches to a same-reaction atom: :o is AtomId(1), bond endpoints [0, :o].
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_remove() {
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:remove :b1}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::from_order(1),
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_remove_error() {
        // Adding then removing the same bond is prohibited — recover-from-lhs cannot reach it.
        let input =
            r##"{:lhs {:atoms ["C" "O"]} :deltas [{:bond {:add [0 1 "1"]}} {:bond {:remove 0}}]}"##;
        let err = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap_err();
        assert!(matches!(err, ParseError::InvalidValue(_)));
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_modify() {
        // lhs :b1 is order 1; modify to order 2 → one ModifyField (old recovered from lhs).
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_add() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None },)
            ))]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_remove() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:remove {:connected {}}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None },)
            ))]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_added_atom_ref() {
        // The constraint ref :o names an atom added in the same reaction (AtomId(1)), resolved
        // against the unified namespace.
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:constraint {:add {:atom [:o {:valence 2}]}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Constraint(ConstraintDelta::Add(Constraint::Atom(
                    AtomId(1),
                    AtomConstraint::Valence(ValueAst::Lit(2)),
                ))),
            ]),
        );
    }

    #[rstest]
    #[case::atom_modify(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##,
        ReactionAst {
            lhs: MoleculeAst::from_edn_str(r##"{:atoms [[:br "Br#c0"]]}"##).unwrap(),
            deltas: Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
            })]),
        }
    )]
    #[case::atom_add_bond_add(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##,
        ReactionAst {
            lhs: MoleculeAst::from_edn_str(r##"{:atoms ["C"]}"##).unwrap(),
            deltas: Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        }
    )]
    fn test_reaction_dsl_from_edn(#[case] input: &str, #[case] expected: ReactionAst) {
        let dsl = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        assert_eq!(dsl.ast(), &expected);
    }

    #[rstest]
    #[case::atom_modify(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##
    )]
    #[case::bond_modify_and_constraint(
        r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}} {:constraint {:add {:connected {}}}}]}"##
    )]
    #[case::atom_add_bond_add(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##
    )]
    fn test_reaction_dsl_from_edn_str_from_edn_parity(#[case] input: &str) {
        let via_tree = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        let via_stream = ReactionDsl::from_edn_str(input).unwrap();
        assert_eq!(via_tree, via_stream);
    }

    #[rstest]
    fn test_reaction_dsl_from_str() {
        let input = r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##;
        let dsl: ReactionDsl = input.parse().unwrap();
        assert_eq!(dsl, ReactionDsl::from_edn_str(input).unwrap());
    }

    /// Shared render metadata: lhs atoms br(0) c(1), bonds b1(0) bx(1); created atom n(2), bond b2(2).
    #[fixture]
    fn meta() -> ReactionMetadata {
        ReactionMetadata {
            lhs: MoleculeMetadata::new()
                .with_atom_id(AtomId(0), "br")
                .with_atom_id(AtomId(1), "c")
                .with_bond_id(BondId(0), "b1")
                .with_bond_id(BondId(1), "bx"),
            ..Default::default()
        }
        .with_atom_id(AtomId(2), "n")
        .with_bond_id(BondId(2), "b2")
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add(vec![Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::N) })], r##"{:atom {:add [:n "N"]}}"##)]
    #[case::remove(vec![Delta::Atom(AtomDelta::Remove { id: AtomId(1), ast: AtomAst::from_element(Element::C) })], "{:atom {:remove :c}}")]
    #[case::modify_field(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) } })], r##"{:atom {:modify [:br "#c-"]}}"##)]
    #[case::modify_set_constraint(vec![Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: None, new: Some(AtomConstraint::valence(4_i64)) })], r##"{:atom {:modify [:br "#v4"]}}"##)]
    #[case::modify_remove_constraint(vec![Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: Some(AtomConstraint::valence(4_i64)), new: None })], r##"{:atom {:modify [:br "#v*"]}}"##)]
    #[case::modify_coalesced(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) } }), Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: Some(AtomConstraint::valence(4_i64)), new: None })], r##"{:atom {:modify [:br "#c-#v*"]}}"##)]
    fn test_render_deltas_atom(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add(vec![Delta::Bond(BondDelta::Add { id: BondId(2), atoms: [AtomId(1), AtomId(2)], ast: BondAst::from_order(1) })], "{:bond {:add {:id :b2 :atoms [:c :n] :type :single}}}")]
    #[case::remove(vec![Delta::Bond(BondDelta::Remove { id: BondId(1), atoms: [AtomId(0), AtomId(1)], ast: BondAst::from_order(1) })], "{:bond {:remove :bx}}")]
    #[case::modify_field(vec![Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) } })], r##"{:bond {:modify [:b1 "2"]}}"##)]
    #[case::modify_constraint(vec![Delta::Bond(BondDelta::ModifyConstraint { id: BondId(0), old: None, new: Some(BondConstraint::Aromatic(BooleanAst::Lit(true))) })], r##"{:bond {:modify [:b1 "#a"]}}"##)]
    fn test_render_deltas_bond(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add_molecule(vec![Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })))], "{:constraint {:add {:connected {}}}}")]
    #[case::add_entity_leaf(vec![Delta::Constraint(ConstraintDelta::Add(Constraint::Atom(AtomId(2), AtomConstraint::Valence(ValueAst::Lit(2)))))], "{:constraint {:add {:atom [:n {:valence 2}]}}}")]
    #[case::remove(vec![Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })))], "{:constraint {:remove {:connected {}}}}")]
    fn test_render_deltas_constraint(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##)]
    #[case::reaction_alias(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add :nu}}] :atom-aliases [:nu "O#h1#c-1"] }"##)]
    #[case::molecule_constraint(r##"{:lhs {:atoms ["C" "N"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}} {:constraint {:add {:connected {}}}}]}"##)]
    #[case::entity_leaf_constraint(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:constraint {:add {:atom [:o {:valence 2}]}}}]}"##)]
    fn test_reaction_dsl_from_edn_to_edn_roundtrip(#[case] input: &str) {
        let dsl = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        let reparsed = ReactionDsl::from_edn(&dsl.to_edn()).unwrap();
        assert_eq!(reparsed, dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##)]
    #[case::add_atom_and_bond(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:add [0 1 "1"]}}]}"##)]
    #[case::molecule_constraint(r##"{:lhs {:atoms ["C" "N"] :bonds [[0 1 "1"]]} :deltas [{:constraint {:add {:connected {}}}}]}"##)]
    #[case::entity_leaf_constraint(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:constraint {:add {:atom [1 {:valence 2}]}}}]}"##)]
    fn test_reaction_ast_to_edn(#[case] input: &str) {
        let ast = ReactionAst::from_edn(&read_string(input).unwrap()).unwrap();
        let reparsed = ReactionAst::from_edn(&ast.to_edn()).unwrap();
        assert_eq!(reparsed, ast);
    }
}
