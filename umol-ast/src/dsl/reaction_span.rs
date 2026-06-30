//! Reaction span DSL: the surface form of `ReactionSpanAst`, where each entity carries its
//! complete before/after value (`EntitySpan`) rather than a delta. Entity ids, bond endpoints,
//! and constraint topology refs are resolved in `into_ast`.

use std::fmt::{self, Display};
use std::str::FromStr;

use indexmap::IndexMap;
use umol_edn::{read_string, DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_graph_core::{EdgeId, Graph};

use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::config::MoleculeDefaults;
use super::constraint::{ConstraintDsl, EntityCounts};
use super::edn_utils::{
    optional_id, pair, parse_vec, read_map, read_vec, required_key, single_key_map, two_atom_refs,
};
use super::error::ParseError;
use super::molecule::{
    parse_atom_aliases, parse_atom_entry, parse_bond_entry, read_atom_aliases, render_atom_value,
    render_bond_entry, resolve_atom_spec, AtomSpecInput, MoleculeMetadata,
};
use super::refs::AtomRef;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::id::{AtomId, BondId};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::{Constraint, ConstraintSpan, EntitySpan, ReactionSpanAst};

/// Surface DSL for a reaction span. Pairs `ReactionSpanAst` with the `MoleculeMetadata` recording
/// its span-frame id↔name bindings; fields private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanDsl {
    ast: ReactionSpanAst,
    metadata: MoleculeMetadata,
}

impl ReactionSpanDsl {
    pub fn from_parts(ast: ReactionSpanAst, metadata: MoleculeMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &ReactionSpanAst {
        &self.ast
    }

    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionSpanAst, MoleculeMetadata) {
        (self.ast, self.metadata)
    }
}

/// Apply `f` to every populated side of an `EntitySpan`.
fn map_span<T, U>(span: &EntitySpan<T>, f: impl Fn(&T) -> U) -> EntitySpan<U> {
    match span {
        EntitySpan::Unchanged(value) => EntitySpan::Unchanged(f(value)),
        EntitySpan::Modified { left, right } => EntitySpan::Modified {
            left: f(left),
            right: f(right),
        },
        EntitySpan::Added(value) => EntitySpan::Added(f(value)),
        EntitySpan::Removed(value) => EntitySpan::Removed(f(value)),
    }
}

impl FromAst<ReactionSpanAst> for ReactionSpanDsl {
    type Ctx = MoleculeDefaults;

    fn from_ast(ast: &ReactionSpanAst, cfg: &Self::Ctx) -> Self {
        let lowered = ReactionSpanAst::from_parts(
            ast.graph().clone(),
            ast.atoms()
                .iter()
                .map(|span| map_span(span, |atom| AtomDsl::from_ast(atom, &cfg.atom).0))
                .collect(),
            ast.bonds()
                .iter()
                .map(|span| map_span(span, |bond| BondDsl::from_ast(bond, &cfg.bond).0))
                .collect(),
            // TODO(I4): overlay span columns through the overlay DSL. Bridged empty so the lib
            // links; the `ReactionSpanDsl` roundtrip drops overlay spans until I4.
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            ast.constraints().to_vec(),
        );
        ReactionSpanDsl::from_parts(lowered, MoleculeMetadata::default())
    }
}

impl IntoAst<ReactionSpanAst> for ReactionSpanDsl {
    type Ctx = MoleculeDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> ReactionSpanAst {
        let ast = self.ast;
        ReactionSpanAst::from_parts(
            ast.graph().clone(),
            ast.atoms()
                .iter()
                .map(|span| map_span(span, |atom| AtomDsl(atom.clone()).into_ast(&cfg.atom)))
                .collect(),
            ast.bonds()
                .iter()
                .map(|span| map_span(span, |bond| BondDsl(bond.clone()).into_ast(&cfg.bond)))
                .collect(),
            // TODO(I4): overlay span columns through the overlay DSL (bridged empty).
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            ast.constraints().to_vec(),
        )
    }
}

/// One molecule-level constraint's span, with its refs still unresolved.
#[derive(Debug, PartialEq)]
pub(crate) enum ConstraintSpanInput {
    Unchanged(ConstraintDsl),
    Added(ConstraintDsl),
    Removed(ConstraintDsl),
}

/// A parsed reaction span before ref resolution: each atom carries an unresolved `AtomSpecInput`
/// (`Bare | Alias`) per `EntitySpan` side, each bond its endpoints + value, plus the `:atom-aliases`
/// table. Aliases, bond endpoints, and constraint refs are resolved in `into_ast`.
#[derive(Debug, PartialEq)]
pub(crate) struct SpanInput {
    atoms: Vec<(Option<String>, EntitySpan<AtomSpecInput>)>,
    bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>)>,
    constraints: Vec<ConstraintSpanInput>,
    atom_aliases: Vec<(String, Box<AtomDsl>)>,
}

const SPAN_VERBS: [&str; 3] = ["add", "modify", "remove"];

/// Split the optional outer `[<id> <body>]` wrapper off a span entry. The id is borrowed.
fn split_span_entry<'a, 'de>(edn: &'a Edn<'de>) -> (Option<&'a str>, &'a Edn<'de>) {
    if let Edn::Vector(v) = edn {
        if v.len() == 2 {
            if let Edn::Keyword(id) = &v[0] {
                return (Some(id.name()), &v[1]);
            }
        }
    }
    (None, edn)
}

/// The verb + payload of a `{:add|:modify|:remove <p>}` wrapper, or `None` for a bare entry. A
/// single-key map whose key is not a verb (e.g. a bare `{:connected …}` constraint) is `None`.
fn verb_wrapper<'a, 'de>(edn: &'a Edn<'de>) -> Option<(&'a str, &'a Edn<'de>)> {
    let Edn::Map(m) = edn else { return None };
    if m.len() != 1 {
        return None;
    }
    let (key, payload) = m.iter().next()?;
    let Edn::Keyword(verb) = key else { return None };
    SPAN_VERBS
        .contains(&verb.name())
        .then_some((verb.name(), payload))
}

fn parse_atom_span_entry(
    edn: &Edn<'_>,
) -> Result<(Option<String>, EntitySpan<AtomSpecInput>), DeError> {
    let (id, body) = split_span_entry(edn);
    let span = match verb_wrapper(body) {
        None => EntitySpan::Unchanged(parse_atom_entry(body)?.spec),
        Some(("add", p)) => EntitySpan::Added(parse_atom_entry(p)?.spec),
        Some(("remove", p)) => EntitySpan::Removed(parse_atom_entry(p)?.spec),
        Some(("modify", p)) => {
            let (left, right) = pair(p, "atom span :modify")?;
            EntitySpan::Modified {
                left: parse_atom_entry(left)?.spec,
                right: parse_atom_entry(right)?.spec,
            }
        }
        Some((verb, _)) => return Err(DeError::Custom(format!("atom span: unexpected verb :{verb}"))),
    };
    Ok((id.map(String::from), span))
}

/// Parse a complete bond-entry payload (`[a b bond]` or the `{:id :atoms :type}` map) and wrap its
/// `BondAst` into the given span side.
#[allow(clippy::type_complexity)]
fn bond_entry_span(
    payload: &Edn<'_>,
    wrap: impl Fn(BondAst) -> EntitySpan<BondAst>,
) -> Result<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>), DeError> {
    let entry = parse_bond_entry(payload)?;
    Ok((entry.id, [entry.first, entry.second], wrap(entry.bond.0)))
}

/// Split a bond `:modify` payload — `[a b X]` or `{:id :atoms [a b] :type X}` — into its id,
/// endpoints, and the raw value `X` (a `[left right]` vector).
fn split_bond_frame<'e>(
    payload: &'e Edn<'e>,
) -> Result<(Option<String>, [AtomRef; 2], &'e Edn<'e>), DeError> {
    match payload {
        Edn::Vector(v) if v.len() == 3 => Ok((
            None,
            [AtomRef::from_edn(&v[0])?, AtomRef::from_edn(&v[1])?],
            &v[2],
        )),
        Edn::Map(m) => Ok((
            optional_id(m)?,
            two_atom_refs(
                parse_vec(required_key(m, "atoms", "bond span")?, ":atoms", |e| {
                    AtomRef::from_edn(e)
                })?,
                "bond span",
            )?,
            required_key(m, "type", "bond span")?,
        )),
        other => Err(DeError::TypeMismatch {
            expected: "bond :modify [a b [left right]] or map",
            got: other.kind(),
            path: vec!["bond span".into()],
        }),
    }
}

#[allow(clippy::type_complexity)]
fn parse_bond_span_entry(
    edn: &Edn<'_>,
) -> Result<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>), DeError> {
    match verb_wrapper(edn) {
        None => bond_entry_span(edn, EntitySpan::Unchanged),
        Some(("add", p)) => bond_entry_span(p, EntitySpan::Added),
        Some(("remove", p)) => bond_entry_span(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let (id, endpoints, value) = split_bond_frame(p)?;
            let (left, right) = pair(value, "bond span :modify")?;
            Ok((
                id,
                endpoints,
                EntitySpan::Modified {
                    left: BondDsl::from_edn(left)?.0,
                    right: BondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!("bond span: unexpected verb :{verb}"))),
    }
}

fn parse_constraint_span_entry(edn: &Edn<'_>) -> Result<ConstraintSpanInput, DeError> {
    match verb_wrapper(edn) {
        None => Ok(ConstraintSpanInput::Unchanged(ConstraintDsl::from_edn(edn)?)),
        Some(("add", p)) => Ok(ConstraintSpanInput::Added(ConstraintDsl::from_edn(p)?)),
        Some(("remove", p)) => Ok(ConstraintSpanInput::Removed(ConstraintDsl::from_edn(p)?)),
        Some((verb, _)) => Err(DeError::Custom(format!(
            "constraint span: unexpected verb :{verb}"
        ))),
    }
}

/// Tree parse of a span map: `:atoms` / `:bonds` / `:constraints` via the entry parsers. A plain
/// molecule map (all entries bare) parses as an all-`Unchanged` span.
fn parse_span_input(edn: &Edn<'_>) -> Result<SpanInput, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "span map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut constraints = Vec::new();
    let mut atom_aliases = Vec::new();
    for (key, value) in m.iter() {
        let Edn::Keyword(key) = key else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: key.kind(),
                path: vec!["span".into()],
            });
        };
        match key.name() {
            "atoms" => atoms = parse_vec(value, ":atoms", parse_atom_span_entry)?,
            "bonds" => bonds = parse_vec(value, ":bonds", parse_bond_span_entry)?,
            "constraints" => {
                constraints = parse_vec(value, ":constraints", parse_constraint_span_entry)?
            }
            "atom-aliases" => atom_aliases = parse_atom_aliases(value)?,
            other => return Err(DeError::Custom(format!("unknown span key :{other}"))),
        }
    }
    Ok(SpanInput {
        atoms,
        bonds,
        constraints,
        atom_aliases,
    })
}

/// Streaming parse of a span map. The span entry grammar is tree-only (it reuses the molecule entry
/// parsers), so each section element is buffered to an `Edn` and dispatched to the tree entry parser.
fn read_span_input(de: &mut EdnStreamDeserializer<'_>) -> Result<SpanInput, EdnError> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut constraints = Vec::new();
    let mut atom_aliases = Vec::new();
    read_map(de, |de, key| {
        match key {
            "atoms" => {
                atoms = read_vec(de, |de| {
                    Ok(parse_atom_span_entry(&read_string(de.read_value_slice()?)?)?)
                })?
            }
            "bonds" => {
                bonds = read_vec(de, |de| {
                    Ok(parse_bond_span_entry(&read_string(de.read_value_slice()?)?)?)
                })?
            }
            "constraints" => {
                constraints = read_vec(de, |de| {
                    Ok(parse_constraint_span_entry(&read_string(de.read_value_slice()?)?)?)
                })?
            }
            "atom-aliases" => atom_aliases = read_atom_aliases(de)?,
            other => return Err(DeError::Custom(format!("unknown span key :{other}")).into()),
        }
        Ok(())
    })?;
    Ok(SpanInput {
        atoms,
        bonds,
        constraints,
        atom_aliases,
    })
}

fn resolve_atom_span(
    span: EntitySpan<AtomSpecInput>,
    aliases: &IndexMap<String, Box<AtomDsl>>,
) -> Result<EntitySpan<AtomAst>, ParseError> {
    Ok(match span {
        EntitySpan::Unchanged(s) => EntitySpan::Unchanged(resolve_atom_spec(s, aliases)?),
        EntitySpan::Added(s) => EntitySpan::Added(resolve_atom_spec(s, aliases)?),
        EntitySpan::Removed(s) => EntitySpan::Removed(resolve_atom_spec(s, aliases)?),
        EntitySpan::Modified { left, right } => EntitySpan::Modified {
            left: resolve_atom_spec(left, aliases)?,
            right: resolve_atom_spec(right, aliases)?,
        },
    })
}

fn resolve_constraint_span(
    input: ConstraintSpanInput,
    counts: &EntityCounts,
    meta: &MoleculeMetadata,
) -> Result<ConstraintSpan, ParseError> {
    Ok(match input {
        ConstraintSpanInput::Unchanged(dsl) => ConstraintSpan::Unchanged(dsl.into_ast(counts, meta)?),
        ConstraintSpanInput::Added(dsl) => ConstraintSpan::Added(dsl.into_ast(counts, meta)?),
        ConstraintSpanInput::Removed(dsl) => ConstraintSpan::Removed(dsl.into_ast(counts, meta)?),
    })
}

impl SpanInput {
    /// Resolve the union-frame span: positions are the union ids (no fresh-id allocation), inline
    /// `:id`s and `:atom-aliases` populate the namespace, atom `AtomSpecInput` sides resolve to
    /// `AtomAst`, bond endpoints and constraint refs resolve against the namespace, and each
    /// projected side must be internally ref-consistent.
    pub(crate) fn into_ast(self) -> Result<(ReactionSpanAst, MoleculeMetadata), ParseError> {
        let atom_count = self.atoms.len();
        let bond_count = self.bonds.len();

        // Namespace: inline ids onto the union positions, then the bijective alias table.
        let mut metadata = MoleculeMetadata::default();
        for (index, (id, _)) in self.atoms.iter().enumerate() {
            if let Some(name) = id {
                if metadata.contains_id(name) {
                    return Err(ParseError::DuplicateId(name.clone()));
                }
                metadata.set_atom_id(AtomId(index as u32), name.clone());
            }
        }
        for (index, (id, _, _)) in self.bonds.iter().enumerate() {
            if let Some(name) = id {
                if metadata.contains_id(name) {
                    return Err(ParseError::DuplicateId(name.clone()));
                }
                metadata.set_bond_id(BondId(index as u32), name.clone());
            }
        }
        let mut aliases: IndexMap<String, Box<AtomDsl>> = IndexMap::new();
        for (name, dsl) in self.atom_aliases {
            if aliases.contains_key(&name) || metadata.contains_id(&name) {
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

        // Resolve atoms (alias → AtomAst), bonds (endpoints + value), constraints.
        let mut atoms: Vec<EntitySpan<AtomAst>> = Vec::with_capacity(atom_count);
        for (_, span) in self.atoms {
            atoms.push(resolve_atom_span(span, &aliases)?);
        }
        let mut bonds: Vec<EntitySpan<BondAst>> = Vec::with_capacity(bond_count);
        let mut endpoints: Vec<[AtomId; 2]> = Vec::with_capacity(bond_count);
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count);
        for (_, [ref_a, ref_b], span) in self.bonds {
            let a = ref_a.into_ast(atom_count, &metadata)?;
            let b = ref_b.into_ast(atom_count, &metadata)?;
            edges.push([a.index() as u32, b.index() as u32]);
            endpoints.push([a, b]);
            bonds.push(span);
        }

        // Per-side ref consistency: a bond present on a side needs both endpoints present there.
        for (span, [a, b]) in bonds.iter().zip(&endpoints) {
            if span.left().is_some()
                && (atoms[a.index()].left().is_none() || atoms[b.index()].left().is_none())
            {
                return Err(ParseError::InvalidValue(
                    "bond present on the left references an atom absent on the left".into(),
                ));
            }
            if span.right().is_some()
                && (atoms[a.index()].right().is_none() || atoms[b.index()].right().is_none())
            {
                return Err(ParseError::InvalidValue(
                    "bond present on the right references an atom absent on the right".into(),
                ));
            }
        }

        let counts = EntityCounts {
            atom_count,
            bond_count,
            dative_bond_count: 0,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
            stereo_atom_count: 0,
            stereo_bond_count: 0,
        };
        let mut constraints: Vec<ConstraintSpan> = Vec::with_capacity(self.constraints.len());
        for input in self.constraints {
            constraints.push(resolve_constraint_span(input, &counts, &metadata)?);
        }

        let graph = Graph::new(atom_count, &edges);
        Ok((
            // TODO(I4): overlay span columns from the parsed DSL (bridged empty).
            ReactionSpanAst::from_parts(
                graph,
                atoms,
                bonds,
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                constraints,
            ),
            metadata,
        ))
    }
}

impl<'de> FromEdn<'de> for ReactionSpanDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let (ast, metadata) = parse_span_input(edn)?
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionSpanDsl::from_parts(ast, metadata))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let span_input = read_span_input(&mut de)?;
        de.expect_eof()?;
        let (ast, metadata) = span_input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionSpanDsl::from_parts(ast, metadata))
    }
}

impl FromStr for ReactionSpanDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ReactionSpanDsl::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

fn render_atom_span_entry(
    id: AtomId,
    span: &EntitySpan<AtomAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let body = match span {
        EntitySpan::Unchanged(a) => render_atom_value(a, meta),
        EntitySpan::Added(a) => single_key_map("add", render_atom_value(a, meta)),
        EntitySpan::Removed(a) => single_key_map("remove", render_atom_value(a, meta)),
        EntitySpan::Modified { left, right } => single_key_map(
            "modify",
            Edn::Vector(vec![render_atom_value(left, meta), render_atom_value(right, meta)].into()),
        ),
    };
    match meta.atom_id(id) {
        Some(name) => {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(name.to_string())), body].into())
        }
        None => body,
    }
}

fn render_bond_span_entry(
    id: BondId,
    endpoints: [AtomId; 2],
    span: &EntitySpan<BondAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |bond: &BondAst| BondDsl::from_ref(bond).to_edn();
    match span {
        EntitySpan::Unchanged(b) => render_bond_entry(id, endpoints, value(b), meta),
        EntitySpan::Added(b) => single_key_map("add", render_bond_entry(id, endpoints, value(b), meta)),
        EntitySpan::Removed(b) => {
            single_key_map("remove", render_bond_entry(id, endpoints, value(b), meta))
        }
        EntitySpan::Modified { left, right } => single_key_map(
            "modify",
            render_bond_entry(
                id,
                endpoints,
                Edn::Vector(vec![value(left), value(right)].into()),
                meta,
            ),
        ),
    }
}

fn render_constraint_span_entry(span: &ConstraintSpan, meta: &MoleculeMetadata) -> Edn<'static> {
    let render = |c: &Constraint| {
        ConstraintDsl::from_ast(c, meta)
            .expect("ConstraintDsl::from_ast is infallible for a well-formed AST")
            .to_edn()
    };
    match span {
        ConstraintSpan::Unchanged(c) => render(c),
        ConstraintSpan::Added(c) => single_key_map("add", render(c)),
        ConstraintSpan::Removed(c) => single_key_map("remove", render(c)),
    }
}

/// Render a `ReactionSpanAst` (+ its `MoleculeMetadata`) to the span EDN map.
fn render_span_edn(ast: &ReactionSpanAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(4);
    let atoms: Vec<Edn<'static>> = ast
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, span)| render_atom_span_entry(AtomId(i as u32), span, meta))
        .collect();
    map.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));

    let bonds: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .enumerate()
        .map(|(j, span)| {
            let [a, b] = ast.graph().edge_endpoints(EdgeId(j as u32));
            render_bond_span_entry(BondId(j as u32), [AtomId::from(a), AtomId::from(b)], span, meta)
        })
        .collect();
    if !bonds.is_empty() {
        map.insert(Edn::keyword("bonds"), Edn::Vector(bonds.into()));
    }

    let constraints: Vec<Edn<'static>> = ast
        .constraints()
        .iter()
        .map(|span| render_constraint_span_entry(span, meta))
        .collect();
    if !constraints.is_empty() {
        map.insert(Edn::keyword("constraints"), Edn::Vector(constraints.into()));
    }

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

impl ToEdn for ReactionSpanDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_span_edn(&self.ast, &self.metadata)
    }
}

impl Display for ReactionSpanDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for ReactionSpanAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        ReactionSpanDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        ReactionSpanDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for ReactionSpanAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for ReactionSpanAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `ReactionSpanAst`: positional form (no `:id` keywords, no aliases) since
/// the AST carries no metadata. For id/alias-bearing output, wrap in [`ReactionSpanDsl`].
impl ToEdn for ReactionSpanAst {
    fn to_edn(&self) -> Edn<'static> {
        render_span_edn(self, &MoleculeMetadata::default())
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use umol_edn::read_string;

    use super::*;
    use crate::ast::boolean::BooleanAst;
    use crate::ast::constraint::{BondConstraint, Constraint, Constraints, MoleculeConstraint};
    use crate::ast::delta::{AtomDelta, BondDelta, ConstraintDelta, Delta, Deltas};
    use crate::ast::edit::BondFieldChange;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::reaction::ReactionAst;
    use crate::ast::value::ValueAst;

    // Modified bond + Unchanged atoms + Unchanged molecule-constraint.
    #[rstest]
    #[case::modify(ReactionAst::new(
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![], vec![],
            Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
        ),
        Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
        })]),
    ))]
    // Unchanged / Removed / Added atoms and bonds + an Added constraint.
    #[case::add_remove(ReactionAst::new(
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove { id: AtomId(1), ast: AtomAst::from_element(Element::O) }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::from_order(1),
            }),
            Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::N) }),
            Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(2)],
                ast: BondAst::from_order(1),
            }),
            Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
            )),
        ]),
    ))]
    fn test_reaction_span_dsl_from_ast(#[case] reaction: ReactionAst) {
        let span = reaction.to_reaction_span().unwrap();
        let cfg = MoleculeDefaults::default();
        assert_eq!(ReactionSpanDsl::from_ast(&span, &cfg).into_ast(&cfg), span);
    }

    #[rstest]
    #[case::unchanged(r#""C""#, (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))))))]
    #[case::add(r#"{:add "O"}"#, (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))))))]
    #[case::remove(r#"{:remove "O"}"#, (None, EntitySpan::Removed(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))))))]
    #[case::modify(r#"{:modify ["C" "N"]}"#, (None, EntitySpan::Modified {
        left: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
        right: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::N)))),
    }))]
    #[case::with_id(r#"[:c "C"]"#, (Some("c".to_string()), EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))))))]
    #[case::alias(r#":nu"#, (None, EntitySpan::Unchanged(AtomSpecInput::Alias("nu".to_string()))))]
    #[case::add_alias(r#"{:add :nu}"#, (None, EntitySpan::Added(AtomSpecInput::Alias("nu".to_string()))))]
    fn test_parse_atom_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, EntitySpan<AtomSpecInput>),
    ) {
        assert_eq!(
            parse_atom_span_entry(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unchanged("[0 1 :single]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondAst::from_order(1))))]
    #[case::add("{:add [0 2 :single]}", (None, [AtomRef::Index(0), AtomRef::Index(2)], EntitySpan::Added(BondAst::from_order(1))))]
    #[case::remove("{:remove [0 1 :single]}", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Removed(BondAst::from_order(1))))]
    #[case::modify("{:modify [0 2 [:single :double]]}", (None, [AtomRef::Index(0), AtomRef::Index(2)], EntitySpan::Modified {
        left: BondAst::from_order(1),
        right: BondAst::from_order(2),
    }))]
    #[case::unchanged_map_id("{:id :b1 :atoms [0 1] :type :single}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondAst::from_order(1))))]
    #[case::modify_map_id("{:modify {:id :b1 :atoms [0 1] :type [:single :double]}}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Modified {
        left: BondAst::from_order(1),
        right: BondAst::from_order(2),
    }))]
    #[case::triple("[0 1 :triple]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondAst::from_order(3))))]
    #[case::aromatic("[0 1 :aromatic]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(
        BondAst::from_order(1).with_constraint(BondConstraint::Aromatic(BooleanAst::Lit(true))),
    )))]
    fn test_parse_bond_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, [AtomRef; 2], EntitySpan<BondAst>),
    ) {
        assert_eq!(
            parse_bond_span_entry(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:connected {}}"#, ConstraintSpanInput::Unchanged(
        ConstraintDsl::from_edn(&read_string(r#"{:connected {}}"#).unwrap()).unwrap(),
    ))]
    #[case::add(r#"{:add {:connected {}}}"#, ConstraintSpanInput::Added(
        ConstraintDsl::from_edn(&read_string(r#"{:connected {}}"#).unwrap()).unwrap(),
    ))]
    #[case::remove(r#"{:remove {:connected {}}}"#, ConstraintSpanInput::Removed(
        ConstraintDsl::from_edn(&read_string(r#"{:connected {}}"#).unwrap()).unwrap(),
    ))]
    fn test_parse_constraint_span_entry(#[case] input: &str, #[case] expected: ConstraintSpanInput) {
        assert_eq!(
            parse_constraint_span_entry(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::full(
        r#"{:atoms ["C" {:add "O"}] :bonds [[0 1 :single]] :constraints [{:connected {}}]}"#,
        SpanInput {
            atoms: vec![
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))))),
                (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))))),
            ],
            bonds: vec![(
                None,
                [AtomRef::Index(0), AtomRef::Index(1)],
                EntitySpan::Unchanged(BondAst::from_order(1)),
            )],
            constraints: vec![ConstraintSpanInput::Unchanged(
                ConstraintDsl::from_edn(&read_string(r#"{:connected {}}"#).unwrap()).unwrap(),
            )],
            atom_aliases: vec![],
        },
    )]
    #[case::plain_molecule(
        r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#,
        SpanInput {
            atoms: vec![
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))))),
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))))),
            ],
            bonds: vec![(
                None,
                [AtomRef::Index(0), AtomRef::Index(1)],
                EntitySpan::Unchanged(BondAst::from_order(1)),
            )],
            constraints: vec![],
            atom_aliases: vec![],
        },
    )]
    #[case::with_aliases(
        r#"{:atoms [:nu {:add "O"}] :atom-aliases [:nu "C#h3"]}"#,
        SpanInput {
            atoms: vec![
                (None, EntitySpan::Unchanged(AtomSpecInput::Alias("nu".to_string()))),
                (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))))),
            ],
            bonds: vec![],
            constraints: vec![],
            atom_aliases: vec![(
                "nu".to_string(),
                Box::new(AtomDsl::from_edn(&read_string(r#""C#h3""#).unwrap()).unwrap()),
            )],
        },
    )]
    fn test_parse_span_input(#[case] input: &str, #[case] expected: SpanInput) {
        assert_eq!(parse_span_input(&read_string(input).unwrap()).unwrap(), expected);
        let mut de = EdnStreamDeserializer::new(input);
        assert_eq!(read_span_input(&mut de).unwrap(), expected);
    }

    #[rstest]
    #[case::alias_and_add(
        r#"{:atoms [:nu {:add "O"}] :bonds [{:add [0 1 :single]}] :atom-aliases [:nu "C"]}"#,
        ReactionSpanAst::from_parts(
            Graph::new(2, &[[0, 1]]),
            vec![
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Added(AtomAst::from_element(Element::O)),
            ],
            vec![EntitySpan::Added(BondAst::from_order(1))],
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            vec![],
        ),
        MoleculeMetadata::new().with_atom_alias("nu", AtomDsl(AtomAst::from_element(Element::C))),
    )]
    fn test_span_input_into_ast(
        #[case] input: &str,
        #[case] expected_ast: ReactionSpanAst,
        #[case] expected_metadata: MoleculeMetadata,
    ) {
        assert_eq!(
            parse_span_input(&read_string(input).unwrap())
                .unwrap()
                .into_ast()
                .unwrap(),
            (expected_ast, expected_metadata),
        );
    }

    #[rstest]
    #[case::unknown_alias(
        r#"{:atoms [:nu]}"#,
        ParseError::InvalidValue("unknown atom alias :nu".to_string()),
    )]
    #[case::unknown_ref(
        r#"{:atoms ["C"] :bonds [[0 5 :single]]}"#,
        ParseError::InvalidRef { kind: "atom", value: "5".to_string() },
    )]
    #[case::left_inconsistent(
        r#"{:atoms ["C" {:add "O"}] :bonds [[0 1 :single]]}"#,
        ParseError::InvalidValue(
            "bond present on the left references an atom absent on the left".to_string(),
        ),
    )]
    fn test_span_input_into_ast_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(
            parse_span_input(&read_string(input).unwrap())
                .unwrap()
                .into_ast()
                .unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::span(
        r#"{:atoms ["C" {:add "O"}] :bonds [{:add [0 1 :single]}]}"#,
        ReactionSpanDsl::from_parts(
            ReactionSpanAst::from_parts(
                Graph::new(2, &[[0, 1]]),
                vec![
                    EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                    EntitySpan::Added(AtomAst::from_element(Element::O)),
                ],
                vec![EntitySpan::Added(BondAst::from_order(1))],
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                vec![],
            ),
            MoleculeMetadata::new(),
        ),
    )]
    #[case::plain_molecule(
        r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#,
        ReactionSpanDsl::from_parts(
            ReactionSpanAst::from_parts(
                Graph::new(2, &[[0, 1]]),
                vec![
                    EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                    EntitySpan::Unchanged(AtomAst::from_element(Element::O)),
                ],
                vec![EntitySpan::Unchanged(BondAst::from_order(1))],
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                Default::default(),
                vec![],
            ),
            MoleculeMetadata::new(),
        ),
    )]
    fn test_reaction_span_dsl_from_edn(#[case] input: &str, #[case] expected: ReactionSpanDsl) {
        assert_eq!(
            ReactionSpanDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
        assert_eq!(ReactionSpanDsl::from_edn_str(input).unwrap(), expected);
        assert_eq!(ReactionSpanDsl::from_str(input).unwrap(), expected);
    }

    #[rstest]
    #[case::unchanged(AtomId(0), EntitySpan::Unchanged(AtomAst::from_element(Element::C)), MoleculeMetadata::new(), r#""C""#)]
    #[case::add(AtomId(0), EntitySpan::Added(AtomAst::from_element(Element::O)), MoleculeMetadata::new(), r#"{:add "O"}"#)]
    #[case::remove(AtomId(0), EntitySpan::Removed(AtomAst::from_element(Element::O)), MoleculeMetadata::new(), r#"{:remove "O"}"#)]
    #[case::modify(AtomId(0), EntitySpan::Modified { left: AtomAst::from_element(Element::C), right: AtomAst::from_element(Element::N) }, MoleculeMetadata::new(), r#"{:modify ["C" "N"]}"#)]
    #[case::with_id(AtomId(0), EntitySpan::Unchanged(AtomAst::from_element(Element::C)), MoleculeMetadata::new().with_atom_id(AtomId(0), "c"), r#"[:c "C"]"#)]
    #[case::alias(AtomId(0), EntitySpan::Unchanged(AtomAst::from_element(Element::C)), MoleculeMetadata::new().with_atom_alias("nu", AtomDsl(AtomAst::from_element(Element::C))), r#":nu"#)]
    fn test_render_atom_span_entry(
        #[case] id: AtomId,
        #[case] span: EntitySpan<AtomAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_atom_span_entry(id, &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(BondAst::from_order(1)), MoleculeMetadata::new(), "[0 1 :single]")]
    #[case::add(EntitySpan::Added(BondAst::from_order(2)), MoleculeMetadata::new(), "{:add [0 1 :double]}")]
    #[case::remove(EntitySpan::Removed(BondAst::from_order(1)), MoleculeMetadata::new(), "{:remove [0 1 :single]}")]
    #[case::modify(EntitySpan::Modified { left: BondAst::from_order(1), right: BondAst::from_order(2) }, MoleculeMetadata::new(), "{:modify [0 1 [:single :double]]}")]
    #[case::with_id(EntitySpan::Unchanged(BondAst::from_order(1)), MoleculeMetadata::new().with_bond_id(BondId(0), "b1"), "{:id :b1 :atoms [0 1] :type :single}")]
    #[case::aromatic(EntitySpan::Unchanged(BondAst::from_order(1).with_constraint(BondConstraint::Aromatic(BooleanAst::Lit(true)))), MoleculeMetadata::new(), "[0 1 :aromatic]")]
    fn test_render_bond_span_entry(
        #[case] span: EntitySpan<BondAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_bond_span_entry(BondId(0), [AtomId(0), AtomId(1)], &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })), "{:connected {}}")]
    #[case::add(ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })), "{:add {:connected {}}}")]
    #[case::remove(ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })), "{:remove {:connected {}}}")]
    fn test_render_constraint_span_entry(#[case] span: ConstraintSpan, #[case] expected: &str) {
        assert_eq!(
            render_constraint_span_entry(&span, &MoleculeMetadata::new()),
            read_string(expected).unwrap(),
        );
    }

    // `ReactionSpanDsl` (ast + metadata) round-trips through the EDN surface (render → reparse).
    #[rstest]
    #[case::plain_molecule(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::modify(r#"{:atoms ["C" "C"] :bonds [{:modify [0 1 [:single :double]]}]}"#)]
    #[case::add_remove(r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :single]}]}"#)]
    #[case::constraint(r#"{:atoms ["C"] :constraints [{:connected {}} {:add {:connected {}}}]}"#)]
    #[case::aliases(r#"{:atoms [:nu {:add "O"}] :atom-aliases [:nu "C"]}"#)]
    fn test_reaction_span_dsl_to_edn(#[case] input: &str) {
        let dsl = ReactionSpanDsl::from_str(input).unwrap();
        assert_eq!(ReactionSpanDsl::from_edn(&dsl.to_edn()).unwrap(), dsl);
    }

    // `ReactionSpanAst` renders positionally and round-trips (metadata discarded on both sides).
    #[rstest]
    #[case::plain_molecule(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::modify(r#"{:atoms ["C" "C"] :bonds [{:modify [0 1 [:single :double]]}]}"#)]
    #[case::add_remove(r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :single]}]}"#)]
    fn test_reaction_span_ast_to_edn(#[case] input: &str) {
        let ast = ReactionSpanAst::from_str(input).unwrap();
        assert_eq!(ReactionSpanAst::from_edn(&ast.to_edn()).unwrap(), ast);
    }

    // An `:add` constraint lands on the right projection only; `:remove` on the left only.
    #[rstest]
    #[case::add(
        r#"{:atoms ["C"] :constraints [{:add {:connected {}}}]}"#,
        vec![],
        vec![Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })],
    )]
    #[case::remove(
        r#"{:atoms ["C"] :constraints [{:remove {:connected {}}}]}"#,
        vec![Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })],
        vec![],
    )]
    fn test_reaction_span_ast_constraint_projection(
        #[case] input: &str,
        #[case] left: Vec<Constraint>,
        #[case] right: Vec<Constraint>,
    ) {
        let ast = ReactionSpanAst::from_str(input).unwrap();
        assert_eq!(ast.left().constraints().as_slice(), left.as_slice());
        assert_eq!(ast.right().constraints().as_slice(), right.as_slice());
    }
}
