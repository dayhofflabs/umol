//! Reaction span DSL: the surface form of `ReactionSpanAst`, where each entity carries its
//! complete before/after value (`EntitySpan`) rather than a delta. Entity ids, bond endpoints,
//! and constraint topology refs are resolved in `into_ast`.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{
    read_string, DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn,
};
use umol_graph_core::{
    EdgeId, FactorOrdering, FixedRelationSet, FixedVarBirelationSet, Graph, NodeId,
    RelationParticipant, VarRelationSet,
};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::config::MoleculeDefaults;
use super::constraint::ConstraintDsl;
use super::dative::DativeBondDsl;
use super::edn_utils::{
    optional_id, pair, parse_vec, read_map, read_vec, required_key, single_key_map, two_atom_refs,
};
use super::error::ParseError;
use super::molecule::{
    parse_aromatic_system_entry, parse_atom_aliases, parse_atom_entry, parse_bond_entry,
    parse_dative_bond_entry, parse_multicenter_bond_entry, parse_noncovalent_bond_entry,
    parse_stereo_atom_entry, parse_stereo_bond_entry, read_atom_aliases, render_aromatic_entry,
    render_atom_value, render_bond_entry, render_dative_entry, render_multicenter_entry,
    render_noncovalent_entry, render_stereo_atom_entry, render_stereo_bond_entry,
    render_stereo_ligand, resolve_atom_spec, AtomSpecInput, MoleculeMetadata,
};
use super::multicenter::MulticenterBondDsl;
use super::namespace::MoleculeNamespace;
use super::noncovalent::NoncovalentBondDsl;
use super::refs::{parse_stereo_ligand, AtomRef, BondRef, StereoLigandRef};
use super::stereo::{StereoAtomDsl, StereoBondDsl};
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::{
    AromaticSystemAst, Constraint, ConstraintSpan, DativeBondAst, EntitySpan, MulticenterBondAst,
    NoncovalentBondAst, ReactionSpanAst, StereoAtomAst, StereoBondAst,
};

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
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => EntitySpan::Modified {
            lhs: f(left),
            rhs: f(right),
        },
        EntitySpan::Added(value) => EntitySpan::Added(f(value)),
        EntitySpan::Removed(value) => EntitySpan::Removed(f(value)),
    }
}

/// Rebuild an overlay relation set with its per-side payload mapped through `f`; participants
/// (the overlay topology) carry through unchanged. One per relation-set shape.
fn map_fixed_var_span<L1, O1, const N1: usize, L2, O2, T, U>(
    set: &FixedVarBirelationSet<L1, O1, N1, L2, O2, EntitySpan<T>>,
    f: impl Fn(&T) -> U,
) -> FixedVarBirelationSet<L1, O1, N1, L2, O2, EntitySpan<U>>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    FixedVarBirelationSet::new(
        set.relation_ids()
            .map(|rid| {
                (
                    *set.participants_1(rid),
                    set.participants_2(rid).to_vec(),
                    map_span(set.data(rid), &f),
                )
            })
            .collect(),
    )
}

fn map_var_span<P, O, T, U>(
    set: &VarRelationSet<P, O, EntitySpan<T>>,
    f: impl Fn(&T) -> U,
) -> VarRelationSet<P, O, EntitySpan<U>>
where
    P: RelationParticipant,
    O: FactorOrdering,
{
    VarRelationSet::new(
        set.relation_ids()
            .map(|rid| (set.participants(rid).to_vec(), map_span(set.data(rid), &f)))
            .collect(),
    )
}

fn map_fixed_span<P, O, const N: usize, T, U>(
    set: &FixedRelationSet<P, O, EntitySpan<T>, N>,
    f: impl Fn(&T) -> U,
) -> FixedRelationSet<P, O, EntitySpan<U>, N>
where
    P: RelationParticipant,
    O: FactorOrdering,
{
    FixedRelationSet::new(
        set.relation_ids()
            .map(|rid| (*set.participants(rid), map_span(set.data(rid), &f)))
            .collect(),
    )
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
            map_fixed_var_span(ast.dative_bonds(), |b| {
                DativeBondDsl::from_ast(b, &cfg.dative_bond).0
            }),
            map_var_span(ast.aromatic_systems(), |s| {
                AromaticSystemDsl::from_ast(s, &cfg.aromatic_system).0
            }),
            map_var_span(ast.multicenter_bonds(), |b| {
                MulticenterBondDsl::from_ast(b, &cfg.multicenter_bond).0
            }),
            map_fixed_span(ast.noncovalent_bonds(), |b| {
                NoncovalentBondDsl::from_ast(b, &cfg.noncovalent_bond).0
            }),
            map_fixed_var_span(ast.stereo_atoms(), |s| {
                StereoAtomDsl::from_ast(s, &cfg.stereo_atom).0
            }),
            map_fixed_var_span(ast.stereo_bonds(), |s| {
                StereoBondDsl::from_ast(s, &cfg.stereo_bond).0
            }),
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
            map_fixed_var_span(ast.dative_bonds(), |b| {
                DativeBondDsl(b.clone()).into_ast(&cfg.dative_bond)
            }),
            map_var_span(ast.aromatic_systems(), |s| {
                AromaticSystemDsl(s.clone()).into_ast(&cfg.aromatic_system)
            }),
            map_var_span(ast.multicenter_bonds(), |b| {
                MulticenterBondDsl(b.clone()).into_ast(&cfg.multicenter_bond)
            }),
            map_fixed_span(ast.noncovalent_bonds(), |b| {
                NoncovalentBondDsl(b.clone()).into_ast(&cfg.noncovalent_bond)
            }),
            map_fixed_var_span(ast.stereo_atoms(), |s| {
                StereoAtomDsl(s.clone()).into_ast(&cfg.stereo_atom)
            }),
            map_fixed_var_span(ast.stereo_bonds(), |s| {
                StereoBondDsl(s.clone()).into_ast(&cfg.stereo_bond)
            }),
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
#[allow(clippy::type_complexity)]
pub(crate) struct SpanInput {
    atoms: Vec<(Option<String>, EntitySpan<AtomSpecInput>)>,
    bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>)>,
    dative_bonds: Vec<(
        Option<String>,
        Vec<AtomRef>,
        AtomRef,
        EntitySpan<DativeBondAst>,
    )>,
    aromatic_systems: Vec<(Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemAst>)>,
    multicenter_bonds: Vec<(Option<String>, Vec<AtomRef>, EntitySpan<MulticenterBondAst>)>,
    noncovalent_bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<NoncovalentBondAst>)>,
    stereo_atoms: Vec<(
        Option<String>,
        AtomRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoAtomAst>,
    )>,
    stereo_bonds: Vec<(
        Option<String>,
        BondRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoBondAst>,
    )>,
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
                lhs: parse_atom_entry(left)?.spec,
                rhs: parse_atom_entry(right)?.spec,
            }
        }
        Some((verb, _)) => {
            return Err(DeError::Custom(format!(
                "atom span: unexpected verb :{verb}"
            )))
        }
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
                    lhs: BondDsl::from_edn(left)?.0,
                    rhs: BondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "bond span: unexpected verb :{verb}"
        ))),
    }
}

fn parse_constraint_span_entry(edn: &Edn<'_>) -> Result<ConstraintSpanInput, DeError> {
    match verb_wrapper(edn) {
        None => Ok(ConstraintSpanInput::Unchanged(ConstraintDsl::from_edn(
            edn,
        )?)),
        Some(("add", p)) => Ok(ConstraintSpanInput::Added(ConstraintDsl::from_edn(p)?)),
        Some(("remove", p)) => Ok(ConstraintSpanInput::Removed(ConstraintDsl::from_edn(p)?)),
        Some((verb, _)) => Err(DeError::Custom(format!(
            "constraint span: unexpected verb :{verb}"
        ))),
    }
}

// Overlay span entries: a bare entry (`Unchanged`) or a `{:add|:modify|:remove <entry>}` wrapper,
// over the shared `parse_<entity>_entry`. For `:add`/`:remove` and the bare form the entry's `:type`
// is one dsl; for `:modify` it is a `[left right]` pair (participants are span-invariant).

#[allow(clippy::type_complexity)]
fn parse_dative_span_entry(
    edn: &Edn<'_>,
) -> Result<
    (
        Option<String>,
        Vec<AtomRef>,
        AtomRef,
        EntitySpan<DativeBondAst>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(DativeBondAst) -> EntitySpan<DativeBondAst>| {
        let e = parse_dative_bond_entry(p)?;
        Ok::<_, DeError>((e.id, e.donors, e.acceptor, wrap(e.bond.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "dative span map",
                    got: p.kind(),
                    path: vec!["dative span".into()],
                });
            };
            let donors = parse_vec(required_key(m, "donors", "dative span")?, ":donors", |e| {
                AtomRef::from_edn(e)
            })?;
            let acceptor = AtomRef::from_edn(required_key(m, "acceptor", "dative span")?)?;
            let (left, right) = pair(
                required_key(m, "type", "dative span")?,
                "dative span :modify",
            )?;
            Ok((
                optional_id(m)?,
                donors,
                acceptor,
                EntitySpan::Modified {
                    lhs: DativeBondDsl::from_edn(left)?.0,
                    rhs: DativeBondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "dative span: unexpected verb :{verb}"
        ))),
    }
}

#[allow(clippy::type_complexity)]
fn parse_aromatic_span_entry(
    edn: &Edn<'_>,
) -> Result<(Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemAst>), DeError> {
    let full = |p: &Edn<'_>, wrap: fn(AromaticSystemAst) -> EntitySpan<AromaticSystemAst>| {
        let e = parse_aromatic_system_entry(p)?;
        Ok::<_, DeError>((e.id, e.atoms, wrap(e.system.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "aromatic span map",
                    got: p.kind(),
                    path: vec!["aromatic span".into()],
                });
            };
            let atoms = parse_vec(required_key(m, "atoms", "aromatic span")?, ":atoms", |e| {
                AtomRef::from_edn(e)
            })?;
            let (left, right) = pair(
                required_key(m, "type", "aromatic span")?,
                "aromatic span :modify",
            )?;
            Ok((
                optional_id(m)?,
                atoms,
                EntitySpan::Modified {
                    lhs: AromaticSystemDsl::from_edn(left)?.0,
                    rhs: AromaticSystemDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "aromatic span: unexpected verb :{verb}"
        ))),
    }
}

// `|e| AtomRef::from_edn(e)` reads redundant but is load-bearing: the bare fn item can't satisfy the
// higher-ranked `for<'a> Fn(&Edn<'a>)` bound `parse_vec` requires.
#[allow(clippy::type_complexity, clippy::redundant_closure)]
fn parse_multicenter_span_entry(
    edn: &Edn<'_>,
) -> Result<(Option<String>, Vec<AtomRef>, EntitySpan<MulticenterBondAst>), DeError> {
    let full = |p: &Edn<'_>, wrap: fn(MulticenterBondAst) -> EntitySpan<MulticenterBondAst>| {
        let e = parse_multicenter_bond_entry(p)?;
        Ok::<_, DeError>((e.id, e.atoms, wrap(e.bond.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "multicenter span map",
                    got: p.kind(),
                    path: vec!["multicenter span".into()],
                });
            };
            let atoms = parse_vec(
                required_key(m, "atoms", "multicenter span")?,
                ":atoms",
                |e| AtomRef::from_edn(e),
            )?;
            let (left, right) = pair(
                required_key(m, "type", "multicenter span")?,
                "multicenter span :modify",
            )?;
            Ok((
                optional_id(m)?,
                atoms,
                EntitySpan::Modified {
                    lhs: MulticenterBondDsl::from_edn(left)?.0,
                    rhs: MulticenterBondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "multicenter span: unexpected verb :{verb}"
        ))),
    }
}

// `|e| AtomRef::from_edn(e)` reads redundant but is load-bearing: the bare fn item can't satisfy the
// higher-ranked `for<'a> Fn(&Edn<'a>)` bound `parse_vec` requires.
#[allow(clippy::type_complexity, clippy::redundant_closure)]
fn parse_noncovalent_span_entry(
    edn: &Edn<'_>,
) -> Result<(Option<String>, [AtomRef; 2], EntitySpan<NoncovalentBondAst>), DeError> {
    let full = |p: &Edn<'_>, wrap: fn(NoncovalentBondAst) -> EntitySpan<NoncovalentBondAst>| {
        let e = parse_noncovalent_bond_entry(p)?;
        Ok::<_, DeError>((e.id, [e.first, e.second], wrap(e.bond.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "noncovalent span map",
                    got: p.kind(),
                    path: vec!["noncovalent span".into()],
                });
            };
            let atoms = two_atom_refs(
                parse_vec(
                    required_key(m, "atoms", "noncovalent span")?,
                    ":atoms",
                    |e| AtomRef::from_edn(e),
                )?,
                "noncovalent span",
            )?;
            let (left, right) = pair(
                required_key(m, "type", "noncovalent span")?,
                "noncovalent span :modify",
            )?;
            Ok((
                optional_id(m)?,
                atoms,
                EntitySpan::Modified {
                    lhs: NoncovalentBondDsl::from_edn(left)?.0,
                    rhs: NoncovalentBondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "noncovalent span: unexpected verb :{verb}"
        ))),
    }
}

#[allow(clippy::type_complexity)]
fn parse_stereo_atom_span_entry(
    edn: &Edn<'_>,
) -> Result<
    (
        Option<String>,
        AtomRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoAtomAst>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(StereoAtomAst) -> EntitySpan<StereoAtomAst>| {
        let e = parse_stereo_atom_entry(p)?;
        Ok::<_, DeError>((e.id, e.site, e.ligands, wrap(e.stereo.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "stereo-atom span map",
                    got: p.kind(),
                    path: vec!["stereo-atom span".into()],
                });
            };
            let site = AtomRef::from_edn(required_key(m, "site", "stereo-atom span")?)?;
            let ligands = parse_vec(
                required_key(m, "ligands", "stereo-atom span")?,
                ":ligands",
                parse_stereo_ligand,
            )?;
            let (left, right) = pair(
                required_key(m, "type", "stereo-atom span")?,
                "stereo-atom span :modify",
            )?;
            Ok((
                optional_id(m)?,
                site,
                ligands,
                EntitySpan::Modified {
                    lhs: StereoAtomDsl::from_edn(left)?.0,
                    rhs: StereoAtomDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "stereo-atom span: unexpected verb :{verb}"
        ))),
    }
}

#[allow(clippy::type_complexity)]
fn parse_stereo_bond_span_entry(
    edn: &Edn<'_>,
) -> Result<
    (
        Option<String>,
        BondRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoBondAst>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(StereoBondAst) -> EntitySpan<StereoBondAst>| {
        let e = parse_stereo_bond_entry(p)?;
        Ok::<_, DeError>((e.id, e.site, e.ligands, wrap(e.stereo.0)))
    };
    match verb_wrapper(edn) {
        None => full(edn, EntitySpan::Unchanged),
        Some(("add", p)) => full(p, EntitySpan::Added),
        Some(("remove", p)) => full(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let Edn::Map(m) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "stereo-bond span map",
                    got: p.kind(),
                    path: vec!["stereo-bond span".into()],
                });
            };
            let site = BondRef::from_edn(required_key(m, "site", "stereo-bond span")?)?;
            let ligands = parse_vec(
                required_key(m, "ligands", "stereo-bond span")?,
                ":ligands",
                parse_stereo_ligand,
            )?;
            let (left, right) = pair(
                required_key(m, "type", "stereo-bond span")?,
                "stereo-bond span :modify",
            )?;
            Ok((
                optional_id(m)?,
                site,
                ligands,
                EntitySpan::Modified {
                    lhs: StereoBondDsl::from_edn(left)?.0,
                    rhs: StereoBondDsl::from_edn(right)?.0,
                },
            ))
        }
        Some((verb, _)) => Err(DeError::Custom(format!(
            "stereo-bond span: unexpected verb :{verb}"
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
    let mut dative_bonds = Vec::new();
    let mut aromatic_systems = Vec::new();
    let mut multicenter_bonds = Vec::new();
    let mut noncovalent_bonds = Vec::new();
    let mut stereo_atoms = Vec::new();
    let mut stereo_bonds = Vec::new();
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
            "dative-bonds" => {
                dative_bonds = parse_vec(value, ":dative-bonds", parse_dative_span_entry)?
            }
            "aromatic-systems" => {
                aromatic_systems = parse_vec(value, ":aromatic-systems", parse_aromatic_span_entry)?
            }
            "multicenter-bonds" => {
                multicenter_bonds =
                    parse_vec(value, ":multicenter-bonds", parse_multicenter_span_entry)?
            }
            "noncovalent-bonds" => {
                noncovalent_bonds =
                    parse_vec(value, ":noncovalent-bonds", parse_noncovalent_span_entry)?
            }
            "stereo-atoms" => {
                stereo_atoms = parse_vec(value, ":stereo-atoms", parse_stereo_atom_span_entry)?
            }
            "stereo-bonds" => {
                stereo_bonds = parse_vec(value, ":stereo-bonds", parse_stereo_bond_span_entry)?
            }
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
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        stereo_atoms,
        stereo_bonds,
        constraints,
        atom_aliases,
    })
}

/// Streaming parse of a span map. The span entry grammar is tree-only (it reuses the molecule entry
/// parsers), so each section element is buffered to an `Edn` and dispatched to the tree entry parser.
fn read_span_input(de: &mut EdnStreamDeserializer<'_>) -> Result<SpanInput, EdnError> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut dative_bonds = Vec::new();
    let mut aromatic_systems = Vec::new();
    let mut multicenter_bonds = Vec::new();
    let mut noncovalent_bonds = Vec::new();
    let mut stereo_atoms = Vec::new();
    let mut stereo_bonds = Vec::new();
    let mut constraints = Vec::new();
    let mut atom_aliases = Vec::new();
    read_map(de, |de, key| {
        match key {
            "atoms" => {
                atoms = read_vec(de, |de| {
                    Ok(parse_atom_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "bonds" => {
                bonds = read_vec(de, |de| {
                    Ok(parse_bond_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "dative-bonds" => {
                dative_bonds = read_vec(de, |de| {
                    Ok(parse_dative_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "aromatic-systems" => {
                aromatic_systems = read_vec(de, |de| {
                    Ok(parse_aromatic_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "multicenter-bonds" => {
                multicenter_bonds = read_vec(de, |de| {
                    Ok(parse_multicenter_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "noncovalent-bonds" => {
                noncovalent_bonds = read_vec(de, |de| {
                    Ok(parse_noncovalent_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "stereo-atoms" => {
                stereo_atoms = read_vec(de, |de| {
                    Ok(parse_stereo_atom_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "stereo-bonds" => {
                stereo_bonds = read_vec(de, |de| {
                    Ok(parse_stereo_bond_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
                })?
            }
            "constraints" => {
                constraints = read_vec(de, |de| {
                    Ok(parse_constraint_span_entry(&read_string(
                        de.read_value_slice()?,
                    )?)?)
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
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        stereo_atoms,
        stereo_bonds,
        constraints,
        atom_aliases,
    })
}

fn resolve_atom_span(
    span: EntitySpan<AtomSpecInput>,
    namespace: &MoleculeNamespace,
) -> Result<EntitySpan<AtomAst>, ParseError> {
    Ok(match span {
        EntitySpan::Unchanged(s) => EntitySpan::Unchanged(resolve_atom_spec(s, namespace)?),
        EntitySpan::Added(s) => EntitySpan::Added(resolve_atom_spec(s, namespace)?),
        EntitySpan::Removed(s) => EntitySpan::Removed(resolve_atom_spec(s, namespace)?),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => EntitySpan::Modified {
            lhs: resolve_atom_spec(left, namespace)?,
            rhs: resolve_atom_spec(right, namespace)?,
        },
    })
}

fn resolve_constraint_span(
    input: ConstraintSpanInput,
    namespace: &MoleculeNamespace,
) -> Result<ConstraintSpan, ParseError> {
    Ok(match input {
        ConstraintSpanInput::Unchanged(dsl) => ConstraintSpan::Unchanged(dsl.into_ast(namespace)?),
        ConstraintSpanInput::Added(dsl) => ConstraintSpan::Added(dsl.into_ast(namespace)?),
        ConstraintSpanInput::Removed(dsl) => ConstraintSpan::Removed(dsl.into_ast(namespace)?),
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

        // The span's namespace: atoms take the union positions as their ids, then the bijective
        // aliases. Every ref resolves against it; `register_*` enforces id-disjointness, and the
        // roundtrip `MoleculeMetadata` is projected from it at the end.
        let mut namespace = MoleculeNamespace::default();
        for (id, _) in self.atoms.iter() {
            namespace.register_atom(id.clone())?;
        }
        for (name, dsl) in self.atom_aliases {
            namespace.register_atom_alias(name, dsl)?;
        }

        // Resolve atoms (alias → AtomAst), bonds (endpoints + value), constraints.
        let mut atoms: Vec<EntitySpan<AtomAst>> = Vec::with_capacity(atom_count);
        for (_, span) in self.atoms {
            atoms.push(resolve_atom_span(span, &namespace)?);
        }
        let mut bonds: Vec<EntitySpan<BondAst>> = Vec::with_capacity(bond_count);
        let mut endpoints: Vec<[AtomId; 2]> = Vec::with_capacity(bond_count);
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count);
        for (id_name, [ref_a, ref_b], span) in self.bonds {
            let a = ref_a.resolve(&namespace)?;
            let b = ref_b.resolve(&namespace)?;
            namespace.register_bond(id_name, a, b)?;
            edges.push([a.index() as u32, b.index() as u32]);
            endpoints.push([a, b]);
            bonds.push(span);
        }

        // Per-side ref consistency: a bond present on a side needs both endpoints present there.
        for (span, [a, b]) in bonds.iter().zip(&endpoints) {
            if span.lhs().is_some()
                && (atoms[a.index()].lhs().is_none() || atoms[b.index()].lhs().is_none())
            {
                return Err(ParseError::InvalidValue(
                    "bond present on the left references an atom absent on the left".into(),
                ));
            }
            if span.rhs().is_some()
                && (atoms[a.index()].rhs().is_none() || atoms[b.index()].rhs().is_none())
            {
                return Err(ParseError::InvalidValue(
                    "bond present on the right references an atom absent on the right".into(),
                ));
            }
        }

        // Overlay per-side ref consistency: an overlay present on a side needs every participant
        // atom present on that side (the bond check above, generalized to overlay participants).
        let side_ok =
            |lhs: bool, rhs: bool, participants: &[AtomId], kind: &str| -> Result<(), ParseError> {
                if lhs
                    && participants
                        .iter()
                        .any(|a| atoms[a.index()].lhs().is_none())
                {
                    return Err(ParseError::InvalidValue(format!(
                        "{kind} present on the left references an atom absent on the left"
                    )));
                }
                if rhs
                    && participants
                        .iter()
                        .any(|a| atoms[a.index()].rhs().is_none())
                {
                    return Err(ParseError::InvalidValue(format!(
                        "{kind} present on the right references an atom absent on the right"
                    )));
                }
                Ok(())
            };

        // Dative bonds: `[acceptor]` fixed side, donors var side.
        let mut dative_entries: Vec<([NodeId; 1], Vec<NodeId>, EntitySpan<DativeBondAst>)> =
            Vec::with_capacity(self.dative_bonds.len());
        for (id, donors, acceptor, span) in self.dative_bonds {
            let acceptor_id = acceptor.resolve(&namespace)?;
            let donor_ids: Vec<AtomId> = donors
                .into_iter()
                .map(|d| d.resolve(&namespace))
                .collect::<Result<_, _>>()?;
            let mut participants = Vec::with_capacity(donor_ids.len() + 1);
            participants.push(acceptor_id);
            participants.extend(donor_ids.iter().copied());
            side_ok(
                span.lhs().is_some(),
                span.rhs().is_some(),
                &participants,
                "dative bond",
            )?;
            namespace.register_dative_bond(id, &donor_ids, acceptor_id)?;
            dative_entries.push((
                [NodeId::from(acceptor_id)],
                donor_ids.iter().map(|&a| NodeId::from(a)).collect(),
                span,
            ));
        }
        let dative_bonds = FixedVarBirelationSet::new(dative_entries);

        // Aromatic systems.
        let mut aromatic_entries: Vec<(Vec<NodeId>, EntitySpan<AromaticSystemAst>)> =
            Vec::with_capacity(self.aromatic_systems.len());
        for (id, atoms_ref, span) in self.aromatic_systems {
            let atom_ids: Vec<AtomId> = atoms_ref
                .into_iter()
                .map(|r| r.resolve(&namespace))
                .collect::<Result<_, _>>()?;
            side_ok(
                span.lhs().is_some(),
                span.rhs().is_some(),
                &atom_ids,
                "aromatic system",
            )?;
            namespace.register_aromatic_system(id, &atom_ids)?;
            aromatic_entries.push((atom_ids.iter().map(|&a| NodeId::from(a)).collect(), span));
        }
        let aromatic_systems = VarRelationSet::new(aromatic_entries);

        // Multicenter bonds.
        let mut multicenter_entries: Vec<(Vec<NodeId>, EntitySpan<MulticenterBondAst>)> =
            Vec::with_capacity(self.multicenter_bonds.len());
        for (id, atoms_ref, span) in self.multicenter_bonds {
            let atom_ids: Vec<AtomId> = atoms_ref
                .into_iter()
                .map(|r| r.resolve(&namespace))
                .collect::<Result<_, _>>()?;
            side_ok(
                span.lhs().is_some(),
                span.rhs().is_some(),
                &atom_ids,
                "multicenter bond",
            )?;
            namespace.register_multicenter_bond(id, &atom_ids)?;
            multicenter_entries.push((atom_ids.iter().map(|&a| NodeId::from(a)).collect(), span));
        }
        let multicenter_bonds = VarRelationSet::new(multicenter_entries);

        // Noncovalent bonds.
        let mut noncovalent_entries: Vec<([NodeId; 2], EntitySpan<NoncovalentBondAst>)> =
            Vec::with_capacity(self.noncovalent_bonds.len());
        for (id, [first, second], span) in self.noncovalent_bonds {
            let a = first.resolve(&namespace)?;
            let b = second.resolve(&namespace)?;
            side_ok(
                span.lhs().is_some(),
                span.rhs().is_some(),
                &[a, b],
                "noncovalent bond",
            )?;
            namespace.register_noncovalent_bond(id, a, b)?;
            noncovalent_entries.push(([NodeId::from(a), NodeId::from(b)], span));
        }
        let noncovalent_bonds = FixedRelationSet::new(noncovalent_entries);

        // Stereo atoms: `[site]` fixed side, ligand frame var side.
        let mut stereo_atom_entries: Vec<(
            [NodeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoAtomAst>,
        )> = Vec::with_capacity(self.stereo_atoms.len());
        for (id, site, ligands, span) in self.stereo_atoms {
            let site_id = site.resolve(&namespace)?;
            let mut participants = vec![site_id];
            let mut ligand_frame = Vec::with_capacity(ligands.len());
            for l in ligands {
                let a = l.atom.resolve(&namespace)?;
                participants.push(a);
                ligand_frame.push(StereoLigand::new(a, l.kind));
            }
            side_ok(
                span.lhs().is_some(),
                span.rhs().is_some(),
                &participants,
                "stereo atom",
            )?;
            namespace.register_stereo_atom(id, site_id, &ligand_frame)?;
            stereo_atom_entries.push(([NodeId::from(site_id)], ligand_frame, span));
        }
        let stereo_atoms = FixedVarBirelationSet::new(stereo_atom_entries);

        // Stereo bonds: `[site]` fixed bond side, ligand frame var side. The site bond and every
        // ligand atom must be present on any side the stereo bond is present.
        let mut stereo_bond_entries: Vec<(
            [EdgeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoBondAst>,
        )> = Vec::with_capacity(self.stereo_bonds.len());
        for (id, site, ligands, span) in self.stereo_bonds {
            let site_id = site.resolve(&namespace)?;
            let mut ligand_atoms = Vec::with_capacity(ligands.len());
            let mut ligand_frame = Vec::with_capacity(ligands.len());
            for l in ligands {
                let a = l.atom.resolve(&namespace)?;
                ligand_atoms.push(a);
                ligand_frame.push(StereoLigand::new(a, l.kind));
            }
            if span.lhs().is_some()
                && (bonds[site_id.index()].lhs().is_none()
                    || ligand_atoms
                        .iter()
                        .any(|a| atoms[a.index()].lhs().is_none()))
            {
                return Err(ParseError::InvalidValue(
                    "stereo bond present on the left references a bond or atom absent on the left"
                        .into(),
                ));
            }
            if span.rhs().is_some()
                && (bonds[site_id.index()].rhs().is_none()
                    || ligand_atoms
                        .iter()
                        .any(|a| atoms[a.index()].rhs().is_none()))
            {
                return Err(ParseError::InvalidValue(
                    "stereo bond present on the right references a bond or atom absent on the right"
                        .into(),
                ));
            }
            namespace.register_stereo_bond(id, site_id, &ligand_frame)?;
            stereo_bond_entries.push(([EdgeId::from(site_id)], ligand_frame, span));
        }
        let stereo_bonds = FixedVarBirelationSet::new(stereo_bond_entries);

        let mut constraints: Vec<ConstraintSpan> = Vec::with_capacity(self.constraints.len());
        for input in self.constraints {
            constraints.push(resolve_constraint_span(input, &namespace)?);
        }

        let graph = Graph::new(atom_count, &edges);
        let metadata = MoleculeMetadata::from(&namespace);
        Ok((
            ReactionSpanAst::from_parts(
                graph,
                atoms,
                bonds,
                dative_bonds,
                aromatic_systems,
                multicenter_bonds,
                noncovalent_bonds,
                stereo_atoms,
                stereo_bonds,
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
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            Edn::Vector(
                vec![
                    render_atom_value(left, meta),
                    render_atom_value(right, meta),
                ]
                .into(),
            ),
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
        EntitySpan::Added(b) => {
            single_key_map("add", render_bond_entry(id, endpoints, value(b), meta))
        }
        EntitySpan::Removed(b) => {
            single_key_map("remove", render_bond_entry(id, endpoints, value(b), meta))
        }
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
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

fn render_dative_span_entry(
    id: DativeBondId,
    donors: &[AtomId],
    acceptor: AtomId,
    span: &EntitySpan<DativeBondAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |b: &DativeBondAst| DativeBondDsl::from_ref(b).to_edn();
    let entry = |type_edn: Edn<'static>| {
        render_dative_entry(id, donors.iter().copied(), acceptor, type_edn, meta)
    };
    match span {
        EntitySpan::Unchanged(b) => entry(value(b)),
        EntitySpan::Added(b) => single_key_map("add", entry(value(b))),
        EntitySpan::Removed(b) => single_key_map("remove", entry(value(b))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
        ),
    }
}

fn render_aromatic_span_entry(
    id: AromaticSystemId,
    atoms: &[AtomId],
    span: &EntitySpan<AromaticSystemAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value =
        |s: &AromaticSystemAst| Edn::Str(Cow::Owned(AromaticSystemDsl::from_ref(s).to_string()));
    let entry =
        |type_edn: Edn<'static>| render_aromatic_entry(id, atoms.iter().copied(), type_edn, meta);
    match span {
        EntitySpan::Unchanged(s) => entry(value(s)),
        EntitySpan::Added(s) => single_key_map("add", entry(value(s))),
        EntitySpan::Removed(s) => single_key_map("remove", entry(value(s))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
        ),
    }
}

fn render_multicenter_span_entry(
    id: MulticenterBondId,
    atoms: &[AtomId],
    span: &EntitySpan<MulticenterBondAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value =
        |b: &MulticenterBondAst| Edn::Str(Cow::Owned(MulticenterBondDsl::from_ref(b).to_string()));
    let entry = |type_edn: Edn<'static>| {
        render_multicenter_entry(id, atoms.iter().copied(), type_edn, meta)
    };
    match span {
        EntitySpan::Unchanged(b) => entry(value(b)),
        EntitySpan::Added(b) => single_key_map("add", entry(value(b))),
        EntitySpan::Removed(b) => single_key_map("remove", entry(value(b))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
        ),
    }
}

fn render_noncovalent_span_entry(
    id: NoncovalentBondId,
    atoms: [AtomId; 2],
    span: &EntitySpan<NoncovalentBondAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |b: &NoncovalentBondAst| NoncovalentBondDsl::from_ref(b).to_edn();
    let entry = |type_edn: Edn<'static>| render_noncovalent_entry(id, atoms, type_edn, meta);
    match span {
        EntitySpan::Unchanged(b) => entry(value(b)),
        EntitySpan::Added(b) => single_key_map("add", entry(value(b))),
        EntitySpan::Removed(b) => single_key_map("remove", entry(value(b))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
        ),
    }
}

fn render_stereo_atom_span_entry(
    id: StereoAtomId,
    site: AtomId,
    ligands: &[StereoLigand],
    span: &EntitySpan<StereoAtomAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |s: &StereoAtomAst| StereoAtomDsl::from_ref(s).to_edn();
    let ligand_edns: Vec<Edn<'static>> = ligands
        .iter()
        .map(|&l| render_stereo_ligand(l, meta))
        .collect();
    let entry = |type_edn: Edn<'static>| {
        render_stereo_atom_entry(id, site, ligand_edns.clone(), type_edn, meta)
    };
    match span {
        EntitySpan::Unchanged(s) => entry(value(s)),
        EntitySpan::Added(s) => single_key_map("add", entry(value(s))),
        EntitySpan::Removed(s) => single_key_map("remove", entry(value(s))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
        ),
    }
}

fn render_stereo_bond_span_entry(
    id: StereoBondId,
    site: BondId,
    ligands: &[StereoLigand],
    span: &EntitySpan<StereoBondAst>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |s: &StereoBondAst| StereoBondDsl::from_ref(s).to_edn();
    let ligand_edns: Vec<Edn<'static>> = ligands
        .iter()
        .map(|&l| render_stereo_ligand(l, meta))
        .collect();
    let entry = |type_edn: Edn<'static>| {
        render_stereo_bond_entry(id, site, ligand_edns.clone(), type_edn, meta)
    };
    match span {
        EntitySpan::Unchanged(s) => entry(value(s)),
        EntitySpan::Added(s) => single_key_map("add", entry(value(s))),
        EntitySpan::Removed(s) => single_key_map("remove", entry(value(s))),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => single_key_map(
            "modify",
            entry(Edn::Vector(vec![value(left), value(right)].into())),
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
            render_bond_span_entry(
                BondId(j as u32),
                [AtomId::from(a), AtomId::from(b)],
                span,
                meta,
            )
        })
        .collect();
    if !bonds.is_empty() {
        map.insert(Edn::keyword("bonds"), Edn::Vector(bonds.into()));
    }

    let dative = ast.dative_bonds();
    let dative_bonds: Vec<Edn<'static>> = dative
        .relation_ids()
        .map(|rid| {
            let acceptor = AtomId::from(dative.participants_1(rid)[0]);
            let donors: Vec<AtomId> = dative
                .participants_2(rid)
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            render_dative_span_entry(
                DativeBondId(rid.index() as u32),
                &donors,
                acceptor,
                dative.data(rid),
                meta,
            )
        })
        .collect();
    if !dative_bonds.is_empty() {
        map.insert(
            Edn::keyword("dative-bonds"),
            Edn::Vector(dative_bonds.into()),
        );
    }

    let aromatic = ast.aromatic_systems();
    let aromatic_systems: Vec<Edn<'static>> = aromatic
        .relation_ids()
        .map(|rid| {
            let atoms: Vec<AtomId> = aromatic
                .participants(rid)
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            render_aromatic_span_entry(
                AromaticSystemId(rid.index() as u32),
                &atoms,
                aromatic.data(rid),
                meta,
            )
        })
        .collect();
    if !aromatic_systems.is_empty() {
        map.insert(
            Edn::keyword("aromatic-systems"),
            Edn::Vector(aromatic_systems.into()),
        );
    }

    let multicenter = ast.multicenter_bonds();
    let multicenter_bonds: Vec<Edn<'static>> = multicenter
        .relation_ids()
        .map(|rid| {
            let atoms: Vec<AtomId> = multicenter
                .participants(rid)
                .iter()
                .map(|&n| AtomId::from(n))
                .collect();
            render_multicenter_span_entry(
                MulticenterBondId(rid.index() as u32),
                &atoms,
                multicenter.data(rid),
                meta,
            )
        })
        .collect();
    if !multicenter_bonds.is_empty() {
        map.insert(
            Edn::keyword("multicenter-bonds"),
            Edn::Vector(multicenter_bonds.into()),
        );
    }

    let noncovalent = ast.noncovalent_bonds();
    let noncovalent_bonds: Vec<Edn<'static>> = noncovalent
        .relation_ids()
        .map(|rid| {
            let [a, b] = noncovalent.participants(rid);
            render_noncovalent_span_entry(
                NoncovalentBondId(rid.index() as u32),
                [AtomId::from(*a), AtomId::from(*b)],
                noncovalent.data(rid),
                meta,
            )
        })
        .collect();
    if !noncovalent_bonds.is_empty() {
        map.insert(
            Edn::keyword("noncovalent-bonds"),
            Edn::Vector(noncovalent_bonds.into()),
        );
    }

    let stereo_a = ast.stereo_atoms();
    let stereo_atoms: Vec<Edn<'static>> = stereo_a
        .relation_ids()
        .map(|rid| {
            let site = AtomId::from(stereo_a.participants_1(rid)[0]);
            render_stereo_atom_span_entry(
                StereoAtomId(rid.index() as u32),
                site,
                stereo_a.participants_2(rid),
                stereo_a.data(rid),
                meta,
            )
        })
        .collect();
    if !stereo_atoms.is_empty() {
        map.insert(
            Edn::keyword("stereo-atoms"),
            Edn::Vector(stereo_atoms.into()),
        );
    }

    let stereo_b = ast.stereo_bonds();
    let stereo_bonds: Vec<Edn<'static>> = stereo_b
        .relation_ids()
        .map(|rid| {
            let site = BondId::from(stereo_b.participants_1(rid)[0]);
            render_stereo_bond_span_entry(
                StereoBondId(rid.index() as u32),
                site,
                stereo_b.participants_2(rid),
                stereo_b.data(rid),
                meta,
            )
        })
        .collect();
    if !stereo_bonds.is_empty() {
        map.insert(
            Edn::keyword("stereo-bonds"),
            Edn::Vector(stereo_bonds.into()),
        );
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
    use std::fs;

    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::boolean::BooleanAst;
    use crate::ast::constraint::{BondConstraint, Constraint, Constraints, MoleculeConstraint};
    use crate::ast::delta::{AtomDelta, BondDelta, ConstraintDelta, Delta, Deltas};
    use crate::ast::edit::BondFieldChange;
    use crate::ast::ligand::StereoLigandKind;
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
        lhs:AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
        rhs:AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::N)))),
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
        lhs:BondAst::from_order(1),
        rhs:BondAst::from_order(2),
    }))]
    #[case::unchanged_map_id("{:id :b1 :atoms [0 1] :type :single}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondAst::from_order(1))))]
    #[case::modify_map_id("{:modify {:id :b1 :atoms [0 1] :type [:single :double]}}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Modified {
        lhs:BondAst::from_order(1),
        rhs:BondAst::from_order(2),
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
    fn test_parse_constraint_span_entry(
        #[case] input: &str,
        #[case] expected: ConstraintSpanInput,
    ) {
        assert_eq!(
            parse_constraint_span_entry(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:donors [0] :acceptor 1 :type "1#R"}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::add(r#"{:add {:donors [0] :acceptor 1 :type "1#R"}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Added(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:donors [0] :acceptor 1 :type "1#R"}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Removed(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:donors [0] :acceptor 1 :type ["1#R" "2#R"]}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Modified {
            lhs:DativeBondDsl::from_str("1#R").unwrap().0,
            rhs:DativeBondDsl::from_str("2#R").unwrap().0,
        },
    ))]
    #[case::with_id(r#"{:id :d1 :donors [0 2] :acceptor 1 :type "1#R"}"#, (
        Some("d1".to_string()), vec![AtomRef::Index(0), AtomRef::Index(2)], AtomRef::Index(1),
        EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    fn test_parse_dative_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            Vec<AtomRef>,
            AtomRef,
            EntitySpan<DativeBondAst>,
        ),
    ) {
        assert_eq!(
            parse_dative_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1 2] :type "*#e6"}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1), AtomRef::Index(2)],
        EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0),
    ))]
    #[case::add(r#"{:add {:atoms [0 1] :type "*#e2"}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Added(AromaticSystemDsl::from_str("*#e2").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:atoms [0 1] :type ["*#e2" "*#e4"]}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Modified {
            lhs:AromaticSystemDsl::from_str("*#e2").unwrap().0,
            rhs:AromaticSystemDsl::from_str("*#e4").unwrap().0,
        },
    ))]
    fn test_parse_aromatic_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemAst>),
    ) {
        assert_eq!(
            parse_aromatic_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1] :type "*#e2"}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Unchanged(MulticenterBondDsl::from_str("*#e2").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:atoms [0 1 2] :type "*#e3"}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1), AtomRef::Index(2)],
        EntitySpan::Removed(MulticenterBondDsl::from_str("*#e3").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:atoms [0 1] :type ["*#e2" "*#e4"]}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Modified {
            lhs:MulticenterBondDsl::from_str("*#e2").unwrap().0,
            rhs:MulticenterBondDsl::from_str("*#e4").unwrap().0,
        },
    ))]
    fn test_parse_multicenter_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, Vec<AtomRef>, EntitySpan<MulticenterBondAst>),
    ) {
        assert_eq!(
            parse_multicenter_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1] :type "Hbd"}"#, (
        None, [AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:atoms [0 1] :type "Hbd"}}"#, (
        None, [AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Removed(NoncovalentBondDsl::from_str("Hbd").unwrap().0),
    ))]
    fn test_parse_noncovalent_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, [AtomRef; 2], EntitySpan<NoncovalentBondAst>),
    ) {
        assert_eq!(
            parse_noncovalent_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:site 0 :ligands [1 2 3 4] :type "Th1"}"#, (
        None, AtomRef::Index(0),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(1) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(2) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(4) },
        ],
        EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0),
    ))]
    #[case::add(r#"{:add {:site 0 :ligands [1 2 [:h 3]] :type "Th1"}}"#, (
        None, AtomRef::Index(0),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(1) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(2) },
            StereoLigandRef { kind: StereoLigandKind::ImplicitHydrogen, atom: AtomRef::Index(3) },
        ],
        EntitySpan::Added(StereoAtomDsl::from_str("Th1").unwrap().0),
    ))]
    fn test_parse_stereo_atom_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            AtomRef,
            Vec<StereoLigandRef>,
            EntitySpan<StereoAtomAst>,
        ),
    ) {
        assert_eq!(
            parse_stereo_atom_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:site 1 :ligands [0 3] :type "Ct1"}"#, (
        None, BondRef::Index(1),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(0) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
        ],
        EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:site 1 :ligands [0 3] :type "Ct1"}}"#, (
        None, BondRef::Index(1),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(0) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
        ],
        EntitySpan::Removed(StereoBondDsl::from_str("Ct1").unwrap().0),
    ))]
    fn test_parse_stereo_bond_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            BondRef,
            Vec<StereoLigandRef>,
            EntitySpan<StereoBondAst>,
        ),
    ) {
        assert_eq!(
            parse_stereo_bond_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
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
            dative_bonds: vec![],
            aromatic_systems: vec![],
            multicenter_bonds: vec![],
            noncovalent_bonds: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
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
            dative_bonds: vec![],
            aromatic_systems: vec![],
            multicenter_bonds: vec![],
            noncovalent_bonds: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
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
            dative_bonds: vec![],
            aromatic_systems: vec![],
            multicenter_bonds: vec![],
            noncovalent_bonds: vec![],
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            constraints: vec![],
            atom_aliases: vec![(
                "nu".to_string(),
                Box::new(AtomDsl::from_edn(&read_string(r#""C#h3""#).unwrap()).unwrap()),
            )],
        },
    )]
    fn test_parse_span_input(#[case] input: &str, #[case] expected: SpanInput) {
        assert_eq!(
            parse_span_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
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
    #[case::dative_overlay(
        r#"{:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1#R"}]}"#,
        ReactionSpanAst::from_parts(
            Graph::new(2, &[]),
            vec![
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Unchanged(AtomAst::from_element(Element::N)),
            ],
            vec![],
            FixedVarBirelationSet::new(vec![(
                [NodeId(1)],
                vec![NodeId(0)],
                EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
            )]),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            vec![],
        ),
        MoleculeMetadata::new().with_dative_bond_id(DativeBondId(0), "d1"),
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
    fn test_span_input_into_ast_structural_bond_ref() {
        // A structural bond ref ({:atoms [0 1]}) names the bond by its endpoints, resolved against
        // the namespace's participant lookup.
        let input = r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#;
        let (ast, _) = parse_span_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.constraints().to_vec(),
            vec![ConstraintSpan::Unchanged(Constraint::Bond(
                BondId(0),
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
            ))]
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
    #[case::duplicate_id(
        r#"{:atoms [[:a "C"] [:a "O"]]}"#,
        ParseError::DuplicateId("a".to_string()),
    )]
    #[case::left_inconsistent(
        r#"{:atoms ["C" {:add "O"}] :bonds [[0 1 :single]]}"#,
        ParseError::InvalidValue(
            "bond present on the left references an atom absent on the left".to_string(),
        ),
    )]
    #[case::overlay_left_inconsistent(
        r#"{:atoms ["C" {:add "N"}] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]}"#,
        ParseError::InvalidValue(
            "dative bond present on the left references an atom absent on the left".to_string(),
        ),
    )]
    #[case::stereo_bond_left_inconsistent(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:add [1 2 "2"]} [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#,
        ParseError::InvalidValue(
            "stereo bond present on the left references a bond or atom absent on the left".to_string(),
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
    #[case::unchanged(
        AtomId(0),
        EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
        MoleculeMetadata::new(),
        r#""C""#
    )]
    #[case::add(
        AtomId(0),
        EntitySpan::Added(AtomAst::from_element(Element::O)),
        MoleculeMetadata::new(),
        r#"{:add "O"}"#
    )]
    #[case::remove(
        AtomId(0),
        EntitySpan::Removed(AtomAst::from_element(Element::O)),
        MoleculeMetadata::new(),
        r#"{:remove "O"}"#
    )]
    #[case::modify(AtomId(0), EntitySpan::Modified { lhs:AtomAst::from_element(Element::C), rhs:AtomAst::from_element(Element::N) }, MoleculeMetadata::new(), r#"{:modify ["C" "N"]}"#)]
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
    #[case::unchanged(
        EntitySpan::Unchanged(BondAst::from_order(1)),
        MoleculeMetadata::new(),
        "[0 1 :single]"
    )]
    #[case::add(
        EntitySpan::Added(BondAst::from_order(2)),
        MoleculeMetadata::new(),
        "{:add [0 1 :double]}"
    )]
    #[case::remove(
        EntitySpan::Removed(BondAst::from_order(1)),
        MoleculeMetadata::new(),
        "{:remove [0 1 :single]}"
    )]
    #[case::modify(EntitySpan::Modified { lhs:BondAst::from_order(1), rhs:BondAst::from_order(2) }, MoleculeMetadata::new(), "{:modify [0 1 [:single :double]]}")]
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

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:donors [0] :acceptor 1 :type "1#R"}"#)]
    #[case::add(EntitySpan::Added(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:add {:donors [0] :acceptor 1 :type "1#R"}}"#)]
    #[case::remove(EntitySpan::Removed(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:donors [0] :acceptor 1 :type "1#R"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:DativeBondDsl::from_str("1#R").unwrap().0, rhs:DativeBondDsl::from_str("2#R").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:donors [0] :acceptor 1 :type ["1#R" "2#R"]}}"#)]
    #[case::with_id(EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new().with_dative_bond_id(DativeBondId(0), "d1"), r#"{:id :d1 :donors [0] :acceptor 1 :type "1#R"}"#)]
    fn test_render_dative_span_entry(
        #[case] span: EntitySpan<DativeBondAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_dative_span_entry(DativeBondId(0), &[AtomId(0)], AtomId(1), &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1 2] :type "*#e6"}"#)]
    #[case::add(EntitySpan::Added(AromaticSystemDsl::from_str("*#e6").unwrap().0), MoleculeMetadata::new(), r#"{:add {:atoms [0 1 2] :type "*#e6"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:AromaticSystemDsl::from_str("*#e6").unwrap().0, rhs:AromaticSystemDsl::from_str("*#e2").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:atoms [0 1 2] :type ["*#e6" "*#e2"]}}"#)]
    #[case::with_id(EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0), MoleculeMetadata::new().with_aromatic_system_id(AromaticSystemId(0), "ar1"), r#"{:id :ar1 :atoms [0 1 2] :type "*#e6"}"#)]
    fn test_render_aromatic_span_entry(
        #[case] span: EntitySpan<AromaticSystemAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_aromatic_span_entry(
                AromaticSystemId(0),
                &[AtomId(0), AtomId(1), AtomId(2)],
                &span,
                &meta
            ),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(MulticenterBondDsl::from_str("*#e2").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1] :type "*#e2"}"#)]
    #[case::remove(EntitySpan::Removed(MulticenterBondDsl::from_str("*#e2").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:atoms [0 1] :type "*#e2"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:MulticenterBondDsl::from_str("*#e2").unwrap().0, rhs:MulticenterBondDsl::from_str("*#e4").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:atoms [0 1] :type ["*#e2" "*#e4"]}}"#)]
    fn test_render_multicenter_span_entry(
        #[case] span: EntitySpan<MulticenterBondAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_multicenter_span_entry(
                MulticenterBondId(0),
                &[AtomId(0), AtomId(1)],
                &span,
                &meta
            ),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1] :type "Hbd"}"#)]
    #[case::add(EntitySpan::Added(NoncovalentBondDsl::from_str("Hbd").unwrap().0), MoleculeMetadata::new(), r#"{:add {:atoms [0 1] :type "Hbd"}}"#)]
    #[case::with_id(EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0), MoleculeMetadata::new().with_noncovalent_bond_id(NoncovalentBondId(0), "nc1"), r#"{:id :nc1 :atoms [0 1] :type "Hbd"}"#)]
    fn test_render_noncovalent_span_entry(
        #[case] span: EntitySpan<NoncovalentBondAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_noncovalent_span_entry(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(1)],
                &span,
                &meta
            ),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0), MoleculeMetadata::new(), r#"{:site 0 :ligands [1 2 3 4] :type :cw}"#)]
    #[case::add(EntitySpan::Added(StereoAtomDsl::from_str("Th1").unwrap().0), MoleculeMetadata::new(), r#"{:add {:site 0 :ligands [1 2 3 4] :type :cw}}"#)]
    #[case::with_id(EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0), MoleculeMetadata::new().with_stereo_atom_id(StereoAtomId(0), "s1"), r#"{:id :s1 :site 0 :ligands [1 2 3 4] :type :cw}"#)]
    fn test_render_stereo_atom_span_entry(
        #[case] span: EntitySpan<StereoAtomAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        assert_eq!(
            render_stereo_atom_span_entry(StereoAtomId(0), AtomId(0), &ligands, &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0), MoleculeMetadata::new(), r#"{:site 1 :ligands [0 3] :type :e}"#)]
    #[case::remove(EntitySpan::Removed(StereoBondDsl::from_str("Ct1").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:site 1 :ligands [0 3] :type :e}}"#)]
    #[case::with_id(EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0), MoleculeMetadata::new().with_stereo_bond_id(StereoBondId(0), "sb1"), r#"{:id :sb1 :site 1 :ligands [0 3] :type :e}"#)]
    fn test_render_stereo_bond_span_entry(
        #[case] span: EntitySpan<StereoBondAst>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
        assert_eq!(
            render_stereo_bond_span_entry(StereoBondId(0), BondId(1), &ligands, &span, &meta),
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
    #[case::dative(
        r#"{:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1#R"}]}"#
    )]
    #[case::dative_add(
        r#"{:atoms ["C" "N"] :dative-bonds [{:add {:donors [0] :acceptor 1 :type "1#R"}}]}"#
    )]
    #[case::dative_modify(r#"{:atoms ["C" "N"] :dative-bonds [{:modify {:donors [0] :acceptor 1 :type ["1#R" "2#R"]}}]}"#)]
    #[case::aromatic(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*#e6"}]}"#)]
    #[case::multicenter(r#"{:atoms ["C" "C"] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]}"#)]
    #[case::noncovalent(
        r#"{:atoms ["N" "H"] :noncovalent-bonds [{:remove {:atoms [0 1] :type "Hbd"}}]}"#
    )]
    #[case::stereo_atom(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#)]
    #[case::stereo_bond(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#)]
    fn test_reaction_span_dsl_to_edn(#[case] input: &str) {
        let dsl = ReactionSpanDsl::from_str(input).unwrap();
        assert_eq!(ReactionSpanDsl::from_edn(&dsl.to_edn()).unwrap(), dsl);
    }

    // `ReactionSpanAst` renders positionally and round-trips (metadata discarded on both sides).
    #[rstest]
    #[case::plain_molecule(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::modify(r#"{:atoms ["C" "C"] :bonds [{:modify [0 1 [:single :double]]}]}"#)]
    #[case::add_remove(r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :single]}]}"#)]
    #[case::dative(r#"{:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]}"#)]
    #[case::dative_modify(r#"{:atoms ["C" "N"] :dative-bonds [{:modify {:donors [0] :acceptor 1 :type ["1#R" "2#R"]}}]}"#)]
    #[case::aromatic(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*#e6"}]}"#)]
    #[case::multicenter(
        r#"{:atoms ["C" "C"] :multicenter-bonds [{:add {:atoms [0 1] :type "*#e2"}}]}"#
    )]
    #[case::noncovalent(r#"{:atoms ["N" "H"] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#)]
    #[case::stereo_atom(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#)]
    #[case::stereo_bond(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#)]
    fn test_reaction_span_ast_to_edn(#[case] input: &str) {
        let ast = ReactionSpanAst::from_str(input).unwrap();
        assert_eq!(ReactionSpanAst::from_edn(&ast.to_edn()).unwrap(), ast);
    }

    /// Every `fuzz_reaction_span` seed must parse and satisfy the fuzz invariant (streaming and tree
    /// parsers agree) — guards the seed corpus against rot as the span DSL evolves.
    #[rstest]
    fn test_fuzz_reaction_span_seeds_valid() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fuzz/seeds/fuzz_reaction_span");
        let mut count = 0;
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let data = fs::read_to_string(&path).unwrap();
            let stream = ReactionSpanDsl::from_edn_str(&data).ok();
            let tree = read_string(&data)
                .ok()
                .and_then(|edn| ReactionSpanDsl::from_edn(&edn).ok());
            assert!(stream.is_some(), "seed {path:?} failed to parse");
            assert_eq!(
                stream, tree,
                "seed {path:?}: streaming and tree parsers disagree"
            );
            count += 1;
        }
        assert_eq!(count, 21);
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
        #[case] lhs: Vec<Constraint>,
        #[case] rhs: Vec<Constraint>,
    ) {
        let ast = ReactionSpanAst::from_str(input).unwrap();
        assert_eq!(ast.lhs().constraints().as_slice(), lhs.as_slice());
        assert_eq!(ast.rhs().constraints().as_slice(), rhs.as_slice());
    }
}
