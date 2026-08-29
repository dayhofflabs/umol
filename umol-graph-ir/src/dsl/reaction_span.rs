//! Reaction span DSL: the surface form of `ReactionSpan`, where each entity carries its
//! complete before/after value (`EntitySpan`) rather than a delta. Entity ids, bond endpoints,
//! and constraint topology refs are resolved in `into_ir`.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{
    read_string, DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn,
};
use umol_graph_core::EdgeId;

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::config::MoleculeDefaults;
use super::constraint::ConstraintDsl;
use super::dative::DativeBondDsl;
use super::edn_utils::{
    optional_id_keyword, pair, parse_vec, read_map, read_vec, required_key, single_key_map,
    two_atom_refs,
};
use super::error::ParseError;
use super::metadata::{MetadataError, MoleculeMetadata};
use super::molecule::{
    parse_aromatic_system_entry, parse_atom_aliases, parse_atom_entry, parse_bond_entry,
    parse_dative_bond_entry, parse_multicenter_bond_entry, parse_noncovalent_bond_entry,
    parse_stereo_atom_entry, parse_stereo_bond_entry, read_atom_aliases, render_aromatic_entry,
    render_atom_value, render_bond_entry, render_dative_entry, render_multicenter_entry,
    render_noncovalent_entry, render_stereo_atom_entry, render_stereo_bond_entry,
    render_stereo_ligand, resolve_atom_spec, AtomSpecInput,
};
use super::multicenter::MulticenterBondDsl;
use super::namespace::MoleculeContext;
use super::noncovalent::NoncovalentBondDsl;
use super::refs::{parse_stereo_ligand, AtomRef, BondRef, StereoLigandRef};
use super::stereo::{StereoAtomDsl, StereoBondDsl};
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::entity::Entity;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::StereoLigand;
use crate::ir::traits::{FromIr, IntoIr};
use crate::ir::{
    AromaticSystemForm, Constraint, ConstraintSpan, DativeBondForm, EntitySpan,
    MulticenterBondForm, NoncovalentBondForm, ReactionSpan, ReactionSpanEntries, StereoAtomForm,
    StereoBondForm,
};

/// Surface DSL for a reaction span. Pairs `ReactionSpan` with the
/// `MoleculeMetadata` recording its span-frame entity-keyword bindings; fields
/// are private so metadata cannot drift onto a different reaction span.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanDsl {
    span: ReactionSpan,
    metadata: MoleculeMetadata,
}

impl ReactionSpanDsl {
    /// Pair a reaction span with coherent surface metadata.
    pub fn new(span: ReactionSpan, metadata: MoleculeMetadata) -> Result<Self, MetadataError> {
        for (entity, _) in metadata.iter_keywords() {
            let contains = match entity {
                Entity::Atom(id) => id.index() < span.atoms().len(),
                Entity::Bond(id) => id.index() < span.bonds().len(),
                Entity::DativeBond(id) => id.index() < span.dative_bonds().count(),
                Entity::AromaticSystem(id) => id.index() < span.aromatic_systems().count(),
                Entity::MulticenterBond(id) => id.index() < span.multicenter_bonds().count(),
                Entity::NoncovalentBond(id) => id.index() < span.noncovalent_bonds().count(),
                Entity::StereoAtom(id) => id.index() < span.stereo_atoms().count(),
                Entity::StereoBond(id) => id.index() < span.stereo_bonds().count(),
            };
            if !contains {
                return Err(MetadataError::EntityOutOfRange(entity));
            }
        }
        Ok(Self::from_parts(span, metadata))
    }

    fn from_parts(span: ReactionSpan, metadata: MoleculeMetadata) -> Self {
        Self { span, metadata }
    }

    pub fn span(&self) -> &ReactionSpan {
        &self.span
    }

    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionSpan, MoleculeMetadata) {
        (self.span, self.metadata)
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

impl FromIr<ReactionSpan> for ReactionSpanDsl {
    type Context = MoleculeDefaults;

    fn from_ir(span: &ReactionSpan, context: &Self::Context) -> Self {
        let lowered = ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: span
                .atoms()
                .iter()
                .map(|entity_span| {
                    map_span(entity_span, |atom| AtomDsl::from_ir(atom, &context.atom).0)
                })
                .collect(),
            bonds: span
                .bonds()
                .iter()
                .enumerate()
                .map(|(index, entity_span)| {
                    let [first, second] = span.graph().edge_endpoints(EdgeId(index as u32));
                    (
                        AtomId::from(first),
                        AtomId::from(second),
                        map_span(entity_span, |bond| BondDsl::from_ir(bond, &context.bond).0),
                    )
                })
                .collect(),
            dative: span
                .dative_bonds()
                .ids()
                .map(|id| {
                    (
                        span.dative_bonds().donors(id).collect(),
                        span.dative_bonds().acceptor(id),
                        map_span(span.dative_bonds().attributes(id), |bond| {
                            DativeBondDsl::from_ir(bond, &context.dative_bond).0
                        }),
                    )
                })
                .collect(),
            aromatic: span
                .aromatic_systems()
                .ids()
                .map(|id| {
                    (
                        span.aromatic_systems().atoms(id).collect(),
                        map_span(span.aromatic_systems().attributes(id), |system| {
                            AromaticSystemDsl::from_ir(system, &context.aromatic_system).0
                        }),
                    )
                })
                .collect(),
            multicenter: span
                .multicenter_bonds()
                .ids()
                .map(|id| {
                    (
                        span.multicenter_bonds().atoms(id).collect(),
                        map_span(span.multicenter_bonds().attributes(id), |bond| {
                            MulticenterBondDsl::from_ir(bond, &context.multicenter_bond).0
                        }),
                    )
                })
                .collect(),
            noncovalent: span
                .noncovalent_bonds()
                .ids()
                .map(|id| {
                    let [first, second] = span.noncovalent_bonds().atoms(id);
                    (
                        first,
                        second,
                        map_span(span.noncovalent_bonds().attributes(id), |bond| {
                            NoncovalentBondDsl::from_ir(bond, &context.noncovalent_bond).0
                        }),
                    )
                })
                .collect(),
            stereo_atoms: span
                .stereo_atoms()
                .ids()
                .map(|id| {
                    (
                        span.stereo_atoms().site(id),
                        span.stereo_atoms().ligands(id).to_vec(),
                        map_span(span.stereo_atoms().attributes(id), |stereo| {
                            StereoAtomDsl::from_ir(stereo, &context.stereo_atom).0
                        }),
                    )
                })
                .collect(),
            stereo_bonds: span
                .stereo_bonds()
                .ids()
                .map(|id| {
                    (
                        span.stereo_bonds().site(id),
                        span.stereo_bonds().ligands(id).to_vec(),
                        map_span(span.stereo_bonds().attributes(id), |stereo| {
                            StereoBondDsl::from_ir(stereo, &context.stereo_bond).0
                        }),
                    )
                })
                .collect(),
            constraints: span.constraints().to_vec(),
        });
        ReactionSpanDsl::from_parts(lowered, MoleculeMetadata::default())
    }
}

impl IntoIr<ReactionSpan> for ReactionSpanDsl {
    type Context = MoleculeDefaults;

    fn into_ir(self, context: &Self::Context) -> ReactionSpan {
        let span = self.span;
        ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: span
                .atoms()
                .iter()
                .map(|entity_span| {
                    map_span(entity_span, |atom| {
                        AtomDsl(atom.clone()).into_ir(&context.atom)
                    })
                })
                .collect(),
            bonds: span
                .bonds()
                .iter()
                .enumerate()
                .map(|(index, entity_span)| {
                    let [first, second] = span.graph().edge_endpoints(EdgeId(index as u32));
                    (
                        AtomId::from(first),
                        AtomId::from(second),
                        map_span(entity_span, |bond| {
                            BondDsl(bond.clone()).into_ir(&context.bond)
                        }),
                    )
                })
                .collect(),
            dative: span
                .dative_bonds()
                .ids()
                .map(|id| {
                    (
                        span.dative_bonds().donors(id).collect(),
                        span.dative_bonds().acceptor(id),
                        map_span(span.dative_bonds().attributes(id), |bond| {
                            DativeBondDsl(bond.clone()).into_ir(&context.dative_bond)
                        }),
                    )
                })
                .collect(),
            aromatic: span
                .aromatic_systems()
                .ids()
                .map(|id| {
                    (
                        span.aromatic_systems().atoms(id).collect(),
                        map_span(span.aromatic_systems().attributes(id), |system| {
                            AromaticSystemDsl(system.clone()).into_ir(&context.aromatic_system)
                        }),
                    )
                })
                .collect(),
            multicenter: span
                .multicenter_bonds()
                .ids()
                .map(|id| {
                    (
                        span.multicenter_bonds().atoms(id).collect(),
                        map_span(span.multicenter_bonds().attributes(id), |bond| {
                            MulticenterBondDsl(bond.clone()).into_ir(&context.multicenter_bond)
                        }),
                    )
                })
                .collect(),
            noncovalent: span
                .noncovalent_bonds()
                .ids()
                .map(|id| {
                    let [first, second] = span.noncovalent_bonds().atoms(id);
                    (
                        first,
                        second,
                        map_span(span.noncovalent_bonds().attributes(id), |bond| {
                            NoncovalentBondDsl(bond.clone()).into_ir(&context.noncovalent_bond)
                        }),
                    )
                })
                .collect(),
            stereo_atoms: span
                .stereo_atoms()
                .ids()
                .map(|id| {
                    (
                        span.stereo_atoms().site(id),
                        span.stereo_atoms().ligands(id).to_vec(),
                        map_span(span.stereo_atoms().attributes(id), |stereo| {
                            StereoAtomDsl(stereo.clone()).into_ir(&context.stereo_atom)
                        }),
                    )
                })
                .collect(),
            stereo_bonds: span
                .stereo_bonds()
                .ids()
                .map(|id| {
                    (
                        span.stereo_bonds().site(id),
                        span.stereo_bonds().ligands(id).to_vec(),
                        map_span(span.stereo_bonds().attributes(id), |stereo| {
                            StereoBondDsl(stereo.clone()).into_ir(&context.stereo_bond)
                        }),
                    )
                })
                .collect(),
            constraints: span.constraints().to_vec(),
        })
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
/// table. Aliases, bond endpoints, and constraint refs are resolved in `into_ir`.
#[derive(Debug, PartialEq)]
#[allow(clippy::type_complexity)]
pub(crate) struct SpanInput {
    atoms: Vec<(Option<String>, EntitySpan<AtomSpecInput>)>,
    bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<BondForm>)>,
    dative_bonds: Vec<(
        Option<String>,
        Vec<AtomRef>,
        AtomRef,
        EntitySpan<DativeBondForm>,
    )>,
    aromatic_systems: Vec<(Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemForm>)>,
    multicenter_bonds: Vec<(
        Option<String>,
        Vec<AtomRef>,
        EntitySpan<MulticenterBondForm>,
    )>,
    noncovalent_bonds: Vec<(
        Option<String>,
        [AtomRef; 2],
        EntitySpan<NoncovalentBondForm>,
    )>,
    stereo_atoms: Vec<(
        Option<String>,
        AtomRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoAtomForm>,
    )>,
    stereo_bonds: Vec<(
        Option<String>,
        BondRef,
        Vec<StereoLigandRef>,
        EntitySpan<StereoBondForm>,
    )>,
    constraints: Vec<ConstraintSpanInput>,
    atom_aliases: Vec<(String, Box<AtomDsl>)>,
}

const SPAN_VERBS: [&str; 3] = ["add", "modify", "remove"];

/// Split the optional outer `[<keyword> <body>]` wrapper off a span entry. The keyword is borrowed.
fn split_span_entry<'a, 'de>(edn: &'a Edn<'de>) -> (Option<&'a str>, &'a Edn<'de>) {
    if let Edn::Vector(v) = edn {
        if v.len() == 2 {
            if let Edn::Keyword(keyword) = &v[0] {
                return (Some(keyword.name()), &v[1]);
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
    let (keyword, body) = split_span_entry(edn);
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
    Ok((keyword.map(String::from), span))
}

/// Parse a complete bond-entry payload (`[a b bond]` or the `{:id :atoms :attrs}` map) and wrap its
/// `BondForm` into the given span side.
#[allow(clippy::type_complexity)]
fn bond_entry_span(
    payload: &Edn<'_>,
    wrap: impl Fn(BondForm) -> EntitySpan<BondForm>,
) -> Result<(Option<String>, [AtomRef; 2], EntitySpan<BondForm>), DeError> {
    let entry = parse_bond_entry(payload)?;
    Ok((
        entry.keyword,
        [entry.first, entry.second],
        wrap(entry.bond.0),
    ))
}

/// Split a bond `:modify` payload — `[a b X]` or `{:id :atoms [a b] :attrs X}` — into its keyword,
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
            optional_id_keyword(m)?,
            two_atom_refs(
                parse_vec(required_key(m, "atoms", "bond span")?, ":atoms", |e| {
                    AtomRef::from_edn(e)
                })?,
                "bond span",
            )?,
            required_key(m, "attrs", "bond span")?,
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
) -> Result<(Option<String>, [AtomRef; 2], EntitySpan<BondForm>), DeError> {
    match verb_wrapper(edn) {
        None => bond_entry_span(edn, EntitySpan::Unchanged),
        Some(("add", p)) => bond_entry_span(p, EntitySpan::Added),
        Some(("remove", p)) => bond_entry_span(p, EntitySpan::Removed),
        Some(("modify", p)) => {
            let (keyword, endpoints, value) = split_bond_frame(p)?;
            let (left, right) = pair(value, "bond span :modify")?;
            Ok((
                keyword,
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
// over the shared `parse_<entity>_entry`. For `:add`/`:remove` and the bare form the entry's `:attrs`
// is one dsl; for `:modify` it is a `[left right]` pair (participants are span-invariant).

#[allow(clippy::type_complexity)]
fn parse_dative_span_entry(
    edn: &Edn<'_>,
) -> Result<
    (
        Option<String>,
        Vec<AtomRef>,
        AtomRef,
        EntitySpan<DativeBondForm>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(DativeBondForm) -> EntitySpan<DativeBondForm>| {
        let e = parse_dative_bond_entry(p)?;
        Ok::<_, DeError>((e.keyword, e.donors, e.acceptor, wrap(e.bond.0)))
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
                required_key(m, "attrs", "dative span")?,
                "dative span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
) -> Result<(Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemForm>), DeError> {
    let full = |p: &Edn<'_>, wrap: fn(AromaticSystemForm) -> EntitySpan<AromaticSystemForm>| {
        let e = parse_aromatic_system_entry(p)?;
        Ok::<_, DeError>((e.keyword, e.atoms, wrap(e.system.0)))
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
                required_key(m, "attrs", "aromatic span")?,
                "aromatic span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
) -> Result<
    (
        Option<String>,
        Vec<AtomRef>,
        EntitySpan<MulticenterBondForm>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(MulticenterBondForm) -> EntitySpan<MulticenterBondForm>| {
        let e = parse_multicenter_bond_entry(p)?;
        Ok::<_, DeError>((e.keyword, e.atoms, wrap(e.bond.0)))
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
                required_key(m, "attrs", "multicenter span")?,
                "multicenter span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
) -> Result<
    (
        Option<String>,
        [AtomRef; 2],
        EntitySpan<NoncovalentBondForm>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(NoncovalentBondForm) -> EntitySpan<NoncovalentBondForm>| {
        let e = parse_noncovalent_bond_entry(p)?;
        Ok::<_, DeError>((e.keyword, [e.first, e.second], wrap(e.bond.0)))
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
                required_key(m, "attrs", "noncovalent span")?,
                "noncovalent span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
        EntitySpan<StereoAtomForm>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(StereoAtomForm) -> EntitySpan<StereoAtomForm>| {
        let e = parse_stereo_atom_entry(p)?;
        Ok::<_, DeError>((e.keyword, e.site, e.ligands, wrap(e.stereo.0)))
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
                required_key(m, "attrs", "stereo-atom span")?,
                "stereo-atom span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
        EntitySpan<StereoBondForm>,
    ),
    DeError,
> {
    let full = |p: &Edn<'_>, wrap: fn(StereoBondForm) -> EntitySpan<StereoBondForm>| {
        let e = parse_stereo_bond_entry(p)?;
        Ok::<_, DeError>((e.keyword, e.site, e.ligands, wrap(e.stereo.0)))
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
                required_key(m, "attrs", "stereo-bond span")?,
                "stereo-bond span :modify",
            )?;
            Ok((
                optional_id_keyword(m)?,
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
    context: &MoleculeContext,
) -> Result<EntitySpan<AtomForm>, ParseError> {
    Ok(match span {
        EntitySpan::Unchanged(s) => EntitySpan::Unchanged(resolve_atom_spec(s, context)?),
        EntitySpan::Added(s) => EntitySpan::Added(resolve_atom_spec(s, context)?),
        EntitySpan::Removed(s) => EntitySpan::Removed(resolve_atom_spec(s, context)?),
        EntitySpan::Modified {
            lhs: left,
            rhs: right,
        } => EntitySpan::Modified {
            lhs: resolve_atom_spec(left, context)?,
            rhs: resolve_atom_spec(right, context)?,
        },
    })
}

fn resolve_constraint_span(
    input: ConstraintSpanInput,
    context: &MoleculeContext,
) -> Result<ConstraintSpan, ParseError> {
    Ok(match input {
        ConstraintSpanInput::Unchanged(dsl) => ConstraintSpan::Unchanged(dsl.into_ir(context)?),
        ConstraintSpanInput::Added(dsl) => ConstraintSpan::Added(dsl.into_ir(context)?),
        ConstraintSpanInput::Removed(dsl) => ConstraintSpan::Removed(dsl.into_ir(context)?),
    })
}

impl SpanInput {
    /// Resolve the union-frame span: positions are the union ids (no fresh-id allocation), inline
    /// entity keywords and `:atom-aliases` populate the context, atom `AtomSpecInput` sides resolve to
    /// `AtomForm`, and participant, site, ligand, and constraint refs resolve against the context.
    /// Checked span construction requires every selected side reference to remain available after
    /// projection; chemistry and other semantic properties are not validated here.
    pub(crate) fn into_ir(self) -> Result<(ReactionSpan, MoleculeMetadata), ParseError> {
        let atom_count = self.atoms.len();
        let bond_count = self.bonds.len();

        // The span context: atoms take the union positions as their ids, then the bijective aliases.
        // Every ref resolves against it; `register_*` enforces keyword disjointness, and the
        // roundtrip `MoleculeMetadata` is moved out at the end.
        let mut context = MoleculeContext::default();
        for (keyword, _) in self.atoms.iter() {
            context.register_atom(keyword.clone())?;
        }
        for (name, dsl) in self.atom_aliases {
            context.register_atom_alias(name, *dsl)?;
        }

        // Resolve atoms (alias → AtomForm), bonds (endpoints + value), constraints.
        let mut atoms: Vec<EntitySpan<AtomForm>> = Vec::with_capacity(atom_count);
        for (_, span) in self.atoms {
            atoms.push(resolve_atom_span(span, &context)?);
        }
        let mut bonds: Vec<(AtomId, AtomId, EntitySpan<BondForm>)> = Vec::with_capacity(bond_count);
        for (keyword, [ref_a, ref_b], span) in self.bonds {
            let a = ref_a.resolve(&context)?;
            let b = ref_b.resolve(&context)?;
            context.register_bond(keyword, a, b)?;
            bonds.push((a, b, span));
        }

        let mut dative = Vec::with_capacity(self.dative_bonds.len());
        for (keyword, donors, acceptor, span) in self.dative_bonds {
            let acceptor_id = acceptor.resolve(&context)?;
            let donor_ids: Vec<AtomId> = donors
                .into_iter()
                .map(|d| d.resolve(&context))
                .collect::<Result<_, _>>()?;
            context.register_dative_bond(keyword, &donor_ids, acceptor_id)?;
            dative.push((donor_ids, acceptor_id, span));
        }

        let mut aromatic = Vec::with_capacity(self.aromatic_systems.len());
        for (keyword, atoms_ref, span) in self.aromatic_systems {
            let atom_ids: Vec<AtomId> = atoms_ref
                .into_iter()
                .map(|r| r.resolve(&context))
                .collect::<Result<_, _>>()?;
            context.register_aromatic_system(keyword, &atom_ids)?;
            aromatic.push((atom_ids, span));
        }

        let mut multicenter = Vec::with_capacity(self.multicenter_bonds.len());
        for (keyword, atoms_ref, span) in self.multicenter_bonds {
            let atom_ids: Vec<AtomId> = atoms_ref
                .into_iter()
                .map(|r| r.resolve(&context))
                .collect::<Result<_, _>>()?;
            context.register_multicenter_bond(keyword, &atom_ids)?;
            multicenter.push((atom_ids, span));
        }

        let mut noncovalent = Vec::with_capacity(self.noncovalent_bonds.len());
        for (keyword, [first, second], span) in self.noncovalent_bonds {
            let a = first.resolve(&context)?;
            let b = second.resolve(&context)?;
            context.register_noncovalent_bond(keyword, a, b)?;
            noncovalent.push((a, b, span));
        }

        let mut stereo_atoms = Vec::with_capacity(self.stereo_atoms.len());
        for (keyword, site, ligands, span) in self.stereo_atoms {
            let site_id = site.resolve(&context)?;
            let mut ligand_frame = Vec::with_capacity(ligands.len());
            for l in ligands {
                let a = l.atom.resolve(&context)?;
                ligand_frame.push(StereoLigand::new(a, l.kind));
            }
            context.register_stereo_atom(keyword, site_id, &ligand_frame)?;
            stereo_atoms.push((site_id, ligand_frame, span));
        }

        let mut stereo_bonds = Vec::with_capacity(self.stereo_bonds.len());
        for (keyword, site, ligands, span) in self.stereo_bonds {
            let site_id = site.resolve(&context)?;
            let mut ligand_frame = Vec::with_capacity(ligands.len());
            for l in ligands {
                let a = l.atom.resolve(&context)?;
                ligand_frame.push(StereoLigand::new(a, l.kind));
            }
            context.register_stereo_bond(keyword, site_id, &ligand_frame)?;
            stereo_bonds.push((site_id, ligand_frame, span));
        }

        let mut constraints: Vec<ConstraintSpan> = Vec::with_capacity(self.constraints.len());
        for input in self.constraints {
            constraints.push(resolve_constraint_span(input, &context)?);
        }

        let metadata = context.into_metadata();
        let span = ReactionSpan::try_from_entries(ReactionSpanEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        })
        .map_err(|error| ParseError::InvalidValue(error.to_string()))?;
        Ok((span, metadata))
    }
}

impl<'de> FromEdn<'de> for ReactionSpanDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let (span, metadata) = parse_span_input(edn)?
            .into_ir()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionSpanDsl::from_parts(span, metadata))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let span_input = read_span_input(&mut de)?;
        de.expect_eof()?;
        let (span, metadata) = span_input
            .into_ir()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionSpanDsl::from_parts(span, metadata))
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
    span: &EntitySpan<AtomForm>,
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
    match meta.keyword(Entity::Atom(id)) {
        Some(name) => {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(name.to_string())), body].into())
        }
        None => body,
    }
}

fn render_bond_span_entry(
    id: BondId,
    endpoints: [AtomId; 2],
    span: &EntitySpan<BondForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |bond: &BondForm| BondDsl::from_ref(bond).to_edn();
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
    span: &EntitySpan<DativeBondForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |b: &DativeBondForm| DativeBondDsl::from_ref(b).to_edn();
    let entry = |attributes_edn: Edn<'static>| {
        render_dative_entry(id, donors.iter().copied(), acceptor, attributes_edn, meta)
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
    span: &EntitySpan<AromaticSystemForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value =
        |s: &AromaticSystemForm| Edn::Str(Cow::Owned(AromaticSystemDsl::from_ref(s).to_string()));
    let entry = |attributes_edn: Edn<'static>| {
        render_aromatic_entry(id, atoms.iter().copied(), attributes_edn, meta)
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

fn render_multicenter_span_entry(
    id: MulticenterBondId,
    atoms: &[AtomId],
    span: &EntitySpan<MulticenterBondForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value =
        |b: &MulticenterBondForm| Edn::Str(Cow::Owned(MulticenterBondDsl::from_ref(b).to_string()));
    let entry = |attributes_edn: Edn<'static>| {
        render_multicenter_entry(id, atoms.iter().copied(), attributes_edn, meta)
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
    span: &EntitySpan<NoncovalentBondForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |b: &NoncovalentBondForm| NoncovalentBondDsl::from_ref(b).to_edn();
    let entry =
        |attributes_edn: Edn<'static>| render_noncovalent_entry(id, atoms, attributes_edn, meta);
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
    span: &EntitySpan<StereoAtomForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |s: &StereoAtomForm| StereoAtomDsl::from_ref(s).to_edn();
    let ligand_edns: Vec<Edn<'static>> = ligands
        .iter()
        .map(|&l| render_stereo_ligand(l, meta))
        .collect();
    let entry = |attributes_edn: Edn<'static>| {
        render_stereo_atom_entry(id, site, ligand_edns.clone(), attributes_edn, meta)
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
    span: &EntitySpan<StereoBondForm>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let value = |s: &StereoBondForm| StereoBondDsl::from_ref(s).to_edn();
    let ligand_edns: Vec<Edn<'static>> = ligands
        .iter()
        .map(|&l| render_stereo_ligand(l, meta))
        .collect();
    let entry = |attributes_edn: Edn<'static>| {
        render_stereo_bond_entry(id, site, ligand_edns.clone(), attributes_edn, meta)
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
        ConstraintDsl::from_ir(c, meta)
            .expect("ConstraintDsl::from_ir is infallible for well-formed graph IR")
            .to_edn()
    };
    match span {
        ConstraintSpan::Unchanged(c) => render(c),
        ConstraintSpan::Added(c) => single_key_map("add", render(c)),
        ConstraintSpan::Removed(c) => single_key_map("remove", render(c)),
    }
}

/// Render a `ReactionSpan` (+ its `MoleculeMetadata`) to the span EDN map.
fn render_span_edn(span: &ReactionSpan, meta: &MoleculeMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(4);
    let atoms: Vec<Edn<'static>> = span
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, entity_span)| render_atom_span_entry(AtomId(i as u32), entity_span, meta))
        .collect();
    map.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));

    let bonds: Vec<Edn<'static>> = span
        .bonds()
        .iter()
        .enumerate()
        .map(|(j, entity_span)| {
            let [a, b] = span.graph().edge_endpoints(EdgeId(j as u32));
            render_bond_span_entry(
                BondId(j as u32),
                [AtomId::from(a), AtomId::from(b)],
                entity_span,
                meta,
            )
        })
        .collect();
    if !bonds.is_empty() {
        map.insert(Edn::keyword("bonds"), Edn::Vector(bonds.into()));
    }

    let dative = span.dative_bonds();
    let dative_bonds: Vec<Edn<'static>> = dative
        .ids()
        .map(|rid| {
            let acceptor = dative.acceptor(rid);
            let donors: Vec<AtomId> = dative.donors(rid).collect();
            render_dative_span_entry(
                DativeBondId(rid.index() as u32),
                &donors,
                acceptor,
                dative.attributes(rid),
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

    let aromatic = span.aromatic_systems();
    let aromatic_systems: Vec<Edn<'static>> = aromatic
        .ids()
        .map(|rid| {
            let atoms: Vec<AtomId> = aromatic.atoms(rid).collect();
            render_aromatic_span_entry(
                AromaticSystemId(rid.index() as u32),
                &atoms,
                aromatic.attributes(rid),
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

    let multicenter = span.multicenter_bonds();
    let multicenter_bonds: Vec<Edn<'static>> = multicenter
        .ids()
        .map(|rid| {
            let atoms: Vec<AtomId> = multicenter.atoms(rid).collect();
            render_multicenter_span_entry(
                MulticenterBondId(rid.index() as u32),
                &atoms,
                multicenter.attributes(rid),
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

    let noncovalent = span.noncovalent_bonds();
    let noncovalent_bonds: Vec<Edn<'static>> = noncovalent
        .ids()
        .map(|rid| {
            let [a, b] = noncovalent.atoms(rid);
            render_noncovalent_span_entry(
                NoncovalentBondId(rid.index() as u32),
                [a, b],
                noncovalent.attributes(rid),
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

    let stereo_a = span.stereo_atoms();
    let stereo_atoms: Vec<Edn<'static>> = stereo_a
        .ids()
        .map(|rid| {
            let site = stereo_a.site(rid);
            render_stereo_atom_span_entry(
                StereoAtomId(rid.index() as u32),
                site,
                stereo_a.ligands(rid),
                stereo_a.attributes(rid),
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

    let stereo_b = span.stereo_bonds();
    let stereo_bonds: Vec<Edn<'static>> = stereo_b
        .ids()
        .map(|rid| {
            let site = stereo_b.site(rid);
            render_stereo_bond_span_entry(
                StereoBondId(rid.index() as u32),
                site,
                stereo_b.ligands(rid),
                stereo_b.attributes(rid),
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

    let constraints: Vec<Edn<'static>> = span
        .constraints()
        .iter()
        .map(|span| render_constraint_span_entry(span, meta))
        .collect();
    if !constraints.is_empty() {
        map.insert(Edn::keyword("constraints"), Edn::Vector(constraints.into()));
    }

    let aliases = meta.iter_atom_aliases();
    if aliases.len() != 0 {
        let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(aliases.len() * 2);
        for (name, dsl) in aliases {
            pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
            pairs.push(dsl.to_edn());
        }
        map.insert(Edn::keyword("atom-aliases"), Edn::Vector(pairs.into()));
    }

    Edn::Map(map)
}

impl ToEdn for ReactionSpanDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_span_edn(&self.span, &self.metadata)
    }
}

impl Display for ReactionSpanDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for ReactionSpan {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        ReactionSpanDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        ReactionSpanDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for ReactionSpan {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for ReactionSpan {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `ReactionSpan`: positional form (no entity keywords or aliases) since
/// the reaction span carries no metadata. For keyword/alias-bearing output, wrap in [`ReactionSpanDsl`].
impl ToEdn for ReactionSpan {
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
    use umol_graph_core::Correspondence;

    use super::*;
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::{BondConstraintForm, Constraint, Constraints, MoleculeConstraint};
    use crate::ir::delta::{AtomDelta, BondDelta, ConstraintDelta, Delta, Deltas};
    use crate::ir::edit::BondFieldChange;
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
    use crate::ir::num::NumForm;
    use crate::ir::reaction::Reaction;
    use crate::ir::{MoleculeCorrespondence, ReactionSpanIntegrityError};

    #[fixture]
    fn populated_reaction_span_dsl() -> ReactionSpanDsl {
        r#"{
            :atoms [[:a "C"] "F" "Cl" "Br" "I"]
            :bonds [
                {:id :b :atoms [0 1] :attrs "2"}
                [0 2 "1"]
                [0 3 "1"]
                [0 4 "1"]
            ]
            :dative-bonds [{:id :d :donors [1] :acceptor 0 :attrs "1#R"}]
            :aromatic-systems [{:id :ar :atoms [0 1] :attrs "*#e2"}]
            :multicenter-bonds [{:id :m :atoms [0 1] :attrs "*#e2"}]
            :noncovalent-bonds [{:id :n :atoms [0 1] :attrs "Hbd"}]
            :stereo-atoms [{:id :sa :site 0 :ligands [1 2 3 4] :attrs "Th1"}]
            :stereo-bonds [{:id :sb :site :b :ligands [2 3 [:h 1] [:lp 1]] :attrs "Ct1"}]
            :atom-aliases [:x "O"]
        }"#
        .parse()
        .unwrap()
    }

    #[rstest]
    #[case::empty(None)]
    #[case::atom(Some(Entity::Atom(AtomId(4))))]
    #[case::bond(Some(Entity::Bond(BondId(3))))]
    #[case::dative_bond(Some(Entity::DativeBond(DativeBondId(0))))]
    #[case::aromatic_system(Some(Entity::AromaticSystem(AromaticSystemId(0))))]
    #[case::multicenter_bond(Some(Entity::MulticenterBond(MulticenterBondId(0))))]
    #[case::noncovalent_bond(Some(Entity::NoncovalentBond(NoncovalentBondId(0))))]
    #[case::stereo_atom(Some(Entity::StereoAtom(StereoAtomId(0))))]
    #[case::stereo_bond(Some(Entity::StereoBond(StereoBondId(0))))]
    fn test_reaction_span_dsl_new(
        populated_reaction_span_dsl: ReactionSpanDsl,
        #[case] entity: Option<Entity>,
    ) {
        let span = populated_reaction_span_dsl.into_parts().0;
        let mut metadata = MoleculeMetadata::new();
        if let Some(entity) = entity {
            metadata.set_keyword(entity, "key").unwrap();
        }

        let actual = ReactionSpanDsl::new(span.clone(), metadata.clone()).unwrap();

        assert_eq!(actual.into_parts(), (span, metadata));
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(5)))]
    #[case::bond(Entity::Bond(BondId(4)))]
    #[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
    #[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
    #[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
    #[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
    fn test_reaction_span_dsl_new_error(
        populated_reaction_span_dsl: ReactionSpanDsl,
        #[case] entity: Entity,
    ) {
        let span = populated_reaction_span_dsl.into_parts().0;
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(entity, "key").unwrap();

        assert_eq!(
            ReactionSpanDsl::new(span, metadata),
            Err(MetadataError::EntityOutOfRange(entity))
        );
    }

    #[rstest]
    fn test_reaction_span_dsl_new_parsed(populated_reaction_span_dsl: ReactionSpanDsl) {
        let expected = populated_reaction_span_dsl.clone();
        let (span, metadata) = populated_reaction_span_dsl.into_parts();

        assert_eq!(ReactionSpanDsl::new(span, metadata).unwrap(), expected);
    }

    // Modified bond + Unchanged atoms + Unchanged molecule-constraint.
    #[rstest]
    #[case::modify(Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
        })]),
    ))]
    // Unchanged / Removed / Added atoms and bonds + an Added constraint.
    #[case::add_remove(Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(1),
                attributes: AtomForm::from_element(Element::O),
            }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: BondForm::from_order(1),
            }),
            Delta::Atom(AtomDelta::Add {
                id: AtomId(2),
                attributes: AtomForm::from_element(Element::N),
            }),
            Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(2)],
                attributes: BondForm::from_order(1),
            }),
            Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
            )),
        ]),
    ))]
    fn test_reaction_span_dsl_from_ir(#[case] reaction: Reaction) {
        let span = reaction.to_reaction_span().unwrap();
        let cfg = MoleculeDefaults::default();
        assert_eq!(ReactionSpanDsl::from_ir(&span, &cfg).into_ir(&cfg), span);
    }

    #[rstest]
    fn test_reaction_span_dsl_from_ir_superimposed() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let rhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(AtomId(0), AtomId(0))], 2, 2).unwrap();
        let correspondence = MoleculeCorrespondence::induce(&lhs, &rhs, atoms).unwrap();
        let expected = ReactionSpan::superimpose(&lhs, &rhs, &correspondence).unwrap();
        let defaults = MoleculeDefaults::default();

        assert_eq!(
            ReactionSpanDsl::from_ir(&expected, &defaults).into_ir(&defaults),
            expected,
        );
    }

    #[rstest]
    #[case::unchanged(r#""C""#, (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::C)))))))]
    #[case::add(r#"{:add "O"}"#, (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::O)))))))]
    #[case::remove(r#"{:remove "O"}"#, (None, EntitySpan::Removed(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::O)))))))]
    #[case::modify(r#"{:modify ["C" "N"]}"#, (None, EntitySpan::Modified {
        lhs:AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::C)))),
        rhs:AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::N)))),
    }))]
    #[case::with_keyword(r#"[:c "C"]"#, (Some("c".to_string()), EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::C)))))))]
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
    #[case::unchanged("[0 1 :single]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondForm::from_order(1))))]
    #[case::add("{:add [0 2 :single]}", (None, [AtomRef::Index(0), AtomRef::Index(2)], EntitySpan::Added(BondForm::from_order(1))))]
    #[case::remove("{:remove [0 1 :single]}", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Removed(BondForm::from_order(1))))]
    #[case::modify("{:modify [0 2 [:single :double]]}", (None, [AtomRef::Index(0), AtomRef::Index(2)], EntitySpan::Modified {
        lhs:BondForm::from_order(1),
        rhs:BondForm::from_order(2),
    }))]
    #[case::unchanged_map_id("{:id :b1 :atoms [0 1] :attrs :single}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondForm::from_order(1))))]
    #[case::modify_map_id("{:modify {:id :b1 :atoms [0 1] :attrs [:single :double]}}", (Some("b1".to_string()), [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Modified {
        lhs:BondForm::from_order(1),
        rhs:BondForm::from_order(2),
    }))]
    #[case::triple("[0 1 :triple]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(BondForm::from_order(3))))]
    #[case::aromatic("[0 1 :aromatic]", (None, [AtomRef::Index(0), AtomRef::Index(1)], EntitySpan::Unchanged(
        BondForm::from_order(1).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    )))]
    fn test_parse_bond_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, [AtomRef; 2], EntitySpan<BondForm>),
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
    #[case::unchanged(r#"{:donors [0] :acceptor 1 :attrs "1#R"}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::add(r#"{:add {:donors [0] :acceptor 1 :attrs "1#R"}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Added(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:donors [0] :acceptor 1 :attrs "1#R"}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Removed(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:donors [0] :acceptor 1 :attrs ["1#R" "2#R"]}}"#, (
        None, vec![AtomRef::Index(0)], AtomRef::Index(1),
        EntitySpan::Modified {
            lhs:DativeBondDsl::from_str("1#R").unwrap().0,
            rhs:DativeBondDsl::from_str("2#R").unwrap().0,
        },
    ))]
    #[case::with_id(r#"{:id :d1 :donors [0 2] :acceptor 1 :attrs "1#R"}"#, (
        Some("d1".to_string()), vec![AtomRef::Index(0), AtomRef::Index(2)], AtomRef::Index(1),
        EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
    ))]
    fn test_parse_dative_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            Vec<AtomRef>,
            AtomRef,
            EntitySpan<DativeBondForm>,
        ),
    ) {
        assert_eq!(
            parse_dative_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1 2] :attrs "*#e6"}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1), AtomRef::Index(2)],
        EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0),
    ))]
    #[case::add(r#"{:add {:atoms [0 1] :attrs "*#e2"}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Added(AromaticSystemDsl::from_str("*#e2").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:atoms [0 1] :attrs ["*#e2" "*#e4"]}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Modified {
            lhs:AromaticSystemDsl::from_str("*#e2").unwrap().0,
            rhs:AromaticSystemDsl::from_str("*#e4").unwrap().0,
        },
    ))]
    fn test_parse_aromatic_span_entry(
        #[case] input: &str,
        #[case] expected: (Option<String>, Vec<AtomRef>, EntitySpan<AromaticSystemForm>),
    ) {
        assert_eq!(
            parse_aromatic_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1] :attrs "*#e2"}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Unchanged(MulticenterBondDsl::from_str("*#e2").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:atoms [0 1 2] :attrs "*#e3"}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1), AtomRef::Index(2)],
        EntitySpan::Removed(MulticenterBondDsl::from_str("*#e3").unwrap().0),
    ))]
    #[case::modify(r#"{:modify {:atoms [0 1] :attrs ["*#e2" "*#e4"]}}"#, (
        None, vec![AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Modified {
            lhs:MulticenterBondDsl::from_str("*#e2").unwrap().0,
            rhs:MulticenterBondDsl::from_str("*#e4").unwrap().0,
        },
    ))]
    fn test_parse_multicenter_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            Vec<AtomRef>,
            EntitySpan<MulticenterBondForm>,
        ),
    ) {
        assert_eq!(
            parse_multicenter_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:atoms [0 1] :attrs "Hbd"}"#, (
        None, [AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:atoms [0 1] :attrs "Hbd"}}"#, (
        None, [AtomRef::Index(0), AtomRef::Index(1)],
        EntitySpan::Removed(NoncovalentBondDsl::from_str("Hbd").unwrap().0),
    ))]
    fn test_parse_noncovalent_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            [AtomRef; 2],
            EntitySpan<NoncovalentBondForm>,
        ),
    ) {
        assert_eq!(
            parse_noncovalent_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:site 0 :ligands [1 2 3 4] :attrs "Th1"}"#, (
        None, AtomRef::Index(0),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(1) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(2) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(4) },
        ],
        EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0),
    ))]
    #[case::add(r#"{:add {:site 0 :ligands [1 2 [:h 3]] :attrs "Th1"}}"#, (
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
            EntitySpan<StereoAtomForm>,
        ),
    ) {
        assert_eq!(
            parse_stereo_atom_span_entry(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::unchanged(r#"{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}"#, (
        None, BondRef::Index(1),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(0) },
            StereoLigandRef { kind: StereoLigandKind::ImplicitHydrogen, atom: AtomRef::Index(1) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
            StereoLigandRef { kind: StereoLigandKind::ImplicitHydrogen, atom: AtomRef::Index(2) },
        ],
        EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0),
    ))]
    #[case::remove(r#"{:remove {:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}}"#, (
        None, BondRef::Index(1),
        vec![
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(0) },
            StereoLigandRef { kind: StereoLigandKind::ImplicitHydrogen, atom: AtomRef::Index(1) },
            StereoLigandRef { kind: StereoLigandKind::Atom, atom: AtomRef::Index(3) },
            StereoLigandRef { kind: StereoLigandKind::ImplicitHydrogen, atom: AtomRef::Index(2) },
        ],
        EntitySpan::Removed(StereoBondDsl::from_str("Ct1").unwrap().0),
    ))]
    fn test_parse_stereo_bond_span_entry(
        #[case] input: &str,
        #[case] expected: (
            Option<String>,
            BondRef,
            Vec<StereoLigandRef>,
            EntitySpan<StereoBondForm>,
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
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::C)))))),
                (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::O)))))),
            ],
            bonds: vec![(
                None,
                [AtomRef::Index(0), AtomRef::Index(1)],
                EntitySpan::Unchanged(BondForm::from_order(1)),
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
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::C)))))),
                (None, EntitySpan::Unchanged(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::O)))))),
            ],
            bonds: vec![(
                None,
                [AtomRef::Index(0), AtomRef::Index(1)],
                EntitySpan::Unchanged(BondForm::from_order(1)),
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
                (None, EntitySpan::Added(AtomSpecInput::Bare(Box::new(AtomDsl(AtomForm::from_element(Element::O)))))),
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
        ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Added(AtomForm::from_element(Element::O)),
            ],
            bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1)))],
            ..Default::default()
        }),
        {
            let mut metadata = MoleculeMetadata::new();
            metadata
                .add_atom_alias("nu", AtomDsl(AtomForm::from_element(Element::C)))
                .unwrap();
            metadata
        },
    )]
    #[case::equivalent_modified(
        r#"{:atoms [{:modify ["C#c1" "C#c{1}"]}]}"#,
        ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: vec![EntitySpan::Modified {
                lhs: AtomForm::from_element(Element::C).with_charge(1_i64),
                rhs: AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])),
            }],
            ..Default::default()
        }),
        MoleculeMetadata::new(),
    )]
    #[case::dative_overlay(
        r#"{:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :attrs "1#R"}]}"#,
        ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::N)),
            ],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0),
            )],
            ..Default::default()
        }),
        {
            let mut metadata = MoleculeMetadata::new();
            metadata
                .set_keyword(Entity::DativeBond(DativeBondId(0)), "d1")
                .unwrap();
            metadata
        },
    )]
    fn test_span_input_into_ir(
        #[case] input: &str,
        #[case] expected_span: ReactionSpan,
        #[case] expected_metadata: MoleculeMetadata,
    ) {
        assert_eq!(
            parse_span_input(&read_string(input).unwrap())
                .unwrap()
                .into_ir()
                .unwrap(),
            (expected_span, expected_metadata),
        );
    }

    #[rstest]
    #[case::atom(r#"{:atoms [[:a "C"]]}"#, Entity::Atom(AtomId(0)), "a")]
    #[case::bond(
        r#"{:atoms ["C" "C"] :bonds [{:id :b :atoms [0 1] :attrs :single}]}"#,
        Entity::Bond(BondId(0)),
        "b"
    )]
    #[case::dative_bond(
        r#"{:atoms ["C" "N"] :dative-bonds [{:id :d :donors [0] :acceptor 1 :attrs "1#R"}]}"#,
        Entity::DativeBond(DativeBondId(0)),
        "d"
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :aromatic-systems [{:id :a :atoms [0 1] :attrs "*#e2"}]}"#,
        Entity::AromaticSystem(AromaticSystemId(0)),
        "a"
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["C" "C"] :multicenter-bonds [{:id :m :atoms [0 1] :attrs "*#e2"}]}"#,
        Entity::MulticenterBond(MulticenterBondId(0)),
        "m"
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["N" "H"] :noncovalent-bonds [{:id :n :atoms [0 1] :attrs "Hbd"}]}"#,
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        "n"
    )]
    #[case::stereo_atom(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:id :s :site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#,
        Entity::StereoAtom(StereoAtomId(0)),
        "s",
    )]
    #[case::stereo_bond(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :s :site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#,
        Entity::StereoBond(StereoBondId(0)),
        "s",
    )]
    fn test_span_input_into_ir_metadata(
        #[case] input: &str,
        #[case] entity: Entity,
        #[case] keyword: &str,
    ) {
        let (_, metadata) = parse_span_input(&read_string(input).unwrap())
            .unwrap()
            .into_ir()
            .unwrap();

        assert_eq!(
            metadata.iter_keywords().collect::<Vec<_>>(),
            vec![(entity, keyword)]
        );
    }

    #[rstest]
    fn test_span_input_into_ir_structural_bond_ref() {
        // A structural bond ref ({:atoms [0 1]}) names the bond by its endpoints, resolved against
        // the context's participant lookup.
        let input = r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#;
        let (span, _) = parse_span_input(&read_string(input).unwrap())
            .unwrap()
            .into_ir()
            .unwrap();
        assert_eq!(
            span.constraints().to_vec(),
            vec![ConstraintSpan::Unchanged(Constraint::Bond(
                BondId(0),
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
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
    #[case::duplicate_keyword(
        r#"{:atoms [[:a "C"] [:a "O"]]}"#,
        ParseError::DuplicateKeyword("a".to_string()),
    )]
    #[case::side_local_bond_mismatch(
        r#"{:atoms ["C" {:add "O"}] :bonds [[0 1 :single]]}"#,
        ParseError::InvalidValue(
            "reaction span lhs is not a valid molecule representation: molecule references unavailable atom 1".to_string()
        ),
    )]
    #[case::side_local_dative_mismatch(
        r#"{:atoms ["C" {:add "N"}] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#R"}]}"#,
        ParseError::InvalidValue(
            "reaction span lhs is not a valid molecule representation: molecule references unavailable atom 1".to_string()
        ),
    )]
    #[case::side_local_stereo_bond_mismatch(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:add [1 2 "2"]} [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#,
        ParseError::InvalidValue(
            "reaction span lhs is not a valid molecule representation: molecule references unavailable bond 2".to_string()
        ),
    )]
    #[case::atom_repeated_virtual(
        r#"{:atoms ["C" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 [:h 0] [:h 0]] :attrs "Th0"}]}"#,
        ParseError::InvalidValue(ReactionSpanIntegrityError::Lhs(
            MoleculeIntegrityError::DuplicateStereoLigand {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                ligand: StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            },
        ).to_string()),
    )]
    #[case::atom_oversized(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "N" "O" "S"] :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4 5 6 7] :attrs "*"}}]}"#,
        ParseError::InvalidValue(ReactionSpanIntegrityError::Rhs(
            MoleculeIntegrityError::StereoFrameDegreeTooLarge {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                degree: 7,
                maximum: 6,
            },
        ).to_string()),
    )]
    #[case::bond_repeated_virtual(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 1]] :attrs "Ct0"}]}"#,
        ParseError::InvalidValue(ReactionSpanIntegrityError::Lhs(
            MoleculeIntegrityError::DuplicateStereoLigand {
                entity: Entity::StereoBond(StereoBondId(0)),
                ligand: StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            },
        ).to_string()),
    )]
    #[case::bond_oversized(
        r#"{:atoms ["C" "C" "F" "Cl" "Br" "I" "N" "O" "S"] :bonds [[0 1 "2"]] :stereo-bonds [{:add {:site 0 :ligands [2 3 4 5 6 7 8] :attrs "*"}}]}"#,
        ParseError::InvalidValue(ReactionSpanIntegrityError::Rhs(
            MoleculeIntegrityError::StereoFrameDegreeTooLarge {
                entity: Entity::StereoBond(StereoBondId(0)),
                degree: 7,
                maximum: 6,
            },
        ).to_string()),
    )]
    fn test_span_input_into_ir_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(
            parse_span_input(&read_string(input).unwrap())
                .unwrap()
                .into_ir()
                .unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::span(
        r#"{:atoms ["C" {:add "O"}] :bonds [{:add [0 1 :single]}]}"#,
        ReactionSpanDsl::new(
            ReactionSpan::from_entries(ReactionSpanEntries {
                atoms: vec![
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Added(AtomForm::from_element(Element::O)),
                ],
                bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1)))],
                ..Default::default()
            }),
            MoleculeMetadata::new(),
        ).unwrap(),
    )]
    #[case::plain_molecule(
        r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#,
        ReactionSpanDsl::new(
            ReactionSpan::from_entries(ReactionSpanEntries {
                atoms: vec![
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
                ],
                bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1)))],
                ..Default::default()
            }),
            MoleculeMetadata::new(),
        ).unwrap(),
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
        EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        MoleculeMetadata::new(),
        r#""C""#
    )]
    #[case::add(
        AtomId(0),
        EntitySpan::Added(AtomForm::from_element(Element::O)),
        MoleculeMetadata::new(),
        r#"{:add "O"}"#
    )]
    #[case::remove(
        AtomId(0),
        EntitySpan::Removed(AtomForm::from_element(Element::O)),
        MoleculeMetadata::new(),
        r#"{:remove "O"}"#
    )]
    #[case::modify(AtomId(0), EntitySpan::Modified { lhs:AtomForm::from_element(Element::C), rhs:AtomForm::from_element(Element::N) }, MoleculeMetadata::new(), r#"{:modify ["C" "N"]}"#)]
    #[case::with_keyword(AtomId(0), EntitySpan::Unchanged(AtomForm::from_element(Element::C)), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(0)), "c").unwrap();
        metadata
    }, r#"[:c "C"]"#)]
    #[case::alias(AtomId(0), EntitySpan::Unchanged(AtomForm::from_element(Element::C)), {
        let mut metadata = MoleculeMetadata::new();
        metadata.add_atom_alias("nu", AtomDsl(AtomForm::from_element(Element::C))).unwrap();
        metadata
    }, r#":nu"#)]
    fn test_render_atom_span_entry(
        #[case] id: AtomId,
        #[case] span: EntitySpan<AtomForm>,
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
        EntitySpan::Unchanged(BondForm::from_order(1)),
        MoleculeMetadata::new(),
        "[0 1 :single]"
    )]
    #[case::add(
        EntitySpan::Added(BondForm::from_order(2)),
        MoleculeMetadata::new(),
        "{:add [0 1 :double]}"
    )]
    #[case::remove(
        EntitySpan::Removed(BondForm::from_order(1)),
        MoleculeMetadata::new(),
        "{:remove [0 1 :single]}"
    )]
    #[case::modify(EntitySpan::Modified { lhs:BondForm::from_order(1), rhs:BondForm::from_order(2) }, MoleculeMetadata::new(), "{:modify [0 1 [:single :double]]}")]
    #[case::with_id(EntitySpan::Unchanged(BondForm::from_order(1)), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Bond(BondId(0)), "b1").unwrap();
        metadata
    }, "{:id :b1 :atoms [0 1] :attrs :single}")]
    #[case::aromatic(EntitySpan::Unchanged(BondForm::from_order(1).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true)))), MoleculeMetadata::new(), "[0 1 :aromatic]")]
    fn test_render_bond_span_entry(
        #[case] span: EntitySpan<BondForm>,
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
    #[case::unchanged(EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:donors [0] :acceptor 1 :attrs "1#R"}"#)]
    #[case::add(EntitySpan::Added(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:add {:donors [0] :acceptor 1 :attrs "1#R"}}"#)]
    #[case::remove(EntitySpan::Removed(DativeBondDsl::from_str("1#R").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:donors [0] :acceptor 1 :attrs "1#R"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:DativeBondDsl::from_str("1#R").unwrap().0, rhs:DativeBondDsl::from_str("2#R").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:donors [0] :acceptor 1 :attrs ["1#R" "2#R"]}}"#)]
    #[case::with_id(EntitySpan::Unchanged(DativeBondDsl::from_str("1#R").unwrap().0), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::DativeBond(DativeBondId(0)), "d1").unwrap();
        metadata
    }, r#"{:id :d1 :donors [0] :acceptor 1 :attrs "1#R"}"#)]
    fn test_render_dative_span_entry(
        #[case] span: EntitySpan<DativeBondForm>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_dative_span_entry(DativeBondId(0), &[AtomId(0)], AtomId(1), &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1 2] :attrs "*#e6"}"#)]
    #[case::add(EntitySpan::Added(AromaticSystemDsl::from_str("*#e6").unwrap().0), MoleculeMetadata::new(), r#"{:add {:atoms [0 1 2] :attrs "*#e6"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:AromaticSystemDsl::from_str("*#e6").unwrap().0, rhs:AromaticSystemDsl::from_str("*#e2").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:atoms [0 1 2] :attrs ["*#e6" "*#e2"]}}"#)]
    #[case::with_id(EntitySpan::Unchanged(AromaticSystemDsl::from_str("*#e6").unwrap().0), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::AromaticSystem(AromaticSystemId(0)), "ar1").unwrap();
        metadata
    }, r#"{:id :ar1 :atoms [0 1 2] :attrs "*#e6"}"#)]
    fn test_render_aromatic_span_entry(
        #[case] span: EntitySpan<AromaticSystemForm>,
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
    #[case::unchanged(EntitySpan::Unchanged(MulticenterBondDsl::from_str("*#e2").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1] :attrs "*#e2"}"#)]
    #[case::remove(EntitySpan::Removed(MulticenterBondDsl::from_str("*#e2").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:atoms [0 1] :attrs "*#e2"}}"#)]
    #[case::modify(EntitySpan::Modified { lhs:MulticenterBondDsl::from_str("*#e2").unwrap().0, rhs:MulticenterBondDsl::from_str("*#e4").unwrap().0 }, MoleculeMetadata::new(), r#"{:modify {:atoms [0 1] :attrs ["*#e2" "*#e4"]}}"#)]
    fn test_render_multicenter_span_entry(
        #[case] span: EntitySpan<MulticenterBondForm>,
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
    #[case::unchanged(EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0), MoleculeMetadata::new(), r#"{:atoms [0 1] :attrs "Hbd"}"#)]
    #[case::add(EntitySpan::Added(NoncovalentBondDsl::from_str("Hbd").unwrap().0), MoleculeMetadata::new(), r#"{:add {:atoms [0 1] :attrs "Hbd"}}"#)]
    #[case::with_id(EntitySpan::Unchanged(NoncovalentBondDsl::from_str("Hbd").unwrap().0), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::NoncovalentBond(NoncovalentBondId(0)), "nc1").unwrap();
        metadata
    }, r#"{:id :nc1 :atoms [0 1] :attrs "Hbd"}"#)]
    fn test_render_noncovalent_span_entry(
        #[case] span: EntitySpan<NoncovalentBondForm>,
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
    #[case::unchanged(EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0), MoleculeMetadata::new(), r#"{:site 0 :ligands [1 2 3 4] :attrs :cw}"#)]
    #[case::add(EntitySpan::Added(StereoAtomDsl::from_str("Th1").unwrap().0), MoleculeMetadata::new(), r#"{:add {:site 0 :ligands [1 2 3 4] :attrs :cw}}"#)]
    #[case::with_id(EntitySpan::Unchanged(StereoAtomDsl::from_str("Th1").unwrap().0), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::StereoAtom(StereoAtomId(0)), "s1").unwrap();
        metadata
    }, r#"{:id :s1 :site 0 :ligands [1 2 3 4] :attrs :cw}"#)]
    fn test_render_stereo_atom_span_entry(
        #[case] span: EntitySpan<StereoAtomForm>,
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
    #[case::unchanged(EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0), MoleculeMetadata::new(), r#"{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs :e}"#)]
    #[case::remove(EntitySpan::Removed(StereoBondDsl::from_str("Ct1").unwrap().0), MoleculeMetadata::new(), r#"{:remove {:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs :e}}"#)]
    #[case::with_id(EntitySpan::Unchanged(StereoBondDsl::from_str("Ct1").unwrap().0), {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::StereoBond(StereoBondId(0)), "sb1").unwrap();
        metadata
    }, r#"{:id :sb1 :site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs :e}"#)]
    fn test_render_stereo_bond_span_entry(
        #[case] span: EntitySpan<StereoBondForm>,
        #[case] meta: MoleculeMetadata,
        #[case] expected: &str,
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ];
        assert_eq!(
            render_stereo_bond_span_entry(StereoBondId(0), BondId(1), &ligands, &span, &meta),
            read_string(expected).unwrap(),
        );
    }

    // `ReactionSpanDsl` (span + metadata) round-trips through the EDN surface (render → reparse).
    #[rstest]
    #[case::plain_molecule(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::modify(r#"{:atoms ["C" "C"] :bonds [{:modify [0 1 [:single :double]]}]}"#)]
    #[case::add_remove(r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :single]}]}"#)]
    #[case::constraint(r#"{:atoms ["C"] :constraints [{:connected {}} {:add {:connected {}}}]}"#)]
    #[case::aliases(r#"{:atoms [:nu {:add "O"}] :atom-aliases [:nu "C"]}"#)]
    #[case::atom_keyword(r#"{:atoms [[:c "C"]]}"#)]
    #[case::bond_keyword(r#"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :attrs :single}]}"#)]
    #[case::dative(
        r#"{:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :attrs "1#R"}]}"#
    )]
    #[case::dative_add(
        r#"{:atoms ["C" "N"] :dative-bonds [{:add {:donors [0] :acceptor 1 :attrs "1#R"}}]}"#
    )]
    #[case::dative_modify(r#"{:atoms ["C" "N"] :dative-bonds [{:modify {:donors [0] :acceptor 1 :attrs ["1#R" "2#R"]}}]}"#)]
    #[case::aromatic_keyword(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:id :a1 :atoms [0 1 2 3 4 5] :attrs "*#e6"}]}"#)]
    #[case::multicenter_keyword(
        r#"{:atoms ["C" "C"] :multicenter-bonds [{:id :m1 :atoms [0 1] :attrs "*#e2"}]}"#
    )]
    #[case::noncovalent(
        r#"{:atoms ["N" "H"] :noncovalent-bonds [{:remove {:id :n1 :atoms [0 1] :attrs "Hbd"}}]}"#
    )]
    #[case::stereo_atom_keyword(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#)]
    #[case::stereo_bond_keyword(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :s1 :site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#)]
    fn test_reaction_span_dsl_to_edn(#[case] input: &str) {
        let dsl = ReactionSpanDsl::from_str(input).unwrap();
        assert_eq!(ReactionSpanDsl::from_edn(&dsl.to_edn()).unwrap(), dsl);
    }

    // `ReactionSpan` renders positionally and round-trips (metadata discarded on both sides).
    #[rstest]
    #[case::plain_molecule(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::modify(r#"{:atoms ["C" "C"] :bonds [{:modify [0 1 [:single :double]]}]}"#)]
    #[case::add_remove(r#"{:atoms ["C" {:remove "O"} {:add "N"}] :bonds [{:remove [0 1 :single]} {:add [0 2 :single]}]}"#)]
    #[case::dative(r#"{:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#R"}]}"#)]
    #[case::dative_modify(r#"{:atoms ["C" "N"] :dative-bonds [{:modify {:donors [0] :acceptor 1 :attrs ["1#R" "2#R"]}}]}"#)]
    #[case::aromatic(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#e6"}]}"#)]
    #[case::multicenter(
        r#"{:atoms ["C" "C"] :multicenter-bonds [{:add {:atoms [0 1] :attrs "*#e2"}}]}"#
    )]
    #[case::noncovalent(r#"{:atoms ["N" "H"] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#)]
    #[case::stereo_atom(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#)]
    #[case::stereo_bond(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#)]
    fn test_reaction_span_to_edn(#[case] input: &str) {
        let span = ReactionSpan::from_str(input).unwrap();
        assert_eq!(ReactionSpan::from_edn(&span.to_edn()).unwrap(), span);
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
        assert_eq!(count, 22);
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
    fn test_reaction_span_constraint_projection(
        #[case] input: &str,
        #[case] lhs: Vec<Constraint>,
        #[case] rhs: Vec<Constraint>,
    ) {
        let span = ReactionSpan::from_str(input).unwrap();
        assert_eq!(span.lhs().constraints().as_slice(), lhs.as_slice());
        assert_eq!(span.rhs().constraints().as_slice(), rhs.as_slice());
    }
}
