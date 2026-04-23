//! Tree-shaped constraint DSLs.
//!
//! Boundary types between the AST `Constraint` tree and its EDN form. Refs in
//! the tree carry either an integer index or a symbolic id; resolution to /
//! from the `AtomIdx` / `BondIdx` / ... on the AST is a separate fallible
//! step that consults the surrounding `Metadata`.

use umol_edn::{DeError, Edn, EdnKeyword, EdnMap, FromEdn, ToEdn};

use super::aromatic::AromaticSystemConstraintDsl;
use super::atom::AtomConstraintDsl;
use super::bond::BondConstraintDsl;
use super::dative::DativeBondConstraintDsl;
use super::error::ParseError;
use super::molecule::{Metadata, MoleculeDsl};
use super::multicenter::MulticenterBondConstraintDsl;
use super::noncovalent::NoncovalentBondConstraintDsl;
use super::value::ValueDsl;
use crate::ast::constraint::{Constraint, Constraints, MoleculeConstraint, SubPatternAnchor};
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::spin::SpinStateAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::dsl::config::MoleculeDefaults;

/// Resolution context for turning DSL refs (`AtomRef`, `BondRef`, ...) into
/// AST indices. Bundles the per-entity counts (for numeric index bounds
/// checking) with a borrowed `Metadata` (for id → index lookup). Used
/// uniformly by all constraint DSLs that contain refs, both on the
/// `FromAst` direction (only `metadata` is read — counts are ignored) and on
/// the `IntoAst` direction (both counts and `metadata` are read).
pub struct ResolveContext<'a> {
    pub atom_count: usize,
    pub bond_count: usize,
    pub dative_bond_count: usize,
    pub aromatic_system_count: usize,
    pub multicenter_bond_count: usize,
    pub noncovalent_bond_count: usize,
    pub metadata: &'a Metadata,
}

impl<'a> ResolveContext<'a> {
    /// Construct a context with only `metadata` populated; all counts zero.
    /// Appropriate for the `FromAst` direction, where counts are not read.
    pub fn for_rendering(metadata: &'a Metadata) -> Self {
        Self {
            atom_count: 0,
            bond_count: 0,
            dative_bond_count: 0,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
            metadata,
        }
    }
}

macro_rules! define_ref {
    ($name:ident, $idx:ident, $field:ident, $kind:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
        }

        impl $name {
            /// Build a ref from an AST index, preferring an id from `metadata`
            /// if one is recorded for this index.
            pub fn from_ast(idx: $idx, metadata: &Metadata) -> Self {
                if let Some(name) = metadata.$field.get(&idx) {
                    Self::Id(name.clone())
                } else {
                    Self::Index(idx.index())
                }
            }

            /// Resolve this ref to an AST index against `metadata`. Fails on
            /// unknown id or out-of-range numeric index.
            pub fn into_ast(self, count: usize, metadata: &Metadata) -> Result<$idx, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($idx::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => metadata
                        .$field
                        .iter()
                        .find(|(_, n)| n.as_str() == name)
                        .map(|(idx, _)| *idx)
                        .ok_or(ParseError::InvalidRef {
                            kind: $kind,
                            value: name,
                        }),
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
                    other => Err(DeError::TypeMismatch {
                        expected: concat!($kind, " ref (int or keyword)"),
                        got: other.kind(),
                        path: Vec::new(),
                    }),
                }
            }
        }

        impl ToEdn for $name {
            fn to_edn(&self) -> Edn<'static> {
                match self {
                    Self::Index(i) => Edn::Int(*i as i64),
                    Self::Id(name) => Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone())),
                }
            }
        }
    };
}

define_ref!(AtomRef, AtomIdx, atom_ids, "atom");
define_ref!(BondRef, BondIdx, bond_ids, "bond");
define_ref!(DativeBondRef, DativeBondIdx, dative_bond_ids, "dative-bond");
define_ref!(
    AromaticSystemRef,
    AromaticSystemIdx,
    aromatic_system_ids,
    "aromatic-system"
);
define_ref!(
    MulticenterBondRef,
    MulticenterBondIdx,
    multicenter_bond_ids,
    "multicenter-bond"
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondIdx,
    noncovalent_bond_ids,
    "noncovalent-bond"
);

// -- MoleculeConstraintDsl -------------------

/// Surface DSL wrapper around `MoleculeConstraint`. EDN form is a single-key
/// map keyed by the variant: `{:charge-sum {...}}`, `{:spin-sum {...}}`,
/// `{:bond-order-sum {...}}`, `{:connected [...]}`, `{:sub-pattern {...}}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraintDsl {
    ChargeSum {
        atoms: Vec<AtomRef>,
        sum: ValueDsl,
    },
    SpinSum {
        atoms: Vec<AtomRef>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Vec<BondRef>,
        sum: ValueDsl,
    },
    Connected(Vec<AtomRef>),
    SubPattern {
        anchor: SubPatternAnchorDsl,
        pattern: Box<MoleculeDsl>,
    },
}

impl FromAst<MoleculeConstraint> for MoleculeConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn from_ast<'a>(c: &MoleculeConstraint, ctx: &Self::Ctx<'a>) -> Result<Self, ParseError> {
        let meta = ctx.metadata;
        Ok(match c {
            MoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum {
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
                sum: ValueDsl::from_ast(sum, &()).unwrap(),
            },
            MoleculeConstraint::SpinSum { atoms, spin } => Self::SpinSum {
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
                spin: spin.clone(),
            },
            MoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum {
                bonds: bonds.iter().map(|&b| BondRef::from_ast(b, meta)).collect(),
                sum: ValueDsl::from_ast(sum, &()).unwrap(),
            },
            MoleculeConstraint::Connected(atoms) => {
                Self::Connected(atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect())
            }
            MoleculeConstraint::SubPattern { anchor, pattern } => {
                let pattern_dsl =
                    MoleculeDsl::from_ast(pattern.as_ref(), &MoleculeDefaults::zeroed())?;
                let anchor_dsl =
                    SubPatternAnchorDsl::from_ast_pair(anchor, meta, pattern_dsl.metadata());
                Self::SubPattern {
                    anchor: anchor_dsl,
                    pattern: Box::new(pattern_dsl),
                }
            }
        })
    }
}

impl IntoAst<MoleculeConstraint> for MoleculeConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn into_ast<'a>(self, ctx: &Self::Ctx<'a>) -> Result<MoleculeConstraint, ParseError> {
        let meta = ctx.metadata;
        Ok(match self {
            Self::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms: atoms
                    .into_iter()
                    .map(|r| r.into_ast(ctx.atom_count, meta))
                    .collect::<Result<_, _>>()?,
                sum: sum.into_ast(&()).unwrap(),
            },
            Self::SpinSum { atoms, spin } => MoleculeConstraint::SpinSum {
                atoms: atoms
                    .into_iter()
                    .map(|r| r.into_ast(ctx.atom_count, meta))
                    .collect::<Result<_, _>>()?,
                spin,
            },
            Self::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds: bonds
                    .into_iter()
                    .map(|r| r.into_ast(ctx.bond_count, meta))
                    .collect::<Result<_, _>>()?,
                sum: sum.into_ast(&()).unwrap(),
            },
            Self::Connected(atoms) => MoleculeConstraint::Connected(
                atoms
                    .into_iter()
                    .map(|r| r.into_ast(ctx.atom_count, meta))
                    .collect::<Result<_, _>>()?,
            ),
            Self::SubPattern { anchor, pattern } => {
                let pattern_ast = (*pattern).into_ast(&MoleculeDefaults::zeroed())?;
                let pattern_ctx = ResolveContext {
                    atom_count: pattern_ast.atom_count(),
                    bond_count: pattern_ast.bond_count(),
                    dative_bond_count: pattern_ast.dative_bond_count(),
                    aromatic_system_count: pattern_ast.aromatic_system_count(),
                    multicenter_bond_count: pattern_ast.multicenter_bond_count(),
                    noncovalent_bond_count: pattern_ast.noncovalent_bond_count(),
                    metadata: ctx.metadata, // pattern metadata is lost after into_ast
                };
                let anchor_ast = anchor.into_ast_pair(ctx, &pattern_ctx)?;
                MoleculeConstraint::SubPattern {
                    anchor: anchor_ast,
                    pattern: Box::new(pattern_ast),
                }
            }
        })
    }
}

impl<'de> FromEdn<'de> for MoleculeConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "molecule-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "molecule-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["molecule-constraint".into()],
            });
        };
        Ok(match key.name() {
            "charge-sum" => {
                let (atoms, sum) = parse_sum_map::<AtomRef>(v, "charge-sum", "atoms")?;
                Self::ChargeSum { atoms, sum }
            }
            "spin-sum" => {
                let m = expect_map(v, "spin-sum")?;
                let atoms_edn = m
                    .get_keyword("atoms")
                    .ok_or_else(|| DeError::MissingField {
                        key: "atoms".into(),
                        path: vec!["spin-sum".into()],
                    })?;
                let spin_edn = m.get_keyword("spin").ok_or_else(|| DeError::MissingField {
                    key: "spin".into(),
                    path: vec!["spin-sum".into()],
                })?;
                Self::SpinSum {
                    atoms: parse_refs::<AtomRef>(atoms_edn)?,
                    spin: parse_spin(spin_edn)?,
                }
            }
            "bond-order-sum" => {
                let (bonds, sum) = parse_sum_map::<BondRef>(v, "bond-order-sum", "bonds")?;
                Self::BondOrderSum { bonds, sum }
            }
            "connected" => Self::Connected(parse_refs::<AtomRef>(v)?),
            "sub-pattern" => {
                let m = expect_map(v, "sub-pattern")?;
                let anchor_edn = m
                    .get_keyword("anchor")
                    .ok_or_else(|| DeError::MissingField {
                        key: "anchor".into(),
                        path: vec!["sub-pattern".into()],
                    })?;
                let pattern_edn =
                    m.get_keyword("pattern")
                        .ok_or_else(|| DeError::MissingField {
                            key: "pattern".into(),
                            path: vec!["sub-pattern".into()],
                        })?;
                Self::SubPattern {
                    anchor: SubPatternAnchorDsl::from_edn(anchor_edn)?,
                    pattern: Box::new(MoleculeDsl::from_edn(pattern_edn)?),
                }
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for MoleculeConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::ChargeSum { atoms, sum } => ("charge-sum", render_sum_map("atoms", atoms, sum)),
            Self::SpinSum { atoms, spin } => {
                let mut m = EdnMap::with_capacity(2);
                m.insert(Edn::keyword("atoms"), render_refs(atoms));
                m.insert(Edn::keyword("spin"), render_spin(spin));
                ("spin-sum", Edn::Map(m))
            }
            Self::BondOrderSum { bonds, sum } => {
                ("bond-order-sum", render_sum_map("bonds", bonds, sum))
            }
            Self::Connected(atoms) => ("connected", render_refs(atoms)),
            Self::SubPattern { anchor, pattern } => {
                let mut m = EdnMap::with_capacity(2);
                m.insert(Edn::keyword("anchor"), anchor.to_edn());
                m.insert(Edn::keyword("pattern"), pattern.to_edn());
                ("sub-pattern", Edn::Map(m))
            }
        };
        let mut outer = EdnMap::with_capacity(1);
        outer.insert(Edn::Keyword(EdnKeyword::owned(key.into())), value);
        Edn::Map(outer)
    }
}

// -- SubPatternAnchorDsl ------------------

/// Surface DSL wrapper around `SubPatternAnchor`. Each vector carries
/// `(target, pattern)` ref pairs. Target-side refs resolve against the outer
/// molecule's `Metadata`; pattern-side refs resolve against the pattern
/// molecule's `Metadata`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SubPatternAnchorDsl {
    pub atoms: Vec<(AtomRef, AtomRef)>,
    pub bonds: Vec<(BondRef, BondRef)>,
    pub dative_bonds: Vec<(DativeBondRef, DativeBondRef)>,
    pub aromatic_systems: Vec<(AromaticSystemRef, AromaticSystemRef)>,
    pub multicenter_bonds: Vec<(MulticenterBondRef, MulticenterBondRef)>,
    pub noncovalent_bonds: Vec<(NoncovalentBondRef, NoncovalentBondRef)>,
}

impl SubPatternAnchorDsl {
    /// Build from an AST anchor. `target_meta` is the outer molecule's
    /// metadata; `pattern_meta` is the pattern molecule's metadata.
    pub fn from_ast_pair(
        anchor: &SubPatternAnchor,
        target_meta: &Metadata,
        pattern_meta: &Metadata,
    ) -> Self {
        Self {
            atoms: anchor
                .atoms()
                .iter()
                .map(|&(t, p)| {
                    (
                        AtomRef::from_ast(t, target_meta),
                        AtomRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            bonds: anchor
                .bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        BondRef::from_ast(t, target_meta),
                        BondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            dative_bonds: anchor
                .dative_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        DativeBondRef::from_ast(t, target_meta),
                        DativeBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            aromatic_systems: anchor
                .aromatic_systems()
                .iter()
                .map(|&(t, p)| {
                    (
                        AromaticSystemRef::from_ast(t, target_meta),
                        AromaticSystemRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            multicenter_bonds: anchor
                .multicenter_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        MulticenterBondRef::from_ast(t, target_meta),
                        MulticenterBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            noncovalent_bonds: anchor
                .noncovalent_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        NoncovalentBondRef::from_ast(t, target_meta),
                        NoncovalentBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
        }
    }

    /// Resolve to an AST anchor. `target_ctx` carries outer-molecule counts +
    /// metadata; `pattern_ctx` carries pattern-molecule counts + metadata.
    pub fn into_ast_pair(
        self,
        target_ctx: &ResolveContext,
        pattern_ctx: &ResolveContext,
    ) -> Result<SubPatternAnchor, ParseError> {
        let mut anchor = SubPatternAnchor::new();
        for (t, p) in self.atoms {
            anchor.push_atom(
                t.into_ast(target_ctx.atom_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.atom_count, pattern_ctx.metadata)?,
            );
        }
        for (t, p) in self.bonds {
            anchor.push_bond(
                t.into_ast(target_ctx.bond_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.bond_count, pattern_ctx.metadata)?,
            );
        }
        for (t, p) in self.dative_bonds {
            anchor.push_dative_bond(
                t.into_ast(target_ctx.dative_bond_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.dative_bond_count, pattern_ctx.metadata)?,
            );
        }
        for (t, p) in self.aromatic_systems {
            anchor.push_aromatic_system(
                t.into_ast(target_ctx.aromatic_system_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.aromatic_system_count, pattern_ctx.metadata)?,
            );
        }
        for (t, p) in self.multicenter_bonds {
            anchor.push_multicenter_bond(
                t.into_ast(target_ctx.multicenter_bond_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.multicenter_bond_count, pattern_ctx.metadata)?,
            );
        }
        for (t, p) in self.noncovalent_bonds {
            anchor.push_noncovalent_bond(
                t.into_ast(target_ctx.noncovalent_bond_count, target_ctx.metadata)?,
                p.into_ast(pattern_ctx.noncovalent_bond_count, pattern_ctx.metadata)?,
            );
        }
        Ok(anchor)
    }
}

impl<'de> FromEdn<'de> for SubPatternAnchorDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let m = expect_map(edn, "sub-pattern-anchor")?;
        let mut out = Self::default();
        for (k, v) in m.iter() {
            let Edn::Keyword(key) = k else {
                return Err(DeError::TypeMismatch {
                    expected: "keyword key",
                    got: k.kind(),
                    path: vec!["sub-pattern-anchor".into()],
                });
            };
            match key.name() {
                "atoms" => out.atoms = parse_ref_pairs::<AtomRef, AtomRef>(v)?,
                "bonds" => out.bonds = parse_ref_pairs::<BondRef, BondRef>(v)?,
                "dative-bonds" => {
                    out.dative_bonds = parse_ref_pairs::<DativeBondRef, DativeBondRef>(v)?
                }
                "aromatic-systems" => {
                    out.aromatic_systems =
                        parse_ref_pairs::<AromaticSystemRef, AromaticSystemRef>(v)?
                }
                "multicenter-bonds" => {
                    out.multicenter_bonds =
                        parse_ref_pairs::<MulticenterBondRef, MulticenterBondRef>(v)?
                }
                "noncovalent-bonds" => {
                    out.noncovalent_bonds =
                        parse_ref_pairs::<NoncovalentBondRef, NoncovalentBondRef>(v)?
                }
                other => {
                    return Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["sub-pattern-anchor".into()],
                    });
                }
            }
        }
        Ok(out)
    }
}

impl ToEdn for SubPatternAnchorDsl {
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(6);
        if !self.atoms.is_empty() {
            m.insert(Edn::keyword("atoms"), render_ref_pairs(&self.atoms));
        }
        if !self.bonds.is_empty() {
            m.insert(Edn::keyword("bonds"), render_ref_pairs(&self.bonds));
        }
        if !self.dative_bonds.is_empty() {
            m.insert(
                Edn::keyword("dative-bonds"),
                render_ref_pairs(&self.dative_bonds),
            );
        }
        if !self.aromatic_systems.is_empty() {
            m.insert(
                Edn::keyword("aromatic-systems"),
                render_ref_pairs(&self.aromatic_systems),
            );
        }
        if !self.multicenter_bonds.is_empty() {
            m.insert(
                Edn::keyword("multicenter-bonds"),
                render_ref_pairs(&self.multicenter_bonds),
            );
        }
        if !self.noncovalent_bonds.is_empty() {
            m.insert(
                Edn::keyword("noncovalent-bonds"),
                render_ref_pairs(&self.noncovalent_bonds),
            );
        }
        Edn::Map(m)
    }
}

// -- Helpers ------------------

fn expect_map<'e>(edn: &'e Edn<'e>, context: &'static str) -> Result<&'e EdnMap<'e>, DeError> {
    match edn {
        Edn::Map(m) => Ok(m),
        other => Err(DeError::TypeMismatch {
            expected: "map",
            got: other.kind(),
            path: vec![context.into()],
        }),
    }
}

fn parse_refs<R>(edn: &Edn<'_>) -> Result<Vec<R>, DeError>
where
    R: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of refs",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    v.iter().map(R::from_edn).collect()
}

fn render_refs<R: ToEdn>(refs: &[R]) -> Edn<'static> {
    Edn::Vector(refs.iter().map(R::to_edn).collect::<Vec<_>>().into())
}

fn parse_ref_pairs<A, B>(edn: &Edn<'_>) -> Result<Vec<(A, B)>, DeError>
where
    A: for<'de> FromEdn<'de>,
    B: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of [target pattern] pairs",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    v.iter()
        .map(|e| {
            let Edn::Vector(pair) = e else {
                return Err(DeError::TypeMismatch {
                    expected: "2-element vector [target pattern]",
                    got: e.kind(),
                    path: Vec::new(),
                });
            };
            if pair.len() != 2 {
                return Err(DeError::Custom(format!(
                    "anchor pair must have 2 elements, got {}",
                    pair.len()
                )));
            }
            Ok((A::from_edn(&pair[0])?, B::from_edn(&pair[1])?))
        })
        .collect()
}

fn render_ref_pairs<A: ToEdn, B: ToEdn>(pairs: &[(A, B)]) -> Edn<'static> {
    Edn::Vector(
        pairs
            .iter()
            .map(|(a, b)| Edn::Vector(vec![a.to_edn(), b.to_edn()].into()))
            .collect::<Vec<_>>()
            .into(),
    )
}

fn parse_sum_map<R>(
    edn: &Edn<'_>,
    context: &'static str,
    refs_key: &'static str,
) -> Result<(Vec<R>, ValueDsl), DeError>
where
    R: for<'de> FromEdn<'de>,
{
    let m = expect_map(edn, context)?;
    let refs_edn = m
        .get_keyword(refs_key)
        .ok_or_else(|| DeError::MissingField {
            key: refs_key.to_string(),
            path: vec![context.into()],
        })?;
    let sum_edn = m.get_keyword("sum").ok_or_else(|| DeError::MissingField {
        key: "sum".into(),
        path: vec![context.into()],
    })?;
    Ok((parse_refs::<R>(refs_edn)?, ValueDsl::from_edn(sum_edn)?))
}

fn render_sum_map<R: ToEdn>(refs_key: &str, refs: &[R], sum: &ValueDsl) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(refs_key.into())),
        render_refs(refs),
    );
    m.insert(Edn::keyword("sum"), sum.to_edn());
    Edn::Map(m)
}

fn parse_spin(edn: &Edn<'_>) -> Result<SpinStateAst, DeError> {
    let m = expect_map(edn, "spin")?;
    let unpaired = m
        .get_keyword("unpaired")
        .ok_or_else(|| DeError::MissingField {
            key: "unpaired".into(),
            path: vec!["spin".into()],
        })?;
    let multiplicity = m
        .get_keyword("multiplicity")
        .ok_or_else(|| DeError::MissingField {
            key: "multiplicity".into(),
            path: vec!["spin".into()],
        })?;
    Ok(SpinStateAst::from_values(
        ValueDsl::from_edn(unpaired)?.into_ast(&()).unwrap(),
        ValueDsl::from_edn(multiplicity)?.into_ast(&()).unwrap(),
    ))
}

fn render_spin(spin: &SpinStateAst) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    m.insert(
        Edn::keyword("unpaired"),
        ValueDsl::from_ast(&spin.unpaired, &()).unwrap().to_edn(),
    );
    m.insert(
        Edn::keyword("multiplicity"),
        ValueDsl::from_ast(&spin.multiplicity, &())
            .unwrap()
            .to_edn(),
    );
    Edn::Map(m)
}

// -- ConstraintDsl ------------------

/// Surface DSL wrapper around `Constraint`. Single-key-map EDN form:
/// `{:atom [<ref> <atom-constraint>]}` for entity leaves;
/// `{:charge-sum {...}}` etc. for molecule-scope leaves (keys flattened from
/// `MoleculeConstraintDsl`); `{:and [...]}` / `{:or [...]}` / `{:not <c>}`
/// for combinators.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstraintDsl {
    Atom(AtomRef, AtomConstraintDsl),
    Bond(BondRef, BondConstraintDsl),
    DativeBond(DativeBondRef, DativeBondConstraintDsl),
    AromaticSystem(AromaticSystemRef, AromaticSystemConstraintDsl),
    MulticenterBond(MulticenterBondRef, MulticenterBondConstraintDsl),
    NoncovalentBond(NoncovalentBondRef, NoncovalentBondConstraintDsl),
    Molecule(MoleculeConstraintDsl),
    And(Vec<ConstraintDsl>),
    Or(Vec<ConstraintDsl>),
    Not(Box<ConstraintDsl>),
}

impl<'de> FromEdn<'de> for ConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let m = expect_map(edn, "constraint")?;
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["constraint".into()],
            });
        };
        Ok(match key.name() {
            "atom" => {
                let (r, c) = parse_entity_leaf::<AtomRef, AtomConstraintDsl>(v, "atom")?;
                Self::Atom(r, c)
            }
            "bond" => {
                let (r, c) = parse_entity_leaf::<BondRef, BondConstraintDsl>(v, "bond")?;
                Self::Bond(r, c)
            }
            "dative-bond" => {
                let (r, c) =
                    parse_entity_leaf::<DativeBondRef, DativeBondConstraintDsl>(v, "dative-bond")?;
                Self::DativeBond(r, c)
            }
            "aromatic-system" => {
                let (r, c) = parse_entity_leaf::<AromaticSystemRef, AromaticSystemConstraintDsl>(
                    v,
                    "aromatic-system",
                )?;
                Self::AromaticSystem(r, c)
            }
            "multicenter-bond" => {
                let (r, c) = parse_entity_leaf::<MulticenterBondRef, MulticenterBondConstraintDsl>(
                    v,
                    "multicenter-bond",
                )?;
                Self::MulticenterBond(r, c)
            }
            "noncovalent-bond" => {
                let (r, c) = parse_entity_leaf::<NoncovalentBondRef, NoncovalentBondConstraintDsl>(
                    v,
                    "noncovalent-bond",
                )?;
                Self::NoncovalentBond(r, c)
            }
            "and" => Self::And(parse_constraint_vec(v, "and")?),
            "or" => Self::Or(parse_constraint_vec(v, "or")?),
            "not" => Self::Not(Box::new(ConstraintDsl::from_edn(v)?)),
            // Molecule-scope keys: delegate to MoleculeConstraintDsl.
            "charge-sum" | "spin-sum" | "bond-order-sum" | "connected" | "sub-pattern" => {
                Self::Molecule(MoleculeConstraintDsl::from_edn(edn)?)
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for ConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Atom(r, c) => entity_leaf_edn("atom", r, c),
            Self::Bond(r, c) => entity_leaf_edn("bond", r, c),
            Self::DativeBond(r, c) => entity_leaf_edn("dative-bond", r, c),
            Self::AromaticSystem(r, c) => entity_leaf_edn("aromatic-system", r, c),
            Self::MulticenterBond(r, c) => entity_leaf_edn("multicenter-bond", r, c),
            Self::NoncovalentBond(r, c) => entity_leaf_edn("noncovalent-bond", r, c),
            Self::Molecule(m) => m.to_edn(),
            Self::And(xs) => combinator_edn("and", xs),
            Self::Or(xs) => combinator_edn("or", xs),
            Self::Not(c) => {
                let mut m = EdnMap::with_capacity(1);
                m.insert(Edn::keyword("not"), c.to_edn());
                Edn::Map(m)
            }
        }
    }
}

impl FromAst<Constraint> for ConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn from_ast<'a>(c: &Constraint, ctx: &Self::Ctx<'a>) -> Result<Self, ParseError> {
        let meta = ctx.metadata;
        Ok(match c {
            Constraint::Atom(idx, c) => Self::Atom(
                AtomRef::from_ast(*idx, meta),
                AtomConstraintDsl::from_ast(c, &()).unwrap(),
            ),
            Constraint::Bond(idx, c) => Self::Bond(
                BondRef::from_ast(*idx, meta),
                BondConstraintDsl::from_ast(c, &()).unwrap(),
            ),
            Constraint::DativeBond(idx, c) => Self::DativeBond(
                DativeBondRef::from_ast(*idx, meta),
                DativeBondConstraintDsl::from_ast(c, ctx).unwrap(),
            ),
            Constraint::AromaticSystem(idx, c) => Self::AromaticSystem(
                AromaticSystemRef::from_ast(*idx, meta),
                AromaticSystemConstraintDsl::from_ast(c, ctx).unwrap(),
            ),
            Constraint::MulticenterBond(idx, c) => Self::MulticenterBond(
                MulticenterBondRef::from_ast(*idx, meta),
                MulticenterBondConstraintDsl::from_ast(c, ctx).unwrap(),
            ),
            Constraint::NoncovalentBond(idx, c) => Self::NoncovalentBond(
                NoncovalentBondRef::from_ast(*idx, meta),
                NoncovalentBondConstraintDsl::from_ast(c, ctx).unwrap(),
            ),
            Constraint::Molecule(m) => Self::Molecule(MoleculeConstraintDsl::from_ast(m, ctx)?),
            Constraint::And(xs) => Self::And(
                xs.iter()
                    .map(|c| ConstraintDsl::from_ast(c, ctx))
                    .collect::<Result<_, _>>()?,
            ),
            Constraint::Or(xs) => Self::Or(
                xs.iter()
                    .map(|c| ConstraintDsl::from_ast(c, ctx))
                    .collect::<Result<_, _>>()?,
            ),
            Constraint::Not(c) => Self::Not(Box::new(ConstraintDsl::from_ast(c, ctx)?)),
        })
    }
}

impl IntoAst<Constraint> for ConstraintDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn into_ast<'a>(self, ctx: &Self::Ctx<'a>) -> Result<Constraint, ParseError> {
        let meta = ctx.metadata;
        Ok(match self {
            Self::Atom(r, c) => {
                Constraint::Atom(r.into_ast(ctx.atom_count, meta)?, c.into_ast(&()).unwrap())
            }
            Self::Bond(r, c) => {
                Constraint::Bond(r.into_ast(ctx.bond_count, meta)?, c.into_ast(&()).unwrap())
            }
            Self::DativeBond(r, c) => {
                Constraint::DativeBond(r.into_ast(ctx.dative_bond_count, meta)?, c.into_ast(ctx)?)
            }
            Self::AromaticSystem(r, c) => Constraint::AromaticSystem(
                r.into_ast(ctx.aromatic_system_count, meta)?,
                c.into_ast(ctx)?,
            ),
            Self::MulticenterBond(r, c) => Constraint::MulticenterBond(
                r.into_ast(ctx.multicenter_bond_count, meta)?,
                c.into_ast(ctx)?,
            ),
            Self::NoncovalentBond(r, c) => Constraint::NoncovalentBond(
                r.into_ast(ctx.noncovalent_bond_count, meta)?,
                c.into_ast(ctx)?,
            ),
            Self::Molecule(m) => Constraint::Molecule(m.into_ast(ctx)?),
            Self::And(xs) => Constraint::And(
                xs.into_iter()
                    .map(|c| c.into_ast(ctx))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Or(xs) => Constraint::Or(
                xs.into_iter()
                    .map(|c| c.into_ast(ctx))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Not(c) => Constraint::Not(Box::new(c.into_ast(ctx)?)),
        })
    }
}

// -- ConstraintsDsl ------------------

/// Surface DSL wrapper around `Constraints` (a flat vec of `Constraint`).
/// EDN form: a vector of `ConstraintDsl` EDN forms.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConstraintsDsl(pub Vec<ConstraintDsl>);

impl<'de> FromEdn<'de> for ConstraintsDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Vector(v) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "vector of constraints",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        Ok(Self(
            v.iter()
                .map(ConstraintDsl::from_edn)
                .collect::<Result<_, _>>()?,
        ))
    }
}

impl ToEdn for ConstraintsDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(self.0.iter().map(|c| c.to_edn()).collect::<Vec<_>>().into())
    }
}

impl FromAst<Constraints> for ConstraintsDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn from_ast<'a>(cs: &Constraints, ctx: &Self::Ctx<'a>) -> Result<Self, ParseError> {
        Ok(Self(
            cs.iter()
                .map(|c| ConstraintDsl::from_ast(c, ctx))
                .collect::<Result<_, _>>()?,
        ))
    }
}

impl IntoAst<Constraints> for ConstraintsDsl {
    type Ctx<'a> = ResolveContext<'a>;
    type Error = ParseError;

    fn into_ast<'a>(self, ctx: &Self::Ctx<'a>) -> Result<Constraints, ParseError> {
        let mut out = Constraints::new();
        for c in self.0 {
            out.push(c.into_ast(ctx)?);
        }
        Ok(out)
    }
}

fn parse_entity_leaf<R, C>(edn: &Edn<'_>, context: &'static str) -> Result<(R, C), DeError>
where
    R: for<'de> FromEdn<'de>,
    C: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "2-element vector [ref constraint]",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    if v.len() != 2 {
        return Err(DeError::Custom(format!(
            "{} entity leaf must have 2 elements, got {}",
            context,
            v.len()
        )));
    }
    Ok((R::from_edn(&v[0])?, C::from_edn(&v[1])?))
}

fn entity_leaf_edn<R: ToEdn, C: ToEdn>(key: &str, r: &R, c: &C) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(key.into())),
        Edn::Vector(vec![r.to_edn(), c.to_edn()].into()),
    );
    Edn::Map(m)
}

fn parse_constraint_vec(
    edn: &Edn<'_>,
    context: &'static str,
) -> Result<Vec<ConstraintDsl>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of constraints",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    v.iter().map(ConstraintDsl::from_edn).collect()
}

fn combinator_edn(key: &str, xs: &[ConstraintDsl]) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(key.into())),
        Edn::Vector(xs.iter().map(|c| c.to_edn()).collect::<Vec<_>>().into()),
    );
    Edn::Map(m)
}

#[cfg(test)]
mod tests {
    use bimap::BiMap;
    use indexmap::IndexMap;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnKeyword};

    use super::*;

    #[fixture]
    fn meta_with_atom_id() -> Metadata {
        let mut atom_ids = IndexMap::new();
        atom_ids.insert(AtomIdx(2), "c1".to_string());
        Metadata {
            atom_ids,
            atom_aliases: BiMap::new(),
            bond_ids: IndexMap::new(),
            dative_bond_ids: IndexMap::new(),
            aromatic_system_ids: IndexMap::new(),
            multicenter_bond_ids: IndexMap::new(),
            noncovalent_bond_ids: IndexMap::new(),
        }
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
    fn test_atom_ref_from_ast_uses_id_when_present(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomIdx(2), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Id("c1".into()));
    }

    #[rstest]
    fn test_atom_ref_from_ast_falls_back_to_index_without_id(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomIdx(4), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Index(4));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_id(meta_with_atom_id: Metadata) {
        let idx = AtomRef::Id("c1".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(idx, AtomIdx(2));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_index(meta_with_atom_id: Metadata) {
        let idx = AtomRef::Index(3).into_ast(5, &meta_with_atom_id).unwrap();
        assert_eq!(idx, AtomIdx(3));
    }

    #[rstest]
    fn test_atom_ref_into_ast_out_of_range_index(meta_with_atom_id: Metadata) {
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
    fn test_atom_ref_into_ast_unknown_id(meta_with_atom_id: Metadata) {
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

    // -- MoleculeConstraintDsl ----------------

    use crate::ast::value::ValueAst;

    fn empty_metadata() -> Metadata {
        Metadata {
            atom_ids: IndexMap::new(),
            atom_aliases: BiMap::new(),
            bond_ids: IndexMap::new(),
            dative_bond_ids: IndexMap::new(),
            aromatic_system_ids: IndexMap::new(),
            multicenter_bond_ids: IndexMap::new(),
            noncovalent_bond_ids: IndexMap::new(),
        }
    }

    fn full_ctx<'a>(meta: &'a Metadata) -> ResolveContext<'a> {
        ResolveContext {
            atom_count: 10,
            bond_count: 10,
            dative_bond_count: 10,
            aromatic_system_count: 10,
            multicenter_bond_count: 10,
            noncovalent_bond_count: 10,
            metadata: meta,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum(
        MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0), AtomIdx(1)],
            sum: ValueAst::Lit(0),
        },
        "{:charge-sum {:atoms [0 1] :sum 0}}"
    )]
    #[case::spin_sum(
        MoleculeConstraint::SpinSum {
            atoms: vec![AtomIdx(0)],
            spin: SpinStateAst::new(1, 2),
        },
        "{:spin-sum {:atoms [0] :spin {:unpaired 1 :multiplicity 2}}}"
    )]
    #[case::bond_order_sum(
        MoleculeConstraint::BondOrderSum {
            bonds: vec![BondIdx(0), BondIdx(1)],
            sum: ValueAst::Lit(4),
        },
        "{:bond-order-sum {:bonds [0 1] :sum 4}}"
    )]
    #[case::connected(
        MoleculeConstraint::Connected(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
        "{:connected [0 1 2]}"
    )]
    fn test_molecule_constraint_dsl_roundtrip(
        #[case] input: MoleculeConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = empty_metadata();
        let render_ctx = ResolveContext::for_rendering(&meta);
        let dsl = MoleculeConstraintDsl::from_ast(&input, &render_ctx).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = MoleculeConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&full_ctx(&meta)).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_molecule_constraint_dsl_rejects_wrong_shape() {
        let err = MoleculeConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_molecule_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = MoleculeConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_molecule_constraint_dsl_charge_sum_rejects_missing_sum() {
        let edn = read_string("{:charge-sum {:atoms [0 1]}}").unwrap();
        let err = MoleculeConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::MissingField { .. }));
    }

    // -- SubPatternAnchorDsl ----------------

    #[rstest]
    fn test_sub_pattern_anchor_dsl_empty_roundtrip() {
        let meta = empty_metadata();
        let anchor = SubPatternAnchor::new();
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        // Empty anchor renders as an empty map.
        assert_eq!(edn, read_string("{}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let ctx = full_ctx(&meta);
        let back = parsed.into_ast_pair(&ctx, &ctx).unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_atoms_roundtrip() {
        let meta = empty_metadata();
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomIdx(3), AtomIdx(0));
        anchor.push_atom(AtomIdx(5), AtomIdx(1));
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("{:atoms [[3 0] [5 1]]}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let ctx = full_ctx(&meta);
        let back = parsed.into_ast_pair(&ctx, &ctx).unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus [[0 0]]}").unwrap();
        let err = SubPatternAnchorDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_rejects_wrong_pair_length() {
        let edn = read_string("{:atoms [[0]]}").unwrap();
        let err = SubPatternAnchorDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    // -- ConstraintDsl ----------------

    use crate::ast::constraint::{AtomConstraint, BondConstraint};

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_leaf(
        Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
        "{:atom [0 {:valence 4}]}"
    )]
    #[case::bond_leaf(
        Constraint::Bond(BondIdx(1), BondConstraint::Aromatic),
        "{:bond [1 :aromatic]}"
    )]
    #[case::molecule_connected(
        Constraint::Molecule(MoleculeConstraint::Connected(vec![AtomIdx(0), AtomIdx(1)])),
        "{:connected [0 1]}"
    )]
    #[case::not(
        Constraint::Not(Box::new(Constraint::Atom(
            AtomIdx(0),
            AtomConstraint::Valence(ValueAst::Lit(3)),
        ))),
        "{:not {:atom [0 {:valence 3}]}}"
    )]
    #[case::and(
        Constraint::And(vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Bond(BondIdx(0), BondConstraint::Aromatic),
        ]),
        "{:and [{:atom [0 {:valence 4}]} {:bond [0 :aromatic]}]}"
    )]
    #[case::or(
        Constraint::Or(vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(3))),
            Constraint::Atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(4))),
        ]),
        "{:or [{:atom [0 {:degree 3}]} {:atom [0 {:degree 4}]}]}"
    )]
    fn test_constraint_dsl_roundtrip(
        #[case] input: Constraint,
        #[case] edn_source: &str,
    ) {
        let meta = empty_metadata();
        let ctx = full_ctx(&meta);
        let dsl = ConstraintDsl::from_ast(&input, &ctx).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = ConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&ctx).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = ConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    // -- ConstraintsDsl ----------------

    #[rstest]
    fn test_constraints_dsl_empty_roundtrip() {
        let meta = empty_metadata();
        let ctx = full_ctx(&meta);
        let cs = Constraints::new();
        let dsl = ConstraintsDsl::from_ast(&cs, &ctx).unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("[]").unwrap());
        let parsed = ConstraintsDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&ctx).unwrap();
        assert_eq!(back, cs);
    }

    #[rstest]
    fn test_constraints_dsl_multi_roundtrip() {
        let meta = empty_metadata();
        let ctx = full_ctx(&meta);
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(
            AtomIdx(0),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        ));
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
            AtomIdx(1),
        ])));
        let dsl = ConstraintsDsl::from_ast(&cs, &ctx).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string("[{:atom [0 {:valence 4}]} {:connected [0 1]}]").unwrap();
        assert_eq!(edn, expected);
        let parsed = ConstraintsDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&ctx).unwrap();
        assert_eq!(back, cs);
    }
}
