//! Surface DSL for cross-entity relational constraints.
//!
//! Mirrors [`crate::ast::constraint::RelationalConstraint`] with
//! `AtomRef`/`BondRef`/`DativeBondRef`/etc. in place of raw `AtomIdx` etc.
//! EDN form is a flat single-key map keyed by the variant — e.g.
//! `{:dative-bond-donor [<bond_ref> <atom_ref>]}`,
//! `{:aromatic-system-contains [<system_ref> <atom_ref>]}`. Payloads are
//! always a 2-element vector: `[<owning_entity_ref> <target>]`, where
//! `<target>` is either a single ref, a vector of refs, an atom constraint,
//! or a 2-element constraint pair.
//!
//! All ref-bearing entity constraints live here; the per-entity
//! `XxxConstraintDsl` types are narrowed to value-only inline variants.

use umol_edn::{DeError, Edn, EdnKeyword, EdnMap, FromEdn, ToEdn};

use super::atom::AtomConstraintDsl;
use super::constraint::{
    AromaticSystemRef, AtomRef, BondRef, DativeBondRef, EntityCounts, MulticenterBondRef,
    NoncovalentBondRef,
};
use super::error::ParseError;
use super::molecule::Metadata;
use crate::ast::constraint::RelationalConstraint;
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL wrapper around [`RelationalConstraint`]. Structural parallel
/// to the AST enum — same 18 variants, with surface refs ([`AtomRef`],
/// [`DativeBondRef`], etc.) in place of raw `*Idx`. Each variant's EDN
/// form is a flat single-key map `{:<entity>-<role> [<owner_ref> <target>]}`.
/// See the AST enum for per-variant semantics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RelationalConstraintDsl {
    // region: Dative bond
    /// EDN: `{:dative-bond-donor [<dative_ref> <atom_ref>]}`.
    DativeBondDonor { bond: DativeBondRef, atom: AtomRef },
    /// EDN: `{:dative-bond-acceptor [<dative_ref> <atom_ref>]}`.
    DativeBondAcceptor { bond: DativeBondRef, atom: AtomRef },
    /// EDN: `{:dative-bond-parallels [<dative_ref> <bond_ref>]}`.
    DativeBondParallels {
        dative: DativeBondRef,
        parallel: BondRef,
    },
    /// EDN: `{:dative-bond-donor-satisfies [<dative_ref> <atom-constraint>]}`.
    DativeBondDonorSatisfies {
        bond: DativeBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:dative-bond-acceptor-satisfies [<dative_ref> <atom-constraint>]}`.
    DativeBondAcceptorSatisfies {
        bond: DativeBondRef,
        predicate: Box<AtomConstraintDsl>,
    },

    // endregion: Dative bond

    // region: Aromatic system
    /// EDN: `{:aromatic-system-atoms [<system_ref> [<atom_ref>+]]}`.
    AromaticSystemAtoms {
        system: AromaticSystemRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:aromatic-system-contains [<system_ref> <atom_ref>]}`.
    AromaticSystemContains {
        system: AromaticSystemRef,
        atom: AtomRef,
    },
    /// EDN: `{:aromatic-system-contains-all [<system_ref> [<atom_ref>+]]}`.
    AromaticSystemContainsAll {
        system: AromaticSystemRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:aromatic-system-all-atoms [<system_ref> <atom-constraint>]}`.
    AromaticSystemAllAtoms {
        system: AromaticSystemRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:aromatic-system-any-atom [<system_ref> <atom-constraint>]}`.
    AromaticSystemAnyAtom {
        system: AromaticSystemRef,
        predicate: Box<AtomConstraintDsl>,
    },

    // endregion: Aromatic system

    // region: Multicenter bond
    /// EDN: `{:multicenter-bond-atoms [<bond_ref> [<atom_ref>+]]}`.
    MulticenterBondAtoms {
        bond: MulticenterBondRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:multicenter-bond-contains [<bond_ref> <atom_ref>]}`.
    MulticenterBondContains {
        bond: MulticenterBondRef,
        atom: AtomRef,
    },
    /// EDN: `{:multicenter-bond-contains-all [<bond_ref> [<atom_ref>+]]}`.
    MulticenterBondContainsAll {
        bond: MulticenterBondRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:multicenter-bond-all-atoms [<bond_ref> <atom-constraint>]}`.
    MulticenterBondAllAtoms {
        bond: MulticenterBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:multicenter-bond-any-atom [<bond_ref> <atom-constraint>]}`.
    MulticenterBondAnyAtom {
        bond: MulticenterBondRef,
        predicate: Box<AtomConstraintDsl>,
    },

    // endregion: Multicenter bond

    // region: Noncovalent bond
    /// EDN: `{:noncovalent-bond-ends [<bond_ref> [<atom_ref> <atom_ref>]]}`.
    NoncovalentBondEnds {
        bond: NoncovalentBondRef,
        atoms: [AtomRef; 2],
    },
    /// EDN: `{:noncovalent-bond-contains [<bond_ref> <atom_ref>]}`.
    NoncovalentBondContains {
        bond: NoncovalentBondRef,
        atom: AtomRef,
    },
    /// EDN: `{:noncovalent-bond-ends-satisfy [<bond_ref>
    /// [<atom-constraint> <atom-constraint>]]}`.
    NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondRef,
        predicates: [Box<AtomConstraintDsl>; 2],
    },
    // endregion: Noncovalent bond
}

/// Top-level EDN keywords for every relational variant. Matches the
/// flat-key convention — every relational form is its own `:<entity>-<role>`
/// keyword, rather than nested under the entity like narrow constraints.
pub(super) const RELATIONAL_KEYS: &[&str] = &[
    "dative-bond-donor",
    "dative-bond-acceptor",
    "dative-bond-parallels",
    "dative-bond-donor-satisfies",
    "dative-bond-acceptor-satisfies",
    "aromatic-system-atoms",
    "aromatic-system-contains",
    "aromatic-system-contains-all",
    "aromatic-system-all-atoms",
    "aromatic-system-any-atom",
    "multicenter-bond-atoms",
    "multicenter-bond-contains",
    "multicenter-bond-contains-all",
    "multicenter-bond-all-atoms",
    "multicenter-bond-any-atom",
    "noncovalent-bond-ends",
    "noncovalent-bond-contains",
    "noncovalent-bond-ends-satisfy",
];

impl<'de> FromEdn<'de> for RelationalConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "relational-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "relational-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["relational-constraint".into()],
            });
        };
        parse_payload(key.name(), v)
    }
}

impl ToEdn for RelationalConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = render_payload(self);
        let mut m = EdnMap::with_capacity(1);
        m.insert(Edn::Keyword(EdnKeyword::owned(key.into())), value);
        Edn::Map(m)
    }
}

impl RelationalConstraintDsl {
    pub(crate) fn from_ast(rel: &RelationalConstraint, meta: &Metadata) -> Self {
        use RelationalConstraint::*;
        match rel {
            DativeBondDonor { bond, atom } => Self::DativeBondDonor {
                bond: DativeBondRef::from_ast(*bond, meta),
                atom: AtomRef::from_ast(*atom, meta),
            },
            DativeBondAcceptor { bond, atom } => Self::DativeBondAcceptor {
                bond: DativeBondRef::from_ast(*bond, meta),
                atom: AtomRef::from_ast(*atom, meta),
            },
            DativeBondParallels { dative, parallel } => Self::DativeBondParallels {
                dative: DativeBondRef::from_ast(*dative, meta),
                parallel: BondRef::from_ast(*parallel, meta),
            },
            DativeBondDonorSatisfies { bond, predicate } => Self::DativeBondDonorSatisfies {
                bond: DativeBondRef::from_ast(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            DativeBondAcceptorSatisfies { bond, predicate } => Self::DativeBondAcceptorSatisfies {
                bond: DativeBondRef::from_ast(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            AromaticSystemAtoms { system, atoms } => Self::AromaticSystemAtoms {
                system: AromaticSystemRef::from_ast(*system, meta),
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
            },
            AromaticSystemContains { system, atom } => Self::AromaticSystemContains {
                system: AromaticSystemRef::from_ast(*system, meta),
                atom: AtomRef::from_ast(*atom, meta),
            },
            AromaticSystemContainsAll { system, atoms } => Self::AromaticSystemContainsAll {
                system: AromaticSystemRef::from_ast(*system, meta),
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
            },
            AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system: AromaticSystemRef::from_ast(*system, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system: AromaticSystemRef::from_ast(*system, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            MulticenterBondAtoms { bond, atoms } => Self::MulticenterBondAtoms {
                bond: MulticenterBondRef::from_ast(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
            },
            MulticenterBondContains { bond, atom } => Self::MulticenterBondContains {
                bond: MulticenterBondRef::from_ast(*bond, meta),
                atom: AtomRef::from_ast(*atom, meta),
            },
            MulticenterBondContainsAll { bond, atoms } => Self::MulticenterBondContainsAll {
                bond: MulticenterBondRef::from_ast(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::from_ast(a, meta)).collect(),
            },
            MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond: MulticenterBondRef::from_ast(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond: MulticenterBondRef::from_ast(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &()).unwrap()),
            },
            NoncovalentBondEnds { bond, atoms } => Self::NoncovalentBondEnds {
                bond: NoncovalentBondRef::from_ast(*bond, meta),
                atoms: [
                    AtomRef::from_ast(atoms[0], meta),
                    AtomRef::from_ast(atoms[1], meta),
                ],
            },
            NoncovalentBondContains { bond, atom } => Self::NoncovalentBondContains {
                bond: NoncovalentBondRef::from_ast(*bond, meta),
                atom: AtomRef::from_ast(*atom, meta),
            },
            NoncovalentBondEndsSatisfy { bond, predicates } => Self::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondRef::from_ast(*bond, meta),
                predicates: [
                    Box::new(AtomConstraintDsl::from_ast(&predicates[0], &()).unwrap()),
                    Box::new(AtomConstraintDsl::from_ast(&predicates[1], &()).unwrap()),
                ],
            },
        }
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<RelationalConstraint, ParseError> {
        use RelationalConstraintDsl::*;
        Ok(match self {
            DativeBondDonor { bond, atom } => RelationalConstraint::DativeBondDonor {
                bond: bond.into_ast(counts.dative_bond_count, meta)?,
                atom: atom.into_ast(counts.atom_count, meta)?,
            },
            DativeBondAcceptor { bond, atom } => RelationalConstraint::DativeBondAcceptor {
                bond: bond.into_ast(counts.dative_bond_count, meta)?,
                atom: atom.into_ast(counts.atom_count, meta)?,
            },
            DativeBondParallels { dative, parallel } => RelationalConstraint::DativeBondParallels {
                dative: dative.into_ast(counts.dative_bond_count, meta)?,
                parallel: parallel.into_ast(counts.bond_count, meta)?,
            },
            DativeBondDonorSatisfies { bond, predicate } => {
                RelationalConstraint::DativeBondDonorSatisfies {
                    bond: bond.into_ast(counts.dative_bond_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            DativeBondAcceptorSatisfies { bond, predicate } => {
                RelationalConstraint::DativeBondAcceptorSatisfies {
                    bond: bond.into_ast(counts.dative_bond_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            AromaticSystemAtoms { system, atoms } => RelationalConstraint::AromaticSystemAtoms {
                system: system.into_ast(counts.aromatic_system_count, meta)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.into_ast(counts.atom_count, meta))
                    .collect::<Result<_, _>>()?,
            },
            AromaticSystemContains { system, atom } => {
                RelationalConstraint::AromaticSystemContains {
                    system: system.into_ast(counts.aromatic_system_count, meta)?,
                    atom: atom.into_ast(counts.atom_count, meta)?,
                }
            }
            AromaticSystemContainsAll { system, atoms } => {
                RelationalConstraint::AromaticSystemContainsAll {
                    system: system.into_ast(counts.aromatic_system_count, meta)?,
                    atoms: atoms
                        .into_iter()
                        .map(|a| a.into_ast(counts.atom_count, meta))
                        .collect::<Result<_, _>>()?,
                }
            }
            AromaticSystemAllAtoms { system, predicate } => {
                RelationalConstraint::AromaticSystemAllAtoms {
                    system: system.into_ast(counts.aromatic_system_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            AromaticSystemAnyAtom { system, predicate } => {
                RelationalConstraint::AromaticSystemAnyAtom {
                    system: system.into_ast(counts.aromatic_system_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            MulticenterBondAtoms { bond, atoms } => RelationalConstraint::MulticenterBondAtoms {
                bond: bond.into_ast(counts.multicenter_bond_count, meta)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.into_ast(counts.atom_count, meta))
                    .collect::<Result<_, _>>()?,
            },
            MulticenterBondContains { bond, atom } => {
                RelationalConstraint::MulticenterBondContains {
                    bond: bond.into_ast(counts.multicenter_bond_count, meta)?,
                    atom: atom.into_ast(counts.atom_count, meta)?,
                }
            }
            MulticenterBondContainsAll { bond, atoms } => {
                RelationalConstraint::MulticenterBondContainsAll {
                    bond: bond.into_ast(counts.multicenter_bond_count, meta)?,
                    atoms: atoms
                        .into_iter()
                        .map(|a| a.into_ast(counts.atom_count, meta))
                        .collect::<Result<_, _>>()?,
                }
            }
            MulticenterBondAllAtoms { bond, predicate } => {
                RelationalConstraint::MulticenterBondAllAtoms {
                    bond: bond.into_ast(counts.multicenter_bond_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            MulticenterBondAnyAtom { bond, predicate } => {
                RelationalConstraint::MulticenterBondAnyAtom {
                    bond: bond.into_ast(counts.multicenter_bond_count, meta)?,
                    predicate: Box::new(predicate.into_ast(&()).unwrap()),
                }
            }
            NoncovalentBondEnds { bond, atoms } => {
                let [a, b] = atoms;
                RelationalConstraint::NoncovalentBondEnds {
                    bond: bond.into_ast(counts.noncovalent_bond_count, meta)?,
                    atoms: [
                        a.into_ast(counts.atom_count, meta)?,
                        b.into_ast(counts.atom_count, meta)?,
                    ],
                }
            }
            NoncovalentBondContains { bond, atom } => {
                RelationalConstraint::NoncovalentBondContains {
                    bond: bond.into_ast(counts.noncovalent_bond_count, meta)?,
                    atom: atom.into_ast(counts.atom_count, meta)?,
                }
            }
            NoncovalentBondEndsSatisfy { bond, predicates } => {
                let [a, b] = predicates;
                RelationalConstraint::NoncovalentBondEndsSatisfy {
                    bond: bond.into_ast(counts.noncovalent_bond_count, meta)?,
                    predicates: [
                        Box::new(a.into_ast(&()).unwrap()),
                        Box::new(b.into_ast(&()).unwrap()),
                    ],
                }
            }
        })
    }
}

// region: EDN payload codec

/// Render the `[owner_ref <target>]` 2-element pair shared by every variant.
fn render_pair(owner: Edn<'static>, target: Edn<'static>) -> Edn<'static> {
    Edn::Vector(vec![owner, target].into())
}

fn render_atom_refs(refs: &[AtomRef]) -> Edn<'static> {
    Edn::Vector(refs.iter().map(AtomRef::to_edn).collect::<Vec<_>>().into())
}

fn render_payload(dsl: &RelationalConstraintDsl) -> (&'static str, Edn<'static>) {
    use RelationalConstraintDsl::*;
    match dsl {
        DativeBondDonor { bond, atom } => (
            "dative-bond-donor",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        DativeBondAcceptor { bond, atom } => (
            "dative-bond-acceptor",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        DativeBondParallels { dative, parallel } => (
            "dative-bond-parallels",
            render_pair(dative.to_edn(), parallel.to_edn()),
        ),
        DativeBondDonorSatisfies { bond, predicate } => (
            "dative-bond-donor-satisfies",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        DativeBondAcceptorSatisfies { bond, predicate } => (
            "dative-bond-acceptor-satisfies",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        AromaticSystemAtoms { system, atoms } => (
            "aromatic-system-atoms",
            render_pair(system.to_edn(), render_atom_refs(atoms)),
        ),
        AromaticSystemContains { system, atom } => (
            "aromatic-system-contains",
            render_pair(system.to_edn(), atom.to_edn()),
        ),
        AromaticSystemContainsAll { system, atoms } => (
            "aromatic-system-contains-all",
            render_pair(system.to_edn(), render_atom_refs(atoms)),
        ),
        AromaticSystemAllAtoms { system, predicate } => (
            "aromatic-system-all-atoms",
            render_pair(system.to_edn(), predicate.to_edn()),
        ),
        AromaticSystemAnyAtom { system, predicate } => (
            "aromatic-system-any-atom",
            render_pair(system.to_edn(), predicate.to_edn()),
        ),
        MulticenterBondAtoms { bond, atoms } => (
            "multicenter-bond-atoms",
            render_pair(bond.to_edn(), render_atom_refs(atoms)),
        ),
        MulticenterBondContains { bond, atom } => (
            "multicenter-bond-contains",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        MulticenterBondContainsAll { bond, atoms } => (
            "multicenter-bond-contains-all",
            render_pair(bond.to_edn(), render_atom_refs(atoms)),
        ),
        MulticenterBondAllAtoms { bond, predicate } => (
            "multicenter-bond-all-atoms",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        MulticenterBondAnyAtom { bond, predicate } => (
            "multicenter-bond-any-atom",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        NoncovalentBondEnds { bond, atoms } => (
            "noncovalent-bond-ends",
            render_pair(
                bond.to_edn(),
                Edn::Vector(vec![atoms[0].to_edn(), atoms[1].to_edn()].into()),
            ),
        ),
        NoncovalentBondContains { bond, atom } => (
            "noncovalent-bond-contains",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        NoncovalentBondEndsSatisfy { bond, predicates } => (
            "noncovalent-bond-ends-satisfy",
            render_pair(
                bond.to_edn(),
                Edn::Vector(vec![predicates[0].to_edn(), predicates[1].to_edn()].into()),
            ),
        ),
    }
}

fn parse_pair<'a, 'de>(
    edn: &'a Edn<'de>,
    key: &str,
) -> Result<(&'a Edn<'de>, &'a Edn<'de>), DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "2-element vector",
            got: edn.kind(),
            path: vec![key.into()],
        });
    };
    if v.len() != 2 {
        return Err(DeError::Custom(format!(
            "{}: expected 2 elements, got {}",
            key,
            v.len()
        )));
    }
    Ok((&v[0], &v[1]))
}

fn parse_atom_refs<'de>(edn: &Edn<'de>, key: &str) -> Result<Vec<AtomRef>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of atom refs",
            got: edn.kind(),
            path: vec![key.into()],
        });
    };
    v.iter().map(AtomRef::from_edn).collect()
}

fn parse_atom_constraint_pair<'de>(
    edn: &Edn<'de>,
    key: &str,
) -> Result<[Box<AtomConstraintDsl>; 2], DeError> {
    let (a, b) = parse_pair(edn, key)?;
    Ok([
        Box::new(AtomConstraintDsl::from_edn(a)?),
        Box::new(AtomConstraintDsl::from_edn(b)?),
    ])
}

fn parse_payload<'de>(key: &str, edn: &Edn<'de>) -> Result<RelationalConstraintDsl, DeError> {
    use RelationalConstraintDsl::*;
    Ok(match key {
        "dative-bond-donor" => {
            let (bond, atom) = parse_pair(edn, key)?;
            DativeBondDonor {
                bond: DativeBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "dative-bond-acceptor" => {
            let (bond, atom) = parse_pair(edn, key)?;
            DativeBondAcceptor {
                bond: DativeBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "dative-bond-parallels" => {
            let (dative, parallel) = parse_pair(edn, key)?;
            DativeBondParallels {
                dative: DativeBondRef::from_edn(dative)?,
                parallel: BondRef::from_edn(parallel)?,
            }
        }
        "dative-bond-donor-satisfies" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            DativeBondDonorSatisfies {
                bond: DativeBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "dative-bond-acceptor-satisfies" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            DativeBondAcceptorSatisfies {
                bond: DativeBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "aromatic-system-atoms" => {
            let (system, atoms) = parse_pair(edn, key)?;
            AromaticSystemAtoms {
                system: AromaticSystemRef::from_edn(system)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "aromatic-system-contains" => {
            let (system, atom) = parse_pair(edn, key)?;
            AromaticSystemContains {
                system: AromaticSystemRef::from_edn(system)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "aromatic-system-contains-all" => {
            let (system, atoms) = parse_pair(edn, key)?;
            AromaticSystemContainsAll {
                system: AromaticSystemRef::from_edn(system)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "aromatic-system-all-atoms" => {
            let (system, predicate) = parse_pair(edn, key)?;
            AromaticSystemAllAtoms {
                system: AromaticSystemRef::from_edn(system)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "aromatic-system-any-atom" => {
            let (system, predicate) = parse_pair(edn, key)?;
            AromaticSystemAnyAtom {
                system: AromaticSystemRef::from_edn(system)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "multicenter-bond-atoms" => {
            let (bond, atoms) = parse_pair(edn, key)?;
            MulticenterBondAtoms {
                bond: MulticenterBondRef::from_edn(bond)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "multicenter-bond-contains" => {
            let (bond, atom) = parse_pair(edn, key)?;
            MulticenterBondContains {
                bond: MulticenterBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "multicenter-bond-contains-all" => {
            let (bond, atoms) = parse_pair(edn, key)?;
            MulticenterBondContainsAll {
                bond: MulticenterBondRef::from_edn(bond)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "multicenter-bond-all-atoms" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            MulticenterBondAllAtoms {
                bond: MulticenterBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "multicenter-bond-any-atom" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            MulticenterBondAnyAtom {
                bond: MulticenterBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "noncovalent-bond-ends" => {
            let (bond, ends) = parse_pair(edn, key)?;
            let (a, b) = parse_pair(ends, key)?;
            NoncovalentBondEnds {
                bond: NoncovalentBondRef::from_edn(bond)?,
                atoms: [AtomRef::from_edn(a)?, AtomRef::from_edn(b)?],
            }
        }
        "noncovalent-bond-contains" => {
            let (bond, atom) = parse_pair(edn, key)?;
            NoncovalentBondContains {
                bond: NoncovalentBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "noncovalent-bond-ends-satisfy" => {
            let (bond, predicates) = parse_pair(edn, key)?;
            NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondRef::from_edn(bond)?,
                predicates: parse_atom_constraint_pair(predicates, key)?,
            }
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["relational-constraint".into()],
            });
        }
    })
}

// endregion: EDN payload codec

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::constraint::AtomConstraint;
    use crate::ast::idx::{
        AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
    };
    use crate::ast::value::ValueAst;

    #[fixture]
    fn meta() -> Metadata {
        Metadata::default()
    }

    #[fixture]
    fn full_counts() -> EntityCounts {
        EntityCounts {
            atom_count: 10,
            bond_count: 10,
            dative_bond_count: 10,
            aromatic_system_count: 10,
            multicenter_bond_count: 10,
            noncovalent_bond_count: 10,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor(RelationalConstraint::DativeBondDonor { bond: DativeBondIdx(0), atom: AtomIdx(2) }, "{:dative-bond-donor [0 2]}")]
    #[case::dative_acceptor(RelationalConstraint::DativeBondAcceptor { bond: DativeBondIdx(1), atom: AtomIdx(3) }, "{:dative-bond-acceptor [1 3]}")]
    #[case::dative_parallels(RelationalConstraint::DativeBondParallels { dative: DativeBondIdx(0), parallel: BondIdx(2) }, "{:dative-bond-parallels [0 2]}")]
    #[case::dative_donor_satisfies(RelationalConstraint::DativeBondDonorSatisfies { bond: DativeBondIdx(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(3))) },
        "{:dative-bond-donor-satisfies [0 {:valence 3}]}")]
    #[case::dative_acceptor_satisfies(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondIdx(0),
        predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(2))) }, "{:dative-bond-acceptor-satisfies [0 {:degree 2}]}")]
    #[case::aromatic_system_atoms(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemIdx(0), atoms: vec![AtomIdx(0), AtomIdx(1)] },
        "{:aromatic-system-atoms [0 [0 1]]}")]
    #[case::aromatic_system_contains(RelationalConstraint::AromaticSystemContains { system: AromaticSystemIdx(0), atom: AtomIdx(2) }, "{:aromatic-system-contains [0 2]}")]
    #[case::aromatic_system_contains_all(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemIdx(0), atoms: vec![AtomIdx(2), AtomIdx(5)] },
        "{:aromatic-system-contains-all [0 [2 5]]}")]
    #[case::aromatic_system_all_atoms(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemIdx(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:aromatic-system-all-atoms [0 {:valence 4}]}")]
    #[case::aromatic_system_any_atom(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemIdx(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:aromatic-system-any-atom [0 {:degree 3}]}")]
    #[case::multicenter_atoms(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondIdx(0), atoms: vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)] },
        "{:multicenter-bond-atoms [0 [0 1 2]]}")]
    #[case::multicenter_contains(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondIdx(0), atom: AtomIdx(3) }, "{:multicenter-bond-contains [0 3]}")]
    #[case::multicenter_contains_all(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondIdx(0), atoms: vec![AtomIdx(0), AtomIdx(1)] },
        "{:multicenter-bond-contains-all [0 [0 1]]}")]
    #[case::multicenter_all_atoms(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondIdx(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:multicenter-bond-all-atoms [0 {:valence 4}]}")]
    #[case::multicenter_any_atom(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondIdx(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:multicenter-bond-any-atom [0 {:degree 3}]}")]
    #[case::noncovalent_ends(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondIdx(0), atoms: [AtomIdx(0), AtomIdx(3)] },
        "{:noncovalent-bond-ends [0 [0 3]]}")]
    #[case::noncovalent_contains(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondIdx(0), atom: AtomIdx(2) }, "{:noncovalent-bond-contains [0 2]}")]
    #[case::noncovalent_ends_satisfy(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondIdx(0),
        predicates: [Box::new(AtomConstraint::Valence(ValueAst::Lit(2))), Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))] },
        "{:noncovalent-bond-ends-satisfy [0 [{:valence 2} {:valence 3}]]}")]
    fn test_relational_constraint_dsl_roundtrip(
        meta: Metadata,
        full_counts: EntityCounts,
        #[case] input: RelationalConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = RelationalConstraintDsl::from_ast(&input, &meta);
        let edn = dsl.clone().to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = RelationalConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&full_counts, &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = RelationalConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_multi_key() {
        let edn = read_string("{:dative-bond-donor [0 1] :dative-bond-acceptor [0 2]}").unwrap();
        let err = RelationalConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_wrong_shape() {
        let err = RelationalConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_wrong_pair_length() {
        let edn = read_string("{:dative-bond-donor [0]}").unwrap();
        let err = RelationalConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_out_of_range_atom(meta: Metadata) {
        let counts = EntityCounts {
            atom_count: 2,
            bond_count: 0,
            dative_bond_count: 1,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
        };
        let dsl = RelationalConstraintDsl::DativeBondDonor {
            bond: DativeBondRef::Index(0),
            atom: AtomRef::Index(99),
        };
        let err = dsl.into_ast(&counts, &meta).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "99".into(),
            }
        );
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_out_of_range_bond(meta: Metadata) {
        let counts = EntityCounts {
            atom_count: 5,
            bond_count: 0,
            dative_bond_count: 5,
            aromatic_system_count: 0,
            multicenter_bond_count: 0,
            noncovalent_bond_count: 0,
        };
        let dsl = RelationalConstraintDsl::DativeBondDonor {
            bond: DativeBondRef::Index(99),
            atom: AtomRef::Index(0),
        };
        let err = dsl.into_ast(&counts, &meta).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "dative-bond",
                value: "99".into(),
            }
        );
    }
}
