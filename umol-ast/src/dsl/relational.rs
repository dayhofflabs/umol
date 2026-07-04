//! Surface DSL for cross-entity relational constraints.
//!
//! Mirrors [`crate::ast::constraint::RelationalConstraint`] with
//! `AtomRef`/`BondRef`/`DativeBondRef`/etc. in place of raw `AtomId` etc.
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
use super::error::ParseError;
use super::molecule::Metadata;
use super::namespace::Namespace;
use super::refs::{
    AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef, NoncovalentBondRef,
    StereoAtomRef, StereoBondRef,
};
use crate::ast::constraint::RelationalConstraint;
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL wrapper around [`RelationalConstraint`]. Structural parallel
/// to the AST enum — same 18 variants, with surface refs ([`AtomRef`],
/// [`DativeBondRef`], etc.) in place of raw `*Idx`. Each variant's EDN
/// form is a flat single-key map `{:<entity>-<role> [<owner_ref> <target>]}`.
/// See the AST enum for per-variant semantics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RelationalConstraintDsl {
    /// EDN: `{:dative-bond-donors [<dative_ref> [<atom_ref>+]]}`.
    DativeBondDonors {
        bond: DativeBondRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:dative-bond-donor [<dative_ref> <atom_ref>]}`.
    DativeBondDonor { bond: DativeBondRef, atom: AtomRef },
    /// EDN: `{:dative-bond-contains-all-donors [<dative_ref> [<atom_ref>+]]}`.
    DativeBondContainsAllDonors {
        bond: DativeBondRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:dative-bond-all-donors [<dative_ref> <atom-constraint>]}`.
    DativeBondAllDonors {
        bond: DativeBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:dative-bond-any-donor [<dative_ref> <atom-constraint>]}`.
    DativeBondAnyDonor {
        bond: DativeBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:dative-bond-acceptor [<dative_ref> <atom_ref>]}`.
    DativeBondAcceptor { bond: DativeBondRef, atom: AtomRef },
    /// EDN: `{:dative-bond-acceptor-satisfies [<dative_ref> <atom-constraint>]}`.
    DativeBondAcceptorSatisfies {
        bond: DativeBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:dative-bond-parallels [<dative_ref> <bond_ref>]}`.
    DativeBondParallels {
        dative: DativeBondRef,
        parallel: BondRef,
    },

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

    /// EDN: `{:stereo-atom-site [<stereo_atom_ref> <atom_ref>]}`.
    StereoAtomSite {
        stereo_atom: StereoAtomRef,
        atom: AtomRef,
    },
    /// EDN: `{:stereo-atom-contains [<stereo_atom_ref> <atom_ref>]}`.
    StereoAtomContains {
        stereo_atom: StereoAtomRef,
        atom: AtomRef,
    },
    /// EDN: `{:stereo-atom-ligands [<stereo_atom_ref> [<atom_ref>+]]}`.
    StereoAtomLigands {
        stereo_atom: StereoAtomRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:stereo-atom-all-ligands [<stereo_atom_ref> <atom-constraint>]}`.
    StereoAtomAllLigands {
        stereo_atom: StereoAtomRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:stereo-atom-any-ligand [<stereo_atom_ref> <atom-constraint>]}`.
    StereoAtomAnyLigand {
        stereo_atom: StereoAtomRef,
        predicate: Box<AtomConstraintDsl>,
    },

    /// EDN: `{:stereo-bond-site [<stereo_bond_ref> <bond_ref>]}`.
    StereoBondSite {
        stereo_bond: StereoBondRef,
        bond: BondRef,
    },
    /// EDN: `{:stereo-bond-contains [<stereo_bond_ref> <atom_ref>]}`.
    StereoBondContains {
        stereo_bond: StereoBondRef,
        atom: AtomRef,
    },
    /// EDN: `{:stereo-bond-ligands [<stereo_bond_ref> [<atom_ref>+]]}`.
    StereoBondLigands {
        stereo_bond: StereoBondRef,
        atoms: Vec<AtomRef>,
    },
    /// EDN: `{:stereo-bond-all-ligands [<stereo_bond_ref> <atom-constraint>]}`.
    StereoBondAllLigands {
        stereo_bond: StereoBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
    /// EDN: `{:stereo-bond-any-ligand [<stereo_bond_ref> <atom-constraint>]}`.
    StereoBondAnyLigand {
        stereo_bond: StereoBondRef,
        predicate: Box<AtomConstraintDsl>,
    },
}

/// Top-level EDN keywords for every relational variant. Matches the
/// flat-key convention — every relational form is its own `:<entity>-<role>`
/// keyword, rather than nested under the entity like narrow constraints.
pub(super) const RELATIONAL_KEYS: &[&str] = &[
    "dative-bond-donors",
    "dative-bond-donor",
    "dative-bond-contains-all-donors",
    "dative-bond-all-donors",
    "dative-bond-any-donor",
    "dative-bond-acceptor",
    "dative-bond-acceptor-satisfies",
    "dative-bond-parallels",
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
    "stereo-atom-site",
    "stereo-atom-contains",
    "stereo-atom-ligands",
    "stereo-atom-all-ligands",
    "stereo-atom-any-ligand",
    "stereo-bond-site",
    "stereo-bond-contains",
    "stereo-bond-ligands",
    "stereo-bond-all-ligands",
    "stereo-bond-any-ligand",
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
    pub(crate) fn from_ast<M: Metadata>(rel: &RelationalConstraint, meta: &M) -> Self {
        use RelationalConstraint::*;
        match rel {
            DativeBondDonors { bond, atoms } => Self::DativeBondDonors {
                bond: DativeBondRef::denote(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            DativeBondDonor { bond, atom } => Self::DativeBondDonor {
                bond: DativeBondRef::denote(*bond, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            DativeBondContainsAllDonors { bond, atoms } => Self::DativeBondContainsAllDonors {
                bond: DativeBondRef::denote(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            DativeBondAllDonors { bond, predicate } => Self::DativeBondAllDonors {
                bond: DativeBondRef::denote(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            DativeBondAnyDonor { bond, predicate } => Self::DativeBondAnyDonor {
                bond: DativeBondRef::denote(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            DativeBondAcceptor { bond, atom } => Self::DativeBondAcceptor {
                bond: DativeBondRef::denote(*bond, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            DativeBondAcceptorSatisfies { bond, predicate } => Self::DativeBondAcceptorSatisfies {
                bond: DativeBondRef::denote(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            DativeBondParallels { dative, parallel } => Self::DativeBondParallels {
                dative: DativeBondRef::denote(*dative, meta),
                parallel: BondRef::denote(*parallel, meta),
            },
            AromaticSystemAtoms { system, atoms } => Self::AromaticSystemAtoms {
                system: AromaticSystemRef::denote(*system, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            AromaticSystemContains { system, atom } => Self::AromaticSystemContains {
                system: AromaticSystemRef::denote(*system, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            AromaticSystemContainsAll { system, atoms } => Self::AromaticSystemContainsAll {
                system: AromaticSystemRef::denote(*system, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system: AromaticSystemRef::denote(*system, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system: AromaticSystemRef::denote(*system, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            MulticenterBondAtoms { bond, atoms } => Self::MulticenterBondAtoms {
                bond: MulticenterBondRef::denote(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            MulticenterBondContains { bond, atom } => Self::MulticenterBondContains {
                bond: MulticenterBondRef::denote(*bond, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            MulticenterBondContainsAll { bond, atoms } => Self::MulticenterBondContainsAll {
                bond: MulticenterBondRef::denote(*bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond: MulticenterBondRef::denote(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond: MulticenterBondRef::denote(*bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            NoncovalentBondEnds { bond, atoms } => Self::NoncovalentBondEnds {
                bond: NoncovalentBondRef::denote(*bond, meta),
                atoms: [
                    AtomRef::denote(atoms[0], meta),
                    AtomRef::denote(atoms[1], meta),
                ],
            },
            NoncovalentBondContains { bond, atom } => Self::NoncovalentBondContains {
                bond: NoncovalentBondRef::denote(*bond, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            NoncovalentBondEndsSatisfy { bond, predicates } => Self::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondRef::denote(*bond, meta),
                predicates: [
                    Box::new(AtomConstraintDsl::from_ast(&predicates[0], &())),
                    Box::new(AtomConstraintDsl::from_ast(&predicates[1], &())),
                ],
            },
            StereoAtomSite { stereo_atom, atom } => Self::StereoAtomSite {
                stereo_atom: StereoAtomRef::denote(*stereo_atom, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            StereoAtomContains { stereo_atom, atom } => Self::StereoAtomContains {
                stereo_atom: StereoAtomRef::denote(*stereo_atom, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            StereoAtomLigands { stereo_atom, atoms } => Self::StereoAtomLigands {
                stereo_atom: StereoAtomRef::denote(*stereo_atom, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAllLigands {
                stereo_atom: StereoAtomRef::denote(*stereo_atom, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand {
                stereo_atom: StereoAtomRef::denote(*stereo_atom, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            StereoBondSite { stereo_bond, bond } => Self::StereoBondSite {
                stereo_bond: StereoBondRef::denote(*stereo_bond, meta),
                bond: BondRef::denote(*bond, meta),
            },
            StereoBondContains { stereo_bond, atom } => Self::StereoBondContains {
                stereo_bond: StereoBondRef::denote(*stereo_bond, meta),
                atom: AtomRef::denote(*atom, meta),
            },
            StereoBondLigands { stereo_bond, atoms } => Self::StereoBondLigands {
                stereo_bond: StereoBondRef::denote(*stereo_bond, meta),
                atoms: atoms.iter().map(|&a| AtomRef::denote(a, meta)).collect(),
            },
            StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => Self::StereoBondAllLigands {
                stereo_bond: StereoBondRef::denote(*stereo_bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
            StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand {
                stereo_bond: StereoBondRef::denote(*stereo_bond, meta),
                predicate: Box::new(AtomConstraintDsl::from_ast(predicate, &())),
            },
        }
    }

    pub(crate) fn into_ast<N: Namespace>(
        self,
        namespace: &N,
    ) -> Result<RelationalConstraint, ParseError> {
        use RelationalConstraintDsl::*;
        Ok(match self {
            DativeBondDonors { bond, atoms } => RelationalConstraint::DativeBondDonors {
                bond: bond.resolve(namespace)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.resolve(namespace))
                    .collect::<Result<_, _>>()?,
            },
            DativeBondDonor { bond, atom } => RelationalConstraint::DativeBondDonor {
                bond: bond.resolve(namespace)?,
                atom: atom.resolve(namespace)?,
            },
            DativeBondContainsAllDonors { bond, atoms } => {
                RelationalConstraint::DativeBondContainsAllDonors {
                    bond: bond.resolve(namespace)?,
                    atoms: atoms
                        .into_iter()
                        .map(|a| a.resolve(namespace))
                        .collect::<Result<_, _>>()?,
                }
            }
            DativeBondAllDonors { bond, predicate } => RelationalConstraint::DativeBondAllDonors {
                bond: bond.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
            DativeBondAnyDonor { bond, predicate } => RelationalConstraint::DativeBondAnyDonor {
                bond: bond.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
            DativeBondAcceptor { bond, atom } => RelationalConstraint::DativeBondAcceptor {
                bond: bond.resolve(namespace)?,
                atom: atom.resolve(namespace)?,
            },
            DativeBondAcceptorSatisfies { bond, predicate } => {
                RelationalConstraint::DativeBondAcceptorSatisfies {
                    bond: bond.resolve(namespace)?,
                    predicate: Box::new(predicate.into_ast(&())),
                }
            }
            DativeBondParallels { dative, parallel } => RelationalConstraint::DativeBondParallels {
                dative: dative.resolve(namespace)?,
                parallel: parallel.resolve(namespace)?,
            },
            AromaticSystemAtoms { system, atoms } => RelationalConstraint::AromaticSystemAtoms {
                system: system.resolve(namespace)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.resolve(namespace))
                    .collect::<Result<_, _>>()?,
            },
            AromaticSystemContains { system, atom } => {
                RelationalConstraint::AromaticSystemContains {
                    system: system.resolve(namespace)?,
                    atom: atom.resolve(namespace)?,
                }
            }
            AromaticSystemContainsAll { system, atoms } => {
                RelationalConstraint::AromaticSystemContainsAll {
                    system: system.resolve(namespace)?,
                    atoms: atoms
                        .into_iter()
                        .map(|a| a.resolve(namespace))
                        .collect::<Result<_, _>>()?,
                }
            }
            AromaticSystemAllAtoms { system, predicate } => {
                RelationalConstraint::AromaticSystemAllAtoms {
                    system: system.resolve(namespace)?,
                    predicate: Box::new(predicate.into_ast(&())),
                }
            }
            AromaticSystemAnyAtom { system, predicate } => {
                RelationalConstraint::AromaticSystemAnyAtom {
                    system: system.resolve(namespace)?,
                    predicate: Box::new(predicate.into_ast(&())),
                }
            }
            MulticenterBondAtoms { bond, atoms } => RelationalConstraint::MulticenterBondAtoms {
                bond: bond.resolve(namespace)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.resolve(namespace))
                    .collect::<Result<_, _>>()?,
            },
            MulticenterBondContains { bond, atom } => {
                RelationalConstraint::MulticenterBondContains {
                    bond: bond.resolve(namespace)?,
                    atom: atom.resolve(namespace)?,
                }
            }
            MulticenterBondContainsAll { bond, atoms } => {
                RelationalConstraint::MulticenterBondContainsAll {
                    bond: bond.resolve(namespace)?,
                    atoms: atoms
                        .into_iter()
                        .map(|a| a.resolve(namespace))
                        .collect::<Result<_, _>>()?,
                }
            }
            MulticenterBondAllAtoms { bond, predicate } => {
                RelationalConstraint::MulticenterBondAllAtoms {
                    bond: bond.resolve(namespace)?,
                    predicate: Box::new(predicate.into_ast(&())),
                }
            }
            MulticenterBondAnyAtom { bond, predicate } => {
                RelationalConstraint::MulticenterBondAnyAtom {
                    bond: bond.resolve(namespace)?,
                    predicate: Box::new(predicate.into_ast(&())),
                }
            }
            NoncovalentBondEnds { bond, atoms } => {
                let [a, b] = atoms;
                RelationalConstraint::NoncovalentBondEnds {
                    bond: bond.resolve(namespace)?,
                    atoms: [a.resolve(namespace)?, b.resolve(namespace)?],
                }
            }
            NoncovalentBondContains { bond, atom } => {
                RelationalConstraint::NoncovalentBondContains {
                    bond: bond.resolve(namespace)?,
                    atom: atom.resolve(namespace)?,
                }
            }
            NoncovalentBondEndsSatisfy { bond, predicates } => {
                let [a, b] = predicates;
                RelationalConstraint::NoncovalentBondEndsSatisfy {
                    bond: bond.resolve(namespace)?,
                    predicates: [Box::new(a.into_ast(&())), Box::new(b.into_ast(&()))],
                }
            }
            StereoAtomSite { stereo_atom, atom } => RelationalConstraint::StereoAtomSite {
                stereo_atom: stereo_atom.resolve(namespace)?,
                atom: atom.resolve(namespace)?,
            },
            StereoAtomContains { stereo_atom, atom } => RelationalConstraint::StereoAtomContains {
                stereo_atom: stereo_atom.resolve(namespace)?,
                atom: atom.resolve(namespace)?,
            },
            StereoAtomLigands { stereo_atom, atoms } => RelationalConstraint::StereoAtomLigands {
                stereo_atom: stereo_atom.resolve(namespace)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.resolve(namespace))
                    .collect::<Result<_, _>>()?,
            },
            StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => RelationalConstraint::StereoAtomAllLigands {
                stereo_atom: stereo_atom.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
            StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => RelationalConstraint::StereoAtomAnyLigand {
                stereo_atom: stereo_atom.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
            StereoBondSite { stereo_bond, bond } => RelationalConstraint::StereoBondSite {
                stereo_bond: stereo_bond.resolve(namespace)?,
                bond: bond.resolve(namespace)?,
            },
            StereoBondContains { stereo_bond, atom } => RelationalConstraint::StereoBondContains {
                stereo_bond: stereo_bond.resolve(namespace)?,
                atom: atom.resolve(namespace)?,
            },
            StereoBondLigands { stereo_bond, atoms } => RelationalConstraint::StereoBondLigands {
                stereo_bond: stereo_bond.resolve(namespace)?,
                atoms: atoms
                    .into_iter()
                    .map(|a| a.resolve(namespace))
                    .collect::<Result<_, _>>()?,
            },
            StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => RelationalConstraint::StereoBondAllLigands {
                stereo_bond: stereo_bond.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
            StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => RelationalConstraint::StereoBondAnyLigand {
                stereo_bond: stereo_bond.resolve(namespace)?,
                predicate: Box::new(predicate.into_ast(&())),
            },
        })
    }
}

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
        DativeBondDonors { bond, atoms } => (
            "dative-bond-donors",
            render_pair(bond.to_edn(), render_atom_refs(atoms)),
        ),
        DativeBondDonor { bond, atom } => (
            "dative-bond-donor",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        DativeBondContainsAllDonors { bond, atoms } => (
            "dative-bond-contains-all-donors",
            render_pair(bond.to_edn(), render_atom_refs(atoms)),
        ),
        DativeBondAllDonors { bond, predicate } => (
            "dative-bond-all-donors",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        DativeBondAnyDonor { bond, predicate } => (
            "dative-bond-any-donor",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        DativeBondAcceptor { bond, atom } => (
            "dative-bond-acceptor",
            render_pair(bond.to_edn(), atom.to_edn()),
        ),
        DativeBondAcceptorSatisfies { bond, predicate } => (
            "dative-bond-acceptor-satisfies",
            render_pair(bond.to_edn(), predicate.to_edn()),
        ),
        DativeBondParallels { dative, parallel } => (
            "dative-bond-parallels",
            render_pair(dative.to_edn(), parallel.to_edn()),
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
        StereoAtomSite { stereo_atom, atom } => (
            "stereo-atom-site",
            render_pair(stereo_atom.to_edn(), atom.to_edn()),
        ),
        StereoAtomContains { stereo_atom, atom } => (
            "stereo-atom-contains",
            render_pair(stereo_atom.to_edn(), atom.to_edn()),
        ),
        StereoAtomLigands { stereo_atom, atoms } => (
            "stereo-atom-ligands",
            render_pair(stereo_atom.to_edn(), render_atom_refs(atoms)),
        ),
        StereoAtomAllLigands {
            stereo_atom,
            predicate,
        } => (
            "stereo-atom-all-ligands",
            render_pair(stereo_atom.to_edn(), predicate.to_edn()),
        ),
        StereoAtomAnyLigand {
            stereo_atom,
            predicate,
        } => (
            "stereo-atom-any-ligand",
            render_pair(stereo_atom.to_edn(), predicate.to_edn()),
        ),
        StereoBondSite { stereo_bond, bond } => (
            "stereo-bond-site",
            render_pair(stereo_bond.to_edn(), bond.to_edn()),
        ),
        StereoBondContains { stereo_bond, atom } => (
            "stereo-bond-contains",
            render_pair(stereo_bond.to_edn(), atom.to_edn()),
        ),
        StereoBondLigands { stereo_bond, atoms } => (
            "stereo-bond-ligands",
            render_pair(stereo_bond.to_edn(), render_atom_refs(atoms)),
        ),
        StereoBondAllLigands {
            stereo_bond,
            predicate,
        } => (
            "stereo-bond-all-ligands",
            render_pair(stereo_bond.to_edn(), predicate.to_edn()),
        ),
        StereoBondAnyLigand {
            stereo_bond,
            predicate,
        } => (
            "stereo-bond-any-ligand",
            render_pair(stereo_bond.to_edn(), predicate.to_edn()),
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
        "dative-bond-donors" => {
            let (bond, atoms) = parse_pair(edn, key)?;
            DativeBondDonors {
                bond: DativeBondRef::from_edn(bond)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "dative-bond-donor" => {
            let (bond, atom) = parse_pair(edn, key)?;
            DativeBondDonor {
                bond: DativeBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "dative-bond-contains-all-donors" => {
            let (bond, atoms) = parse_pair(edn, key)?;
            DativeBondContainsAllDonors {
                bond: DativeBondRef::from_edn(bond)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "dative-bond-all-donors" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            DativeBondAllDonors {
                bond: DativeBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "dative-bond-any-donor" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            DativeBondAnyDonor {
                bond: DativeBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "dative-bond-acceptor" => {
            let (bond, atom) = parse_pair(edn, key)?;
            DativeBondAcceptor {
                bond: DativeBondRef::from_edn(bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "dative-bond-acceptor-satisfies" => {
            let (bond, predicate) = parse_pair(edn, key)?;
            DativeBondAcceptorSatisfies {
                bond: DativeBondRef::from_edn(bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "dative-bond-parallels" => {
            let (dative, parallel) = parse_pair(edn, key)?;
            DativeBondParallels {
                dative: DativeBondRef::from_edn(dative)?,
                parallel: BondRef::from_edn(parallel)?,
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
        "stereo-atom-site" => {
            let (stereo_atom, atom) = parse_pair(edn, key)?;
            StereoAtomSite {
                stereo_atom: StereoAtomRef::from_edn(stereo_atom)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "stereo-atom-contains" => {
            let (stereo_atom, atom) = parse_pair(edn, key)?;
            StereoAtomContains {
                stereo_atom: StereoAtomRef::from_edn(stereo_atom)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "stereo-atom-ligands" => {
            let (stereo_atom, atoms) = parse_pair(edn, key)?;
            StereoAtomLigands {
                stereo_atom: StereoAtomRef::from_edn(stereo_atom)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "stereo-atom-all-ligands" => {
            let (stereo_atom, predicate) = parse_pair(edn, key)?;
            StereoAtomAllLigands {
                stereo_atom: StereoAtomRef::from_edn(stereo_atom)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "stereo-atom-any-ligand" => {
            let (stereo_atom, predicate) = parse_pair(edn, key)?;
            StereoAtomAnyLigand {
                stereo_atom: StereoAtomRef::from_edn(stereo_atom)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "stereo-bond-site" => {
            let (stereo_bond, bond) = parse_pair(edn, key)?;
            StereoBondSite {
                stereo_bond: StereoBondRef::from_edn(stereo_bond)?,
                bond: BondRef::from_edn(bond)?,
            }
        }
        "stereo-bond-contains" => {
            let (stereo_bond, atom) = parse_pair(edn, key)?;
            StereoBondContains {
                stereo_bond: StereoBondRef::from_edn(stereo_bond)?,
                atom: AtomRef::from_edn(atom)?,
            }
        }
        "stereo-bond-ligands" => {
            let (stereo_bond, atoms) = parse_pair(edn, key)?;
            StereoBondLigands {
                stereo_bond: StereoBondRef::from_edn(stereo_bond)?,
                atoms: parse_atom_refs(atoms, key)?,
            }
        }
        "stereo-bond-all-ligands" => {
            let (stereo_bond, predicate) = parse_pair(edn, key)?;
            StereoBondAllLigands {
                stereo_bond: StereoBondRef::from_edn(stereo_bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
            }
        }
        "stereo-bond-any-ligand" => {
            let (stereo_bond, predicate) = parse_pair(edn, key)?;
            StereoBondAnyLigand {
                stereo_bond: StereoBondRef::from_edn(stereo_bond)?,
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?),
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::super::molecule::MoleculeMetadata;
    use super::super::namespace::MoleculeNamespace;
    use super::*;
    use crate::ast::constraint::AtomConstraint;
    use crate::ast::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
    use crate::ast::value::ValueAst;

    #[fixture]
    fn meta() -> MoleculeMetadata {
        MoleculeMetadata::default()
    }

    /// A namespace with ten entities of each kind, so index refs up to 9 resolve.
    #[fixture]
    fn full_namespace() -> MoleculeNamespace {
        let mut namespace = MoleculeNamespace::default();
        for _ in 0..10 {
            namespace.register_atom(None).unwrap();
        }
        for _ in 0..10 {
            namespace.register_bond(None, AtomId(0), AtomId(1)).unwrap();
            namespace
                .register_dative_bond(None, &[AtomId(0)], AtomId(1))
                .unwrap();
            namespace
                .register_aromatic_system(None, &[AtomId(0)])
                .unwrap();
            namespace
                .register_multicenter_bond(None, &[AtomId(0)])
                .unwrap();
            namespace
                .register_noncovalent_bond(None, AtomId(0), AtomId(1))
                .unwrap();
            namespace
                .register_stereo_atom(None, AtomId(0), &[])
                .unwrap();
            namespace
                .register_stereo_bond(None, BondId(0), &[])
                .unwrap();
        }
        namespace
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(2) }, "{:dative-bond-donor [0 2]}")]
    #[case::dative_acceptor(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(1), atom: AtomId(3) }, "{:dative-bond-acceptor [1 3]}")]
    #[case::dative_parallels(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(2) }, "{:dative-bond-parallels [0 2]}")]
    #[case::dative_donors(RelationalConstraint::DativeBondDonors { bond: DativeBondId(0), atoms: vec![AtomId(1), AtomId(2)] },
        "{:dative-bond-donors [0 [1 2]]}")]
    #[case::dative_contains_all_donors(RelationalConstraint::DativeBondContainsAllDonors { bond: DativeBondId(0), atoms: vec![AtomId(1), AtomId(2)] },
        "{:dative-bond-contains-all-donors [0 [1 2]]}")]
    #[case::dative_all_donors(RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(3))) },
        "{:dative-bond-all-donors [0 {:valence 3}]}")]
    #[case::dative_any_donor(RelationalConstraint::DativeBondAnyDonor { bond: DativeBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(3))) },
        "{:dative-bond-any-donor [0 {:valence 3}]}")]
    #[case::dative_acceptor_satisfies(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondId(0),
        predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(2))) }, "{:dative-bond-acceptor-satisfies [0 {:degree 2}]}")]
    #[case::aromatic_system_atoms(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0), atoms: vec![AtomId(0), AtomId(1)] },
        "{:aromatic-system-atoms [0 [0 1]]}")]
    #[case::aromatic_system_contains(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(2) }, "{:aromatic-system-contains [0 2]}")]
    #[case::aromatic_system_contains_all(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0), atoms: vec![AtomId(2), AtomId(5)] },
        "{:aromatic-system-contains-all [0 [2 5]]}")]
    #[case::aromatic_system_all_atoms(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:aromatic-system-all-atoms [0 {:valence 4}]}")]
    #[case::aromatic_system_any_atom(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemId(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:aromatic-system-any-atom [0 {:degree 3}]}")]
    #[case::multicenter_atoms(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0), atoms: vec![AtomId(0), AtomId(1), AtomId(2)] },
        "{:multicenter-bond-atoms [0 [0 1 2]]}")]
    #[case::multicenter_contains(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(3) }, "{:multicenter-bond-contains [0 3]}")]
    #[case::multicenter_contains_all(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0), atoms: vec![AtomId(0), AtomId(1)] },
        "{:multicenter-bond-contains-all [0 [0 1]]}")]
    #[case::multicenter_all_atoms(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:multicenter-bond-all-atoms [0 {:valence 4}]}")]
    #[case::multicenter_any_atom(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:multicenter-bond-any-atom [0 {:degree 3}]}")]
    #[case::noncovalent_ends(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(3)] },
        "{:noncovalent-bond-ends [0 [0 3]]}")]
    #[case::noncovalent_contains(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(2) }, "{:noncovalent-bond-contains [0 2]}")]
    #[case::noncovalent_ends_satisfy(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
        predicates: [Box::new(AtomConstraint::Valence(ValueAst::Lit(2))), Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))] },
        "{:noncovalent-bond-ends-satisfy [0 [{:valence 2} {:valence 3}]]}")]
    #[case::stereo_atom_site(RelationalConstraint::StereoAtomSite { stereo_atom: StereoAtomId(0), atom: AtomId(2) },
        "{:stereo-atom-site [0 2]}")]
    #[case::stereo_atom_contains(RelationalConstraint::StereoAtomContains { stereo_atom: StereoAtomId(0), atom: AtomId(3) },
        "{:stereo-atom-contains [0 3]}")]
    #[case::stereo_atom_ligands(RelationalConstraint::StereoAtomLigands { stereo_atom: StereoAtomId(0), atoms: vec![AtomId(0), AtomId(1)] },
        "{:stereo-atom-ligands [0 [0 1]]}")]
    #[case::stereo_atom_all_ligands(RelationalConstraint::StereoAtomAllLigands { stereo_atom: StereoAtomId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:stereo-atom-all-ligands [0 {:valence 4}]}")]
    #[case::stereo_atom_any_ligand(RelationalConstraint::StereoAtomAnyLigand { stereo_atom: StereoAtomId(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:stereo-atom-any-ligand [0 {:degree 3}]}")]
    #[case::stereo_bond_site(RelationalConstraint::StereoBondSite { stereo_bond: StereoBondId(0), bond: BondId(2) },
        "{:stereo-bond-site [0 2]}")]
    #[case::stereo_bond_contains(RelationalConstraint::StereoBondContains { stereo_bond: StereoBondId(0), atom: AtomId(3) },
        "{:stereo-bond-contains [0 3]}")]
    #[case::stereo_bond_ligands(RelationalConstraint::StereoBondLigands { stereo_bond: StereoBondId(0), atoms: vec![AtomId(0), AtomId(1)] },
        "{:stereo-bond-ligands [0 [0 1]]}")]
    #[case::stereo_bond_all_ligands(RelationalConstraint::StereoBondAllLigands { stereo_bond: StereoBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) },
        "{:stereo-bond-all-ligands [0 {:valence 4}]}")]
    #[case::stereo_bond_any_ligand(RelationalConstraint::StereoBondAnyLigand { stereo_bond: StereoBondId(0), predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) },
        "{:stereo-bond-any-ligand [0 {:degree 3}]}")]
    fn test_relational_constraint_dsl_roundtrip(
        meta: MoleculeMetadata,
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[case] input: RelationalConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = RelationalConstraintDsl::from_ast(&input, &meta);
        let edn = dsl.clone().to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = RelationalConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&namespace).unwrap();
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
    fn test_relational_constraint_dsl_rejects_out_of_range_atom() {
        let mut namespace = MoleculeNamespace::default();
        namespace.register_atom(None).unwrap();
        namespace.register_atom(None).unwrap();
        namespace
            .register_dative_bond(None, &[AtomId(0)], AtomId(1))
            .unwrap();
        let dsl = RelationalConstraintDsl::DativeBondDonor {
            bond: DativeBondRef::Index(0),
            atom: AtomRef::Index(99),
        };
        let err = dsl.into_ast(&namespace).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "99".into(),
            }
        );
    }

    #[rstest]
    fn test_relational_constraint_dsl_rejects_out_of_range_bond() {
        let mut namespace = MoleculeNamespace::default();
        for _ in 0..5 {
            namespace.register_atom(None).unwrap();
        }
        for _ in 0..5 {
            namespace
                .register_dative_bond(None, &[AtomId(0)], AtomId(1))
                .unwrap();
        }
        let dsl = RelationalConstraintDsl::DativeBondDonor {
            bond: DativeBondRef::Index(99),
            atom: AtomRef::Index(0),
        };
        let err = dsl.into_ast(&namespace).unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "dative-bond",
                value: "99".into(),
            }
        );
    }
}
