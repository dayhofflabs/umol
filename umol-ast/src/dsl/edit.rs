//! Surface encoding for handles in standalone edit documents.

#![expect(
    dead_code,
    reason = "edit-family codecs are assembled by the aggregate EditsDsl root"
)]

use umol_edn::{DeError, Edn, EdnMap, EdnMapHelper, FromEdn, ToEdn};

use super::aromatic::{AromaticSystemDsl, AromaticSystemUpdateDsl};
use super::atom::{AtomDsl, AtomUpdateDsl};
use super::bond::{BondDsl, BondUpdateDsl};
use super::config::MoleculeDefaults;
use super::constraint::ConstraintDsl;
use super::dative::{DativeBondDsl, DativeBondUpdateDsl};
use super::edn_utils::{parse_single_key_map, single_key_map};
use super::metadata::MoleculeMetadata;
use super::multicenter::{MulticenterBondDsl, MulticenterBondUpdateDsl};
use super::namespace::Namespace;
use super::noncovalent::{NoncovalentBondDsl, NoncovalentBondUpdateDsl};
use crate::ast::aromatic::AromaticSystemUpdate;
use crate::ast::atom::AtomUpdate;
use crate::ast::bond::BondUpdate;
use crate::ast::constraint::{
    AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, Constraint,
    DativeBondConstraintAst, MulticenterBondConstraintAst, NoncovalentBondConstraintAst,
};
use crate::ast::dative::DativeBondUpdate;
use crate::ast::edit::{
    AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange,
    BondHandle, DativeBondFieldChange, DativeBondHandle, Edit, Edits, MulticenterBondFieldChange,
    MulticenterBondHandle, NoncovalentBondFieldChange, NoncovalentBondHandle, StereoAtomHandle,
    StereoBondHandle,
};
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;
use crate::ast::multicenter::MulticenterBondUpdate;
use crate::ast::noncovalent::NoncovalentBondUpdate;
use crate::ast::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};
use crate::ast::traits::{FromAst, IntoAst, Lattice};

/// Surface form shared by every typed handle in a standalone edit document.
///
/// A bare integer identifies an entity in the initial host. `{:new n}` identifies the `n`th
/// same-kind entity created earlier in the edit sequence.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum EditHandleDsl {
    Id(u32),
    New(usize),
}

impl<'de> FromEdn<'de> for EditHandleDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Int(_) => u32::from_edn(edn).map(Self::Id),
            Edn::Map(map) => {
                let mut helper = EdnMapHelper::new(map);
                let index = helper.required("new")?;
                helper.finalize()?;
                if map.len() != 1 {
                    return Err(DeError::Custom(
                        "edit handle map keys must be keywords".to_string(),
                    ));
                }
                Ok(Self::New(index))
            }
            other => Err(DeError::TypeMismatch {
                expected: "edit handle (non-negative integer or {:new n} map)",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for EditHandleDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Id(index) => index.to_edn(),
            Self::New(index) => single_key_map("new", index.to_edn()),
        }
    }
}

macro_rules! impl_typed_handle_conversion {
    ($handle:ident, $id:ident) => {
        impl From<$handle> for EditHandleDsl {
            fn from(handle: $handle) -> Self {
                match handle {
                    $handle::Id(id) => Self::Id(id.0),
                    $handle::New(index) => Self::New(index),
                }
            }
        }

        impl From<EditHandleDsl> for $handle {
            fn from(handle: EditHandleDsl) -> Self {
                match handle {
                    EditHandleDsl::Id(index) => Self::Id($id(index)),
                    EditHandleDsl::New(index) => Self::New(index),
                }
            }
        }

        impl<'de> FromEdn<'de> for $handle {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                EditHandleDsl::from_edn(edn).map(Self::from)
            }
        }

        impl ToEdn for $handle {
            fn to_edn(&self) -> Edn<'static> {
                EditHandleDsl::from(self.clone()).to_edn()
            }
        }
    };
}

impl_typed_handle_conversion!(AtomHandle, AtomId);
impl_typed_handle_conversion!(BondHandle, BondId);
impl_typed_handle_conversion!(DativeBondHandle, DativeBondId);
impl_typed_handle_conversion!(AromaticSystemHandle, AromaticSystemId);
impl_typed_handle_conversion!(MulticenterBondHandle, MulticenterBondId);
impl_typed_handle_conversion!(NoncovalentBondHandle, NoncovalentBondId);
impl_typed_handle_conversion!(StereoAtomHandle, StereoAtomId);
impl_typed_handle_conversion!(StereoBondHandle, StereoBondId);

#[derive(Clone, Debug, PartialEq, Eq)]
enum EditInput {
    AtomAdd(AtomDsl),
    AtomRemove(AtomHandle),
    AtomModify {
        id: AtomHandle,
        expect: AtomUpdate,
        update: AtomUpdate,
    },
    BondAdd {
        atoms: [AtomHandle; 2],
        ast: BondDsl,
    },
    BondRemove(BondHandle),
    BondModify {
        id: BondHandle,
        expect: BondUpdate,
        update: BondUpdate,
    },
    DativeBondAdd {
        donors: Vec<AtomHandle>,
        acceptor: AtomHandle,
        ast: DativeBondDsl,
    },
    DativeBondsRemove(Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondDsl)>),
    DativeBondModify {
        id: DativeBondHandle,
        expect: DativeBondUpdate,
        update: DativeBondUpdate,
    },
    AromaticSystemAdd {
        atoms: Vec<AtomHandle>,
        ast: AromaticSystemDsl,
    },
    AromaticSystemsRemove(Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemDsl)>),
    AromaticSystemModify {
        id: AromaticSystemHandle,
        expect: AromaticSystemUpdate,
        update: AromaticSystemUpdate,
    },
    MulticenterBondAdd {
        atoms: Vec<AtomHandle>,
        ast: MulticenterBondDsl,
    },
    MulticenterBondsRemove(Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondDsl)>),
    MulticenterBondModify {
        id: MulticenterBondHandle,
        expect: MulticenterBondUpdate,
        update: MulticenterBondUpdate,
    },
    NoncovalentBondAdd {
        atoms: [AtomHandle; 2],
        ast: NoncovalentBondDsl,
    },
    NoncovalentBondsRemove(Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondDsl)>),
    NoncovalentBondModify {
        id: NoncovalentBondHandle,
        expect: NoncovalentBondUpdate,
        update: NoncovalentBondUpdate,
    },
    TopologyRemove {
        atoms: Vec<AtomHandle>,
        bonds: Vec<BondHandle>,
    },
    ConstraintAdd(Constraint),
    ConstraintRemove(Constraint),
}

impl<'de> FromEdn<'de> for EditInput {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let (entity, body) = parse_single_key_map(edn, "edit")?;
        match entity {
            "atom" => parse_atom_edit(body),
            "bond" => parse_bond_edit(body),
            "dative-bond" => parse_dative_bond_edit(body),
            "dative-bonds" => parse_dative_bonds_edit(body),
            "aromatic-system" => parse_aromatic_system_edit(body),
            "aromatic-systems" => parse_aromatic_systems_edit(body),
            "multicenter-bond" => parse_multicenter_bond_edit(body),
            "multicenter-bonds" => parse_multicenter_bonds_edit(body),
            "noncovalent-bond" => parse_noncovalent_bond_edit(body),
            "noncovalent-bonds" => parse_noncovalent_bonds_edit(body),
            "topology" => parse_topology_edit(body),
            "constraint" => parse_constraint_edit(body),
            other => Err(DeError::Custom(format!("unknown edit :{other}"))),
        }
    }
}

impl ToEdn for EditInput {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::AtomAdd(ast) => edit_map("atom", "add", ast.to_edn()),
            Self::AtomRemove(id) => edit_map("atom", "remove", id.to_edn()),
            Self::AtomModify { id, expect, update } => edit_map(
                "atom",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    AtomUpdateDsl(expect.clone()).to_edn(),
                    AtomUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::BondAdd {
                atoms: [first, second],
                ast,
            } => edit_map(
                "bond",
                "add",
                Edn::Vector(vec![first.to_edn(), second.to_edn(), ast.to_edn()].into()),
            ),
            Self::BondRemove(id) => edit_map("bond", "remove", id.to_edn()),
            Self::BondModify { id, expect, update } => edit_map(
                "bond",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    BondUpdateDsl(expect.clone()).to_edn(),
                    BondUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::DativeBondAdd {
                donors,
                acceptor,
                ast,
            } => edit_map(
                "dative-bond",
                "add",
                dative_entry_edn(None, donors, acceptor, ast.to_edn()),
            ),
            Self::DativeBondsRemove(removes) => edit_map(
                "dative-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, ast)| {
                            let (acceptor, donors) = atoms
                                .split_last()
                                .expect("dative edit always has an acceptor");
                            dative_entry_edn(Some(id.to_edn()), donors, acceptor, ast.to_edn())
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::DativeBondModify { id, expect, update } => edit_map(
                "dative-bond",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    DativeBondUpdateDsl(expect.clone()).to_edn(),
                    DativeBondUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::AromaticSystemAdd { atoms, ast } => edit_map(
                "aromatic-system",
                "add",
                relation_entry_edn(None, atoms, ast.to_edn()),
            ),
            Self::AromaticSystemsRemove(removes) => edit_map(
                "aromatic-systems",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, ast)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, ast.to_edn())
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::AromaticSystemModify { id, expect, update } => edit_map(
                "aromatic-system",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    AromaticSystemUpdateDsl(expect.clone()).to_edn(),
                    AromaticSystemUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::MulticenterBondAdd { atoms, ast } => edit_map(
                "multicenter-bond",
                "add",
                relation_entry_edn(None, atoms, ast.to_edn()),
            ),
            Self::MulticenterBondsRemove(removes) => edit_map(
                "multicenter-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, ast)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, ast.to_edn())
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::MulticenterBondModify { id, expect, update } => edit_map(
                "multicenter-bond",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    MulticenterBondUpdateDsl(expect.clone()).to_edn(),
                    MulticenterBondUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::NoncovalentBondAdd { atoms, ast } => edit_map(
                "noncovalent-bond",
                "add",
                relation_entry_edn(None, atoms, ast.to_edn()),
            ),
            Self::NoncovalentBondsRemove(removes) => edit_map(
                "noncovalent-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, ast)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, ast.to_edn())
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::NoncovalentBondModify { id, expect, update } => edit_map(
                "noncovalent-bond",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    NoncovalentBondUpdateDsl(expect.clone()).to_edn(),
                    NoncovalentBondUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::TopologyRemove { atoms, bonds } => {
                let mut removal = EdnMap::with_capacity(2);
                removal.insert(
                    Edn::keyword("atoms"),
                    Edn::Vector(atoms.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into()),
                );
                removal.insert(
                    Edn::keyword("bonds"),
                    Edn::Vector(bonds.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into()),
                );
                edit_map("topology", "remove", Edn::Map(removal))
            }
            Self::ConstraintAdd(constraint) => {
                edit_map("constraint", "add", render_constraint(constraint).to_edn())
            }
            Self::ConstraintRemove(constraint) => edit_map(
                "constraint",
                "remove",
                render_constraint(constraint).to_edn(),
            ),
        }
    }
}

impl EditInput {
    fn append_to(self, edits: &mut Edits, defaults: &MoleculeDefaults) -> Result<(), DeError> {
        match self {
            Self::AtomAdd(ast) => {
                edits.add_atom(ast.into_ast(&defaults.atom));
            }
            Self::AtomRemove(id) => edits.remove_atom(id),
            Self::AtomModify { id, expect, update } => {
                append_atom_modify(edits, id, expect, update)?;
            }
            Self::BondAdd {
                atoms: [first, second],
                ast,
            } => {
                edits.add_bond(first, second, ast.into_ast(&defaults.bond));
            }
            Self::BondRemove(id) => edits.remove_bond(id),
            Self::BondModify { id, expect, update } => {
                append_bond_modify(edits, id, expect, update)?;
            }
            Self::DativeBondAdd {
                mut donors,
                acceptor,
                ast,
            } => {
                donors.push(acceptor);
                edits.add_dative_bond(donors, ast.into_ast(&defaults.dative_bond));
            }
            Self::DativeBondsRemove(removes) => edits.remove_dative_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, ast)| (id, atoms, ast.into_ast(&defaults.dative_bond)))
                    .collect(),
            ),
            Self::DativeBondModify { id, expect, update } => {
                append_dative_bond_modify(edits, id, expect, update)?;
            }
            Self::AromaticSystemAdd { atoms, ast } => {
                edits.add_aromatic_system(atoms, ast.into_ast(&defaults.aromatic_system));
            }
            Self::AromaticSystemsRemove(removes) => edits.remove_aromatic_systems(
                removes
                    .into_iter()
                    .map(|(id, atoms, ast)| (id, atoms, ast.into_ast(&defaults.aromatic_system)))
                    .collect(),
            ),
            Self::AromaticSystemModify { id, expect, update } => {
                append_aromatic_system_modify(edits, id, expect, update)?;
            }
            Self::MulticenterBondAdd { atoms, ast } => {
                edits.add_multicenter_bond(atoms, ast.into_ast(&defaults.multicenter_bond));
            }
            Self::MulticenterBondsRemove(removes) => edits.remove_multicenter_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, ast)| (id, atoms, ast.into_ast(&defaults.multicenter_bond)))
                    .collect(),
            ),
            Self::MulticenterBondModify { id, expect, update } => {
                append_multicenter_bond_modify(edits, id, expect, update)?;
            }
            Self::NoncovalentBondAdd { atoms, ast } => {
                edits.add_noncovalent_bond(atoms, ast.into_ast(&defaults.noncovalent_bond));
            }
            Self::NoncovalentBondsRemove(removes) => edits.remove_noncovalent_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, ast)| (id, atoms, ast.into_ast(&defaults.noncovalent_bond)))
                    .collect(),
            ),
            Self::NoncovalentBondModify { id, expect, update } => {
                append_noncovalent_bond_modify(edits, id, expect, update)?;
            }
            Self::TopologyRemove { atoms, bonds } => edits.remove_topology(atoms, bonds),
            Self::ConstraintAdd(constraint) => edits.add_molecule_constraint(constraint),
            Self::ConstraintRemove(constraint) => edits.remove_molecule_constraint(constraint),
        }
        Ok(())
    }

    fn from_edit(edit: &Edit, defaults: &MoleculeDefaults) -> Result<Option<Vec<Self>>, DeError> {
        let inputs = match edit {
            Edit::AddAtoms { atoms } => atoms
                .iter()
                .map(|ast| Self::AtomAdd(AtomDsl::from_ast(ast, &defaults.atom)))
                .collect(),
            Edit::AddBonds { bonds } => bonds
                .iter()
                .map(|bond| Self::BondAdd {
                    atoms: bond.endpoints.clone(),
                    ast: BondDsl::from_ast(&bond.ast, &defaults.bond),
                })
                .collect(),
            Edit::RemoveTopology { atoms, bonds } if atoms.len() == 1 && bonds.is_empty() => {
                vec![Self::AtomRemove(atoms[0].clone())]
            }
            Edit::RemoveTopology { atoms, bonds } if atoms.is_empty() && bonds.len() == 1 => {
                vec![Self::BondRemove(bonds[0].clone())]
            }
            Edit::RemoveTopology { atoms, bonds } => vec![Self::TopologyRemove {
                atoms: atoms.clone(),
                bonds: bonds.clone(),
            }],
            Edit::ModifyAtomField { id, change } => {
                let (expect, update) = atom_field_updates(change);
                vec![Self::AtomModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyBondField { id, change } => {
                let (expect, update) = bond_field_updates(change);
                vec![Self::BondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyAtomConstraint { id, old, new } => {
                let (expect, update) = atom_constraint_updates(old, new)?;
                vec![Self::AtomModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyBondConstraint { id, old, new } => {
                let (expect, update) = bond_constraint_updates(old, new)?;
                vec![Self::BondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddDativeBond { atoms, ast } => {
                let (acceptor, donors) = atoms.split_last().ok_or_else(|| {
                    DeError::Custom("dative-bond addition has no acceptor".to_string())
                })?;
                vec![Self::DativeBondAdd {
                    donors: donors.to_vec(),
                    acceptor: acceptor.clone(),
                    ast: DativeBondDsl::from_ast(ast, &defaults.dative_bond),
                }]
            }
            Edit::RemoveDativeBonds { removes } => vec![Self::DativeBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        if atoms.is_empty() {
                            return Err(DeError::Custom(
                                "dative-bond removal has no acceptor".to_string(),
                            ));
                        }
                        Ok((
                            id.clone(),
                            atoms.clone(),
                            DativeBondDsl::from_ast(ast, &defaults.dative_bond),
                        ))
                    })
                    .collect::<Result<_, _>>()?,
            )],
            Edit::ModifyDativeBondField { id, change } => {
                let (expect, update) = dative_bond_field_updates(change);
                vec![Self::DativeBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyDativeBondConstraint { id, old, new } => {
                let (expect, update) = dative_bond_constraint_updates(old, new)?;
                vec![Self::DativeBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddAromaticSystem { atoms, ast } => vec![Self::AromaticSystemAdd {
                atoms: atoms.clone(),
                ast: AromaticSystemDsl::from_ast(ast, &defaults.aromatic_system),
            }],
            Edit::RemoveAromaticSystems { removes } => vec![Self::AromaticSystemsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            AromaticSystemDsl::from_ast(ast, &defaults.aromatic_system),
                        )
                    })
                    .collect(),
            )],
            Edit::ModifyAromaticSystemField { id, change } => {
                let (expect, update) = aromatic_system_field_updates(change);
                vec![Self::AromaticSystemModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyAromaticSystemConstraint { id, old, new } => {
                let (expect, update) = aromatic_system_constraint_updates(old, new)?;
                vec![Self::AromaticSystemModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddMulticenterBond { atoms, ast } => vec![Self::MulticenterBondAdd {
                atoms: atoms.clone(),
                ast: MulticenterBondDsl::from_ast(ast, &defaults.multicenter_bond),
            }],
            Edit::RemoveMulticenterBonds { removes } => vec![Self::MulticenterBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            MulticenterBondDsl::from_ast(ast, &defaults.multicenter_bond),
                        )
                    })
                    .collect(),
            )],
            Edit::ModifyMulticenterBondField { id, change } => {
                let (expect, update) = multicenter_bond_field_updates(change);
                vec![Self::MulticenterBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyMulticenterBondConstraint { id, old, new } => {
                let (expect, update) = multicenter_bond_constraint_updates(old, new)?;
                vec![Self::MulticenterBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddNoncovalentBond { atoms, ast } => vec![Self::NoncovalentBondAdd {
                atoms: atoms.clone(),
                ast: NoncovalentBondDsl::from_ast(ast, &defaults.noncovalent_bond),
            }],
            Edit::RemoveNoncovalentBonds { removes } => vec![Self::NoncovalentBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            NoncovalentBondDsl::from_ast(ast, &defaults.noncovalent_bond),
                        )
                    })
                    .collect(),
            )],
            Edit::ModifyNoncovalentBondField { id, change } => {
                let (expect, update) = noncovalent_bond_field_updates(change);
                vec![Self::NoncovalentBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyNoncovalentBondConstraint { id, old, new } => {
                let (expect, update) = noncovalent_bond_constraint_updates(old, new)?;
                vec![Self::NoncovalentBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddMoleculeConstraint { constraint } => {
                vec![Self::ConstraintAdd(constraint.clone())]
            }
            Edit::RemoveMoleculeConstraint { constraint } => {
                vec![Self::ConstraintRemove(constraint.clone())]
            }
            _ => return Ok(None),
        };
        Ok(Some(inputs))
    }
}

fn parse_atom_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "atom edit")?;
    match op {
        "add" => Ok(EditInput::AtomAdd(AtomDsl::from_edn(payload)?)),
        "remove" => Ok(EditInput::AtomRemove(AtomHandle::from_edn(payload)?)),
        "modify" => {
            let (id, expect, update) = parse_atom_checked_update(payload)?;
            validate_atom_update_pair(&expect, &update)?;
            Ok(EditInput::AtomModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!("unknown atom edit op :{other}"))),
    }
}

fn parse_bond_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "bond edit")?;
    match op {
        "add" => {
            let Edn::Vector(parts) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "bond :add [first second dsl]",
                    got: payload.kind(),
                    path: vec!["bond edit".to_string()],
                });
            };
            if parts.len() != 3 {
                return Err(DeError::Custom(format!(
                    "bond :add expects [first second dsl], got {} elements",
                    parts.len()
                )));
            }
            Ok(EditInput::BondAdd {
                atoms: [
                    AtomHandle::from_edn(&parts[0])?,
                    AtomHandle::from_edn(&parts[1])?,
                ],
                ast: BondDsl::from_edn(&parts[2])?,
            })
        }
        "remove" => Ok(EditInput::BondRemove(BondHandle::from_edn(payload)?)),
        "modify" => {
            let (id, expect, update) = parse_bond_checked_update(payload)?;
            validate_bond_update_pair(&expect, &update)?;
            Ok(EditInput::BondModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!("unknown bond edit op :{other}"))),
    }
}

fn parse_dative_bond_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "dative-bond edit")?;
    match op {
        "add" => {
            let (donors, acceptor, ast) = parse_dative_bond_addition(payload)?;
            Ok(EditInput::DativeBondAdd {
                donors,
                acceptor,
                ast,
            })
        }
        "modify" => {
            let (id, expect, update) = parse_dative_bond_checked_update(payload)?;
            validate_dative_bond_update_pair(&expect, &update)?;
            Ok(EditInput::DativeBondModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown dative-bond edit op :{other}"
        ))),
    }
}

fn parse_dative_bonds_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "dative-bonds edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown dative-bonds edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of dative-bond removals",
            got: payload.kind(),
            path: vec!["dative-bonds edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_dative_bond_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::DativeBondsRemove(removes))
}

fn parse_aromatic_system_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "aromatic-system edit")?;
    match op {
        "add" => {
            let (atoms, ast) = parse_aromatic_system_addition(payload)?;
            Ok(EditInput::AromaticSystemAdd { atoms, ast })
        }
        "modify" => {
            let (id, expect, update) = parse_aromatic_system_checked_update(payload)?;
            validate_aromatic_system_update_pair(&expect, &update)?;
            Ok(EditInput::AromaticSystemModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown aromatic-system edit op :{other}"
        ))),
    }
}

fn parse_aromatic_systems_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "aromatic-systems edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown aromatic-systems edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of aromatic-system removals",
            got: payload.kind(),
            path: vec!["aromatic-systems edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_aromatic_system_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::AromaticSystemsRemove(removes))
}

fn parse_multicenter_bond_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "multicenter-bond edit")?;
    match op {
        "add" => {
            let (atoms, ast) = parse_multicenter_bond_addition(payload)?;
            Ok(EditInput::MulticenterBondAdd { atoms, ast })
        }
        "modify" => {
            let (id, expect, update) = parse_multicenter_bond_checked_update(payload)?;
            validate_multicenter_bond_update_pair(&expect, &update)?;
            Ok(EditInput::MulticenterBondModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown multicenter-bond edit op :{other}"
        ))),
    }
}

fn parse_multicenter_bonds_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "multicenter-bonds edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown multicenter-bonds edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of multicenter-bond removals",
            got: payload.kind(),
            path: vec!["multicenter-bonds edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_multicenter_bond_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::MulticenterBondsRemove(removes))
}

fn parse_noncovalent_bond_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "noncovalent-bond edit")?;
    match op {
        "add" => {
            let (atoms, ast) = parse_noncovalent_bond_addition(payload)?;
            Ok(EditInput::NoncovalentBondAdd { atoms, ast })
        }
        "modify" => {
            let (id, expect, update) = parse_noncovalent_bond_checked_update(payload)?;
            validate_noncovalent_bond_update_pair(&expect, &update)?;
            Ok(EditInput::NoncovalentBondModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown noncovalent-bond edit op :{other}"
        ))),
    }
}

fn parse_noncovalent_bonds_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "noncovalent-bonds edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown noncovalent-bonds edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of noncovalent-bond removals",
            got: payload.kind(),
            path: vec!["noncovalent-bonds edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_noncovalent_bond_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::NoncovalentBondsRemove(removes))
}

fn parse_dative_bond_addition(
    edn: &Edn<'_>,
) -> Result<(Vec<AtomHandle>, AtomHandle, DativeBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "dative-bond addition map",
            got: edn.kind(),
            path: vec!["dative-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let donors = helper.required("donors")?;
    let acceptor = helper.required("acceptor")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((donors, acceptor, ast))
}

fn parse_dative_bond_removal(
    edn: &Edn<'_>,
) -> Result<(DativeBondHandle, Vec<AtomHandle>, DativeBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "dative-bond removal map",
            got: edn.kind(),
            path: vec!["dative-bonds edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let mut donors: Vec<AtomHandle> = helper.required("donors")?;
    let acceptor = helper.required("acceptor")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    donors.push(acceptor);
    Ok((id, donors, ast))
}

fn parse_aromatic_system_addition(
    edn: &Edn<'_>,
) -> Result<(Vec<AtomHandle>, AromaticSystemDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "aromatic-system addition map",
            got: edn.kind(),
            path: vec!["aromatic-system edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((atoms, ast))
}

fn parse_aromatic_system_removal(
    edn: &Edn<'_>,
) -> Result<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "aromatic-system removal map",
            got: edn.kind(),
            path: vec!["aromatic-systems edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((id, atoms, ast))
}

fn parse_multicenter_bond_addition(
    edn: &Edn<'_>,
) -> Result<(Vec<AtomHandle>, MulticenterBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "multicenter-bond addition map",
            got: edn.kind(),
            path: vec!["multicenter-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((atoms, ast))
}

fn parse_multicenter_bond_removal(
    edn: &Edn<'_>,
) -> Result<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "multicenter-bond removal map",
            got: edn.kind(),
            path: vec!["multicenter-bonds edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((id, atoms, ast))
}

fn parse_noncovalent_bond_addition(
    edn: &Edn<'_>,
) -> Result<([AtomHandle; 2], NoncovalentBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "noncovalent-bond addition map",
            got: edn.kind(),
            path: vec!["noncovalent-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((atoms, ast))
}

fn parse_noncovalent_bond_removal(
    edn: &Edn<'_>,
) -> Result<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondDsl), DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "noncovalent-bond removal map",
            got: edn.kind(),
            path: vec!["noncovalent-bonds edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let atoms = helper.required("atoms")?;
    let ast = helper.required("type")?;
    helper.finalize()?;
    Ok((id, atoms, ast))
}

fn parse_topology_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "topology edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!("unknown topology edit op :{op}")));
    }
    let Edn::Map(map) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "topology :remove map",
            got: payload.kind(),
            path: vec!["topology edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let atoms = helper.required("atoms")?;
    let bonds = helper.required("bonds")?;
    helper.finalize()?;
    Ok(EditInput::TopologyRemove { atoms, bonds })
}

fn parse_constraint_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "constraint edit")?;
    let constraint = ConstraintDsl::from_edn(payload)?
        .into_ast(&EditConstraintNamespace)
        .map_err(|error| DeError::subgrammar("standalone-edit-constraint", error))?;
    match op {
        "add" => Ok(EditInput::ConstraintAdd(constraint)),
        "remove" => Ok(EditInput::ConstraintRemove(constraint)),
        other => Err(DeError::Custom(format!(
            "unknown constraint edit op :{other}"
        ))),
    }
}

fn parse_atom_checked_update(
    edn: &Edn<'_>,
) -> Result<(AtomHandle, AtomUpdate, AtomUpdate), DeError> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "atom :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["atom edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "atom :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "atom :modify changes map",
            got: parts[1].kind(),
            path: vec!["atom edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: AtomUpdateDsl = helper.required("expect")?;
    let update: AtomUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((AtomHandle::from_edn(&parts[0])?, expect.0, update.0))
}

fn parse_bond_checked_update(
    edn: &Edn<'_>,
) -> Result<(BondHandle, BondUpdate, BondUpdate), DeError> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "bond :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["bond edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "bond :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "bond :modify changes map",
            got: parts[1].kind(),
            path: vec!["bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: BondUpdateDsl = helper.required("expect")?;
    let update: BondUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((BondHandle::from_edn(&parts[0])?, expect.0, update.0))
}

fn parse_dative_bond_checked_update(
    edn: &Edn<'_>,
) -> Result<(DativeBondHandle, DativeBondUpdate, DativeBondUpdate), DeError> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "dative-bond :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["dative-bond edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "dative-bond :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "dative-bond :modify changes map",
            got: parts[1].kind(),
            path: vec!["dative-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: DativeBondUpdateDsl = helper.required("expect")?;
    let update: DativeBondUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((DativeBondHandle::from_edn(&parts[0])?, expect.0, update.0))
}

fn parse_aromatic_system_checked_update(
    edn: &Edn<'_>,
) -> Result<
    (
        AromaticSystemHandle,
        AromaticSystemUpdate,
        AromaticSystemUpdate,
    ),
    DeError,
> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "aromatic-system :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["aromatic-system edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "aromatic-system :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "aromatic-system :modify changes map",
            got: parts[1].kind(),
            path: vec!["aromatic-system edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: AromaticSystemUpdateDsl = helper.required("expect")?;
    let update: AromaticSystemUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((
        AromaticSystemHandle::from_edn(&parts[0])?,
        expect.0,
        update.0,
    ))
}

fn parse_multicenter_bond_checked_update(
    edn: &Edn<'_>,
) -> Result<
    (
        MulticenterBondHandle,
        MulticenterBondUpdate,
        MulticenterBondUpdate,
    ),
    DeError,
> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "multicenter-bond :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["multicenter-bond edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "multicenter-bond :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "multicenter-bond :modify changes map",
            got: parts[1].kind(),
            path: vec!["multicenter-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: MulticenterBondUpdateDsl = helper.required("expect")?;
    let update: MulticenterBondUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((
        MulticenterBondHandle::from_edn(&parts[0])?,
        expect.0,
        update.0,
    ))
}

fn parse_noncovalent_bond_checked_update(
    edn: &Edn<'_>,
) -> Result<
    (
        NoncovalentBondHandle,
        NoncovalentBondUpdate,
        NoncovalentBondUpdate,
    ),
    DeError,
> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "noncovalent-bond :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["noncovalent-bond edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "noncovalent-bond :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "noncovalent-bond :modify changes map",
            got: parts[1].kind(),
            path: vec!["noncovalent-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: NoncovalentBondUpdateDsl = helper.required("expect")?;
    let update: NoncovalentBondUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((
        NoncovalentBondHandle::from_edn(&parts[0])?,
        expect.0,
        update.0,
    ))
}

fn validate_atom_update_pair(expect: &AtomUpdate, update: &AtomUpdate) -> Result<(), DeError> {
    let fields_match = expect.element.is_some() == update.element.is_some()
        && expect.isotope_mass.is_some() == update.isotope_mass.is_some()
        && expect.charge.is_some() == update.charge.is_some()
        && expect.implicit_hydrogens.is_some() == update.implicit_hydrogens.is_some()
        && expect.lone_pairs.is_some() == update.lone_pairs.is_some()
        && expect.unpaired_electrons.count.is_some() == update.unpaired_electrons.count.is_some()
        && expect.unpaired_electrons.multiplicity.is_some()
            == update.unpaired_electrons.multiplicity.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(AtomConstraintAst::key)
        .eq(update.constraints.iter().map(AtomConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "atom :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    if expect.unpaired_electrons.count.is_some() != expect.unpaired_electrons.multiplicity.is_some()
    {
        return Err(DeError::Custom(
            "atom :modify unpaired-electron changes require both #u and #s".to_string(),
        ));
    }
    Ok(())
}

fn validate_bond_update_pair(expect: &BondUpdate, update: &BondUpdate) -> Result<(), DeError> {
    let fields_match = expect.order.is_some() == update.order.is_some()
        && expect.charge.is_some() == update.charge.is_some()
        && expect.unpaired_electrons.count.is_some() == update.unpaired_electrons.count.is_some()
        && expect.unpaired_electrons.multiplicity.is_some()
            == update.unpaired_electrons.multiplicity.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(BondConstraintAst::key)
        .eq(update.constraints.iter().map(BondConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    if expect.unpaired_electrons.count.is_some() != expect.unpaired_electrons.multiplicity.is_some()
    {
        return Err(DeError::Custom(
            "bond :modify unpaired-electron changes require both #u and #s".to_string(),
        ));
    }
    Ok(())
}

fn validate_dative_bond_update_pair(
    expect: &DativeBondUpdate,
    update: &DativeBondUpdate,
) -> Result<(), DeError> {
    let fields_match = expect.order.is_some() == update.order.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(DativeBondConstraintAst::key)
        .eq(update.constraints.iter().map(DativeBondConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "dative-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    Ok(())
}

fn validate_aromatic_system_update_pair(
    expect: &AromaticSystemUpdate,
    update: &AromaticSystemUpdate,
) -> Result<(), DeError> {
    let fields_match = expect.electrons.is_some() == update.electrons.is_some()
        && expect.charge.is_some() == update.charge.is_some()
        && expect.unpaired_electrons.count.is_some() == update.unpaired_electrons.count.is_some()
        && expect.unpaired_electrons.multiplicity.is_some()
            == update.unpaired_electrons.multiplicity.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(AromaticSystemConstraintAst::key)
        .eq(update
            .constraints
            .iter()
            .map(AromaticSystemConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "aromatic-system :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    validate_complete_unpaired_electrons(
        &expect.unpaired_electrons,
        "aromatic-system :modify unpaired-electron changes require both #u and #s",
    )
}

fn validate_multicenter_bond_update_pair(
    expect: &MulticenterBondUpdate,
    update: &MulticenterBondUpdate,
) -> Result<(), DeError> {
    let fields_match = expect.electrons.is_some() == update.electrons.is_some()
        && expect.charge.is_some() == update.charge.is_some()
        && expect.unpaired_electrons.count.is_some() == update.unpaired_electrons.count.is_some()
        && expect.unpaired_electrons.multiplicity.is_some()
            == update.unpaired_electrons.multiplicity.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(MulticenterBondConstraintAst::key)
        .eq(update
            .constraints
            .iter()
            .map(MulticenterBondConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "multicenter-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    validate_complete_unpaired_electrons(
        &expect.unpaired_electrons,
        "multicenter-bond :modify unpaired-electron changes require both #u and #s",
    )
}

fn validate_noncovalent_bond_update_pair(
    expect: &NoncovalentBondUpdate,
    update: &NoncovalentBondUpdate,
) -> Result<(), DeError> {
    let fields_match = expect.kind.is_some() == update.kind.is_some();
    let constraints_match = expect
        .constraints
        .iter()
        .map(NoncovalentBondConstraintAst::key)
        .eq(update
            .constraints
            .iter()
            .map(NoncovalentBondConstraintAst::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "noncovalent-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    Ok(())
}

fn validate_complete_unpaired_electrons(
    update: &UnpairedElectronsUpdate,
    message: &'static str,
) -> Result<(), DeError> {
    if update.count.is_some() != update.multiplicity.is_some() {
        return Err(DeError::Custom(message.to_string()));
    }
    Ok(())
}

fn append_atom_modify(
    edits: &mut Edits,
    id: AtomHandle,
    expect: AtomUpdate,
    update: AtomUpdate,
) -> Result<(), DeError> {
    validate_atom_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.element, update.element) {
        edits.push(Edit::ModifyAtomField {
            id: id.clone(),
            change: AtomFieldChange::Element { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.isotope_mass, update.isotope_mass) {
        edits.push(Edit::ModifyAtomField {
            id: id.clone(),
            change: AtomFieldChange::IsotopeMass { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.charge, update.charge) {
        edits.push(Edit::ModifyAtomField {
            id: id.clone(),
            change: AtomFieldChange::Charge { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.implicit_hydrogens, update.implicit_hydrogens) {
        edits.push(Edit::ModifyAtomField {
            id: id.clone(),
            change: AtomFieldChange::ImplicitHydrogens { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.lone_pairs, update.lone_pairs) {
        edits.push(Edit::ModifyAtomField {
            id: id.clone(),
            change: AtomFieldChange::LonePairs { old, new },
        });
    }
    append_atom_unpaired_electrons(
        edits,
        id.clone(),
        expect.unpaired_electrons,
        update.unpaired_electrons,
    );
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyAtomConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_bond_modify(
    edits: &mut Edits,
    id: BondHandle,
    expect: BondUpdate,
    update: BondUpdate,
) -> Result<(), DeError> {
    validate_bond_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.order, update.order) {
        edits.push(Edit::ModifyBondField {
            id: id.clone(),
            change: BondFieldChange::Order { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.charge, update.charge) {
        edits.push(Edit::ModifyBondField {
            id: id.clone(),
            change: BondFieldChange::Charge { old, new },
        });
    }
    append_bond_unpaired_electrons(
        edits,
        id.clone(),
        expect.unpaired_electrons,
        update.unpaired_electrons,
    );
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyBondConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_dative_bond_modify(
    edits: &mut Edits,
    id: DativeBondHandle,
    expect: DativeBondUpdate,
    update: DativeBondUpdate,
) -> Result<(), DeError> {
    validate_dative_bond_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.order, update.order) {
        edits.push(Edit::ModifyDativeBondField {
            id: id.clone(),
            change: DativeBondFieldChange::Order { old, new },
        });
    }
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyDativeBondConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_aromatic_system_modify(
    edits: &mut Edits,
    id: AromaticSystemHandle,
    expect: AromaticSystemUpdate,
    update: AromaticSystemUpdate,
) -> Result<(), DeError> {
    validate_aromatic_system_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.electrons, update.electrons) {
        edits.push(Edit::ModifyAromaticSystemField {
            id: id.clone(),
            change: AromaticSystemFieldChange::Electrons { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.charge, update.charge) {
        edits.push(Edit::ModifyAromaticSystemField {
            id: id.clone(),
            change: AromaticSystemFieldChange::Charge { old, new },
        });
    }
    append_aromatic_system_unpaired_electrons(
        edits,
        id.clone(),
        expect.unpaired_electrons,
        update.unpaired_electrons,
    );
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyAromaticSystemConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_multicenter_bond_modify(
    edits: &mut Edits,
    id: MulticenterBondHandle,
    expect: MulticenterBondUpdate,
    update: MulticenterBondUpdate,
) -> Result<(), DeError> {
    validate_multicenter_bond_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.electrons, update.electrons) {
        edits.push(Edit::ModifyMulticenterBondField {
            id: id.clone(),
            change: MulticenterBondFieldChange::Electrons { old, new },
        });
    }
    if let (Some(old), Some(new)) = (expect.charge, update.charge) {
        edits.push(Edit::ModifyMulticenterBondField {
            id: id.clone(),
            change: MulticenterBondFieldChange::Charge { old, new },
        });
    }
    append_multicenter_bond_unpaired_electrons(
        edits,
        id.clone(),
        expect.unpaired_electrons,
        update.unpaired_electrons,
    );
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyMulticenterBondConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_noncovalent_bond_modify(
    edits: &mut Edits,
    id: NoncovalentBondHandle,
    expect: NoncovalentBondUpdate,
    update: NoncovalentBondUpdate,
) -> Result<(), DeError> {
    validate_noncovalent_bond_update_pair(&expect, &update)?;
    if let (Some(old), Some(new)) = (expect.kind, update.kind) {
        edits.push(Edit::ModifyNoncovalentBondField {
            id: id.clone(),
            change: NoncovalentBondFieldChange::Kind { old, new },
        });
    }
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyNoncovalentBondConstraint {
            id: id.clone(),
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_atom_unpaired_electrons(
    edits: &mut Edits,
    id: AtomHandle,
    expect: UnpairedElectronsUpdate,
    update: UnpairedElectronsUpdate,
) {
    if let (Some(old_count), Some(old_multiplicity), Some(new_count), Some(new_multiplicity)) = (
        expect.count,
        expect.multiplicity,
        update.count,
        update.multiplicity,
    ) {
        edits.push(Edit::ModifyAtomField {
            id,
            change: AtomFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsAst {
                    count: new_count,
                    multiplicity: new_multiplicity,
                },
            },
        });
    }
}

fn append_bond_unpaired_electrons(
    edits: &mut Edits,
    id: BondHandle,
    expect: UnpairedElectronsUpdate,
    update: UnpairedElectronsUpdate,
) {
    if let (Some(old_count), Some(old_multiplicity), Some(new_count), Some(new_multiplicity)) = (
        expect.count,
        expect.multiplicity,
        update.count,
        update.multiplicity,
    ) {
        edits.push(Edit::ModifyBondField {
            id,
            change: BondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsAst {
                    count: new_count,
                    multiplicity: new_multiplicity,
                },
            },
        });
    }
}

fn append_aromatic_system_unpaired_electrons(
    edits: &mut Edits,
    id: AromaticSystemHandle,
    expect: UnpairedElectronsUpdate,
    update: UnpairedElectronsUpdate,
) {
    if let (Some(old_count), Some(old_multiplicity), Some(new_count), Some(new_multiplicity)) = (
        expect.count,
        expect.multiplicity,
        update.count,
        update.multiplicity,
    ) {
        edits.push(Edit::ModifyAromaticSystemField {
            id,
            change: AromaticSystemFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsAst {
                    count: new_count,
                    multiplicity: new_multiplicity,
                },
            },
        });
    }
}

fn append_multicenter_bond_unpaired_electrons(
    edits: &mut Edits,
    id: MulticenterBondHandle,
    expect: UnpairedElectronsUpdate,
    update: UnpairedElectronsUpdate,
) {
    if let (Some(old_count), Some(old_multiplicity), Some(new_count), Some(new_multiplicity)) = (
        expect.count,
        expect.multiplicity,
        update.count,
        update.multiplicity,
    ) {
        edits.push(Edit::ModifyMulticenterBondField {
            id,
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsAst {
                    count: new_count,
                    multiplicity: new_multiplicity,
                },
            },
        });
    }
}

fn atom_field_updates(change: &AtomFieldChange) -> (AtomUpdate, AtomUpdate) {
    let mut expect = AtomUpdate::default();
    let mut update = AtomUpdate::default();
    match change {
        AtomFieldChange::Element { old, new } => {
            expect.element = Some(old.clone());
            update.element = Some(new.clone());
        }
        AtomFieldChange::IsotopeMass { old, new } => {
            expect.isotope_mass = Some(old.clone());
            update.isotope_mass = Some(new.clone());
        }
        AtomFieldChange::Charge { old, new } => {
            expect.charge = Some(old.clone());
            update.charge = Some(new.clone());
        }
        AtomFieldChange::ImplicitHydrogens { old, new } => {
            expect.implicit_hydrogens = Some(old.clone());
            update.implicit_hydrogens = Some(new.clone());
        }
        AtomFieldChange::LonePairs { old, new } => {
            expect.lone_pairs = Some(old.clone());
            update.lone_pairs = Some(new.clone());
        }
        AtomFieldChange::UnpairedElectrons { old, new } => {
            expect.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(old.count.clone()),
                multiplicity: Some(old.multiplicity.clone()),
            };
            update.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(new.count.clone()),
                multiplicity: Some(new.multiplicity.clone()),
            };
        }
    }
    (expect, update)
}

fn bond_field_updates(change: &BondFieldChange) -> (BondUpdate, BondUpdate) {
    let mut expect = BondUpdate::default();
    let mut update = BondUpdate::default();
    match change {
        BondFieldChange::Order { old, new } => {
            expect.order = Some(old.clone());
            update.order = Some(new.clone());
        }
        BondFieldChange::Charge { old, new } => {
            expect.charge = Some(old.clone());
            update.charge = Some(new.clone());
        }
        BondFieldChange::UnpairedElectrons { old, new } => {
            expect.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(old.count.clone()),
                multiplicity: Some(old.multiplicity.clone()),
            };
            update.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(new.count.clone()),
                multiplicity: Some(new.multiplicity.clone()),
            };
        }
    }
    (expect, update)
}

fn dative_bond_field_updates(
    change: &DativeBondFieldChange,
) -> (DativeBondUpdate, DativeBondUpdate) {
    let mut expect = DativeBondUpdate::default();
    let mut update = DativeBondUpdate::default();
    match change {
        DativeBondFieldChange::Order { old, new } => {
            expect.order = Some(old.clone());
            update.order = Some(new.clone());
        }
    }
    (expect, update)
}

fn aromatic_system_field_updates(
    change: &AromaticSystemFieldChange,
) -> (AromaticSystemUpdate, AromaticSystemUpdate) {
    let mut expect = AromaticSystemUpdate::default();
    let mut update = AromaticSystemUpdate::default();
    match change {
        AromaticSystemFieldChange::Electrons { old, new } => {
            expect.electrons = Some(old.clone());
            update.electrons = Some(new.clone());
        }
        AromaticSystemFieldChange::Charge { old, new } => {
            expect.charge = Some(old.clone());
            update.charge = Some(new.clone());
        }
        AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
            expect.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(old.count.clone()),
                multiplicity: Some(old.multiplicity.clone()),
            };
            update.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(new.count.clone()),
                multiplicity: Some(new.multiplicity.clone()),
            };
        }
    }
    (expect, update)
}

fn multicenter_bond_field_updates(
    change: &MulticenterBondFieldChange,
) -> (MulticenterBondUpdate, MulticenterBondUpdate) {
    let mut expect = MulticenterBondUpdate::default();
    let mut update = MulticenterBondUpdate::default();
    match change {
        MulticenterBondFieldChange::Electrons { old, new } => {
            expect.electrons = Some(old.clone());
            update.electrons = Some(new.clone());
        }
        MulticenterBondFieldChange::Charge { old, new } => {
            expect.charge = Some(old.clone());
            update.charge = Some(new.clone());
        }
        MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
            expect.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(old.count.clone()),
                multiplicity: Some(old.multiplicity.clone()),
            };
            update.unpaired_electrons = UnpairedElectronsUpdate {
                count: Some(new.count.clone()),
                multiplicity: Some(new.multiplicity.clone()),
            };
        }
    }
    (expect, update)
}

fn noncovalent_bond_field_updates(
    change: &NoncovalentBondFieldChange,
) -> (NoncovalentBondUpdate, NoncovalentBondUpdate) {
    let mut expect = NoncovalentBondUpdate::default();
    let mut update = NoncovalentBondUpdate::default();
    match change {
        NoncovalentBondFieldChange::Kind { old, new } => {
            expect.kind = Some(old.clone());
            update.kind = Some(new.clone());
        }
    }
    (expect, update)
}

fn atom_constraint_updates(
    old: &Option<AtomConstraintAst>,
    new: &Option<AtomConstraintAst>,
) -> Result<(AtomUpdate, AtomUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "atom constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = AtomUpdate::default();
    let mut update = AtomUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn bond_constraint_updates(
    old: &Option<BondConstraintAst>,
    new: &Option<BondConstraintAst>,
) -> Result<(BondUpdate, BondUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "bond constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = BondUpdate::default();
    let mut update = BondUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn dative_bond_constraint_updates(
    old: &Option<DativeBondConstraintAst>,
    new: &Option<DativeBondConstraintAst>,
) -> Result<(DativeBondUpdate, DativeBondUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "dative-bond constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = DativeBondUpdate::default();
    let mut update = DativeBondUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn aromatic_system_constraint_updates(
    old: &Option<AromaticSystemConstraintAst>,
    new: &Option<AromaticSystemConstraintAst>,
) -> Result<(AromaticSystemUpdate, AromaticSystemUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "aromatic-system constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = AromaticSystemUpdate::default();
    let mut update = AromaticSystemUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn multicenter_bond_constraint_updates(
    old: &Option<MulticenterBondConstraintAst>,
    new: &Option<MulticenterBondConstraintAst>,
) -> Result<(MulticenterBondUpdate, MulticenterBondUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "multicenter-bond constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = MulticenterBondUpdate::default();
    let mut update = MulticenterBondUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn noncovalent_bond_constraint_updates(
    old: &Option<NoncovalentBondConstraintAst>,
    new: &Option<NoncovalentBondConstraintAst>,
) -> Result<(NoncovalentBondUpdate, NoncovalentBondUpdate), DeError> {
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "noncovalent-bond constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = NoncovalentBondUpdate::default();
    let mut update = NoncovalentBondUpdate::default();
    match (old, new) {
        (Some(old), Some(new)) => {
            expect.constraints.set(old.clone());
            update.constraints.set(new.clone());
        }
        (Some(old), None) => {
            expect.constraints.set(old.clone());
            update.constraints.set(old.as_undetermined());
        }
        (None, Some(new)) => {
            expect.constraints.set(new.as_undetermined());
            update.constraints.set(new.clone());
        }
        (None, None) => unreachable!(),
    }
    Ok((expect, update))
}

fn checked_update_edn(
    handle: Edn<'static>,
    expect: Edn<'static>,
    update: Edn<'static>,
) -> Edn<'static> {
    let mut changes = EdnMap::with_capacity(2);
    changes.insert(Edn::keyword("expect"), expect);
    changes.insert(Edn::keyword("update"), update);
    Edn::Vector(vec![handle, Edn::Map(changes)].into())
}

fn dative_entry_edn(
    id: Option<Edn<'static>>,
    donors: &[AtomHandle],
    acceptor: &AtomHandle,
    type_edn: Edn<'static>,
) -> Edn<'static> {
    let mut entry = EdnMap::with_capacity(4);
    if let Some(id) = id {
        entry.insert(Edn::keyword("id"), id);
    }
    entry.insert(
        Edn::keyword("donors"),
        Edn::Vector(donors.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into()),
    );
    entry.insert(Edn::keyword("acceptor"), acceptor.to_edn());
    entry.insert(Edn::keyword("type"), type_edn);
    Edn::Map(entry)
}

fn relation_entry_edn(
    id: Option<Edn<'static>>,
    atoms: &[AtomHandle],
    type_edn: Edn<'static>,
) -> Edn<'static> {
    let mut entry = EdnMap::with_capacity(3);
    if let Some(id) = id {
        entry.insert(Edn::keyword("id"), id);
    }
    entry.insert(
        Edn::keyword("atoms"),
        Edn::Vector(atoms.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into()),
    );
    entry.insert(Edn::keyword("type"), type_edn);
    Edn::Map(entry)
}

fn edit_map(entity: &str, operation: &str, payload: Edn<'static>) -> Edn<'static> {
    single_key_map(entity, single_key_map(operation, payload))
}

fn render_constraint(constraint: &Constraint) -> ConstraintDsl {
    ConstraintDsl::from_ast(constraint, &MoleculeMetadata::new())
        .expect("anonymous positional constraint rendering is infallible")
}

struct EditConstraintNamespace;

impl EditConstraintNamespace {
    fn positional_count() -> usize {
        usize::try_from(u64::from(u32::MAX) + 1).unwrap_or(usize::MAX)
    }
}

impl Namespace for EditConstraintNamespace {
    fn atom_count(&self) -> usize {
        Self::positional_count()
    }

    fn bond_count(&self) -> usize {
        Self::positional_count()
    }

    fn dative_bond_count(&self) -> usize {
        Self::positional_count()
    }

    fn aromatic_system_count(&self) -> usize {
        Self::positional_count()
    }

    fn multicenter_bond_count(&self) -> usize {
        Self::positional_count()
    }

    fn noncovalent_bond_count(&self) -> usize {
        Self::positional_count()
    }

    fn stereo_atom_count(&self) -> usize {
        Self::positional_count()
    }

    fn stereo_bond_count(&self) -> usize {
        Self::positional_count()
    }

    fn find_atom_by_keyword(&self, _keyword: &str) -> Option<AtomId> {
        None
    }

    fn find_bond_by_keyword(&self, _keyword: &str) -> Option<BondId> {
        None
    }

    fn find_dative_bond_by_keyword(&self, _keyword: &str) -> Option<DativeBondId> {
        None
    }

    fn find_aromatic_system_by_keyword(&self, _keyword: &str) -> Option<AromaticSystemId> {
        None
    }

    fn find_multicenter_bond_by_keyword(&self, _keyword: &str) -> Option<MulticenterBondId> {
        None
    }

    fn find_noncovalent_bond_by_keyword(&self, _keyword: &str) -> Option<NoncovalentBondId> {
        None
    }

    fn find_stereo_atom_by_keyword(&self, _keyword: &str) -> Option<StereoAtomId> {
        None
    }

    fn find_stereo_bond_by_keyword(&self, _keyword: &str) -> Option<StereoBondId> {
        None
    }

    fn find_bond_by_participants(&self, _first: AtomId, _second: AtomId) -> Option<BondId> {
        None
    }

    fn find_dative_bond_by_participants(
        &self,
        _donors: &[AtomId],
        _acceptor: AtomId,
    ) -> Option<DativeBondId> {
        None
    }

    fn find_aromatic_system_by_participants(&self, _atoms: &[AtomId]) -> Option<AromaticSystemId> {
        None
    }

    fn find_multicenter_bond_by_participants(
        &self,
        _atoms: &[AtomId],
    ) -> Option<MulticenterBondId> {
        None
    }

    fn find_noncovalent_bond_by_participants(
        &self,
        _first: AtomId,
        _second: AtomId,
    ) -> Option<NoncovalentBondId> {
        None
    }

    fn find_stereo_atom_by_participants(
        &self,
        _site: AtomId,
        _ligands: &[StereoLigand],
    ) -> Option<StereoAtomId> {
        None
    }

    fn find_stereo_bond_by_participants(
        &self,
        _site: BondId,
        _ligands: &[StereoLigand],
    ) -> Option<StereoBondId> {
        None
    }

    fn contains_keyword(&self, _keyword: &str) -> bool {
        false
    }

    fn find_atom_alias(&self, _name: &str) -> Option<&AtomDsl> {
        None
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::{read_string, EdnError};

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::{AtomAst, ElementAst, IsotopeMassAst};
    use crate::ast::bond::BondAst;
    use crate::ast::boolean::BooleanAst;
    use crate::ast::constraint::{
        AromaticSystemConstraintAst, AtomConstraintsAst, BondConstraintsAst,
        DativeBondConstraintAst, MoleculeConstraint, MulticenterBondConstraintAst,
        NoncovalentBondConstraintAst, RingMembershipAst, RingScope,
    };
    use crate::ast::dative::DativeBondAst;
    use crate::ast::edit::AddBond;
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{
        NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst,
    };
    use crate::ast::value::ValueAst;
    use crate::mol_dsl;

    #[rstest]
    #[case::id("7", EditHandleDsl::Id(7))]
    #[case::new("{:new 7}", EditHandleDsl::New(7))]
    fn test_edit_handle_dsl_edn_roundtrip(#[case] input: &str, #[case] expected: EditHandleDsl) {
        let parsed = EditHandleDsl::from_edn_str(input).unwrap();

        assert_eq!(parsed, expected);
        assert_eq!(parsed.to_edn(), read_string(input).unwrap());
    }

    #[rstest]
    #[case::negative_id(
        "-1",
        EdnError::De(DeError::OutOfRange {
            value: "-1".to_string(),
            target: "u32",
            path: Vec::new(),
        }),
    )]
    #[case::keyword(
        ":carbon",
        EdnError::De(DeError::TypeMismatch {
            expected: "edit handle (non-negative integer or {:new n} map)",
            got: "keyword",
            path: Vec::new(),
        }),
    )]
    #[case::structural_ref(
        "{:atoms [0 1]}",
        EdnError::De(DeError::MissingField {
            key: "new".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::empty_map(
        "{}",
        EdnError::De(DeError::MissingField {
            key: "new".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::negative_new(
        "{:new -1}",
        EdnError::De(DeError::OutOfRange {
            value: "-1".to_string(),
            target: "usize",
            path: Vec::new(),
        }),
    )]
    #[case::non_integer_new(
        "{:new :first}",
        EdnError::De(DeError::TypeMismatch {
            expected: "int",
            got: "keyword",
            path: Vec::new(),
        }),
    )]
    #[case::extra_keyword(
        "{:new 0 :atoms [0 1]}",
        EdnError::De(DeError::UnknownField {
            key: "atoms".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::extra_non_keyword(
        "{:new 0 \"atoms\" [0 1]}",
        EdnError::De(DeError::Custom(
            "edit handle map keys must be keywords".to_string(),
        )),
    )]
    fn test_edit_handle_dsl_from_edn_error(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(EditHandleDsl::from_edn_str(input), Err(expected));
    }

    #[rstest]
    #[case::id(EditHandleDsl::Id(7))]
    #[case::new(EditHandleDsl::New(7))]
    fn test_typed_edit_handle_edn_roundtrip(#[case] handle: EditHandleDsl) {
        match handle {
            EditHandleDsl::Id(index) => {
                assert_eq!(
                    AtomHandle::from_edn_str("7"),
                    Ok(AtomHandle::Id(AtomId(index)))
                );
                assert_eq!(
                    AtomHandle::Id(AtomId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    BondHandle::from_edn_str("7"),
                    Ok(BondHandle::Id(BondId(index)))
                );
                assert_eq!(
                    BondHandle::Id(BondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    DativeBondHandle::from_edn_str("7"),
                    Ok(DativeBondHandle::Id(DativeBondId(index)))
                );
                assert_eq!(
                    DativeBondHandle::Id(DativeBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    AromaticSystemHandle::from_edn_str("7"),
                    Ok(AromaticSystemHandle::Id(AromaticSystemId(index)))
                );
                assert_eq!(
                    AromaticSystemHandle::Id(AromaticSystemId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    MulticenterBondHandle::from_edn_str("7"),
                    Ok(MulticenterBondHandle::Id(MulticenterBondId(index)))
                );
                assert_eq!(
                    MulticenterBondHandle::Id(MulticenterBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    NoncovalentBondHandle::from_edn_str("7"),
                    Ok(NoncovalentBondHandle::Id(NoncovalentBondId(index)))
                );
                assert_eq!(
                    NoncovalentBondHandle::Id(NoncovalentBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    StereoAtomHandle::from_edn_str("7"),
                    Ok(StereoAtomHandle::Id(StereoAtomId(index)))
                );
                assert_eq!(
                    StereoAtomHandle::Id(StereoAtomId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    StereoBondHandle::from_edn_str("7"),
                    Ok(StereoBondHandle::Id(StereoBondId(index)))
                );
                assert_eq!(
                    StereoBondHandle::Id(StereoBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
            }
            EditHandleDsl::New(index) => {
                assert_eq!(
                    AtomHandle::from_edn_str("{:new 7}"),
                    Ok(AtomHandle::New(index))
                );
                assert_eq!(
                    AtomHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    BondHandle::from_edn_str("{:new 7}"),
                    Ok(BondHandle::New(index))
                );
                assert_eq!(
                    BondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    DativeBondHandle::from_edn_str("{:new 7}"),
                    Ok(DativeBondHandle::New(index))
                );
                assert_eq!(
                    DativeBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    AromaticSystemHandle::from_edn_str("{:new 7}"),
                    Ok(AromaticSystemHandle::New(index))
                );
                assert_eq!(
                    AromaticSystemHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    MulticenterBondHandle::from_edn_str("{:new 7}"),
                    Ok(MulticenterBondHandle::New(index))
                );
                assert_eq!(
                    MulticenterBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    NoncovalentBondHandle::from_edn_str("{:new 7}"),
                    Ok(NoncovalentBondHandle::New(index))
                );
                assert_eq!(
                    NoncovalentBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    StereoAtomHandle::from_edn_str("{:new 7}"),
                    Ok(StereoAtomHandle::New(index))
                );
                assert_eq!(
                    StereoAtomHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    StereoBondHandle::from_edn_str("{:new 7}"),
                    Ok(StereoBondHandle::New(index))
                );
                assert_eq!(
                    StereoBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
            }
        }
    }

    #[rstest]
    #[case::atom_add(
        r#"{:atom {:add "C"}}"#,
        MoleculeDefaults::ground(),
        Edits::from_iter([Edit::AddAtoms {
            atoms: vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                isotope_mass: IsotopeMassAst::Natural,
                charge: ValueAst::Lit(0),
                implicit_hydrogens: ValueAst::Lit(0),
                lone_pairs: ValueAst::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::closed_shell(),
                constraints: AtomConstraintsAst::new(),
            }],
        }]),
    )]
    #[case::atom_remove(
        "{:atom {:remove {:new 1}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::New(1)],
            bonds: Vec::new(),
        }]),
    )]
    #[case::atom_field(
        r##"{:atom {:modify [0 {:expect "#c0" :update "#c-"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(0),
                new: ValueAst::Lit(-1),
            },
        }]),
    )]
    #[case::atom_constraint_add(
        r##"{:atom {:modify [0 {:expect "#v*" :update "#v4"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintAst::valence(4_i64)),
        }]),
    )]
    #[case::atom_constraint_remove(
        r##"{:atom {:modify [0 {:expect "#v4" :update "#v*"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: Some(AtomConstraintAst::valence(4_i64)),
            new: None,
        }]),
    )]
    #[case::bond_add(
        "{:bond {:add [0 {:new 0} :single]}}",
        MoleculeDefaults::ground(),
        Edits::from_iter([Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                ast: BondAst {
                    order: ValueAst::Lit(1),
                    charge: ValueAst::Lit(0),
                    unpaired_electrons: UnpairedElectronsAst::closed_shell(),
                    constraints: BondConstraintsAst::new(),
                },
            }],
        }]),
    )]
    #[case::bond_remove(
        "{:bond {:remove 1}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![BondHandle::Id(BondId(1))],
        }]),
    )]
    #[case::bond_field(
        r#"{:bond {:modify [0 {:expect "1" :update "2"}]}}"#,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Order {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        }]),
    )]
    #[case::bond_constraint_add(
        r##"{:bond {:modify [0 {:expect "#R(6)*" :update "#R(6)"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: None,
            new: Some(BondConstraintAst::RingMembership(RingMembershipAst::new(
                RingScope::Size(6),
                ValueAst::Lit(1),
            ))),
        }]),
    )]
    #[case::topology_remove(
        "{:topology {:remove {:atoms [0 {:new 0}] :bonds [1 {:new 1}]}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            bonds: vec![BondHandle::Id(BondId(1)), BondHandle::New(1)],
        }]),
    )]
    #[case::constraint_add(
        "{:constraint {:add {:connected {}}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        }]),
    )]
    #[case::constraint_remove(
        "{:constraint {:remove {:connected {}}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::RemoveMoleculeConstraint {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        }]),
    )]
    #[case::constraint_positional(
        "{:constraint {:add {:atom [2 {:valence 4}]}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::Atom(AtomId(2), AtomConstraintAst::valence(4_i64)),
        }]),
    )]
    fn test_edit_input_append_to(
        #[case] input: &str,
        #[case] defaults: MoleculeDefaults,
        #[case] expected: Edits,
    ) {
        let parsed = EditInput::from_edn_str(input).unwrap();
        let mut edits = Edits::new();

        parsed.append_to(&mut edits, &defaults).unwrap();

        assert_eq!(edits, expected);
    }

    #[rstest]
    #[case::atom_add(
        Edit::AddAtoms {
            atoms: vec![AtomAst::from_element(Element::C).into_ground()],
        },
        MoleculeDefaults::ground(),
        r#"{:atom {:add "C"}}"#,
    )]
    #[case::bond_add(
        Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                ast: BondAst::from_order(1).into_ground(),
            }],
        },
        MoleculeDefaults::ground(),
        "{:bond {:add [0 {:new 0} :single]}}",
    )]
    #[case::atom_remove(
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        MoleculeDefaults::new(),
        "{:atom {:remove 0}}",
    )]
    #[case::bond_remove(
        Edit::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![BondHandle::New(0)],
        },
        MoleculeDefaults::new(),
        "{:bond {:remove {:new 0}}}",
    )]
    #[case::topology_remove(
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: vec![BondHandle::Id(BondId(0))],
        },
        MoleculeDefaults::new(),
        "{:topology {:remove {:atoms [0] :bonds [0]}}}",
    )]
    #[case::atom_field(
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(0),
                new: ValueAst::Lit(-1),
            },
        },
        MoleculeDefaults::new(),
        r##"{:atom {:modify [0 {:expect "#c0" :update "#c-"}]}}"##,
    )]
    #[case::atom_constraint_add(
        Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintAst::valence(4_i64)),
        },
        MoleculeDefaults::new(),
        r##"{:atom {:modify [0 {:expect "#v*" :update "#v4"}]}}"##,
    )]
    #[case::constraint_add(
        Edit::AddMoleculeConstraint {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        },
        MoleculeDefaults::new(),
        "{:constraint {:add {:connected {}}}}",
    )]
    fn test_edit_input_from_edit(
        #[case] edit: Edit,
        #[case] defaults: MoleculeDefaults,
        #[case] expected: &str,
    ) {
        let rendered = EditInput::from_edit(&edit, &defaults)
            .unwrap()
            .unwrap()
            .into_iter()
            .map(|input| input.to_edn())
            .collect::<Vec<_>>();

        assert_eq!(rendered, vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_add(
        r#"{:dative-bond {:add {:donors [0 {:new 0}] :acceptor 2 :type :single}}}"#,
        Edit::AddDativeBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            ast: DativeBondAst::from_order(1),
        },
    )]
    #[case::dative_remove(
        r#"{:dative-bonds {:remove [{:id 0 :donors [1] :acceptor {:new 2} :type :single} {:id {:new 0} :donors [{:new 1}] :acceptor 3 :type :double}]}}"#,
        Edit::RemoveDativeBonds { removes: vec![
            (DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(1)), AtomHandle::New(2)], DativeBondAst::from_order(1)),
            (DativeBondHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(3))], DativeBondAst::from_order(2)),
        ] },
    )]
    #[case::dative_field(
        r#"{:dative-bond {:modify [{:new 0} {:expect "1" :update "2"}]}}"#,
        Edit::ModifyDativeBondField {
            id: DativeBondHandle::New(0),
            change: DativeBondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
        },
    )]
    #[case::dative_constraint(
        r##"{:dative-bond {:modify [0 {:expect "#a*" :update "#a"}]}}"##,
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::Id(DativeBondId(0)),
            old: None,
            new: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        },
    )]
    #[case::dative_constraint_remove(
        r##"{:dative-bond {:modify [{:new 1} {:expect "#a" :update "#a*"}]}}"##,
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::New(1),
            old: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
            new: None,
        },
    )]
    #[case::aromatic_add(
        r#"{:aromatic-system {:add {:atoms [0 {:new 0} 2] :type "[1,1,1]"}}}"#,
        Edit::AddAromaticSystem {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            ast: AromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
    )]
    #[case::aromatic_remove(
        r#"{:aromatic-systems {:remove [{:id 0 :atoms [0 {:new 0}] :type "[1,1]"} {:id {:new 0} :atoms [{:new 1} 2] :type "[2,2]"}]}}"#,
        Edit::RemoveAromaticSystems { removes: vec![
            (AromaticSystemHandle::Id(AromaticSystemId(0)), vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], AromaticSystemAst::from_electrons(vec![1, 1])),
            (AromaticSystemHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(2))], AromaticSystemAst::from_electrons(vec![2, 2])),
        ] },
    )]
    #[case::aromatic_field(
        r##"{:aromatic-system {:modify [{:new 0} {:expect "#c0" :update "#c-"}]}}"##,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::New(0),
            change: AromaticSystemFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
        },
    )]
    #[case::aromatic_electrons(
        r#"{:aromatic-system {:modify [1 {:expect "[1,1]" :update "[2,0]"}]}}"#,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(1)),
            change: AromaticSystemFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1]), new: ElectronCountsAst::Lit(vec![2, 0]) },
        },
    )]
    #[case::aromatic_unpaired_electrons(
        r##"{:aromatic-system {:modify [{:new 1} {:expect "#u0#s" :update "#u2#s3"}]}}"##,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::New(1),
            change: AromaticSystemFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(1) },
                new: UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) },
            },
        },
    )]
    #[case::aromatic_constraint(
        r##"{:aromatic-system {:modify [0 {:expect "#e*" :update "#e6"}]}}"##,
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            old: None,
            new: Some(AromaticSystemConstraintAst::electron_count(ValueAst::Lit(6))),
        },
    )]
    #[case::aromatic_constraint_remove(
        r##"{:aromatic-system {:modify [{:new 2} {:expect "#e6" :update "#e*"}]}}"##,
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::New(2),
            old: Some(AromaticSystemConstraintAst::electron_count(ValueAst::Lit(6))),
            new: None,
        },
    )]
    #[case::multicenter_add(
        r#"{:multicenter-bond {:add {:atoms [0 {:new 0} 2] :type "[1,1,0]"}}}"#,
        Edit::AddMulticenterBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            ast: MulticenterBondAst::from_electrons(vec![1, 1, 0]),
        },
    )]
    #[case::multicenter_remove(
        r#"{:multicenter-bonds {:remove [{:id 0 :atoms [0 {:new 0}] :type "[1,1]"} {:id {:new 0} :atoms [{:new 1} 2] :type "[2,0]"}]}}"#,
        Edit::RemoveMulticenterBonds { removes: vec![
            (MulticenterBondHandle::Id(MulticenterBondId(0)), vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], MulticenterBondAst::from_electrons(vec![1, 1])),
            (MulticenterBondHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(2))], MulticenterBondAst::from_electrons(vec![2, 0])),
        ] },
    )]
    #[case::multicenter_field(
        r#"{:multicenter-bond {:modify [{:new 0} {:expect "[1,1]" :update "[2,0]"}]}}"#,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::New(0),
            change: MulticenterBondFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1]), new: ElectronCountsAst::Lit(vec![2, 0]) },
        },
    )]
    #[case::multicenter_charge(
        r##"{:multicenter-bond {:modify [1 {:expect "#c0" :update "#c+"}]}}"##,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(1)),
            change: MulticenterBondFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(1) },
        },
    )]
    #[case::multicenter_unpaired_electrons(
        r##"{:multicenter-bond {:modify [{:new 1} {:expect "#u0#s" :update "#u2#s3"}]}}"##,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::New(1),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(1) },
                new: UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) },
            },
        },
    )]
    #[case::multicenter_constraint(
        r##"{:multicenter-bond {:modify [0 {:expect "#e*" :update "#e2"}]}}"##,
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            old: None,
            new: Some(MulticenterBondConstraintAst::electron_count(ValueAst::Lit(2))),
        },
    )]
    #[case::multicenter_constraint_remove(
        r##"{:multicenter-bond {:modify [{:new 2} {:expect "#e2" :update "#e*"}]}}"##,
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::New(2),
            old: Some(MulticenterBondConstraintAst::electron_count(ValueAst::Lit(2))),
            new: None,
        },
    )]
    #[case::noncovalent_add(
        r#"{:noncovalent-bond {:add {:atoms [0 {:new 0}] :type "Hbd"}}}"#,
        Edit::AddNoncovalentBond {
            atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        },
    )]
    #[case::noncovalent_remove(
        r#"{:noncovalent-bonds {:remove [{:id 0 :atoms [0 {:new 0}] :type "Hbd"} {:id {:new 0} :atoms [{:new 1} 2] :type "Ion"}]}}"#,
        Edit::RemoveNoncovalentBonds { removes: vec![
            (NoncovalentBondHandle::Id(NoncovalentBondId(0)), [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)),
            (NoncovalentBondHandle::New(0), [AtomHandle::New(1), AtomHandle::Id(AtomId(2))], NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic)),
        ] },
    )]
    #[case::noncovalent_field(
        r#"{:noncovalent-bond {:modify [{:new 0} {:expect "Hbd" :update "Ion"}]}}"#,
        Edit::ModifyNoncovalentBondField {
            id: NoncovalentBondHandle::New(0),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
            },
        },
    )]
    #[case::noncovalent_constraint(
        r##"{:noncovalent-bond {:modify [0 {:expect "#I*" :update "#I"}]}}"##,
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            old: None,
            new: Some(NoncovalentBondConstraintAst::intramolecular(true)),
        },
    )]
    #[case::noncovalent_constraint_remove(
        r##"{:noncovalent-bond {:modify [{:new 1} {:expect "#I" :update "#I*"}]}}"##,
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::New(1),
            old: Some(NoncovalentBondConstraintAst::intramolecular(true)),
            new: None,
        },
    )]
    fn test_overlay_edit_input_roundtrip(#[case] input: &str, #[case] expected: Edit) {
        let parsed = EditInput::from_edn_str(input).unwrap();
        let mut edits = Edits::new();

        parsed.append_to(&mut edits, &MoleculeDefaults::new()).unwrap();
        let rendered = EditInput::from_edit(&expected, &MoleculeDefaults::new())
            .unwrap()
            .unwrap()
            .into_iter()
            .map(|input| input.to_edn())
            .collect::<Vec<_>>();

        assert_eq!(edits, Edits::from_iter([expected]));
        assert_eq!(rendered, vec![read_string(input).unwrap()]);
    }

    #[rstest]
    #[case::aromatic_add(
        r#"{:aromatic-system {:add {:atoms [0 1] :type "[1,1]"}}}"#,
        Edit::AddAromaticSystem {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            ast: AromaticSystemAst::from_electrons(vec![1, 1]).into_ground(),
        },
    )]
    #[case::aromatic_remove(
        r#"{:aromatic-systems {:remove [{:id 0 :atoms [0 1] :type "[1,1]"}]}}"#,
        Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                AromaticSystemAst::from_electrons(vec![1, 1]).into_ground(),
            )],
        },
    )]
    #[case::multicenter_add(
        r#"{:multicenter-bond {:add {:atoms [0 1] :type "[1,1]"}}}"#,
        Edit::AddMulticenterBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            ast: MulticenterBondAst::from_electrons(vec![1, 1]).into_ground(),
        },
    )]
    #[case::multicenter_remove(
        r#"{:multicenter-bonds {:remove [{:id 0 :atoms [0 1] :type "[1,1]"}]}}"#,
        Edit::RemoveMulticenterBonds {
            removes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                MulticenterBondAst::from_electrons(vec![1, 1]).into_ground(),
            )],
        },
    )]
    fn test_overlay_edit_input_ground_defaults(#[case] input: &str, #[case] expected: Edit) {
        let mut edits = Edits::new();

        EditInput::from_edn_str(input)
            .unwrap()
            .append_to(&mut edits, &MoleculeDefaults::ground())
            .unwrap();
        let rendered = EditInput::from_edit(&expected, &MoleculeDefaults::ground())
            .unwrap()
            .unwrap()
            .into_iter()
            .map(|input| input.to_edn())
            .collect::<Vec<_>>();

        assert_eq!(edits, Edits::from_iter([expected]));
        assert_eq!(rendered, vec![read_string(input).unwrap()]);
    }

    #[rstest]
    #[case::atom_field(
        r##"{:atom {:modify [0 {:expect "#c0" :update "#h0"}]}}"##,
        EdnError::De(DeError::Custom(
            "atom :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::atom_constraint(
        r##"{:atom {:modify [0 {:expect "#v4" :update "#d4"}]}}"##,
        EdnError::De(DeError::Custom(
            "atom :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::atom_spin(
        r##"{:atom {:modify [0 {:expect "#u2" :update "#u0"}]}}"##,
        EdnError::De(DeError::Custom(
            "atom :modify unpaired-electron changes require both #u and #s".to_string(),
        )),
    )]
    #[case::bond_field(
        r##"{:bond {:modify [0 {:expect "1" :update "#c0"}]}}"##,
        EdnError::De(DeError::Custom(
            "bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::dative_field(
        r##"{:dative-bond {:modify [0 {:expect "1" :update "#a"}]}}"##,
        EdnError::De(DeError::Custom(
            "dative-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::aromatic_field(
        r##"{:aromatic-system {:modify [0 {:expect "#c0" :update "[1,1]"}]}}"##,
        EdnError::De(DeError::Custom(
            "aromatic-system :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::aromatic_spin(
        r##"{:aromatic-system {:modify [0 {:expect "#u0" :update "#u2"}]}}"##,
        EdnError::De(DeError::Custom(
            "aromatic-system :modify unpaired-electron changes require both #u and #s".to_string(),
        )),
    )]
    #[case::multicenter_field(
        r##"{:multicenter-bond {:modify [0 {:expect "#c0" :update "[1,1]"}]}}"##,
        EdnError::De(DeError::Custom(
            "multicenter-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::multicenter_spin(
        r##"{:multicenter-bond {:modify [0 {:expect "#u0" :update "#u2"}]}}"##,
        EdnError::De(DeError::Custom(
            "multicenter-bond :modify unpaired-electron changes require both #u and #s".to_string(),
        )),
    )]
    #[case::noncovalent_field(
        r##"{:noncovalent-bond {:modify [0 {:expect "Hbd" :update "#I"}]}}"##,
        EdnError::De(DeError::Custom(
            "noncovalent-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        )),
    )]
    #[case::constraint_keyword(
        "{:constraint {:add {:atom [:carbon {:valence 4}]}}}",
        EdnError::De(DeError::Subgrammar {
            grammar: "standalone-edit-constraint",
            message: "invalid atom ref: carbon".to_string(),
            path: Vec::new(),
        }),
    )]
    fn test_edit_input_from_edn_error(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(EditInput::from_edn_str(input), Err(expected));
    }

    #[rstest]
    #[case::incident_bond(
        mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
        "{:topology {:remove {:atoms [0] :bonds [0]}}}",
        mol_dsl!(r#"{:atoms ["N"]}"#),
    )]
    fn test_edit_input_append_to_topology(
        #[case] molecule: MoleculeAst,
        #[case] input: &str,
        #[case] expected: MoleculeAst,
    ) {
        let mut edits = Edits::new();
        EditInput::from_edn_str(input)
            .unwrap()
            .append_to(&mut edits, &MoleculeDefaults::new())
            .unwrap();

        assert_eq!(molecule.apply(edits), Ok(expected));
    }
}
