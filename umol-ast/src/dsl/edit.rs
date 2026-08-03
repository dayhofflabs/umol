//! Surface encoding for handles in standalone edit documents.

#![expect(
    dead_code,
    reason = "edit-family codecs are assembled by the aggregate EditsDsl root"
)]

use umol_edn::{DeError, Edn, EdnMap, EdnMapHelper, FromEdn, ToEdn};

use super::atom::{AtomDsl, AtomUpdateDsl};
use super::bond::{BondDsl, BondUpdateDsl};
use super::config::MoleculeDefaults;
use super::constraint::ConstraintDsl;
use super::edn_utils::{parse_single_key_map, single_key_map};
use super::metadata::MoleculeMetadata;
use super::namespace::Namespace;
use crate::ast::atom::AtomUpdate;
use crate::ast::bond::BondUpdate;
use crate::ast::constraint::{AtomConstraintAst, BondConstraintAst, Constraint};
use crate::ast::edit::{
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    DativeBondHandle, Edit, Edits, MulticenterBondHandle, NoncovalentBondHandle, StereoAtomHandle,
    StereoBondHandle,
};
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::StereoLigand;
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
    use crate::ast::atom::{AtomAst, ElementAst, IsotopeMassAst};
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{
        AtomConstraintsAst, BondConstraintsAst, MoleculeConstraint, RingMembershipAst, RingScope,
    };
    use crate::ast::edit::AddBond;
    use crate::ast::molecule::MoleculeAst;
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
