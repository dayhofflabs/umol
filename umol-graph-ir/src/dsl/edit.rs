//! Surface encoding for standalone edit documents.

use std::convert::Infallible;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnMap, EdnMapHelper, FromEdn, ToEdn};

use super::aromatic::{AromaticSystemConstraintDsl, AromaticSystemDsl, AromaticSystemUpdateDsl};
use super::atom::{AtomConstraintDsl, AtomDsl, AtomUpdateDsl};
use super::bond::{BondConstraintDsl, BondDsl, BondUpdateDsl};
use super::config::MoleculeDefaults;
use super::constraint::{expect_map, parse_unpaired_electrons, render_unpaired_electrons};
use super::dative::{DativeBondConstraintDsl, DativeBondDsl, DativeBondUpdateDsl};
use super::edn_utils::{parse_single_key_map, single_key_map};
use super::error::ParseError;
use super::multicenter::{
    MulticenterBondConstraintDsl, MulticenterBondDsl, MulticenterBondUpdateDsl,
};
use super::noncovalent::{
    NoncovalentBondConstraintDsl, NoncovalentBondDsl, NoncovalentBondUpdateDsl,
};
use super::num::NumDsl;
use super::stereo::{
    StereoAtomConstraintDsl, StereoAtomDsl, StereoAtomUpdateDsl, StereoBondConstraintDsl,
    StereoBondDsl, StereoBondUpdateDsl,
};
use crate::ir::aromatic::AromaticSystemUpdate;
use crate::ir::atom::AtomUpdate;
use crate::ir::bond::BondUpdate;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
    DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
    NoncovalentBondConstraintForm, RelationalConstraint, StereoAtomConstraintForm,
    StereoBondConstraintForm,
};
use crate::ir::dative::DativeBondUpdate;
use crate::ir::edit::{
    AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange,
    BondHandle, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit, Edits, EntityHandle,
    MulticenterBondFieldChange, MulticenterBondHandle, NoncovalentBondFieldChange,
    NoncovalentBondHandle, StereoAtomFieldChange, StereoAtomHandle, StereoBondFieldChange,
    StereoBondHandle,
};
use crate::ir::entity::Entity;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::StereoLigandKind;
use crate::ir::multicenter::MulticenterBondUpdate;
use crate::ir::noncovalent::NoncovalentBondUpdate;
use crate::ir::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use crate::ir::stereo::{
    StereoAtomUpdate, StereoBondUpdate, StereoConfigurationForm, StereoConfigurationUpdate,
    StereoKind,
};
use crate::ir::traits::{FromIr, IntoIr, Lattice};

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

type StereoLigandInput = (AtomHandle, StereoLigandKind);
type StereoAtomAdditionInput = (AtomHandle, Vec<StereoLigandInput>, StereoAtomDsl);
type StereoAtomRemovalInput = (
    StereoAtomHandle,
    AtomHandle,
    Vec<StereoLigandInput>,
    StereoAtomDsl,
);
type StereoBondAdditionInput = (BondHandle, Vec<StereoLigandInput>, StereoBondDsl);
type StereoBondRemovalInput = (
    StereoBondHandle,
    BondHandle,
    Vec<StereoLigandInput>,
    StereoBondDsl,
);

/// Ordered standalone surface form for a batch of host-specific molecule edits.
///
/// Parsing validates each checked update and recorded removal shape before the batch can be
/// converted to [`Edits`]. Full entity definitions and recorded removal state are interpreted under
/// [`MoleculeDefaults`]; partial `:expect` and `:update` values are not defaulted.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EditsDsl {
    inputs: Vec<EditInput>,
}

impl FromStr for EditsDsl {
    type Err = ParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(input).map_err(|error| ParseError::EdnParse(error.to_string()))
    }
}

impl Display for EditsDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for EditsDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Vector(entries) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "vector of edits",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        Ok(Self {
            inputs: entries
                .iter()
                .map(EditInput::from_edn)
                .collect::<Result<_, _>>()?,
        })
    }
}

impl ToEdn for EditsDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(
            self.inputs
                .iter()
                .map(ToEdn::to_edn)
                .collect::<Vec<_>>()
                .into(),
        )
    }
}

impl FromIr<Edits> for EditsDsl {
    type Context = MoleculeDefaults;

    fn from_ir(edits: &Edits, context: &Self::Context) -> Self {
        let mut inputs = Vec::new();
        for edit in edits.iter() {
            inputs.extend(
                EditInput::from_edit(edit, context)
                    .expect("Edit variants satisfy their representational invariants"),
            );
        }
        Self { inputs }
    }
}

impl IntoIr<Edits> for EditsDsl {
    type Context = MoleculeDefaults;

    fn into_ir(self, context: &Self::Context) -> Edits {
        let mut edits = Edits::new();
        for input in self.inputs {
            input
                .append_to(&mut edits, context)
                .expect("EditsDsl stores only validated edit inputs");
        }
        edits
    }
}

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
        attributes: BondDsl,
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
        attributes: DativeBondDsl,
    },
    DativeBondsRemove(Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondDsl)>),
    DativeBondModify {
        id: DativeBondHandle,
        expect: DativeBondUpdate,
        update: DativeBondUpdate,
    },
    AromaticSystemAdd {
        atoms: Vec<AtomHandle>,
        attributes: AromaticSystemDsl,
    },
    AromaticSystemsRemove(Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemDsl)>),
    AromaticSystemModify {
        id: AromaticSystemHandle,
        expect: AromaticSystemUpdate,
        update: AromaticSystemUpdate,
    },
    MulticenterBondAdd {
        atoms: Vec<AtomHandle>,
        attributes: MulticenterBondDsl,
    },
    MulticenterBondsRemove(Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondDsl)>),
    MulticenterBondModify {
        id: MulticenterBondHandle,
        expect: MulticenterBondUpdate,
        update: MulticenterBondUpdate,
    },
    NoncovalentBondAdd {
        atoms: [AtomHandle; 2],
        attributes: NoncovalentBondDsl,
    },
    NoncovalentBondsRemove(Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondDsl)>),
    NoncovalentBondModify {
        id: NoncovalentBondHandle,
        expect: NoncovalentBondUpdate,
        update: NoncovalentBondUpdate,
    },
    StereoAtomAdd {
        site: AtomHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoAtomDsl,
    },
    StereoAtomsRemove(Vec<StereoAtomRemovalInput>),
    StereoAtomModify {
        id: StereoAtomHandle,
        expect: StereoAtomUpdate,
        update: StereoAtomUpdate,
    },
    StereoBondAdd {
        site: BondHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoBondDsl,
    },
    StereoBondsRemove(Vec<StereoBondRemovalInput>),
    StereoBondModify {
        id: StereoBondHandle,
        expect: StereoBondUpdate,
        update: StereoBondUpdate,
    },
    TopologyRemove {
        atoms: Vec<AtomHandle>,
        bonds: Vec<BondHandle>,
    },
    ConstraintAdd(ConstraintEdit),
    ConstraintRemove(ConstraintEdit),
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
            "stereo-atom" => parse_stereo_atom_edit(body),
            "stereo-atoms" => parse_stereo_atoms_edit(body),
            "stereo-bond" => parse_stereo_bond_edit(body),
            "stereo-bonds" => parse_stereo_bonds_edit(body),
            "topology" => parse_topology_edit(body),
            "constraint" => parse_constraint_edit(body),
            other => Err(DeError::Custom(format!("unknown edit :{other}"))),
        }
    }
}

impl ToEdn for EditInput {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::AtomAdd(attributes) => edit_map("atom", "add", attributes.to_edn()),
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
                attributes,
            } => edit_map(
                "bond",
                "add",
                Edn::Vector(vec![first.to_edn(), second.to_edn(), attributes.to_edn()].into()),
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
                attributes,
            } => edit_map(
                "dative-bond",
                "add",
                dative_entry_edn(None, donors, acceptor, attributes.to_edn()),
            ),
            Self::DativeBondsRemove(removes) => edit_map(
                "dative-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, attributes)| {
                            let (acceptor, donors) = atoms
                                .split_last()
                                .expect("dative edit always has an acceptor");
                            dative_entry_edn(
                                Some(id.to_edn()),
                                donors,
                                acceptor,
                                attributes.to_edn(),
                            )
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
            Self::AromaticSystemAdd { atoms, attributes } => edit_map(
                "aromatic-system",
                "add",
                relation_entry_edn(None, atoms, attributes.to_edn()),
            ),
            Self::AromaticSystemsRemove(removes) => edit_map(
                "aromatic-systems",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, attributes)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, attributes.to_edn())
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
            Self::MulticenterBondAdd { atoms, attributes } => edit_map(
                "multicenter-bond",
                "add",
                relation_entry_edn(None, atoms, attributes.to_edn()),
            ),
            Self::MulticenterBondsRemove(removes) => edit_map(
                "multicenter-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, attributes)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, attributes.to_edn())
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
            Self::NoncovalentBondAdd { atoms, attributes } => edit_map(
                "noncovalent-bond",
                "add",
                relation_entry_edn(None, atoms, attributes.to_edn()),
            ),
            Self::NoncovalentBondsRemove(removes) => edit_map(
                "noncovalent-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, atoms, attributes)| {
                            relation_entry_edn(Some(id.to_edn()), atoms, attributes.to_edn())
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
            Self::StereoAtomAdd {
                site,
                ligands,
                attributes,
            } => edit_map(
                "stereo-atom",
                "add",
                stereo_entry_edn(None, site.to_edn(), ligands, attributes.to_edn()),
            ),
            Self::StereoAtomsRemove(removes) => edit_map(
                "stereo-atoms",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, site, ligands, attributes)| {
                            stereo_entry_edn(
                                Some(id.to_edn()),
                                site.to_edn(),
                                ligands,
                                attributes.to_edn(),
                            )
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::StereoAtomModify { id, expect, update } => edit_map(
                "stereo-atom",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    StereoAtomUpdateDsl(expect.clone()).to_edn(),
                    StereoAtomUpdateDsl(update.clone()).to_edn(),
                ),
            ),
            Self::StereoBondAdd {
                site,
                ligands,
                attributes,
            } => edit_map(
                "stereo-bond",
                "add",
                stereo_entry_edn(None, site.to_edn(), ligands, attributes.to_edn()),
            ),
            Self::StereoBondsRemove(removes) => edit_map(
                "stereo-bonds",
                "remove",
                Edn::Vector(
                    removes
                        .iter()
                        .map(|(id, site, ligands, attributes)| {
                            stereo_entry_edn(
                                Some(id.to_edn()),
                                site.to_edn(),
                                ligands,
                                attributes.to_edn(),
                            )
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ),
            Self::StereoBondModify { id, expect, update } => edit_map(
                "stereo-bond",
                "modify",
                checked_update_edn(
                    id.to_edn(),
                    StereoBondUpdateDsl(expect.clone()).to_edn(),
                    StereoBondUpdateDsl(update.clone()).to_edn(),
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
            Self::ConstraintAdd(constraint) => edit_map(
                "constraint",
                "add",
                ConstraintEditDsl::from_edit(constraint.clone()).to_edn(),
            ),
            Self::ConstraintRemove(constraint) => edit_map(
                "constraint",
                "remove",
                ConstraintEditDsl::from_edit(constraint.clone()).to_edn(),
            ),
        }
    }
}

impl EditInput {
    fn append_to(self, edits: &mut Edits, defaults: &MoleculeDefaults) -> Result<(), DeError> {
        match self {
            Self::AtomAdd(attributes) => {
                edits.add_atom(attributes.into_ir(&defaults.atom));
            }
            Self::AtomRemove(id) => edits.remove_atom(id),
            Self::AtomModify { id, expect, update } => {
                append_atom_modify(edits, id, expect, update)?;
            }
            Self::BondAdd {
                atoms: [first, second],
                attributes,
            } => {
                edits.add_bond(first, second, attributes.into_ir(&defaults.bond));
            }
            Self::BondRemove(id) => edits.remove_bond(id),
            Self::BondModify { id, expect, update } => {
                append_bond_modify(edits, id, expect, update)?;
            }
            Self::DativeBondAdd {
                mut donors,
                acceptor,
                attributes,
            } => {
                donors.push(acceptor);
                edits.add_dative_bond(donors, attributes.into_ir(&defaults.dative_bond));
            }
            Self::DativeBondsRemove(removes) => edits.remove_dative_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, attributes)| {
                        (id, atoms, attributes.into_ir(&defaults.dative_bond))
                    })
                    .collect(),
            ),
            Self::DativeBondModify { id, expect, update } => {
                append_dative_bond_modify(edits, id, expect, update)?;
            }
            Self::AromaticSystemAdd { atoms, attributes } => {
                edits.add_aromatic_system(atoms, attributes.into_ir(&defaults.aromatic_system));
            }
            Self::AromaticSystemsRemove(removes) => edits.remove_aromatic_systems(
                removes
                    .into_iter()
                    .map(|(id, atoms, attributes)| {
                        (id, atoms, attributes.into_ir(&defaults.aromatic_system))
                    })
                    .collect(),
            ),
            Self::AromaticSystemModify { id, expect, update } => {
                append_aromatic_system_modify(edits, id, expect, update)?;
            }
            Self::MulticenterBondAdd { atoms, attributes } => {
                edits.add_multicenter_bond(atoms, attributes.into_ir(&defaults.multicenter_bond));
            }
            Self::MulticenterBondsRemove(removes) => edits.remove_multicenter_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, attributes)| {
                        (id, atoms, attributes.into_ir(&defaults.multicenter_bond))
                    })
                    .collect(),
            ),
            Self::MulticenterBondModify { id, expect, update } => {
                append_multicenter_bond_modify(edits, id, expect, update)?;
            }
            Self::NoncovalentBondAdd { atoms, attributes } => {
                edits.add_noncovalent_bond(atoms, attributes.into_ir(&defaults.noncovalent_bond));
            }
            Self::NoncovalentBondsRemove(removes) => edits.remove_noncovalent_bonds(
                removes
                    .into_iter()
                    .map(|(id, atoms, attributes)| {
                        (id, atoms, attributes.into_ir(&defaults.noncovalent_bond))
                    })
                    .collect(),
            ),
            Self::NoncovalentBondModify { id, expect, update } => {
                append_noncovalent_bond_modify(edits, id, expect, update)?;
            }
            Self::StereoAtomAdd {
                site,
                ligands,
                attributes,
            } => {
                edits.add_stereo_atom(site, ligands, attributes.into_ir(&defaults.stereo_atom));
            }
            Self::StereoAtomsRemove(removes) => edits.remove_stereo_atoms(
                removes
                    .into_iter()
                    .map(|(id, site, ligands, attributes)| {
                        (id, site, ligands, attributes.into_ir(&defaults.stereo_atom))
                    })
                    .collect(),
            ),
            Self::StereoAtomModify { id, expect, update } => {
                append_stereo_atom_modify(edits, id, expect, update)?;
            }
            Self::StereoBondAdd {
                site,
                ligands,
                attributes,
            } => {
                edits.add_stereo_bond(site, ligands, attributes.into_ir(&defaults.stereo_bond));
            }
            Self::StereoBondsRemove(removes) => edits.remove_stereo_bonds(
                removes
                    .into_iter()
                    .map(|(id, site, ligands, attributes)| {
                        (id, site, ligands, attributes.into_ir(&defaults.stereo_bond))
                    })
                    .collect(),
            ),
            Self::StereoBondModify { id, expect, update } => {
                append_stereo_bond_modify(edits, id, expect, update)?;
            }
            Self::TopologyRemove { atoms, bonds } => edits.remove_topology(atoms, bonds),
            Self::ConstraintAdd(constraint) => edits.add_molecule_constraint(constraint),
            Self::ConstraintRemove(constraint) => edits.remove_molecule_constraint(constraint),
        }
        Ok(())
    }

    fn from_edit(edit: &Edit, defaults: &MoleculeDefaults) -> Result<Vec<Self>, DeError> {
        let inputs = match edit {
            Edit::AddAtoms { atoms } => atoms
                .iter()
                .map(|attributes| Self::AtomAdd(AtomDsl::from_ir(attributes, &defaults.atom)))
                .collect(),
            Edit::AddBonds { bonds } => bonds
                .iter()
                .map(|bond| Self::BondAdd {
                    atoms: bond.endpoints.clone(),
                    attributes: BondDsl::from_ir(&bond.attributes, &defaults.bond),
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
            Edit::AddDativeBond { atoms, attributes } => {
                let (acceptor, donors) = atoms.split_last().ok_or_else(|| {
                    DeError::Custom("dative-bond addition has no acceptor".to_string())
                })?;
                vec![Self::DativeBondAdd {
                    donors: donors.to_vec(),
                    acceptor: acceptor.clone(),
                    attributes: DativeBondDsl::from_ir(attributes, &defaults.dative_bond),
                }]
            }
            Edit::RemoveDativeBonds { removes } => vec![Self::DativeBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        if atoms.is_empty() {
                            return Err(DeError::Custom(
                                "dative-bond removal has no acceptor".to_string(),
                            ));
                        }
                        Ok((
                            id.clone(),
                            atoms.clone(),
                            DativeBondDsl::from_ir(attributes, &defaults.dative_bond),
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
            Edit::AddAromaticSystem { atoms, attributes } => vec![Self::AromaticSystemAdd {
                atoms: atoms.clone(),
                attributes: AromaticSystemDsl::from_ir(attributes, &defaults.aromatic_system),
            }],
            Edit::RemoveAromaticSystems { removes } => vec![Self::AromaticSystemsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            AromaticSystemDsl::from_ir(attributes, &defaults.aromatic_system),
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
            Edit::AddMulticenterBond { atoms, attributes } => vec![Self::MulticenterBondAdd {
                atoms: atoms.clone(),
                attributes: MulticenterBondDsl::from_ir(attributes, &defaults.multicenter_bond),
            }],
            Edit::RemoveMulticenterBonds { removes } => vec![Self::MulticenterBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            MulticenterBondDsl::from_ir(attributes, &defaults.multicenter_bond),
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
            Edit::AddNoncovalentBond { atoms, attributes } => vec![Self::NoncovalentBondAdd {
                atoms: atoms.clone(),
                attributes: NoncovalentBondDsl::from_ir(attributes, &defaults.noncovalent_bond),
            }],
            Edit::RemoveNoncovalentBonds { removes } => vec![Self::NoncovalentBondsRemove(
                removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.clone(),
                            atoms.clone(),
                            NoncovalentBondDsl::from_ir(attributes, &defaults.noncovalent_bond),
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
            Edit::AddStereoAtom {
                site,
                ligands,
                attributes,
            } => vec![Self::StereoAtomAdd {
                site: site.clone(),
                ligands: ligands.clone(),
                attributes: StereoAtomDsl::from_ir(attributes, &defaults.stereo_atom),
            }],
            Edit::RemoveStereoAtoms { removes } => vec![Self::StereoAtomsRemove(
                removes
                    .iter()
                    .map(|(id, site, ligands, attributes)| {
                        (
                            id.clone(),
                            site.clone(),
                            ligands.clone(),
                            StereoAtomDsl::from_ir(attributes, &defaults.stereo_atom),
                        )
                    })
                    .collect(),
            )],
            Edit::ModifyStereoAtomField { id, change } => {
                let (expect, update) = stereo_atom_field_updates(change);
                vec![Self::StereoAtomModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyStereoAtomConstraint { id, kind, old, new } => {
                let (expect, update) = stereo_atom_constraint_updates(*kind, old, new)?;
                vec![Self::StereoAtomModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::AddStereoBond {
                site,
                ligands,
                attributes,
            } => vec![Self::StereoBondAdd {
                site: site.clone(),
                ligands: ligands.clone(),
                attributes: StereoBondDsl::from_ir(attributes, &defaults.stereo_bond),
            }],
            Edit::RemoveStereoBonds { removes } => vec![Self::StereoBondsRemove(
                removes
                    .iter()
                    .map(|(id, site, ligands, attributes)| {
                        (
                            id.clone(),
                            site.clone(),
                            ligands.clone(),
                            StereoBondDsl::from_ir(attributes, &defaults.stereo_bond),
                        )
                    })
                    .collect(),
            )],
            Edit::ModifyStereoBondField { id, change } => {
                let (expect, update) = stereo_bond_field_updates(change);
                vec![Self::StereoBondModify {
                    id: id.clone(),
                    expect,
                    update,
                }]
            }
            Edit::ModifyStereoBondConstraint { id, kind, old, new } => {
                let (expect, update) = stereo_bond_constraint_updates(*kind, old, new)?;
                vec![Self::StereoBondModify {
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
        };
        Ok(inputs)
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
                attributes: BondDsl::from_edn(&parts[2])?,
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
            let (donors, acceptor, attributes) = parse_dative_bond_addition(payload)?;
            Ok(EditInput::DativeBondAdd {
                donors,
                acceptor,
                attributes,
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
            let (atoms, attributes) = parse_aromatic_system_addition(payload)?;
            Ok(EditInput::AromaticSystemAdd { atoms, attributes })
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
            let (atoms, attributes) = parse_multicenter_bond_addition(payload)?;
            Ok(EditInput::MulticenterBondAdd { atoms, attributes })
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
            let (atoms, attributes) = parse_noncovalent_bond_addition(payload)?;
            Ok(EditInput::NoncovalentBondAdd { atoms, attributes })
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((donors, acceptor, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    donors.push(acceptor);
    Ok((id, donors, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((atoms, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((id, atoms, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((atoms, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((id, atoms, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((atoms, attributes))
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
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((id, atoms, attributes))
}

fn parse_stereo_atom_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-atom edit")?;
    match op {
        "add" => {
            let (site, ligands, attributes) = parse_stereo_atom_addition(payload)?;
            Ok(EditInput::StereoAtomAdd {
                site,
                ligands,
                attributes,
            })
        }
        "modify" => {
            let (id, expect, update) = parse_stereo_atom_checked_update(payload)?;
            validate_stereo_atom_update_pair(&expect, &update)?;
            Ok(EditInput::StereoAtomModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown stereo-atom edit op :{other}"
        ))),
    }
}

fn parse_stereo_atoms_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-atoms edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown stereo-atoms edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of stereo-atom removals",
            got: payload.kind(),
            path: vec!["stereo-atoms edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_stereo_atom_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::StereoAtomsRemove(removes))
}

fn parse_stereo_bond_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-bond edit")?;
    match op {
        "add" => {
            let (site, ligands, attributes) = parse_stereo_bond_addition(payload)?;
            Ok(EditInput::StereoBondAdd {
                site,
                ligands,
                attributes,
            })
        }
        "modify" => {
            let (id, expect, update) = parse_stereo_bond_checked_update(payload)?;
            validate_stereo_bond_update_pair(&expect, &update)?;
            Ok(EditInput::StereoBondModify { id, expect, update })
        }
        other => Err(DeError::Custom(format!(
            "unknown stereo-bond edit op :{other}"
        ))),
    }
}

fn parse_stereo_bonds_edit(edn: &Edn<'_>) -> Result<EditInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-bonds edit")?;
    if op != "remove" {
        return Err(DeError::Custom(format!(
            "unknown stereo-bonds edit op :{op}"
        )));
    }
    let Edn::Vector(entries) = payload else {
        return Err(DeError::TypeMismatch {
            expected: "vector of stereo-bond removals",
            got: payload.kind(),
            path: vec!["stereo-bonds edit".to_string()],
        });
    };
    let removes = entries
        .iter()
        .map(parse_stereo_bond_removal)
        .collect::<Result<_, _>>()?;
    Ok(EditInput::StereoBondsRemove(removes))
}

fn parse_stereo_atom_addition(edn: &Edn<'_>) -> Result<StereoAtomAdditionInput, DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-atom addition map",
            got: edn.kind(),
            path: vec!["stereo-atom edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let site = helper.required("site")?;
    let ligands: Vec<Edn<'_>> = helper.required("ligands")?;
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((site, parse_stereo_ligands(&ligands)?, attributes))
}

fn parse_stereo_atom_removal(edn: &Edn<'_>) -> Result<StereoAtomRemovalInput, DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-atom removal map",
            got: edn.kind(),
            path: vec!["stereo-atoms edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let site = helper.required("site")?;
    let ligands: Vec<Edn<'_>> = helper.required("ligands")?;
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((id, site, parse_stereo_ligands(&ligands)?, attributes))
}

fn parse_stereo_bond_addition(edn: &Edn<'_>) -> Result<StereoBondAdditionInput, DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-bond addition map",
            got: edn.kind(),
            path: vec!["stereo-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let site = helper.required("site")?;
    let ligands: Vec<Edn<'_>> = helper.required("ligands")?;
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((site, parse_stereo_ligands(&ligands)?, attributes))
}

fn parse_stereo_bond_removal(edn: &Edn<'_>) -> Result<StereoBondRemovalInput, DeError> {
    let Edn::Map(map) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-bond removal map",
            got: edn.kind(),
            path: vec!["stereo-bonds edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(map);
    let id = helper.required("id")?;
    let site = helper.required("site")?;
    let ligands: Vec<Edn<'_>> = helper.required("ligands")?;
    let attributes = helper.required("attrs")?;
    helper.finalize()?;
    Ok((id, site, parse_stereo_ligands(&ligands)?, attributes))
}

fn parse_stereo_ligands(ligands: &[Edn<'_>]) -> Result<Vec<StereoLigandInput>, DeError> {
    ligands.iter().map(parse_stereo_ligand).collect()
}

fn parse_stereo_ligand(edn: &Edn<'_>) -> Result<StereoLigandInput, DeError> {
    match edn {
        Edn::Vector(parts) if parts.len() == 2 => {
            let Edn::Keyword(tag) = &parts[0] else {
                return Err(DeError::TypeMismatch {
                    expected: "stereo ligand kind keyword",
                    got: parts[0].kind(),
                    path: vec!["stereo ligand".to_string()],
                });
            };
            let kind = match tag.name() {
                "h" => StereoLigandKind::ImplicitHydrogen,
                "lp" => StereoLigandKind::LonePair,
                other => {
                    return Err(DeError::Custom(format!(
                        "unknown stereo ligand kind :{other}"
                    )));
                }
            };
            Ok((AtomHandle::from_edn(&parts[1])?, kind))
        }
        Edn::Vector(_) => Err(DeError::Custom(
            "stereo ligand vector expects [kind atom-handle]".to_string(),
        )),
        _ => Ok((AtomHandle::from_edn(edn)?, StereoLigandKind::Atom)),
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
    let constraint = ConstraintEditDsl::from_edn(payload)?.into_edit();
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
    Ok((
        AtomHandle::from_edn(&parts[0])?,
        expect.into_ir(&()),
        update.into_ir(&()),
    ))
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
    Ok((
        BondHandle::from_edn(&parts[0])?,
        expect.into_ir(&()),
        update.into_ir(&()),
    ))
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
    Ok((
        DativeBondHandle::from_edn(&parts[0])?,
        expect.into_ir(&()),
        update.into_ir(&()),
    ))
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
        expect.into_ir(&()),
        update.into_ir(&()),
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
        expect.into_ir(&()),
        update.into_ir(&()),
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
        expect.into_ir(&()),
        update.into_ir(&()),
    ))
}

fn parse_stereo_atom_checked_update(
    edn: &Edn<'_>,
) -> Result<(StereoAtomHandle, StereoAtomUpdate, StereoAtomUpdate), DeError> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-atom :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["stereo-atom edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "stereo-atom :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-atom :modify changes map",
            got: parts[1].kind(),
            path: vec!["stereo-atom edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: StereoAtomUpdateDsl = helper.required("expect")?;
    let update: StereoAtomUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((
        StereoAtomHandle::from_edn(&parts[0])?,
        expect.into_ir(&()),
        update.into_ir(&()),
    ))
}

fn parse_stereo_bond_checked_update(
    edn: &Edn<'_>,
) -> Result<(StereoBondHandle, StereoBondUpdate, StereoBondUpdate), DeError> {
    let Edn::Vector(parts) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-bond :modify [handle {:expect dsl :update dsl}]",
            got: edn.kind(),
            path: vec!["stereo-bond edit".to_string()],
        });
    };
    if parts.len() != 2 {
        return Err(DeError::Custom(format!(
            "stereo-bond :modify expects [handle changes], got {} elements",
            parts.len()
        )));
    }
    let Edn::Map(changes) = &parts[1] else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-bond :modify changes map",
            got: parts[1].kind(),
            path: vec!["stereo-bond edit".to_string()],
        });
    };
    let mut helper = EdnMapHelper::new(changes);
    let expect: StereoBondUpdateDsl = helper.required("expect")?;
    let update: StereoBondUpdateDsl = helper.required("update")?;
    helper.finalize()?;
    Ok((
        StereoBondHandle::from_edn(&parts[0])?,
        expect.into_ir(&()),
        update.into_ir(&()),
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
        .map(AtomConstraintForm::key)
        .eq(update.constraints.iter().map(AtomConstraintForm::key));
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
        .map(BondConstraintForm::key)
        .eq(update.constraints.iter().map(BondConstraintForm::key));
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
        .map(DativeBondConstraintForm::key)
        .eq(update.constraints.iter().map(DativeBondConstraintForm::key));
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
        .map(AromaticSystemConstraintForm::key)
        .eq(update
            .constraints
            .iter()
            .map(AromaticSystemConstraintForm::key));
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
        .map(MulticenterBondConstraintForm::key)
        .eq(update
            .constraints
            .iter()
            .map(MulticenterBondConstraintForm::key));
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
        .map(NoncovalentBondConstraintForm::key)
        .eq(update
            .constraints
            .iter()
            .map(NoncovalentBondConstraintForm::key));
    if !fields_match || !constraints_match {
        return Err(DeError::Custom(
            "noncovalent-bond :modify :expect and :update must address the same fields and constraints"
                .to_string(),
        ));
    }
    Ok(())
}

fn validate_stereo_atom_update_pair(
    expect: &StereoAtomUpdate,
    update: &StereoAtomUpdate,
) -> Result<(), DeError> {
    let expect_changes_configuration = !matches!(
        expect.configuration,
        StereoConfigurationUpdate::Unchanged
            | StereoConfigurationUpdate::Kinded { coset: None, .. }
    );
    let update_changes_configuration = !matches!(
        update.configuration,
        StereoConfigurationUpdate::Unchanged
            | StereoConfigurationUpdate::Kinded { coset: None, .. }
    );
    let constraints_match = expect
        .constraints
        .iter()
        .map(StereoAtomConstraintForm::key)
        .eq(update.constraints.iter().map(StereoAtomConstraintForm::key));
    if expect_changes_configuration != update_changes_configuration || !constraints_match {
        return Err(DeError::Custom(
            "stereo-atom :modify :expect and :update must address the same field and constraints"
                .to_string(),
        ));
    }
    if !expect.constraints.is_empty() {
        let (Some(expect_kind), Some(update_kind)) =
            (expect.configuration.kind(), update.configuration.kind())
        else {
            return Err(DeError::Custom(
                "stereo-atom constraint changes require a stereo kind in both :expect and :update"
                    .to_string(),
            ));
        };
        if expect_kind != update_kind {
            return Err(DeError::Custom(
                "stereo-atom constraint changes require the same stereo kind in :expect and :update"
                    .to_string(),
            ));
        }
    }
    Ok(())
}

fn validate_stereo_bond_update_pair(
    expect: &StereoBondUpdate,
    update: &StereoBondUpdate,
) -> Result<(), DeError> {
    let expect_changes_configuration = !matches!(
        expect.configuration,
        StereoConfigurationUpdate::Unchanged
            | StereoConfigurationUpdate::Kinded { coset: None, .. }
    );
    let update_changes_configuration = !matches!(
        update.configuration,
        StereoConfigurationUpdate::Unchanged
            | StereoConfigurationUpdate::Kinded { coset: None, .. }
    );
    let constraints_match = expect
        .constraints
        .iter()
        .map(StereoBondConstraintForm::key)
        .eq(update.constraints.iter().map(StereoBondConstraintForm::key));
    if expect_changes_configuration != update_changes_configuration || !constraints_match {
        return Err(DeError::Custom(
            "stereo-bond :modify :expect and :update must address the same field and constraints"
                .to_string(),
        ));
    }
    if !expect.constraints.is_empty() {
        let (Some(expect_kind), Some(update_kind)) =
            (expect.configuration.kind(), update.configuration.kind())
        else {
            return Err(DeError::Custom(
                "stereo-bond constraint changes require a stereo kind in both :expect and :update"
                    .to_string(),
            ));
        };
        if expect_kind != update_kind {
            return Err(DeError::Custom(
                "stereo-bond constraint changes require the same stereo kind in :expect and :update"
                    .to_string(),
            ));
        }
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

fn append_stereo_atom_modify(
    edits: &mut Edits,
    id: StereoAtomHandle,
    expect: StereoAtomUpdate,
    update: StereoAtomUpdate,
) -> Result<(), DeError> {
    validate_stereo_atom_update_pair(&expect, &update)?;
    let kind = expect.configuration.kind();
    let old_configuration = match expect.configuration {
        StereoConfigurationUpdate::Unchanged
        | StereoConfigurationUpdate::Kinded { coset: None, .. } => None,
        StereoConfigurationUpdate::Undetermined => Some(StereoConfigurationForm::Undetermined),
        StereoConfigurationUpdate::Kinded {
            kind,
            coset: Some(coset),
        } => Some(StereoConfigurationForm::kinded(kind, coset)),
    };
    let new_configuration = match update.configuration {
        StereoConfigurationUpdate::Unchanged
        | StereoConfigurationUpdate::Kinded { coset: None, .. } => None,
        StereoConfigurationUpdate::Undetermined => Some(StereoConfigurationForm::Undetermined),
        StereoConfigurationUpdate::Kinded {
            kind,
            coset: Some(coset),
        } => Some(StereoConfigurationForm::kinded(kind, coset)),
    };
    if let (Some(old), Some(new)) = (old_configuration, new_configuration) {
        edits.push(Edit::ModifyStereoAtomField {
            id: id.clone(),
            change: StereoAtomFieldChange::Configuration { old, new },
        });
    }
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyStereoAtomConstraint {
            id: id.clone(),
            kind,
            old: (!old.is_undetermined()).then(|| old.clone()),
            new: (!new.is_undetermined()).then(|| new.clone()),
        });
    }
    Ok(())
}

fn append_stereo_bond_modify(
    edits: &mut Edits,
    id: StereoBondHandle,
    expect: StereoBondUpdate,
    update: StereoBondUpdate,
) -> Result<(), DeError> {
    validate_stereo_bond_update_pair(&expect, &update)?;
    let kind = expect.configuration.kind();
    let old_configuration = match expect.configuration {
        StereoConfigurationUpdate::Unchanged
        | StereoConfigurationUpdate::Kinded { coset: None, .. } => None,
        StereoConfigurationUpdate::Undetermined => Some(StereoConfigurationForm::Undetermined),
        StereoConfigurationUpdate::Kinded {
            kind,
            coset: Some(coset),
        } => Some(StereoConfigurationForm::kinded(kind, coset)),
    };
    let new_configuration = match update.configuration {
        StereoConfigurationUpdate::Unchanged
        | StereoConfigurationUpdate::Kinded { coset: None, .. } => None,
        StereoConfigurationUpdate::Undetermined => Some(StereoConfigurationForm::Undetermined),
        StereoConfigurationUpdate::Kinded {
            kind,
            coset: Some(coset),
        } => Some(StereoConfigurationForm::kinded(kind, coset)),
    };
    if let (Some(old), Some(new)) = (old_configuration, new_configuration) {
        edits.push(Edit::ModifyStereoBondField {
            id: id.clone(),
            change: StereoBondFieldChange::Configuration { old, new },
        });
    }
    for (old, new) in expect.constraints.iter().zip(update.constraints.iter()) {
        edits.push(Edit::ModifyStereoBondConstraint {
            id: id.clone(),
            kind,
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
                old: UnpairedElectronsForm {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsForm {
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
                old: UnpairedElectronsForm {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsForm {
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
                old: UnpairedElectronsForm {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsForm {
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
                old: UnpairedElectronsForm {
                    count: old_count,
                    multiplicity: old_multiplicity,
                },
                new: UnpairedElectronsForm {
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

fn stereo_atom_field_updates(
    change: &StereoAtomFieldChange,
) -> (StereoAtomUpdate, StereoAtomUpdate) {
    let StereoAtomFieldChange::Configuration { old, new } = change;
    let expect = StereoAtomUpdate {
        configuration: match old {
            StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
            StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                kind: *kind,
                coset: Some(coset.clone()),
            },
        },
        ..Default::default()
    };
    let update = StereoAtomUpdate {
        configuration: match new {
            StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
            StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                kind: *kind,
                coset: Some(coset.clone()),
            },
        },
        ..Default::default()
    };
    (expect, update)
}

fn stereo_bond_field_updates(
    change: &StereoBondFieldChange,
) -> (StereoBondUpdate, StereoBondUpdate) {
    let StereoBondFieldChange::Configuration { old, new } = change;
    let expect = StereoBondUpdate {
        configuration: match old {
            StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
            StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                kind: *kind,
                coset: Some(coset.clone()),
            },
        },
        ..Default::default()
    };
    let update = StereoBondUpdate {
        configuration: match new {
            StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
            StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                kind: *kind,
                coset: Some(coset.clone()),
            },
        },
        ..Default::default()
    };
    (expect, update)
}

fn atom_constraint_updates(
    old: &Option<AtomConstraintForm>,
    new: &Option<AtomConstraintForm>,
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
    old: &Option<BondConstraintForm>,
    new: &Option<BondConstraintForm>,
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
    old: &Option<DativeBondConstraintForm>,
    new: &Option<DativeBondConstraintForm>,
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
    old: &Option<AromaticSystemConstraintForm>,
    new: &Option<AromaticSystemConstraintForm>,
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
    old: &Option<MulticenterBondConstraintForm>,
    new: &Option<MulticenterBondConstraintForm>,
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
    old: &Option<NoncovalentBondConstraintForm>,
    new: &Option<NoncovalentBondConstraintForm>,
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

fn stereo_atom_constraint_updates(
    kind: Option<StereoKind>,
    old: &Option<StereoAtomConstraintForm>,
    new: &Option<StereoAtomConstraintForm>,
) -> Result<(StereoAtomUpdate, StereoAtomUpdate), DeError> {
    let Some(kind) = kind else {
        return Err(DeError::Custom(
            "stereo-atom constraint edit requires a stereo kind".to_string(),
        ));
    };
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "stereo-atom constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = StereoAtomUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind, coset: None },
        ..Default::default()
    };
    let mut update = StereoAtomUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind, coset: None },
        ..Default::default()
    };
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

fn stereo_bond_constraint_updates(
    kind: Option<StereoKind>,
    old: &Option<StereoBondConstraintForm>,
    new: &Option<StereoBondConstraintForm>,
) -> Result<(StereoBondUpdate, StereoBondUpdate), DeError> {
    let Some(kind) = kind else {
        return Err(DeError::Custom(
            "stereo-bond constraint edit requires a stereo kind".to_string(),
        ));
    };
    let key_matches = match (old, new) {
        (Some(old), Some(new)) => old.key() == new.key(),
        (Some(_), None) | (None, Some(_)) => true,
        (None, None) => false,
    };
    if !key_matches {
        return Err(DeError::Custom(
            "stereo-bond constraint edit must address one constraint key".to_string(),
        ));
    }
    let mut expect = StereoBondUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind, coset: None },
        ..Default::default()
    };
    let mut update = StereoBondUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind, coset: None },
        ..Default::default()
    };
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
    attributes_edn: Edn<'static>,
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
    entry.insert(Edn::keyword("attrs"), attributes_edn);
    Edn::Map(entry)
}

fn relation_entry_edn(
    id: Option<Edn<'static>>,
    atoms: &[AtomHandle],
    attributes_edn: Edn<'static>,
) -> Edn<'static> {
    let mut entry = EdnMap::with_capacity(3);
    if let Some(id) = id {
        entry.insert(Edn::keyword("id"), id);
    }
    entry.insert(
        Edn::keyword("atoms"),
        Edn::Vector(atoms.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into()),
    );
    entry.insert(Edn::keyword("attrs"), attributes_edn);
    Edn::Map(entry)
}

fn stereo_entry_edn(
    id: Option<Edn<'static>>,
    site: Edn<'static>,
    ligands: &[(AtomHandle, StereoLigandKind)],
    attributes_edn: Edn<'static>,
) -> Edn<'static> {
    let mut entry = EdnMap::with_capacity(4);
    if let Some(id) = id {
        entry.insert(Edn::keyword("id"), id);
    }
    entry.insert(Edn::keyword("site"), site);
    entry.insert(
        Edn::keyword("ligands"),
        Edn::Vector(
            ligands
                .iter()
                .map(stereo_ligand_edn)
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    entry.insert(Edn::keyword("attrs"), attributes_edn);
    Edn::Map(entry)
}

fn stereo_ligand_edn(ligand: &(AtomHandle, StereoLigandKind)) -> Edn<'static> {
    let (atom, kind) = ligand;
    match kind {
        StereoLigandKind::Atom => atom.to_edn(),
        StereoLigandKind::ImplicitHydrogen => {
            Edn::Vector(vec![Edn::keyword("h"), atom.to_edn()].into())
        }
        StereoLigandKind::LonePair => Edn::Vector(vec![Edn::keyword("lp"), atom.to_edn()].into()),
    }
}

fn edit_map(entity: &str, operation: &str, payload: Edn<'static>) -> Edn<'static> {
    single_key_map(entity, single_key_map(operation, payload))
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ConstraintEditDsl {
    constraint: Constraint,
    handles: ConstraintHandles,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct ConstraintHandles {
    atoms: Vec<AtomHandle>,
    bonds: Vec<BondHandle>,
    dative_bonds: Vec<DativeBondHandle>,
    aromatic_systems: Vec<AromaticSystemHandle>,
    multicenter_bonds: Vec<MulticenterBondHandle>,
    noncovalent_bonds: Vec<NoncovalentBondHandle>,
    stereo_atoms: Vec<StereoAtomHandle>,
    stereo_bonds: Vec<StereoBondHandle>,
}

impl<'de> FromEdn<'de> for ConstraintEditDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let mut handles = ConstraintHandles::default();
        let constraint = parse_constraint(edn, &mut handles)?;
        Ok(Self {
            constraint,
            handles,
        })
    }
}

impl ToEdn for ConstraintEditDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_constraint(&self.constraint, &self.handles)
    }
}

impl ConstraintEditDsl {
    pub(super) fn from_edit(edit: ConstraintEdit) -> Self {
        let mut handles = ConstraintHandles::default();
        let ConstraintHandles {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        } = &mut handles;
        let constraint = edit
            .resolve(
                |handle| Ok::<_, Infallible>(AtomId::from(intern(atoms, handle))),
                |handle| Ok::<_, Infallible>(BondId::from(intern(bonds, handle))),
                |handle| Ok::<_, Infallible>(DativeBondId::from(intern(dative_bonds, handle))),
                |handle| {
                    Ok::<_, Infallible>(AromaticSystemId::from(intern(aromatic_systems, handle)))
                },
                |handle| {
                    Ok::<_, Infallible>(MulticenterBondId::from(intern(multicenter_bonds, handle)))
                },
                |handle| {
                    Ok::<_, Infallible>(NoncovalentBondId::from(intern(noncovalent_bonds, handle)))
                },
                |handle| Ok::<_, Infallible>(StereoAtomId::from(intern(stereo_atoms, handle))),
                |handle| Ok::<_, Infallible>(StereoBondId::from(intern(stereo_bonds, handle))),
            )
            .expect("normalizing constraint-edit handles is infallible");
        Self {
            constraint,
            handles,
        }
    }

    pub(super) fn into_edit(self) -> ConstraintEdit {
        let Self {
            constraint,
            handles,
        } = self;
        ConstraintEdit::new(constraint, |entity| handles.get(entity))
            .expect("constraint-edit DSL indices have kind-correct handles")
    }
}

impl ConstraintHandles {
    fn atom(&mut self, handle: AtomHandle) -> AtomId {
        AtomId::from(intern(&mut self.atoms, handle))
    }

    fn bond(&mut self, handle: BondHandle) -> BondId {
        BondId::from(intern(&mut self.bonds, handle))
    }

    fn dative_bond(&mut self, handle: DativeBondHandle) -> DativeBondId {
        DativeBondId::from(intern(&mut self.dative_bonds, handle))
    }

    fn aromatic_system(&mut self, handle: AromaticSystemHandle) -> AromaticSystemId {
        AromaticSystemId::from(intern(&mut self.aromatic_systems, handle))
    }

    fn multicenter_bond(&mut self, handle: MulticenterBondHandle) -> MulticenterBondId {
        MulticenterBondId::from(intern(&mut self.multicenter_bonds, handle))
    }

    fn noncovalent_bond(&mut self, handle: NoncovalentBondHandle) -> NoncovalentBondId {
        NoncovalentBondId::from(intern(&mut self.noncovalent_bonds, handle))
    }

    fn stereo_atom(&mut self, handle: StereoAtomHandle) -> StereoAtomId {
        StereoAtomId::from(intern(&mut self.stereo_atoms, handle))
    }

    fn stereo_bond(&mut self, handle: StereoBondHandle) -> StereoBondId {
        StereoBondId::from(intern(&mut self.stereo_bonds, handle))
    }

    fn get(&self, entity: Entity) -> Option<EntityHandle> {
        match entity {
            Entity::Atom(id) => self.atoms.get(id.index()).cloned().map(EntityHandle::Atom),
            Entity::Bond(id) => self.bonds.get(id.index()).cloned().map(EntityHandle::Bond),
            Entity::DativeBond(id) => self
                .dative_bonds
                .get(id.index())
                .cloned()
                .map(EntityHandle::DativeBond),
            Entity::AromaticSystem(id) => self
                .aromatic_systems
                .get(id.index())
                .cloned()
                .map(EntityHandle::AromaticSystem),
            Entity::MulticenterBond(id) => self
                .multicenter_bonds
                .get(id.index())
                .cloned()
                .map(EntityHandle::MulticenterBond),
            Entity::NoncovalentBond(id) => self
                .noncovalent_bonds
                .get(id.index())
                .cloned()
                .map(EntityHandle::NoncovalentBond),
            Entity::StereoAtom(id) => self
                .stereo_atoms
                .get(id.index())
                .cloned()
                .map(EntityHandle::StereoAtom),
            Entity::StereoBond(id) => self
                .stereo_bonds
                .get(id.index())
                .cloned()
                .map(EntityHandle::StereoBond),
        }
    }
}

fn intern<H: Eq>(handles: &mut Vec<H>, handle: H) -> usize {
    handles
        .iter()
        .position(|candidate| candidate == &handle)
        .unwrap_or_else(|| {
            handles.push(handle);
            handles.len() - 1
        })
}

fn parse_constraint(edn: &Edn<'_>, handles: &mut ConstraintHandles) -> Result<Constraint, DeError> {
    let (key, payload) = parse_single_key_map(edn, "constraint")?;
    Ok(match key {
        "atom" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::Atom(
                handles.atom(AtomHandle::from_edn(handle)?),
                AtomConstraintDsl::from_edn(constraint)?.into_ir(&()),
            )
        }
        "bond" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::Bond(
                handles.bond(BondHandle::from_edn(handle)?),
                BondConstraintDsl::from_edn(constraint)?.into_ir(&()),
            )
        }
        "dative-bond" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::DativeBond(
                handles.dative_bond(DativeBondHandle::from_edn(handle)?),
                DativeBondConstraintDsl::from_edn(constraint)?.into_ir(),
            )
        }
        "aromatic-system" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::AromaticSystem(
                handles.aromatic_system(AromaticSystemHandle::from_edn(handle)?),
                AromaticSystemConstraintDsl::from_edn(constraint)?.into_ir(),
            )
        }
        "multicenter-bond" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::MulticenterBond(
                handles.multicenter_bond(MulticenterBondHandle::from_edn(handle)?),
                MulticenterBondConstraintDsl::from_edn(constraint)?.into_ir(),
            )
        }
        "noncovalent-bond" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            Constraint::NoncovalentBond(
                handles.noncovalent_bond(NoncovalentBondHandle::from_edn(handle)?),
                NoncovalentBondConstraintDsl::from_edn(constraint)?.into_ir(),
            )
        }
        "stereo-atom" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            let StereoAtomConstraintDsl(kind, constraint) =
                StereoAtomConstraintDsl::from_edn(constraint)?;
            Constraint::StereoAtom(
                handles.stereo_atom(StereoAtomHandle::from_edn(handle)?),
                kind,
                constraint,
            )
        }
        "stereo-bond" => {
            let (handle, constraint) = parse_pair(payload, key)?;
            let StereoBondConstraintDsl(kind, constraint) =
                StereoBondConstraintDsl::from_edn(constraint)?;
            Constraint::StereoBond(
                handles.stereo_bond(StereoBondHandle::from_edn(handle)?),
                kind,
                constraint,
            )
        }
        "and" => Constraint::And(parse_constraint_vec(payload, handles)?),
        "or" => Constraint::Or(parse_constraint_vec(payload, handles)?),
        "not" => Constraint::Not(Box::new(parse_constraint(payload, handles)?)),
        "charge-sum" | "unpaired-electron-coupling" | "bond-order-sum" | "connected" => {
            Constraint::Molecule(parse_molecule_constraint(key, payload, handles)?)
        }
        _ => Constraint::Relational(parse_relational_constraint(key, payload, handles)?),
    })
}

fn render_constraint(constraint: &Constraint, handles: &ConstraintHandles) -> Edn<'static> {
    match constraint {
        Constraint::Atom(id, constraint) => entity_leaf_edn(
            "atom",
            handles.atoms[id.index()].to_edn(),
            AtomConstraintDsl::from_ir(constraint, &()).to_edn(),
        ),
        Constraint::Bond(id, constraint) => entity_leaf_edn(
            "bond",
            handles.bonds[id.index()].to_edn(),
            BondConstraintDsl::from_ir(constraint, &()).to_edn(),
        ),
        Constraint::DativeBond(id, constraint) => entity_leaf_edn(
            "dative-bond",
            handles.dative_bonds[id.index()].to_edn(),
            DativeBondConstraintDsl::from_ir(constraint).to_edn(),
        ),
        Constraint::AromaticSystem(id, constraint) => entity_leaf_edn(
            "aromatic-system",
            handles.aromatic_systems[id.index()].to_edn(),
            AromaticSystemConstraintDsl::from_ir(constraint).to_edn(),
        ),
        Constraint::MulticenterBond(id, constraint) => entity_leaf_edn(
            "multicenter-bond",
            handles.multicenter_bonds[id.index()].to_edn(),
            MulticenterBondConstraintDsl::from_ir(constraint).to_edn(),
        ),
        Constraint::NoncovalentBond(id, constraint) => entity_leaf_edn(
            "noncovalent-bond",
            handles.noncovalent_bonds[id.index()].to_edn(),
            NoncovalentBondConstraintDsl::from_ir(constraint).to_edn(),
        ),
        Constraint::StereoAtom(id, kind, constraint) => entity_leaf_edn(
            "stereo-atom",
            handles.stereo_atoms[id.index()].to_edn(),
            StereoAtomConstraintDsl::from_ir(constraint, kind).to_edn(),
        ),
        Constraint::StereoBond(id, kind, constraint) => entity_leaf_edn(
            "stereo-bond",
            handles.stereo_bonds[id.index()].to_edn(),
            StereoBondConstraintDsl::from_ir(constraint, kind).to_edn(),
        ),
        Constraint::Relational(constraint) => render_relational_constraint(constraint, handles),
        Constraint::Molecule(constraint) => render_molecule_constraint(constraint, handles),
        Constraint::And(constraints) => combinator_edn("and", constraints, handles),
        Constraint::Or(constraints) => combinator_edn("or", constraints, handles),
        Constraint::Not(constraint) => {
            single_key_map("not", render_constraint(constraint, handles))
        }
    }
}

fn parse_constraint_vec(
    edn: &Edn<'_>,
    handles: &mut ConstraintHandles,
) -> Result<Vec<Constraint>, DeError> {
    let Edn::Vector(values) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of constraints",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    values
        .iter()
        .map(|constraint| parse_constraint(constraint, handles))
        .collect()
}

fn combinator_edn(
    key: &str,
    constraints: &[Constraint],
    handles: &ConstraintHandles,
) -> Edn<'static> {
    single_key_map(
        key,
        Edn::Vector(
            constraints
                .iter()
                .map(|constraint| render_constraint(constraint, handles))
                .collect::<Vec<_>>()
                .into(),
        ),
    )
}

fn entity_leaf_edn(key: &str, handle: Edn<'static>, constraint: Edn<'static>) -> Edn<'static> {
    single_key_map(key, Edn::Vector(vec![handle, constraint].into()))
}

fn parse_pair<'a>(edn: &'a Edn<'_>, context: &str) -> Result<(&'a Edn<'a>, &'a Edn<'a>), DeError> {
    let Edn::Vector(values) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "2-element vector",
            got: edn.kind(),
            path: vec![context.to_string()],
        });
    };
    if values.len() != 2 {
        return Err(DeError::Custom(format!(
            "{context}: expected 2 elements, got {}",
            values.len()
        )));
    }
    Ok((&values[0], &values[1]))
}

fn parse_handle_vec<H>(edn: &Edn<'_>, context: &str) -> Result<Vec<H>, DeError>
where
    H: for<'de> FromEdn<'de>,
{
    let Edn::Vector(values) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of edit handles",
            got: edn.kind(),
            path: vec![context.to_string()],
        });
    };
    values.iter().map(H::from_edn).collect()
}

fn render_handle_vec<H: ToEdn>(handles: &[H]) -> Edn<'static> {
    Edn::Vector(handles.iter().map(ToEdn::to_edn).collect::<Vec<_>>().into())
}

fn parse_molecule_constraint(
    key: &str,
    payload: &Edn<'_>,
    handles: &mut ConstraintHandles,
) -> Result<MoleculeConstraint, DeError> {
    let map = expect_map(payload, "molecule constraint")?;
    Ok(match key {
        "charge-sum" => {
            let mut helper = EdnMapHelper::new(map);
            let atoms = helper
                .optional::<Vec<AtomHandle>>("atoms")?
                .map(|atoms| atoms.into_iter().map(|atom| handles.atom(atom)).collect());
            let sum: NumDsl = helper.required("sum")?;
            helper.finalize()?;
            MoleculeConstraint::ChargeSum {
                atoms,
                sum: sum.into_ir(&()),
            }
        }
        "unpaired-electron-coupling" => {
            let mut helper = EdnMapHelper::new(map);
            let atoms = helper
                .optional::<Vec<AtomHandle>>("atoms")?
                .map(|atoms| atoms.into_iter().map(|atom| handles.atom(atom)).collect());
            let unpaired_electrons = parse_unpaired_electrons(
                map.get_keyword("unpaired-electrons")
                    .ok_or_else(|| DeError::MissingField {
                        key: "unpaired-electrons".to_string(),
                        path: vec!["unpaired-electron-coupling".to_string()],
                    })?,
            )?;
            helper.optional::<Edn<'_>>("unpaired-electrons")?;
            helper.finalize()?;
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            }
        }
        "bond-order-sum" => {
            let mut helper = EdnMapHelper::new(map);
            let bonds = helper
                .optional::<Vec<BondHandle>>("bonds")?
                .map(|bonds| bonds.into_iter().map(|bond| handles.bond(bond)).collect());
            let sum: NumDsl = helper.required("sum")?;
            helper.finalize()?;
            MoleculeConstraint::BondOrderSum {
                bonds,
                sum: sum.into_ir(&()),
            }
        }
        "connected" => {
            let mut helper = EdnMapHelper::new(map);
            let atoms = helper
                .optional::<Vec<AtomHandle>>("atoms")?
                .map(|atoms| atoms.into_iter().map(|atom| handles.atom(atom)).collect());
            helper.finalize()?;
            MoleculeConstraint::Connected { atoms }
        }
        _ => unreachable!("molecule-constraint key was classified by parse_constraint"),
    })
}

fn render_molecule_constraint(
    constraint: &MoleculeConstraint,
    handles: &ConstraintHandles,
) -> Edn<'static> {
    match constraint {
        MoleculeConstraint::ChargeSum { atoms, sum } => {
            let mut map = EdnMap::with_capacity(2);
            if let Some(atoms) = atoms {
                map.insert(
                    Edn::keyword("atoms"),
                    render_handle_vec(
                        &atoms
                            .iter()
                            .map(|id| handles.atoms[id.index()].clone())
                            .collect::<Vec<_>>(),
                    ),
                );
            }
            map.insert(Edn::keyword("sum"), NumDsl::from_ir(sum, &()).to_edn());
            single_key_map("charge-sum", Edn::Map(map))
        }
        MoleculeConstraint::UnpairedElectronCoupling {
            atoms,
            unpaired_electrons,
        } => {
            let mut map = EdnMap::with_capacity(2);
            if let Some(atoms) = atoms {
                map.insert(
                    Edn::keyword("atoms"),
                    render_handle_vec(
                        &atoms
                            .iter()
                            .map(|id| handles.atoms[id.index()].clone())
                            .collect::<Vec<_>>(),
                    ),
                );
            }
            map.insert(
                Edn::keyword("unpaired-electrons"),
                render_unpaired_electrons(unpaired_electrons),
            );
            single_key_map("unpaired-electron-coupling", Edn::Map(map))
        }
        MoleculeConstraint::BondOrderSum { bonds, sum } => {
            let mut map = EdnMap::with_capacity(2);
            if let Some(bonds) = bonds {
                map.insert(
                    Edn::keyword("bonds"),
                    render_handle_vec(
                        &bonds
                            .iter()
                            .map(|id| handles.bonds[id.index()].clone())
                            .collect::<Vec<_>>(),
                    ),
                );
            }
            map.insert(Edn::keyword("sum"), NumDsl::from_ir(sum, &()).to_edn());
            single_key_map("bond-order-sum", Edn::Map(map))
        }
        MoleculeConstraint::Connected { atoms } => {
            let mut map = EdnMap::with_capacity(1);
            if let Some(atoms) = atoms {
                map.insert(
                    Edn::keyword("atoms"),
                    render_handle_vec(
                        &atoms
                            .iter()
                            .map(|id| handles.atoms[id.index()].clone())
                            .collect::<Vec<_>>(),
                    ),
                );
            }
            single_key_map("connected", Edn::Map(map))
        }
    }
}
fn parse_relational_constraint(
    key: &str,
    payload: &Edn<'_>,
    handles: &mut ConstraintHandles,
) -> Result<RelationalConstraint, DeError> {
    use RelationalConstraint as R;

    let constraint = match key {
        "dative-bond-donors" => {
            let (bond, atoms) = parse_pair(payload, key)?;
            R::DativeBondDonors {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "dative-bond-donor" => {
            let (bond, atom) = parse_pair(payload, key)?;
            R::DativeBondDonor {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "dative-bond-contains-all-donors" => {
            let (bond, atoms) = parse_pair(payload, key)?;
            R::DativeBondContainsAllDonors {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "dative-bond-all-donors" => {
            let (bond, predicate) = parse_pair(payload, key)?;
            R::DativeBondAllDonors {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "dative-bond-any-donor" => {
            let (bond, predicate) = parse_pair(payload, key)?;
            R::DativeBondAnyDonor {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "dative-bond-acceptor" => {
            let (bond, atom) = parse_pair(payload, key)?;
            R::DativeBondAcceptor {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "dative-bond-acceptor-satisfies" => {
            let (bond, predicate) = parse_pair(payload, key)?;
            R::DativeBondAcceptorSatisfies {
                bond: handles.dative_bond(DativeBondHandle::from_edn(bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "dative-bond-parallels" => {
            let (dative, parallel) = parse_pair(payload, key)?;
            R::DativeBondParallels {
                dative: handles.dative_bond(DativeBondHandle::from_edn(dative)?),
                parallel: handles.bond(BondHandle::from_edn(parallel)?),
            }
        }
        "aromatic-system-atoms" => {
            let (system, atoms) = parse_pair(payload, key)?;
            R::AromaticSystemAtoms {
                system: handles.aromatic_system(AromaticSystemHandle::from_edn(system)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "aromatic-system-contains" => {
            let (system, atom) = parse_pair(payload, key)?;
            R::AromaticSystemContains {
                system: handles.aromatic_system(AromaticSystemHandle::from_edn(system)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "aromatic-system-contains-all" => {
            let (system, atoms) = parse_pair(payload, key)?;
            R::AromaticSystemContainsAll {
                system: handles.aromatic_system(AromaticSystemHandle::from_edn(system)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "aromatic-system-all-atoms" => {
            let (system, predicate) = parse_pair(payload, key)?;
            R::AromaticSystemAllAtoms {
                system: handles.aromatic_system(AromaticSystemHandle::from_edn(system)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "aromatic-system-any-atom" => {
            let (system, predicate) = parse_pair(payload, key)?;
            R::AromaticSystemAnyAtom {
                system: handles.aromatic_system(AromaticSystemHandle::from_edn(system)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "multicenter-bond-atoms" => {
            let (bond, atoms) = parse_pair(payload, key)?;
            R::MulticenterBondAtoms {
                bond: handles.multicenter_bond(MulticenterBondHandle::from_edn(bond)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "multicenter-bond-contains" => {
            let (bond, atom) = parse_pair(payload, key)?;
            R::MulticenterBondContains {
                bond: handles.multicenter_bond(MulticenterBondHandle::from_edn(bond)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "multicenter-bond-contains-all" => {
            let (bond, atoms) = parse_pair(payload, key)?;
            R::MulticenterBondContainsAll {
                bond: handles.multicenter_bond(MulticenterBondHandle::from_edn(bond)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "multicenter-bond-all-atoms" => {
            let (bond, predicate) = parse_pair(payload, key)?;
            R::MulticenterBondAllAtoms {
                bond: handles.multicenter_bond(MulticenterBondHandle::from_edn(bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "multicenter-bond-any-atom" => {
            let (bond, predicate) = parse_pair(payload, key)?;
            R::MulticenterBondAnyAtom {
                bond: handles.multicenter_bond(MulticenterBondHandle::from_edn(bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "noncovalent-bond-ends" => {
            let (bond, atoms) = parse_pair(payload, key)?;
            let (first, second) = parse_pair(atoms, key)?;
            R::NoncovalentBondEnds {
                bond: handles.noncovalent_bond(NoncovalentBondHandle::from_edn(bond)?),
                atoms: [
                    handles.atom(AtomHandle::from_edn(first)?),
                    handles.atom(AtomHandle::from_edn(second)?),
                ],
            }
        }
        "noncovalent-bond-contains" => {
            let (bond, atom) = parse_pair(payload, key)?;
            R::NoncovalentBondContains {
                bond: handles.noncovalent_bond(NoncovalentBondHandle::from_edn(bond)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "noncovalent-bond-ends-satisfy" => {
            let (bond, predicates) = parse_pair(payload, key)?;
            let (first, second) = parse_pair(predicates, key)?;
            R::NoncovalentBondEndsSatisfy {
                bond: handles.noncovalent_bond(NoncovalentBondHandle::from_edn(bond)?),
                predicates: [
                    Box::new(AtomConstraintDsl::from_edn(first)?.into_ir(&())),
                    Box::new(AtomConstraintDsl::from_edn(second)?.into_ir(&())),
                ],
            }
        }
        "stereo-atom-site" => {
            let (stereo_atom, atom) = parse_pair(payload, key)?;
            R::StereoAtomSite {
                stereo_atom: handles.stereo_atom(StereoAtomHandle::from_edn(stereo_atom)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "stereo-atom-contains" => {
            let (stereo_atom, atom) = parse_pair(payload, key)?;
            R::StereoAtomContains {
                stereo_atom: handles.stereo_atom(StereoAtomHandle::from_edn(stereo_atom)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "stereo-atom-ligands" => {
            let (stereo_atom, atoms) = parse_pair(payload, key)?;
            R::StereoAtomLigands {
                stereo_atom: handles.stereo_atom(StereoAtomHandle::from_edn(stereo_atom)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "stereo-atom-all-ligands" => {
            let (stereo_atom, predicate) = parse_pair(payload, key)?;
            R::StereoAtomAllLigands {
                stereo_atom: handles.stereo_atom(StereoAtomHandle::from_edn(stereo_atom)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "stereo-atom-any-ligand" => {
            let (stereo_atom, predicate) = parse_pair(payload, key)?;
            R::StereoAtomAnyLigand {
                stereo_atom: handles.stereo_atom(StereoAtomHandle::from_edn(stereo_atom)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "stereo-bond-site" => {
            let (stereo_bond, bond) = parse_pair(payload, key)?;
            R::StereoBondSite {
                stereo_bond: handles.stereo_bond(StereoBondHandle::from_edn(stereo_bond)?),
                bond: handles.bond(BondHandle::from_edn(bond)?),
            }
        }
        "stereo-bond-contains" => {
            let (stereo_bond, atom) = parse_pair(payload, key)?;
            R::StereoBondContains {
                stereo_bond: handles.stereo_bond(StereoBondHandle::from_edn(stereo_bond)?),
                atom: handles.atom(AtomHandle::from_edn(atom)?),
            }
        }
        "stereo-bond-ligands" => {
            let (stereo_bond, atoms) = parse_pair(payload, key)?;
            R::StereoBondLigands {
                stereo_bond: handles.stereo_bond(StereoBondHandle::from_edn(stereo_bond)?),
                atoms: parse_handle_vec::<AtomHandle>(atoms, key)?
                    .into_iter()
                    .map(|atom| handles.atom(atom))
                    .collect(),
            }
        }
        "stereo-bond-all-ligands" => {
            let (stereo_bond, predicate) = parse_pair(payload, key)?;
            R::StereoBondAllLigands {
                stereo_bond: handles.stereo_bond(StereoBondHandle::from_edn(stereo_bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        "stereo-bond-any-ligand" => {
            let (stereo_bond, predicate) = parse_pair(payload, key)?;
            R::StereoBondAnyLigand {
                stereo_bond: handles.stereo_bond(StereoBondHandle::from_edn(stereo_bond)?),
                predicate: Box::new(AtomConstraintDsl::from_edn(predicate)?.into_ir(&())),
            }
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["constraint edit".to_string()],
            });
        }
    };
    Ok(constraint)
}

fn render_relational_constraint(
    constraint: &RelationalConstraint,
    handles: &ConstraintHandles,
) -> Edn<'static> {
    use RelationalConstraint as R;

    let (key, payload) = match constraint {
        R::DativeBondDonors { bond, atoms } => (
            "dative-bond-donors",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::DativeBondDonor { bond, atom } => (
            "dative-bond-donor",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::DativeBondContainsAllDonors { bond, atoms } => (
            "dative-bond-contains-all-donors",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::DativeBondAllDonors { bond, predicate } => (
            "dative-bond-all-donors",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::DativeBondAnyDonor { bond, predicate } => (
            "dative-bond-any-donor",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::DativeBondAcceptor { bond, atom } => (
            "dative-bond-acceptor",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::DativeBondAcceptorSatisfies { bond, predicate } => (
            "dative-bond-acceptor-satisfies",
            relation_pair(
                handles.dative_bonds[bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::DativeBondParallels { dative, parallel } => (
            "dative-bond-parallels",
            relation_pair(
                handles.dative_bonds[dative.index()].to_edn(),
                handles.bonds[parallel.index()].to_edn(),
            ),
        ),
        R::AromaticSystemAtoms { system, atoms } => (
            "aromatic-system-atoms",
            relation_pair(
                handles.aromatic_systems[system.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::AromaticSystemContains { system, atom } => (
            "aromatic-system-contains",
            relation_pair(
                handles.aromatic_systems[system.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::AromaticSystemContainsAll { system, atoms } => (
            "aromatic-system-contains-all",
            relation_pair(
                handles.aromatic_systems[system.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::AromaticSystemAllAtoms { system, predicate } => (
            "aromatic-system-all-atoms",
            relation_pair(
                handles.aromatic_systems[system.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::AromaticSystemAnyAtom { system, predicate } => (
            "aromatic-system-any-atom",
            relation_pair(
                handles.aromatic_systems[system.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::MulticenterBondAtoms { bond, atoms } => (
            "multicenter-bond-atoms",
            relation_pair(
                handles.multicenter_bonds[bond.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::MulticenterBondContains { bond, atom } => (
            "multicenter-bond-contains",
            relation_pair(
                handles.multicenter_bonds[bond.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::MulticenterBondContainsAll { bond, atoms } => (
            "multicenter-bond-contains-all",
            relation_pair(
                handles.multicenter_bonds[bond.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::MulticenterBondAllAtoms { bond, predicate } => (
            "multicenter-bond-all-atoms",
            relation_pair(
                handles.multicenter_bonds[bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::MulticenterBondAnyAtom { bond, predicate } => (
            "multicenter-bond-any-atom",
            relation_pair(
                handles.multicenter_bonds[bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::NoncovalentBondEnds { bond, atoms } => (
            "noncovalent-bond-ends",
            relation_pair(
                handles.noncovalent_bonds[bond.index()].to_edn(),
                Edn::Vector(
                    vec![
                        handles.atoms[atoms[0].index()].to_edn(),
                        handles.atoms[atoms[1].index()].to_edn(),
                    ]
                    .into(),
                ),
            ),
        ),
        R::NoncovalentBondContains { bond, atom } => (
            "noncovalent-bond-contains",
            relation_pair(
                handles.noncovalent_bonds[bond.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::NoncovalentBondEndsSatisfy { bond, predicates } => (
            "noncovalent-bond-ends-satisfy",
            relation_pair(
                handles.noncovalent_bonds[bond.index()].to_edn(),
                Edn::Vector(
                    vec![
                        AtomConstraintDsl::from_ir(&predicates[0], &()).to_edn(),
                        AtomConstraintDsl::from_ir(&predicates[1], &()).to_edn(),
                    ]
                    .into(),
                ),
            ),
        ),
        R::StereoAtomSite { stereo_atom, atom } => (
            "stereo-atom-site",
            relation_pair(
                handles.stereo_atoms[stereo_atom.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::StereoAtomContains { stereo_atom, atom } => (
            "stereo-atom-contains",
            relation_pair(
                handles.stereo_atoms[stereo_atom.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::StereoAtomLigands { stereo_atom, atoms } => (
            "stereo-atom-ligands",
            relation_pair(
                handles.stereo_atoms[stereo_atom.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::StereoAtomAllLigands {
            stereo_atom,
            predicate,
        } => (
            "stereo-atom-all-ligands",
            relation_pair(
                handles.stereo_atoms[stereo_atom.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::StereoAtomAnyLigand {
            stereo_atom,
            predicate,
        } => (
            "stereo-atom-any-ligand",
            relation_pair(
                handles.stereo_atoms[stereo_atom.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::StereoBondSite { stereo_bond, bond } => (
            "stereo-bond-site",
            relation_pair(
                handles.stereo_bonds[stereo_bond.index()].to_edn(),
                handles.bonds[bond.index()].to_edn(),
            ),
        ),
        R::StereoBondContains { stereo_bond, atom } => (
            "stereo-bond-contains",
            relation_pair(
                handles.stereo_bonds[stereo_bond.index()].to_edn(),
                handles.atoms[atom.index()].to_edn(),
            ),
        ),
        R::StereoBondLigands { stereo_bond, atoms } => (
            "stereo-bond-ligands",
            relation_pair(
                handles.stereo_bonds[stereo_bond.index()].to_edn(),
                render_atom_ids(atoms, handles),
            ),
        ),
        R::StereoBondAllLigands {
            stereo_bond,
            predicate,
        } => (
            "stereo-bond-all-ligands",
            relation_pair(
                handles.stereo_bonds[stereo_bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
        R::StereoBondAnyLigand {
            stereo_bond,
            predicate,
        } => (
            "stereo-bond-any-ligand",
            relation_pair(
                handles.stereo_bonds[stereo_bond.index()].to_edn(),
                AtomConstraintDsl::from_ir(predicate, &()).to_edn(),
            ),
        ),
    };
    single_key_map(key, payload)
}

fn relation_pair(owner: Edn<'static>, target: Edn<'static>) -> Edn<'static> {
    Edn::Vector(vec![owner, target].into())
}

fn render_atom_ids(atoms: &[AtomId], handles: &ConstraintHandles) -> Edn<'static> {
    Edn::Vector(
        atoms
            .iter()
            .map(|atom| handles.atoms[atom.index()].to_edn())
            .collect::<Vec<_>>()
            .into(),
    )
}
#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::{read_string, EdnError};

    use super::*;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::{AtomForm, ElementForm, IsotopeMassForm};
    use crate::ir::bond::BondForm;
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::{
        AromaticSystemConstraintForm, AtomConstraintsForm, BondConstraintsForm, Constraint,
        DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
        NoncovalentBondConstraintForm, RingMembershipForm, RingScope, StereoAtomConstraintForm,
        StereoBondConstraintForm, StereogenicityForm,
    };
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::AddBond;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::molecule::Molecule;
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
    };
    use crate::ir::num::NumForm;
    use crate::ir::stereo::{
        StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoKind, Stereogenicity,
    };
    use crate::mol_dsl;

    #[rstest]
    #[case::empty("[]", MoleculeDefaults::new(), Edits::new())]
    #[case::construction(
        r#"[{:atom {:add "C#h3"}} {:bond {:add [0 {:new 0} :single]}}]"#,
        MoleculeDefaults::new(),
        Edits::from_iter([
            Edit::AddAtoms {
                atoms: vec![AtomForm {
                    element: ElementForm::Lit(Element::C),
                    implicit_hydrogens: NumForm::Lit(3),
                    ..Default::default()
                }],
            },
            Edit::AddBonds {
                bonds: vec![AddBond {
                    endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                    attributes: BondForm::from_order(1),
                }],
            },
        ]),
    )]
    #[case::checked_update(
        r##"[{:atom {:modify [0 {:expect "#h3" :update "#h2"}]}}]"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::ImplicitHydrogens {
                old: NumForm::Lit(3),
                new: NumForm::Lit(2),
            },
        }]),
    )]
    #[case::ground_default(
        r#"[{:atom {:add "O"}}]"#,
        MoleculeDefaults::concrete(),
        Edits::from_iter([Edit::AddAtoms {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::O),
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(0),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                constraints: AtomConstraintsForm::new(),
            }],
        }]),
    )]
    fn test_edits_dsl_roundtrip(
        #[case] input: &str,
        #[case] defaults: MoleculeDefaults,
        #[case] expected: Edits,
    ) {
        let dsl = EditsDsl::from_str(input).unwrap();
        let rendered = dsl.to_edn();
        let displayed = dsl.to_string();
        let edits = dsl.into_ir(&defaults);
        let rebuilt = EditsDsl::from_ir(&expected, &defaults);

        assert_eq!(edits, expected);
        assert_eq!(rendered, read_string(input).unwrap());
        assert_eq!(displayed, rendered.to_string());
        assert_eq!(rebuilt.to_edn(), read_string(input).unwrap());
    }

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
        MoleculeDefaults::concrete(),
        Edits::from_iter([Edit::AddAtoms {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::C),
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(0),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                constraints: AtomConstraintsForm::new(),
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
                old: NumForm::Lit(0),
                new: NumForm::Lit(-1),
            },
        }]),
    )]
    #[case::atom_constraint_add(
        r##"{:atom {:modify [0 {:expect "#v*" :update "#v4"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintForm::valence(4_i64)),
        }]),
    )]
    #[case::atom_constraint_remove(
        r##"{:atom {:modify [0 {:expect "#v4" :update "#v*"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: Some(AtomConstraintForm::valence(4_i64)),
            new: None,
        }]),
    )]
    #[case::bond_add(
        "{:bond {:add [0 {:new 0} :single]}}",
        MoleculeDefaults::concrete(),
        Edits::from_iter([Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                attributes: BondForm {
                    order: NumForm::Lit(1),
                    charge: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                    constraints: BondConstraintsForm::new(),
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
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }]),
    )]
    #[case::bond_constraint_add(
        r##"{:bond {:modify [0 {:expect "#R(6)*" :update "#R(6)"}]}}"##,
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: None,
            new: Some(BondConstraintForm::RingMembership(RingMembershipForm::new(
                RingScope::Size(6),
                NumForm::Lit(1),
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
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }).into(),
        }]),
    )]
    #[case::constraint_remove(
        "{:constraint {:remove {:connected {}}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::RemoveMoleculeConstraint {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }).into(),
        }]),
    )]
    #[case::constraint_positional(
        "{:constraint {:add {:atom [2 {:valence 4}]}}}",
        MoleculeDefaults::new(),
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::Atom(AtomId(2), AtomConstraintForm::valence(4_i64)).into(),
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
            atoms: vec![AtomForm::from_element(Element::C).into_concrete()],
        },
        MoleculeDefaults::concrete(),
        r#"{:atom {:add "C"}}"#,
    )]
    #[case::bond_add(
        Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                attributes: BondForm::from_order(1).into_concrete(),
            }],
        },
        MoleculeDefaults::concrete(),
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
                old: NumForm::Lit(0),
                new: NumForm::Lit(-1),
            },
        },
        MoleculeDefaults::new(),
        r##"{:atom {:modify [0 {:expect "#c0" :update "#c-"}]}}"##,
    )]
    #[case::atom_constraint_add(
        Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintForm::valence(4_i64)),
        },
        MoleculeDefaults::new(),
        r##"{:atom {:modify [0 {:expect "#v*" :update "#v4"}]}}"##,
    )]
    #[case::constraint_add(
        Edit::AddMoleculeConstraint {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }).into(),
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
            .into_iter()
            .map(|input| input.to_edn())
            .collect::<Vec<_>>();

        assert_eq!(rendered, vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_add(
        r#"{:dative-bond {:add {:donors [0 {:new 0}] :acceptor 2 :attrs :single}}}"#,
        Edit::AddDativeBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            attributes: DativeBondForm::from_order(1),
        },
    )]
    #[case::dative_remove(
        r#"{:dative-bonds {:remove [{:id 0 :donors [1] :acceptor {:new 2} :attrs :single} {:id {:new 0} :donors [{:new 1}] :acceptor 3 :attrs :double}]}}"#,
        Edit::RemoveDativeBonds { removes: vec![
            (DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(1)), AtomHandle::New(2)], DativeBondForm::from_order(1)),
            (DativeBondHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(3))], DativeBondForm::from_order(2)),
        ] },
    )]
    #[case::dative_field(
        r#"{:dative-bond {:modify [{:new 0} {:expect "1" :update "2"}]}}"#,
        Edit::ModifyDativeBondField {
            id: DativeBondHandle::New(0),
            change: DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
        },
    )]
    #[case::dative_constraint(
        r##"{:dative-bond {:modify [0 {:expect "#a*" :update "#a"}]}}"##,
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::Id(DativeBondId(0)),
            old: None,
            new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        },
    )]
    #[case::dative_constraint_remove(
        r##"{:dative-bond {:modify [{:new 1} {:expect "#a" :update "#a*"}]}}"##,
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::New(1),
            old: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            new: None,
        },
    )]
    #[case::aromatic_add(
        r#"{:aromatic-system {:add {:atoms [0 {:new 0} 2] :attrs "[1,1,1]"}}}"#,
        Edit::AddAromaticSystem {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            attributes: AromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
    )]
    #[case::aromatic_remove(
        r#"{:aromatic-systems {:remove [{:id 0 :atoms [0 {:new 0}] :attrs "[1,1]"} {:id {:new 0} :atoms [{:new 1} 2] :attrs "[2,2]"}]}}"#,
        Edit::RemoveAromaticSystems { removes: vec![
            (AromaticSystemHandle::Id(AromaticSystemId(0)), vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], AromaticSystemForm::from_electrons(vec![1, 1])),
            (AromaticSystemHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(2))], AromaticSystemForm::from_electrons(vec![2, 2])),
        ] },
    )]
    #[case::aromatic_field(
        r##"{:aromatic-system {:modify [{:new 0} {:expect "#c0" :update "#c-"}]}}"##,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::New(0),
            change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Lit(-1) },
        },
    )]
    #[case::aromatic_electrons(
        r#"{:aromatic-system {:modify [1 {:expect "[1,1]" :update "[2,0]"}]}}"#,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(1)),
            change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1]), new: ElectronCountsForm::Lit(vec![2, 0]) },
        },
    )]
    #[case::aromatic_unpaired_electrons(
        r##"{:aromatic-system {:modify [{:new 1} {:expect "#u0#s" :update "#u2#s3"}]}}"##,
        Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::New(1),
            change: AromaticSystemFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(1) },
                new: UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) },
            },
        },
    )]
    #[case::aromatic_constraint(
        r##"{:aromatic-system {:modify [0 {:expect "#e*" :update "#e6"}]}}"##,
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            old: None,
            new: Some(AromaticSystemConstraintForm::electron_count(NumForm::Lit(6))),
        },
    )]
    #[case::aromatic_constraint_remove(
        r##"{:aromatic-system {:modify [{:new 2} {:expect "#e6" :update "#e*"}]}}"##,
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::New(2),
            old: Some(AromaticSystemConstraintForm::electron_count(NumForm::Lit(6))),
            new: None,
        },
    )]
    #[case::multicenter_add(
        r#"{:multicenter-bond {:add {:atoms [0 {:new 0} 2] :attrs "[1,1,0]"}}}"#,
        Edit::AddMulticenterBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0), AtomHandle::Id(AtomId(2))],
            attributes: MulticenterBondForm::from_electrons(vec![1, 1, 0]),
        },
    )]
    #[case::multicenter_remove(
        r#"{:multicenter-bonds {:remove [{:id 0 :atoms [0 {:new 0}] :attrs "[1,1]"} {:id {:new 0} :atoms [{:new 1} 2] :attrs "[2,0]"}]}}"#,
        Edit::RemoveMulticenterBonds { removes: vec![
            (MulticenterBondHandle::Id(MulticenterBondId(0)), vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], MulticenterBondForm::from_electrons(vec![1, 1])),
            (MulticenterBondHandle::New(0), vec![AtomHandle::New(1), AtomHandle::Id(AtomId(2))], MulticenterBondForm::from_electrons(vec![2, 0])),
        ] },
    )]
    #[case::multicenter_field(
        r#"{:multicenter-bond {:modify [{:new 0} {:expect "[1,1]" :update "[2,0]"}]}}"#,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::New(0),
            change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1]), new: ElectronCountsForm::Lit(vec![2, 0]) },
        },
    )]
    #[case::multicenter_charge(
        r##"{:multicenter-bond {:modify [1 {:expect "#c0" :update "#c+"}]}}"##,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(1)),
            change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Lit(1) },
        },
    )]
    #[case::multicenter_unpaired_electrons(
        r##"{:multicenter-bond {:modify [{:new 1} {:expect "#u0#s" :update "#u2#s3"}]}}"##,
        Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::New(1),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(1) },
                new: UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(3) },
            },
        },
    )]
    #[case::multicenter_constraint(
        r##"{:multicenter-bond {:modify [0 {:expect "#e*" :update "#e2"}]}}"##,
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            old: None,
            new: Some(MulticenterBondConstraintForm::electron_count(NumForm::Lit(2))),
        },
    )]
    #[case::multicenter_constraint_remove(
        r##"{:multicenter-bond {:modify [{:new 2} {:expect "#e2" :update "#e*"}]}}"##,
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::New(2),
            old: Some(MulticenterBondConstraintForm::electron_count(NumForm::Lit(2))),
            new: None,
        },
    )]
    #[case::noncovalent_add(
        r#"{:noncovalent-bond {:add {:atoms [0 {:new 0}] :attrs "Hbd"}}}"#,
        Edit::AddNoncovalentBond {
            atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        },
    )]
    #[case::noncovalent_remove(
        r#"{:noncovalent-bonds {:remove [{:id 0 :atoms [0 {:new 0}] :attrs "Hbd"} {:id {:new 0} :atoms [{:new 1} 2] :attrs "Ion"}]}}"#,
        Edit::RemoveNoncovalentBonds { removes: vec![
            (NoncovalentBondHandle::Id(NoncovalentBondId(0)), [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)], NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)),
            (NoncovalentBondHandle::New(0), [AtomHandle::New(1), AtomHandle::Id(AtomId(2))], NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic)),
        ] },
    )]
    #[case::noncovalent_field(
        r#"{:noncovalent-bond {:modify [{:new 0} {:expect "Hbd" :update "Ion"}]}}"#,
        Edit::ModifyNoncovalentBondField {
            id: NoncovalentBondHandle::New(0),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            },
        },
    )]
    #[case::noncovalent_constraint(
        r##"{:noncovalent-bond {:modify [0 {:expect "#I*" :update "#I"}]}}"##,
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            old: None,
            new: Some(NoncovalentBondConstraintForm::intramolecular(true)),
        },
    )]
    #[case::noncovalent_constraint_remove(
        r##"{:noncovalent-bond {:modify [{:new 1} {:expect "#I" :update "#I*"}]}}"##,
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::New(1),
            old: Some(NoncovalentBondConstraintForm::intramolecular(true)),
            new: None,
        },
    )]
    #[case::stereo_atom_add(
        r#"{:stereo-atom {:add {:site {:new 0} :ligands [0 [:h {:new 1}] [:lp 2] {:new 3}] :attrs :ccw}}}"#,
        Edit::AddStereoAtom {
            site: AtomHandle::New(0),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::New(1), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::LonePair),
                (AtomHandle::New(3), StereoLigandKind::Atom),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        },
    )]
    #[case::stereo_atom_remove(
        r#"{:stereo-atoms {:remove [{:id 0 :site 1 :ligands [2 [:h 3] [:lp {:new 0}] 4] :attrs :cw} {:id {:new 1} :site {:new 2} :ligands [5 6 7 8] :attrs :ccw}]}}"#,
        Edit::RemoveStereoAtoms { removes: vec![
            (
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(1)),
                vec![
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::ImplicitHydrogen),
                    (AtomHandle::New(0), StereoLigandKind::LonePair),
                    (AtomHandle::Id(AtomId(4)), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
            ),
            (
                StereoAtomHandle::New(1),
                AtomHandle::New(2),
                vec![
                    (AtomHandle::Id(AtomId(5)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(6)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(7)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(8)), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
            ),
        ] },
    )]
    #[case::stereo_atom_field(
        r#"{:stereo-atom {:modify [{:new 0} {:expect "Th0" :update "Th1"}]}}"#,
        Edit::ModifyStereoAtomField {
            id: StereoAtomHandle::New(0),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32),
            },
        },
    )]
    #[case::stereo_atom_constraint_add(
        r##"{:stereo-atom {:modify [0 {:expect "Th#g*" :update "Th#g/"}]}}"##,
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        },
    )]
    #[case::stereo_atom_constraint_remove(
        r##"{:stereo-atom {:modify [{:new 1} {:expect "Th#g/" :update "Th#g*"}]}}"##,
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::New(1),
            kind: Some(StereoKind::Tetrahedral),
            old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
            new: None,
        },
    )]
    #[case::stereo_bond_add(
        r#"{:stereo-bond {:add {:site {:new 0} :ligands [0 [:h {:new 1}] [:lp 2] {:new 3}] :attrs :z}}}"#,
        Edit::AddStereoBond {
            site: BondHandle::New(0),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::New(1), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::LonePair),
                (AtomHandle::New(3), StereoLigandKind::Atom),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32),
        },
    )]
    #[case::stereo_bond_remove(
        r#"{:stereo-bonds {:remove [{:id 0 :site 1 :ligands [2 3 4 5] :attrs :e} {:id {:new 1} :site {:new 2} :ligands [[:h 6] 7 [:lp {:new 0}] 8] :attrs :z}]}}"#,
        Edit::RemoveStereoBonds { removes: vec![
            (
                StereoBondHandle::Id(StereoBondId(0)),
                BondHandle::Id(BondId(1)),
                vec![
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(4)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(5)), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
            (
                StereoBondHandle::New(1),
                BondHandle::New(2),
                vec![
                    (AtomHandle::Id(AtomId(6)), StereoLigandKind::ImplicitHydrogen),
                    (AtomHandle::Id(AtomId(7)), StereoLigandKind::Atom),
                    (AtomHandle::New(0), StereoLigandKind::LonePair),
                    (AtomHandle::Id(AtomId(8)), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 0_u32),
            ),
        ] },
    )]
    #[case::stereo_bond_field_clear(
        r#"{:stereo-bond {:modify [0 {:expect "Ct1" :update "*"}]}}"#,
        Edit::ModifyStereoBondField {
            id: StereoBondHandle::Id(StereoBondId(0)),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32),
                new: StereoConfigurationForm::Undetermined,
            },
        },
    )]
    #[case::stereo_bond_constraint_replace(
        r##"{:stereo-bond {:modify [{:new 0} {:expect "Ct#g/" :update "Ct#g="}]}}"##,
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::New(0),
            kind: Some(StereoKind::CisTrans),
            old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
            new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        },
    )]
    #[case::stereo_bond_constraint_add(
        r##"{:stereo-bond {:modify [0 {:expect "Ct#g*" :update "Ct#g/"}]}}"##,
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::Id(StereoBondId(0)),
            kind: Some(StereoKind::CisTrans),
            old: None,
            new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        },
    )]
    #[case::stereo_bond_constraint_remove(
        r##"{:stereo-bond {:modify [{:new 1} {:expect "Ct#g/" :update "Ct#g*"}]}}"##,
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::New(1),
            kind: Some(StereoKind::CisTrans),
            old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
            new: None,
        },
    )]
    fn test_overlay_edit_input_roundtrip(#[case] input: &str, #[case] expected: Edit) {
        let parsed = EditInput::from_edn_str(input).unwrap();
        let mut edits = Edits::new();

        parsed.append_to(&mut edits, &MoleculeDefaults::new()).unwrap();
        let rendered = EditInput::from_edit(&expected, &MoleculeDefaults::new())
            .unwrap()
            .into_iter()
            .map(|input| input.to_edn())
            .collect::<Vec<_>>();

        assert_eq!(edits, Edits::from_iter([expected]));
        assert_eq!(rendered, vec![read_string(input).unwrap()]);
    }

    #[rstest]
    #[case::aromatic_add(
        r#"{:aromatic-system {:add {:atoms [0 1] :attrs "[1,1]"}}}"#,
        Edit::AddAromaticSystem {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: AromaticSystemForm::from_electrons(vec![1, 1]).into_concrete(),
        },
    )]
    #[case::aromatic_remove(
        r#"{:aromatic-systems {:remove [{:id 0 :atoms [0 1] :attrs "[1,1]"}]}}"#,
        Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                AromaticSystemForm::from_electrons(vec![1, 1]).into_concrete(),
            )],
        },
    )]
    #[case::multicenter_add(
        r#"{:multicenter-bond {:add {:atoms [0 1] :attrs "[1,1]"}}}"#,
        Edit::AddMulticenterBond {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: MulticenterBondForm::from_electrons(vec![1, 1]).into_concrete(),
        },
    )]
    #[case::multicenter_remove(
        r#"{:multicenter-bonds {:remove [{:id 0 :atoms [0 1] :attrs "[1,1]"}]}}"#,
        Edit::RemoveMulticenterBonds {
            removes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                MulticenterBondForm::from_electrons(vec![1, 1]).into_concrete(),
            )],
        },
    )]
    fn test_overlay_edit_input_ground_defaults(#[case] input: &str, #[case] expected: Edit) {
        let mut edits = Edits::new();

        EditInput::from_edn_str(input)
            .unwrap()
            .append_to(&mut edits, &MoleculeDefaults::concrete())
            .unwrap();
        let rendered = EditInput::from_edit(&expected, &MoleculeDefaults::concrete())
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
    #[case::stereo_atom_field(
        r#"{:stereo-atom {:modify [0 {:expect "Th0" :update "Th"}]}}"#,
        EdnError::De(DeError::Custom(
            "stereo-atom :modify :expect and :update must address the same field and constraints"
                .to_string(),
        )),
    )]
    #[case::stereo_atom_constraint_kind(
        r##"{:stereo-atom {:modify [0 {:expect "Th#g/" :update "Ct#g="}]}}"##,
        EdnError::De(DeError::Custom(
            "stereo-atom constraint changes require the same stereo kind in :expect and :update"
                .to_string(),
        )),
    )]
    #[case::stereo_bond_ligand_kind(
        r#"{:stereo-bond {:add {:site 0 :ligands [[:x 1]] :attrs :z}}}"#,
        EdnError::De(DeError::Custom(
            "unknown stereo ligand kind :x".to_string(),
        )),
    )]
    #[case::constraint_keyword(
        "{:constraint {:add {:atom [:carbon {:valence 4}]}}}",
        EdnError::De(DeError::TypeMismatch {
            expected: "edit handle (non-negative integer or {:new n} map)",
            got: "keyword",
            path: Vec::new(),
        }),
    )]
    fn test_edit_input_from_edn_error(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(EditInput::from_edn_str(input), Err(expected));
    }

    #[rstest]
    #[case::stereo_atom(
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            kind: None,
            old: None,
            new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            ))),
        },
        DeError::Custom("stereo-atom constraint edit requires a stereo kind".to_string()),
    )]
    #[case::stereo_bond(
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::Id(StereoBondId(0)),
            kind: None,
            old: None,
            new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            ))),
        },
        DeError::Custom("stereo-bond constraint edit requires a stereo kind".to_string()),
    )]
    fn test_edit_input_from_edit_requires_stereo_kind(
        #[case] edit: Edit,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            EditInput::from_edit(&edit, &MoleculeDefaults::new()),
            Err(expected)
        );
    }

    #[rstest]
    #[case::incident_bond(
        mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
        "{:topology {:remove {:atoms [0] :bonds [0]}}}",
        mol_dsl!(r#"{:atoms ["N"]}"#),
    )]
    fn test_edit_input_append_to_topology(
        #[case] molecule: Molecule,
        #[case] input: &str,
        #[case] expected: Molecule,
    ) {
        let mut edits = Edits::new();
        EditInput::from_edn_str(input)
            .unwrap()
            .append_to(&mut edits, &MoleculeDefaults::new())
            .unwrap();

        assert_eq!(molecule.apply(edits), Ok(expected));
    }
    #[rustfmt::skip]
    #[rstest]
    #[case::entity_leaf(
        "{:atom [2 {:valence 4}]}")]
    #[case::logical_repeated_handle(
        "{:and [{:atom [{:new 0} {:valence 4}]} {:not {:atom [{:new 0} {:degree 3}]}}]}")]
    #[case::entity_kinds(
        "{:and [{:bond [{:new 0} {:aromatic true}]} {:dative-bond [{:new 1} {:aromatic false}]} {:aromatic-system [{:new 2} {:electron-count 6}]} {:multicenter-bond [{:new 3} {:electron-count 2}]} {:noncovalent-bond [{:new 4} {:intramolecular true}]} {:stereo-atom [{:new 5} [:tetrahedral {:stereogenicity {:relation :stereogenic}}]]} {:stereo-bond [{:new 6} [:cis-trans {:stereogenicity {:relation :stereogenic}}]]}]}")]
    #[case::relational_kinds(
        "{:and [{:dative-bond-parallels [{:new 0} {:new 0}]} {:aromatic-system-contains [{:new 0} {:new 0}]} {:multicenter-bond-contains [{:new 0} {:new 0}]} {:noncovalent-bond-contains [{:new 0} {:new 0}]} {:stereo-atom-site [{:new 0} {:new 0}]} {:stereo-bond-site [{:new 0} {:new 0}]}]}")]
    #[case::quantified_predicate(
        "{:aromatic-system-all-atoms [{:new 0} {:valence 3}]}")]
    #[case::atom_subset(
        "{:charge-sum {:atoms [1 {:new 0}] :sum 0}}")]
    #[case::bond_subset(
        "{:bond-order-sum {:bonds [2 {:new 0}] :sum 3}}")]
    #[case::whole_molecule(
        "{:connected {}}")]
    fn test_constraint_edit_dsl_roundtrip(#[case] input: &str) {
        let parsed = ConstraintEditDsl::from_edn_str(input).unwrap();
        let rebuilt = ConstraintEditDsl::from_edit(parsed.clone().into_edit());

        assert_eq!(parsed.to_edn(), read_string(input).unwrap());
        assert_eq!(rebuilt, parsed);
    }

    #[rstest]
    #[case::repeated(
        "{:and [{:atom [{:new 2} {:valence 4}]} {:atom [{:new 2} {:degree 3}]}]}",
        ConstraintHandles {
            atoms: vec![AtomHandle::New(2)],
            ..ConstraintHandles::default()
        },
    )]
    #[case::quantified(
        "{:aromatic-system-all-atoms [{:new 1} {:valence 3}]}",
        ConstraintHandles {
            aromatic_systems: vec![AromaticSystemHandle::New(1)],
            ..ConstraintHandles::default()
        },
    )]
    fn test_constraint_edit_dsl_handles(#[case] input: &str, #[case] expected: ConstraintHandles) {
        assert_eq!(
            ConstraintEditDsl::from_edn_str(input).unwrap().handles,
            expected
        );
    }

    #[rstest]
    #[case::keyword(
        "{:atom [:carbon {:valence 4}]}",
        DeError::TypeMismatch {
            expected: "edit handle (non-negative integer or {:new n} map)",
            got: "keyword",
            path: Vec::new(),
        },
    )]
    #[case::structural(
        "{:bond [{:atoms [0 1]} {:aromatic true}]}",
        DeError::MissingField {
            key: "new".to_string(),
            path: Vec::new(),
        },
    )]
    fn test_constraint_edit_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            ConstraintEditDsl::from_edn_str(input),
            Err(EdnError::De(expected))
        );
    }
}
