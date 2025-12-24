//! Edit definitions for GraphIR.

use serde::{Deserialize, Serialize};
use strum::{AsRefStr, EnumDiscriminants, EnumIter};

#[derive(
    Debug, Default, Clone, PartialEq, Serialize, Deserialize, EnumDiscriminants, EnumIter, AsRefStr,
)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Edit {
    #[default]
    NoOp,
}
// SetAtomCharge {
//     atom: usize,
//     charge: i32,
// },
//     SetAtomExplicitHCount {
//         atom: usize,
//         count: u32,
//     },
//     SetAtomImplicitHCount {
//         atom: usize,
//         count: u32,
//     },
//     SetAtomImplicitH {
//         atom: usize,
//         implicit: bool,
//     },
//     SetAtomAromaticFlag {
//         atom: usize,
//         aromatic: Option<bool>,
//     },
//     SetAtomChirality {
//         atom: usize,
//         chirality: Option<Chirality>,
//     },
//     SetAtomClass {
//         atom: usize,
//         class: Option<u32>,
//     },
//     SetAtomUnpairedECount {
//         atom: usize,
//         count: u32,
//     },
//     SetBondOrder {
//         bond: usize,
//         order: BondOrder,
//     },
//     SetBondSymbol {
//         bond: usize,
//         symbol: BondSymbol,
//     },
//     SetBondDirection {
//         bond: usize,
//         direction: Option<BondDir>,
//     },
//     SetBondStereo {
//         bond: usize,
//         stereo: Option<BondStereo>,
//     },
//     AddAtom {
//         atom: usize,
//         element: Element,
//     },
//     RemoveAtom {
//         atom: usize,
//     },
//     AddBond {
//         bond: usize,
//         atoms: (usize, usize),
//         order: BondOrder,
//     },
//     RemoveBond {
//         bond: usize,
//     },
//     RetargetBond {
//         bond: usize,
//         atoms: (usize, usize),
//     },
//     SetBondRing {
//         bond: usize,
//         ring: Option<u32>,
//     },
// }

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct EditList {
    pub edits: Vec<Edit>,
}

impl EditList {
    pub fn new() -> Self {
        Self { edits: Vec::new() }
    }

    pub fn push(&mut self, edit: Edit) {
        self.edits.push(edit);
    }

    pub fn extend<I: IntoIterator<Item = Edit>>(&mut self, edits: I) {
        self.edits.extend(edits);
    }

    pub fn append_list(&mut self, other: &mut EditList) {
        self.edits.append(&mut other.edits);
    }

    pub fn iter(&self) -> impl Iterator<Item = &Edit> {
        self.edits.iter()
    }

    pub fn into_vec(self) -> Vec<Edit> {
        self.edits
    }
}

impl IntoIterator for EditList {
    type Item = Edit;
    type IntoIter = std::vec::IntoIter<Edit>;

    fn into_iter(self) -> Self::IntoIter {
        self.edits.into_iter()
    }
}
