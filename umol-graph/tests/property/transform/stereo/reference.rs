//! Independent binary masks, domain membership, source-frame parity, and joint XOR.
//! This reference does not call perception, validation, or a transformer.

use std::collections::BTreeSet;

use umol_graph::ops::stereo::StereoSite;
use umol_graph_ir::ir::{
    AtomConstraintKey, BondConstraintKey, Constraint, MoleculeEntries, StereoAtomConstraintForm,
    StereoAtomId, StereoBondConstraintForm, StereoBondId, StereoCoset, StereoKind, Stereogenicity,
    StereogenicityForm,
};

use super::strategies::Site;

pub(super) fn selected(
    info: &Site,
    index: usize,
    sites: u8,
    kinds: u8,
    classes: Option<u8>,
) -> bool {
    let kind = if info.kind == StereoKind::Tetrahedral {
        1
    } else {
        2
    };
    let class = match info.class {
        Stereogenicity::Symmetric => 1,
        Stereogenicity::Prochiral => 2,
        Stereogenicity::Stereogenic => 4,
    };
    sites & (1 << index) != 0
        && kinds & kind != 0
        && classes.is_none_or(|classes| classes & class != 0)
}

pub(super) fn assignments(
    sites: &[Site],
    selected: &BTreeSet<usize>,
    order: &[usize],
    joint: Option<(usize, usize)>,
) -> Vec<Vec<StereoCoset>> {
    let mut results = Vec::new();
    for bits in 0..(1usize << selected.len()) {
        let mut candidate: Vec<_> = sites.iter().map(|site| site.coset.clone()).collect();
        let mut values = vec![0u32; sites.len()];
        let mut admitted = true;
        for (bit, &site) in selected.iter().enumerate() {
            let value = ((bits >> bit) & 1) as u32;
            values[site] = value;
            admitted &= sites[site].admitted & (1 << value) != 0;
            admitted &= sites[site].asserted & (1 << (value ^ sites[site].parity)) != 0;
            candidate[site] = StereoCoset::Lit(value);
        }
        if let Some((left, right)) = joint {
            admitted &=
                (values[left] ^ sites[left].parity) != (values[right] ^ sites[right].parity);
        }
        if admitted {
            let key: Vec<_> = order
                .iter()
                .filter(|index| selected.contains(index))
                .map(|&index| values[index])
                .collect();
            results.push((key, candidate));
        }
    }
    results.sort_by(|left, right| left.0.cmp(&right.0));
    results.into_iter().map(|(_, values)| values).collect()
}

/// Filter the generated entries and number retained records independently of IR compaction.
pub(super) fn removed_entries(
    source: &MoleculeEntries,
    sites: &[Site],
    removed: &BTreeSet<StereoSite>,
) -> MoleculeEntries {
    let mut entries = source.clone();
    entries
        .stereo_atoms
        .retain(|(site, _, _)| !removed.contains(&StereoSite::Atom(*site)));
    entries
        .stereo_bonds
        .retain(|(site, _, _)| !removed.contains(&StereoSite::Bond(*site)));
    for site in removed {
        match site {
            StereoSite::Atom(atom) => {
                entries.atoms[atom.index()]
                    .constraints
                    .remove(AtomConstraintKey::TetrahedralStereo);
            }
            StereoSite::Bond(bond) => {
                entries.bonds[bond.index()]
                    .2
                    .constraints
                    .remove(BondConstraintKey::CisTransStereo);
            }
        }
    }
    let mut constraints = Vec::new();
    let mut atoms = 0;
    let mut bonds = 0;
    for info in sites.iter().filter(|info| !removed.contains(&info.site)) {
        match info.site {
            StereoSite::Atom(_) => {
                constraints.push(Constraint::StereoAtom(
                    StereoAtomId(atoms),
                    info.kind,
                    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                ));
                atoms += 1;
            }
            StereoSite::Bond(_) => {
                constraints.push(Constraint::StereoBond(
                    StereoBondId(bonds),
                    info.kind,
                    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                ));
                bonds += 1;
            }
        }
    }
    constraints.push(source.constraints.iter().last().unwrap().clone());
    entries.constraints = constraints.into();
    entries
}
