//! Finite explicit Th/Ct sites with independently recorded ligand classes and domains.
//! Addition and classification cases have undetermined cosets. Absolute cases vary finite
//! domains and source assertions. All cases retain an unrelated partial atom.

use std::collections::BTreeSet;
use std::ops::RangeInclusive;

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::stereo::StereoSite;
use umol_graph_core::{GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondConstraintForm, BondConstraintKey,
    BondForm, BondId, CisTransStereoForm, Constraint, IsotopeMassForm, Molecule, MoleculeEntries,
    MoleculeRemapping, NumForm, RingScope, StereoAtomConstraintForm, StereoAtomForm, StereoAtomId,
    StereoBondConstraintForm, StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind, Stereogenicity, StereogenicityForm, TetrahedralStereoForm,
};

#[derive(Clone, Copy, Debug)]
pub(super) enum Domain {
    Addition,
    Undetermined,
    Absolute,
    AbsoluteNonempty,
}

#[derive(Clone, Debug)]
pub(super) struct Site {
    pub site: StereoSite,
    pub kind: StereoKind,
    pub natural: [AtomId; 4],
    pub class: Stereogenicity,
    pub coset: StereoCoset,
    pub admitted: u8,
    pub asserted: u8,
    pub parity: u32,
    pub record: Option<usize>,
}

#[derive(Clone, Debug)]
pub(super) struct StereoCase {
    pub source: Molecule,
    pub entries: MoleculeEntries,
    pub sites: Vec<Site>,
    pub images: Vec<AtomId>,
    pub order: Vec<usize>,
}

pub(super) fn stereo_cases(
    domain: Domain,
    count: RangeInclusive<usize>,
) -> impl Strategy<Value = StereoCase> {
    let labels = prop_oneof![
        prop::array::uniform4(0usize..5),
        Just([1, 2, 3, 4]),
        Just([0, 0, 1, 2]),
        Just([0, 0, 0, 0]),
    ];
    prop::collection::vec(
        (
            any::<bool>(),
            labels,
            prop::array::uniform4(any::<u16>()),
            0u8..8,
            0u8..7,
            0u8..6,
            any::<bool>(),
            any::<bool>(),
            any::<u16>(),
        ),
        count,
    )
    .prop_flat_map(move |specs| {
        let atoms = 1 + specs
            .iter()
            .map(|spec| if spec.0 { 5 } else { 6 })
            .sum::<usize>();
        (Just(specs), prop::collection::vec(any::<u16>(), atoms))
    })
    .prop_map(move |(specs, keys)| {
        let mut entries = MoleculeEntries::default();
        let mut sites = Vec::new();
        let mut order_keys = Vec::new();
        let elements = [Element::H, Element::F, Element::Cl, Element::Br, Element::I];
        for (
            th,
            labels,
            frame_keys,
            ct_action,
            domain_index,
            assertion,
            stored,
            isotope,
            order_key,
        ) in specs
        {
            let base = entries.atoms.len() as u32;
            let mut center = AtomForm::from_element(Element::C).into_concrete();
            center
                .constraints
                .set(AtomConstraintForm::Degree(NumForm::Lit(if th {
                    4
                } else {
                    3
                })));
            if isotope {
                center.isotope_mass = IsotopeMassForm::Lit(13);
            }
            entries.atoms.push(center.clone());
            if !th {
                entries.atoms.push(center);
            }
            let natural = [0, 1, 2, 3].map(|i| AtomId(base + if th { 1 } else { 2 } + i));
            entries.atoms.extend(
                labels.map(|label| AtomForm::from_element(elements[label]).into_concrete()),
            );
            let site = if th {
                entries.bonds.extend(natural.map(|ligand| {
                    (
                        AtomId(base),
                        ligand,
                        BondForm::from_order(1).into_concrete(),
                    )
                }));
                StereoSite::Atom(AtomId(base))
            } else {
                let bond = BondId(entries.bonds.len() as u32);
                let attributes = BondForm::from_order(2)
                    .into_concrete()
                    .with_constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 0));
                entries
                    .bonds
                    .push((AtomId(base), AtomId(base + 1), attributes));
                entries
                    .bonds
                    .extend(natural.into_iter().enumerate().map(|(i, ligand)| {
                        (
                            AtomId(base + u32::from(i >= 2)),
                            ligand,
                            BondForm::from_order(1).into_concrete(),
                        )
                    }));
                StereoSite::Bond(bond)
            };
            let kind = if th {
                StereoKind::Tetrahedral
            } else {
                StereoKind::CisTrans
            };
            let distinct = labels.into_iter().collect::<BTreeSet<_>>().len();
            let class = if th {
                match distinct {
                    4 => Stereogenicity::Stereogenic,
                    3 => Stereogenicity::Prochiral,
                    _ => Stereogenicity::Symmetric,
                }
            } else if labels[0] != labels[1] && labels[2] != labels[3] {
                Stereogenicity::Stereogenic
            } else {
                Stereogenicity::Symmetric
            };
            let mut image = [0usize, 1, 2, 3];
            if th {
                image.sort_by_key(|&i| (frame_keys[i], i));
            } else {
                if ct_action & 1 != 0 {
                    image.swap(0, 1);
                }
                if ct_action & 2 != 0 {
                    image.swap(2, 3);
                }
                if ct_action & 4 != 0 {
                    image.rotate_left(2);
                }
            }
            let parity = (0..4)
                .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                .filter(|&(i, j)| image[i] > image[j])
                .count() as u32
                % 2;
            let (coset, admitted) = if matches!(domain, Domain::Absolute | Domain::AbsoluteNonempty)
            {
                let index = if matches!(domain, Domain::AbsoluteNonempty) {
                    domain_index % 6
                } else {
                    domain_index
                };
                match index {
                    0 => (StereoCoset::Undetermined, 3),
                    1 => (StereoCoset::Lit(0), 1),
                    2 => (StereoCoset::Lit(1), 2),
                    3 => (StereoCoset::lit_set([0, 1]), 3),
                    4 => (StereoCoset::lit_set([0]), 1),
                    5 => (StereoCoset::lit_set([1]), 2),
                    _ => (StereoCoset::lit_set([]), 0),
                }
            } else {
                (StereoCoset::Undetermined, 3)
            };
            let assertion = if matches!(domain, Domain::Addition) {
                assertion % 2
            } else {
                assertion
            };
            let (asserted_coset, asserted) = match assertion {
                0 => (None, 3),
                1 => (Some(StereoCoset::Undetermined), 3),
                2 => (Some(StereoCoset::Lit(0)), 1),
                3 => (Some(StereoCoset::Lit(1)), 2),
                4 => (None, 0),
                _ => (Some(StereoCoset::lit_set([0, 1])), 3),
            };
            match site {
                StereoSite::Atom(atom) => {
                    if let Some(coset) = asserted_coset {
                        entries.atoms[atom.index()].constraints.set(
                            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                                coset,
                            )),
                        );
                    } else if asserted == 0 {
                        entries.atoms[atom.index()].constraints.set(
                            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo),
                        );
                    }
                }
                StereoSite::Bond(bond) => {
                    if let Some(coset) = asserted_coset {
                        entries.bonds[bond.index()].2.constraints.set(
                            BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset)),
                        );
                    } else if asserted == 0 {
                        entries.bonds[bond.index()].2.constraints.set(
                            BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo),
                        );
                    }
                }
            }
            let record = if !matches!(domain, Domain::Addition) || stored {
                let ligands = image
                    .map(|i| StereoLigand::new(natural[i], StereoLigandKind::Atom))
                    .to_vec();
                match site {
                    StereoSite::Atom(atom) => {
                        let id = entries.stereo_atoms.len();
                        entries.stereo_atoms.push((
                            atom,
                            ligands,
                            StereoAtomForm::new(kind, coset.clone()),
                        ));
                        entries.constraints.push(Constraint::StereoAtom(
                            StereoAtomId(id as u32),
                            kind,
                            StereoAtomConstraintForm::Stereogenicity(
                                StereogenicityForm::Undetermined,
                            ),
                        ));
                        Some(id)
                    }
                    StereoSite::Bond(bond) => {
                        let id = entries.stereo_bonds.len();
                        entries.stereo_bonds.push((
                            bond,
                            ligands,
                            StereoBondForm::new(kind, coset.clone()),
                        ));
                        entries.constraints.push(Constraint::StereoBond(
                            StereoBondId(id as u32),
                            kind,
                            StereoBondConstraintForm::Stereogenicity(
                                StereogenicityForm::Undetermined,
                            ),
                        ));
                        Some(id)
                    }
                }
            } else {
                None
            };
            sites.push(Site {
                site,
                kind,
                natural,
                class,
                coset,
                admitted,
                asserted,
                parity,
                record,
            });
            order_keys.push(order_key);
        }
        let noise = AtomId(entries.atoms.len() as u32);
        let mut partial = AtomForm::from_element(Element::Ne).into_concrete();
        partial.charge = NumForm::var("retained_charge");
        entries.atoms.push(partial);
        entries.constraints.push(Constraint::Atom(
            noise,
            AtomConstraintForm::Degree(NumForm::Lit(0)),
        ));
        let mut sorted: Vec<_> = (0..keys.len()).collect();
        sorted.sort_by_key(|&i| (keys[i], i));
        if sorted.len() > 1 && sorted.iter().copied().eq(0..keys.len()) {
            sorted.swap(0, 1);
        }
        let mut images = vec![AtomId(0); sorted.len()];
        for (target, source) in sorted.into_iter().enumerate() {
            images[source] = AtomId(target as u32);
        }
        let mut order: Vec<_> = (0..sites.len()).collect();
        order.sort_by_key(|&i| (order_keys[i], i));
        StereoCase {
            source: Molecule::from_entries(entries.clone()),
            entries,
            sites,
            images,
            order,
        }
    })
}

/// Renumber ids and transport the assertions' implicit, sorted ligand frames.
pub(super) fn renumber(molecule: Molecule, case: &StereoCase) -> Molecule {
    let mapping = MoleculeRemapping::new(
        GraphRemapping::new(
            Remapping::new(case.images.iter().map(|id| NodeId(id.0)).collect()).unwrap(),
            Remapping::identity(molecule.bonds().count()),
        ),
        Remapping::identity(0),
        Remapping::identity(0),
        Remapping::identity(0),
        Remapping::identity(0),
        Remapping::identity(molecule.stereo_atoms().count()),
        Remapping::identity(molecule.stereo_bonds().count()),
    );
    let mut result = molecule.remap(&mapping);
    for info in &case.sites {
        let mapped = info.natural.map(|id| case.images[id.index()]);
        let parity = match info.site {
            StereoSite::Atom(_) => {
                (0..4)
                    .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                    .filter(|&(i, j)| mapped[i] > mapped[j])
                    .count() as u32
                    % 2
            }
            StereoSite::Bond(_) => {
                u32::from(mapped[0] > mapped[1]) ^ u32::from(mapped[2] > mapped[3])
            }
        };
        match info.site {
            StereoSite::Atom(atom) => {
                let mut view = result.atom_mut(case.images[atom.index()]);
                let attributes = view.attributes_mut();
                if let Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                    coset,
                ))) = attributes
                    .constraints
                    .get(AtomConstraintKey::TetrahedralStereo)
                    .cloned()
                {
                    attributes
                        .constraints
                        .set(AtomConstraintForm::TetrahedralStereo(
                            TetrahedralStereoForm::Stereo(flip_coset(coset, parity)),
                        ));
                }
            }
            StereoSite::Bond(bond) => {
                let mut view = result.bond_mut(bond);
                let attributes = view.attributes_mut();
                if let Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) =
                    attributes
                        .constraints
                        .get(BondConstraintKey::CisTransStereo)
                        .cloned()
                {
                    attributes
                        .constraints
                        .set(BondConstraintForm::CisTransStereo(
                            CisTransStereoForm::Stereo(flip_coset(coset, parity)),
                        ));
                }
            }
        }
    }
    result
}

fn flip_coset(coset: StereoCoset, parity: u32) -> StereoCoset {
    match coset {
        StereoCoset::Lit(value) => StereoCoset::Lit(value ^ parity),
        StereoCoset::LitSet(values) => {
            StereoCoset::lit_set(values.into_iter().map(|value| value ^ parity))
        }
        other => other,
    }
}
