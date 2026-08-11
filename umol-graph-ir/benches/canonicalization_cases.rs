use super::*;

pub(super) struct CorpusCase {
    pub(super) name: &'static str,
    pub(super) molecule: Molecule,
}

pub(super) const LEVELS: [IncidenceLevel; 3] = [
    IncidenceLevel::Topology,
    IncidenceLevel::Constitution,
    IncidenceLevel::Full,
];

pub(super) fn level_name(level: IncidenceLevel) -> &'static str {
    match level {
        IncidenceLevel::Topology => "topology",
        IncidenceLevel::Constitution => "constitution",
        IncidenceLevel::Full => "full",
    }
}

fn atom(element: Element) -> AtomForm {
    AtomForm::from_element(element)
}

fn bond(first: u32, second: u32, order: u8) -> (AtomId, AtomId, BondForm) {
    (AtomId(first), AtomId(second), BondForm::from_order(order))
}

fn ligand(atom: u32) -> StereoLigand {
    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
}

fn implicit_hydrogen(site: u32) -> StereoLigand {
    StereoLigand::new(AtomId(site), StereoLigandKind::ImplicitHydrogen)
}

fn carbon_graph(atom_count: usize, edges: &[(u32, u32)]) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); atom_count],
        bonds: edges
            .iter()
            .map(|&(first, second)| bond(first, second, 1))
            .collect(),
        ..Default::default()
    })
}

fn ordinary_naphthalene() -> Molecule {
    carbon_graph(
        10,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 4),
        ],
    )
}

fn disconnected_rings() -> Molecule {
    carbon_graph(
        12,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ],
    )
}

fn overlay_heavy() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::O,
            Element::F,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 2),
            bond(2, 3, 1),
            bond(3, 0, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(2, 7, 1),
        ],
        dative: vec![(
            vec![AtomId(4), AtomId(5)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(4), AtomId(5)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            AtomId(6),
            AtomId(7),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(0), ligand(2), ligand(4), ligand(5)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![ligand(0), ligand(4), ligand(3), ligand(6)],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

fn tetrahedral_stereo() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(atom)
            .collect(),
        bonds: vec![bond(0, 1, 1), bond(0, 2, 1), bond(0, 3, 1), bond(0, 4, 1)],
        stereo_atoms: vec![(
            AtomId(0),
            vec![ligand(1), ligand(2), ligand(3), ligand(4)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

fn meso_dichlorobutane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(0, 2, 1),
            bond(1, 3, 1),
            bond(0, 4, 1),
            bond(1, 5, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![ligand(1), ligand(2), ligand(4), implicit_hydrogen(0)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![ligand(0), ligand(3), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

fn para_stereo_trichloropentane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(3, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(3, 7, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(1),
                vec![ligand(0), ligand(2), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(2),
                vec![ligand(1), ligand(3), ligand(6), implicit_hydrogen(2)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(3),
                vec![ligand(2), ligand(4), ligand(7), implicit_hydrogen(3)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

pub(super) fn corpus() -> [CorpusCase; 6] {
    [
        CorpusCase {
            name: "ordinary_naphthalene",
            molecule: ordinary_naphthalene(),
        },
        CorpusCase {
            name: "disconnected_rings",
            molecule: disconnected_rings(),
        },
        CorpusCase {
            name: "overlay_heavy",
            molecule: overlay_heavy(),
        },
        CorpusCase {
            name: "tetrahedral_stereo",
            molecule: tetrahedral_stereo(),
        },
        CorpusCase {
            name: "meso_dichlorobutane",
            molecule: meso_dichlorobutane(),
        },
        CorpusCase {
            name: "para_stereo_trichloropentane",
            molecule: para_stereo_trichloropentane(),
        },
    ]
}
