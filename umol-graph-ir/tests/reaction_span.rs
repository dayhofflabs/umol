use rstest::{fixture, rstest};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomForm, AtomId, BondForm, BondId, DativeBondForm,
    DativeBondId, EntitySpan, MulticenterBondForm, MulticenterBondId, NoncovalentBondForm,
    NoncovalentBondId, ReactionSpan, ReactionSpanEntries, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoLigand, StereoLigandKind,
};

#[fixture]
fn reaction_span() -> ReactionSpan {
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::default()),
            EntitySpan::Unchanged(AtomForm::default()),
            EntitySpan::Unchanged(AtomForm::default()),
        ],
        bonds: vec![(
            AtomId(0),
            AtomId(1),
            EntitySpan::Unchanged(BondForm::default()),
        )],
        dative: vec![(
            vec![AtomId(1)],
            AtomId(0),
            EntitySpan::Unchanged(DativeBondForm::default()),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::default()),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(MulticenterBondForm::default()),
        )],
        noncovalent: vec![(
            [AtomId(0), AtomId(2)],
            EntitySpan::Unchanged(NoncovalentBondForm::default()),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            EntitySpan::Unchanged(StereoAtomForm::default()),
        )],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            EntitySpan::Unchanged(StereoBondForm::default()),
        )],
        constraints: Vec::new(),
    })
}

#[rstest]
fn test_reaction_span_dative_bonds(reaction_span: ReactionSpan) {
    let id = DativeBondId(0);

    assert_eq!(reaction_span.dative_bonds().acceptor(id), AtomId(0));
    assert_eq!(
        reaction_span.dative_bonds().donors(id).collect::<Vec<_>>(),
        [AtomId(1)]
    );
    assert_eq!(
        reaction_span.dative_bonds().attributes(id),
        &EntitySpan::Unchanged(DativeBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_aromatic_systems(reaction_span: ReactionSpan) {
    let id = AromaticSystemId(0);

    assert_eq!(
        reaction_span
            .aromatic_systems()
            .atoms(id)
            .collect::<Vec<_>>(),
        [AtomId(0), AtomId(1), AtomId(2)],
    );
    assert_eq!(
        reaction_span.aromatic_systems().attributes(id),
        &EntitySpan::Unchanged(AromaticSystemForm::default()),
    );
}

#[rstest]
fn test_reaction_span_multicenter_bonds(reaction_span: ReactionSpan) {
    let id = MulticenterBondId(0);

    assert_eq!(
        reaction_span
            .multicenter_bonds()
            .atoms(id)
            .collect::<Vec<_>>(),
        [AtomId(0), AtomId(1), AtomId(2)],
    );
    assert_eq!(
        reaction_span.multicenter_bonds().attributes(id),
        &EntitySpan::Unchanged(MulticenterBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_noncovalent_bonds(reaction_span: ReactionSpan) {
    let id = NoncovalentBondId(0);

    assert_eq!(
        reaction_span.noncovalent_bonds().atoms(id),
        [AtomId(0), AtomId(2)],
    );
    assert_eq!(
        reaction_span.noncovalent_bonds().attributes(id),
        &EntitySpan::Unchanged(NoncovalentBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_stereo_atoms(reaction_span: ReactionSpan) {
    let id = StereoAtomId(0);
    let ligand = StereoLigand::new(AtomId(1), StereoLigandKind::Atom);

    assert_eq!(reaction_span.stereo_atoms().site(id), AtomId(0));
    assert_eq!(reaction_span.stereo_atoms().ligands(id), [ligand]);
    assert_eq!(
        reaction_span.stereo_atoms().attributes(id),
        &EntitySpan::Unchanged(StereoAtomForm::default()),
    );
}

#[rstest]
fn test_reaction_span_stereo_bonds(reaction_span: ReactionSpan) {
    let id = StereoBondId(0);
    let ligands = [
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
    ];

    assert_eq!(reaction_span.stereo_bonds().site(id), BondId(0));
    assert_eq!(reaction_span.stereo_bonds().ligands(id), ligands);
    assert_eq!(
        reaction_span.stereo_bonds().attributes(id),
        &EntitySpan::Unchanged(StereoBondForm::default()),
    );
}
