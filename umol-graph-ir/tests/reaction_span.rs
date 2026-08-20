use rstest::{fixture, rstest};
use umol_graph_core::{EdgeId, NodeId, RelationId};
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomForm, AtomId, BondForm, BondId, DativeBondForm, EntitySpan,
    MulticenterBondForm, NoncovalentBondForm, ReactionSpan, ReactionSpanEntries, StereoAtomForm,
    StereoBondForm, StereoLigand, StereoLigandKind,
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
            AtomId(0),
            AtomId(2),
            EntitySpan::Unchanged(NoncovalentBondForm::default()),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
            EntitySpan::Unchanged(StereoAtomForm::default()),
        )],
        stereo_bonds: vec![(
            BondId(0),
            vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
            EntitySpan::Unchanged(StereoBondForm::default()),
        )],
        constraints: Vec::new(),
    })
}

#[rstest]
fn test_reaction_span_dative_bonds(reaction_span: ReactionSpan) {
    let id = RelationId(0);

    assert_eq!(
        reaction_span.dative_bonds().participants_1(id),
        &[NodeId(0)]
    );
    assert_eq!(
        reaction_span.dative_bonds().participants_2(id),
        &[NodeId(1)]
    );
    assert_eq!(
        reaction_span.dative_bonds().data(id),
        &EntitySpan::Unchanged(DativeBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_aromatic_systems(reaction_span: ReactionSpan) {
    let id = RelationId(0);

    assert_eq!(
        reaction_span.aromatic_systems().participants(id),
        &[NodeId(0), NodeId(1), NodeId(2)],
    );
    assert_eq!(
        reaction_span.aromatic_systems().data(id),
        &EntitySpan::Unchanged(AromaticSystemForm::default()),
    );
}

#[rstest]
fn test_reaction_span_multicenter_bonds(reaction_span: ReactionSpan) {
    let id = RelationId(0);

    assert_eq!(
        reaction_span.multicenter_bonds().participants(id),
        &[NodeId(0), NodeId(1), NodeId(2)],
    );
    assert_eq!(
        reaction_span.multicenter_bonds().data(id),
        &EntitySpan::Unchanged(MulticenterBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_noncovalent_bonds(reaction_span: ReactionSpan) {
    let id = RelationId(0);

    assert_eq!(
        reaction_span.noncovalent_bonds().participants(id),
        &[NodeId(0), NodeId(2)],
    );
    assert_eq!(
        reaction_span.noncovalent_bonds().data(id),
        &EntitySpan::Unchanged(NoncovalentBondForm::default()),
    );
}

#[rstest]
fn test_reaction_span_stereo_atoms(reaction_span: ReactionSpan) {
    let id = RelationId(0);
    let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);

    assert_eq!(
        reaction_span.stereo_atoms().participants_1(id),
        &[NodeId(0)]
    );
    assert_eq!(reaction_span.stereo_atoms().participants_2(id), &[ligand]);
    assert_eq!(
        reaction_span.stereo_atoms().data(id),
        &EntitySpan::Unchanged(StereoAtomForm::default()),
    );
}

#[rstest]
fn test_reaction_span_stereo_bonds(reaction_span: ReactionSpan) {
    let id = RelationId(0);
    let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);

    assert_eq!(
        reaction_span.stereo_bonds().participants_1(id),
        &[EdgeId(0)]
    );
    assert_eq!(reaction_span.stereo_bonds().participants_2(id), &[ligand]);
    assert_eq!(
        reaction_span.stereo_bonds().data(id),
        &EntitySpan::Unchanged(StereoBondForm::default()),
    );
}
