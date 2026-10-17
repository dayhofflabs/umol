//! Sequential reaction composition (A;B): the reaction whose application equals applying A then
//! B. Per overlap of A's product `R_A` with B's reactant `L_B`, the composite is built in one
//! id space and `normalize`d; overlaps with no `B.apply(A.apply(H))` witness (the DPO gluing
//! conditions) are rejected.

use std::ops::ControlFlow;

use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm, Correspondence, EdgeId, EmbeddingKind, GraphCorrespondence,
    NodeId,
};

use super::correspondence::MoleculeCorrespondence;
use super::id::{AtomId, BondId, StereoAtomId, StereoBondId};
use super::reaction::Reaction;
use super::traits::{Lattice, Normalize};

impl Reaction {
    /// Sequential composites of `self` (A) then `other` (B): one per admissible overlap of A's
    /// product with B's reactant, including the empty overlap.
    /// `compose(A,B).apply(H)` equals `B.apply(A.apply(H))`.
    pub fn compose(
        &self,
        other: &Reaction,
        algorithm: CommonSubgraphEnumerationAlgorithm,
    ) -> Vec<Reaction> {
        compose_all(self, other, algorithm).unwrap_or_default()
    }
}

/// The sequential composite for one overlap of A's product with `b`'s reactant: glue `a_inverse.lhs`
/// (`R_A`) and `b.lhs` (`L_B`) over `overlap`, apply `A⁻¹` and `B` at the glue, and read off the
/// composite `L_c → R_c` from the two application witnesses. `None` when the overlap is
/// inadmissible (`⊥` meet) or either application dangles.
fn compose_overlap(
    a_inverse: &Reaction,
    b: &Reaction,
    overlap: &GraphCorrespondence,
) -> Option<Reaction> {
    let (glue, correspondence) = a_inverse.lhs().tracked_meet_pushout(b.lhs(), overlap)?;
    let (lhs, left) = a_inverse
        .tracked_apply_at(&glue, &correspondence.left)
        .ok()??;
    let (rhs, right) = b.tracked_apply_at(&glue, &correspondence.right).ok()??;
    let mut correspondence = left
        .reverse()
        .compose(&right)
        .expect("both applications share the pushout host");
    let atom_survives = |&(left, right): &(StereoAtomId, StereoAtomId)| {
        !lhs.stereo_atom(left)
            .attributes()
            .configuration
            .kind()
            .zip(rhs.stereo_atom(right).attributes().configuration.kind())
            .is_some_and(|(left, right)| left != right)
    };
    let bond_survives = |&(left, right): &(StereoBondId, StereoBondId)| {
        !lhs.stereo_bond(left)
            .attributes()
            .configuration
            .kind()
            .zip(rhs.stereo_bond(right).attributes().configuration.kind())
            .is_some_and(|(left, right)| left != right)
    };
    if !correspondence
        .stereo_atoms()
        .matched_pairs()
        .iter()
        .all(atom_survives)
        || !correspondence
            .stereo_bonds()
            .matched_pairs()
            .iter()
            .all(bond_survives)
    {
        correspondence = MoleculeCorrespondence::new(
            correspondence.atoms().clone(),
            correspondence.bonds().clone(),
            correspondence.dative_bonds().clone(),
            correspondence.aromatic_systems().clone(),
            correspondence.multicenter_bonds().clone(),
            correspondence.noncovalent_bonds().clone(),
            Correspondence::new(
                correspondence
                    .stereo_atoms()
                    .matched_pairs()
                    .iter()
                    .copied()
                    .filter(atom_survives)
                    .collect(),
                lhs.stereo_atoms().count(),
                rhs.stereo_atoms().count(),
            )
            .expect("removing pairs preserves the partial bijection"),
            Correspondence::new(
                correspondence
                    .stereo_bonds()
                    .matched_pairs()
                    .iter()
                    .copied()
                    .filter(bond_survives)
                    .collect(),
                lhs.stereo_bonds().count(),
                rhs.stereo_bonds().count(),
            )
            .expect("removing pairs preserves the partial bijection"),
        );
    }
    let deltas = lhs.difference_to(&rhs, &correspondence)?;
    Some(Reaction::new(lhs, deltas.normalize().ok()?))
}

fn compose_all(
    a: &Reaction,
    b: &Reaction,
    algorithm: CommonSubgraphEnumerationAlgorithm,
) -> Option<Vec<Reaction>> {
    let span_a = a.to_reaction_span().ok()?;
    let r_a = span_a.rhs();
    let l_b = b.lhs();

    let mut node_match = |ra: NodeId, lb: NodeId| {
        r_a.atom(AtomId::from(ra))
            .attributes()
            .meet(l_b.atom(AtomId::from(lb)).attributes())
            .is_some()
    };
    let mut edge_match = |re: EdgeId, le: EdgeId| {
        r_a.bond(BondId::from(re))
            .attributes()
            .meet(l_b.bond(BondId::from(le)).attributes())
            .is_some()
    };
    let a_inverse = a.reverse().ok()?;

    // Every overlap of R_A with L_B — the *complete* common-subgraph visitation, not just the
    // maximal ones: each distinct (incl. partial and empty) overlap is a distinct sequential
    // composite, so completeness (`seq ⊆ composed`) requires all of them. Monomorphism (not induced):
    // an R_A bond absent in the matched L_B region stays as R_A context (the R1 case), rather than
    // forcing the overlap edge-for-edge. `compose_overlap` builds the composite by gluing over the
    // overlap.
    let mut results = Vec::new();
    let _: ControlFlow<()> = r_a.raw_graph().visit_common_subgraphs(
        l_b.raw_graph(),
        &mut node_match,
        &mut edge_match,
        EmbeddingKind::Monomorphism,
        algorithm,
        |pairs| {
            let nodes = Correspondence::new(
                pairs.to_vec(),
                r_a.raw_graph().node_count(),
                l_b.raw_graph().node_count(),
            )
            .expect("common-subgraph node pairs form a valid correspondence");
            let overlap = GraphCorrespondence::induce(r_a.raw_graph(), l_b.raw_graph(), nodes)
                .expect("a common-subgraph pairing induces a unique graph correspondence");
            if let Some(composite) = compose_overlap(&a_inverse, b, &overlap) {
                results.push(composite);
            }
            ControlFlow::Continue(())
        },
    );
    Some(results)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{
        Correspondence, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
    };

    use super::super::aromatic::AromaticSystemForm;
    use super::super::atom::AtomForm;
    use super::super::bond::BondForm;
    use super::super::constraint::Constraints;
    use super::super::delta::{
        AromaticSystemDelta, AtomDelta, BondDelta, Delta, Deltas, NoncovalentBondDelta,
        StereoAtomDelta, StereoBondDelta,
    };
    use super::super::edit::{
        AtomFieldChange, BondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
        StereoBondFieldChange,
    };
    use super::super::id::{AromaticSystemId, NoncovalentBondId, StereoAtomId, StereoBondId};
    use super::super::ligand::{StereoLigand, StereoLigandKind};
    use super::super::molecule::{Molecule, MoleculeEntries};
    use super::super::noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
    };
    use super::super::num::NumForm;
    use super::super::stereo::{
        StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
    };
    use super::super::substructure::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
    use super::*;

    const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
        match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    };

    #[rstest]
    #[case::kindless_intermediate(
        false,
        StereoAtomForm::default(),
        StereoBondForm::default(),
        StereoAtomForm::new(StereoKind::SquarePlanar, 1_u32),
        StereoBondForm::new(StereoKind::Axial, 1_u32)
    )]
    #[case::explicit_same_kind(
        true,
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        StereoBondForm::new(StereoKind::CisTrans, 0_u32)
    )]
    fn test_compose_overlap_stereo_replacement(
        #[case] explicit: bool,
        #[case] middle_atom: StereoAtomForm,
        #[case] middle_bond: StereoBondForm,
        #[case] new_atom: StereoAtomForm,
        #[case] new_bond: StereoBondForm,
    ) {
        let atom_ligands: Vec<_> = [1, 2, 3, 4]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let bond_ligands: Vec<_> = [2, 3, 4, 5]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: (0..6)
                .flat_map(|a| {
                    ((a + 1)..6).map(move |b| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                })
                .collect(),
            ..Default::default()
        };
        let old_atom = StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32);
        let old_bond = StereoBondForm::new(StereoKind::CisTrans, 0_u32);
        let lhs = Molecule::from_entries(MoleculeEntries {
            stereo_atoms: vec![(AtomId(0), atom_ligands.clone(), old_atom.clone())],
            stereo_bonds: vec![(BondId(0), bond_ligands.clone(), old_bond.clone())],
            ..entries.clone()
        });
        let remove_atom = Delta::StereoAtom(StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: atom_ligands.clone(),
            attributes: old_atom.clone(),
        });
        let remove_bond = Delta::StereoBond(StereoBondDelta::Remove {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: bond_ligands.clone(),
            attributes: old_bond.clone(),
        });
        let first = Reaction::new(
            lhs.clone(),
            if explicit {
                Deltas::from_iter([
                    remove_atom.clone(),
                    remove_bond.clone(),
                    Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(1),
                        site: AtomId(0),
                        ligands: atom_ligands.clone(),
                        attributes: middle_atom.clone(),
                    }),
                    Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(1),
                        site: BondId(0),
                        ligands: bond_ligands.clone(),
                        attributes: middle_bond.clone(),
                    }),
                ])
            } else {
                Deltas::from_iter([
                    Delta::StereoAtom(StereoAtomDelta::ModifyField {
                        id: StereoAtomId(0),
                        change: StereoAtomFieldChange::Configuration {
                            old: old_atom.configuration,
                            new: middle_atom.configuration.clone(),
                        },
                    }),
                    Delta::StereoBond(StereoBondDelta::ModifyField {
                        id: StereoBondId(0),
                        change: StereoBondFieldChange::Configuration {
                            old: old_bond.configuration,
                            new: middle_bond.configuration.clone(),
                        },
                    }),
                ])
            },
        );
        let second = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                stereo_atoms: vec![(AtomId(0), atom_ligands.clone(), middle_atom.clone())],
                stereo_bonds: vec![(BondId(0), bond_ligands.clone(), middle_bond.clone())],
                ..entries
            }),
            Deltas::from_iter([
                Delta::StereoAtom(StereoAtomDelta::ModifyField {
                    id: StereoAtomId(0),
                    change: StereoAtomFieldChange::Configuration {
                        old: middle_atom.configuration,
                        new: new_atom.configuration.clone(),
                    },
                }),
                Delta::StereoBond(StereoBondDelta::ModifyField {
                    id: StereoBondId(0),
                    change: StereoBondFieldChange::Configuration {
                        old: middle_bond.configuration,
                        new: new_bond.configuration.clone(),
                    },
                }),
            ]),
        );
        let expected = Reaction::new(
            lhs,
            Deltas::from_iter([
                remove_atom,
                remove_bond,
                Delta::StereoAtom(StereoAtomDelta::Add {
                    id: StereoAtomId(1),
                    site: AtomId(0),
                    ligands: atom_ligands,
                    attributes: new_atom,
                }),
                Delta::StereoBond(StereoBondDelta::Add {
                    id: StereoBondId(1),
                    site: BondId(0),
                    ligands: bond_ligands,
                    attributes: new_bond,
                }),
            ])
            .normalize()
            .unwrap(),
        );
        let overlap =
            GraphCorrespondence::new(Correspondence::identity(6), Correspondence::identity(15));
        assert_eq!(
            compose_overlap(&first.reverse().unwrap(), &second, &overlap),
            Some(expected)
        );
    }

    // C-O order 1→2 then 2→3; the single overlap fuses to 1→3.
    #[rstest]
    #[case::fuse(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(3) },
            })]),
        )]
    )]
    // A appends an O bonded to C (O is A-created); B raises that C-O 1→2. The composite appends the
    // O already at order 2 (create-then-modify fuses across the seam).
    #[case::created_atom(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(1),
                }),
            ]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(2),
                }),
            ]),
        )]
    )]
    // A appends a C to N (R_A = N-C); B's reactant N-C-O maps the A-created C onto the middle atom,
    // whose bond to the extra O is a boundary bond on an A-created atom — unrealizable, rejected.
    #[case::inadmissible(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::N)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(1),
                }),
            ]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ], bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(1), AtomId(2), BondForm::from_order(1)),
                ], ..Default::default() }),
            Deltas::new(),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![]
    )]
    // A raises C-N 1→2 and adds a hydrogen bond across the pair (a created overlay); B raises 2→3.
    // The full overlap fuses the bond to 1→3 and carries the noncovalent bond at id 0.
    #[case::overlay(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's lhs carries a hydrogen bond it never touches (only a covalent-order edit); B raises 2→3.
    // The composite carries the noncovalent bond (class ①) and fuses the order to 1→3.
    #[case::carry(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(3) },
            })]),
        )]
    )]
    // A removes its carried hydrogen bond. The composite carries it (class ①) and re-anchors A's
    // Remove delta onto composite noncovalent id 0.
    #[case::remove_carried(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // Both A's product and B's reactant carry the hydrogen bond on the overlap; B retypes it. The
    // overlap-region overlay corresponds (no fresh id), so B's modify re-anchors onto A's bond.
    #[case::correspondence(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // B's reactant requires a hydrogen bond on the overlap that A's product does not supply — the
    // overlap has no overlay correspondent, so it is skipped and compose yields nothing.
    #[case::required_absent(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], noncovalent: vec![([AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![]
    )]
    // A carries an aromatic system (a positioned overlay) it never touches; the composite carries it
    // (class ①, identity participants) and fuses the order.
    #[case::aromatic_carry(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(3) },
            })]),
        )]
    )]
    // A modifies its hydrogen bond and B raises the covalent order. The composite carries both
    // changes across the shared atoms.
    #[case::overlay_modify(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // A removes its hydrogen bond and B raises the covalent order. The composite carries the
    // hydrogen bond in its lhs and removes it.
    #[case::overlay_remove(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A creates a hydrogen bond and B raises the covalent order. The composite creates the
    // hydrogen bond at id 0.
    #[case::overlay_add(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's only edit removes its aromatic system (a positional kind, no atom/bond delta) — exercises
    // the aromatic carry path. The overlay removal participates in the overlap.
    #[case::aromatic_remove(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                id: AromaticSystemId(0),
                atoms: vec![AtomId(0), AtomId(1)],
                attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::AromaticSystem(AromaticSystemDelta::Remove {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
                }),
            ]),
        )]
    )]
    fn test_reaction_compose(
        #[case] a: Reaction,
        #[case] b: Reaction,
        #[case] algorithm: CommonSubgraphEnumerationAlgorithm,
        #[case] expected: Vec<Reaction>,
    ) {
        // Complete overlap enumeration returns a composite per overlap (partial and empty
        // included), so each case pins that its specific composite is *among* the results; the
        // whole set's soundness / completeness / well-formedness is the property suite's job.
        let composites = a.compose(&b, algorithm);
        for e in &expected {
            assert!(
                composites.contains(e),
                "expected composite absent from compose result: {e:?}",
            );
        }
    }

    // Disjoint reactants (C-C, N-N — no matchable atom) have only the empty overlap, so `compose`
    // is exactly the disjoint sum A ⊔ B (ids concatenated, both bond modifies relabeled:
    // B's bond 0 → 1) — an exact-set assertion.
    #[rstest]
    #[case::disjoint_sum(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::N), AtomForm::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::N),
                ], bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(2), AtomId(3), BondForm::from_order(1)),
                ], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(1),
                    change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
                }),
            ]),
        )]
    )]
    #[case::deletion_only(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            })]),
        ),
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::N).with_charge(0)],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(0),
                    new: NumForm::Lit(1),
                },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::N).with_charge(0),
                    AtomForm::from_element(Element::C),
                ],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::C),
                }),
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                }),
            ]),
        )]
    )]
    fn test_reaction_compose_disjoint(
        #[case] a: Reaction,
        #[case] b: Reaction,
        #[case] algorithm: CommonSubgraphEnumerationAlgorithm,
        #[case] expected: Vec<Reaction>,
    ) {
        assert_eq!(a.compose(&b, algorithm), expected);
    }

    #[rstest]
    fn test_reaction_compose_apply_equivalence() {
        // compose(A,B).apply(H) == B.apply(A.apply(H)): C-O 1→2 then 2→3 on host C-O order 1.
        let a = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        );
        let b = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(3),
                },
            })]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });

        let composed: Vec<Molecule> = a
            .compose(
                &b,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            )
            .iter()
            .flat_map(|c| c.apply(&host, MATCH_CONFIG).unwrap().map(Result::unwrap))
            .collect();
        let sequential: Vec<Molecule> = a
            .apply(&host, MATCH_CONFIG)
            .unwrap()
            .map(Result::unwrap)
            .flat_map(|intermediate| {
                b.apply(&intermediate, MATCH_CONFIG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .collect();

        let product = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(3))],
            ..Default::default()
        });
        assert_eq!(composed, vec![product.clone()]);
        assert_eq!(sequential, vec![product]);
    }

    // compose_overlap builds the span composite for one overlap: A (C–O order 1→2) then B (2→3), fused
    // over the shared overlap; applied at the host it reproduces B(A(host)).
    #[rstest]
    #[case::order_fuse(1, 2, 3)]
    fn test_compose_overlap(#[case] start: i64, #[case] mid: i64, #[case] end: i64) {
        let a = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(start as u8))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(start),
                    new: NumForm::Lit(mid),
                },
            })]),
        );
        let b = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(mid as u8))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(mid),
                    new: NumForm::Lit(end),
                },
            })]),
        );
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let a_inverse = a.reverse().unwrap();
        let composite = compose_overlap(&a_inverse, &b, &overlap).expect("admissible composite");

        let alg = SubgraphIsomorphismAlgorithm::Vf2;
        let host = a.lhs().clone();
        let intermediate = a
            .apply(
                &host,
                SubstructureMatchConfig {
                    subgraph_isomorphism_algorithm: alg,
                    ..MATCH_CONFIG
                },
            )
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        let sequential = b
            .apply(
                &intermediate,
                SubstructureMatchConfig {
                    subgraph_isomorphism_algorithm: alg,
                    ..MATCH_CONFIG
                },
            )
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        let composed = composite
            .apply(
                &host,
                SubstructureMatchConfig {
                    subgraph_isomorphism_algorithm: alg,
                    ..MATCH_CONFIG
                },
            )
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(composed, sequential);
    }

    // A tetrahedral center inverted by A then inverted back by B. `meet_pushout` keeps A⁻¹'s (self's)
    // ligand frame, so A⁻¹ applies at the glue untouched; B states the same center in ligand order
    // `b_ligands`, so its stereo delta must be re-framed into the glue frame before `apply_at`. The net
    // composite is a stereo no-op, reproducing `B(A(host)) = host`. `same_frame` is the control (B's
    // frame already matches the glue, no reframe); `swapped_frame` forces the reframe.
    #[rstest]
    #[case::same_frame([1, 2, 3, 4], 1, 0)]
    #[case::swapped_frame([2, 1, 3, 4], 0, 1)]
    fn test_reaction_compose_stereo(
        #[case] b_ligands: [u32; 4],
        #[case] b_old: u32,
        #[case] b_new: u32,
    ) {
        let a = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                    AtomForm::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            })]),
        );
        let b = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                    AtomForm::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    b_ligands
                        .iter()
                        .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, b_old),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(b_old),
                    ),
                    new: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(b_new),
                    ),
                },
            })]),
        );
        // Invert-then-invert is a net stereo no-op; whatever ligand frame B states its center in, the
        // reframe carries B's delta into the glue frame so the composite folds to A's reactant with no
        // deltas — a frame-invariant result.
        let expected = Reaction::new(a.lhs().clone(), Deltas::new());
        assert!(a
            .compose(
                &b,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking
            )
            .contains(&expected));
    }
}
