//! Sequential reaction composition (A;B): the reaction whose application equals applying A then
//! B. Per overlap of A's product `R_A` with B's reactant `L_B`, the composite is built in one
//! id space and `canonicalize`d; overlaps with no `B.apply(A.apply(H))` witness (the DPO gluing
//! conditions) are rejected.

use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm, EdgeId, EmbeddingKind, GraphCorrespondence, NodeId,
};

use super::id::{AtomId, BondId};
use super::reaction::ReactionAst;
use super::traits::{Canonicalize, Lattice};

impl ReactionAst {
    /// Sequential composites of `self` (A) then `other` (B): one per admissible overlap of A's
    /// product with B's reactant, including the empty overlap.
    /// `compose(A,B).apply(H)` equals `B.apply(A.apply(H))`.
    pub fn compose(
        &self,
        other: &ReactionAst,
        algorithm: CommonSubgraphEnumerationAlgorithm,
    ) -> Vec<ReactionAst> {
        compose_all(self, other, algorithm).unwrap_or_default()
    }
}

/// The sequential composite for one overlap of A's product with `b`'s reactant: glue `a_inverse.lhs`
/// (`R_A`) and `b.lhs` (`L_B`) over `overlap`, apply `A⁻¹` and `B` at the glue, and read off the
/// composite `L_c → R_c` (`from_sides` diffs the two applied glues). `None` when the overlap is
/// inadmissible (`⊥` meet) or either application dangles.
fn compose_overlap(
    a_inverse: &ReactionAst,
    b: &ReactionAst,
    overlap: &GraphCorrespondence,
) -> Option<ReactionAst> {
    let glue = a_inverse.lhs.meet_pushout(&b.lhs, overlap)?;
    let derivation_a = a_inverse.apply_at(&glue.object, &glue.left).ok()?;
    let derivation_b = b.apply_at(&glue.object, &glue.right).ok()?;
    let correspondence = derivation_a
        .atom_map()
        .reverse()
        .compose(derivation_b.atom_map());
    let composite = ReactionAst::from_sides(
        derivation_a.rhs().clone(),
        derivation_b.rhs().clone(),
        correspondence,
    );
    Some(ReactionAst::new(
        composite.lhs,
        composite.deltas.canonicalize().ok()?,
    ))
}

fn compose_all(
    a: &ReactionAst,
    b: &ReactionAst,
    algorithm: CommonSubgraphEnumerationAlgorithm,
) -> Option<Vec<ReactionAst>> {
    let span_a = a.to_reaction_span().ok()?;
    let r_a = span_a.rhs();
    let l_b = &b.lhs;

    let mut node_match = |ra: NodeId, lb: NodeId| {
        r_a.atom(AtomId::from(ra))
            .ast
            .meet(l_b.atom(AtomId::from(lb)).ast)
            .is_some()
    };
    let mut edge_match = |re: EdgeId, le: EdgeId| {
        r_a.bond(BondId::from(re))
            .ast
            .meet(l_b.bond(BondId::from(le)).ast)
            .is_some()
    };
    // Every overlap of R_A with L_B — the *complete* common-subgraph enumeration, not just the
    // maximal ones: each distinct (incl. partial and empty) overlap is a distinct sequential
    // composite, so completeness (`seq ⊆ composed`) requires all of them. Monomorphism (not induced):
    // an R_A bond absent in the matched L_B region stays as R_A context (the R1 case), rather than
    // forcing the overlap edge-for-edge. `compose_overlap` builds the composite by gluing over the
    // overlap.
    let overlaps = r_a.raw_graph().enumerate_common_subgraphs(
        l_b.raw_graph(),
        &mut node_match,
        &mut edge_match,
        EmbeddingKind::Monomorphism,
        algorithm,
    );
    let a_inverse = a.reverse().ok()?;

    let mut results = Vec::new();
    for overlap in overlaps {
        if let Some(composite) = compose_overlap(&a_inverse, b, &overlap) {
            results.push(composite);
        }
    }
    Some(results)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{Correspondence, SubgraphIsomorphismAlgorithm};

    use super::super::aromatic::AromaticSystemAst;
    use super::super::atom::AtomAst;
    use super::super::bond::BondAst;
    use super::super::constraint::Constraints;
    use super::super::delta::{
        AromaticSystemDelta, AtomDelta, BondDelta, Delta, Deltas, NoncovalentBondDelta,
        StereoAtomDelta,
    };
    use super::super::edit::{
        AtomFieldChange, BondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    };
    use super::super::id::{AromaticSystemId, NoncovalentBondId, StereoAtomId};
    use super::super::ligand::{StereoLigand, StereoLigandKind};
    use super::super::molecule::{MoleculeAst, MoleculeParts};
    use super::super::noncovalent::{
        NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst,
    };
    use super::super::stereo::{StereoAtomAst, StereoConfigurationAst, StereoCoset, StereoKind};
    use super::super::substructure::SubstructureMatchAlgorithm;
    use super::super::value::ValueAst;
    use super::*;

    // C-O order 1→2 then 2→3; the single overlap fuses to 1→3.
    #[rstest]
    #[case::fuse(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A appends an O bonded to C (O is A-created); B raises that C-O 1→2. The composite appends the
    // O already at order 2 (create-then-modify fuses across the seam).
    #[case::created_atom(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(2),
                }),
            ]),
        )]
    )]
    // A appends a C to N (R_A = N-C); B's reactant N-C-O maps the A-created C onto the middle atom,
    // whose bond to the extra O is a boundary bond on an A-created atom — unrealizable, rejected.
    #[case::inadmissible(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::N)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ], bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(1), AtomId(2), BondAst::from_order(1)),
                ], ..Default::default() }),
            Deltas::new(),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![]
    )]
    // A raises C-N 1→2 and adds a hydrogen bond across the pair (a created overlay); B raises 2→3.
    // The full overlap fuses the bond to 1→3 and carries the noncovalent bond at id 0.
    #[case::overlay(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's lhs carries a hydrogen bond it never touches (only a covalent-order edit); B raises 2→3.
    // The composite carries the noncovalent bond (class ①) and fuses the order to 1→3.
    #[case::carry(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A removes its carried hydrogen bond. The composite carries it (class ①) and re-anchors A's
    // Remove delta onto composite noncovalent id 0.
    #[case::remove_carried(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // Both A's product and B's reactant carry the hydrogen bond on the overlap; B retypes it. The
    // overlap-region overlay corresponds (no fresh id), so B's modify re-anchors onto A's bond.
    #[case::correspondence(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // B's reactant requires a hydrogen bond on the overlap that A's product does not supply — the
    // overlap has no overlay correspondent, so it is skipped and compose yields nothing.
    #[case::required_absent(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], noncovalent: vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![]
    )]
    // A carries an aromatic system (a positional family) it never touches; the composite carries it
    // (class ①, identity participants) and fuses the order.
    #[case::aromatic_carry(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A modifies its hydrogen bond and B raises the covalent order. The composite carries both
    // changes across the shared atoms.
    #[case::overlay_modify(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // A removes its hydrogen bond and B raises the covalent order. The composite carries the
    // hydrogen bond in its lhs and removes it.
    #[case::overlay_remove(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A creates a hydrogen bond and B raises the covalent order. The composite creates the
    // hydrogen bond at id 0.
    #[case::overlay_add(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's only edit removes its aromatic system (a positional kind, no atom/bond delta) — exercises
    // the aromatic carry path. The overlay removal participates in the overlap.
    #[case::aromatic_remove(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                id: AromaticSystemId(0),
                atoms: vec![AtomId(0), AtomId(1)],
                ast: AromaticSystemAst::from_electrons(vec![1, 2]),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::AromaticSystem(AromaticSystemDelta::Remove {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    ast: AromaticSystemAst::from_electrons(vec![1, 2]),
                }),
            ]),
        )]
    )]
    fn test_reaction_ast_compose(
        #[case] a: ReactionAst,
        #[case] b: ReactionAst,
        #[case] algorithm: CommonSubgraphEnumerationAlgorithm,
        #[case] expected: Vec<ReactionAst>,
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
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::N)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::N),
                ], bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(2), AtomId(3), BondAst::from_order(1)),
                ], ..Default::default() }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(1),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
            ]),
        )]
    )]
    #[case::deletion_only(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C)],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::N).with_charge(0)],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(0),
                    new: ValueAst::Lit(1),
                },
            })]),
        ),
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::N).with_charge(0),
                    AtomAst::from_element(Element::C),
                ],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: ValueAst::Lit(0),
                        new: ValueAst::Lit(1),
                    },
                }),
            ]),
        )]
    )]
    fn test_reaction_ast_compose_disjoint(
        #[case] a: ReactionAst,
        #[case] b: ReactionAst,
        #[case] algorithm: CommonSubgraphEnumerationAlgorithm,
        #[case] expected: Vec<ReactionAst>,
    ) {
        assert_eq!(a.compose(&b, algorithm), expected);
    }

    #[rstest]
    fn test_reaction_ast_compose_apply_equivalence() {
        // compose(A,B).apply(H) == B.apply(A.apply(H)): C-O 1→2 then 2→3 on host C-O order 1.
        let a = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
        let b = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(2),
                    new: ValueAst::Lit(3),
                },
            })]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });

        let composed: Vec<MoleculeAst> = a
            .compose(
                &b,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            )
            .iter()
            .flat_map(|c| {
                c.apply(
                    &host,
                    SubstructureMatchAlgorithm::GraphAndOverlays,
                    SubgraphIsomorphismAlgorithm::Vf2,
                )
                .unwrap()
                .map(Result::unwrap)
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let sequential: Vec<MoleculeAst> = a
            .apply(
                &host,
                SubstructureMatchAlgorithm::GraphAndOverlays,
                SubgraphIsomorphismAlgorithm::Vf2,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .flat_map(|intermediate| {
                b.apply(
                    &intermediate,
                    SubstructureMatchAlgorithm::GraphAndOverlays,
                    SubgraphIsomorphismAlgorithm::Vf2,
                )
                .unwrap()
                .map(Result::unwrap)
                .map(|derivation| derivation.rhs().clone())
                .collect::<Vec<_>>()
            })
            .collect();

        let product = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(3))],
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
        let a = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(start as u8))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(start),
                    new: ValueAst::Lit(mid),
                },
            })]),
        );
        let b = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(mid as u8))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(mid),
                    new: ValueAst::Lit(end),
                },
            })]),
        );
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1),
        );
        let a_inverse = a.reverse().unwrap();
        let composite = compose_overlap(&a_inverse, &b, &overlap).expect("admissible composite");

        let alg = SubgraphIsomorphismAlgorithm::Vf2;
        let host = a.lhs.clone();
        let intermediate = a
            .apply(&host, SubstructureMatchAlgorithm::GraphAndOverlays, alg)
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .rhs()
            .clone();
        let sequential = b
            .apply(
                &intermediate,
                SubstructureMatchAlgorithm::GraphAndOverlays,
                alg,
            )
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .rhs()
            .clone();
        let composed = composite
            .apply(&host, SubstructureMatchAlgorithm::GraphAndOverlays, alg)
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .rhs()
            .clone();
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
    fn test_reaction_ast_compose_stereo(
        #[case] b_ligands: [u32; 4],
        #[case] b_old: u32,
        #[case] b_new: u32,
    ) {
        let a = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::F),
                    AtomAst::from_element(Element::Cl),
                    AtomAst::from_element(Element::Br),
                    AtomAst::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(0), AtomId(2), BondAst::from_order(1)),
                    (AtomId(0), AtomId(3), BondAst::from_order(1)),
                    (AtomId(0), AtomId(4), BondAst::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            })]),
        );
        let b = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::F),
                    AtomAst::from_element(Element::Cl),
                    AtomAst::from_element(Element::Br),
                    AtomAst::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(0), AtomId(2), BondAst::from_order(1)),
                    (AtomId(0), AtomId(3), BondAst::from_order(1)),
                    (AtomId(0), AtomId(4), BondAst::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    b_ligands
                        .iter()
                        .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                        .collect(),
                    StereoAtomAst::new(StereoKind::Tetrahedral, b_old),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(b_old),
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(b_new),
                    ),
                },
            })]),
        );
        // Invert-then-invert is a net stereo no-op; whatever ligand frame B states its center in, the
        // reframe carries B's delta into the glue frame so the composite folds to A's reactant with no
        // deltas — a frame-invariant result.
        let expected = ReactionAst::new(a.lhs.clone(), Deltas::new());
        assert!(a
            .compose(
                &b,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking
            )
            .contains(&expected));
    }
}
