//! The attributed pushout of two molecules: glue `self` and `other` on a shared subgraph, meeting
//! atom / bond / overlay data and combining molecule constraints where they coincide — the umol-graph-ir
//! attribute layer over graph-core's structural `pushout`. A child of `molecule` so it reaches the
//! private overlay relation-sets directly, without exposing raw accessors.

use umol_graph_core::{
    Correspondence, EdgeId, FixedVarBirelationSet, GraphCorrespondence, NodeId, Ordered,
    RelationData, RelationParticipant, Remapping, Unordered, VarRelationSet,
};

use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::Constraints;
use super::super::correspondence::MoleculeCorrespondence;
use super::super::dative::DativeBondAst;
use super::super::id::{AtomId, BondId};
use super::super::ligand::StereoLigand;
use super::super::noncovalent::NoncovalentBondAst;
use super::super::remap::IdRemapping;
use super::super::stereo::{StereoAtomAst, StereoBondAst};
use super::super::traits::Lattice;
use super::{MoleculeAst, MoleculeEntries};

/// The attributed pushout of two molecules over a graph `overlap`: `self` and `other` glued on their
/// shared subgraph, with atom / bond data `meet`-combined where they coincide. `object` keeps `self`'s
/// ids; `left` / `right` embed each side into it.
#[allow(dead_code)]
pub struct MoleculePushout {
    pub object: MoleculeAst,
    pub left: MoleculeCorrespondence,
    pub right: MoleculeCorrespondence,
}

impl MoleculeAst {
    /// Glue `self` (left) and `other` (right) over `overlap` (a common subgraph — its edges the
    /// coincident bonds), meeting atom / bond / overlay data at coincident entities and combining the
    /// two molecule-constraint sets; `None` when any coincident `meet` is `⊥` (the overlap is
    /// inadmissible). Stereo overlays keep `self`'s ligand frame; `other`'s coincident cosets are
    /// aligned to it (`transform_frame`) before the pushout, so they `meet` in a shared frame.
    #[allow(dead_code)]
    pub fn meet_pushout(
        &self,
        other: &MoleculeAst,
        overlap: &GraphCorrespondence,
    ) -> Option<MoleculePushout> {
        let po = self.raw_graph().pushout(other.raw_graph(), overlap);

        let mut atoms: Vec<AtomAst> = Vec::with_capacity(po.object.node_count());
        for node in 0..po.object.node_count() as u32 {
            let object = NodeId(node);
            let atom = match (
                po.left.nodes().left_of(object),
                po.right.nodes().left_of(object),
            ) {
                (Some(l), Some(r)) => self
                    .atom(AtomId::from(l))
                    .ast
                    .meet(other.atom(AtomId::from(r)).ast)?,
                (Some(l), None) => self.atom(AtomId::from(l)).ast.clone(),
                (None, Some(r)) => other.atom(AtomId::from(r)).ast.clone(),
                (None, None) => unreachable!("a glued node originates from a side"),
            };
            atoms.push(atom);
        }

        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::with_capacity(po.object.edge_count());
        for edge in 0..po.object.edge_count() as u32 {
            let object = EdgeId(edge);
            let [u, v] = po.object.edge_endpoints(object);
            let bond = match (
                po.left.edges().left_of(object),
                po.right.edges().left_of(object),
            ) {
                (Some(l), Some(r)) => self
                    .bond(BondId::from(l))
                    .ast
                    .meet(other.bond(BondId::from(r)).ast)?,
                (Some(l), None) => self.bond(BondId::from(l)).ast.clone(),
                (None, Some(r)) => other.bond(BondId::from(r)).ast.clone(),
                (None, None) => unreachable!("a glued edge originates from a side"),
            };
            bonds.push((AtomId::from(u), AtomId::from(v), bond));
        }

        // Overlays glue over the same pushout: relabel `other`'s participants into the glue space
        // (`self` already keeps its ids), then merge coinciding overlays by `meet`; non-coinciding
        // ones are appended (context). `⊥` on any coincident meet makes the whole glue inadmissible.
        let participant_remapping = Remapping::new(
            (0..other.raw_graph().node_count())
                .map(|index| {
                    po.right
                        .nodes()
                        .right_of(NodeId(index as u32))
                        .expect("right total on other nodes")
                })
                .collect(),
            (0..other.raw_graph().edge_count())
                .map(|index| {
                    po.right
                        .edges()
                        .right_of(EdgeId(index as u32))
                        .expect("right total on other edges")
                })
                .collect(),
        );

        let aromatic = glue_var_overlays(
            &self.aromatic_systems,
            &other.aromatic_systems,
            &participant_remapping,
        )?;
        let multicenter = glue_var_overlays(
            &self.multicenter_bonds,
            &other.multicenter_bonds,
            &participant_remapping,
        )?;

        let dative_glue = other.dative_bonds.apply_remapping(&participant_remapping);
        let dative_merged = self.dative_bonds.pushout(&dative_glue, |a, b| a.meet(b))?;
        let dative_object = &dative_merged.object;
        let dative: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = dative_object
            .relation_ids()
            .map(|id| {
                (
                    dative_object
                        .participants_2(id)
                        .iter()
                        .map(|&n| AtomId::from(n))
                        .collect(),
                    AtomId::from(dative_object.participants_1(id)[0]),
                    dative_object.data(id).clone(),
                )
            })
            .collect();

        let noncovalent_glue = other
            .noncovalent_bonds
            .apply_remapping(&participant_remapping);
        let noncovalent_merged = self
            .noncovalent_bonds
            .pushout(&noncovalent_glue, |a, b| a.meet(b))?;
        let noncovalent_object = &noncovalent_merged.object;
        let noncovalent: Vec<(AtomId, AtomId, NoncovalentBondAst)> = noncovalent_object
            .relation_ids()
            .map(|id| {
                let &[u, v] = noncovalent_object.participants(id);
                (
                    AtomId::from(u),
                    AtomId::from(v),
                    noncovalent_object.data(id).clone(),
                )
            })
            .collect();

        // Stereo overlays differ: ligand order is meaningful (the coset is frame-relative), but a pure
        // id-remap preserves the sequence and leaves the coset untouched. So the glue keeps `self`'s
        // ligand frame as-is; `other`'s overlays are relabeled into the glue space and, where one
        // coincides with a `self` site (same site + ligand multiset), its coset is aligned to `self`'s
        // frame (`transform_frame`) before the full-participant relation `pushout` `meet`s the two
        // (`⊥ → None`). `other`-only sites keep their own (relabeled) frame. A same-site/different-ligand
        // collision leaves two overlays on one site — over-coordination, rejected by the `has_conflict`
        // gate below.
        let remapped_stereo_atoms = other.stereo_atoms.apply_remapping(&participant_remapping);

        let stereo_atom_right = FixedVarBirelationSet::new(stereo_glue_entries(
            &self.stereo_atoms,
            &remapped_stereo_atoms,
            |d, before, after| d.transform_frame(before, after),
        )?);
        let stereo_atom_merged = self
            .stereo_atoms
            .pushout(&stereo_atom_right, |a, b| a.meet(b))?;
        let sa_object = &stereo_atom_merged.object;
        let stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)> = sa_object
            .relation_ids()
            .map(|id| {
                (
                    AtomId::from(sa_object.participants_1(id)[0]),
                    sa_object.participants_2(id).to_vec(),
                    sa_object.data(id).clone(),
                )
            })
            .collect();

        let remapped_stereo_bonds = other.stereo_bonds.apply_remapping(&participant_remapping);
        let stereo_bond_right = FixedVarBirelationSet::new(stereo_glue_entries(
            &self.stereo_bonds,
            &remapped_stereo_bonds,
            |d, before, after| d.transform_frame(before, after),
        )?);
        let stereo_bond_merged = self
            .stereo_bonds
            .pushout(&stereo_bond_right, |a, b| a.meet(b))?;
        let sb_object = &stereo_bond_merged.object;
        let stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> = sb_object
            .relation_ids()
            .map(|id| {
                (
                    BondId::from(sb_object.participants_1(id)[0]),
                    sb_object.participants_2(id).to_vec(),
                    sb_object.data(id).clone(),
                )
            })
            .collect();

        let mut object = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints: Constraints::new(),
        });

        let atom_correspondence = |nodes: &Correspondence<NodeId>| {
            Correspondence::new(
                nodes
                    .matched_pairs()
                    .iter()
                    .map(|&(left, right)| (AtomId::from(left), AtomId::from(right)))
                    .collect(),
                nodes.left_count(),
                nodes.right_count(),
            )
            .expect("graph pushout preserves atom correspondence invariants")
        };
        let left =
            MoleculeCorrespondence::induce(self, &object, atom_correspondence(po.left.nodes()))?;
        let right =
            MoleculeCorrespondence::induce(other, &object, atom_correspondence(po.right.nodes()))?;

        // Molecule-level constraints: `self`'s hold in the glue as-is (it keeps `self`'s ids); `other`'s
        // are re-anchored through the `right` embedding. Conjunction, deduplicated.
        let remapping = IdRemapping::new(
            right.atoms().matched_pairs().iter().copied().collect(),
            right.bonds().matched_pairs().iter().copied().collect(),
            right
                .dative_bonds()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            right
                .aromatic_systems()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            right
                .multicenter_bonds()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            right
                .noncovalent_bonds()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            right
                .stereo_atoms()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            right
                .stereo_bonds()
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
        );
        let mut constraints = self.constraints.clone();
        for c in other.constraints.iter() {
            let remapped = c.clone().remap(&remapping);
            if !constraints.iter().any(|existing| existing == &remapped) {
                constraints.push(remapped);
            }
        }
        object.constraints = constraints;

        // Emit-compliance: the glue is a generated molecule, so it must satisfy every per-entity
        // structural invariant — gluing can collide bonds/overlays a well-formed input never would
        // (parallel bonds, overlapping systems, two stereo centers on one site, …). The per-entity
        // `has_conflict` primitives are the shared gates (also consulted by the validator and
        // `apply_at`); enforced per generating op pending a single central emit gate.
        if object.bonds().has_conflict()
            || object.dative_bonds().has_conflict()
            || object.aromatic_systems().has_conflict()
            || object.multicenter_bonds().has_conflict()
            || object.noncovalent_bonds().has_conflict()
            || object.stereo_atoms().has_conflict()
            || object.stereo_bonds().has_conflict()
        {
            return None;
        }

        Some(MoleculePushout {
            object,
            left,
            right,
        })
    }
}

/// Glue two `VarRelationSet` overlay families after relabeling `right` into the glue id space.
/// Coinciding overlays merge by `meet`; non-coinciding overlays are appended as context. `None` if
/// any coincident meet is `⊥`.
fn glue_var_overlays<D: Lattice + RelationData>(
    left: &VarRelationSet<NodeId, Unordered, D>,
    right: &VarRelationSet<NodeId, Unordered, D>,
    remapping: &Remapping,
) -> Option<Vec<(Vec<AtomId>, D)>> {
    let right_glue = right.apply_remapping(remapping);
    let merged = left.pushout(&right_glue, |a, b| a.meet(b))?;
    Some(
        merged
            .object
            .relation_ids()
            .map(|id| {
                (
                    merged
                        .object
                        .participants(id)
                        .iter()
                        .map(|&n| AtomId::from(n))
                        .collect(),
                    merged.object.data(id).clone(),
                )
            })
            .collect(),
    )
}

/// Align already-remapped right-side stereo entries to coincident left-side ligand frames, carrying
/// each coset through `transform_frame`. Right-only entries retain their remapped frame. A right
/// entry whose site collides with a left entry under a different ligand set remains distinct;
/// `meet_pushout` rejects the resulting over-coordination through its conflict gate.
#[allow(clippy::type_complexity)]
fn stereo_glue_entries<S, D>(
    left: &FixedVarBirelationSet<S, Ordered, 1, StereoLigand, Ordered, D>,
    right: &FixedVarBirelationSet<S, Ordered, 1, StereoLigand, Ordered, D>,
    transform: impl Fn(&D, &[StereoLigand], &[StereoLigand]) -> Option<D>,
) -> Option<Vec<([S; 1], Vec<StereoLigand>, D)>>
where
    S: RelationParticipant,
    D: Clone,
{
    right
        .relation_ids()
        .map(|id| {
            let site = right.participants_1(id)[0];
            let remapped_frame = right.participants_2(id).to_vec();
            match left.find_by_participants(&[site], &remapped_frame) {
                Some(hit) => {
                    let target = left.participants_2(hit).to_vec();
                    let data = transform(right.data(id), &remapped_frame, &target)?;
                    Some(([site], target, data))
                }
                None => Some(([site], remapped_frame, right.data(id).clone())),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::Correspondence;

    use super::super::super::aromatic::AromaticSystemAst;
    use super::super::super::constraint::{AtomConstraintAst, Constraint};
    use super::super::super::ligand::StereoLigandKind;
    use super::super::super::multicenter::MulticenterBondAst;
    use super::super::super::stereo::StereoKind;
    use super::*;

    // A single shared atom (node 0 ↔ node 0), no shared bond; each side has one unmatched atom.
    #[fixture]
    fn overlap() -> GraphCorrespondence {
        GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        )
    }

    // meet_pushout glues over the shared atom: `object` keeps left's ids, appends right's unmatched
    // atom, and the shared atom carries the two sides' meet (either order → the more specific).
    #[rstest]
    #[case::top_meets_element(AtomAst::default(), AtomAst::from_element(Element::C), Element::C)]
    #[case::element_meets_top(AtomAst::from_element(Element::C), AtomAst::default(), Element::C)]
    fn test_molecule_ast_meet_pushout(
        overlap: GraphCorrespondence,
        #[case] left_shared: AtomAst,
        #[case] right_shared: AtomAst,
        #[case] shared_element: Element,
    ) {
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![left_shared, AtomAst::from_element(Element::N)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![right_shared, AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(shared_element),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
            ],
            ..Default::default()
        });
        assert_eq!(
            left.meet_pushout(&right, &overlap)
                .expect("admissible glue")
                .object,
            expected,
        );
    }

    // The glue is inadmissible (`None`) when it would be malformed: a coincident-atom meet is `⊥`
    // (`carbon_nitrogen` / `oxygen_nitrogen`), or an emit-compliance invariant fails — here two aromatic
    // systems that share glue atom 0 (`[0,1]` from left, `[0,2]` from right's unmatched atom), which
    // the
    // `has_conflict` gate rejects.
    #[rstest]
    #[case::carbon_nitrogen(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::oxygen_nitrogen(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::N)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::aromatic_overlap(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default())],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default())],
            constraints: Constraints::new(),
            ..Default::default()
        }),
    )]
    fn test_molecule_ast_meet_pushout_inadmissible(
        overlap: GraphCorrespondence,
        #[case] left: MoleculeAst,
        #[case] right: MoleculeAst,
    ) {
        assert!(left.meet_pushout(&right, &overlap).is_none());
    }

    // The crossing node correspondence reverses each right-side participant pair. Electron counts
    // follow their atoms into the canonical glue order; coincident overlays meet and disjoint
    // overlays remain as context.
    #[rstest]
    fn test_molecule_ast_meet_pushout_overlays() {
        let full_overlap = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(1)),
                    (NodeId(1), NodeId(0)),
                    (NodeId(2), NodeId(3)),
                    (NodeId(3), NodeId(2)),
                ],
                4,
                4,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 4],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemAst::from_electrons(vec![1, 2]),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1)],
                MulticenterBondAst::from_electrons(vec![7, 11]),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 4],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemAst::from_electrons(vec![2, 1]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    AromaticSystemAst::from_electrons(vec![5, 3]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    MulticenterBondAst::from_electrons(vec![11, 7]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    MulticenterBondAst::from_electrons(vec![17, 13]),
                ),
            ],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 4],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemAst::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    AromaticSystemAst::from_electrons(vec![3, 5]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    MulticenterBondAst::from_electrons(vec![7, 11]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    MulticenterBondAst::from_electrons(vec![13, 17]),
                ),
            ],
            constraints: Constraints::new(),
            ..Default::default()
        });
        assert_eq!(
            left.meet_pushout(&right, &full_overlap)
                .expect("admissible")
                .object,
            expected,
        );
    }

    // Molecule-level constraints are carried and re-anchored: `self`'s stay put; `other`'s are remapped
    // through the embedding — `other`'s atom 1 becomes glue atom 2.
    #[rstest]
    #[case::atom_valences(4, 2)]
    fn test_molecule_ast_meet_pushout_constraints(
        overlap: GraphCorrespondence,
        #[case] left_valence: i64,
        #[case] right_valence: i64,
    ) {
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            constraints: Constraints::from(vec![Constraint::Atom(
                AtomId(0),
                AtomConstraintAst::valence(left_valence),
            )]),
            ..Default::default()
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            constraints: Constraints::from(vec![Constraint::Atom(
                AtomId(1),
                AtomConstraintAst::valence(right_valence),
            )]),
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
            ],
            constraints: Constraints::from(vec![
                Constraint::Atom(AtomId(0), AtomConstraintAst::valence(left_valence)),
                Constraint::Atom(AtomId(2), AtomConstraintAst::valence(right_valence)),
            ]),
            ..Default::default()
        });
        assert_eq!(
            left.meet_pushout(&right, &overlap)
                .expect("admissible")
                .object,
            expected,
        );
    }

    // One tetrahedral center on atom 0 (ligands {1,2,3,4}) shared by both molecules, but `other`'s
    // ligand frame swaps 1 and 2 — a transposition, which flips the coset. `meet_pushout` canonicalizes
    // both frames (`transform_frame`) before the meet, so `agree` (opposite coset in the swapped frame
    // = same physical configuration) folds to `self`, while `contradict` (same coset = opposite
    // configuration) is `⊥`.
    #[rstest]
    #[case::agree([2, 1, 3, 4], 1, true)]
    #[case::contradict([2, 1, 3, 4], 0, false)]
    #[case::ligand_set([2, 1, 3, 5], 1, false)]
    fn test_molecule_ast_meet_pushout_stereo_atom(
        #[case] other_ligands: [u32; 4],
        #[case] other_coset: u32,
        #[case] admissible: bool,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let self_mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
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
        });
        let other_mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            stereo_atoms: vec![(
                AtomId(0),
                other_ligands
                    .into_iter()
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, other_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::new((0..6u32).map(|i| (NodeId(i), NodeId(i))).collect(), 6, 6)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let expected = admissible.then(|| self_mol.clone());
        assert_eq!(
            self_mol
                .meet_pushout(&other_mol, &overlap)
                .map(|po| po.object),
            expected,
        );
    }

    #[rstest]
    #[case::agree(
        [
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        ],
        1,
        true,
    )]
    #[case::contradict(
        [
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        ],
        0,
        false,
    )]
    #[case::ligand_set(
        [
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        ],
        1,
        false,
    )]
    fn test_molecule_ast_meet_pushout_stereo_bond(
        #[case] other_ligands: [StereoLigand; 4],
        #[case] other_coset: u32,
        #[case] admissible: bool,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
        ];
        let bonds = vec![(AtomId(0), AtomId(1), BondAst::from_order(2))];
        let self_mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondAst::new(StereoKind::CisTrans, 0u32),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let other_mol = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                other_ligands.to_vec(),
                StereoBondAst::new(StereoKind::CisTrans, other_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::new((0..4u32).map(|i| (NodeId(i), NodeId(i))).collect(), 4, 4)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let expected = admissible.then(|| self_mol.clone());
        assert_eq!(
            self_mol
                .meet_pushout(&other_mol, &overlap)
                .map(|po| po.object),
            expected,
        );
    }
}
