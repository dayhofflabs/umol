//! The attributed pushout of two molecules: glue `self` and `other` on a shared subgraph, meeting
//! atom / bond data (and, later, overlays and constraints) where they coincide — the umol-ast
//! attribute layer over graph-core's structural `pushout`. A child of `molecule` so it reaches the
//! private overlay relation-sets directly, without exposing raw accessors.

use umol_graph_core::{
    EdgeId, FactorOrdering, FixedRelationSet, FixedVarBirelationSet, GraphCorrespondence, NodeId,
    ParticipantPosition, Unordered, VarRelationSet,
};

use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::Constraints;
use super::super::correspondence::MoleculeCorrespondence;
use super::super::dative::DativeBondAst;
use super::super::id::{AtomId, BondId};
use super::super::noncovalent::NoncovalentBondAst;
use super::super::traits::Lattice;
use super::MoleculeAst;

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
    /// coincident bonds), meeting atom and bond data at coincident entities; `None` when any
    /// coincident `meet` is `⊥` (the overlap is inadmissible). Node/edge layer only — overlays and
    /// molecule constraints are not yet carried.
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
        let map_node = |n: NodeId| po.right.nodes().right_of(n).expect("right total on other");

        let aromatic = glue_var_overlays(
            &self.aromatic_systems,
            &other.aromatic_systems,
            map_node,
            |ast, sigma| ast.permute(sigma),
        )?;
        let multicenter = glue_var_overlays(
            &self.multicenter_bonds,
            &other.multicenter_bonds,
            map_node,
            |ast, sigma| ast.permute(sigma),
        )?;

        let dative_glue = FixedVarBirelationSet::new(
            other
                .dative_bonds
                .relation_ids()
                .map(|id| {
                    let acceptor = [map_node(other.dative_bonds.participants_1(id)[0])];
                    let donors: Vec<NodeId> = other
                        .dative_bonds
                        .participants_2(id)
                        .iter()
                        .map(|&n| map_node(n))
                        .collect();
                    (acceptor, donors, other.dative_bonds.data(id).clone())
                })
                .collect(),
        );
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

        let noncovalent_glue = FixedRelationSet::new(
            other
                .noncovalent_bonds
                .relation_ids()
                .map(|id| {
                    let &[u, v] = other.noncovalent_bonds.participants(id);
                    (
                        [map_node(u), map_node(v)],
                        other.noncovalent_bonds.data(id).clone(),
                    )
                })
                .collect(),
        );
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

        let object = MoleculeAst::from_parts(
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            Vec::new(),
            Vec::new(),
            Constraints::new(),
        );

        let left = MoleculeCorrespondence::induce(self, &object, po.left.nodes().clone());
        let right = MoleculeCorrespondence::induce(other, &object, po.right.nodes().clone());
        Some(MoleculePushout {
            object,
            left,
            right,
        })
    }
}

/// Glue two `VarRelationSet` overlay families over the pushout node relabeling `map_node` (`right` →
/// glue; `left` already lives in the glue id space). Each `right` datum's electron ordering is
/// re-indexed to the canonicalized participant order, then coinciding overlays merge by `meet` and
/// non-coinciding ones are appended (context). `None` if any coincident meet is `⊥`.
fn glue_var_overlays<D: Lattice + Clone>(
    left: &VarRelationSet<NodeId, Unordered, D>,
    right: &VarRelationSet<NodeId, Unordered, D>,
    map_node: impl Fn(NodeId) -> NodeId,
    mut permute: impl FnMut(&mut D, &[ParticipantPosition]),
) -> Option<Vec<(Vec<AtomId>, D)>> {
    let right_glue = VarRelationSet::new(
        right
            .relation_ids()
            .map(|id| {
                let mut members: Vec<NodeId> = right
                    .participants(id)
                    .iter()
                    .map(|&n| map_node(n))
                    .collect();
                let sigma = Unordered::canonicalize_positions(&mut members);
                let mut data = right.data(id).clone();
                permute(&mut data, &sigma);
                (members, data)
            })
            .collect(),
    );
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

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::Correspondence;

    use super::super::super::aromatic::AromaticSystemAst;
    use super::*;

    // A single shared atom (node 0 ↔ node 0), no shared bond; each side has one exposed atom.
    #[fixture]
    fn overlap() -> GraphCorrespondence {
        GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(0))], 2, 2),
            Correspondence::new(vec![], 1, 1),
        )
    }

    // meet_pushout glues over the shared atom: `object` keeps left's ids, appends right's exposed
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
        let left = MoleculeAst::from_atoms_and_bonds(
            vec![left_shared, AtomAst::from_element(Element::N)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let right = MoleculeAst::from_atoms_and_bonds(
            vec![right_shared, AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let expected = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(shared_element),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
            ],
        );
        assert_eq!(
            left.meet_pushout(&right, &overlap)
                .expect("admissible glue")
                .object,
            expected,
        );
    }

    // A conflicting shared atom (meet = ⊥) makes the whole glue inadmissible.
    #[rstest]
    #[case::carbon_nitrogen(Element::C, Element::N)]
    #[case::oxygen_nitrogen(Element::O, Element::N)]
    fn test_molecule_ast_meet_pushout_inadmissible(
        overlap: GraphCorrespondence,
        #[case] left_element: Element,
        #[case] right_element: Element,
    ) {
        let left = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(left_element),
                AtomAst::from_element(Element::N),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let right = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(right_element),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        assert!(left.meet_pushout(&right, &overlap).is_none());
    }

    // Overlays glue over the same pushout: a system shared on the same atoms coincides (met to one),
    // one on other atoms is context (kept) — both sides' overlays survive.
    #[rstest]
    fn test_molecule_ast_meet_pushout_overlays() {
        let full_overlap = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(0)),
                    (NodeId(1), NodeId(1)),
                    (NodeId(2), NodeId(2)),
                ],
                3,
                3,
            ),
            Correspondence::new(vec![(EdgeId(0), EdgeId(0)), (EdgeId(1), EdgeId(1))], 2, 2),
        );
        // left aromatic on {0,1}; right aromatic on {0,1} (coincides) and {1,2} (context).
        let left = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 3],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            vec![],
            vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let right = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 3],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            vec![],
            vec![
                (vec![AtomId(0), AtomId(1)], AromaticSystemAst::default()),
                (vec![AtomId(1), AtomId(2)], AromaticSystemAst::default()),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let expected = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 3],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            vec![],
            vec![
                (vec![AtomId(0), AtomId(1)], AromaticSystemAst::default()),
                (vec![AtomId(1), AtomId(2)], AromaticSystemAst::default()),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        assert_eq!(
            left.meet_pushout(&right, &full_overlap)
                .expect("admissible")
                .object,
            expected,
        );
    }
}
