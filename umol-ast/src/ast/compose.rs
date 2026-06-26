//! Sequential reaction composition (A;B): the reaction whose application equals applying A then
//! B. Per overlap of A's product `R_A` with B's reactant `L_B`, the composite is built in one
//! frame and `canonicalize`d; overlaps with no `B.apply(A.apply(H))` witness (the DPO gluing
//! conditions) are rejected. See doc 132 §W3 for the frame algebra.

use std::collections::{HashMap, HashSet};

use umol_graph_core::{CommonSubgraphEnumerationAlgorithm, EdgeId, NodeId};

use super::atom::AtomAst;
use super::bond::BondAst;
use super::delta::{remap_delta, AtomDelta, BondDelta, Delta, Deltas};
use super::id::{AtomId, BondId};
use super::molecule::MoleculeAst;
use super::reaction::ReactionAst;
use super::traits::{Canonicalize, Lattice};

/// Which overlaps `compose` keeps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompositionScope {
    /// Keep only overlaps touching A's reaction center (sequential composites that actually
    /// chain A into B).
    RcAnchored,
    /// Keep every overlap, including the empty one (the free rule-algebra sum).
    Full,
}

impl ReactionAst {
    /// Sequential composites of `self` (A) then `other` (B): one per admissible overlap of A's
    /// product with B's reactant. `compose(A,B).apply(H)` equals `B.apply(A.apply(H))`.
    pub fn compose(&self, other: &ReactionAst, scope: CompositionScope) -> Vec<ReactionAst> {
        compose_all(self, other, scope).unwrap_or_default()
    }
}

fn created_atom_ids(deltas: &Deltas) -> Vec<AtomId> {
    let mut ids: Vec<AtomId> = deltas
        .iter()
        .filter_map(|d| match d {
            Delta::Atom(AtomDelta::Add { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    ids.sort();
    ids
}

fn created_bond_ids(deltas: &Deltas) -> Vec<BondId> {
    let mut ids: Vec<BondId> = deltas
        .iter()
        .filter_map(|d| match d {
            Delta::Bond(BondDelta::Add { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    ids.sort();
    ids
}

fn compose_all(
    a: &ReactionAst,
    b: &ReactionAst,
    scope: CompositionScope,
) -> Option<Vec<ReactionAst>> {
    let da = a.deltas.clone().canonicalize().ok()?;
    let db = b.deltas.clone().canonicalize().ok()?;
    let span_a = a.to_reaction_span().ok()?;
    let r_a = span_a.right();
    let l_b = &b.lhs;
    let n_a = a.lhs.atoms().count();
    let m_a = a.lhs.bonds().count();
    let n_b = l_b.atoms().count();
    let m_b = l_b.bonds().count();

    // R_A id ⇒ A-frame id: the span union frame is `lhs_A` in place then A-created appended, and
    // `right()` keeps the survivors in that order, so the k-th survivor's A-frame index is the
    // A-frame id of R_A atom/bond k.
    let ra_atom_aframe: Vec<usize> = span_a
        .atoms()
        .iter()
        .enumerate()
        .filter(|(_, change)| change.right().is_some())
        .map(|(aframe, _)| aframe)
        .collect();
    let mut ra_bond_aframe: Vec<usize> = Vec::new();
    for (aframe, change) in span_a.bonds().iter().enumerate() {
        if change.right().is_none() {
            continue;
        }
        let [x, y] = span_a.graph().edge_endpoints(EdgeId(aframe as u32));
        if span_a.atoms()[x.index()].right().is_some() && span_a.atoms()[y.index()].right().is_some()
        {
            ra_bond_aframe.push(aframe);
        }
    }
    let a_created_atoms = span_a.atoms().len() - n_a;
    let a_created_bonds = span_a.bonds().len() - m_a;

    // A-created delta id ⇒ A-frame rank (the span appends them sorted by id).
    let a_atom_rank: HashMap<AtomId, usize> = created_atom_ids(&da)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let a_bond_rank: HashMap<BondId, usize> = created_bond_ids(&da)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();

    // Reaction center, projected onto R_A: A-created elements (all) plus `lhs_A` elements A
    // modifies. A-frame id of a changed element → which R_A atoms/bonds it is.
    let mut rc_aframe_atoms: HashSet<usize> = HashSet::new();
    let mut rc_aframe_bonds: HashSet<usize> = HashSet::new();
    for delta in da.iter() {
        match delta {
            Delta::Atom(AtomDelta::SetField { id, .. })
            | Delta::Atom(AtomDelta::SetConstraint { id, .. }) => {
                rc_aframe_atoms.insert(id.index());
            }
            Delta::Atom(AtomDelta::Add { id, .. }) => {
                rc_aframe_atoms.insert(n_a + a_atom_rank[id]);
            }
            Delta::Bond(BondDelta::SetField { id, .. })
            | Delta::Bond(BondDelta::SetConstraint { id, .. }) => {
                rc_aframe_bonds.insert(id.index());
            }
            Delta::Bond(BondDelta::Add { id, .. }) => {
                rc_aframe_bonds.insert(m_a + a_bond_rank[id]);
            }
            _ => {}
        }
    }
    let rc_ra_atoms: HashSet<AtomId> = (0..ra_atom_aframe.len() as u32)
        .map(AtomId)
        .filter(|k| rc_aframe_atoms.contains(&ra_atom_aframe[k.index()]))
        .collect();
    let rc_ra_bonds: HashSet<BondId> = (0..ra_bond_aframe.len() as u32)
        .map(BondId)
        .filter(|k| rc_aframe_bonds.contains(&ra_bond_aframe[k.index()]))
        .collect();

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
    let overlaps = r_a.raw_graph().enumerate_common_subgraphs(
        l_b.raw_graph(),
        &mut node_match,
        &mut edge_match,
        CommonSubgraphEnumerationAlgorithm::BronKerbosch,
    );

    let db_created_atom_rank: HashMap<AtomId, usize> = created_atom_ids(&db)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let db_created_bond_rank: HashMap<BondId, usize> = created_bond_ids(&db)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let db_removed_atoms: HashSet<AtomId> = db
        .iter()
        .filter_map(|d| match d {
            Delta::Atom(AtomDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    let db_removed_bonds: HashSet<BondId> = db
        .iter()
        .filter_map(|d| match d {
            Delta::Bond(BondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();

    let mut results = Vec::new();
    for overlap in overlaps {
        let mapping = overlap.mapping();
        let lb_to_ra: HashMap<AtomId, AtomId> = mapping
            .iter()
            .map(|&(ra, lb)| (AtomId::from(lb), AtomId::from(ra)))
            .collect();
        let ra_to_lb: HashMap<AtomId, AtomId> = mapping
            .iter()
            .map(|&(ra, lb)| (AtomId::from(ra), AtomId::from(lb)))
            .collect();
        let overlap_lb: HashSet<AtomId> = lb_to_ra.keys().copied().collect();
        let overlap_ra: HashSet<AtomId> = ra_to_lb.keys().copied().collect();

        if matches!(scope, CompositionScope::RcAnchored) {
            let touches_atom = overlap_ra.iter().any(|ra| rc_ra_atoms.contains(ra));
            let touches_bond = rc_ra_bonds.iter().any(|&rb| {
                let [x, y] = r_a.raw_graph().edge_endpoints(EdgeId(rb.0));
                overlap_ra.contains(&AtomId::from(x)) && overlap_ra.contains(&AtomId::from(y))
            });
            if !touches_atom && !touches_bond {
                continue;
            }
        }

        let is_ra_created = |ra: AtomId| ra_atom_aframe[ra.index()] >= n_a;

        // Extra (non-overlap) L_B atoms → composite class 2 (`n_a..n_a+e`).
        let lb_extra: Vec<AtomId> = (0..n_b as u32)
            .map(AtomId)
            .filter(|x| !overlap_lb.contains(x))
            .collect();
        let e = lb_extra.len();
        let lb_extra_comp: HashMap<AtomId, usize> = lb_extra
            .iter()
            .enumerate()
            .map(|(rank, &x)| (x, n_a + rank))
            .collect();

        // Classify L_B bonds; context bonds (incident to an extra atom) become composite class 2.
        // Pushout-complement: a context bond whose overlap endpoint is A-created has no place in
        // `lhs_C` → the overlap is inadmissible.
        let mut context_bonds: Vec<BondId> = Vec::new();
        let mut admissible = true;
        for j in 0..m_b as u32 {
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(j));
            let (u, v) = (AtomId::from(u), AtomId::from(v));
            let u_overlap = overlap_lb.contains(&u);
            let v_overlap = overlap_lb.contains(&v);
            if u_overlap && v_overlap {
                continue;
            }
            if (u_overlap && is_ra_created(lb_to_ra[&u]))
                || (v_overlap && is_ra_created(lb_to_ra[&v]))
            {
                admissible = false;
                break;
            }
            context_bonds.push(BondId(j));
        }
        if !admissible {
            continue;
        }
        let f = context_bonds.len();
        let context_comp: HashMap<BondId, usize> = context_bonds
            .iter()
            .enumerate()
            .map(|(rank, &x)| (x, m_a + rank))
            .collect();

        // Combined-frame dangling: if B deletes a shared (overlap) atom, every R_A bond incident
        // to its image must be an overlap bond B also deletes; an A-product bond B cannot see
        // would dangle.
        let mut dangling = false;
        for &u in &db_removed_atoms {
            if !overlap_lb.contains(&u) {
                continue;
            }
            let ru = lb_to_ra[&u];
            for rb in r_a.atom(ru).bond_ids() {
                let [x, y] = r_a.raw_graph().edge_endpoints(EdgeId(rb.0));
                let other = if AtomId::from(x) == ru {
                    AtomId::from(y)
                } else {
                    AtomId::from(x)
                };
                if !overlap_ra.contains(&other) {
                    dangling = true;
                    break;
                }
                let w = ra_to_lb[&other];
                match l_b.raw_graph().find_edge(NodeId::from(u), NodeId::from(w)) {
                    Some(le) if db_removed_bonds.contains(&BondId::from(le)) => {}
                    _ => {
                        dangling = true;
                        break;
                    }
                }
            }
            if dangling {
                break;
            }
        }
        if dangling {
            continue;
        }

        // Composite id maps. A-frame id → composite: `lhs_A`/`lhs_A`-bonds keep their id; A-created
        // shift past the appended extras (atoms by `e`, bonds by `f`).
        let aframe_atom_comp = |aframe: usize| {
            if aframe < n_a {
                aframe
            } else {
                n_a + e + (aframe - n_a)
            }
        };
        let aframe_bond_comp = |aframe: usize| {
            if aframe < m_a {
                aframe
            } else {
                m_a + f + (aframe - m_a)
            }
        };
        let ra_atom_comp = |ra: AtomId| AtomId(aframe_atom_comp(ra_atom_aframe[ra.index()]) as u32);
        let ra_bond_comp = |rb: BondId| BondId(aframe_bond_comp(ra_bond_aframe[rb.index()]) as u32);

        let mut da_atom: HashMap<AtomId, AtomId> = (0..n_a as u32).map(|i| (AtomId(i), AtomId(i))).collect();
        for (&id, &rank) in &a_atom_rank {
            da_atom.insert(id, AtomId((n_a + e + rank) as u32));
        }
        let mut da_bond: HashMap<BondId, BondId> = (0..m_a as u32).map(|j| (BondId(j), BondId(j))).collect();
        for (&id, &rank) in &a_bond_rank {
            da_bond.insert(id, BondId((m_a + f + rank) as u32));
        }

        let b_atom_base = n_a + e + a_created_atoms;
        let b_bond_base = m_a + f + a_created_bonds;
        let mut db_atom: HashMap<AtomId, AtomId> = HashMap::new();
        for j in 0..n_b as u32 {
            let lb = AtomId(j);
            let comp = match lb_to_ra.get(&lb) {
                Some(&ra) => ra_atom_comp(ra),
                None => AtomId(lb_extra_comp[&lb] as u32),
            };
            db_atom.insert(lb, comp);
        }
        for (&id, &rank) in &db_created_atom_rank {
            db_atom.insert(id, AtomId((b_atom_base + rank) as u32));
        }
        let mut db_bond: HashMap<BondId, BondId> = HashMap::new();
        for j in 0..m_b as u32 {
            let lb = BondId(j);
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(j));
            let (u, v) = (AtomId::from(u), AtomId::from(v));
            let comp = if overlap_lb.contains(&u) && overlap_lb.contains(&v) {
                let re = r_a
                    .raw_graph()
                    .find_edge(NodeId::from(lb_to_ra[&u]), NodeId::from(lb_to_ra[&v]))
                    .expect("an induced overlap bond exists in R_A");
                ra_bond_comp(BondId::from(re))
            } else {
                BondId(context_comp[&lb] as u32)
            };
            db_bond.insert(lb, comp);
        }
        for (&id, &rank) in &db_created_bond_rank {
            db_bond.insert(id, BondId((b_bond_base + rank) as u32));
        }

        // lhs_C = lhs_A + extra L_B atoms, with lhs_A bonds + the L_B context bonds.
        let mut lc_atoms: Vec<AtomAst> = Vec::with_capacity(n_a + e);
        for i in 0..n_a as u32 {
            lc_atoms.push(a.lhs.atom(AtomId(i)).ast.clone());
        }
        for &x in &lb_extra {
            lc_atoms.push(l_b.atom(x).ast.clone());
        }
        let mut lc_bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::new();
        for j in 0..m_a as u32 {
            let [x, y] = a.lhs.raw_graph().edge_endpoints(EdgeId(j));
            lc_bonds.push((
                AtomId::from(x),
                AtomId::from(y),
                a.lhs.bond(BondId(j)).ast.clone(),
            ));
        }
        for &cb in &context_bonds {
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(cb.0));
            lc_bonds.push((
                db_atom[&AtomId::from(u)],
                db_atom[&AtomId::from(v)],
                l_b.bond(cb).ast.clone(),
            ));
        }
        let lhs_c = MoleculeAst::from_atoms_and_bonds(lc_atoms, lc_bonds);

        let mut deltas: Vec<Delta> = Vec::with_capacity(da.len() + db.len());
        for delta in da.iter() {
            deltas.push(remap_delta(delta.clone(), &da_atom, &da_bond));
        }
        for delta in db.iter() {
            deltas.push(remap_delta(delta.clone(), &db_atom, &db_bond));
        }
        let Ok(deltas) = Deltas::from_iter(deltas).canonicalize() else {
            continue;
        };
        results.push(ReactionAst::new(lhs_c, deltas));
    }
    Some(results)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::SubgraphIsomorphismAlgorithm;

    use super::super::edit::BondFieldChange;
    use super::super::value::ValueAst;
    use super::*;

    fn carbon_oxygen(order: u8) -> MoleculeAst {
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(order))],
        )
    }

    fn bond_order_rule(from: u8, to: u8) -> ReactionAst {
        ReactionAst::new(
            carbon_oxygen(from),
            Deltas::from_iter([Delta::Bond(BondDelta::SetField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(from as i64),
                    new: ValueAst::Lit(to as i64),
                },
            })]),
        )
    }

    #[rstest]
    fn test_reaction_ast_compose() {
        // A: C-O order 1→2, B: C-O order 2→3. The single overlap fuses to order 1→3.
        let a = bond_order_rule(1, 2);
        let b = bond_order_rule(2, 3);
        assert_eq!(
            a.compose(&b, CompositionScope::Full),
            vec![ReactionAst::new(
                carbon_oxygen(1),
                Deltas::from_iter([Delta::Bond(BondDelta::SetField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit(1),
                        new: ValueAst::Lit(3),
                    },
                })]),
            )],
        );
    }

    #[rstest]
    fn test_reaction_ast_compose_apply_equivalence() {
        // compose(A,B).apply(H) == B.apply(A.apply(H)).
        let a = bond_order_rule(1, 2);
        let b = bond_order_rule(2, 3);
        let host = carbon_oxygen(1);

        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|c| c.apply(&host, SubgraphIsomorphismAlgorithm::Vf2))
            .collect();
        let sequential: Vec<MoleculeAst> = a
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .flat_map(|intermediate| {
                b.apply(&intermediate, SubgraphIsomorphismAlgorithm::Vf2)
                    .collect::<Vec<_>>()
            })
            .collect();

        assert_eq!(composed, vec![carbon_oxygen(3)]);
        assert_eq!(sequential, vec![carbon_oxygen(3)]);
    }

    #[rstest]
    fn test_reaction_ast_compose_created_overlap() {
        // A: C → append an O bonded (order 1) — O is A-created. B: C-O order 1→2. The composite
        // appends the O already at order 2 (the create-then-modify fuses across the seam).
        let a = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
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
        );
        let b = bond_order_rule(1, 2);
        assert_eq!(
            a.compose(&b, CompositionScope::Full),
            vec![ReactionAst::new(
                MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
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
            )],
        );
    }

    #[rstest]
    fn test_reaction_ast_compose_inadmissible() {
        // A: N → append a C bonded (C is A-created), so R_A = N-C. B's reactant is N-C-O, so the
        // sole overlap maps the A-created C onto the middle atom, whose bond to the extra O is a
        // boundary bond on an A-created atom — no `B.apply(A.apply(H))` realizes it. Rejected.
        let a = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::N)], vec![]),
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
        );
        let b = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(1), AtomId(2), BondAst::from_order(1)),
                ],
            ),
            Deltas::new(),
        );
        assert_eq!(a.compose(&b, CompositionScope::Full), vec![]);
    }
}
