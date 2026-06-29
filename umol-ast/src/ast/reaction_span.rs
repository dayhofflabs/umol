//! Reaction span AST: the superimposed `L ∪_K R` graph encoding a reaction's DPO rule span.
//!
//! Materialized superimposed graph carrying, per atom/bond, both its before and after state plus a
//! membership tag. The DPO span `L ←K─ R` is read off the tags — `K = Unchanged ∪ Modified`,
//! `L = K ∪ Removed`, `R = K ∪ Added` — and `right()` / `left()` project the two sides back to
//! a `MoleculeAst`. `Modified` (a preserved entity relabeled across the reaction) is the
//! relabeling-DPO reading: the entity persists in `K`, its label resolved per side.

// TODO: Add overlays. Molecule-level constraints and overlays not represented here yet,
// dropped on conversion from ReactionAst.

use std::collections::{BTreeMap, HashMap, HashSet};

use umol_graph_core::{EdgeId, Graph};

use super::atom::AtomAst;
use super::bond::BondAst;
use super::delta::{
    apply_atom_change, apply_bond_change, remap_delta, AtomDelta, BondDelta, Delta, Deltas,
    EntityFold, EntitySpan,
};
use super::error::Contradiction;
use super::id::{AtomId, BondId};
use super::molecule::MoleculeAst;
use super::reaction::ReactionAst;
use super::traits::Canonicalize;

/// The superimposed reaction graph — the reaction's DPO rule span, materialized. The union
/// topology is the `lhs` frame (deleted entities kept as nodes/edges) with created entities
/// appended; `atoms` / `bonds` are indexed parallel to the graph's nodes / edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst {
    graph: Graph,
    atoms: Vec<EntitySpan<AtomAst>>,
    bonds: Vec<EntitySpan<BondAst>>,
}

impl ReactionSpanAst {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn atoms(&self) -> &[EntitySpan<AtomAst>] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[EntitySpan<BondAst>] {
        &self.bonds
    }

    /// The left-hand (reactant) molecule: every entity present on the left, in a compacted
    /// frame (created entities dropped).
    pub fn left(&self) -> MoleculeAst {
        self.project(|atom| atom.left(), |bond| bond.left())
    }

    /// The right-hand (product) molecule: every entity present on the right, in a compacted
    /// frame (deleted entities dropped).
    pub fn right(&self) -> MoleculeAst {
        self.project(|atom| atom.right(), |bond| bond.right())
    }

    /// Recover the operational `ReactionAst` from the span — the inverse of
    /// `ReactionAst::to_reaction_span`, up to delta normal form. `lhs = left()` (which preserves
    /// the original lhs frame); each entity's `EntitySpan` yields its delta, a `Modified` one
    /// via an AST-diff of its left/right values.
    pub fn to_reaction(&self) -> ReactionAst {
        let mut deltas = AtomDelta::deltas_from_states(&self.atoms, |_| ());
        deltas.extend(BondDelta::deltas_from_states(&self.bonds, |edge| {
            let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
            [AtomId::from(a), AtomId::from(b)]
        }));
        ReactionAst::new(self.left(), Deltas::from_iter(deltas))
    }

    /// Project one side to a `MoleculeAst`. `atom_side` / `bond_side` pick the left or right
    /// value of each entity; absent entities are dropped and the survivors are renumbered.
    fn project(
        &self,
        atom_side: impl Fn(&EntitySpan<AtomAst>) -> Option<&AtomAst>,
        bond_side: impl Fn(&EntitySpan<BondAst>) -> Option<&BondAst>,
    ) -> MoleculeAst {
        let mut compacted: Vec<Option<AtomId>> = vec![None; self.atoms.len()];
        let mut atoms: Vec<AtomAst> = Vec::new();
        for (node, change) in self.atoms.iter().enumerate() {
            if let Some(ast) = atom_side(change) {
                compacted[node] = Some(AtomId(atoms.len() as u32));
                atoms.push(ast.clone());
            }
        }
        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::new();
        for (edge, change) in self.bonds.iter().enumerate() {
            if let Some(ast) = bond_side(change) {
                let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
                if let (Some(a), Some(b)) = (compacted[a.index()], compacted[b.index()]) {
                    bonds.push((a, b, ast.clone()));
                }
            }
        }
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }
}

impl ReactionAst {
    /// Materialize the superimposed reaction span. Canonicalizes the deltas, then annotates
    /// each `lhs` entity (in its own frame) with its before/after state — `Removed` /
    /// `Added` / `Modified` / `Unchanged` — appending created entities. A `Modified` entity's
    /// right value is its left value with the entity's field/constraint changes applied.
    /// `Err(Contradiction)` if the deltas are inconsistent (or inconsistent with `lhs`).
    pub fn to_reaction_span(&self) -> Result<ReactionSpanAst, Contradiction> {
        let deltas = self.deltas.clone().canonicalize()?;
        let lhs = &self.lhs;
        let atom_count = lhs.atoms().count();
        let bond_count = lhs.bonds().count();

        let mut removed_atoms: HashMap<AtomId, AtomAst> = HashMap::new();
        let mut added_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut atom_changes: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut removed_bonds: HashMap<BondId, BondAst> = HashMap::new();
        let mut added_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut bond_changes: HashMap<BondId, Vec<BondDelta>> = HashMap::new();

        for delta in deltas.iter() {
            match delta {
                Delta::Atom(atom) => match atom {
                    AtomDelta::Remove { id, ast } => {
                        removed_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::Add { id, ast } => {
                        added_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::ModifyField { id, .. } | AtomDelta::ModifyConstraint { id, .. } => {
                        atom_changes.entry(*id).or_default().push(atom.clone());
                    }
                },
                Delta::Bond(bond) => match bond {
                    BondDelta::Remove { id, ast, .. } => {
                        removed_bonds.insert(*id, ast.clone());
                    }
                    BondDelta::Add { id, atoms, ast } => {
                        added_bonds.insert(*id, (*atoms, ast.clone()));
                    }
                    BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. } => {
                        bond_changes.entry(*id).or_default().push(bond.clone());
                    }
                },
                // Molecule-level constraints are carried by the operational form, not the span.
                Delta::Constraint(_) => {}
            }
        }

        // Union node index per frame atom id: lhs atoms keep their id, created atoms append.
        let mut atom_index: HashMap<AtomId, usize> =
            HashMap::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            atom_index.insert(AtomId(node as u32), node);
        }
        for (offset, &id) in added_atoms.keys().enumerate() {
            atom_index.insert(id, atom_count + offset);
        }

        let mut atoms: Vec<EntitySpan<AtomAst>> =
            Vec::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            let id = AtomId(node as u32);
            if let Some(ast) = removed_atoms.get(&id) {
                atoms.push(EntitySpan::Removed(ast.clone()));
            } else if let Some(changes) = atom_changes.get(&id) {
                let left = lhs.atom(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_atom_change(&mut right, change)?;
                }
                atoms.push(EntitySpan::Modified { left, right });
            } else {
                atoms.push(EntitySpan::Unchanged(lhs.atom(id).ast.clone()));
            }
        }
        for ast in added_atoms.into_values() {
            atoms.push(EntitySpan::Added(ast));
        }

        let mut bonds: Vec<EntitySpan<BondAst>> =
            Vec::with_capacity(bond_count + added_bonds.len());
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count + added_bonds.len());
        for edge in 0..bond_count {
            let id = BondId(edge as u32);
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(edge as u32));
            edges.push([a.0, b.0]);
            if let Some(ast) = removed_bonds.get(&id) {
                bonds.push(EntitySpan::Removed(ast.clone()));
            } else if let Some(changes) = bond_changes.get(&id) {
                let left = lhs.bond(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_bond_change(&mut right, change)?;
                }
                bonds.push(EntitySpan::Modified { left, right });
            } else {
                bonds.push(EntitySpan::Unchanged(lhs.bond(id).ast.clone()));
            }
        }
        for (atoms, ast) in added_bonds.into_values() {
            edges.push([atom_index[&atoms[0]] as u32, atom_index[&atoms[1]] as u32]);
            bonds.push(EntitySpan::Added(ast));
        }

        let graph = Graph::new(atoms.len(), &edges);
        Ok(ReactionSpanAst {
            graph,
            atoms,
            bonds,
        })
    }

    /// The reverse reaction: the product becomes the reactant and every delta is inverted and
    /// re-anchored to the product's (compacted) frame. `reverse().to_reaction_span()` swaps the
    /// sides of `self`'s span. `Err(Contradiction)` if the deltas are inconsistent.
    pub fn reverse(&self) -> Result<ReactionAst, Contradiction> {
        let deltas = self.deltas.clone().canonicalize()?;
        let new_lhs = self.to_reaction_span()?.right();
        let atom_count = self.lhs.atoms().count();
        let bond_count = self.lhs.bonds().count();

        let mut removed_atoms: Vec<AtomId> = Vec::new();
        let mut created_atoms: Vec<AtomId> = Vec::new();
        let mut removed_bonds: Vec<BondId> = Vec::new();
        let mut created_bonds: Vec<BondId> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::Atom(AtomDelta::Remove { id, .. }) => removed_atoms.push(*id),
                Delta::Atom(AtomDelta::Add { id, .. }) => created_atoms.push(*id),
                Delta::Bond(BondDelta::Remove { id, .. }) => removed_bonds.push(*id),
                Delta::Bond(BondDelta::Add { id, .. }) => created_bonds.push(*id),
                _ => {}
            }
        }
        created_atoms.sort();
        created_bonds.sort();

        // Forward → reverse-frame maps, matching `right()`'s compaction: survivors take ids in
        // union order (lhs in place, created appended); deleted entities become created in the
        // reverse and take fresh ids after the survivors.
        let removed_atom_set: HashSet<AtomId> = removed_atoms.iter().copied().collect();
        let removed_bond_set: HashSet<BondId> = removed_bonds.iter().copied().collect();
        let survivor_atoms = (0..atom_count as u32)
            .map(AtomId)
            .filter(|id| !removed_atom_set.contains(id))
            .chain(created_atoms.iter().copied());
        let rev_atom: HashMap<AtomId, AtomId> = survivor_atoms
            .chain(removed_atoms.iter().copied())
            .enumerate()
            .map(|(rev, id)| (id, AtomId(rev as u32)))
            .collect();
        let survivor_bonds = (0..bond_count as u32)
            .map(BondId)
            .filter(|id| !removed_bond_set.contains(id))
            .chain(created_bonds.iter().copied());
        let rev_bond: HashMap<BondId, BondId> = survivor_bonds
            .chain(removed_bonds.iter().copied())
            .enumerate()
            .map(|(rev, id)| (id, BondId(rev as u32)))
            .collect();

        let reversed: Deltas = deltas
            .iter()
            .map(|delta| remap_delta(delta.clone().inverse(), &rev_atom, &rev_bond))
            .collect();
        Ok(ReactionAst::new(new_lhs, reversed))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::delta::Deltas;
    use super::super::edit::BondFieldChange;
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    fn test_reaction_ast_to_reaction_span() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
        let span = reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [EntitySpan::Modified {
                left: BondAst::from_order(1),
                right: BondAst::from_order(2),
            }],
        );
    }

    #[rstest]
    fn test_reaction_span_ast_right(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Removed(AtomAst::from_element(Element::O)),
                EntitySpan::Added(AtomAst::from_element(Element::N)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [
                EntitySpan::Removed(BondAst::from_order(1)),
                EntitySpan::Added(BondAst::from_order(1)),
            ],
        );
        assert_eq!(
            span.right(),
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_left(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.left(),
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
        );
    }

    #[rstest]
    #[case::order_change(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
    )]
    #[case::substitution(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
    )]
    fn test_reaction_ast_reverse(
        #[case] forward: ReactionAst,
        #[case] expected_reactant: MoleculeAst,
        #[case] expected_product: MoleculeAst,
    ) {
        // The reverse reaction's reactant is the forward product; its product is the forward
        // reactant.
        let span = forward.reverse().unwrap().to_reaction_span().unwrap();
        assert_eq!(span.left(), expected_reactant);
        assert_eq!(span.right(), expected_product);
    }

    // C-O with atom 1 (O) and its bond removed, replaced by a new N (atom 2) bonded to C.
    #[fixture]
    fn substitution_reaction() -> ReactionAst {
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        )
    }
}
