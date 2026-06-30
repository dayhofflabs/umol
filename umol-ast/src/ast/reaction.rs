//! Reaction AST: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations live in `reaction_span.rs`).

use std::collections::{BTreeMap, HashMap, HashSet};

use umol_graph_core::SubgraphIsomorphismAlgorithm;

use super::atom::AtomAst;
use super::bond::BondAst;
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, DativeBondDelta, Delta, Deltas, MulticenterBondDelta,
    NoncovalentBondDelta,
};
use super::edit::{
    AddBond, AromaticSystemRef, AtomRef, BondRef, DativeBondRef, Edit, MulticenterBondRef,
    NoncovalentBondRef,
};
use super::embedding::MoleculeEmbedding;
use super::error::{ApplyError, Contradiction};
use super::id::{AtomId, BondId};
use super::molecule::MoleculeAst;
use super::substructure::SubstructureMatchAlgorithm;
use super::traits::Canonicalize;

/// A reaction as one full molecule state (`lhs`) plus one resolved delta (`deltas`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionAst {
    pub lhs: MoleculeAst,
    pub deltas: Deltas,
}

impl ReactionAst {
    pub fn new(lhs: MoleculeAst, deltas: Deltas) -> Self {
        Self { lhs, deltas }
    }

    /// Apply the reaction at one match `m` of `lhs` into a host (`m.ast()`), producing the
    /// transformed host. DPO: a deleted host atom must carry no localized bond the rule does not
    /// also delete (else `ApplyError::Dangling`). Created atoms/bonds are appended, preserved
    /// entities are mutated in place, deleted entities are removed (the host renumbers).
    /// Molecule-level constraints are not applied (deferred with the span's overlay scope).
    pub fn apply_at(&self, m: &MoleculeEmbedding) -> Result<MoleculeAst, ApplyError> {
        let deltas = self.deltas.clone().canonicalize()?;
        let host = m.ast();

        let mut created_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut created_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut sets: Vec<Edit> = Vec::new();
        let mut remove_atoms: Vec<AtomRef> = Vec::new();
        let mut remove_bonds: Vec<BondRef> = Vec::new();
        let mut removed_host_atoms: Vec<AtomId> = Vec::new();
        let mut removed_host_bonds: HashSet<BondId> = HashSet::new();

        for delta in deltas.iter() {
            match delta {
                Delta::Atom(AtomDelta::Add { id, ast }) => {
                    created_atoms.insert(*id, ast.clone());
                }
                Delta::Atom(AtomDelta::Remove { id, .. }) => {
                    let host_atom = m.host_atom(*id);
                    removed_host_atoms.push(host_atom);
                    remove_atoms.push(AtomRef::Id(host_atom));
                }
                Delta::Atom(AtomDelta::ModifyField { id, change }) => {
                    sets.push(Edit::ModifyAtomField {
                        id: AtomRef::Id(m.host_atom(*id)),
                        change: change.clone(),
                    })
                }
                Delta::Atom(AtomDelta::ModifyConstraint { id, old, new }) => {
                    sets.push(Edit::ModifyAtomConstraint {
                        id: AtomRef::Id(m.host_atom(*id)),
                        old: old.clone(),
                        new: new.clone(),
                    })
                }
                Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                    created_bonds.insert(*id, (*atoms, ast.clone()));
                }
                Delta::Bond(BondDelta::Remove { id, .. }) => {
                    let host_bond = m.host_bond(*id);
                    removed_host_bonds.insert(host_bond);
                    remove_bonds.push(BondRef::Id(host_bond));
                }
                Delta::Bond(BondDelta::ModifyField { id, change }) => {
                    sets.push(Edit::ModifyBondField {
                        id: BondRef::Id(m.host_bond(*id)),
                        change: change.clone(),
                    })
                }
                Delta::Bond(BondDelta::ModifyConstraint { id, old, new }) => {
                    sets.push(Edit::ModifyBondConstraint {
                        id: BondRef::Id(m.host_bond(*id)),
                        old: old.clone(),
                        new: new.clone(),
                    })
                }
                Delta::DativeBond(d) => match d {
                    DativeBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyDativeBondField {
                            id: DativeBondRef::Id(m.host_dative_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    DativeBondDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyDativeBondConstraint {
                            id: DativeBondRef::Id(m.host_dative_bond(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    DativeBondDelta::Add { .. } | DativeBondDelta::Remove { .. } => {}
                },
                Delta::AromaticSystem(a) => match a {
                    AromaticSystemDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyAromaticSystemField {
                            id: AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                            change: change.clone(),
                        })
                    }
                    AromaticSystemDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyAromaticSystemConstraint {
                            id: AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    AromaticSystemDelta::Add { .. } | AromaticSystemDelta::Remove { .. } => {}
                },
                Delta::MulticenterBond(mc) => match mc {
                    MulticenterBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyMulticenterBondField {
                            id: MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    MulticenterBondDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyMulticenterBondConstraint {
                            id: MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    MulticenterBondDelta::Add { .. } | MulticenterBondDelta::Remove { .. } => {}
                },
                Delta::NoncovalentBond(nc) => match nc {
                    NoncovalentBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyNoncovalentBondField {
                            id: NoncovalentBondRef::Id(m.host_noncovalent_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    // `NoncovalentBondConstraint` is uninhabited — no `Edit` variant, no-op.
                    NoncovalentBondDelta::ModifyConstraint { .. } => {}
                    NoncovalentBondDelta::Add { .. } | NoncovalentBondDelta::Remove { .. } => {}
                },
                // Molecule-level constraints are deferred (see `to_reaction_span`).
                Delta::Constraint(_) => {}
            }
        }

        // DPO gluing condition: a deleted host atom keeps no bond the rule does not delete.
        for &host_atom in &removed_host_atoms {
            for bond in host.atom(host_atom).bond_ids() {
                if !removed_host_bonds.contains(&bond) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
        }

        // `AddAtoms` is the first edit, so created atoms take `New(0..k)` in ascending id order.
        let new_index: HashMap<AtomId, usize> = created_atoms
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();
        let atom_ref = |id: AtomId| match new_index.get(&id) {
            Some(&index) => AtomRef::New(index),
            None => AtomRef::Id(m.host_atom(id)),
        };

        // Overlay create/remove need `atom_ref` (created participants resolve to `New`), so they
        // are lowered in a second pass: adds after the topology adds, removes before
        // `RemoveTopology`. Dative `atoms` is `[donors…, acceptor]` (acceptor last, per transact).
        let mut overlay_adds: Vec<Edit> = Vec::new();
        let mut overlay_removes: Vec<Edit> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::DativeBond(DativeBondDelta::Add {
                    donors, acceptor, ast, ..
                }) => {
                    let mut atoms: Vec<AtomRef> = donors.iter().map(|a| atom_ref(*a)).collect();
                    atoms.push(atom_ref(*acceptor));
                    overlay_adds.push(Edit::AddDativeBond {
                        atoms,
                        ast: ast.clone(),
                    });
                }
                Delta::DativeBond(DativeBondDelta::Remove {
                    id, donors, acceptor, ast,
                }) => {
                    let mut atoms: Vec<AtomRef> = donors.iter().map(|a| atom_ref(*a)).collect();
                    atoms.push(atom_ref(*acceptor));
                    overlay_removes.push(Edit::RemoveDativeBond {
                        id: DativeBondRef::Id(m.host_dative_bond(*id)),
                        atoms,
                        ast: ast.clone(),
                    });
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddAromaticSystem {
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, ast }) => {
                    overlay_removes.push(Edit::RemoveAromaticSystem {
                        id: AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddMulticenterBond {
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, ast }) => {
                    overlay_removes.push(Edit::RemoveMulticenterBond {
                        id: MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddNoncovalentBond {
                        atoms: [atom_ref(atoms[0]), atom_ref(atoms[1])],
                        ast: ast.clone(),
                    });
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, ast }) => {
                    overlay_removes.push(Edit::RemoveNoncovalentBond {
                        id: NoncovalentBondRef::Id(m.host_noncovalent_bond(*id)),
                        atoms: [atom_ref(atoms[0]), atom_ref(atoms[1])],
                        ast: ast.clone(),
                    });
                }
                _ => {}
            }
        }

        let mut edits: Vec<Edit> = Vec::new();
        if !created_atoms.is_empty() {
            edits.push(Edit::AddAtoms {
                atoms: created_atoms.values().cloned().collect(),
            });
        }
        if !created_bonds.is_empty() {
            edits.push(Edit::AddBonds {
                bonds: created_bonds
                    .values()
                    .map(|(atoms, ast)| AddBond {
                        endpoints: [atom_ref(atoms[0]), atom_ref(atoms[1])],
                        ast: ast.clone(),
                    })
                    .collect(),
            });
        }
        edits.extend(overlay_adds);
        edits.extend(sets);
        edits.extend(overlay_removes);
        if !remove_atoms.is_empty() || !remove_bonds.is_empty() {
            edits.push(Edit::RemoveTopology {
                atoms: remove_atoms,
                bonds: remove_bonds,
            });
        }

        let mut builder = host.edit();
        builder.transact(edits)?;
        Ok(builder.build())
    }

    /// Every product of applying the reaction to `host`: one per injective match of `lhs` into
    /// `host` (via `subiso`) that satisfies the DPO gluing condition. Matches that dangle are
    /// skipped.
    pub fn apply<'h>(
        &'h self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> impl Iterator<Item = MoleculeAst> + 'h {
        self.lhs
            .substructure_matches(host, SubstructureMatchAlgorithm::GraphAndOverlays, subiso)
            .into_iter()
            .filter_map(move |m| self.apply_at(&m).ok())
    }
}

impl Canonicalize for ReactionAst {
    /// Value-level in a fixed atom id space: `deltas` are canonicalized;
    /// `lhs` is passed through (`MoleculeAst` has no whole-molecule canonical form — its
    /// equality is structural). Equality up to atom renumbering is a separate `umol-graph`
    /// operation.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            lhs: self.lhs,
            deltas: self.deltas.canonicalize()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::edit::{AtomFieldChange, BondFieldChange};
    use super::super::value::ValueAst;
    use super::*;

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(id),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(old),
                new: ValueAst::Lit(new),
            },
        })
    }

    #[rstest]
    fn test_reaction_ast_canonicalize() {
        // The delta chain fuses; the lhs is passed through unchanged.
        let reaction = ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)]),
        );
        assert_eq!(
            reaction.canonicalize().unwrap(),
            ReactionAst::new(
                MoleculeAst::default(),
                Deltas::from_iter([charge_set(0, 0, 2)])
            ),
        );
    }

    #[rstest]
    fn test_reaction_ast_apply_at() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
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
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let embedding = MoleculeEmbedding::from_match(
            &host,
            &reaction.lhs,
            vec![AtomId(0), AtomId(1)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(
            reaction.apply_at(&embedding).unwrap(),
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
        );
    }

    #[rstest]
    fn test_reaction_ast_apply_at_error() {
        // The rule deletes a lone atom; its host image still carries an undeleted bond → dangling.
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        );
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let embedding = MoleculeEmbedding::from_match(
            &host,
            &reaction.lhs,
            vec![AtomId(0)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(
            reaction.apply_at(&embedding),
            Err(ApplyError::Dangling {
                host_atom: AtomId(0),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_apply() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
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
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .collect();
        assert_eq!(
            products,
            vec![MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            )],
        );
    }
}
