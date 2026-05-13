//! DPO (Double Pushout) graph rewriting on `MoleculeAst`.

use std::collections::{HashMap, HashSet};

use super::super::error::RewriteError;
use super::super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::reaction::{Assignment, ReactionRuleAst};
use super::MoleculeAst;

fn map_participants(
    participants: impl Iterator<Item = AtomId>,
    atom_map: &HashMap<AtomId, AtomId>,
) -> Option<HashSet<AtomId>> {
    participants.map(|a| atom_map.get(&a).copied()).collect()
}

impl MoleculeAst {
    /// Apply a DPO reaction rule to this molecule given a match assignment.
    ///
    /// The four phases (add R\K, modify K, remove L\K overlays, remove L\K
    /// topology) are executed on a `MoleculeBuilder` in the order that
    /// preserves index stability for K atoms during addition.
    pub fn apply_rule(
        &self,
        rule: &ReactionRuleAst,
        assignment: &Assignment,
    ) -> Result<MoleculeAst, RewriteError> {
        let lhs = &rule.lhs;
        let rhs = &rule.rhs;

        let l_to_g: HashMap<AtomId, AtomId> = assignment.atoms.iter().copied().collect();
        let l_to_r: HashMap<AtomId, AtomId> = rule.atom_map.iter().copied().collect();
        let r_to_l: HashMap<AtomId, AtomId> =
            rule.atom_map.iter().map(|&(l, r)| (r, l)).collect();

        let k_l: HashSet<AtomId> = l_to_r.keys().copied().collect();
        let l_atoms_g: HashSet<AtomId> = l_to_g.values().copied().collect();

        // Gluing condition: every neighbor of a deleted atom must be in L
        for l_atom in (0..lhs.atoms().count()).map(AtomId::from) {
            if k_l.contains(&l_atom) {
                continue;
            }
            let g_atom = *l_to_g
                .get(&l_atom)
                .ok_or(RewriteError::UnmappedLhsAtom(l_atom))?;
            for n in self.neighbors(g_atom) {
                if !l_atoms_g.contains(&n.atom) {
                    return Err(RewriteError::DanglingEdge {
                        atom: g_atom,
                        neighbor: n.atom,
                    });
                }
            }
        }

        // R-atom → G-atom map (K atoms filled now, R\K atoms filled in Phase 1)
        let mut r_to_g: HashMap<AtomId, AtomId> = HashMap::new();
        for &(l_atom, r_atom) in &rule.atom_map {
            let g_atom = *l_to_g
                .get(&l_atom)
                .ok_or(RewriteError::UnmappedLhsAtom(l_atom))?;
            r_to_g.insert(r_atom, g_atom);
        }

        let mut builder = self.edit();

        // Phase 1: add R\K atoms
        for r_atom in (0..rhs.atoms().count()).map(AtomId::from) {
            if r_to_l.contains_key(&r_atom) {
                continue;
            }
            let new_g = builder.add_atom(rhs[r_atom].clone());
            r_to_g.insert(r_atom, new_g);
        }

        // Phase 1: add R\K bonds
        for bv in rhs.bonds().iter() {
            let in_k = r_to_l.contains_key(&bv.atom_ids()[0])
                && r_to_l.contains_key(&bv.atom_ids()[1])
                && lhs
                    .bonds().connecting_id(r_to_l[&bv.atom_ids()[0]], r_to_l[&bv.atom_ids()[1]])
                    .is_some();
            if !in_k {
                builder.add_bond(r_to_g[&bv.atom_ids()[0]], r_to_g[&bv.atom_ids()[1]], bv.ast.clone());
            }
        }

        // Phase 1: add R\K dative bonds
        for dv in rhs.dative_bonds().iter() {
            let r_parts: Vec<AtomId> = dv.atom_ids().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomId> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.dative_bonds().connecting_id(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_donors: Vec<AtomId> = dv.donor_ids().map(|d| r_to_g[&d]).collect();
                let g_acceptor = r_to_g[&dv.acceptor_id];
                builder.add_dative_bond(g_donors, g_acceptor, dv.ast.clone());
            }
        }

        // Phase 1: add R\K aromatic systems
        for av in rhs.aromatic_systems().iter() {
            let r_parts: Vec<AtomId> = av.atom_ids().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomId> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.aromatic_systems().connecting_id(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_parts: Vec<AtomId> = r_parts.iter().map(|a| r_to_g[a]).collect();
                builder.add_aromatic_system(g_parts, av.ast.clone());
            }
        }

        // Phase 1: add R\K multicenter bonds
        for mv in rhs.multicenter_bonds().iter() {
            let r_parts: Vec<AtomId> = mv.atom_ids().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomId> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.multicenter_bonds().connecting_id(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_parts: Vec<AtomId> = r_parts.iter().map(|a| r_to_g[a]).collect();
                builder.add_multicenter_bond(g_parts, mv.ast.clone());
            }
        }

        // Phase 1: add R\K noncovalent bonds
        for nv in rhs.noncovalent_bonds().iter() {
            let [a, b] = nv.atom_ids();
            let in_k = r_to_l.contains_key(&a)
                && r_to_l.contains_key(&b)
                && lhs
                    .noncovalent_bonds().connecting_id(r_to_l[&a], r_to_l[&b])
                    .is_some();
            if !in_k {
                builder.add_noncovalent_bond(
                    [r_to_g[&a], r_to_g[&b]],
                    nv.ast.clone(),
                );
            }
        }

        // Phase 2: update K atom attributes from R
        for &(_, r_atom) in &rule.atom_map {
            let g_atom = r_to_g[&r_atom];
            *builder.atom_mut(g_atom).ast = rhs[r_atom].clone();
        }

        // Phase 2: update K bond attributes from R
        for bv in rhs.bonds().iter() {
            if !r_to_l.contains_key(&bv.atom_ids()[0]) || !r_to_l.contains_key(&bv.atom_ids()[1]) {
                continue;
            }
            let l_src = r_to_l[&bv.atom_ids()[0]];
            let l_tgt = r_to_l[&bv.atom_ids()[1]];
            if lhs.bonds().connecting_id(l_src, l_tgt).is_none() {
                continue;
            }
            if let Some(g_bond) = self.bonds().connecting_id(r_to_g[&bv.atom_ids()[0]], r_to_g[&bv.atom_ids()[1]]) {
                *builder.bond_mut(g_bond).ast = bv.ast.clone();
            }
        }

        // Phase 3: remove L\K dative bonds
        let mut remove_dative: Vec<DativeBondId> = Vec::new();
        for dv in lhs.dative_bonds().iter() {
            let l_parts: Vec<AtomId> = dv.atom_ids().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts: Option<HashSet<AtomId>> =
                    l_parts.iter().map(|a| l_to_r.get(a).copied()).collect();
                r_parts.is_some_and(|rp| rhs.dative_bonds().connecting_id(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts: Option<HashSet<AtomId>> =
                    l_parts.iter().map(|a| l_to_g.get(a).copied()).collect();
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.dative_bonds().connecting_id(gp.iter().copied()) {
                        remove_dative.push(g_idx);
                    }
                }
            }
        }
        if !remove_dative.is_empty() {
            builder.remove_dative_bonds(&remove_dative);
        }

        // Phase 3: remove L\K aromatic systems
        let mut remove_aromatic: Vec<AromaticSystemId> = Vec::new();
        for av in lhs.aromatic_systems().iter() {
            let l_parts: Vec<AtomId> = av.atom_ids().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts = map_participants(l_parts.iter().copied(), &l_to_r);
                r_parts.is_some_and(|rp| rhs.aromatic_systems().connecting_id(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts = map_participants(l_parts.iter().copied(), &l_to_g);
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.aromatic_systems().connecting_id(gp.iter().copied()) {
                        remove_aromatic.push(g_idx);
                    }
                }
            }
        }
        if !remove_aromatic.is_empty() {
            builder.remove_aromatic_systems(&remove_aromatic);
        }

        // Phase 3: remove L\K multicenter bonds
        let mut remove_multicenter: Vec<MulticenterBondId> = Vec::new();
        for mv in lhs.multicenter_bonds().iter() {
            let l_parts: Vec<AtomId> = mv.atom_ids().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts = map_participants(l_parts.iter().copied(), &l_to_r);
                r_parts.is_some_and(|rp| rhs.multicenter_bonds().connecting_id(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts = map_participants(l_parts.iter().copied(), &l_to_g);
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.multicenter_bonds().connecting_id(gp.iter().copied()) {
                        remove_multicenter.push(g_idx);
                    }
                }
            }
        }
        if !remove_multicenter.is_empty() {
            builder.remove_multicenter_bonds(&remove_multicenter);
        }

        // Phase 3: remove L\K noncovalent bonds
        let mut remove_noncovalent: Vec<NoncovalentBondId> = Vec::new();
        for nv in lhs.noncovalent_bonds().iter() {
            let [a, b] = nv.atom_ids();
            let in_k = k_l.contains(&a)
                && k_l.contains(&b)
                && rhs
                    .noncovalent_bonds().connecting_id(l_to_r[&a], l_to_r[&b])
                    .is_some();
            if !in_k {
                if let (Some(&ga), Some(&gb)) = (l_to_g.get(&a), l_to_g.get(&b))
                {
                    if let Some(g_idx) = self.noncovalent_bonds().connecting_id(ga, gb) {
                        remove_noncovalent.push(g_idx);
                    }
                }
            }
        }
        if !remove_noncovalent.is_empty() {
            builder.remove_noncovalent_bonds(&remove_noncovalent);
        }

        // Phase 4: remove L\K atoms and bonds from the topological graph
        let remove_atoms_g: Vec<AtomId> = (0..lhs.atoms().count())
            .map(AtomId::from)
            .filter(|a| !k_l.contains(a))
            .filter_map(|a| l_to_g.get(&a).copied())
            .collect();

        let mut remove_bonds_g: Vec<BondId> = Vec::new();
        for bv in lhs.bonds().iter() {
            let in_k = k_l.contains(&bv.atom_ids()[0])
                && k_l.contains(&bv.atom_ids()[1])
                && rhs
                    .bonds().connecting_id(l_to_r[&bv.atom_ids()[0]], l_to_r[&bv.atom_ids()[1]])
                    .is_some();
            if !in_k {
                if let (Some(&gs), Some(&gt_)) = (l_to_g.get(&bv.atom_ids()[0]), l_to_g.get(&bv.atom_ids()[1])) {
                    if let Some(g_bond) = self.bonds().connecting_id(gs, gt_) {
                        remove_bonds_g.push(g_bond);
                    }
                }
            }
        }

        builder.remove(&remove_atoms_g, &remove_bonds_g);
        Ok(builder.build())
    }
}
