//! DPO (Double Pushout) graph rewriting on `MoleculeAst`.

use std::collections::{HashMap, HashSet};

use super::super::error::RewriteError;
use super::super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::super::reaction::{Assignment, ReactionRuleAst};
use super::MoleculeAst;

fn map_participants(
    participants: impl Iterator<Item = AtomIdx>,
    atom_map: &HashMap<AtomIdx, AtomIdx>,
) -> Option<HashSet<AtomIdx>> {
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

        let l_to_g: HashMap<AtomIdx, AtomIdx> = assignment.atoms.iter().copied().collect();
        let l_to_r: HashMap<AtomIdx, AtomIdx> = rule.atom_map.iter().copied().collect();
        let r_to_l: HashMap<AtomIdx, AtomIdx> =
            rule.atom_map.iter().map(|&(l, r)| (r, l)).collect();

        let k_l: HashSet<AtomIdx> = l_to_r.keys().copied().collect();
        let l_atoms_g: HashSet<AtomIdx> = l_to_g.values().copied().collect();

        // Gluing condition: every neighbor of a deleted atom must be in L
        for l_atom in (0..lhs.atoms().count()).map(AtomIdx::from) {
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
        let mut r_to_g: HashMap<AtomIdx, AtomIdx> = HashMap::new();
        for &(l_atom, r_atom) in &rule.atom_map {
            let g_atom = *l_to_g
                .get(&l_atom)
                .ok_or(RewriteError::UnmappedLhsAtom(l_atom))?;
            r_to_g.insert(r_atom, g_atom);
        }

        let mut builder = self.edit();

        // Phase 1: add R\K atoms
        for r_atom in (0..rhs.atoms().count()).map(AtomIdx::from) {
            if r_to_l.contains_key(&r_atom) {
                continue;
            }
            let new_g = builder.add_atom(rhs[r_atom].clone());
            r_to_g.insert(r_atom, new_g);
        }

        // Phase 1: add R\K bonds
        for bv in rhs.bonds().iter() {
            let in_k = r_to_l.contains_key(&bv.src)
                && r_to_l.contains_key(&bv.tgt)
                && lhs
                    .connecting_bond(r_to_l[&bv.src], r_to_l[&bv.tgt])
                    .is_some();
            if !in_k {
                builder.add_bond(r_to_g[&bv.src], r_to_g[&bv.tgt], bv.data.clone());
            }
        }

        // Phase 1: add R\K dative bonds
        for dv in rhs.dative_bonds().iter() {
            let r_parts: Vec<AtomIdx> = dv.atoms().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomIdx> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.connecting_dative_bond(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_donors: Vec<AtomIdx> = dv.donors().map(|d| r_to_g[&d]).collect();
                let g_acceptor = r_to_g[&dv.acceptor];
                builder.add_dative_bond(g_donors, g_acceptor, dv.data.clone());
            }
        }

        // Phase 1: add R\K aromatic systems
        for av in rhs.aromatic_systems().iter() {
            let r_parts: Vec<AtomIdx> = av.atoms().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomIdx> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.connecting_aromatic_system(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_parts: Vec<AtomIdx> = r_parts.iter().map(|a| r_to_g[a]).collect();
                builder.add_aromatic_system(g_parts, av.data.clone());
            }
        }

        // Phase 1: add R\K multicenter bonds
        for mv in rhs.multicenter_bonds().iter() {
            let r_parts: Vec<AtomIdx> = mv.atoms().collect();
            let all_k = r_parts.iter().all(|a| r_to_l.contains_key(a));
            let in_k = all_k && {
                let l_parts: HashSet<AtomIdx> = r_parts.iter().map(|a| r_to_l[a]).collect();
                lhs.connecting_multicenter_bond(l_parts.iter().copied()).is_some()
            };
            if !in_k {
                let g_parts: Vec<AtomIdx> = r_parts.iter().map(|a| r_to_g[a]).collect();
                builder.add_multicenter_bond(g_parts, mv.data.clone());
            }
        }

        // Phase 1: add R\K noncovalent bonds
        for nv in rhs.noncovalent_bonds().iter() {
            let in_k = r_to_l.contains_key(&nv.atoms[0])
                && r_to_l.contains_key(&nv.atoms[1])
                && lhs
                    .connecting_noncovalent_bond(r_to_l[&nv.atoms[0]], r_to_l[&nv.atoms[1]])
                    .is_some();
            if !in_k {
                builder.add_noncovalent_bond(
                    [r_to_g[&nv.atoms[0]], r_to_g[&nv.atoms[1]]],
                    nv.data.clone(),
                );
            }
        }

        // Phase 2: update K atom attributes from R
        for &(_, r_atom) in &rule.atom_map {
            let g_atom = r_to_g[&r_atom];
            *builder.atom_mut(g_atom) = rhs[r_atom].clone();
        }

        // Phase 2: update K bond attributes from R
        for bv in rhs.bonds().iter() {
            if !r_to_l.contains_key(&bv.src) || !r_to_l.contains_key(&bv.tgt) {
                continue;
            }
            let l_src = r_to_l[&bv.src];
            let l_tgt = r_to_l[&bv.tgt];
            if lhs.connecting_bond(l_src, l_tgt).is_none() {
                continue;
            }
            if let Some(g_bond) = self.connecting_bond(r_to_g[&bv.src], r_to_g[&bv.tgt]) {
                *builder.bond_mut(g_bond) = bv.data.clone();
            }
        }

        // Phase 3: remove L\K dative bonds
        let mut remove_dative: Vec<DativeBondIdx> = Vec::new();
        for dv in lhs.dative_bonds().iter() {
            let l_parts: Vec<AtomIdx> = dv.atoms().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts: Option<HashSet<AtomIdx>> =
                    l_parts.iter().map(|a| l_to_r.get(a).copied()).collect();
                r_parts.is_some_and(|rp| rhs.connecting_dative_bond(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts: Option<HashSet<AtomIdx>> =
                    l_parts.iter().map(|a| l_to_g.get(a).copied()).collect();
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.connecting_dative_bond(gp.iter().copied()) {
                        remove_dative.push(g_idx);
                    }
                }
            }
        }
        if !remove_dative.is_empty() {
            builder.remove_dative_bonds(&remove_dative);
        }

        // Phase 3: remove L\K aromatic systems
        let mut remove_aromatic: Vec<AromaticSystemIdx> = Vec::new();
        for av in lhs.aromatic_systems().iter() {
            let l_parts: Vec<AtomIdx> = av.atoms().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts = map_participants(l_parts.iter().copied(), &l_to_r);
                r_parts.is_some_and(|rp| rhs.connecting_aromatic_system(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts = map_participants(l_parts.iter().copied(), &l_to_g);
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.connecting_aromatic_system(gp.iter().copied()) {
                        remove_aromatic.push(g_idx);
                    }
                }
            }
        }
        if !remove_aromatic.is_empty() {
            builder.remove_aromatic_systems(&remove_aromatic);
        }

        // Phase 3: remove L\K multicenter bonds
        let mut remove_multicenter: Vec<MulticenterBondIdx> = Vec::new();
        for mv in lhs.multicenter_bonds().iter() {
            let l_parts: Vec<AtomIdx> = mv.atoms().collect();
            let all_k = l_parts.iter().all(|a| k_l.contains(a));
            let in_k = all_k && {
                let r_parts = map_participants(l_parts.iter().copied(), &l_to_r);
                r_parts.is_some_and(|rp| rhs.connecting_multicenter_bond(rp.iter().copied()).is_some())
            };
            if !in_k {
                let g_parts = map_participants(l_parts.iter().copied(), &l_to_g);
                if let Some(gp) = g_parts {
                    if let Some(g_idx) = self.connecting_multicenter_bond(gp.iter().copied()) {
                        remove_multicenter.push(g_idx);
                    }
                }
            }
        }
        if !remove_multicenter.is_empty() {
            builder.remove_multicenter_bonds(&remove_multicenter);
        }

        // Phase 3: remove L\K noncovalent bonds
        let mut remove_noncovalent: Vec<NoncovalentBondIdx> = Vec::new();
        for nv in lhs.noncovalent_bonds().iter() {
            let in_k = k_l.contains(&nv.atoms[0])
                && k_l.contains(&nv.atoms[1])
                && rhs
                    .connecting_noncovalent_bond(l_to_r[&nv.atoms[0]], l_to_r[&nv.atoms[1]])
                    .is_some();
            if !in_k {
                if let (Some(&ga), Some(&gb)) = (l_to_g.get(&nv.atoms[0]), l_to_g.get(&nv.atoms[1]))
                {
                    if let Some(g_idx) = self.connecting_noncovalent_bond(ga, gb) {
                        remove_noncovalent.push(g_idx);
                    }
                }
            }
        }
        if !remove_noncovalent.is_empty() {
            builder.remove_noncovalent_bonds(&remove_noncovalent);
        }

        // Phase 4: remove L\K atoms and bonds from the topological graph
        let remove_atoms_g: Vec<AtomIdx> = (0..lhs.atoms().count())
            .map(AtomIdx::from)
            .filter(|a| !k_l.contains(a))
            .filter_map(|a| l_to_g.get(&a).copied())
            .collect();

        let mut remove_bonds_g: Vec<BondIdx> = Vec::new();
        for bv in lhs.bonds().iter() {
            let in_k = k_l.contains(&bv.src)
                && k_l.contains(&bv.tgt)
                && rhs
                    .connecting_bond(l_to_r[&bv.src], l_to_r[&bv.tgt])
                    .is_some();
            if !in_k {
                if let (Some(&gs), Some(&gt_)) = (l_to_g.get(&bv.src), l_to_g.get(&bv.tgt)) {
                    if let Some(g_bond) = self.connecting_bond(gs, gt_) {
                        remove_bonds_g.push(g_bond);
                    }
                }
            }
        }

        builder.remove(&remove_atoms_g, &remove_bonds_g);
        Ok(builder.build())
    }
}
