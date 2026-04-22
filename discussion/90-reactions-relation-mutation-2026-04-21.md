# Reactions and Relation Mutation

## Structural Model

- **Topology** = atoms + localized bonds (sigma scaffold) in the `Graph`
- **Overlays** = aromatic systems, dative bonds, multicenter bonds, noncovalent bonds — relations outside the topological graph

Example: benzene has 6 C atoms and 6 single localized bonds (sigma scaffold) in the topological graph. The pi system is represented as an aromatic system relation over the same 6 atoms with 6 electrons. The aromatic system is not part of the topological graph.

## DPO Rewriting Phases

Given a rule L ← K → R and a match m: L → G on a target molecule G:

1. **Add R \ K**: new atoms, localized bonds, and overlay relations (K atoms at original indices)
2. **Modify K attributes**: attribute changes on preserved entities
3. **Remove overlay relations**: overlay relations in L \ K (aromatic systems, dative bonds, multicenter bonds, noncovalent bonds whose participant atoms are all preserved)
4. **Remove L \ K atoms/bonds**: topological graph changes via `remove(atoms, bonds)`; overlays whose participant atoms are removed are cleaned up automatically

Phases 1–2 operate at original indices. Phase 4 remaps all remaining indices.

```rust
let mut builder = target.edit();

// Phase 1
let new_atom = builder.add_atom(r_atom_data);
builder.add_bond(k_atom_original, new_atom, bond_data);

// Phase 2
builder.atom_mut(k_atom).element = new_element;

// Phase 3
builder.remove_aromatic_system(system_idx);

// Phase 4
let _remap = builder.remove(&l_minus_k_atoms, &l_minus_k_bonds);

let result = builder.build();
```

## Aromatic System Cases

### S_E Ar: Bromination of Benzene

C₆H₆ → C₆H₅Br (substitution at one carbon).

| | L (matched) | K (interface) | R (replacement) |
|---|---|---|---|
| Atoms | 6 C + 1 H | 6 C | 6 C + 1 Br |
| Localized bonds | 6 C-C single + 1 C-H | 6 C-C single | 6 C-C single + 1 C-Br |
| Aromatic system | 6 C, 6e⁻ | 6 C, 6e⁻ | 6 C, 6e⁻ |

- L \ K = {H atom, C-H bond}
- R \ K = {Br atom, C-Br bond}
- Aromatic system is in K → preserved automatically

Phases used: 1 (add Br, C-Br bond), 4 (remove H, C-H bond). No relation operations.

### Hydrogenation: Benzene → Cyclohexane

Full hydrogenation destroys aromaticity without changing the sigma scaffold.

| | L | K | R |
|---|---|---|---|
| Atoms | 6 C | 6 C | 6 C |
| Localized bonds | 6 C-C single | 6 C-C single | 6 C-C single |
| Aromatic system | 6 C, 6e⁻ | — | — |

- L \ K = {aromatic system} — overlay only, no topology change
- R \ K = {} (plus explicit H atoms if modeled)

Phases used: 3 (`remove_aromatic_system`). No atoms or bonds removed from topology.

Without `remove_aromatic_system`, this transformation is inexpressible: `remove([], [])` does nothing, and the aromatic system survives because all its participant atoms are preserved.

### Partial Reduction: Naphthalene → Tetralin

One ring of naphthalene is hydrogenated; the other retains aromaticity.

Representation depends on whether naphthalene has one aromatic system (10 atoms, 10e⁻) or two (6 atoms each, sharing 2 bridgehead atoms, 6e⁻ each).

**Single-system representation:**

| | L | K | R |
|---|---|---|---|
| Atoms | 10 C | 10 C | 10 C |
| Localized bonds | 11 C-C single | 11 C-C single | 11 C-C single |
| Aromatic system | 10 C, 10e⁻ | — | 6 C (ring B), 6e⁻ |

Phases used: 1 (`add_aromatic_system` for ring B), 3 (`remove_aromatic_system` for the 10-atom system).

**Two-system representation:**

Phases used: 3 (`remove_aromatic_system` for ring A only). Ring B's aromatic system is in K, preserved.

## Dative Bond Cases

### Amine Oxide Reduction: R₃N→O → R₃N

| | L | K | R |
|---|---|---|---|
| Atoms | N, O | N | N |
| Localized bonds | N-O sigma (if present) | — | — |
| Dative bond | N→O | — | — |

- L \ K = {O atom, N-O bond if present, dative bond}
- Removing O drops the dative bond automatically (endpoint removed)

Phases used: 4 (remove O, N-O bond). Dative bond cleaned up by atom removal. No explicit relation removal needed.

### Ligand Exchange: M←L₁ → M←L₂

| | L | K | R |
|---|---|---|---|
| Atoms | M, L₁ atoms | M | M, L₂ atoms |
| Localized bonds | M-L₁ sigma (if any) | — | M-L₂ sigma (if any) |
| Dative bond | L₁→M | — | L₂→M |

- L \ K = {L₁ atoms, L₁ bonds, L₁→M dative bond}
- R \ K = {L₂ atoms, L₂ bonds, L₂→M dative bond}
- Removing L₁ atoms cleans up the dative bond automatically

Phases used: 1 (add L₂, `add_dative_bond`), 4 (remove L₁ atoms/bonds). No explicit relation removal needed.

### Lewis Acid-Base Complex Formation: BF₃ + NH₃ → F₃B←NH₃

Both B and N exist before and after. Dative bond forms between preserved atoms.

| | L | K | R |
|---|---|---|---|
| Atoms | B, N | B, N | B, N |
| Dative bond | — | — | N→B |

- L \ K = {} (nothing removed)
- R \ K = {dative bond N→B, possibly B-N sigma bond}

Phases used: 1 (`add_dative_bond`). Only addition, no removal.

### Ligand Dissociation (No Atom Removal)

A dative bond breaks while both atoms remain in the molecule (e.g., in a crystal or supramolecular complex).

| | L | K | R |
|---|---|---|---|
| Atoms | M, L atoms | M, L atoms | M, L atoms |
| Dative bond | L→M | — | — |

- L \ K = {dative bond only}
- All atoms preserved

Phases used: 3 (`remove_dative_bond`). Without it, the transformation is inexpressible.

## Multicenter Bond Cases

### Diborane Bridge Opening: B₂H₆ → 2 BH₃

Diborane has two 3-center-2-electron (3c-2e) B-H-B bridging bonds.

If the bridge opens while all atoms remain in the molecule (e.g., Lewis base cleaves the bridge):

| | L | K | R |
|---|---|---|---|
| Atoms | 2 B, 2 H_bridge | 2 B, 2 H_bridge | 2 B, 2 H_bridge |
| Multicenter bonds | 2 × (B, H, B) 3c-2e | — | — |
| Localized bonds | — | — | 2 × B-H terminal |

Phases used: 3 (`remove_multicenter_bond` × 2), 1 (add terminal B-H bonds). Without `remove_multicenter_bond`, the transformation is inexpressible.

If the molecule actually splits (atoms removed from one fragment), multicenter bond is cleaned up by atom removal.

### Hapticity Change: η⁵ → η¹ in Metallocene

A cyclopentadienyl ring changes from η⁵ coordination (multicenter bond: 5 C + Fe) to η¹ coordination (single Fe-C sigma or dative bond).

| | L | K | R |
|---|---|---|---|
| Atoms | Fe, 5 C | Fe, 5 C | Fe, 5 C |
| Multicenter bond | (Fe, C₁, C₂, C₃, C₄, C₅) | — | — |
| Localized/dative bond | — | — | Fe-C₁ sigma or C₁→Fe dative |

Phases used: 3 (`remove_multicenter_bond`), 1 (add localized or dative bond). All atoms preserved.

## Noncovalent Bond Cases

### Hydrogen Bond Breaking via Proton Transfer

D-H···A → D⁻ + H-A. The proton transfers from donor to acceptor, breaking the hydrogen bond.

| | L | K | R |
|---|---|---|---|
| Atoms | D, H, A | D, H, A | D, H, A |
| Localized bonds | D-H | — | H-A |
| Noncovalent bond | H-bond (D, A) | — | — |

- All atoms preserved; noncovalent bond endpoints are atoms, not bonds
- `remove([], [D-H bond])` does not clean up the noncovalent bond — both endpoint atoms survive

Phases used: 1 (add H-A bond), 3 (`remove_noncovalent_bond`), 4 (remove D-H bond).

### Halogen/Chalcogen Bond Disruption

Same pattern as hydrogen bond breaking: if both atoms remain in the molecule but the noncovalent interaction is broken by electronic or geometric changes, explicit removal is needed.

## Summary

| Case | Atoms removed? | Overlay change | Needs explicit removal? |
|---|---|---|---|
| S_E Ar (bromination) | H removed | aromatic preserved | No |
| Hydrogenation | none | aromatic removed | **Yes** |
| Partial reduction (naphthalene) | none | aromatic removed + added | **Yes** |
| Amine oxide reduction | O removed | dative cleaned up | No |
| Ligand exchange | L₁ removed | dative cleaned up | No |
| Lewis complex formation | none | dative added | No |
| Ligand dissociation (no atom loss) | none | dative removed | **Yes** |
| Diborane bridge opening | none | multicenter removed | **Yes** |
| Hapticity change (η⁵ → η¹) | none | multicenter removed | **Yes** |
| H-bond breaking (proton transfer) | none | noncovalent removed | **Yes** |

**Pattern:** explicit relation removal is needed whenever an overlay relation is destroyed but all its participant atoms survive. This occurs for all four relation types in realistic chemistry.

## Required API: Full Operations Inventory

### Already present

**MoleculeBuilder — structural addition (Phase 1):**
- `add_atom(AtomAst) -> AtomIdx`
- `add_bond(src, tgt, BondAst) -> BondIdx`
- `add_dative_bond(donor, acceptor, DativeBondAst) -> DativeBondIdx`
- `add_aromatic_system(atoms, AromaticSystemAst) -> AromaticSystemIdx`
- `add_multicenter_bond(atoms, MulticenterBondAst) -> MulticenterBondIdx`
- `add_noncovalent_bond(ends, NoncovalentBondAst) -> NoncovalentBondIdx`

**MoleculeBuilder — constraint addition (Phase 1):**
- `push_atom_constraint`, `push_bond_constraint`
- `push_dative_bond_constraint`, `push_aromatic_system_constraint`
- `push_multicenter_bond_constraint`, `push_noncovalent_bond_constraint`
- `push_molecule_constraint`

**MoleculeBuilder — topological removal (Phase 4):**
- `remove(atoms, bonds) -> IdxRemapping` — removes atoms and localized bonds from the graph, auto-cleans overlay relations whose participant atoms are removed, remaps constraints, returns full index remapping.

**MoleculeAst — read operations for gluing condition check:**
- `neighbors()`, `connecting_bond()`
- `dative_bonds_incident()`, `aromatic_systems_incident()`
- `multicenter_bonds_incident()`, `noncovalent_bonds_incident()`

### Missing on MoleculeBuilder

**Attribute mutation (Phase 2):**

Attribute mutation exists on `MoleculeAst` (`atom_mut`, `bond_mut`, `dative_bond_mut`, `aromatic_system_mut`, `multicenter_bond_mut`, `noncovalent_bond_mut`) but not on `MoleculeBuilder`. For DPO, K-entity attributes must be modified as part of the rewrite, not on the source molecule before `edit()`.

- `atom_mut(AtomIdx) -> &mut AtomAst`
- `bond_mut(BondIdx) -> &mut BondAst`
- `dative_bond_mut(DativeBondIdx) -> &mut DativeBondAst`
- `aromatic_system_mut(AromaticSystemIdx) -> &mut AromaticSystemAst`
- `multicenter_bond_mut(MulticenterBondIdx) -> &mut MulticenterBondAst`
- `noncovalent_bond_mut(NoncovalentBondIdx) -> &mut NoncovalentBondAst`

**Constraint access (Phase 2):**

`Constraints` has `remove_*`, `retain_*`, and `clear` methods, but `MoleculeBuilder` only exposes `push_*`. DPO needs to remove L\K constraints and modify K constraints.

- `constraints_mut(&mut self) -> &mut Constraints`

**Overlay relation removal (Phase 3):**

Batched per relation kind. Each method takes a slice of indices, removes them, shifts remaining indices within that kind, and remaps associated constraints. Does not affect atom, bond, or other relation-kind indices.

- `remove_dative_bonds(&mut self, indices: &[DativeBondIdx])`
- `remove_aromatic_systems(&mut self, indices: &[AromaticSystemIdx])`
- `remove_multicenter_bonds(&mut self, indices: &[MulticenterBondIdx])`
- `remove_noncovalent_bonds(&mut self, indices: &[NoncovalentBondIdx])`

### Missing on MoleculeAst

**Induced subgraph (separate from DPO):**

- `induced_subgraph(atoms: &[AtomIdx]) -> MoleculeSubgraph` — extracts the fragment induced by the given atoms, including all localized bonds between them, all overlay relations whose participants are all in the set, and remapped constraints. Returns new-to-old index maps for all six entity types.

### Complete DPO Sequence

```rust
// Gluing condition check (on source MoleculeAst)
for &atom in &l_minus_k_atoms {
    for n in target.neighbors(atom) {
        assert!(l_atoms.contains(&n.atom), "dangling edge");
    }
    // same for dative_bonds_incident, aromatic_systems_incident, etc.
}

let mut builder = target.edit();

// Phase 1: add R \ K (K atoms at original indices)
let new_atom = builder.add_atom(r_atom_data);
builder.add_bond(k_atom, new_atom, bond_data);
builder.add_aromatic_system(r_system_atoms, system_data);
builder.push_atom_constraint(new_atom, constraint);

// Phase 2: modify K attributes
builder.atom_mut(k_atom).charge = new_charge;
builder.bond_mut(k_bond).order = new_order;
builder.constraints_mut().remove_atom(k_atom);

// Phase 3: remove overlay relations from L \ K
builder.remove_aromatic_systems(&l_minus_k_aromatic_systems);
builder.remove_dative_bonds(&l_minus_k_dative_bonds);

// Phase 4: remove L \ K atoms/bonds from topology
let remap = builder.remove(&l_minus_k_atoms, &l_minus_k_bonds);

let result = builder.build();
```
