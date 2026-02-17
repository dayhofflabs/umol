# Molecule Validation Scheme

## Prompt

I'd like to talk through a validation approach for molecular graphs that works in phases. Input is a molecular graph notation like CTab or SMILES that is parsed into a list of atoms (incl aromaticity and chiral flags) and bonds (incl bond order, aromaticity, and bond direction flags). The validation / resolution approach consists of topology, valence, aromaticity, and stereo phases. 1. Topology checks for disconnected components, parallel edges, self-loops (easy). 2. Valence compares the atomic features to a set of predefined validated atomic valence states (which include charge, number of unpaired electrons, lone pairs, and bond sum (valence). 3. Aromaticity constructs representations from the atomic valence states and aromatic flags for atoms and bonds. 4. Stereo does the same for the chiral and bond stereo flags. Now my question. Aromatic bonds do not allow to define a bond order so that one needs to defer the bond sum check (total number of bonds is not known at that phase). But if the aromaticity detection fails, that means that one would need to backtrack. Is there a way of avoiding backtracking? Or is it actually acceptable? One thought I had was to extend the atomic valence state to include the number of electrons donated to an aromatic system, which would make it a bit cleaner (also specify which elements we allow to be involved in aromaticity). But this would not entirely avoid backtracking either, I'm afraid. If we need backtracking, what is a conceptually clean way of doing that?

## Response

The backtracking problem disappears if you reframe phases as narrowing a **set of candidate states** rather than committing to a single assignment.

### The core trick: separate σ-skeleton from π-system

For aromatic atoms, you *can* check a partial valence constraint at phase 2 — you just check against the σ-skeleton. Every aromatic bond contributes exactly one σ-bond regardless of Kekulé assignment. What's unknown is only the π-electron contribution per atom. So define your aromatic valence states as:

```
(charge, unpaired_e, lone_pairs, σ_bond_sum, π_contribution)
```

where `π_contribution ∈ {0, 1, 2}` depending on atom type. For example:

- Aromatic C (neutral): σ_bonds=3, π=1
- Pyridine N: σ_bonds=2, π=1  
- Pyrrole N: σ_bonds=2 (+ one H), π=2
- Furan O: σ_bonds=2, π=2

Phase 2 now enumerates the **set of valid aromatic valence states** per atom (usually 1–2 candidates). Phase 3 selects from those candidates by checking whether a consistent assignment satisfying 4n+2 per ring exists. If the set becomes empty → invalid molecule. No backtracking, because you never committed.

### The pipeline type

Each phase is conceptually:

```
Phase: Set<CandidateState> → Set<CandidateState>
```

Composition is just `phase4 ∘ phase3 ∘ phase2 ∘ phase1`, and `∅` at any point is a short-circuit failure. This is exactly the nondeterminism monad (`Vec<T>` / `List` monad) — or if you want to track *why* something was eliminated, use `Vec<Result<T, ValidationError>>`.

In Rust terms, each atom carries a `SmallVec<[ValenceState; 4]>` of surviving candidates through the pipeline. Aromaticity resolution becomes constraint propagation over these small domains — which for typical organic molecules is trivially small.

### Why this is clean

- No backtracking: you defer *commitment*, not *checking*
- Valence phase still does real work (rejects impossible states early)
- Aromaticity phase has a well-defined contract: select from pre-validated candidates
- The aromatic valence states also naturally encode your "allowed aromatic elements" — if an element has no aromatic valence state defined, it can't be aromatic

The only subtlety: if you want to propagate *which* Kekulé assignment was chosen back for downstream use (stereo phase might care about double bond placement), phase 3 should output that as part of its refined state.

## Reply

That's good. Now if I want to go beyond organic chemistry alone and include transition metals, do I need to separate further into other types of contributions (delta bonds, dative bonds (think Cp^-)?

## Response

Honest answer: for TMs, the clean σ/π decomposition breaks down not because you need more orbital symmetry types, but because the **electron counting model itself changes**.

### What actually differs

In organic chemistry, every bond is a shared-electron-pair bond, and the valence state is just bookkeeping over those pairs. For TMs:

- **Dative bonds** break the symmetry of contribution — donor gives 2e, acceptor gives 0. Your valence state needs to track electrons *received* vs *shared* separately.
- **Hapticity** (η⁵-Cp⁻, η⁶-arene) means a single "bond" in the graph-theoretic sense is really a multi-center interaction. This is a *representation* problem before it's a validation problem.
- **δ bonds** are honestly a non-issue for validation — they only appear in metal-metal multiple bonds and can be folded into bond order.

### The real fork in the road

You have two choices:

**1. Extend the valence state formalism** with a TM-specific electron counting scheme (CBC/Green's method or ionic formalism):

```
(oxidation_state, d_electron_count, coordination_number, total_e_count)
```

Validation becomes: does the total electron count satisfy 18e (or a known exception)? This is clean but requires bond type annotations (covalent vs dative) in the input, which SMILES/CTab handle poorly.

**2. Punt on TM valence validation.** This is what RDKit effectively does, and it's not unreasonable — the combinatorial space of valid TM states is enormous, the input representations are lossy, and false negatives are worse than no check.

### My suggestion

Keep your phase architecture but make TM validation a **separate, opt-in phase** with a different contract. Organic validation narrows candidate sets toward a unique assignment. TM validation is better framed as *consistency checking* against a permissive rule set (electron count in [12,18], coordination number in allowed range for geometry, etc.) rather than exhaustive enumeration of valid states.

The Cp⁻ case specifically: if your graph representation uses individual M–C edges for haptic bonds, you need a "haptic bond group" annotation so the valence phase knows those 5 edges collectively donate a certain number of electrons rather than each being an independent 1e or 2e bond.

## Reply

What's the Green's formalism?

## Response

Malcolm Green's Covalent Bond Classification (CBC) method. Every metal-ligand interaction is classified into three types based on the number of electrons the ligand contributes to the bond:

- **X-type**: ligand contributes 1e, metal contributes 1e (normal covalent bond). Cl⁻, alkyl, H.
- **L-type**: ligand contributes 2e, metal contributes 0e (dative). CO, PR₃, NH₃.
- **Z-type**: ligand contributes 0e, metal contributes 2e (reverse dative). BF₃, Lewis acids.

A complex is then classified as `[MLₗXₓZᵤ]`, and the electron count is:

```
total_e = metal_d_electrons + 2l + x
```

where the metal's d-electron count comes from the "equivalent neutral class" (you formally remove all ligands homolytically for X, heterolytically for L and Z, and count what's left on the metal).

The elegance is that it's **oxidation-state-free** — you don't need to assign formal charges to get the electron count right. The classification is a property of the bond, not of the ligand in isolation (Cp can be L₂X or other decompositions depending on hapticity).

For your validation framework: each edge in the molecular graph gets typed as L, X, or Z, each contributing a known electron count. This maps naturally onto your candidate-set approach — a TM atom's valid states become `{(l, x, z) : 2l + x + metal_e ∈ {16, 18}}` modulo known exceptions. The combinatorial blowup you'd fear is actually modest because L/X/Z typing is usually unambiguous from the ligand identity.

## Reply

Can you point me to literature?

## Response

The key papers:

1. **The original**: M.L.H. Green, "A new approach to the formal classification of covalent compounds of the elements," *J. Organomet. Chem.* **500**, 127–148 (1995). doi:10.1016/0022-328X(95)00508-N

2. **The pedagogical treatment** (probably the most useful for your purposes): M.L.H. Green & G. Parkin, "Application of the Covalent Bond Classification Method for the Teaching of Inorganic Chemistry," *J. Chem. Educ.* **91**, 807–816 (2014). doi:10.1021/ed400504f

3. **Extension to 3-center 2-electron bonds**: J.C. Green, M.L.H. Green & G. Parkin, *Chem. Commun.* (2012), and later expanded as a chapter in *Structure and Bonding* **171** (Springer, 2016). doi:10.1007/430_2015_206

The 2014 J. Chem. Ed. paper is the one I'd start with — it lays out the full `[MLₗXₓZᵤ]` formalism with worked examples and electron counting rules in a way that maps most directly onto what you'd need to encode as validation rules.

## Reply

Regarding punting: RDKit is not a model for me. It is its tendency to take the easy way out and lack of principled impls that gave me endless problems having to work around its bugs and shortcomings.

## Response

Fair enough. Then let me reconsider the TM validation properly.

The candidate-set architecture actually scales to this. The key insight is that CBC gives you a **finite, enumerable** set of valid `[MLₗXₓZᵤ]` configurations per metal, and the constraints are well-defined:

1. Each edge is typed L, X, or Z (determined by ligand identity — this is usually unambiguous)
2. Valid states satisfy: `2l + x + metal_group_electrons ∈ {12, 14, 16, 18}` (not just 18 — you need the known exceptions)
3. Coordination number and d-electron count must be consistent with known geometries

The combinatorial explosion people fear is a myth for practical chemistry. For a given metal center, once you know the ligand set, the number of valid CBC assignments is small — usually 1, occasionally 2–3 when a ligand is ambiguous (e.g., NO can be L or X depending on bending).

The hard part is **haptic ligands**, because they break the 1-edge-per-interaction assumption. But this is a representation problem you solve once: a haptic group annotation (η⁵-Cp is `L₂X`, η⁶-arene is `L₃`) collapses the multi-edge subgraph into a single CBC contribution. Your graph needs a hyperedge or group annotation for this — individual M–C edges without it are simply not enough information to validate against.

So the honest answer to your earlier question: no, you don't need δ-bond or dative-bond edge types as separate primitives. You need **CBC-typed edges** (L/X/Z) and **haptic group annotations**, and then the same candidate-set narrowing pipeline works for TMs as for organic molecules, just with a different rule set per phase.