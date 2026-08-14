# Assignment-level aromatic selection

Status: Proposed. 2026-08-13.

## Defect

`AromaticityResolver::select` (the doc 194 S4c/S4d joint selection) accepts aromatic systems
per candidate **member set**, accumulated across the joint valence assignments it enumerates.
Two different assignments of one flexible atom can each produce a valid system with different
members — for quinoline `c1ccc2ccccc2n1`, the pyridinic split (N `#h0`, contributing 1) yields
the 10-atom fused system while the pyrrolic split (N `#h1`, contributing 2) breaks the
azine ring's count and yields the benzo 6-ring alone. Both land in `per_system` under distinct
keys; the sequential-narrowing `retain` checks only that each option's restriction forms are
still offered — the fusion carbons carry the same form in both — so **both systems are
accepted**, they overlap on the fusion atoms, and the commit trips the editor invariant:

```
invalid molecule editor state: aromatic systems: overlap on atom AtomId(3)
```

294 of the 9151 OpenSMILES corpus molecules panic this way; every one is a fused
heteroaromatic with a bare aromatic nitrogen (quinoline, benzimidazole, acridine cores). The
defect is latent since the joint selection landed: no lib or conformance input ingests a
bare-`n` fused system (naphthalene has no flexible atom; indole's `[nH]` is pinned), and it
surfaced only when the substructure bench's corpus loader stopped discarding the panic
(2026-08-13). The `claimed` set exists in `select` but gates only the unclaimed-carrier check,
not acceptance.

The correction: selection chooses among **assignments** — each already a consistent partition
— never mixing systems across them.

## Design

Settled 2026-08-13 except the items marked open.

**Assignment.** One completion choice per flexible atom of a component (see factorization
below), together with the aromatic systems perception finds under that narrowing. The
recorded restriction covers system members only, as today: non-member flexible atoms stay
plural and fall through to the resolver's finalization tie-break.

**Validity.** An assignment is valid iff:

1. its systems are consistent with every stored aromatic system (a stored system must
   reappear with the same member set; an assignment whose partition conflicts with a stored
   system is invalid — the S4d commit-side suppression of already-stored member sets is
   unchanged, one level down);
2. under the `Error` value of `aromatic_valence_failure`, it leaves no unclaimed carrier —
   no atom whose every remaining disjunct requires aromaticity but which no system of the
   assignment claims. Under a non-`Error` policy this condition does not invalidate; the
   tolerated carriers keep their assertions, as in the S4h discharge semantics.

For quinoline this is decisive: the pyrrolic assignment leaves the azine-ring carbons as
unclaimed all-aromatic carriers, is eliminated under the default policy, and the pyridinic
assignment is the unique survivor — no tie-break consulted.

**Selection among valid assignments.**

1. **Structural order (`AromaticityTieBreak`, new policy).** Reachable only under a
   non-`Error` failure policy, where partial-fit assignments survive validity. The ordering
   is lexicographic by: unrealized-carrier count ascending (full 4n+2 fit first), then
   claimed-atom count descending (coverage). Under `Error` every valid assignment has zero
   unrealized carriers and the coverage comparison is the only structural component that can
   fire.
2. **Value order.** Assignments surviving the structural order compare by the existing
   `ValenceTieBreak` key, member-wise over the union of their restrictions, lexicographic by
   atom id — the S4b comparison lifted from per-system options to assignments. `Strict`
   leaves a surviving tie plural: the molecule is `Underdetermined` and the report projects
   the difference onto the differing atoms' candidate sets (bare-`n` imidazole `c1cncn1` is
   the witness: the two tautomeric assignments tie structurally, `MostSaturated` picks one,
   `Strict` reports both nitrogen splits).
3. **Canonical partition pick.** Assignments identical in restriction but differing in
   partition (degenerate symmetry) leave no chemical observable apart — the pick is
   representation canonicalization, not policy: lexicographically smallest partition, member
   sets sorted, partitions compared as sorted lists of sorted sets. Applies under every
   policy, `Strict` included, exactly as canonical labeling picks among automorphic
   labelings.

The winner's systems are accepted as a whole; its restrictions narrow the carrier;
`tie_breaks` records atoms narrowed by the value key, as today.

**Component factorization.** The 4n+2 constraint couples atoms only within a connected
component of the aromatic-candidate graph: build the perception ring set once, keep the
candidate rings (every member aromatic-capable), and take connected components over shared
atoms. Enumeration, validity, selection, and the bound are all per component; component
results combine freely (no joint constraint spans components — biphenyl's rings, or any
multi-fragment molecule). The assignment count drops from the product over the molecule to
the maximum over components.

**Search.** Within a component, replace the flat odometer with a depth-first search over
flexible atoms ordered ring-by-ring, pruning a partial assignment when the violation is
certain:

- a candidate ring with all members assigned whose π sum misses every applicable 4n+2 value,
  when the assignment cannot become valid without it (under `Error`: some fully-assigned
  all-aromatic atom has no other candidate ring left);
- interval bounds: each unassigned member carries a minimum and maximum π contribution over
  its remaining disjuncts; a ring whose reachable sum interval excludes every 4n+2 value is
  settled early.

Pruning must be sound — it may never eliminate a valid assignment — and the executable
specification is equivalence: on components within the bound, the pruned search returns
exactly the flat enumeration's valid-assignment set (property test).

**Bound.** `MAX_JOINT_ASSIGNMENTS` (4096) becomes per-component. Exceeding it returns
`Underdetermined` with the report naming the component's atoms — a stated bound, never
sampling. Symmetry reduction is not available at resolution (graph symmetry is
validation/canonicalization machinery) and is not used.

## Model surface

```rust
// umol-graph::ops::model
/// Disposal policy for structurally distinct valid aromatic assignments.
pub enum AromaticityTieBreak {
    /// No structural preference: structurally distinct survivors stay plural.
    Strict,
    /// Full 4n+2 fit first (unrealized carriers ascending), then coverage
    /// (claimed atoms descending).
    MostAromatic,
}

pub struct AromaticityModel {
    pub scope: ElementScope,
    pub rule: AromaticityRule,
    pub tie_break: AromaticityTieBreak,   // new field
}
```

The policy is a modeling choice and sits in the model envelope beside the valence tie-break,
per the doc 194 envelope design. Breaking for `AromaticityModel` constructors, presets, and
the umol-py bindings.

Open (naming and defaults, for review):

- the variant name `MostAromatic` (parallel to `MostSaturated`; "aromatic" here means
  fit-then-coverage, not a scope claim);
- the default: `Strict` is the neutral choice and is inert under `Error` policies (validity
  already forces full fit; only same-coverage ties remain, which `Strict` leaves plural);
  whether the `daylight`/`mdl`/`permissive` presets opt into `MostAromatic`;
- whether the conformance matrix (194 S5c2) gains this axis — under `Error` policies the
  cells would be identical, so the current position is no.

## Corpus

New resolution-conformance category exercising fused and coupled aromatics; each case feeds
the S5c2 matrix cells on regeneration.

| case | pins |
| --- | --- |
| naphthalene, anthracene, tetracene | linear fusion, single maximal system, no flexible atom |
| phenanthrene | angular fusion |
| biphenyl | ring coupling through a non-ring bond: two components |
| coronene | superring: inner/outer circuits under the fused-ring limits |
| azulene | non-alternant envelope count |
| quinoline, isoquinoline | the defect's minimal repro; carrier check eliminates the pyrrolic assignment |
| benzimidazole, `CCCSc1ccc2c(c1)[nH]c(n2)NC(=O)OC` | corpus panic representatives |
| acridine | aza-linear fusion |
| carbazole, dibenzofuran | five-ring bridges |
| purine | fused bare-`n` with a genuine tautomeric tie: value tie-break witness |
| pteridine | several flexible atoms in one component |
| porphine | macrocycle, 2×NH + 2×N: tautomer pair, component-size stress |

## Plan

- A0: component factorization in `select`'s enumeration — candidate-ring components, per-component
  bound; behavior-preserving for single-component inputs. Additive. [dep: none]
- A1: assignment-level validity and selection replacing the `per_system` acceptance:
  validity (stored-system consistency, carrier elimination under `Error`), value tie-break
  over assignments, canonical partition pick. Fixes the overlap panic; quinoline family
  resolves, bare-`n` imidazole under `Strict` goes honestly underdetermined. Breaking
  (outputs for the panicking family change from panic to resolution). [dep: A0]
- A2: `AromaticityTieBreak` in the model envelope with the structural order under
  non-`Error` policies; umol-py bindings; preset decision from the open list. Breaking.
  [dep: A1]
- A3: pruned depth-first search replacing the flat per-component enumeration, with the
  enumeration-equivalence property test. Additive (same results, larger reachable inputs).
  [dep: A1]
- A4: the corpus category plus acceptance tests for quinoline/isoquinoline/purine through
  `ingest_smiles`; the substructure and fingerprint bench corpus loaders run panic-free.
  [dep: A1; before 194 S5d so the one regeneration covers the new category]

The critical path is A0–A1; A2 and A3 are independent of each other. 194 S5d gains the
dependency on A4.
