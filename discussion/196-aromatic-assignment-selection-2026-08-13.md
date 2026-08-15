# 196 — Assignment-level aromatic selection

Status: Proposed
Date: 2026-08-13
Relates: [174](174-aromatic-hydrogen-resolution-2026-07-31.md),
[194](194-constraint-assertion-semantics-2026-08-10.md)

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

**Closed-world admission (`asserted_complete`).** Settled 2026-08-13. Resolution converts
existing constraints into entities: in resolution, the absence of a constraint is read as
actual absence, not indeterminacy. An atom with no aromatic mark in either placement — its
own `#a` assertion (the SMILES placement) or an incident bond's `#a` mark (the MOL
placement) — admits only non-aromatic completions; there is no aromatize-or-not choice for
unmarked rings, and perception-from-silence stays with the aromatizer transform. The
reading gets its own keyed accessor on the constraints views, `asserted_complete(key)`:
the asserted side under resolution's closed-world claim, where the absence cell of an
entity-creating overlay key closes to its definite negative (`NotAromatic`, `NotStereo`).
`asserted(key)` keeps the open-world reading for matching, where a pattern constrains only
what it mentions (recorded in doc 195). Admission's aromatic-evidence check becomes
`asserted_complete(AromaticValence)` over both placements. The four readings stay four
flat accessors — each consumer picks exactly one statically, and the two closures are
justified by different claims, so no parameterized accessor unifies them.

**Validity.** An assignment is valid iff:

1. its systems are consistent with every stored aromatic system (a stored system must
   reappear with the same member set; an assignment whose partition conflicts with a stored
   system is invalid — the S4d commit-side suppression of already-stored member sets is
   unchanged, one level down). Stored systems are facts in resolution (settled 2026-08-13):
   when no valid assignment remains because of a stored-system conflict, the outcome is
   `Contradictory` naming the stored system under `aromatic_system_failure: Error`; under a
   non-`Error` value the component accepts no new systems and tolerated carriers keep
   their assertions. The `ReplaceEntity`/`RemoveConstraint` machinery stays with the plan
   path and never enters selection;
2. under the `Error` value of `aromatic_valence_failure`, it leaves no unclaimed carrier —
   no atom whose every remaining disjunct requires aromaticity but which no system of the
   assignment claims. Under a non-`Error` policy this condition does not invalidate; the
   tolerated carriers keep their assertions, as in the S4h discharge semantics.

For quinoline this is decisive: the pyrrolic assignment leaves the azine-ring carbons as
unclaimed all-aromatic carriers, is eliminated under the default policy, and the pyridinic
assignment is the unique survivor — no tie-break consulted.

**Selection among valid assignments.**

1. **Structural order (`AromaticityTieBreak`, new policy).** Reachable only under a
   non-`Error` failure policy, where partial-fit assignments survive validity. One
   criterion (corrected 2026-08-13: under closed-world admission, claimed atoms and
   unrealized carriers partition the evidenced set, so "full fit first" and "coverage" are
   complements — the earlier two-component phrasing was redundant): claimed-atom count
   descending, then the member sets themselves (lexicographic) — so the value order below
   never compares assignments restricting different atoms. Ring size and system
   granularity carry no preference here: what counts as a system is the perception's
   decision under `AromaticityRule` and `RingLimits`, one decomposition per assignment;
   selection never re-partitions. Under `Error` the order is fully inert: validity forces
   zero unrealized carriers, and closed-world admission fixes the claimable atoms —
   marked atoms must be claimed, unmarked atoms never can be — so valid assignments
   cannot differ in coverage or member sets.
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

The winner's systems are accepted as a whole; its restrictions narrow the carrier.
`tie_breaks` generalizes from "the value key acted here" to "a configured policy acted
here" (settled 2026-08-13): it records atoms narrowed by the value key, as today, and the
chosen systems' member atoms when the structural order decided between valid assignments.
Validity elimination is a criterion, not a preference, and never records.

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

**Bound.** `MAX_JOINT_ASSIGNMENTS` (4096) becomes per-component and is renamed
`MAX_ASSIGNMENTS` (settled 2026-08-13): an assignment is by definition one choice per
flexible atom of the component, so the "joint" qualifier distinguishes nothing. It bounds
the assignment
enumeration and is distinct from `RingLimits`, which bounds the perception's ring
enumeration (a marked macrocycle beyond the ring limits produces no candidate rings and
falls to the carrier check). Exceeding the assignment bound makes the whole molecule
`Underdetermined` (settled 2026-08-13): nothing commits, per the single-commit doctrine —
resolved sibling components included — and the report carries the over-bound component's
atoms, since resolved components' entries are singleton and the report lists plural
survivors. A stated bound, never sampling. Symmetry reduction is not available at
resolution (graph symmetry is validation/canonicalization machinery) and is not used.

## Model surface

```rust
// umol-graph::ops::model
/// Disposal policy for structurally distinct valid aromatic assignments.
#[derive(Default)]
pub enum AromaticityTieBreak {
    /// No structural preference: structurally distinct survivors stay plural.
    /// The default, as for `ValenceTieBreak`: any other default biases the
    /// scheme.
    #[default]
    Strict,
    /// Maximal claimed-atom coverage; the perception decides what a system
    /// is, the policy only orders how much of the evidence is realized.
    MaxCoverage,
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

Naming and defaults, settled 2026-08-13:

- the variant name is `MaxCoverage` (`MostAromatic` was rejected — "aromatic atom" reads as
  aromaticity being the atom's property);
- the default is `Strict`, as for the valence tie-break: any other default would bias the
  scheme;
- the `daylight`/`mdl`/`permissive` presets do **not** opt into `MaxCoverage`; they inherit
  the default. The presets reproduce a format convention's reading, and on this axis there
  is no reference behaviour to reproduce: every reference toolkit pins implicit hydrogens at
  parse time from its fixed valence table, so perception runs over one already-fixed valence
  state and no set of assignments arises. Where a fixed reading is impossible — bare `n` in
  a pyrrole position — RDKit refuses the input rather than guessing (doc 174). The
  valence-level half of that convention is already carried by `ValenceModel::smiles()` with
  `MostSaturated`. Under `Error` failure policies the structural order is inert regardless,
  so opting in would change nothing except in Keep-mode, where it would bias exactly the
  case the caller is inspecting;
- the conformance matrix (194 S5c2) does not gain this axis: under `Error` policies the
  cells would be identical.

Contingent, not open: whether A5's class distribution argues for a shipped bulk suite after
all, and if so the category names for it.

## Corpus

Assignment selection sits on the main resolution path, so it needs validation at a scale no
hand-curated suite reaches. The redistribution question decides the shape (settled
2026-08-13): the bulk corpus is a **local development instrument**, not a shipped suite.
What ships is a curated set promoted from what the instrument finds.

**Why not vendored.** The ChEMBL FTP licence (CC BY-SA 3.0) covers "the data content in
ChEMBL", and VEHICLe is not ChEMBL content — the ChEMBL repository states only that the
regids will be integrated "at some point in the future", and the VEHICLe FTP README carries
no licence statement, just column definitions and a UCB contact. Sitting in a sibling
directory of the same FTP host is not a grant. Even the favourable reading is share-alike,
which would put a copyleft obligation and a separate licence file inside an MIT/Apache
repository. UK/EU sui generis database right is the reason the "exhaustive enumerations are
uncopyrightable facts" argument does not dispose of the question; that right runs 15 years
from publication and the file is dated 2010-04-12, so it has most likely lapsed — "most
likely" being the wrong footing for a redistribution decision. Downloading for local testing
raises none of this.

**The instrument.** A script fetches VEHICLe (24 867 enumerated 5- and 6-membered mono- and
bicyclic aromatic ring systems over C/N/O/S; Pitt et al., *J. Med. Chem.* 2009,
[10.1021/jm801513z](https://doi.org/10.1021/jm801513z); ChEMBL FTP
`pub/databases/chembl/VEHICLe/`) into a gitignored working directory. Nothing downloaded is
committed. It is the defect's own chemical space: 4 221 rows carry a bare aromatic nitrogen,
1 006 carry `[nH]`, and 1 544 tautomer clusters group rows that must resolve alike.

**What it measures.** Each row resolves under both valence tie-breaks, and the outcome pair
classifies it:

| `Strict` | `MostSaturated` | class |
| --- | --- | --- |
| `Determined` | `Determined`, equal | one resolution, no discretion |
| `Underdetermined` | `Determined` | plural under `Strict`, tie-break picks |
| `Underdetermined` | `Underdetermined` | tie survives the key |
| `Contradictory` | `Contradictory` | fails at every level |
| anything else | | hierarchy violation — `MostSaturated` only ever narrows further, so this combination is impossible and indicts the implementation |

The class distribution answers whether perception and selection work well: contradictions
concentrated in chemically implausible rows is health, contradictions scattered through
ordinary heteroaromatics is not. Two invariants are checked on every row regardless of
class: every emitted system fits 4n+2, and no two systems overlap — the invariant this
document's defect violates. Tautomer clusters are the second oracle: rows in one cluster
should share a class, and under `MostSaturated` should pick consistently.

**What ships.** From the instrument's output, a curated set is promoted into the doc 194
S5c2 resolution conformance suite, where the four-cell matrix covers each case across both
candidate sources: the panic representatives (quinoline, isoquinoline, benzimidazole), at
least one tautomer cluster carried in full so the `Strict`-plural and `MostSaturated`-pick
behaviours are both pinned, purine and pteridine for fused bare-`n` ties, and whatever the
distribution surfaces as unexpected. Alongside them, hand-written cases covering what
VEHICLe structurally cannot supply — it stops at mono- and bicyclics: polycyclic aromatics
(acenes, phenanthrene, coronene, azulene), large macrocycles (porphine, phthalocyanine),
and coupled or bridged systems (biphenyl, carbazole, dibenzofuran, acridine).

The instrument is then expected to be deleted rather than kept: it is a one-time validation
of the algorithm, and the promoted cases carry the regression coverage. A dedicated shipped
`aromaticity_perception` suite (raw + cleaned + snapshots, shaped like the SMILES-parsing
suite of doc 047) stays available as an option if the instrument shows the curated set
cannot cover the behaviour, and would then need a corpus whose redistribution is settled —
an independent enumeration of the same space, generated in-repo, being the obvious candidate.

## Plan

Stages keep the A0–A6 identities referenced from doc 194 and the status table. Every stage
ends green; a subitem may go red only where marked breaking, with its caller migration in
the same stage. Subitems are ordered foundation-first: graph-ir views, then admission, then
selection, then the model envelope and bindings, then the suites.

### A0 — closed-world reading and admission (breaking)

- A0a: `asserted_complete(key)` on all eight constraints views
  (`umol-graph-ir/src/ir/view/constraints.rs` and the per-entity derivation modules).
  Signature per view mirrors `asserted`: `pub fn asserted_complete(&self, key: K) ->
  Option<Form>` (owned, since closure can synthesize). Contract: the stored assertion when
  present; otherwise the same absence cells as `derived_complete` — atom `AromaticValence`
  → `NotAromatic`, `TetrahedralStereo` → `NotStereo`, `MulticenterValence` →
  `NotMulticenter`, dative pair counts → zero, bond `Aromatic` → `Lit(false)`,
  `CisTransStereo` → `NotStereo`; topology keys return the assertion or `None` (no absence
  cell, as on the derived side). One placement merge: atom `AromaticValence` reads the
  atom's own assertion, else an incident bond asserting `#a` as `Aromatic(Undetermined)`,
  else the closure — the two placements are one evidence notion. Determination-style table
  tests per view (present / absent / bond-adjacent rows). Additive. [dep: none]
  **Done 2026-08-13:** the accessor on all eight views (atom and bond with their absence
  cells and the atom-side placement merge; dative with the aromatic cell; the five
  no-absence-cell families as the cloned assertion, documented as such). Twenty-five table
  rows pin the contract, including the two semantic edges: bond-adjacent evidence reads
  `Aromatic(Undetermined)`, and an atom or bond inside a stored system without its own
  assertion still reads the negative — the accessor never reads relations. The
  nomenclature guide's constraint-readings entry now lists all four accessors with the
  consumer table (matching → `asserted`, admission → `asserted_complete`). Normalized to
  the facade rule on review: the three nontrivial closures live beneath the views as
  `atom_`/`bond_`/`dative_bond_asserted_complete_constraint` in the derivation modules,
  and every view method is a one-line delegation like its `asserted` sibling. graph-ir
  6063 green, clippy zero.
- A0b: counts admission on the closed-world reading
  (`umol-graph/src/ops/valence/counts.rs`). The aromatic-evidence check becomes
  `asserted_complete(AromaticValence)`; unmarked atoms take only the non-aromatic
  enumeration; field-ground atoms whose evidence is present but contribution undetermined
  are admitted with candidates varying only in the table's aromatic valences.
  `candidate_states` table tests migrate: unmarked rows lose aromatic candidates,
  ground-evidenced rows are new. Breaking. [dep: A0a]
- A0c: atom-typing admission on the same reading
  (`umol-graph/src/ops/valence/atom_typing.rs`). Unmarked atoms no longer admit aromatic
  registry rows (the mixed-candidate case dissolves); field-ground evidenced atoms with
  undetermined contribution are admitted with the aromatic rows. `admit` table tests
  migrate. Breaking. [dep: A0a]
- A0d: caller migration and the stage pin. Un-ignore
  `test_resolver_resolve_stages::case_2_aromaticity_bond_marks`; sweep lib and pytest
  expectations for inputs that previously acquired aromatic candidates without marks.
  Stage green: full lib suite, pytest. [dep: A0b, A0c]

### A1 — component factorization (additive)

- A1a: candidate-ring components (`umol-graph/src/ops/resolve/aromaticity.rs`): from the
  perception ring set, keep rings whose members are all aromatic-capable, connect over
  shared atoms, emit disjoint atom-set components. Unit tests: fused pair = one component,
  biphenyl-coupled rings = two, no candidate rings = none. Additive. [dep: A0d]
- A1b: per-component enumeration: the odometer runs per component;
  `MAX_JOINT_ASSIGNMENTS` → `MAX_ASSIGNMENTS`, now per component. Behavior-preserving for
  single-component inputs — the existing suite is the preservation test; new test: two
  independent flexible rings whose assignment product exceeds the bound while each
  component is under it. Additive. [dep: A1a]

### A2 — assignment-level selection (breaking; fixes the panic)

- A2a: the perception decomposition contract: pin that `find_systems` returns one
  decomposition per input molecule (deterministic, disjoint systems) with a test per rule;
  the selection relies on it. Additive. [dep: none]
- A2b: assignment enumeration replacing the `per_system` acceptance
  (`umol-graph/src/ops/resolve/aromaticity.rs`): per assignment, collect the member
  restriction and the partition; validity = stored-system consistency plus carrier
  elimination under `Error`; the value tie-break compares member-wise (sets equal under
  `Error` by construction); the canonical partition pick takes restriction-identical
  survivors; `tie_breaks` records value-key uses. Interim until A3: under a non-`Error`
  policy, structurally distinct survivors stay plural. Selection can no longer emit
  overlapping systems — asserted by test on the quinoline family. Breaking (the panicking
  family changes from panic to resolution). [dep: A1b, A2a]
- A2c: caller migration and acceptance. Quinoline and isoquinoline resolve through
  `ingest_smiles` with full-molecule expectations; bare-`n` imidazole under `Strict` is
  `Underdetermined` with both nitrogen splits in the report, and the existing
  `MostSaturated` imidazole case stays green; expectation deltas swept. Stage green: full
  lib suite, pytest. [dep: A2b]

### A3 — the aromaticity tie-break envelope (breaking)

- A3a: `AromaticityTieBreak` (`Strict` default, `MaxCoverage`) and the
  `AromaticityModel.tie_break` field (`umol-graph/src/ops/model.rs`), constructors and
  presets inheriting `Strict`; every `AromaticityModel` struct-literal site across the
  workspace migrates. Model table tests. Breaking. [dep: A2c]
- A3b: the structural order in selection under non-`Error` policies: claimed-atom count
  descending, then member sets lexicographic; `tie_breaks` gains the chosen systems'
  members when the structural order decided. Keep-mode tests: the tolerated-carrier family
  under `Strict` (plural) versus `MaxCoverage` (picked, recorded). Additive given A3a.
  [dep: A3a]
- A3c: umol-py: `AromaticityTieBreak` class, the third `AromaticityModel` field, package
  exports and inventory, pytest coverage mirroring the rust model tests. Breaking.
  [dep: A3a]

### A4 — pruned search (additive, deferrable)

- A4a: depth-first assignment search with the two sound prunes (completed-ring 4n+2 misses
  when no valid completion remains; interval bounds on unassigned members) replacing the
  flat per-component enumeration. All suites green unchanged. [dep: A2c]
- A4b: the enumeration-equivalence property: on components within the bound, the pruned
  search returns exactly the flat enumeration's valid-assignment set; placed per
  `docs/development/property-tests.md` (a new umol-graph property target if none exists).
  Additive. [dep: A4a]

### A5 — the local instrument (additive, never committed data)

- A5a: fetch script for the VEHICLe CSV into a gitignored working directory, with the
  provenance note; the `.gitignore` entry lands with it. Additive. [dep: none]
- A5b: the instrument run: an `#[ignore]`d test that skips when the local CSV is absent,
  classifies every row by the outcome pair under both valence tie-breaks, asserts 4n+2 and
  no-overlap on every emitted system, checks tautomer-cluster consistency, and writes the
  class distribution and failure manifest for review. The distribution is the A5
  deliverable, recorded in this document. Additive. [dep: A2c, A5a]

### A6 — promotions (additive; gates 194 S5d)

- A6a: curated promotions into the doc 194 S5c2 resolution conformance suite: the panic
  representatives, one full tautomer cluster, purine and pteridine, and the hand-written
  cases VEHICLe cannot supply (acenes, phenanthrene, coronene, azulene; porphine,
  phthalocyanine; biphenyl, carbazole, dibenzofuran, acridine). Additive; snapshots
  regenerate at S5d. [dep: A5b]
- A6b: the corpus loaders as acceptance: the substructure and fingerprint bench fixtures
  ingest the OpenSMILES corpus panic-free (bench `--test` mode); the instrument is deleted
  unless A5's distribution argued for a shipped bulk suite, and that decision is recorded
  here. Stage green; 194 S5d unblocks. [dep: A6a]

The critical path is A0 → A1 → A2; A2c is the defect fix's acceptance. A3 (Keep-mode
policy) and A4 (scale) are independent of each other and deferrable behind the critical
path; A5–A6 gate 194 S5d and need only A2. The doc 194 S5d dependency stays on A6.
