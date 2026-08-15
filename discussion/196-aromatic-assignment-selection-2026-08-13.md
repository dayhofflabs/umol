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

1. **Structural order (`AromaticityTieBreak`, new policy).** One
   criterion (corrected 2026-08-13: under closed-world admission, claimed atoms and
   unrealized carriers partition the evidenced set, so "full fit first" and "coverage" are
   complements — the earlier two-component phrasing was redundant; extended 2026-08-15
   with the electron component on the pteridine evidence): claimed-atom count
   descending, then the systems' electron total ascending, then the member sets
   themselves (lexicographic) — so the value order below
   never compares assignments restricting different atoms. The electron component acts
   whenever the rule admits more than one total for the same members (bare-`n` pteridine:
   all-pyridinic 10 and all-pyrrolic 14 both satisfy 4n+2, the mixed sums fail; minimal
   synthesis picks 10); the coverage component is reachable only under a non-`Error`
   failure policy, where partial-fit assignments survive validity. Ring size and system
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

**Component factorization.** The aromaticity rule's acceptance — whether a set of atoms is
perceived as aromatic — couples atoms only within a connected component of the
aromatic-candidate graph: build the perception ring set once, keep the
candidate rings (every member aromatic-capable), and take connected components over shared
atoms. Enumeration, validity, selection, and the bound are all per component; component
results combine freely (no joint constraint spans components — biphenyl's rings, or any
multi-fragment molecule). The assignment count drops from the product over the molecule to
the maximum over components.

**Search.** Within a component, replace the flat enumeration with a depth-first search over
flexible atoms ordered ring-by-ring, pruning a partial assignment when the violation is
certain:

- a candidate ring with all members assigned that the rule rejects, when the assignment
  cannot become valid without it (under `Error`: some fully-assigned all-aromatic atom has
  no other candidate ring left);
- interval bounds where the rule admits them: each unassigned member carries a minimum and
  maximum contribution over its remaining disjuncts, and a ring whose reachable
  contribution range the rule can never accept is settled early (the Hückel rule prunes by
  the 4n+2 residues; a rule without a usable bound simply prunes less — soundness, not
  completeness, is the requirement).

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
    /// Claimed-atom count descending, then electron total ascending, then
    /// member sets lexicographic: realize as much of the evidence as the
    /// rule allows while synthesizing as little as possible beyond the
    /// written structure. The perception decides what a system is; the
    /// policy only orders the valid assignments.
    MinElectronCount,
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

Naming and defaults, settled 2026-08-13; revised 2026-08-15 on the pteridine evidence
(see the A6a open item and its resolution below):

- the variant is `MinElectronCount` (2026-08-15, replacing `MaxAtomCount`, itself renamed
  from `MaxCoverage` 2026-08-14 — "coverage" duplicated the claimed-atom count under a
  second word, and `atom_count` is the established count vocabulary; `MostAromatic` was
  rejected earlier — "aromatic atom" reads as aromaticity being the atom's property).
  Coverage-only carried no meaning where the two orderings differ, so the axis keeps ONE
  non-`Strict` variant with the full ordering — claimed-atom count descending, electron
  total ascending, member sets lexicographic — named for its semantically active
  component, as `MostSaturated` names its leading key. The electron component is the
  active one under `Error` policies; coverage discriminates only in Keep-mode;
- the default is `Strict`, as for the valence tie-break: any other default would bias the
  scheme. The neutrality doctrine holds at the type level; conventions live in presets;
- the `daylight`/`mdl`/`permissive` presets **adopt `MinElectronCount`** (2026-08-15,
  reversing the 2026-08-13 decision). The earlier rationale — "no reference behaviour to
  reproduce, since reference toolkits pin hydrogens at parse time" — is contradicted by
  the data: ChEBI and ChemSpider store pteridine as bare-`n` `c1cnc2ncncc2n1` (PubChem
  holds the Kekulé form), and the reference reading of such input is minimal-H — RDKit
  gives bare `n` zero implicit hydrogens and either kekulizes on that basis or refuses
  (H0 or bust, never H1; bare purine cannot kekulize — nine π atoms admit no perfect
  matching — and is refused; ChemSpider accordingly writes purine with explicit `[nH]`).
  `MinElectronCount` reproduces that convention as the usable superset: add exactly the
  hydrogens the rule forces, never more — pteridine resolves all-pyridinic like every
  reference, bare purine (which references refuse) resolves to the unique one-`[nH]`
  assignment, the original doc 174 goal. SMILES and MOL inputs read alike (no
  interpretation difference when hydrogens are missing), and `permissive` adopts it too —
  its purpose is processing arbitrary real-world input;
- the conformance matrix (194 S5c2) still does not gain this axis as a cell dimension;
  the cells inherit the presets, so `MinElectronCount` lands in every cell at the S5d
  regeneration.

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
class: every emitted system re-validates under the model's aromaticity rule, and no two
systems overlap — the latter being the invariant this document's defect violates. Tautomer clusters are the second oracle: rows in one cluster
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
  ground-evidenced rows are new. Breaking. [dep: A0a] **Done 2026-08-13:** the evidence
  check reads stored-system membership from the derived side plus
  `asserted_complete(AromaticValence)` — extensionally equal on the positive side for
  counts (its per-atom boolean already made candidates uniform), so no existing row
  migrated; the semantic change is the ground-evidenced branch in `admitted_completions`:
  a field-ground atom whose closed-world evidence is `Aromatic(Undetermined)` and which is
  in no stored system enters the carrier with one candidate per table aromatic valence.
  The in-system guard is incidence (`is_in_aromatic_system`), not the derived reading —
  the derived bond-adjacency arm would otherwise exclude exactly the bond-marked atoms the
  branch serves. Only the bond-carried placement reaches the branch: an atom-carried `#a+`
  makes the constraint store non-ground and the existing path admits it. Admit tests:
  the six-atom bond-marked ring enters with `a1` candidates; contribution-asserted,
  in-system, and unmarked ground atoms stay complete. The A0d stage pin came early:
  `case_2_aromaticity_bond_marks` un-ignored and green — the MDL explicit-hydrogens shape
  resolves end-to-end on the counts model. Lib 910 green, pytest 1306, clippy zero.
- A0c: atom-typing admission on the same reading
  (`umol-graph/src/ops/valence/atom_typing.rs`). Unmarked atoms no longer admit aromatic
  registry rows (the mixed-candidate case dissolves); field-ground evidenced atoms with
  undetermined contribution are admitted with the aromatic rows. `admit` table tests
  migrate. Breaking. [dep: A0a] **Done 2026-08-13:** the row filter's open-world
  `(None, None) ⇒ compatible` arm — exactly where unmarked atoms admitted aromatic rows —
  now reads `asserted_complete(key)`; present evidence (asserted, derived, or both)
  evaluates as before, and ring keys stay asserted-only. The ground-evidenced branch
  admits only aromatic rows, field-compatibility filtering among them (a ground `N#h1`
  takes the pyrrolic `h1·a2` row; the pyridinic `h0·a1` row is field-incompatible).
  One test migrated as planned: the plural-admission input gains its `#a+` mark, since
  unmarked plurality was the dissolved case; new pins cover the dissolution (an unmarked
  atom against an aromatic-rows-only registry is `NoMatch`) and the ground-evidenced
  bond-marked pair. `classify_molecule_atom` is untouched — classification is the
  validator's conformance reading, not resolution's. Lib 912 green, pytest 1306, clippy
  zero. One recorded asymmetry: with aromatic evidence unrealizable by any state,
  atom-typing errors at admission while counts defers to the discharge `⊥` — each source's
  native surfacing; alignment, if wanted, is an A2-adjacent question.
- A0d: caller migration and the stage pin. Un-ignore
  `test_resolver_resolve_stages::case_2_aromaticity_bond_marks`; sweep lib and pytest
  expectations for inputs that previously acquired aromatic candidates without marks.
  Stage green: full lib suite, pytest. [dep: A0b, A0c] **Done 2026-08-13**, absorbed by
  the earlier subitems: the un-ignore landed with A0b (the counts fixture powers the
  case), and the expectation sweep reduced to the one planned migration at A0c (the
  plural-admission input gaining its `#a+` mark) — no other lib or pytest expectation
  depended on unmarked aromatic candidates. Stage closed green: lib 916, pytest 1306,
  clippy zero.

### A1 — component factorization (additive)

- A1a: candidate-ring components (`umol-graph/src/ops/resolve/aromaticity.rs`): from the
  perception ring set, keep rings whose members are all aromatic-capable, connect over
  shared atoms, emit disjoint atom-set components. Unit tests: fused pair = one component,
  biphenyl-coupled rings = two, no candidate rings = none. Additive. [dep: A0d]
  **Done 2026-08-13:** `candidate_components(rings, capable)` — iterative merge over
  shared atoms, deterministic sorted output, capability as a caller-supplied predicate so
  the A1b wiring decides the evidence composition. Four table cases: the fused pair one
  component, coupled rings two, a chain none, and a capability split leaving the capable
  ring alone. Carries a targeted dead-code allow until A1b wires it. The
  "odometer" identifiers and doc phrasing are purged (`assignment_indices`; "index-vector
  enumeration"). Lib 916 green, clippy zero.
- A1b: per-component enumeration: the index-vector enumeration runs per component;
  `MAX_JOINT_ASSIGNMENTS` → `MAX_ASSIGNMENTS`, now per component. Behavior-preserving for
  single-component inputs — the existing suite is the preservation test; new test: two
  independent flexible rings whose assignment product exceeds the bound while each
  component is under it. Additive. [dep: A1a] **Done 2026-08-13:** the perception exposes
  `candidate_rings` (the exact `find_systems` construction, now shared); `select` derives
  the components with capability = any carrier disjunct contributing, else the stored
  contribution — the fallback chain extracted as `stored_contribution`, shared by the
  enumeration closure. The bound checks per component before any enumeration; each
  component enumerates its own flexible atoms with other components held at their first
  disjunct, and accumulation keeps only systems inside the component, so cross-component
  duplicates never enter `per_system`. Flexible atoms outside every component are not
  enumerated — no candidate ring can claim them, and they fall to the finalization
  tie-break. Preservation: the full suite unchanged. Acceptance: five bridged triazines —
  fifteen flexible nitrogens, assignment product 2¹⁵ over the whole molecule, 2³ per
  component — flip from bound-underdetermined to `Determined` with five systems
  (`test_ingest_smiles_components`), confirmed end-to-end from python. Lib 917 green,
  pytest 1306, clippy zero (one displaced `#[allow(clippy::complexity)]` restored to
  `find_systems`).

### A2 — assignment-level selection (breaking; fixes the panic)

- A2a: the perception decomposition contract: pin that `find_systems` returns one
  decomposition per input molecule (deterministic, disjoint systems) with a test per rule;
  the selection relies on it. Additive. [dep: none] **Done 2026-08-13:** the contract
  holds as implemented and needed no production change — Hückel merges overlapping
  passing candidates to disjoint unions through union-find, Hmo and Clar are
  component- and independent-set-based, every rule sorts its member lists, and the shared
  first-atom sort in `find_systems` is therefore a total deterministic order (the
  `HashMap` inside the Hückel merge affects only pre-sort intermediate order). Pinned by
  a four-case table over the three rules: naphthalene's overlapping candidates merge to
  one ten-atom system under Hückel, coupled rings come out in ascending order, Hmo at
  threshold zero accepts the plain ring, Clar selects the first sextet; each case asserts
  exact member sets, pairwise disjointness, sorted lists, and run-to-run equality.
  Lib 921 green, clippy zero.
- A2b: assignment enumeration replacing the `per_system` acceptance
  (`umol-graph/src/ops/resolve/aromaticity.rs`): per assignment, collect the member
  restriction and the partition; validity = stored-system consistency plus carrier
  elimination under `Error`; the value tie-break compares member-wise (sets equal under
  `Error` by construction); the canonical partition pick takes restriction-identical
  survivors; `tie_breaks` records value-key uses. Interim until A3: under a non-`Error`
  policy, structurally distinct survivors stay plural. Selection can no longer emit
  overlapping systems — asserted by test on the quinoline family. Breaking (the panicking
  family changes from panic to resolution). [dep: A1b, A2a]
  **Done 2026-08-14:** `select` collects one assignment per index choice — the partition
  (perceived systems inside the component, sorted by member list) and the restriction
  (chosen forms of flexible member atoms) — deduplicated, then filtered by validity:
  every stored system touching the component must reappear with the same member set
  (elimination to empty names the first unreproduced stored system under
  `aromatic_system_failure: Error`; a non-`Error` policy leaves the component inert), and
  under `aromatic_valence_failure: Error` an assignment may not leave a component atom
  whose every disjunct requires aromaticity unclaimed. Selection: a unique surviving
  restriction accepts directly; several compare by the value key member-wise when their
  domains coincide, recording key-narrowed atoms in `tie_breaks`; under `Strict`, a key
  tie, or differing domains the members stay plural and nothing is accepted;
  restriction-identical survivors take the lexicographically smallest partition. `claimed`
  covers every valid survivor's members, so tolerated and tied carriers pass the
  unclaimed-carrier check unchanged. The sequential per-system narrowing and its
  `retain` are gone. Tests: the quinoline select case pins one 10-atom system with the
  pyridinic nitrogen narrowed and no tie-break record (the pyrrolic assignment dies on
  the azine carriers); the stored-conflict table pins `AromaticSystemFailure` under
  `Error` and the inert component under `Keep`. End-to-end: quinoline and isoquinoline
  `from_smiles` yield one fused system each (previously the overlap panic). Lib 924
  green, pytest 1306, clippy zero; the resolution conformance suite stays in its known
  S5c2 stale-snapshot state (653 key renames, regeneration at S5d).
- A2c: caller migration and acceptance. Quinoline and isoquinoline resolve through
  `ingest_smiles` with full-molecule expectations; bare-`n` imidazole under `Strict` is
  `Underdetermined` with both nitrogen splits in the report, and the existing
  `MostSaturated` imidazole case stays green; expectation deltas swept. Stage green: full
  lib suite, pytest. [dep: A2b]
  **Done 2026-08-14:** quinoline and isoquinoline pinned in
  `test_ingest_smiles_aromatic_nitrogen` with full-molecule expectations — one ten-atom
  system, electrons `[1;10]`, fusion carbons `#h0`, nitrogen `#h0#n`. New
  `test_ingest_smiles_with_tie_break`: bare-`n` imidazole under `Strict` (via
  `ValenceModel { tie_break: Strict, ..smiles() }`) errors `Underdetermined` with the
  report carrying exactly the two nitrogen entries (`#h0#n#v2#a` / `#h#n0#v2#a2` each) and
  no tie-break records; the `MostSaturated` imidazole case stays green beside it. The
  expectation sweep was empty — no other lib or pytest expectation moved (A2b left both
  suites green). Stage closed green: lib 927, pytest 1306, clippy zero. A2 complete; the
  panic family's acceptance is pinned end-to-end.

### A3 — the aromaticity tie-break envelope (breaking)

- A3a: `AromaticityTieBreak` (`Strict` default, `MaxAtomCount`) and the
  `AromaticityModel.tie_break` field (`umol-graph/src/ops/model.rs`), constructors and
  presets inheriting `Strict`; every `AromaticityModel` struct-literal site across the
  workspace migrates. Model table tests. Breaking. [dep: A2c]
  **Done 2026-08-14:** the enum and field landed per the model-surface spec; the three
  presets set `Strict` explicitly, matching the `ValenceModel` constructor style. Forty-five
  struct-literal sites migrated across umol-graph (ingest, parse, aromaticity, resolve,
  validate, model tests) and umol-py's rust side, `tie_break` last in declaration order.
  umol-py interim until A3c: the pyclass `AromaticityModel` does not carry the field yet;
  `to_rust` fills `Strict`. Model tests: presets pin `Strict` in their full-literal
  expectations, `test_aromaticity_model_eq_difference` gains the `tie_break` case,
  `test_aromaticity_tie_break_default` pins the default variant. Lib 929, umol-py 1620,
  pytest 1306, workspace clippy zero.
- A3b: the structural order in selection under non-`Error` policies: claimed-atom count
  descending, then member sets lexicographic; `tie_breaks` gains the chosen systems'
  members when the structural order decided. Keep-mode tests: the tolerated-carrier family
  under `Strict` (plural) versus `MaxAtomCount` (picked, recorded). Additive given A3a.
  [dep: A3a]
  **Done 2026-08-14:** `AromaticityResolver` carries the model's `tie_break`; under
  `MaxAtomCount` the structural stage runs on the valid survivors between validity and the
  value order — key `(claimed-atom count descending, member-set lists lexicographically
  smallest)` — and the winner's member atoms enter `tie_breaks` on acceptance when the
  stage actually decided (more than one structure present). Under the model's `Strict`
  the stage never runs, keeping the A2b behavior; under the `Error` failure policy it is
  inert by construction, so no gate is needed. Tests: the tolerated-carrier pair on the
  flexible-nitrogen five-ring under `Keep` — `MaxAtomCount` accepts the full ring, narrows
  the nitrogen, and records all five members; the model-`Strict` twin passes the state
  through unchanged. Lib 935 green, clippy zero.
- A3c: umol-py: `AromaticityTieBreak` class, the third `AromaticityModel` field, package
  exports and inventory, pytest coverage mirroring the rust model tests. Breaking.
  [dep: A3a]
  **Done 2026-08-14:** the pyclass enum mirrors `ValenceTieBreak`'s shape (repr,
  `from_rust`/`to_rust`); the py `AromaticityModel` carries the third field with the
  `Strict` constructor default matching `ValenceModel`'s, and `to_rust` maps the field
  instead of filling `Strict`. Registered in the module, exported in `__init__.py` and the
  `test_import` inventory. Rust-side py tests: the twelve model literals carry the field,
  the hmo rows pin the `MaxAtomCount` mapping in both directions and its repr, and the
  `ChemistryModel` repr expectations embed the new segment. pytest mirrors the valence
  tie-break coverage: equality, reprs, constructor default, the explicit-field
  constructor, the mutation row, and the embedded `ChemistryModel` repr. umol-py rust
  1620, pytest 1311, clippy zero. A3 stage complete.

### A4 — pruned search (additive, deferrable)

- A4a: depth-first assignment search with the two sound prunes (rule-rejected completed
  rings when no valid completion remains; rule-supplied interval bounds on unassigned
  members) replacing the flat per-component enumeration. All suites green unchanged. [dep: A2c]
  **Done 2026-08-15** (design settled in session): the two prune families collapse into
  one criterion — claim candidates settle via intervals, and a partial assignment is cut
  when some aromatic-only atom has every containing candidate settled-rejected. Claim
  candidates count fused unions (settled: rings-only would wrongly eliminate
  azulene-class assignments where every ring fails but a union passes); the perception
  exposes `claim_candidates` (rings plus, under Hückel, the fused-combination enumeration
  within the ring limits). The rule bound is `accepts_range(&self, members: &[(u32, u32)])`
  — per-member reachable contribution ranges, what the search holds — as a method on each
  rule struct beside `find_from_rings`: Hückel tests the summed interval for a 4n+2 value,
  Hmo and Clar accept every range (no usable bound, never settle). A member fixed or
  narrowed to no contribution settles its candidates directly. The DFS is inline iterative
  backtracking in `select` (graph-core traversal has no visitor surface and the decision
  tree is not a materialized graph); flexible atoms order ring-by-ring; leaves run
  `find_systems` unchanged; pruning is active only under `aromatic_valence_failure: Error`,
  the criterion's policy. One recorded semantics adjustment: the validity filter now runs
  carrier elimination before the stored-system naming loop — the search prunes by the
  carrier criterion, so stored-failure naming must act on the carrier-valid set to stay
  enumeration-independent; observable only when both eliminations empty a component,
  which now surfaces as the failure family (unclaimed carrier via the global check)
  rather than the stored family — no suite expectation moved. Tests: `accepts_range`
  tables per rule (seven Hückel rows incl. span and below-two edges), `claim_candidates`
  table (Hückel rings-plus-union vs Hmo/Clar rings-only). Lib 947 green unchanged,
  pytest 1311, clippy zero; quinoline, imidazole, and the bridged triazines confirmed
  end-to-end.
- A4b: the enumeration-equivalence property: on components within the bound, the pruned
  search returns exactly the flat enumeration's valid-assignment set; placed per
  `docs/development/property-tests.md` (a new umol-graph property target if none exists).
  Additive. [dep: A4a]
  **Done 2026-08-15:** new `umol-graph` property target (`tests/property.rs`, feature-gated
  `proptest` mirroring graph-ir). The executable form is search-independence at the public
  boundary: `select` compared by `==` against a definition-level flat-enumeration selection
  (`exhaustive_select` in `tests/property/resolve.rs`), sharing only perception
  (`find_systems`) and the value-key comparison (`compare_by_key`) with production.
  Operational domain (`strategies::select_scenario`): one- and two-ring (fused/coupled)
  Hückel skeletons, sizes five and six, every ring atom in the carrier with literal
  contributions from a five-entry completion pool, at most three flexible atoms, no stored
  systems, both failure policies and both tie-breaks on each axis. The semantic property is
  stated in `select`'s rustdoc (`# Semantic properties`). Domain notes: stored-system
  consistency and non-carrier stored contributions are outside the generated domain (pinned
  by the A2b unit tests); the assignment bound is unreachable within it. Green at 2000
  cases (3 s); default run 0.4 s. Lib 947, clippy zero. A4 stage complete.

### A5 — the local instrument (additive, never committed data)

- A5a: fetch script for the VEHICLe CSV into a gitignored working directory, with the
  provenance note; the `.gitignore` entry lands with it. Additive. [dep: none]
  **Done 2026-08-15:** `scripts/fetch-vehicle.sh` stages `VEHICLe.csv` and the FTP README
  into `materials/aromaticity/vehicle/` (staging directory settled in session; the whole
  `materials/` tree is already covered by the repository `.gitignore`, so no new entry was
  needed). The script header carries the provenance (Pitt et al., doi:10.1021/jm801513z;
  ChEMBL FTP source of record) and the never-commit/never-redistribute note. The live CSV
  has 24 867 rows and one column beyond the README's list (`Pgood`).
- A5b: the instrument run: an `#[ignore]`d test that skips when the local CSV is absent,
  classifies every row by the outcome pair under both valence tie-breaks, asserts
  rule-acceptance and no-overlap on every emitted system, checks tautomer-cluster consistency, and writes the
  class distribution and failure manifest for review. The distribution is the A5
  deliverable, recorded in this document. Additive. [dep: A2c, A5a]
  **Done 2026-08-15:** `umol-graph/tests/vehicle.rs` (`--release --test vehicle --
  --ignored`, ~3.6 s; skips with a fetch hint when the CSV is absent). Each row ingests
  under `ValenceModel { tie_break, ..smiles() }` for both tie-breaks; rule-acceptance is
  re-validated by running perception `derive` on every resolved molecule (each stored
  system reassessed from its electron contributions under the daylight model) and overlap
  checked directly; outputs land beside the corpus as `distribution.txt`/`manifest.txt`.
  **The distribution (24 867 rows):** determined 24 521 (98.61 %), tie-break picks 346
  (1.39 %); zero tie-survives, zero contradictory, zero parse failures, zero hierarchy
  violations, zero invariant violations. The 294-row panic family resolves inside those
  classes. Only 346 of the 4 221 bare-`n` rows need the value key at all — the ring sums
  force a unique valid assignment for the rest. Tautomer oracle over the 1 544 clusters:
  five class splits, all the same benign shape — an aromatic-written member (flexible,
  tie-break picks) clustered with a Kekulé-written member (no aromatic marks, determined),
  i.e. different input evidence, not an inconsistency; canonical-form counts over resolved
  `MostSaturated` members: 1 370 two-member clusters keep two distinct forms (distinct
  pinned tautomers), 11 converge to one form, larger clusters analogous
  ((3,3)=124, (4,4)=24, (5,5)=3, (6,6)=2, (3,1)=8, (3,2)=1, (4,1)=1). The health
  criterion is met: no contradictions anywhere in the corpus and both defect-class
  invariants hold on every emitted system. A5 complete.

### A6 — promotions (additive; gates 194 S5d)

- A6a: curated promotions into the doc 194 S5c2 resolution conformance suite: the panic
  representatives, one full tautomer cluster, purine and pteridine, and the hand-written
  cases VEHICLe cannot supply (acenes, phenanthrene, coronene, azulene; porphine,
  phthalocyanine; biphenyl, carbazole, dibenzofuran, acridine). Additive; snapshots
  regenerate at S5d. [dep: A5b]
  **Done 2026-08-15:** fifteen data files. Thirteen under `data/aromatic/`: quinoline,
  isoquinoline, benzimidazole (bare), purine, pteridine, anthracene, azulene, coronene,
  carbazole, dibenzofuran, acridine, porphine, phthalocyanine (biphenyl and phenanthrene
  already present); the full tautomer cluster 1390 under `data/vehicle/` as
  `s16534.edn`/`s19784.edn` (the aromatic-written/Kekulé-written pair, category named for
  provenance). Inputs are the verbatim SMILES raise (generated via
  `Smiles::parse → try_into_ir`, self-contained, no config overrides); the porphine,
  phthalocyanine, and pteridine aromatic transliterations were cross-checked against the
  PubChem kekulized references on atoms, bonds, element multisets, and degree sequences.
  All fifteen execute panic-free; snapshots regenerate at S5d as planned. Cell outcomes at
  promotion time: quinoline, isoquinoline, coronene, azulene, and the Kekulé cluster
  member fully determined; purine, benzimidazole, and the aromatic cluster member
  `Strict`-plural with `MostSaturated` picking; porphine contradictory (aromaticity
  inconsistency) and phthalocyanine contradictory (discharge) in all four cells —
  recorded as-is, the macrocycle families are new behavior to assess.
  **Surfaced by the pteridine promotion, settled 2026-08-15:** on fused polyaza systems
  where more than one total passes the Hückel rule, the `MostSaturated` value key
  preferred the maximally-hydrogenated assignment — bare pteridine resolved to all four
  nitrogens pyrrolic (electrons sum 14 = 4·3+2), where every reference toolkit reads
  bare `n` pyridinic (sum 10); ChEBI and ChemSpider store exactly this bare-`n` form.
  Resolution: the structural order gained the electron-total-ascending component and the
  axis's non-`Strict` variant became `MinElectronCount`, adopted by the presets — see the
  revised "Naming and defaults". `Strict` continues to report the plurality honestly.
  Implemented 2026-08-15: variant, ordering, and presets landed with the equivalence
  property's definition-level reference updated in step (2000 cases green); pinned by the
  six-atom tetrazine select case (totals 6 and 10 both pass, the electron component picks
  6 and records the members with no value key consulted) and the end-to-end pteridine
  ingest expectation (all-pyridinic, electrons `[1;10]`, every nitrogen `#h0`). The
  VEHICLe instrument re-run under the new presets: **determined 24 867 (100.00 %)** —
  every former picks-class row was total-count ambiguity, now settled structurally, so
  both tie-breaks agree on the entire corpus; zero cluster class splits remain, and the
  formerly converging clusters (11+8+1) now hold pairwise-distinct canonical forms — bare
  members resolve minimal instead of collapsing onto their `[nH]` siblings. Invariants
  zero throughout. Lib 949, umol-py 1620, pytest 1311, property 2000 cases, workspace
  clippy zero.
- A6b: the corpus loaders as acceptance: the substructure and fingerprint bench fixtures
  ingest the OpenSMILES corpus panic-free (bench `--test` mode); the instrument is deleted
  unless A5's distribution argued for a shipped bulk suite, and that decision is recorded
  here. Stage green; 194 S5d unblocks. [dep: A6a]

The critical path is A0 → A1 → A2; A2c is the defect fix's acceptance. A3 (Keep-mode
policy) and A4 (scale) are independent of each other and deferrable behind the critical
path; A5–A6 gate 194 S5d and need only A2. The doc 194 S5d dependency stays on A6.
