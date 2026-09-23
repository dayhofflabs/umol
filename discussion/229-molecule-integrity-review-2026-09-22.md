# 229 — Molecule integrity review

Status: Proposed
Date: 2026-09-22
Relates: [213](213-editor-overlay-storage-2026-08-27.md),
[215](215-integrity-minimization-2026-08-28.md),
[review guide](../docs/development/code-reviews.md),
[integrity guide](../docs/development/integrity.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Review result

The review and independent refutation are complete. **The molecule gate's accepted
domain matches doc 215; its execution still repeats work.** Fourteen findings
survived: one downstream consumer panic, ten local efficiency findings, two incorrect
documentation claims, and one focused regression gap. None calls for tightening or
weakening integrity, or introducing a checking framework. Implementation has not
started; this record contains no staged implementation plan.

Reviewed commit: `1fd0c3aa7fbda4e47ec122d0222dd5829ab772a9`.
Three review agents covered contracts/efficiency, names/structure/documentation,
and tests/generators; a fourth independently challenged their evidence and
normative premises. All worked in detached worktrees at that commit. No tracked
source changed. The user's relocation of the unchanged reference-check block in
the primary molecule.rs was preserved and was not reported as a defect.

Scope: Molecule construction, reference checks, the authoritative integrity gate,
editor publication, checked mutation, and direct shared-check consumers in
Reaction and ReactionSpan. The admission audit also traced downstream identity,
frame, symmetry, and remapping consumers. This is not a full audit of those
subsystems or Python. Doc 213's proposed editor lifecycle is a future consumer,
not missing implemented behavior.

The reviewers consulted the living guides, applicable skills, the status index,
docs 211/214/215 including 215's final corrections, and the whitepaper's design
intent. Each defense checked whether the governing documents deliberately deferred
or required the observed behavior. No such deferral defeats the findings below.

## Semantic admission audit

The governing discussion is **215 — Integrity closure and minimization**, including
its final 2026-08-29 corrections. The living integrity guide's admission test requires
a concrete downstream failure or repeated prerequisite; convenience at construction
is insufficient. The audit traced every current predicate to a consumer, rather than
treating the guide's inventory as proof of its own necessity.

**The accepted/rejected domain matches that contract. The execution is not yet
minimal.** No chemistry, satisfiability, normalization, or operation-specific rule
was found in the molecule gate. Repeated enforcement is identified separately below;
removing it need not admit any additional molecules.

### What each condition buys downstream

All source locations refer to the reviewed commit. Paths are under
umol-graph-ir/src/ir unless qualified. Each row states a necessary representation
condition, not a claim that every current execution of its check is necessary.

| Condition / error | Concrete downstream dependency |
| --- | --- |
| InvalidReference | Raw bond endpoints must be checked before graph-core Graph::build_csr indexes them (graph.rs:573–579). Relation/site/ligand IDs feed entity views and orbit arrays; nested constraint IDs feed Constraint::map and dense remapping (constraint/molecule.rs:321–368, molecule/remap.rs:120–126). Virtual ligand anchors also need valid IDs. |
| ElectronCountLengthMismatch | Aromatic and multicenter incidence reads counts by participant position (incidence.rs:270–292). Short vectors can panic; longer vectors have unmatched values and cannot reframe with the participant permutation (electrons.rs:50–54). |
| DuplicateParticipant | Actual-atom frames need a unique positional action, and simple incidence matching must distinguish occurrences. This covers covalent/noncovalent self-loops, repeated dative donors, donor/acceptor overlap, and repeated aromatic/multicenter participants (incidence.rs:250–300; molecule/pushout.rs:195–247). Stereo actual-ligand uniqueness is already implied by complete-ligand uniqueness; atom-site exclusion remains a separate condition. |
| BondsParallel | BondViews::of_id promises a singular bond for its unordered endpoint pair (view/bond.rs:74–79). Parallel rows make that lookup arbitrary. |
| DativeBondsIdentical | DativeBonds::coincident_id identifies one entity by the complete acceptor plus donor multiset (dative.rs:184–190). Distinct complete keys may share participants. |
| NoncovalentBondsParallel | NoncovalentBonds::coincident_id identifies one entity by endpoints, regardless of attribute kind (noncovalent.rs:154–163). |
| AromaticSystemsOverlap | AtomView::aromatic_system_id selects one incident owner, and aromatic_valence uses that owner (view/atom.rs:279–283,381–393). Overlap would require different membership semantics. |
| MulticenterBondsIdentical | MulticenterBonds::coincident_id uses the complete participant set as a singular identity (multicenter.rs:169–176). |
| StereoAtomSitesDuplicate | StereoAtomViews::at_id and AtomView::stereo_atom promise one stereo entity at an atom site (view/stereo.rs:146–149; view/atom.rs:320–322). |
| StereoBondSitesDuplicate | StereoBondViews::at_id and BondView::stereo_bond promise one stereo entity at a bond site (view/stereo.rs:472–474; view/bond.rs:183–185). |
| DuplicateStereoLigand | Reframing and pushout derive one permutation between frames (view/stereo.rs:288,612; molecule/pushout.rs:268,290). Repeated complete (anchor, kind) values do not determine that action uniquely. Different virtual kinds at one anchor remain distinct. |
| StereoFrameDegreeTooLarge | Kindless frames also undergo representative-action selection; Permutation::between asserts the six-position bound (stereo.rs:228–229; umol-perm/src/permutation.rs:146). |
| StereoLigandIncidenceMismatch | Atom ligands denote neighbors or virtual values anchored at their site: symmetry keeps virtual positions with the fixed site (symmetry.rs:447–467), and depiction interprets actual ligands as incident wedge bonds (umol-io/src/depict/molecule.rs:248–256). Bond frames have exactly two consecutive endpoint pairs, even without a kind; stereo action preserves or exchanges these blocks (stereo.rs:187–188; umol-perm/src/coset.rs:187–205), and layout relies on their endpoint association (umol-io/src/layout/coordgen.rs:72–82). |
| StereoKindSiteMismatch | Atom and bond sites select different frame-action groups (stereo.rs:181–188,957). Arity alone cannot distinguish degree-four geometries. This concerns frame interpretation, not chemical stereogenicity. |
| StereoLigandArity | Stored configurations and top-level stereo constraint wrappers must each fit the frame they address. Canonicalization applies kind-degree permutations to that frame (canonicalize.rs:2752–2759,3602–3603). A wrapper's kind need not equal the stored kind. |
| StereoCosetOutOfRange | Every literal, set member, and nested variable-domain member must name a coset of its kind. Otherwise canon_coset can mistake an invalid same-cardinality domain for the full universe (stereo.rs:1812–1815), or group action reports an unrelated failure. |
| StereoPermutationDegree | Inline symmetry/fluxionality actions and expression Apply terms compose with the frame or kind action (stereo.rs:811–817,1786–1788). Permutation::compose asserts equal degree. |
| StereoLigandPositionOutOfRange | Topicity transport applies the inverse frame action to positions; evaluation indexes position orbits (stereo.rs:831–837; symmetry.rs:419–425). An absent position has no interpretation. |

Doc 215's final correction matters: its earlier donor/acceptor-overlap permission
was withdrawn because it breaks incidence matching. Two bounded public graph-core
probes confirmed the dependency: a self-loop query falsely matches a loop-free
triangle three ways; a two-node donor/acceptor incidence graph with parallel,
differently labelled edges loses even its identity match (zero instead of one).
These counterfactual graphs are prevented by molecule integrity, not examples of
malformed molecules that currently pass publication.

The gate still admits kindless stereo, non-ground forms, and interpretable but
contradictory assertions, including empty coset sets. It does not check electron
totals/signs, valence, aromaticity, physical stereogenicity, model conformance, or
constraint truth. A same-degree expression action outside a kind's action group
remains an operation-level semantic failure, as settled in 214:196–204; it does not
justify an additional integrity check. Removed EntityCountMismatch and
incidence-wide DativeBondsParallel checks have not returned.

## Confirmed findings

Source locations below refer to the pinned commit. The main files are
[molecule.rs](../umol-graph-ir/src/ir/molecule.rs),
[integrity.rs](../umol-graph-ir/src/ir/molecule/integrity.rs), and
[editor.rs](../umol-graph-ir/src/ir/molecule/editor.rs).
Local efficiency findings use the review guide's **Priority order** (lines 16–24),
which explicitly admits suboptimal implementation under a sound design. Neither
allocations nor dynamic dispatch are forbidden in general.

### F1 — Unnecessary dynamic dispatch in reference checks

**Confirmed; initial lead L1.** The reference helpers in molecule.rs:148–307 and
integrity.rs:378–397 take &dyn Fn(Entity) -> bool. Their three callers supply
concrete closures over entry, molecule, or span sizes. Optimized production LLVM
retains indirect predicate calls. A generic borrowed Fn supports the same recursion.

**Defense:** one compiled traversal limits code size; some calls may be optimized
away, and benefit depends on constraint prevalence. The borrowed trait object does
not itself allocate. **Proposal:** use static dispatch through the existing
helpers. The isolated comparison below supports this choice, not a molecule-wide
speedup claim. No trait, public checker, or runtime dispatch object is needed.

### F2 — Dative identity keys are cloned on success

**Confirmed; L2.** integrity.rs:137–143 builds a donor BTreeSet and clones it into
the identity collection, retaining the original only for a duplicate diagnostic.

**Defense:** the clone keeps sorted error data immediately available.
**Proposal:** move the key and reconstruct the sorted donor diagnostic only on
failure. Preserve complete donor-plus-acceptor uniqueness and the selected error.
The final correction in 215 requires that uniqueness, not the clone. Empty sets
need not allocate; the issue is avoidable copying for nonempty donor sets.

### F3 — Stereo frames are copied for read-only checks

**Confirmed; L3.** integrity.rs:196/220 calls ligand_frame, whose implementation in
view/stereo.rs:717–721 collects a Vec. StereoAtoms::ligands and StereoBonds::ligands
already return the identical stored slice (stereo.rs:76/304).

**Defense:** the view method is convenient and valid frames are small; malformed
frames can be larger. **Proposal:** borrow the existing slices, preserving frame
order and all checks. No new accessor or public surface is necessary. Nonempty
frames avoid a temporary allocation; the gate never needs ownership of the copy.

### F4 — Multicenter diagnostics allocate before an error exists

**Confirmed; L4.** integrity.rs:168–172 collects participants into a Vec and then
a BTreeSet; the Vec is otherwise used only by the duplicate error.

**Defense:** unlike the identity set, the Vec preserves the supplied participant
order. **Proposal:** construct the identity set directly, then collect the original
sequence only on failure. Preserve diagnostic order; do not replace it with the
sorted identity key or claim that all identity-storage allocation disappears.

### F5 — Stereo-bond uniqueness is checked twice

**Confirmed.** check_stereo_bond_entry (integrity.rs:348–360) calls
check_stereo_frame, then separately hashes actual-atom IDs. Complete StereoLigand
equality includes both atom_id and kind. Two Atom-kind ligands with the same ID
would already have failed the first check with DuplicateStereoLigand.

**Defense:** the shared participant helper makes the actual-atom invariant visible.
It is correct, but its failure branch is unreachable after frame success.
**Proposal:** remove this second pass, retaining configuration and inline-constraint
checks. This is class-3 revalidation under the review guide (lines 44–56); the
same proof applies when Reaction uses the local stereo checker.

### F6 — Stereo-atom checking needs site exclusion, not another set

**Confirmed.** After check_stereo_frame, check_stereo_atom_entry
(integrity.rs:329–345) only needs to reject an actual ligand equal to the site.
Actual-ligand uniqueness has already been established.

**Defense:** the existing shared helper correctly combines both conditions.
**Proposal:** scan for an Atom-kind ligand equal to the site, returning the existing
DuplicateParticipant error. Allow virtual ligands anchored there. Keep frame errors
first. Basis: the review guide's class-3 revalidation and local efficiency criteria;
the exact-error comparison below supplements the proof.

### F7 — Bounded stereo uniqueness does not need a HashSet

**Confirmed; L5 narrowed to stereo.** check_stereo_frame (integrity.rs:308–327)
rejects lengths above MAX_DEGREE, currently six, before checking uniqueness. Testing
each ligand against its preceding slice needs at most 15 equality comparisons,
returns the same first duplicate, and uses no temporary collection.

**Defense:** hashing scales well for larger domains. That advantage does not apply
under this hard bound. **Proposal:** use the prefix scan here and retain the degree
guard first. Do not extend this conclusion to unbounded relation participants.
The bound and complete-value uniqueness are specified in data-types.md:894–903;
benchmark and exact-error evidence are below.

### F8 — Dative insertion documentation promises nonexistent sorting

**Confirmed; documentation.** editor.rs:666–668 says add_dative_bond sorts donors
through Unordered canonicalization. Insertion, storage conversion, and publication
preserve the supplied order. A public-API reproduction confirms donors [2,0] remain
[2,0] after try_build in a three-atom molecule with acceptor 1.

**Defense:** the comment predates the storage migration, and current dative payloads
are frame-invariant. It is nevertheless an individually false promise, not an
absent new documentation convention. **Proposal:** correct the comment, not storage.
Basis: code-reviews.md:77–94 and the frame-preservation contract in
data-types.md:894–897, retained by completed docs 211/214.

### F9 — Recursive stereo-expression failures lack retained boundary cases

**Confirmed; test coverage only.** check_term (integrity.rs:609–633) checks variable
domains, nested terms, and expression permutation degrees. Existing malformed
molecule cases cover literal cosets and constraint permutations, not malformed
expressions through try_from_entries. The scratch cases below reject correctly.

**Defense:** expression parsing/algebra and final error variants have other tests;
serialization generators intentionally omit operator terms. Those do not exercise
this independent integrity traversal. **Proposal:** add a compact public-constructor
table for malformed terms and an accepted raw term. No generator framework or
all-publishers test matrix is needed. Basis: integrity.md:101–102 and 174–178,
and the data-type-contracts skill's boundary-test requirements. Docs 214/215 do
not defer this contract.

### F10 — Assembled covalent endpoints are checked again

**Confirmed.** integrity.rs:123 repeats endpoint bounds checking after Graph
construction has already established it. Graph::new and Graph::add_edge build CSR before
returning or replacing storage; restoration checks endpoint bounds. Molecule's
graph/atom-table parallelism is a producer guarantee, the same reason 215 removed
EntityCountMismatch.

**Defense:** a complete gate might catch an internal publisher bug. Doc 215:121–129
explicitly rejects that as the runtime boundary for structural parallelism.
**Proposal:** remove only this assembled-graph endpoint recheck. Keep raw-entry
endpoint preflight: otherwise CSR construction can panic before returning the typed
InvalidReference error. Overlay and constraint references remain checked.

### F11 — Closed-frame consumers repeat kinded arity checks

**Confirmed.** symmetry.rs:128–130,190,228 checks kind-degree agreement after
stereo_center reads the stored frame/configuration from a published Molecule.
canonicalize.rs:3533 repeats it too: its sole caller remaps stored ligand IDs
without changing length, then immediately expects success because integrity
established arity (3602–3603).

**Defense:** these guards support malformed-frame fallbacks in otherwise reusable
helpers. Both paths are private and receive only closed-source frames; no such
fallback is reachable. **Proposal:** remove the redundant arity branches. Preserve
legitimate undetermined/nonliteral handling and action/normalization failures.
Basis: integrity.md's closed-container enforcement rule and 215:96–98.

### F12 — Symmetry repeats complete-ligand uniqueness

**Confirmed.** grade_center and has_oriented_center (symmetry.rs:190,228) call
all_distinct, allocating a HashSet for stored frames whose complete-value uniqueness
is already guaranteed. This is separate from meaningful equality of ligand classes
under molecular symmetry.

**Defense:** rejecting malformed frames avoids ambiguous transport. Publication
already rejects precisely those frames, and stereo_center only copies them.
**Proposal:** remove these two defenses, retaining generator compatibility and
ligand-class comparisons. Same closed-container basis as F11; no deferral in 214/215.

### F13 — Coset normalization documentation assigns range checks to tier 2

**Confirmed; documentation.** canon_coset's comment (stereo.rs:1795–1798) says
coset range checking belongs to the validator. It is a tier-1 aggregate condition
(integrity.md, StereoCosetOutOfRange), enforced recursively by the current gate.

**Defense:** open standalone forms need not run that aggregate check during
normalization. That is true, but does not make it tier 2. **Proposal:** correct the
classification; do not add checking to normalization or change accepted forms.

### F14 — Graph symmetry panics on an admitted kindless stereo frame

**Confirmed; consumer correctness.** graph_symmetry reaches stereo_center
(symmetry.rs:236–258), which calls view.kind(). That accessor expects a determined
kind (view/stereo.rs:187–191). A five-atom star with four distinct actual ligands
and StereoAtomForm::default() passes try_from_entries, then graph_symmetry panics
with "stereo view has a concrete kind". The public reproduction used
ConstitutionColoring::entity_only, Nauty, and one refinement iteration.

**Defense:** symmetry needs a geometry to grade orientation. However, its existing
nonliteral handling already permits no orientation contribution; graph_symmetry has
no declared determined-kind precondition. Kindless frames are expressly admitted
by 214 and the integrity guide. **Proposal:** handle absent kind in this consumer,
as it handles unresolved cosets. Do not reject these molecules at publication.
Add a focused public regression. This is a consumer contract defect, not a missing
integrity check or a complete audit of stereo-view getters.

## Nomenclature and visibility migrations

These are migrations, not additional correctness findings. The nomenclature guide
(lines 3–7) and review guide (90–94) distinguish later conventions from historical
defects. The validate_* reference helper names existed at cd64b530f (2026-08-05);
the Integrity check glossary entry appears at 1257445b0 (2026-08-08).

- Rename the five molecule validate_* reference helpers to check_*: they check
  representation references, not constraint satisfaction. Apply the same terminology
  to the directly related span entry-reference helper when touched. Name length alone
  is not a finding; an exhaustive constraint match remains appropriate.
- Without relocating code, make these four helpers private: molecule's
  validate_entry_references; integrity's check_stereo_frame, require_reference,
  and require_references. Workspace-wide searches found no external consumers.
  Shared entry/kind stereo helpers genuinely have Reaction consumers; the shared
  constraint-reference traversal genuinely has a ReactionSpan consumer.

Keep the existing module boundaries for the bounded optimization. The 686-line
integrity module is coherent; neither a new checker object nor a forced split
improves the identified problems. Consolidating reference helpers can be considered
separately if it improves locality, with visibility justified by actual callers.

## Open alternatives and rejected claims

| Question | Verdict and reason |
| --- | --- |
| Replace unbounded/global HashSets or BTreeSets (rest of L5) | **Open.** Sorted vectors and dense membership arrays have different sorting, initialization, and sparse-domain costs. The six-ligand experiment does not select these algorithms. Retain current choices pending a relevant paired comparison. |
| Fuse the two constraint traversals (L6) | **Open.** Reference checking currently precedes stereo-domain checking for the whole constraint tree (integrity.rs:246–250,442–469). Plain leaf-by-leaf fusion changes error priority and must still protect indexing. Saving domain errors while checking remaining references adds complexity. No equivalent simplification has been established. |
| Narrow constructor reference preflight (L7) | **Open.** Graph::build_csr indexes bond endpoints (graph.rs:576–579), so deleting preflight can panic. Relation incidence construction stores/sorts opaque IDs, so not every early relation check is needed for storage safety. Narrowing the pass still changes cross-entity error priority. Preserve the full editor publication gate; do not add validity flags or witness types. |
| Remove the result gate from trusted remapping | **Open policy clarification.** molecule/remap.rs:129 calls try_from_arcs, which invokes the whole gate. Doc 215:96–98 says trusted transformations do not recheck outputs, but data-types.md:148–170 retains checked/asserted publication through the same gate, and 215 S2a explicitly removed only the source check. The retained output check is factual; calling it a settled defect would ignore that conflict. No unchecked constructor is proposed here. |
| Full candidate checks are inherently wrong | **Refuted.** Completed 215:198–221 and data-types.md:153–170 expressly settle candidate-and-authoritative-gate mutation. Incremental integrity requires a separate preservation argument. |
| The integrity module needs a new checker framework or forced split | **Refuted.** Its size and shared consumers do not establish a structural defect under code-reviews.md:63–68. Local functions and existing storage access suffice. |

## Evidence and limits

All compilation and sampling ran in isolated worktrees on rustc 1.96.0,
LLVM 22.1.2, aarch64-apple-darwin. No Python, workspace-wide, or MSRV gate ran.
These probes inform local implementation choices; they do not measure total
integrity time, total allocations, code coverage percentage, or real workload frequency.

### Callback dispatch

The unchanged library was compiled with:

```sh
cargo rustc --offline -p umol-graph-ir --lib --release \
  --target-dir /Users/dr/.cargo-target/229-integrity-asm -- --emit=llvm-ir
```

Production LLVM retains indirect calls through %contains.1.val inside
validate_constraint_references. A separate scratch executable copied the exact
reference-helper bodies into dynamic and generic variants; the only body/signature
change was &dyn Fn to &impl Fn. It used the original public Constraint/Entity types
and a runtime eight-kind bounds predicate, each bound 4096.

Each tree is And of the stated number of DativeBondDonors leaves, with four donor
IDs and one bond ID each. A late-error variant appends an invalid AtomId(4096)
through DativeBondAcceptor. Both variants returned exactly the same results.
Construction was outside timing; black_box covered inputs/output. Medians of
11 alternating samples, 20,000 iterations each, nanoseconds per tree:

| Leaves | Success: dynamic / generic | Late error: dynamic / generic |
| ---: | ---: | ---: |
| 1 | 14.4 / 6.8 | 24.6 / 15.0 |
| 32 | 421.7 / 151.7 | 425.1 / 157.4 |
| 256 | 3318.4 / 1142.7 | 3263.7 / 1136.8 |

This supports static dispatch for this traversal. It is not a measured replacement
of the full molecule gate and does not establish a whole-operation speedup.

### Stereo uniqueness

A std-only model compared the current HashSet checks against ordered-prefix
complete-ligand equality plus actual-ligand/site exclusion. It enumerated every
frame of length 0..6 over three atom IDs and three ligand kinds, for a bond and
three atom-site choices: **2,391,484 exact Result agreements**. The unchanged
maximum-degree guard remains before this algorithm.

Timing compares only complete-frame uniqueness, not the combined F5/F6 changes.
An optimized single-run comparison used 500,000 iterations per case and black_box.
Distinct frames used successive atom IDs with Atom kind; early/late duplicates
repeated the first ligand at the second/final position. Nanoseconds per frame:

| Length | Distinct: hash / scan | Early duplicate: hash / scan | Late duplicate: hash / scan |
| ---: | ---: | ---: | ---: |
| 2 | 65.0 / 2.6 | 47.2 / 2.0 | 48.6 / 2.0 |
| 4 | 199.6 / 7.4 | 48.8 / 1.9 | 196.9 / 5.0 |
| 6 | 261.3 / 12.8 | 49.1 / 2.0 | 242.4 / 10.6 |

This supports the bounded scan and removal of the implied checks. It neither
selects an unbounded relation algorithm nor substitutes for production regressions.

### Existing tests and additional boundary probes

Focused commands used --offline and
CARGO_TARGET_DIR=/Users/dr/.cargo-target/229-integrity-tests:

| Command after cargo test -p umol-graph-ir | Result |
| --- | --- |
| --lib test_molecule_try_from_entries | 76 passed. |
| --lib test_molecule_try_modify | 15 passed. |
| --features proptest --test property molecule::publication | 5 passed, PROPTEST_CASES=256 each. |
| --features proptest --test property test_molecule_try_from_entries_rejects_stereo_coset_out_of_range | 1 passed, PROPTEST_CASES=256. |

Existing constructor cases cover every MoleculeIntegrityError variant, and checked
mutation cases verify exact errors and unchanged sources. Publication properties
cover remap, extract, combine, split, and successful apply. This does not imply
complete branch coverage.

The actual molecule_entries_with_constraints_strategy was sampled 2,048 times
using Proptest 1.9.0's deterministic runner, from an unchanged copy of the strategy
module in a scratch crate using public production APIs. All samples published.
Atom counts 0..5 occurred [370,333,315,358,342,330] times; presence of each entity
kind in AB/DAMN/SS order was [1678,876,599,530,504,818,62,494]. Stored stereo kinds
were CisTrans 633, Octahedral 1, SquarePlanar 29, Tetrahedral 34, and
TrigonalBipyramidal 4. Cosets were Undetermined 332, Lit 369, LitSet 0, and Term 0.
The generator deliberately chooses only Undetermined/Lit here; this is supporting
evidence for F9, not a separate generator defect or a model of real chemistry.

Six public-constructor probes used a five-atom star with tetrahedral site 0 and
ligands [1,2,3,4]. LitSet({0,2}), Term(Var with domain {0,2}), Term(Lit(2)),
Term(LitSet({0,2})), and Term(Swap(Mirror(Lit(2)))) returned
StereoCosetOutOfRange for StereoAtom(0), Tetrahedral, coset 2, count 2.
Term(Apply(Var, identity permutation of degree 3)) returned StereoPermutationDegree,
expected 4, actual 3. The six exact assertions passed without a production change.

## Proposed design and next action

Fix the confirmed kindless-symmetry panic in the consumer, preserving the admitted
domain. The bounded optimization then preserves accepted inputs, existing
diagnostics, and publication boundaries: generic reference predicates, borrowed
stereo frames, error-only diagnostic collection, one bounded uniqueness scan,
separate site exclusion, and removal of already-implied endpoint/closed-frame
checks. Include the precise naming/visibility migrations, the two documentation
corrections, and focused expression and symmetry regressions. No new public symbols
are needed.

Keep the two-pass constraint walk, constructor preflight, remapping result gate,
and unbounded set algorithms unchanged in that work. Their alternatives remain
separately identified above; they do not block the local improvements. A later
implementation should measure the combined effect on the actual gate without conflating it with storage
construction or editor rebuilding. These decisions and evidence are durable here;
review scratch files are disposable.
