# 229 — Aggregate integrity review

Status: In Progress
Date: 2026-09-22
Relates: [213](213-editor-overlay-storage-2026-08-27.md),
[215](215-integrity-minimization-2026-08-28.md),
[review guide](../docs/development/code-reviews.md),
[integrity guide](../docs/development/integrity.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Review result

The Molecule gate review and independent refutation are complete. **The molecule
gate's accepted domain matches doc 215; its execution still repeats work.**
Sixteen Molecule findings survived: one downstream consumer panic, twelve local
efficiency findings, two incorrect documentation claims, and one focused
regression gap. None calls for tightening or
weakening integrity, or introducing a checking framework. S0a has established
the benchmark baseline; no integrity-check behavior has changed.

The separate Reaction and ReactionSpan review at `33ecfe5e03d714d2ac801015d20fa00243ea3a82`
confirmed five Reaction publication gaps under the then-current contract, one
Reaction integrity overreach, two incorrect ReactionSpan doc claims, and two
bounded execution costs. No new
ReactionSpan admission gap was found.
The [verdicts](#adversarial-verdicts) distinguish these from the earlier Molecule
cycle and from deliberately deferred reaction materializability.
A findings-blind Reaction rereview at the same commit independently confirmed
the publication gaps and resolved a challenge to the kind-change verdict; its
additional dispositions are recorded below.
A subsequent focused review found that removal-incidence agreement can be
treated as an expected-old precondition, but did not establish that moving its
check out of Reaction construction is worthwhile now. The existing constructor
check stays for this work. R1 still requires a constructor fix for removals of
newly added bonds.
A second findings-blind boundary review confirmed R2–R5 as local representation
gaps and found that bare `Deltas::normalize` can erase a same-id `Add` followed
by `Remove` without checking the removal. That cancellation must be fixed
independently of Reaction construction.

The [small-relation comparison](#small-relation-integrity-checks) covers 2–10
participants and IDs below 100. Two-word bitmaps remain faster when sorted keys
also use inline storage. For aromatic systems, replacing the global membership
HashSet provides most of the gain. This supports temporary bitmap checking for
small namespaces. The decision below uses that path for at most 128 atoms and
the unrestricted sorted-key path above that size. No permanent graph-core index
is required for these gains.

Molecule review commit: `1fd0c3aa7fbda4e47ec122d0222dd5829ab772a9`.
Three review agents covered contracts/efficiency, names/structure/documentation,
and tests/generators; a fourth independently challenged their evidence and
normative premises. All worked in detached worktrees at that commit. No tracked
source changed. The user's relocation of the unchanged reference-check block in
the primary molecule.rs was preserved and was not reported as a defect.

Scope: Molecule construction, reference checks, the authoritative integrity gate,
editor publication, checked mutation, and direct shared-check consumers. The
Reaction and ReactionSpan review covers their construction, publication, and
projection gates, with downstream operations traced where they expose a gate
failure. The Molecule admission audit also traced downstream identity, frame,
symmetry, and remapping consumers. Neither cycle is a full audit of reaction
application or Python. Doc 213's proposed editor lifecycle is a future consumer,
not missing implemented behavior.

The reviewers consulted the living guides, applicable skills, the status index,
docs 211/214/215, and the whitepaper's design intent. Each defense checked whether
the governing documents deliberately deferred or required the observed behavior.
No such deferral defeats the findings below.

## Semantic admission audit

The governing discussion is **215 — Integrity closure and minimization**.
The living integrity guide's admission test requires a concrete downstream failure
or repeated prerequisite; convenience at construction is insufficient. The audit
traced every current predicate to a consumer, rather than treating the guide's
inventory as proof of its own necessity.

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
| DuplicateParticipant | Actual-atom frames need a unique positional action, and simple incidence matching must distinguish occurrences. This covers localized and noncovalent bond self-loops, repeated dative donors, donor/acceptor overlap, and repeated aromatic/multicenter participants (incidence.rs:250–300; molecule/pushout.rs:195–247). Stereo actual-ligand uniqueness is already implied by complete-ligand uniqueness; atom-site exclusion remains a separate condition. |
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

Doc 215 rejects donor/acceptor overlap because it breaks incidence matching.
Two bounded public graph-core probes confirmed the dependency: a self-loop query
falsely matches a loop-free triangle three ways; a two-node donor/acceptor incidence graph with parallel,
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

### F2 — Dative identity checks build unnecessary collections

**Confirmed; L2.** integrity.rs:135–143 builds a participant HashSet, then a donor
BTreeSet, then clones the tree into the identity table. The original tree is retained
only for a possible sorted-donor diagnostic. No tree operation is needed after
constructing this order-independent key.

**Defense:** BTreeSet expresses set identity directly, and the separate participant
check reports an early duplicate without sorting the rest of the row.
**Unrestricted baseline:** sort one copied donor Vec, detect adjacent duplicates, then binary-search
for the acceptor. If donors repeat, use the existing original-order check only on
that failure path to select the same first duplicate. This preserves donors-before-
acceptor precedence. Move the unique sorted Vec into the identity table; reconstruct
sorted diagnostic donors only on an identity collision. Stored frames are untouched.
This eliminates both the successful-path participant HashSet and the trees, while
retaining one Vec per nonempty identity key. The measured benefit and slower early
duplicate rejection are recorded below. Doc 215 requires the identity and distinct
participants, not these collection types.

For at most 128 atoms, use an exact two-word donor bitmap as the identity key;
the unrestricted sorted-key proposal above is the fallback. Both preserve the
same diagnostics without per-row participant HashSets or trees.

### F3 — Stereo frames are copied for read-only checks

**Confirmed; L3.** integrity.rs:196/220 calls ligand_frame, whose implementation in
view/stereo.rs:717–721 collects a Vec. StereoAtoms::ligands and StereoBonds::ligands
already return the identical stored slice (stereo.rs:76/304).

**Defense:** the view method is convenient and valid frames are small; malformed
frames can be larger. **Proposal:** borrow the existing slices, preserving frame
order and all checks. No new accessor or public surface is necessary. Nonempty
frames avoid a temporary allocation; the gate never needs ownership of the copy.

### F4 — Multicenter identity checks build unnecessary collections

**Confirmed; L4.** integrity.rs:166–172 allocates a participant HashSet, an
original-order Vec for a possible diagnostic, and a BTreeSet identity key.

**Defense:** the original-order Vec preserves error data, and the tree erases frame
order for identity. **Unrestricted baseline:** use one sorted Vec for duplicate detection and
identity, with the same error-only original-order duplicate check as F2. On identity
collision, collect the original participant order for the diagnostic. Empty rows
remain supported; two empty rows still have identical identity. Preserve identity
checking before electron-count shape checking. This removes the per-row HashSet
and tree, not the outer identity table or its owned vector keys. The early-error
tradeoff is the same as F2.

For at most 128 atoms, use an exact two-word participant bitmap as the identity
key instead; the sorted-key proposal is the unrestricted fallback. Keep the
outer identity table and original-order diagnostics in either case.

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
**Proposal:** scan for an Atom-kind ligand equal to the site, returning
DuplicateAtom after the approved rename. Allow virtual ligands anchored there.
Keep frame errors first. Basis: the review guide's class-3 revalidation and local efficiency criteria;
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

### F10 — Assembled localized-bond endpoints are checked again

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
ConstitutionColoring::full also calls view.kind() while coloring stereo entities,
before graph_symmetry reaches stereo_center.

**Defense:** symmetry needs a geometry to grade orientation. However, its existing
nonliteral handling already permits no orientation contribution; graph_symmetry has
no declared determined-kind precondition. Kindless frames are expressly admitted
by 214 and the integrity guide. **Proposal:** handle absent kind in this consumer,
as it handles unresolved cosets. Do not reject these molecules at publication.
Add a focused public regression. This is a consumer contract defect, not a missing
integrity check or a complete audit of stereo-view getters.

### F15 — Aromatic participant checks allocate a fresh HashSet per system

**Confirmed.** integrity.rs:150 calls check_unique_participants, creating a local
HashSet for each system before inserting its atoms into the global membership set.

**Defense:** hashing has expected linear cost and finds the first repeated atom
without reading the rest of a malformed row. **Unrestricted baseline:** check adjacent duplicates
in a sorted temporary Vec; on failure use the original-order helper to select the
same error. Drop the Vec before the existing original-order membership loop.
Within-row duplicates still precede overlap, which still precedes electron-count
shape errors. Unlike dative/multicenter, this creates no retained identity key.

The comparison below supports this local change under code-reviews.md:16–24.
It retains one allocation per nonempty row and changes local expected-linear hashing
to O(k log k) sorting; it is not a universal speed claim for arbitrary row sizes.
Retaining the largest scratch buffer or enlarging the global table is unnecessary.
No contract or deferral in 214/215 requires a per-row hash table.

For at most 128 atoms, use bitmaps for both row uniqueness and global membership.
The unrestricted baseline applies above that size. Finish row uniqueness before
overlap checking; recover the first overlapping atom from the input order.

### F16 — Completed family checks retain their temporary tables

**Confirmed; lifetime.** The seven family-wide HashSets in integrity.rs:119–245
share the function's scope, retaining earlier identity keys while later families
are checked. None has a consumer after its own loop.

**Defense:** one scope avoids nested blocks, and an optimizer may affect allocation
behavior. **Proposal:** scope each independent family loop with its own table so
its storage is released on completion. This follows the review guide's local
efficiency criterion without adding helpers, visibility, or shared scratch state.
No numerical whole-gate peak-memory reduction was measured.

## Allocation choices by entity

This inventory concerns integrity-check scratch, not molecule storage. Empty
collections need not allocate; growing a table may allocate more than once.

| Entity | Current successful-path collections | Proposed disposition |
| --- | --- | --- |
| Atom | No atom-form checking collections | None needed. |
| Localized bond | One family-wide HashSet of endpoint pairs | Replace with lookup in the graph's sorted adjacency; no extra collection. |
| Dative bond | Family identity table; per-row participant HashSet, donor tree and tree clone | F2: bitmap key at ≤128 atoms, sorted key above; reserve the outer identity table from relation count and remove row HashSets/trees. |
| Aromatic system | Global membership HashSet; fresh participant HashSet per row | F15: row/global bitmaps at ≤128 atoms, sorted row plus global HashSet above. |
| Multicenter bond | Family identity table; per-row participant HashSet, diagnostic Vec and tree | F4: bitmap key at ≤128 atoms, sorted key above; reserve the outer identity table from relation count and construct original-order diagnostics only on failure. |
| Noncovalent bond | One family-wide HashSet of fixed endpoint pairs | Sort exact packed endpoint keys in an inline SmallVec; no heap for up to 16 bonds. |
| Stereo atom/bond | Family site HashSet; copied ligand Vec and two row HashSets | F3/F5–F7: borrow the frame, use the bounded uniqueness scan, and check atom-site exclusion separately. Reserve the stereo-atom site table; use existing site incidence for stereo bonds. |

Reference/constraint traversal and electron-count shape checks create no such
collections. First use of a stereo kind can initialize a shared coset space; that
is persistent algebra storage, not repeated per-entity scratch. F16 shortens the
lifetimes of the seven family tables without changing their algorithms.

The existing coincidence methods are not allocation-free substitutes. Their
adapters copy queries, and graph-core sorts copied query/candidate frames
(dative.rs:184–189; relation/fixed_var.rs:272–285; relation.rs:193–200).
Repeated lookup among many relations sharing an anchor can require quadratic
candidate scanning. Multicenter lookup also needs a participant anchor, so it
cannot alone detect two empty identical rows. The stored-key comparison below
examines a possible storage improvement; no storage/API change has been adopted.

## Nomenclature and visibility migrations

These are migrations, not additional correctness findings. The nomenclature guide
(lines 3–7) and review guide (90–94) distinguish later conventions from historical
defects. The validate_* reference helper names existed at cd64b530f (2026-08-05);
the Integrity check glossary entry appears at 1257445b0 (2026-08-08).

MoleculeIntegrityError renames:

| Previous | Current |
| --- | --- |
| BondsParallel | ParallelBonds |
| NoncovalentBondsParallel | ParallelNoncovalentBonds |
| DativeBondsIdentical | IdenticalDativeBonds |
| MulticenterBondsIdentical | IdenticalMulticenterBonds |
| StereoAtomSitesDuplicate | DuplicateStereoAtomSites |
| StereoBondSitesDuplicate | DuplicateStereoBondSites |
| DuplicateParticipant | DuplicateAtom |

The first six rows move the trailing adjective to the front, preserving the
existing terms and plurals. DuplicateAtom names the actual failure: an AtomId
occurs twice in one entity's atom references, including a stereo-atom site reused
as an actual atom ligand. Repeated complete stereo ligands already produce
DuplicateStereoLigand. Other variants retain their names. Payloads, diagnostics,
and rejection semantics are unchanged. The review evidence above uses the
pinned names.

- Move the five molecule validate_* reference functions out of `molecule.rs`
  into `molecule::integrity` and name them check_*: they check representation
  references, not constraint satisfaction. Apply the same terminology to the
  directly related span entry-reference function when touched. Keep the
  exhaustive constraint match. Only the entry and constraint checks with
  consumers outside `molecule::integrity` need crate-private access.
- Move the local stereo-form checks shared with Reaction out of
  `molecule::integrity` into a small `stereo::integrity` module. Its private error
  vocabulary maps into each aggregate's own public integrity error; it owns no
  aggregate gate. Keep `Molecule::check_integrity` in `molecule::integrity` and
  `Reaction::check_integrity` in `reaction::integrity`. Only entry points used
  across these modules need crate-private visibility; the rest stay private.
  No public checker object or exported helper API is needed.

## Open alternatives and rejected claims

| Question | Verdict and reason |
| --- | --- |
| Replace family-wide identity/membership/site HashSets | Use a two-word bitmap for aromatic membership at ≤128 atoms, sorted graph adjacency for localized-bond duplicate detection, and inline sorted keys for noncovalent bonds. Reserve the dative/multicenter identity and stereo-atom site HashSets from their relation counts; use existing site incidence for stereo bonds. Keep the aromatic HashSet above 128 atoms. Additional inline identity and site-bitmap paths are not selected. |
| Reuse participant scratch or combine aromatic ownership checks | **Not selected.** Reusing a HashSet can retain and repeatedly clear the largest table. A reused sorted Vec also retains the largest row. An aromatic owner map combines checks but enlarges every global entry and needs delayed-overlap handling. Short-lived row storage avoids that retention; bounded bitmaps avoid its allocation entirely when IDs fit. |
| Fuse the two constraint traversals (L6) | **Open.** Reference checking currently precedes stereo-domain checking for the whole constraint tree (integrity.rs:246–250,442–469). Plain leaf-by-leaf fusion changes error priority and must still protect indexing. Saving domain errors while checking remaining references adds complexity. No equivalent simplification has been established. |
| Narrow constructor reference preflight (L7) | **Open.** Graph::build_csr indexes bond endpoints (graph.rs:576–579), so deleting preflight can panic. Relation incidence construction stores/sorts opaque IDs, so not every early relation check is needed for storage safety. Narrowing the pass still changes cross-entity error priority. Preserve the full editor publication gate; do not add validity flags or witness types. |
| Remove the result gate from trusted remapping | **Open policy clarification.** molecule/remap.rs:129 calls try_from_arcs, which invokes the whole gate. Doc 215:96–98 says trusted transformations do not recheck outputs, but data-types.md:148–170 retains checked/asserted publication through the same gate, and 215 S2a explicitly removed only the source check. The retained output check is factual; calling it a settled defect would ignore that conflict. No unchecked constructor is proposed here. |
| Full candidate checks are inherently wrong | **Refuted.** Completed 215:198–221 and data-types.md:153–170 expressly settle candidate-and-authoritative-gate mutation. Incremental integrity requires a separate preservation argument. |
| The integrity module needs a new checker framework or forced split | **Refuted.** Its size and shared consumers do not establish a structural defect under code-reviews.md:63–68. Local functions and existing storage access suffice. |

## Evidence and limits

The initial review compilation and sampling ran in isolated worktrees on rustc 1.96.0,
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

### Localized-bond duplicate detection

A release-mode comparison used the production Graph CSR, with construction outside
the timed region. Each kernel checked self-loops and parallel endpoint pairs; the
common reference checks and other molecule integrity work were excluded. Nine
rotated-order timing samples of 20,000 calls each report medians with the normal
allocator; a separate instrumented build measured allocations. An exhaustive
11,111-graph check (all endpoint sequences of length 0–4 on four atoms) found no
accepted/rejected disagreement. The adjacency and local-band variants can report
a different first error when several defects coexist; exact error selection is not
required here.

The fixtures were paths on 8, 24, or 96 atoms; a 64-atom path plus 32 chords
`(i, i+2)` for `i=0..31`; a 96-leaf star; and a 96-edge path through atoms
`0, 20, ..., 1920` with the intervening atoms isolated. A separate late-error
fixture repeated `(0, 1)` after the 95 local bonds.

Nanoseconds per successful bond-family check, on arm64 with rustc 1.96.0:

| Method | 7-bond chain | 23-bond chain | 95-bond chain | 95 local bonds / 64 atoms | 96-bond star | 96 bonds / 1,921 atoms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Current pair HashSet | 265 | 1,193 | 4,957 | 4,968 | 5,012 | 4,944 |
| HashSet of 16-bit `(lower ID, difference)` keys | 157 | 834 | 3,251 | 3,271 | 3,399 | — |
| Stack triangular pair bitmap | 20 | 53 | 157 | 171 | 207 | — |
| Scan all sorted neighbor lists | 17 | 51 | 211 | 180 | 225 | 3,302 |
| Per-bond lookup in sorted neighbors | 20 | 65 | 267 | 342 | 751 | 270 |
| 16-neighbor bit band, adjacency for long bonds | 10 | 29 | 113 | 112 | 659 | — |

The 16-bit key is exact only when the atom namespace fits 128 IDs. It reduces
HashSet memory but retains 2–6 allocations per check; at 95 bonds its peak
requested heap was 592 bytes, versus 1,744 for endpoint-pair keys. The stack
triangular bitmap reserves 1,016 bytes and the 16-neighbor band 256 bytes; both
need a separate path for larger atom namespaces. A heap triangular bitmap needed
230,520 bytes for the sparse 1,921-atom fixture despite only 96 bonds. The full
neighbor-list scan also pays for isolated atoms. Packing local endpoints is useful,
but neither locality nor a small bond count bounds the atom namespace.

Per-bond adjacency lookup performs one binary search in an already-sorted
degree-sized slice and then checks the next entry. It needs no hashing or
scratch initialization. A full adjacency scan visits all atoms and both
directions of each bond; the triangular bitmap clears up to 127 words before
setting bits. These are logical operation counts, not measured CPU instructions.

**Selected:** reuse the graph's sorted adjacency for each localized bond. Find
the lower endpoint's neighbor group and reject when it contains a second entry
for the upper endpoint. This uses one algorithm for all graph sizes, allocates
nothing, and took 0.02–0.76 µs across these successful fixtures. It keeps the
same validity boundary and reports an actual duplicate endpoint pair; the first
error among multiple defects may change. The bit-band and triangular bitmap are
measured alternatives if a later whole-gate profile justifies a specialized
small-graph path. Across the six successful fixtures, this is 6.7–18.6 times
faster than the current pair HashSet: roughly an order-of-magnitude improvement
for this check. The timings do not measure the complete molecule gate.

### Noncovalent-bond duplicate detection

Noncovalent bonds are rows of a FixedRelationSet, not graph edges. Its incidence
index sorts relation IDs by atom, not by the other endpoint. Scanning an anchor's
incidence therefore repeats work when many relations share that atom. The generic
coincidence lookup also copies and sorts the query and candidate pairs; it is not
a cheap duplicate check.

A release-mode comparison used the production FixedRelationSet, with construction
outside timing. All methods iterated stored rows directly; view construction,
reference checks, and other integrity work were excluded. Seven rotated-order
samples report medians with the normal allocator; a separate instrumented build
counted allocations. All 69,905 sequences of 0–4 ordered endpoint pairs over
four atoms agreed on acceptance, including reversed frames, loops, and parallel
pairs. The expected common case is fewer than ten mostly disjoint noncovalent
bonds; larger fixtures test the unbounded path, not a claimed frequency.

Nanoseconds per successful noncovalent-bond check, on arm64 with rustc 1.96.0:

| Relations / shape | Current pair HashSet | Direct incidence scan | Inline sorted u64 keys | Heap sorted u64 keys |
| --- | ---: | ---: | ---: | ---: |
| 1 | 46 | 4 | 5 | 16 |
| 4 disjoint | 205 | 22 | 13 | 25 |
| 8 disjoint | 489 | 55 | 19 | 32 |
| 10 disjoint | 516 | 82 | 26 | 32 |
| 95 local / 64 atoms | 5,125 | 1,848 | 680 | 522 |
| 96-leaf star | 5,174 | 7,547 | 338 | 194 |
| 96 spread / 1,921 atoms | 5,145 | 1,642 | 340 | 196 |

The inline method uses SmallVec<[u64; 16]>, already a graph-IR dependency, with
capacity set to the relation count. It makes no heap allocation through 16 rows
and one allocation above that; at 95 rows the measured peak requested heap was
760 bytes, versus six allocations and 1,744 bytes for the HashSet. The 1,000-
and 5,000-row stress fixtures took 8.0/45.0 µs with inline sorted keys, versus
66.6/279.0 µs with the HashSet. Packing `(lower ID, upper ID)` into u64 covers
every AtomId without an atom-count cutoff; using their difference saves no space
in this width.

**Selected:** collect one exact packed key per noncovalent bond, sort the inline
buffer, and reject adjacent equal keys. This keeps one path for the common small
case and the larger boundary, without adding an index to relation storage.
Different first-error selection is acceptable when several defects coexist.
The timings isolate this check, not the complete molecule gate.

### Identity and participant collection kernels

A second bounded comparison tested the actual collection/control-flow shapes of
the dative, multicenter, and aromatic loops, using u32 for the equivalent AtomId
wrapper. Source files for these checks were unchanged between the reviewed commit
and HEAD 33ecfe5e0. Production code was not modified. Row reference checks and
electron-count checks are outside these kernels and must retain their current
positions when implementing the proposals.

The alternatives were: current code; tree keys moved without eager diagnostic
copies; sorted Vec keys with separate participant hashing; reused participant
HashSet; and combined sorted-key uniqueness/identity. Aromatic comparisons covered
current code, reused HashSet, reused sorted Vec, owner HashMap, and the selected
short-lived sorted Vec. All use standard collections, with no new dependency.

**324,672 exact modeled-result comparisons passed.** These are twelve alternative
comparisons for each of 1,456 one-row cases and 25,600 two-row cases: frames of
length 0..5 or 0..3 respectively over three atom IDs, with four acceptor choices.
They cover empty keys, reordered identities, first-duplicate selection, acceptor
conflicts, and duplicate-before-overlap precedence. Dative identity diagnostics
remain sorted; multicenter diagnostics retain input order. This bounded evidence
supports the direct equivalence argument; it does not replace production regressions.

Successful fixtures contain 64 disjoint rows of lengths 1, 4, 16, 64, or 256,
deterministically shuffled before measurement. A skewed fixture has one 4,096-member
row followed by 1,024 singleton rows. These shapes test algorithm costs, not a
claimed distribution of real molecules. Separate optimized rustc binaries measured
timing without instrumentation and allocations with a counting System allocator.
Seven rotated-order timing samples report medians; fixture construction is excluded.

Selected results per **64-row pass**, current → proposed. Time is microseconds;
allocation calls include reallocations. Peak is live requested allocation bytes,
excluding input storage, allocator metadata, size-class rounding, and transient
old/new-buffer overlap inside realloc; it is not process RSS.

| Family | Members per row | Time, µs | Allocation calls | Peak requested bytes |
| --- | ---: | ---: | ---: | ---: |
| Dative | 4 | 21.7 → 7.0 | 326 → 70 | 9,600 → 7,264 |
| Dative | 64 | 277.7 → 34.4 | 1,350 → 70 | 35,952 → 20,944 |
| Multicenter | 4 | 18.8 → 5.7 | 326 → 70 | 8,024 → 5,728 |
| Multicenter | 64 | 252.6 → 32.6 | 966 → 70 | 34,952 → 19,592 |
| Aromatic | 4 | 24.5 → 12.5 | 136 → 72 | 3,856 → 3,856 |
| Aromatic | 64 | 349.8 → 176.9 | 396 → 76 | 61,456 → 61,456 |

For the 64-member dative fixture, merely moving the tree costs 249.7 µs;
a sorted key with separate participant hashing costs 176.3 µs. Combining the two
checks costs 34.4 µs. Thus the recommendation removes the tree and redundant
successful-path collection, rather than only its clone. One owned Vec per nonempty
dative/multicenter key remains; the outer identity table remains too.

**Error-path cost:** in a single malformed dative row with an early duplicate,
current → combined time is 49 ns → 1.20 µs for 256 donors and 47 ns → 26.0 µs
for 4,096 donors. A late duplicate at 4,096 donors costs 176.3 → 185.7 µs.
Sorting the full row delays early rejection. The proposal favors successful
publication while preserving exact errors; retaining the separate uniqueness scan
would preserve early rejection at the higher successful-path cost measured above.
These error timings measure dative kernels, not every entity family.

For aromatic systems, dropping the sorted row before global membership insertion
avoids retaining the largest buffer. In the skewed fixture, peak requested bytes
are 61,456 for current and selected code, versus 77,840 for reused Vec, 102,424
for reused HashSet, and 110,608 for owner map. Selected time is 187.7 µs versus
305.0 µs current; calls fall only from 1,048 to 1,037, but cumulative requested
bytes fall from 192,664 to 102,476. Owner maps can be faster, but their larger
global entries and delayed-overlap logic are not needed for this local improvement.

These results justify the proposed collection choices. They do not establish a
whole-gate speedup or quantify F16's combined-family peak-memory reduction.

### Small-relation integrity checks

This is the primary comparison for choosing local check implementations. It
uses the expected sizes: 2, 4, 6, or 10 participants per relation, with 1, 4, or
8 relations independently, all atom IDs below 100. Larger participant counts in
the stored-key experiments below are boundary cases, not the expected workload.
For dative bonds, participant count here means donors; the acceptor is additional.

The alternatives reproduce the current collection/check shapes, use sorted
SmallVec keys, or use two-word bitmap keys. SmallVec 1.15.2 has inline capacity
10 with union/const_generics enabled: no tested row spills to the heap. Its
elements are a transparent ID newtype with derived Hash/Eq/Ord matching AtomId;
primitive-u32 slice hashing can differ and is not used for the sorted keys here.
Inline capacity 10 is a measured layout choice, not a lower bound on sorted-key
size or a settled production choice.

The kernels include dative donor uniqueness, acceptor exclusion, and full
identity; multicenter participant uniqueness and identity; and aromatic row
uniqueness and cross-row overlap. Reference checks and electron-count shape
checks are excluded identically and must retain their original positions.
Timing includes constructing and dropping all check scratch; no key is prebuilt.
It is not whole-molecule construction or whole-gate timing.

Sorted rows are copied inline, sorted, and scanned for adjacent duplicates.
On failure only, an original-order prefix scan selects the first duplicate.
Bitmap construction tests each bit before setting it. Both preserve donor-first
then acceptor checking, sorted dative identity diagnostics, and original-order
multicenter identity diagnostics. Aromatic checking finishes row uniqueness
before testing overlap; on overlap, it scans input order for the same first atom.

**7,944 oracle cases and all 49 benchmark fixtures passed** for every applicable
kernel in both runs. The independent oracle uses ordered prefix scans and direct
set-membership comparisons. Exhaustion covers two rows of lengths 0–2 over
[0,63,64,99], including all dative acceptor choices; longer explicit cases check
competing errors, reordered identities, and diagnostic order. Bitmap indexing
assumes the established ID bounds; these bounds are not a new admission limit.

Selected successful passes with **four relations**, nanoseconds per complete
family kernel; seven rotated-order medians on the compiler/platform above:

| Family | Participants per relation | Current | Inline sorted | Bitmap |
| --- | ---: | ---: | ---: | ---: |
| Dative | 2 donors | 657 | 304 | 289 |
| Dative | 4 donors | 1,125 | 338 | 296 |
| Multicenter | 4 atoms | 1,082 | 300 | 226 |
| Multicenter | 6 atoms | 1,271 | 344 | 229 |
| Aromatic | 6 atoms | 1,045 | 597 | 28 |
| Aromatic | 10 atoms | 2,048 | 1,111 | 78 |

The bitmap was faster across all measured valid shapes, but the additional gain
over inline sorting is modest for short dative/multicenter rows. At four
participants and 1/4/8 relations, dative sorted → bitmap times were
74→57 / 338→296 / 776→694 ns; multicenter 71→50 / 300→226 / 681→540 ns.
Removing the current per-row HashSets and trees provides most of their gain.

For aromatic systems, the sorted column retains the global HashSet. A matched
control uses sorted rows with bitmap global membership. With four ten-member
systems, times are **1,111 ns sorted+HashSet → 146 ns sorted+bitmap → 78 ns
bitmap+bitmap**. Global membership is the larger cost; comparing row keys alone
would miss it. After row uniqueness, the bitmap path performs two-word overlap
and union instead of hashing every member into the global table.

Memory was measured in a separate allocator-instrumented binary. Selected
four-relation cases, current / inline sorted / bitmap:

| Family and row size | Allocation/reallocation calls | Retained scratch bytes | Peak requested heap bytes |
| --- | --- | --- | --- |
| Dative, 4 donors | 22 / 2 / 2 | 544 / 512 / 256 | 692 / 700 / 316 |
| Multicenter, 4 atoms | 22 / 2 / 2 | 480 / 448 / 192 | 556 / 604 / 220 |
| Aromatic, 10 atoms | 17 / 5 / 0 | 376 / 376 / 16 | 496 / 496 / 0 |

Retained bytes include the resulting table's inline header and live requested
heap capacity. Peak counts heap only, excluding input frames, stack temporaries,
allocator metadata/rounding, and allocator-internal realloc overlap. The aromatic
sorted+bitmap control also has zero allocations and 16 retained bytes; its row
buffer is temporary. Bitmap row/global values are each 16 bytes; SmallVec is
48 bytes, or 56 with a dative acceptor, versus 16/24 for bitmap keys.
Fewer allocations do not imply lower memory: for one four-donor relation, inline
sorting increases retained scratch from 244 to 284 bytes because of the reserved
inline capacity; the bitmap uses 156 bytes. The table above also exposes higher
peak heap use for inline sorting during table growth.

Error timings use four four-member rows. An immediate duplicate costs
17/17/16 ns for sorted dative/multicenter/aromatic checks and 4.6/4.7/5.6 ns for
bitmaps. Later duplicates, identity collisions, acceptor conflicts, overlap,
and competing-error cases also retained exact results; no early-error regression
was observed in this bounded comparison. That does not extend to arbitrary sizes.

Valid frames use disjoint slices of the permutation 37i mod 100, independently
shuffled; acceptors are outside their donor row. Thus all eight ten-member
aromatic rows fit without overlap. An independent reviewer checked the source,
oracle, fixture validity, hashing behavior, timing, and memory accounting.
The result supports temporary bitmap checks for a fitting namespace, with a
general path for other IDs. It does not require persistent stored keys or a
graph-core API change. Measurements assume concrete key types; a larger tagged
key representation would need its own memory accounting.

### Stored order-independent keys

**For the expected common case of molecule IDs below 100, an exact two-word
bitmap is a strong alternative to sorted Vec keys.** It occupies 16 bytes,
needs no key allocation, and compares two words without sorting or collision
fallback. It reduced both construction work and memory in the identity-table
experiment below. It did not win every prepared comparison. No storage change
or dependency has been adopted.

#### Stored-key construction and repeated comparisons

These earlier representation probes use primitive u32 keys. They isolate key
costs and do not replace the AtomId-equivalent combined-check comparison above.
The 16/64-member cases explore size effects beyond the usual relation sizes.

The molecule's IDs are contiguous; an individual relation's participants need
not be. Fixtures therefore use 4, 16, or 64 distinct IDs chosen as 37i mod 100,
independently shuffled for storage and query. Unequal queries replace one ID by
an absent ID still below 100, preserving cardinality. The bitmap is two inline
u64 words and can represent any unique-participant subset of IDs 0..127.
The original ordered frame remains present for every alternative.

Compared: current candidate copy/sort; stored sorted Vec; dense Vec of words;
Roaring; 64/128-bit fingerprints alone or with sorted Vec; and the inline bitmap.
Fingerprints pair cardinality with a wrapping sum of per-ID XXH3 hashes of
little-endian u32 bytes. A mismatch proves inequality; a match requires exact
sorted comparison. Both widths passed forced-collision fallback checks. No
collision-rate guarantee is inferred from this experimental aggregation.

Seven rotated-order medians used an optimized standalone executable with
black-boxed inputs/results. Repeated inputs were hot and single-threaded.
Allocator counts came from a separate binary. These are kernels, not complete
molecule operations or measured CPU instruction counts. Independent refutation
checked the code, exactness, and accounting.

Selected results, nanoseconds per operation. “Fresh equal” includes building
the query key; “prepared” excludes it. Build/drop includes destroying the key.
Bytes include inline headers and live requested heap capacity, additional to the
original frame; allocator metadata, rounding, and query storage are excluded.

| Participants | Key | Build/drop | Prepared equal | Prepared unequal, late difference | Fresh equal | Retained bytes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 4 | Sorted Vec | 17.9 | 1.9 | 2.3 | 29.2 | 40 |
| 4 | Fingerprint64 + Vec | 22.6 | 1.9 | 0.7 | 33.7 | 56 |
| 4 | Inline bitmap | 8.2 | 3.5 | 3.5 | 8.3 | 16 |
| 16 | Sorted Vec | 55.9 | 3.1 | 3.1 | 76.5 | 88 |
| 16 | Fingerprint64 + Vec | 75.8 | 3.2 | 0.7 | 93.9 | 104 |
| 16 | Inline bitmap | 29.6 | 3.3 | 3.4 | 32.9 | 16 |
| 64 | Sorted Vec | 253.5 | 6.8 | 7.1 | 274.0 | 280 |
| 64 | Fingerprint64 + Vec | 325.3 | 7.1 | 0.7 | 343.5 | 296 |
| 64 | Inline bitmap | 121.7 | 2.9 | 3.4 | 131.7 | 16 |

Current prepared equality costs 19/58/263 ns at these sizes and allocates one
candidate Vec. Fingerprint64 alone occupies 16 prototype bytes and cheaply
rejects differences, but equal/colliding candidates still incur that copy/sort.
Fresh unequal queries cost about 7/21/77 ns for fingerprint64, versus roughly
8–14/33/88–103 ns for the inline bitmap. Thus the bitmap's advantage is exactness,
construction, and memory together, not universally fastest negative filtering.
Fingerprint128 adds construction/storage cost without a measured rejection benefit.

Roaring occupies 192/216/312 bytes and takes 77/210/798 ns for fresh equality.
The heap-backed dense bitmap occupies 40 bytes and takes 28/58/194 ns. These
results favor the inline representation over those bitmap containers in this
bounded namespace. Whole-key replacement costs 9/32/85 ns for inline bitmap,
16/59/276 ns for sorted Vec, and 20/73/349 ns for fingerprint64+Vec. This rebuilds
the key and drops the old one; it does not measure relation or incidence updates.

**Integrity identity tables.** Pairwise equality alone does not predict gate
performance: the gate already hashes keys. A separate pass inserted 64 distinct
identities with n participants (37i+r) mod 100. Cold passes build owned keys from
raw frames. Warm passes construct a fresh HashSet borrowing precomputed keys.
All use the ordinary default HashSet hasher; fingerprint keys hash only their
cardinality/digest and retain exact equality. Times are microseconds per pass:

| Participants per row | Cold sorted | Cold fingerprint+Vec | Cold bitmap | Warm sorted | Warm fingerprint+Vec | Warm bitmap |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 5.44 | 6.46 | 4.92 | 3.89 | 3.95 | 3.86 |
| 16 | 9.77 | 9.55 | 5.48 | 6.06 | 3.92 | 3.86 |
| 64 | 30.59 | 24.60 | 10.95 | 16.19 | 3.91 | 3.97 |

Cold bitmap tables make 6 allocation/reallocation calls and retain 2,232 bytes
at every size. Sorted keys make 70 calls and retain 4,280/7,352/19,640 bytes;
fingerprint+Vec also makes 70 calls and retains 6,328/9,400/21,688 bytes. These
totals include table capacity and inline header, not the input frames. Warm
tables each retain 1,208 bytes, **in addition to their precomputed keys**: bitmap
keys occupy 1,024 bytes for 64 rows; sorted keys occupy 2,560/5,632/17,920 bytes,
and fingerprint+Vec 3,584/6,656/18,944 bytes. An outer Vec adds its 24-byte header.

The table pass checks whole-identity uniqueness, not participant uniqueness or
other integrity predicates. A separate checked bitmap builder tests each bit
before setting it and reports the first duplicate in input order. It needs no
second pass or allocation: valid inputs cost 8.5/23.3/123.3 ns; an immediate
repeat at the second position costs 2.1–4.0 ns. This shows that key construction
and participant-uniqueness checking can share one pass. It is not a claim that
checking is free: small differences also reflect compiler/harness behavior.
Reference/range preflight must precede indexed word access, as it does before
current participant checks. Dative acceptor exclusion remains a separate role
check, which can test the donor bitmap.

#### Operation counts and storage placement

For C equal-length candidates and H fingerprint matches (including collisions):

| Retained representation | Candidate work after shared query preparation |
| --- | --- |
| None: current | C copies, allocations, and sorts, then exact comparisons. |
| Sorted frames | No candidate allocations or sorts; exact comparisons remain. |
| Fingerprints only | C fixed-size filters; H candidate copies, allocations, and sorts. |
| Fingerprints + sorted frames | C fixed-size filters; H exact comparisons; no candidate allocation or sorting. |
| Inline bitmap | Two-word exact equality; no candidate allocation, sorting, or fallback. |

Bitmap construction zeroes two words and sets one bit per participant;
uniqueness adds one test per bit. Fingerprints hash every participant; sorted
keys copy and sort. A fingerprint query's sorted copy can be delayed until its
first digest match and then reused. Prepared probes build it once eagerly;
fresh-query probes skip it on mismatch. Vec equality may compare several values
per instruction, so these logical counts are not scalar instruction counts.

Standalone key layouts are not compulsory. Graph-core already stores variable
frames flat with offsets: sorted copies could reuse those offsets, avoiding
per-row Vec headers/allocations. Existing row lengths could supply fingerprint
cardinality, leaving only 8/16 digest bytes. Packed sorted storage was not measured.
Four sorted u32 values and the inline bitmap both have 16-byte payloads; the
40-byte Vec result is not an information lower bound.

A persistent key belongs with its owning storage, which must maintain it on
addition, participant replacement/restoration, reference remapping/compaction,
and algebraic construction. Payload mutation and participant permutation need
not rebuild it. Gate-local keys avoid persistent maintenance but repeat construction
at each gate. The integrity choice above is gate-local; a persistent graph-core
index is a separate, unselected design.

Generic graph-core coincidence preserves **multisets of complete participant
values**, separately for each factor. Plain atom bitmaps erase multiplicity and
non-reference data, so they cannot replace that generic contract. Sorted keys
and fingerprints with exact multiset fallback can. Dative keys retain acceptor
versus donors; stereo sites, ligand kinds, and factor boundaries remain distinct.
Keys compare IDs within one namespace or after explicit transport. Fixed graph
endpoint pairs already have an exact two-u32 key and do not need a larger bitmap.

#### Larger-ID boundary cases

IDs below 100 are the user's expected common case, not a new limit on admitted
molecules. Larger namespaces require another representation. The first probe
therefore remains useful as boundary evidence, rather than the primary workload.
It compared four nearby IDs [8,10,12,14]; four spread IDs
[8,262152,524296,1048568]; 64 IDs 8+16381i; 4,096 consecutive IDs starting at 257;
and 4,096 IDs 8+251i. Equal queries reorder the same IDs; negative queries replace
the minimum/maximum by its predecessor/successor without changing cardinality.

| Shape | Key | Prepared equal, ns | Build/drop, µs | Retained bytes |
| --- | --- | ---: | ---: | ---: |
| 4 spread IDs | Sorted Vec | 1.9 | 0.017 | 40 |
| 4 spread IDs | Dense bitmap | 3,235 | 1.22 | 131,096 |
| 4 spread IDs | Roaring | 25.1 | 0.119 | 216 |
| 4,096 consecutive | Sorted Vec | 435 | 25.5 | 16,408 |
| 4,096 consecutive | Dense bitmap | 14.7 | 4.6 | 576 |
| 4,096 consecutive | Roaring | 3.7 | 43.6 | 200 |
| 4,096 spread IDs | Sorted Vec | 455 | 25.5 | 16,408 |
| 4,096 spread IDs | Dense bitmap | 3,266 | 6.6 | 128,512 |
| 4,096 spread IDs | Roaring | 226 | 45.5 | 16,536 |

Roaring's consecutive set uses one run; its spread set uses sixteen arrays.
Fresh equal consecutive queries still cost 43.0 µs with Roaring, versus 26.7 µs
sorted and 4.7 µs dense: compressed comparison does not eliminate key construction.
A single inline u64 handles the four-nearby-ID case in 8 bytes, building in 1.9 ns
and comparing in 1.7 ns, but cannot represent the far-ID fixture.

The scratch crate used roaring 0.11.5 and xxhash-rust 0.8.18, release thin LTO and
one codegen unit, on the compiler/platform stated above. Roaring construction
includes sorting, from_sorted_iter, and optimize; dense construction includes
maximum-ID discovery and zeroing every represented word. Memory comes from the
allocator ledger, not serialization size or library statistics. See the
[Roaring API](https://docs.rs/roaring/0.11.5/roaring/bitmap/struct.RoaringBitmap.html)
and [format specification](https://github.com/RoaringBitmap/RoaringFormatSpec).
All exactness checks passed, including empty sets, reordered equality, forced
fingerprint collisions, word boundaries 63/64, and bitmap multiplicity counterexamples.
No whole-gate speedup or automatic choice of key representation is claimed.

### Retained family sets and large-namespace membership

A follow-on comparison used an ignored scratch executable linked to the current
graph-core relation storage, compiled in release mode with rustc 1.96.0 on arm64.
Seven timed samples per case report the median. Inputs and relation incidence were
built outside the timed region. A separate allocator-instrumented build counted
allocation/reallocation calls and requested heap bytes for one check. These are
gate-local kernels, not a whole-Molecule timing or a workload-frequency estimate.
The executable checked both distinct and duplicate identities/sites, the 63/64/127
bit boundaries, and aromatic membership at atom ID 128. Production code was not
changed.

The identity cases have four distinct participants per row, IDs below 100, and
distinct complete keys. All alternatives build the proposed exact bitmap key;
this isolates the outer table after F2/F4, not the current BTreeSet-based gate.
The inline alternative stores keys in a SmallVec with capacity 16 and scans
earlier keys; the reserved alternative uses
HashSet::with_capacity(relation_count). The growing HashSet starts empty.
Times are nanoseconds per bitmap-key building and uniqueness pass:

| Family, rows | Growing HashSet | Reserved HashSet | Inline scan | Allocation calls: growing / reserved / inline |
| --- | ---: | ---: | ---: | --- |
| Dative, 8 | 644 | 290 | 47 | 3 / 1 / 0 |
| Multicenter, 8 | 412 | 214 | 50 | 3 / 1 / 0 |
| Dative, 64 | 5,227 | 1,920 | 1,458 | 6 / 1 / 2 |
| Multicenter, 64 | 3,461 | 1,468 | 1,819 | 6 / 1 / 2 |

At eight rows, the reserved tables request 536 bytes for dative and 280 for
multicenter, versus 948 and 500 for the growing tables. The inline capacity
itself occupies 512 bytes for dative pair keys or 256 bytes for multicenter
keys, plus the SmallVec header; the allocation ledger excludes stack storage.
A single relation needs no cross-row identity table at all.
Reserving the current table is a small, general improvement. The inline scan is
substantially faster for small families, but its quadratic comparison count calls
for a large-family path. The crossover and whole-gate benefit are not established.

The stereo cases use 100 possible atom/bond sites, four ligand anchors per
relation, and distinct sites. The scratch relation sets use NodeId for those
anchors, which contributes the same node references as StereoLigand; view-wrapper
costs are excluded. The incidence alternatives query the existing graph-core index;
the atom query filters out relations using the atom only as a ligand and checks
earlier IDs. The bond index directly records site-bond incidence. A two-word
bitmap covers the same bounded namespace. Times are nanoseconds per site-family
pass; each incidence and bitmap pass made no heap allocation:

| Entries | Current HashSet | Reserved HashSet | Atom incidence | Bond incidence | Bitmap |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 156 | 64 | 42 | 15 | 4 |
| 8 | 359 | 122 | 105 | 42 | 9 |
| 16 | 695 | 200 | 274 | 108 | 16 |
| 64 | 2,934 | 718 | 1,485 | 741 | 65 |

For eight entries, the current site HashSet made three allocation/reallocation
calls; reserving made one. The existing bond incidence index is a simple
allocation-free replacement at any site ID. Atom incidence also needs no new
storage or namespace branch, but is slower than a reserved HashSet at 16 and 64
entries. The bitmap is fastest locally, but needs another path above 128 atom or
bond IDs. One stereo entry needs no site-uniqueness table.

For aromatic membership above 128 atoms, a fresh dynamic Vec<u64> bitmap was
compared with the current HashSet. Only global membership was timed; row
uniqueness is common to both alternatives. Distinct members were spread across
the indicated atom namespace:

| Atoms / aromatic members | HashSet ns, calls, requested bytes | Dynamic bitmap ns, calls, requested bytes |
| --- | --- | --- |
| 129 / 8 | 371, 3, 164 | 24, 1, 24 |
| 10,000 / 8 | 366, 3, 164 | 40, 1, 1,256 |
| 10,000 / 128 | 5,628, 7, 2,596 | 130, 1, 1,256 |

The dynamic bitmap is faster in these kernels. Its memory follows the whole atom
namespace, however: the 10,000-atom, eight-member case requests 1,256 bytes
instead of 164. The sparse large-molecule fallback therefore remains a HashSet.
With at most one aromatic system, the global overlap check can be skipped after
within-row uniqueness. No additional namespace-dependent policy is justified by
these isolated timings.

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

## Reaction and reaction-span gates

Three review agents examined contracts, structure/efficiency, and tests/generators
in detached worktrees at `33ecfe5e03d714d2ac801015d20fa00243ea3a82`; a fourth
independently challenged their findings and normative citations. The agents
read the development guides, docs 214/215, the status index, and the whitepaper.
All tracked source worktrees remained clean. Exact public-API probes reproduced
the five Reaction failures below. A separate sample of 512 comprehensive
reactions and 512 independent span entries passed publication and projection
checks; all entity families appeared. Those samples do not establish complete
coverage. No production code changed.

### Adversarial verdicts

**Reaction integrity.** Under the contract reviewed at this commit, five
omissions contradicted [data-types](../docs/development/data-types.md)
(lines 16–34, 85–90) and [doc 215](215-integrity-minimization-2026-08-28.md)
(lines 50–66). R2–R5 remain publication gaps. R1's observed behavior remains a
defect, but the later R10 review questioned whether its constructor check is
needed. The check remains for this work. Product materializability, other
old/new value continuity, and constraint satisfiability remain deferred.

- **R1 — Confirmed, high: created-bond removal incidence.** In
  `reaction/integrity.rs:298–307`, a bond removal whose id is absent from the
  lhs skips comparison with a same-id addition. Adding bond 0 on atoms 0–1 and
  removing bond 0 on atoms 0–2 passes `Reaction::try_new`; normalization silently
  cancels both because `delta.rs:1895–1931` ignores the removal payload. The
  strongest defense is that delta old/new continuity is lazy. That does not
  justify silently discarding the recorded incidence. Check the added bond's
  endpoints at construction under the retained R10 boundary, and make bare
  delta normalization verify the removal before cancellation.
- **R2 — Confirmed, medium-high: positional electron-count lengths.** Aromatic
  and multicenter Add, Remove, and ModifyField deltas can carry two literal
  counts for a three-atom frame; `try_new` accepts them while even identity
  frame transport fails. `reaction/integrity.rs:114–183` checks references but
  not these payload shapes. A prospective product may be nonmaterializable,
  but its delta's positional values still need one meaning in the recorded
  frame ([integrity guide](../docs/development/integrity.md), lines 76–77;
  [data-types](../docs/development/data-types.md), lines 139–146). Add the direct
  `ReactionIntegrityError::ElectronCountLengthMismatch` variant with the
  Molecule error's `{ entity, participants, electron_counts }` fields.
- **R3 — Confirmed, medium: stereo Remove configuration domain.** A removal
  with valid site and ligands but an out-of-range tetrahedral coset passes
  `try_new`; later preconditions and span conversion reject it, and identity
  reframing fails. `reaction/integrity.rs:250–289` has no Remove payload check.
  Deferred old-value continuity does not make the coset denote a configuration
  ([data-types](../docs/development/data-types.md), lines 117–123). The symmetric
  stereo-bond path has the same omission.
- **R4 — Confirmed, medium: stereo ModifyField configuration domain.** The
  same out-of-range coset passes `try_new` as a new field value because
  `reaction/integrity.rs:278–286` checks only old/new kind continuity; the
  stereo-bond branch does the same. Application checks the domain later, but
  identity reframing already fails on the published Reaction. This is a
  representation-domain failure, not prospective-product consistency
  ([data-types](../docs/development/data-types.md), lines 117–123).
- **R5 — Confirmed, medium: stereo constraint positions.** Either a stereo
  ModifyConstraint or a top-level stereo ConstraintDelta with a tetrahedral
  topicity pair at positions 0 and 99 passes `try_new` and preconditions, then
  fails span conversion and identity reframing. `reaction/integrity.rs:258–273,
  400–445` checks references and site kind, not the payload against its owning
  ligand frame. Satisfiability is deferred, but an out-of-frame position is a
  tier-1 error ([integrity guide](../docs/development/integrity.md), line 103).

**Reaction integrity overreach.**

- **R11 — Confirmed: determined old/new stereo-kind disagreement.**
  `reaction/integrity.rs:74–85,278–286` rejects a ModifyField whose individually
  valid old and new configurations assert different kinds. The delta is still
  interpretable: `StereoAtomDelta` and `StereoBondDelta` transport the two values
  independently, and delta normalization can retain the change. Execution and
  two-sided materialization cannot preserve one stereo entity across a kind
  change, but those operations already reject it through
  `stereo_delta_domains_are_valid` and ReactionSpan publication. Construction
  also accepts an undetermined-to-Axial change against a Tetrahedral lhs, which
  has the same materialization failure. Checking only determined old/new kinds
  therefore establishes no general Reaction invariant used by ordinary
  operations. Remove this constructor rejection and its Reaction error variant;
  retain validation of each value against its frame and retain the distinct
  ReactionSpan `StereoKindModified` check. No foldability check moves into
  Reaction integrity.

**ReactionSpan documentation.** Both claims fail on the checked, valid span in
`reaction_span.rs:5133–5203`, which puts an Added atom before Unchanged atoms.

- **R6 — Confirmed:** the type rustdoc at `reaction_span.rs:60–62` describes
  lhs anchoring as a type invariant. That is true for `superimpose` output, but
  raw construction permits any valid dense union order
  ([data-types](../docs/development/data-types.md), lines 1064–1074).
- **R7 — Confirmed:** `to_reaction` rustdoc at `reaction_span.rs:1328–1332`
  omits lhs anchoring from its conditions for an exact span roundtrip. The
  preceding idempotence claim remains sound; the reordered example has no
  constraints or redundant Modified entries yet reanchors on roundtrip
  ([data-types](../docs/development/data-types.md), lines 1082–1100).

**Removal-incidence boundary.**

- **R10 — Keep the constructor check for now.** A `Remove`'s site and
  participants are an expected-old structural condition on the entity named by
  its id. The review did not show that checking agreement at construction
  unlocks broad downstream simplification. Moving the existing check changes
  the failure boundary and requires guarding each path that rewrites or discards
  the condition; that change is deferred. Complete the current constructor
  check for newly added bonds (R1). Bare `Deltas::normalize` has a separate
  cancellation defect and must check same-collection `Add`/`Remove` pairs.
  The copying measurements below describe the current cost, not a decision to
  change the constructor boundary.

**Bounded local efficiency.** The [review guide](../docs/development/code-reviews.md)
(lines 16–24) records local costs under a sound design below correctness and
layering. These measurements compare controlled scratch inputs, not production
frequency or peak memory.

- **R8 — Confirmed:** `projected_ids` hashes dense `0..N` IDs in eight maps.
  The prior 80-atom comparison below reduced construction from 61.9 to 38.2 µs
  with indexed vectors, including a nonempty-constraint case. Allocation calls
  increased from 187 to 199 in the unchanged case, so this is a time finding.
  The generic maps are simple and preserve absent-after-present translation.
- **R9 — Confirmed:** `project_entries` and `to_reaction` build eight
  correspondences solely for top-level constraint transport even when there
  are no constraints. The measured no-constraint projection guard reduced
  80-atom construction from 62.3 to 36.2 µs and 187 to 99 allocations; a
  separate `to_reaction` guard reduced it from 43.5 to 30.2 µs. The unconditional
  path is simpler and remains necessary with selected constraints; both side
  checks stayed in the probe.

The R10 constructor-copy measurements remain in
[Existing execution-cost evidence](#existing-execution-cost-evidence). They do
not by themselves justify changing the check's placement.

**Refuted and open.** A missing focused test for lhs/add and add/add id
collisions is useful future coverage, not a defect: the guard is correct, and
the focused-regression rule postdates it. The [review guide](../docs/development/code-reviews.md)
(lines 90–94) treats such historical conventions as migrations. Span side
checks at construction are required. Later `lhs`/`rhs` repeat a known side
check, but [data-types](../docs/development/data-types.md) (lines 169–172)
explicitly directs projection through the Molecule constructor; bypassing that
route needs a separate preservation argument. The large-file split is already
tracked in [168](168-api-hygiene-2026-07-27.md). No new naming, visibility, or
generator defect survived that cycle.

### Findings-blind Reaction rereview (2026-09-23)

Two new reviewers examined Reaction construction and its downstream consumers
in separate detached worktrees at `33ecfe5e03d714d2ac801015d20fa00243ea3a82`.
They received the review rules, guides, doc 215, and the whitepaper, but neither
this document nor the earlier findings. A third reviewer refuted their pooled
claims. All three worktrees remained clean; this pass traced concrete code-path
cases without running new tests or timings.

- **R1 behavior and R2–R5 publication gaps confirmed independently.** A same-id
  bond addition and removal with
  different endpoints is accepted and then silently canceled by normalization.
  Malformed electron-count, stereo-configuration, and stereo-constraint payloads
  are accepted even though identity frame transport cannot interpret them. The
  strongest defense is that product materialization and old/new continuity are
  deliberately lazy. That defense does not permit normalization to erase R1's
  recorded condition; R2–R5 concern payloads with no local interpretation.
- **R11 confirmed after challenge.** Both initial reviewers defended the
  constructor's `StereoKindModified` rejection. Refutation found that
  Tetrahedral and SquarePlanar atom configurations have the same degree-four
  parent action group, and `StereoAtomDelta::reframe_by` transports old and new
  independently (`umol-perm/src/class.rs:60–64,89–93`;
  `ir/delta.rs:1212–1220`). The kind mismatch can therefore be an interpretable
  raw delta. Application and span conversion have their own fallible boundaries.
  Conversely, an undetermined old value with a conflicting determined new kind
  already passes construction but fails application. The guide's claim that
  different kinds necessarily require different action groups is too broad;
  retain the R11 proposal and correct that rationale when implementing it.
- **Collision diagnostic narrowed.** The pinned integrity guide explicitly
  assigns duplicate additions to `InvalidReference`, so a separate variant is
  not required by the existing constructor contract. But the public rustdoc and
  Display say the colliding ID is *unavailable* (`ir/reaction/integrity.rs:50–53`),
  and Python and the DSL expose that false message. The proposed
  `DuplicateReference` split remains the chosen clearer error design; this
  rereview alone establishes the diagnostic defect, not the need for another
  variant.
- **Whole application-check removal refuted.**
  `stereo_delta_domains_are_valid` repeats some coset-range checks already
  established for the closed lhs and stereo additions, but also rejects empty
  sets/domains and term actions outside a kind's group, which construction
  deliberately permits (`ir/reaction.rs:257–327`;
  `ir/molecule/integrity.rs:591–668`). Its removal and modification compatibility
  checks are application-specific. Removing the whole predicate would change
  behavior; this review proposes no separate optimization of its repeated parts.

### Removal-incidence boundary review (2026-09-23)

Two fresh reviewers traced Reaction publishers and removal consumers at
`33ecfe5e03d714d2ac801015d20fa00243ea3a82`; a third refuted their claims.
They did not receive this document or the earlier findings. The question was
whether a `Remove`'s recorded incidence must agree with the lhs or same-reaction
`Add` at construction, independently of the existing rule in doc 215. The review
was read-only and ran no new tests.

**Review finding: removal incidence is an expected-old structural condition.**
The id names the entity; its recorded site and participants constrain the old
state, like the removal's old attributes or a `ModifyField` old value. A
contradictory condition could be stored and rejected when an operation needs to
align it with the owner. Editor removals already compare incidence and
attributes together at execution
(`ir/molecule/transact.rs:624–640`). Reference existence, unique added ids, and
locally interpretable payloads remain construction concerns.

No infallible public Reaction operation was found to require owner-incidence
agreement. `normalize_reaction_deltas` is shared by normalization, span
conversion, and application; `reframe_reaction_deltas` is the separate fallible
transport path (`ir/reaction.rs:381–545,716–869`). Thus the integrity guide's
concern about most operations repeating the same prerequisite is not established
for this check. The `FrameTransport` laws apply to compatible actions, and the
`Reframe` result laws to satisfiable receivers (`ir/traits.rs:175–215`).

The constructor check remains in place for now. This preserves its existing
failure boundary while other integrity work proceeds; the review did not
establish a broad downstream simplification that requires the check there.
Moving it later would require checking each path before it rewrites or discards
the recorded condition. Reaction normalization currently overwrites a
mismatched dative acceptor or stereo site with the owner's value; bond
materialization uses lhs endpoints by id (`ir/reaction.rs:427–545`;
`ir/reaction_span.rs:2310–2557`). The separate transport path would also need
to preserve or decline incompatible local frames.

Public `Deltas::normalize` independently folds a created entity's `Add` and
`Remove` to nothing without reading the removal's incidence or old attributes
(`ir/delta.rs:1895–1931,2866–2914,3018–3057`). Cancellation must require
matching structured incidence and `normalized_eq` of the removal's old value
with the value after intervening modifications, after transport into the Add
frame. Otherwise it returns `Contradiction`. Compatible local-frame
permutations still cancel. This fixes standalone deltas and the R1 case even
when delta normalization runs before Reaction construction.

### R2–R5 and R10 boundary rereview (2026-09-23)

Two findings-blind reviewers independently traced construction and consumers
at `33ecfe5e03d714d2ac801015d20fa00243ea3a82`; a third refuted their
claims. They read the review rules, guides, doc 215, and whitepaper without
this document. The pass was read-only and ran no new tests or timings.

- **R2 — Confirmed.** A literal electron-count vector missing a member cannot
  assign one value to each position of its local or owning aromatic/multicenter
  frame. `Reaction::try_new` accepts it in `Add`, `Remove`, and either side of an
  Electrons `ModifyField` (`ir/reaction/integrity.rs:87–110,114–183`). The
  defense that later transport or product construction can fail addresses
  execution, not this missing positional meaning. Old/new agreement stays lazy.
- **R3 — Confirmed.** A stereo `Remove` may agree with its owner yet carry a
  coset outside its own declared kind (`ir/reaction/integrity.rs:250–289`). The
  defense that a removal's old value can be an unsatisfied precondition does
  not make an out-of-domain index denote a configuration.
- **R4 — Confirmed.** A stereo `ModifyField` may carry an out-of-domain old or
  new configuration even when their kinds agree; the current check compares
  kinds only (`ir/reaction/integrity.rs:74–85,278–286`). Matching the old value
  against the owner stays lazy, but both values must denote configurations.
- **R5, entity constraint — Confirmed.** `ModifyConstraint` checks an optional
  kind's site admissibility but not either constraint's positions or permutation
  degree against the owner's ligand frame (`ir/reaction/integrity.rs:258–273`).
  An impossible constraint may be retained; a position outside the frame or an
  action of another degree has no local meaning.
- **R5, top-level constraint — Confirmed separately.** Recursive
  `ConstraintDelta` reference checks likewise admit an out-of-frame stereo
  position or action (`ir/reaction/integrity.rs:400–445`). Deferred constraint
  satisfaction does not supply the missing frame position.
- **R10 — Boundary change deferred.** The current constructor follows the
  committed integrity inventory and completed doc 215. No earlier public
  Reaction consumer was found to require owner agreement before fallible
  normalization, transport, span conversion, or application. This questions
  the placement but does not justify moving it during the current work.
  Separately, `Deltas::normalize` can erase a mismatched same-id `Add`/`Remove`
  before a Reaction exists; the constructor cannot cover that public operation
  (`ir/delta.rs:1895–1931,2866–2914,3018–3057`).

Doc 215 is Completed and records no R2–R5 deferral. Its constructor rule for
R10 remains the policy for this work.
The earlier public-API reproducers for R2–R5 remain the executable evidence;
this rereview tested check placement and producer/consumer premises.
The integrity guide must be revised with the implementation; completed doc 215
records the earlier constructor-boundary decision.

### Construction contract and current gaps

| Type | Open input and publication | Infallible result and deferred property |
| --- | --- | --- |
| `Reaction` | A closed lhs `Molecule` plus independently assembled `Deltas`; the target publication contract establishes references, unique additions, local payload domains, and removal-incidence agreement with lhs or a same-reaction `Add`. R1–R5 identify gaps and R11 an overreach in the current gate. `new` asserts the same gate as `try_new`. | Stored deltas need not be executable or materializable as a two-sided span. `to_reaction_span` may return `Contradiction`; host applicability and chemistry remain later checks. |
| `ReactionSpan` | `ReactionSpanEntries` may use any dense union order. `try_from_entries` checks union references, both projected molecules, and modified stereo-kind continuity. `from_entries` asserts the same contract. | `lhs`, `rhs`, `to_reaction`, and `correspondence` are infallible. Explicit equivalent `Modified` entries survive construction. No side entity or constraint is silently dropped. |

The checked Rust constructors return `ReactionIntegrityError` and
`ReactionSpanIntegrityError`, respectively; their asserted counterparts may
panic on malformed independently assembled input. Python construction maps
these failures to `ValueError`. No host-relative DPO condition, delta normal
form, satisfiability, or chemistry check was found in either eager gate.

The target `ReactionIntegrityError` has direct `InvalidReference`,
`DuplicateReference`, and `ElectronCountLengthMismatch` variants.
`IncidenceMismatch` remains a constructor error while R10 stays in the gate.
`InvalidReference { entity }` identifies an unavailable entity ID;
`DuplicateReference { entity }` identifies an `Add` ID already used by the lhs
or an earlier `Add`. Both carry the typed entity reference. The current
`InvalidReference` also reports collisions, with an inaccurate "unavailable"
diagnostic.

Replace `StereoIntegrityError(MoleculeIntegrityError)` with direct
`DuplicateAtom`, `DuplicateStereoLigand`, `StereoFrameDegreeTooLarge`,
`StereoKindSiteMismatch`, `StereoLigandArity`, `StereoCosetOutOfRange`,
`StereoPermutationDegree`, and `StereoLigandPositionOutOfRange` variants, using
the same fields as their Molecule counterparts. These are the eight failures
reachable through the shared local stereo checks; the wrapper currently admits
unrelated Molecule errors as well. Remove Reaction's `StereoKindModified`; keep
the distinct ReactionSpan variant.

Molecule and Reaction own separate aggregate checks: Molecule establishes its
stored sites and incidence, while Reaction resolves lhs/Add references and
checks delta frames and source incidence. The local stereo rules over a site,
ligand frame, and form are the same in both contexts. Keep one implementation of
those rules, but do not make their public error representation or placement under
Molecule define Reaction's error hierarchy. The internal helper return type and
module placement are implementation choices; no public helper type is needed.

### Existing execution-cost evidence

In `reaction/integrity.rs:88–110`, `source_frames(lhs)` runs before the first
delta pass. It clones all six overlay families into a `HashMap`, including each
variable participant or ligand list (`:576–624`). A later pass indexes the map
only for overlay removals (`:291–398`). The four delta passes keep a useful
ordering: collect all additions before checking references, then check stereo
and removal incidence. Combining those passes is not required by this finding.
For 20 dative removals, `unordered_ids` sorts two new vectors per removal
(`:689–694`). Fixed stereo-bond blocks likewise allocate four two-element
vectors per comparison (`:703–713`). These are secondary costs on removal
paths, not a reason to add a general indexing framework.

In `reaction_span.rs:158–213`, checked publication validates union references,
constructs the span, then constructs and discards both side molecules through
`Molecule::try_from_entries`. Both side checks are required: `lhs` and `rhs`
promise valid molecules. `project_entries` (`:1642–1880`) creates eight dense-id
maps, then unconditionally builds eight `Correspondence` values solely to map
the top-level constraint spans. `Correspondence::new` checks each map with two
`BTreeSet`s (`umol-graph-core/src/correspondence.rs:96–135`). With an empty
top-level constraint list, the correspondence is unused. `lhs`/`rhs` call this
projection again and reconstruct a checked molecule (`:1311–1318,1887`).
`to_reaction` also constructs eight ID maps for its delta conversion and an
eight-part correspondence used only by constraint deltas, then calls
`Reaction::new(self.lhs(), deltas)` (`:1333–1636`). Its correspondence is
unnecessary when there are no constraint deltas; the ID maps still serve the
delta conversion.

The existing `projected_ids` mapping intentionally places absent side entities
after the valid dense prefix. A retained reference to one then reaches the
side-molecule integrity error instead of being dropped or panicking. An
optimization must preserve that behavior for interleaved union IDs.

### Measurements and recommendation

The temporary probe used public graph-IR constructors in an isolated scratch
crate, release mode, with inputs prepared outside the timed region. Times are
medians of five batches; allocation counts are medians of seven isolated calls
and include reallocations. The probe used 12- and 80-atom chains, an
80-atom span with alternating retained/added/removed atoms and side-specific
bonds, and an 80-atom molecule with 20 valid dative overlays. These are
controlled costs, not a claim about production workload frequency or peak
memory. The standalone scratch lock resolved dependencies separately from the
workspace lock.

| Operation | Current | Temporary no-constraint paths |
| --- | ---: | ---: |
| 12-atom span construction | 9.5 µs, 101 allocations | 7.1 µs, 81 allocations |
| 80-atom unchanged span construction | 62.3 µs, 187 allocations | 36.2 µs, 99 allocations |
| 80-atom interleaved span construction | 58.7 µs, 184 allocations | 31.5 µs, 95 allocations |
| 80-atom unchanged `lhs()` | 27.6 µs, 87 allocations | 14.2 µs, 43 allocations |
| 80-atom interleaved `lhs()` | 25.4 µs, 85 allocations | 12.2 µs, 41 allocations |
| 80-atom unchanged `to_reaction()` | 43.5 µs, 135 allocations | 17.7 µs, 47 allocations |
| 80-atom interleaved `to_reaction()` | 48.4 µs, 144 allocations | 23.9 µs, 56 allocations |

The temporary projection change returned the already assembled side entries
with empty constraints before constructing their correspondence. A second
temporary guard skipped `to_reaction`'s correspondence when the span had no
constraints. On its own, the second guard changed the unchanged 80-atom
`to_reaction()` from 43.5 µs to 30.2 µs; the table shows both guards together.
Both molecule integrity checks remained in place. The changes were reverted.
For a future implementation, skip each correspondence when it has no selected
constraints or constraint deltas to transport; the measured
`self.constraints.is_empty()` guards are safe common cases. Retain the existing
path when those values are present. Do not cache projected molecules or
introduce a separate projection representation for this saving.

For `Reaction::try_new`, an 80-atom lhs without overlays and empty deltas took
1.8 µs with zero allocations. Adding 20 dative overlays while keeping deltas
empty took 3.8 µs and 24 allocations. Supplying 20 matching removals took
5.8 µs and 64 allocations. The additional no-delta work comes from the eager
frame copy. A temporary guard that built lhs source frames only when an overlay
removal occurred changed the 20-dative, empty-delta case from 3.68 to 2.09 µs
and from 24 to zero allocations. The 20-removal case remained about 5.6 µs
and 64 allocations. That guard skips copying the lhs frames, but the current
`added_entity_and_frame` still clones frames carried by `Add` deltas even when
there is no removal.

One dative removal from the same 20-overlay lhs took 3.96 µs and 26
allocations: 24 for the global source-frame copy and two for sorting the
compared donor lists. A separate dative-only scratch variant read the lhs
frame through its view instead of building the global map. It took 2.25 µs
and two allocations for one removal; 20 removals changed from 5.71 to 3.70 µs
and from 64 to 40 allocations. That probe alone does not measure borrowed
source lookup for the other five overlay kinds or for forward-referenced
additions. It does show that the no-removal guard alone leaves the isolated
removal cost in place. The remaining constructor checks retain their error
precedence and forward references to additions. Keep that ordering while
replacing the separate created-ID set and copied source-frame map with one
borrowed `Add` map keyed by `Entity`. Use it for duplicate detection, reference
availability, and created-entity source lookup; read existing source frames
through lhs views. This keeps the constructor check without cloning every lhs
frame. These measurements alone do not settle the eventual R10 boundary.

### Dense projection maps and repeated checks

`projected_ids` maps every dense union id through a `HashMap`, although the key
space is exactly `0..count`. A scratch copy replaced only these projection maps
with private typed indexing over `Vec<Id>`. It retained the present-id prefix,
absent-id suffix, and the checked side-molecule construction. A second variant
also skipped correspondences when no top-level constraint needed them. Medians
for the same 80-atom unchanged span, in µs per complete operation:

| Operation | Current | Dense maps | Dense maps and empty-constraint guards | Plus trusted later projection |
| --- | ---: | ---: | ---: | ---: |
| `try_from_entries` | 61.9 | 38.2 | 28.4 | — |
| `lhs` | 26.9 | 15.7 | 10.9 | 7.6 |
| `to_reaction` | 42.7 | 21.1 | 11.3 | 6.9 |

The dense-map-only variant was faster but did not reduce allocation calls:
80-atom construction used 199 rather than 187, while requesting slightly fewer
bytes. Combining it with the empty-constraint guards reduced construction to
99 allocations; adding trusted later projection did not change construction.
The trusted path reduced `lhs` from 43 to 37 allocations after the other two
changes. These counts include reallocations, not peak live memory.

Three top-level constraint spans, including atom references on both sides,
prevent the empty-constraint shortcut. For the 80-atom unchanged span,
construction changed from 61.3 to 38.5 µs with dense maps; allocation calls
rose from 195 to 205. With dense maps and a private trusted projection, `lhs`
changed from 26.5 to 11.8 µs and `to_reaction` from 42.6 to 16.6 µs. The
interleaved union fixture changed from 57.8 to 35.1 µs for construction and
48.1 to 20.6 µs for `to_reaction`. No production source or public API changed
in these comparisons.

Dense projection maps and the empty-constraint guards have direct,
representation-preserving explanations and improve complete operations.
The additional trusted projection saved about 3–4 µs per later 80-atom
projection after those changes, but requires a private molecule-construction
path that relies on every ReactionSpan publisher and transformation preserving
both side invariants. Keep that boundary separate from the initial local
optimization; do not bypass the molecule checks at span publication. The
`Reaction::new` check in `to_reaction` also remains until a preservation proof
and a measured need justify changing its internal path.

### Verification

The earlier cost investigation ran 97 focused `reaction::` properties with one
pre-existing ignored case and forty exact span-constructor unit cases. This
adversarial cycle used exact public-API reproducers for R1–R5 and sampled the
existing generators with repository-locked proptest 1.9.0; it did not rerun the
full property target. Existing tests cover other invalid references, overlay
incidence, span side integrity, and modified stereo kinds, but not R1–R5.
No tests or production code changed.

## Proposed changes and next action

Fix the confirmed kindless-symmetry panic in the consumer, preserving the admitted
domain. The bounded optimization then preserves accepted inputs, typed
diagnostics, and publication boundaries: generic reference predicates, borrowed
stereo frames, error-only diagnostic collection, one bounded uniqueness scan,
separate site exclusion, combined participant-uniqueness/identity checking,
and family-local table lifetimes. Use temporary bitmap keys and aromatic
membership when the molecule has at most 128 atoms; sorted keys and an aromatic
membership HashSet cover larger namespaces.
Remove already-implied endpoint/closed-frame checks. Include the precise
naming/visibility migrations, the two documentation
corrections, and focused expression and symmetry regressions. No new public
Molecule symbols are needed.

Keep the two-pass constraint walk, constructor preflight, and remapping result
gate. Reserve the dative/multicenter identity and stereo-atom site HashSets from
their known relation counts. Detect duplicate stereo-bond sites through the
existing bond-site incidence index, checking for an earlier relation at the site.
Replace the localized-bond HashSet with adjacency lookup and the noncovalent-bond
HashSet with inline sorted keys; their first-error selection may change. Other
alternatives remain separately identified above; they do not block the local
improvements. A later implementation should measure the combined effect on the
actual gate without conflating it with storage construction or editor rebuilding.
These decisions and evidence are durable here;
review scratch files are disposable.

For F2/F4/F15, choose the gate-local bitmap path when the molecule has at most
128 atoms. Valid atom references then fit the two words; retain reference checks
before indexing them. Above 128 atoms, use the unrestricted sorted-key path for
dative/multicenter identity and aromatic row uniqueness, with an aromatic global
HashSet. Select the path once per family, so the compact bitmap key is not widened
by a per-key tagged representation. Keep the outer dative/multicenter identity
HashSets with reserved capacity and existing error order. Do not add an inline
identity-table branch, a stereo-site bitmap branch, or a dynamic aromatic bitmap
above 128 atoms on the basis of isolated kernel timings. No permanent
relation-storage key is part of this local change.

For Reaction, close R2–R5 at `try_new`/`new` by establishing the domain of each
frame-bearing electron or stereo payload against its owning or explicit local
frame. Keep the R10 removal-incidence check in the constructor. Close R1 by
checking a removed bond's endpoints against a same-reaction `Add` when the id
does not exist on the lhs. Build one borrowed `Add` map keyed by `Entity` for
duplicate detection, reference availability, and created-entity source lookup;
read existing source frames from lhs views. Preserve the current check ordering
and error precedence. Public `Deltas::normalize` must independently check
the structured incidence and old attributes of a same-id `Add`/`Remove` pair
before canceling it. Transport a compatible removal into the Add frame, apply
intervening modifications, and compare its old attributes using
`normalized_eq`; report `Contradiction` on disagreement. Do not change the
normal form for a matching pair. Keep other delta continuity, product
materializability, satisfiability, DPO, and chemistry with their existing
first consumers. Report an `Add` ID collision as
`ReactionIntegrityError::DuplicateReference { entity }`, leaving
`InvalidReference` for unavailable IDs. Use the direct
`ReactionIntegrityError::ElectronCountLengthMismatch` variant with
`{ entity, participants, electron_counts }`; update the integrity guide with
the revised Reaction error inventory. Exact constructor failures and
successful identity frame transport should exercise the corrected gate.
Replace the broad stereo wrapper with the eight direct variants above while
reusing local checks; exact error assertions should cover the affected
constructor and DSL paths.

For R11, remove Reaction's `StereoKindModified` constructor error and its
ModifyField old/new comparison. Preserve the application precondition and the
ReactionSpan kind-continuity check. Exact cases should accept a locally
well-formed different-kind ModifyField in `Reaction::try_new`, reject its
application as `InconsistentReaction`, and reject two-sided conversion as
`Contradiction`. The same checks must still distinguish out-of-domain payloads
under R3/R4 from valid but unexecutable changes.

For R1 and delta cancellation, cover localized bonds and all six overlay
kinds, lhs and forward-`Add` sources, compatible local permutations, and
mismatches rejected by their respective constructor or normalization boundary.
Cover bare `Deltas::normalize` with matching and mismatching created
`Add`/`Remove` pairs, including an intervening modification. Update the
integrity guide for the completed constructor check; do not move R10's failure
boundary or change the corresponding Rust/Python constructor error contract.

For ReactionSpan, correct R6/R7 rustdoc, use dense maps for projected union IDs,
and omit correspondences when there are no selected top-level constraints or
constraint deltas to transport. Keep both side-molecule checks at span
publication and preserve current constructor/error semantics. Compare complete
operations with empty and nonempty constraints and interleaved union IDs.
Trusted later projections and micro-optimizing removal-frame comparison remain
separate decisions.

## Staged implementation plan

Each subitem includes focused exact cases at its public boundary and leaves the
affected crate green. Internal moves are green only after their callers move in
the same subitem. Public error changes are breaking subitems: migrate every
Rust and Python caller, test, and rustdoc reference before ending that subitem.
Run focused tests and checks during stages, including a Python binding compile
when a public error changes; reserve full workspace tests, feature suites,
Python tests, and the Rust 1.87 gate for S6. Benchmark inputs are constructed
outside the timed region. Record before/after operation times and allocation
counts here so that scratch experiments may be deleted.

### S0 — Baselines and reference-check locality

- **S0a — umol-graph-ir/benches/integrity.rs; additive (green) [dep: none].**
  Register controlled Molecule construction, Reaction construction, and
  ReactionSpan construction/projection cases with empty and nonempty constraints
  and with zero, one, and several removals. Run the existing exact publication
  cases to establish that benchmark inputs satisfy their intended boundaries;
  record baseline time and allocation protocol before changing the gates.
  **Complete 2026-09-23.**
- **S0b — molecule::integrity, molecule.rs, reaction_span.rs; internal move
  (green) [dep: S0a].** Move the entry and recursive constraint-reference
  traversals from molecule.rs into molecule::integrity. Rename the five
  validate_* functions to check_*, use a generic borrowed predicate instead
  of dynamic dispatch (F1), and expose only the two traversals needed by
  constructors outside that module. Rename ReactionSpan's related
  validate_reaction_span_entries function to check_reaction_span_entries.
  Keep reference-first error order and check the existing exact
  missing-reference cases, including nested constraints and ReactionSpan.
  **Complete 2026-09-23.**
- **S0c — MoleculeIntegrityError and callers; breaking, restored green
  [dep: S0b].** Apply the seven approved variant renames in this record,
  including DuplicateAtom. Migrate constructors, editor errors, DSL and
  Python bindings, rustdoc, and exact error assertions together. Check
  compilation of affected crates and the exact error cases.
  **Complete 2026-09-23.**

#### S0a baseline (2026-09-23)

The permanent integrity benchmark has 36 cases. Molecule inputs are 12- and
80-atom chains plus an 80-atom case with 20 dative and 10 of each other overlay
kind. Reaction inputs use that overlay-rich lhs with 0, 1, or 20 matching dative
removals. ReactionSpan inputs use unchanged 12/80-atom chains or an 80-atom
interleaved union with 20 added and 20 removed atoms. Every family has empty
and one-constraint cases. Benchmark test mode successfully published every
input and projection. Existing focused constructor cases also passed: 116
Molecule/ReactionSpan try_from_entries cases, 31 Reaction try_new cases, and
six external ReactionSpan view cases. The bench passed strict Clippy, and
workspace formatting passed.

Selected baselines on arm64, rustc 1.96.0:

| Complete operation | Time, µs | Allocations | Requested bytes |
| --- | ---: | ---: | ---: |
| Molecule 12-atom chain, no constraint | 1.20 | 23 | 2,132 |
| Molecule 80-atom chain, no constraint | 6.73 | 27 | 14,732 |
| Molecule 80-atom overlays, no constraint | 26.10 | 330 | 40,316 |
| Molecule 80-atom overlays, one constraint | 26.39 | 330 | 40,316 |
| Reaction, no constraint, 0 removals | 6.02 | 66 | 11,220 |
| Reaction, no constraint, 1 removal | 6.00 | 68 | 11,228 |
| Reaction, no constraint, 20 removals | 7.17 | 106 | 11,380 |
| Reaction, one constraint, 20 removals | 7.53 | 106 | 11,380 |
| ReactionSpan 80-atom unchanged construction, no constraint | 60.04 | 187 | 214,966 |
| ReactionSpan 80-atom unchanged construction, one constraint | 61.68 | 189 | 215,414 |
| ReactionSpan 80-atom interleaved construction, no constraint | 55.42 | 183 | 127,158 |
| ReactionSpan 80-atom interleaved construction, one constraint | 55.55 | 185 | 127,606 |
| ReactionSpan 80-atom interleaved lhs projection, no constraint | 26.26 | 85 | 54,555 |
| ReactionSpan 80-atom interleaved to_reaction, no constraint | 49.13 | 145 | 112,654 |
| ReactionSpan 80-atom interleaved to_reaction, one constraint | 49.03 | 146 | 112,878 |

Times are Criterion point estimates from a 0.5-second warmup, 1-second
measurement, and 20 samples per case; input cloning is outside the timed
operation. Two projection outliers were rerun in isolation and are omitted from
this selected table. Allocation measurements used a separate release-mode
System counting allocator with the workspace lockfile: warm once, build or clone
the input before resetting counters, then count allocation and reallocation
calls during one complete operation. Values are medians of seven calls.
Requested bytes sum successful allocation and reallocation sizes; they are
neither peak live memory nor process RSS. Timing and allocation binaries are
separate, and scratch probe files are disposable. Reuse this protocol and these
fixtures for S2, S4, and S5 comparisons; compare complete operations, not only
isolated kernels.

#### S0b result (2026-09-23)

The five reference traversals now live in molecule::integrity; only
check_entry_references and check_constraint_references are crate-private.
ReactionSpan's entry preflight uses the check_* name. The borrowed predicate is
statically dispatched in both entry and closed-molecule reference checks.
The 116 focused try_from_entries cases and strict crate Clippy passed;
formatting and diff checks passed.

Complete-operation timings did not show a consistent material gain: the
80-atom overlay-rich Molecule case with one constraint changed from 26.39 to
25.84 µs, while 80-atom interleaved ReactionSpan construction with one
constraint stayed at 55.55 versus 55.63 µs. The no-constraint span case moved
from 55.42 to 56.12 µs. F1 remains a local removal of avoidable dispatch, not
an end-to-end speedup claim.

#### S0c result (2026-09-23)

All seven MoleculeIntegrityError variants now use the approved names in the
integrity gate, its Rust consumers, exact error assertions, and the living
integrity inventory. Their fields, display messages, and rejection paths did
not change. The Python binding converts MoleculeIntegrityError generically, so
it had no variant-specific code to migrate. Graph-IR's 6,816 active library
tests and IO's 4,137 library tests passed; strict Clippy passed for graph-IR,
IO, and umol-py with Python 3.13 active. Graph-IR rustdoc passed with warnings
denied.

### S1 — Stereo integrity and admitted inputs

- **S1a — stereo::integrity, molecule::integrity, reaction::integrity; breaking,
  restored green [dep: S0c].** Move only shared local stereo-frame, kind,
  configuration, and constraint checks from molecule::integrity into
  stereo::integrity. Give that module a crate-internal local error mapped into
  each aggregate's own public error variants; replace Reaction's broad
  StereoIntegrityError wrapper with its eight direct variants and migrate
  callers. The aggregate gates stay in their respective integrity modules.
  Check exact existing stereo errors and add the accepted and malformed
  nested-term cases (F9). Add no public helper or checker type.
  **Complete 2026-09-23.**
- **S1b — symmetry.rs and coloring.rs; green [dep: S1a].** Handle a published
  kindless stereo frame without calling the determined-kind accessor (F14),
  including the full-coloring path. Keep it admitted by Molecule construction
  and add the public graph_symmetry regression; retain unresolved and
  nonliteral behavior.
  **Complete 2026-09-23.**

#### S1a result (2026-09-23)

The local frame, kind, configuration, and constraint checks now live in
stereo::integrity. Molecule and Reaction map its eight failures into separate
public error enums; Reaction no longer wraps MoleculeIntegrityError for stereo
delta payloads. Reference and site-incidence checks remain in the aggregate
gates. New Molecule::try_from_entries cases reject malformed nested terms and
preserve an accepted raw term. The 6,820 active graph-IR library tests, strict
crate Clippy, rustdoc with warnings denied, and the Python binding crate check
passed.

#### S1b result (2026-09-23)

Graph symmetry now treats a stereo atom or bond without a determined kind as
having no orientation contribution. Full constitution coloring hashes the
optional kind instead of calling an asserted accessor. Public Molecule
construction still admits both frames; exact atom- and bond-frame regressions
check their resulting orbits and achirality under entity-only and full coloring.
The 6,822 active graph-IR library tests and strict crate Clippy passed.

### S2 — Delta and Reaction correctness

- **S2a — delta.rs; green [dep: S1a].** In the generic and both stereo
  created-entity folds, cancel Add followed by Remove only after structured
  incidence agrees and the removal's old attributes, transported into the Add
  frame, are normalized_eq to the state after intervening modifications.
  Return Contradiction otherwise. Exact cases cover atoms, bonds, six
  overlays, reordered compatible frames, changed sites/factors, and changed
  old values; run the existing delta normal-form properties.
  **Complete 2026-09-23.**
- **S2b — reaction::integrity; breaking, restored green [dep: S2a].** Replace
  the created-ID set and copied source-frame map with one borrowed Add map
  keyed by Entity. Read existing frames through lhs views; preserve pass order.
  Report duplicate additions as DuplicateReference while retaining
  InvalidReference for unavailable IDs; migrate error consumers. Complete R1
  by checking a created bond's recorded removal endpoints against its Add.
  Keep R10 at construction. Check lhs and forward-Add sources, duplicate
  additions, exact IncidenceMismatch, and the constructor benchmark.
  **Complete 2026-09-23.** One borrowed Add map now serves uniqueness,
  availability, and added-source lookup; lhs frames are read only for removals.
  The constructor rejects added-ID collisions as DuplicateReference and checks
  created-bond removal endpoints. Exact cases cover lhs and forward-Add
  sources, reordered frames, mismatches, and first-pass error precedence.
  On the same Criterion 80-atom constructor fixture, the committed pre-S2b
  code measured 5.86, 5.97, and 8.32 µs for zero, one, and 20 removals
  without constraints. S2b measured 0.025, 0.062, and 1.07 µs respectively.
  Inputs were cloned outside the timed body in both runs. The 20-removal
  baseline had outliers; the constrained fixture measured 7.30 versus
  1.07 µs for that case. Eliminating the eager lhs frame copy improves the
  sparse cases and does not trade that gain for slower dense removals. The
  earlier scratch probe used a different harness and is not a numeric
  baseline for this comparison.
  Crate tests, Reaction property tests, strict crate Clippy, and rustdoc passed.
- **S2c — reaction::integrity; breaking, restored green [dep: S2b].** Add
  ElectronCountLengthMismatch for R2 and check literal aromatic/multicenter
  counts in Add, Remove, and both ModifyField sides against their frame.
  Migrate error consumers and check exact lengths and errors, retaining
  undetermined values.
  **Complete 2026-09-23.** The constructor now checks local Add/Remove frames and
  the lhs or same-reaction Add frame for ModifyField. Exact cases, crate tests,
  Reaction properties, strict crate Clippy, rustdoc, and the Python binding
  check passed. The Python error mapping already handles new integrity variants.
- **S2d — reaction::integrity; green [dep: S1a, S2b].** Use the shared local
  stereo integrity checks for R3–R5: Remove configurations, both ModifyField
  configurations, entity constraint changes, and top-level ConstraintDelta
  positions/actions against the owning frame. Exact tests distinguish local
  domain failures from deferred old-value and product failures; run Reaction
  identity-frame and publication properties.
  **Complete 2026-09-23.** The constructor checks these payloads against the
  lhs or a same-reaction Add frame. Exact cases cover each payload and preserve
  deferred old-value and product failures. An application-precondition generator
  now uses a locally valid but incompatible stereo kind. Crate tests, Reaction
  properties, strict crate Clippy, and rustdoc passed.
- **S2e — reaction::integrity and callers; breaking, restored green
  [dep: S2d].** Remove the constructor's determined-old/new stereo-kind
  comparison and ReactionIntegrityError::StereoKindModified (R11). Keep
  application and ReactionSpan checks. Exact cases show constructor
  acceptance followed by the existing application/span failures.
  **Complete 2026-09-23.** Reaction construction accepts individually valid
  old/new configurations of different kinds. Exact atom and bond cases show
  `try_new` acceptance, `check_preconditions` returning `InconsistentReaction`,
  and `to_reaction_span` returning `Contradiction`. The ReactionSpan check remains.
  Crate tests, Reaction properties, the Python binding check, strict crate
  Clippy, and rustdoc passed.

#### S2a result (2026-09-23)

Created-entity cancellation now compares the removal's incidence and old
attributes with the state after intervening changes. Overlay and stereo
removals transport their old attributes into the Add frame; incompatible frames
and values return Contradiction. The 27 focused cases, 127 delta unit cases,
delta normal-form property, and strict graph-IR Clippy passed.

### S3 — Documentation boundaries

- **S3a — editor.rs, stereo.rs, reaction_span.rs, integrity guide; green
  [dep: S1a, S2e].** Correct donor-order and coset-tier claims (F8, F13),
  ReactionSpan anchoring and roundtrip claims (R6, R7), and the completed
  Molecule/Reaction error inventories. Run rustdoc with warnings denied and
  the existing order/roundtrip cases; do not change the represented values.
  **Complete 2026-09-23.** Corrected the donor-order and coset-range comments,
  the raw-versus-lhs-anchored span description, and the exact roundtrip condition.
  The Molecule and Reaction integrity-guide inventories already match their
  current error enums. Focused dative-order, coset, and span-roundtrip cases,
  ReactionSpan properties, and rustdoc with warnings denied passed.

### S4 — Molecule gate allocations and repeated checks

- **S4a — stereo::integrity and molecule::integrity; green [dep: S1a, S0a].**
  Borrow stored ligand slices (F3), use the bounded prefix uniqueness scan
  (F7), remove the implied second stereo-bond atom set (F5), and use a site
  exclusion scan for stereo atoms (F6). Preserve degree and first-error
  order; compare exact malformed-frame cases and complete gate timings.
  **Complete 2026-09-23.** The Molecule gate borrows both stored ligand frames.
  The shared stereo check scans the bounded frame prefix; stereo atoms scan
  only for an Atom-kind ligand equal to the site, and stereo bonds no longer
  repeat the implied atom-uniqueness check. Existing exact malformed-frame and
  first-error cases passed, as did all 6,897 active graph-IR library tests,
  strict crate Clippy, rustdoc with warnings denied, and formatting.

  In the 80-atom overlay-rich complete Molecule construction benchmark, the
  no-constraint case changed from 25.369 to 22.953 µs and the one-constraint
  case from 25.480 to 22.916 µs (Criterion point estimates, 20 samples, same
  run protocol). The S0a allocation probe changed from 330 calls and 40,316
  requested bytes to 270 calls and 37,876 bytes in both cases. The unchanged
  80-atom chain allocation count remained 27 calls and 14,732 bytes.
- **S4b — molecule::integrity; green [dep: S0b, S0a].** Remove the assembled
  graph's redundant endpoint recheck while retaining raw-entry preflight
  (F10). Use sorted adjacency for localized-bond duplicates and inline
  sorted keys for noncovalent duplicates. Check parallel/self-loop errors
  and compare construction time and allocations.
- **S4c — molecule::integrity; green [dep: S0a, S0c].** Replace dative and
  multicenter per-row HashSets/trees with exact two-word keys for molecules of
  at most 128 atoms and sorted keys above that bound (F2, F4). Reserve the
  retained family identity tables. Check duplicate participants, identical relations,
  empty rows, threshold edges, and first-error precedence; compare complete
  gate costs.
- **S4d — molecule::integrity; green [dep: S0a, S4c].** Use row and global
  bitmaps for molecules of at most 128 atoms and a sorted row with
  global HashSet above that bound (F15). Check within-row duplicates before
  overlap and electron-count errors; cover the 128/129 boundary and compare
  complete gate costs.
- **S4e — molecule::integrity; green [dep: S4a, S4b, S4c, S4d].** Scope each
  family table to its loop (F16), reserve the stereo-atom site table, and use
  existing site incidence for stereo-bond duplicates. Preserve exact
  duplicate-site behavior; measure allocation counts and peak-live memory
  separately from operation time.
- **S4f — symmetry.rs and canonicalize.rs; green [dep: S1b, S4a].** Remove
  the closed-frame arity and uniqueness rechecks (F11, F12), retaining
  nonliteral, action, and normalization failure paths. Run focused
  symmetry/canonicalization cases and their existing public properties.

### S5 — ReactionSpan allocations

- **S5a — reaction_span.rs; green [dep: S0a, S3a].** Skip correspondence
  construction only when no selected top-level constraints or constraint
  deltas need transport (R9). Compare exact empty/nonempty constraint
  projections and complete operation costs.
- **S5b — reaction_span.rs; green [dep: S5a].** Replace dense union-ID
  HashMaps with indexed vectors while preserving present-prefix and absent-ID
  translation (R8). Check interleaved union ordering, dangling references,
  side projections, and exact roundtrips; compare time and allocations for
  complete construction and projection.

### S6 — Final gate and closeout

- **S6a — workspace and discussion record; green [dep: S2e, S3a, S4e, S4f, S5b].**
  Run formatting, focused and full workspace tests, feature-gated property
  suites, strict Clippy, rustdoc with warnings denied, and Python tests using
  the repository Python 3.13 environment. Run the Rust 1.87 CI gate once.
  Re-run complete-operation benchmarks only after the final changes, record
  decisions and remaining limits here, review the full diff against scope,
  then synchronize this document and the status index.

The correctness work runs through S0–S3. S4 and S5 follow it; both are part of
this plan and must be complete before S6 closes.
