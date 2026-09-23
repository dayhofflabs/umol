# 230 — graph-IR delta module review cycle

Status: Proposed
Date: 2026-09-23
Relates: [184](184-deltas-and-edits-2026-08-04.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[209](209-normalization-canonical-semantics-2026-08-25.md),
[229](229-aggregate-integrity-review-2026-09-22.md),
[134](134-reaction-application-overlays-2026-06-26.md),
[168](168-api-hygiene-2026-07-27.md),
[227](227-repository-structure-hygiene-2026-09-10.md),
[206](206-umol-perm-review-2026-08-21.md),
[code reviews guide](../docs/development/code-reviews.md),
[data type contracts guide](../docs/development/data-types.md),
[integrity guide](../docs/development/integrity.md),
[nomenclature guide](../docs/development/nomenclature.md),
[property tests guide](../docs/development/property-tests.md)

## Purpose

Review-cycle record for the graph-IR delta module, produced under the `review-cycle` skill and
the code-reviews guide. Five area reviews (construction/integrity/fallibility, nomenclature,
module structure and visibility, tests and generators, documentation) ran against pinned commit
`5fbdf5399`, followed by one refutation pass over the pooled findings. This document records the
surviving findings with both arguments, the refuted findings with their dismissing citations, the
migrations, the consolidated open items, and the proposed design. It contains no staged
implementation plan; accepted findings enter the ordinary design and planning lifecycle after
triage.

## Review result

Forty findings were pushed (construction 5, nomenclature 9, structure 3, tests 16, documentation
7). After refutation: 18 confirmed (17 distinct defects; one documentation finding duplicates a
nomenclature finding), 11 reduced, 2 open questions, 9 refuted, of which 4 were reclassified as
migrations. Every reproduction supporting a surviving finding was re-run by the refutation agent
and reproduced with the recorded output.

**Check placement is fully compliant.** Every fallible path in the module body classifies as
tier-1 integrity or first-requiring-operation validation; there are no defensive re-validations,
no construction-time enforcement of unrequired properties, no inherited fallibility, and no
panic constructs. The S2a cancellation logic landed at the pinned commit is compliant, and the
direction of its frame transport was confirmed with non-involutive probes.

**The defects concentrate in `Normalize for Deltas`.** Three confirmed semantic defects: a
cancelled created atom leaves a surviving created bond or overlay referencing it, so
`Reaction::normalize` and `Reaction::reframe` publish a `Reaction` that fails its own integrity
gate and `to_reaction_span`/`canonicalize` panic on the original; the normal form never reduces
carried payload values, so `normalized_eq` and `framed_eq` on `Deltas` and `Reaction` distinguish
encodings of one value and `reframe` disagrees with `canonicalize` in the canonical frame; and the
generic folds accept unsatisfiable payloads that the stereo folds reject. A fourth, reduced:
mixed stereo `kind` stamps on one entity make the normal form depend on input order.

**The S2a change has no generative evidence.** After S2a every generated same-id `Add`/`Remove`
pair in the delta property suite normalizes to `Contradiction`, and the reaction generators never
produce such a pair, so the cancellation branch is reached by exact cases only. The bare-`Deltas`
normal-form property samples no overlay or stereo kinds.

**Documentation carries six inaccurate statements**, including a field doc that says
normalization does not read a value that normalization reads, derives, and rewrites. Missing
contract sections are a migration, not a defect. **Identifiers are compliant**; the nomenclature
defects are in rustdoc prose and test names. **Module size, inline tests, and `pub(crate)`
markers are migrations** tracked by docs 168 and 227; a complete split proposal with its
visibility map is recorded below. One structural design question is open: the uneven `id`
surface across the eight `*Delta` enums.

## Scope and method

Target: `umol-graph-ir/src/ir/delta.rs` (5157 lines: the eight entity `*Delta` enums,
`ConstraintDelta`, `Delta`, `Deltas` with `Normalize for Deltas`, `EntitySpan`, `ConstraintSpan`,
the `EntityOp`/`EntityFold` fold machinery, the `FrameTransport` impls, the `apply_*_change`
family, the stereo folds, and the inline unit-test module), `umol-graph-ir/tests/property/delta.rs`
with the generators it draws from `tests/property/strategies.rs`, and the module's re-exports in
`ir.rs`. Consumers (`reaction.rs`, `reaction/integrity.rs`, `reaction_span.rs`, `compose.rs`,
`canonicalize.rs`, `molecule/transact.rs`, `umol-py/src/delta.rs`) were traced for check
ownership, synonym checks, and defense arguments, not reviewed. The pinned commit is doc 229's
stage S2a and is in scope.

Each review agent and the refutation agent worked in its own detached worktree at the pinned
commit; no build, test, or edit touched the primary checkout, no tracked file changed, and the
worktrees were removed after this document was written. Every agent read the code-reviews guide,
the review-cycle skill, the living guides and skills for its area, the status index, docs 184,
214, 211, 209 (top notice respected), 229, 134, 135, 212, 215, 176, 161, 156, and the
whitepaper. Every finding carries a violation argument with its citation and a defense argument
with a deferral check naming the documents consulted. The refutation agent verified every
citation and factual premise at the pinned commit, re-ran the reproductions, and issued graded
verdicts; ties favor the original.

Evidence at the pinned commit: `cargo test -p umol-graph-ir --lib delta` 325 passed
(`ir::delta::tests` 127 cases, including the 27 focused S2a cases doc 229 records);
`--features proptest --test property delta` 30 passed. Line coverage of `src/ir/delta.rs` over
the full lib and property runs (`cargo llvm-cov` 0.9.1): 91.4 % lines, 87.1 % regions, 94.0 %
functions (refutation re-measurement 91.3 % / 199 missed lines), above the 80 % baseline.
Convention dating by read-only git: the code-reviews guide has one commit (2026-08-31); the
`# Semantic properties` policy entered `property-tests.md` 2026-08-06 and the contract-section
paragraph entered `data-types.md` 2026-08-31; the AGENTS.md visibility, `module.rs`, and comment
rules landed 2026-09-03; the test-writing skill's no-behavior-in-names and no-helper-constructor
rules existed from 2026-05-09 in the `.claude/skills` catalog; `delta.rs` dates from 2026-06-25
with 85 commits.

## Generator statistics

Sampled with `TestRunner::deterministic()`, 10,000 instances per strategy; the refutation pass
re-sampled `deltas_strategy`, `comprehensive_reaction_strategy`, and `atom_form_strategy` and
matched these figures. Instrumentation was discarded with the worktrees.

| Strategy | Sampled distribution |
| --- | --- |
| `deltas_strategy` | length 0..7 uniform, empty 12.3 %; kinds Atom 34.0 % / Bond 35.0 % / Constraint 34.0 %, six overlay kinds and both stereo kinds 0 %; variants ≈ 25 % each; ≥1 same-id group 33.4 %; same-id Add+Remove pair 5.74 % (574), of which `normalize` Ok 0, Err 574, cancellation reached 0; `normalize` Ok 71.5 % / `Contradiction` 28.5 %; of Ok, empty 18.9 %; Ok samples that fold a same-id group 6.3 % of all samples (81 % of grouped samples are `Contradiction`); `ModifyField` identity 15.0 %; `ModifyConstraint` identity 25.8 %, `(None, None)` 24.9 % |
| `atom_delta_strategy`, `bond_delta_strategy` | variants ≈ 25 % each; `ModifyField` identity 13.9 % / 14.9 %; `ModifyConstraint` `(None, None)` 24.9 % / 26.3 % |
| `constraint_delta_strategy` | six values (Add/Remove × ChargeSum 0/1/2), 16–17 % each |
| atom, bond, dative, aromatic, multicenter, noncovalent form strategies | non-normalized 0.00 %; `normalize` Err 0.00 %; pairs `normalized_eq` 0 / 0 / 2.9 / 0 / 0.02 / 11.2 % (each equal to the `==` rate); `==` vs `normalized_eq` disagreement 0.00 % |
| `stereo_atom_form_strategy` | kinds Tetrahedral 24.6 / SquarePlanar 25.1 / TrigonalBipyramidal 25.4 / Octahedral 24.9 %; cosets `Lit` 49.6 % / `Undetermined` 50.4 % (`LitSet`, `Term` 0); constraint counts 0–3 balanced; pairs `normalized_eq` 0.47 % |
| `stereo_bond_form_strategy` | CisTrans 100 %; cosets 49.9 / 50.1 %; pairs `normalized_eq` 2.3 % |
| `comprehensive_reaction_strategy` (as consumed by `test_delta_inverse`) | 67,193 deltas; kinds Atom 50.7 / Bond 20.7 / Noncovalent 7.8 / Dative 6.5 / StereoBond 4.6 / Aromatic 4.5 / Multicenter 4.5 / StereoAtom 0.85 / Constraint 0 %; `ModifyConstraint` only for DativeBond (324; re-sample 291); Noncovalent `ModifyField` 0; stereo Add/Remove kinds Tetrahedral 391 / SquarePlanar 103 / TrigonalBipyramidal 6, degrees 4 and 5 only; zero-delta reactions 3.3 %; same-id groups in 1.0 % of reactions (all dative); same-id Add+Remove pairs 0; `Deltas::normalize` Err 0 % |

Direction probe: all seven S2a reordered-frame cases use a single transposition, for which
`between(remove, add)` equals `between(add, remove)`; a 3-cycle on aromatic or multicenter
participants distinguishes the directions and confirms the shipped direction is correct.

## Confirmed and reduced findings

Each finding records the claim, the strongest case for the original (the defense), and the
refutation verdict. Line references are at commit `5fbdf5399`. Original finding ids are given in
parentheses for audit against the refutation record.

### Normalization semantics

**F1 — `Deltas::normalize` cancels a created atom that a surviving created delta still
references** (`delta.rs:1934-1938`, `3339-3435`; confirmed; FA1). `[Atom Add 1, Bond Add 0 on
[0, 1], Atom Remove 1]` normalizes to `Ok([Bond Add 0 on [0, 1]])`. `Reaction::try_new` accepts
the input; `Reaction::normalize` and `Reaction::reframe` publish a `Reaction` whose parts fail
`Reaction::try_new` with `InvalidReference { Atom(1) }`; `to_reaction_span` and `canonicalize` on
the original panic at `reaction_span.rs:2369` ("no entry found for key"); `check_preconditions`
passes and `apply` reports `CorrespondenceMismatch`. The aromatic twin behaves the same. Claim:
data-types guide §"Closed containers and the minimum eager contract" (every publisher invokes the
gate or establishes preservation by construction); integrity guide lines 57–59 and the
`InvalidReference` row ("a missing id would panic or select no source frame"); the module's own
guard already rejects a created entity on a net-removed lhs atom (`:3358`). Defense: doc 214
lines 235–236 state the cancellation rule without a reference condition; doc 229 §S2a specifies
incidence and old-attribute agreement only, which the code implements; integrity guide lines 111
and 177 and doc 229 lines 1573–1576 keep mutual consistency and materializability with the first
consumer. Verdict: the deferral covers the check, not publication of a value that violates the
closed-container contract; the first consumer does not report the deferred failure, it panics,
which the data-types guide permits only where the producer establishes the property. Owner:
`Normalize for Deltas`.

**F2 — the delta normal form leaves carried payload values unreduced** (`delta.rs:1943-1946`,
`1953-2016`, `2992-2999`, `3082-3093`, `3147-3154`, `3237-3248`; only `ConstraintDelta` payloads
pass through `normalize()?` at `3437-3446`; confirmed; FA2). Charge `LitSet({2})` against `Lit(2)`:
`normalize(a) != normalize(b)`, `a.normalized_eq(&b)` false; for reactions `ra`/`rb` differing only
in that encoding, `normalized_eq` and `framed_eq` are false while `canonical_eq` is true;
`reframe(ra).deltas` carries `LitSet({2})`, `canonicalize(ra).deltas` carries `Lit(2)`, so
`reframe(ra) != canonicalize(ra)` while `reframe(rb) == canonicalize(ra)`. `Molecule::normalize`
(`molecule.rs:1954-1970`) reduces every form. Claim: doc 214 line 496 ("`Normalize for Reaction`
reduces the lhs and deltas in their stored id and participant frames"); nomenclature guide
§Normalize (folds value expressions and normalizes set representations for "fixed-frame
transformation values such as `Deltas`"); §Equality ladder (`normalized_eq` is equality of normal
forms); data-types guide §"Raw representation and semantic normal form" lines 306–310 (an
operation deriving a delta normal form documents the form and its roundtrip under the semantic
equality). Defense: doc 211 lines 603–609 describe the reduction as the id-keyed fold; the fold's
decisions are encoding-tolerant (`normalized_eq` throughout), so consistency is unaffected;
`canonical_eq` answers correctly; doc 209's "retains its existing fixed-frame normal form" lies in
its superseded portion. Verdict: doc 214 is the latest governing statement and requires
reduction; doc 211's fold description is not a decision to leave payloads unreduced. Owner: the
rebuild step of `fold_created`, `fold_preserved`, `fold_stereo_atom_group`, and
`fold_stereo_bond_group`.

**F3 — the generic folds accept unsatisfiable payload values that the stereo folds reject**
(`apply_field` at `delta.rs:1560-1567` stores `new` unchecked; `fuse_field` at `1629-1637` fuses
through an unsatisfiable link because `normalized_eq` equates two contradictions; stereo paths use
`normalize()?` at `2928`, `2933`, `3063`, `3218`; confirmed; FA3). A charge of `LitSet({})` in an
`Add`, a `ModifyField`, or a chain `0 → ⊥ → 1` normalizes to `Ok` on the generic path (the chain
to `0 → 1`); the structurally identical stereo cases return `Contradiction`. Claim: `Normalize`
trait doc (`traits.rs:149-152`, `161-165`: `Err` means unsatisfiable); nomenclature guide
§Normalize line 1214; `Molecule::normalize` rejects the same atom, so `Reaction::normalize` accepts
in its deltas what it rejects in its lhs. Defense: the settled `normalized_eq` rule equates
unsatisfiable values (`traits.rs:168`; doc 214 lines 258–260), so fusing through `⊥` follows from
it; reaction integrity keeps satisfiability lazy (integrity guide lines 144–147); the stereo
strictness may be incidental. Verdict: the contract is on `Normalize` itself and the two fold
families implement one trait with opposite outcomes. Resolved by the F2 owner (payloads through
`normalize()?`), after which the stereo folds need no separate `normalize()?` comparison.

**F4 — mixed stereo `kind` stamps make the normal form depend on input order** (`delta.rs:3022`,
`3031`, `3089`; twins `3177`, `3186`, `3244`; reduced; FA4). Two `ModifyConstraint` deltas on
different constraint keys of one stereo entity with `kind` stamps `Tetrahedral` and
`SquarePlanar`: one input order stamps both outputs `Tetrahedral`, the other stamps both
`SquarePlanar`; the same multiset has two normal forms under derived `==` and `normalized_eq` is
false. Claim: `delta.rs:3308-3310` ("each entity's fold is deterministic over input order, so the
result is a unique normal form"); doc 214 line 497; `traits.rs:171-173`. Defense: doc 134 §6 I6e
(line 410) records that the fold threads `kind` through the constraint group; the in-crate
producers stamp one kind per entity, so mixed stamps arise only from hand-assembled or DSL input.
Verdict: reduced. Survives: for differing stamps the selection is first-seen in input order and
is written onto every output. Did not survive: the `None → Some(..)` derivation as a defect,
which is the recorded design (the false field doc is F5). Owner: `fold_stereo_atom_group` and
`fold_stereo_bond_group`.

### Documentation accuracy

**F5 — the `kind` field doc says normalize does not read it; normalize reads, derives, and
rewrites it** (`delta.rs:776-779`, cross-reference `:902`; confirmed; FE1, also raised by areas A
and B). Text: "Not read by apply/normalize/diff." Readers: `fold_stereo_atom_group` at `3022`,
`3031`, `3089` and its twin; `FrameTransport for StereoAtomDelta` gates admissibility on it at
`1232` (`1296` for bonds); `reaction.rs:1595-1601` forwards it into the lowered edit. Probe: a
`ModifyConstraint` with `kind: None` beside a tetrahedral configuration change normalizes to
`kind: Some(Tetrahedral)`. Defense: doc 134 §6 I6e is the origin ("serialization-only",
"threads it through"), and `apply_stereo_atom_change` does ignore it. Verdict: the sentence
contradicts both the code and the recorded decision.

**F6 — the property doc on `test_deltas_normalize_idempotent` claims an order-independence
check the test does not perform** (`tests/property/delta.rs:261-269`; confirmed; FE2). Text: "the
confluence check — the normal form is unique, independent of the order the deltas arrived in";
the body asserts `N(N(x)) == N(x)` on one input order. Within-entity order is semantic by the
type's own contract (`delta.rs:3252-3255`): the chain `[0→1, 1→2]` normalizes and `[1→2, 0→1]` is
`Contradiction`. Claim: property-tests guide §"Specification and evidence" and §"Test source
should point back to the property"; code-reviews guide §Documentation (a failing property must
identify its semantic reason). Defense: the executable assertion is the stated idempotence law;
"order" could be read as cross-entity order. Verdict: confirmed; complementary to F16.

**F7 — "Different entities are independent" overstates acceptance** (`delta.rs:3308`; behavior
at `3339-3436`; reduced; FE3). A created bond alone normalizes `Ok`; with `Atom Remove 0` in the
same collection it is `Contradiction`. The rejection is stated only in line comments (`:3358`,
`:3410-3411`). Claim: code-reviews guide §Documentation; property-tests guide §"Public
documentation shape" lines 503–509. Defense: "independent" describes fold outputs, which is true;
an explicit statement of every failure condition is the `# Errors` migration. Verdict: reduced to
the overstatement; the wording should state whatever scope the F1 fix settles.

**F8 — the `EntityOp` doc describes a two-kind abstraction that has six implementors**
(`delta.rs:1812-1813`; confirmed; FE4). "abstracting `AtomDelta`/`BondDelta` … (`()` for an atom,
its two ids for a bond)" against `Atoms` types `(Vec<AtomId>, AtomId)`, `Vec<AtomId>`, and
`[AtomId; 2]` on the four overlay impls landed 2026-06-30. Defense: the parenthetical reads as
examples; `type Atoms` is visible on every impl. Verdict: the head clause names two of six
implementors as the abstraction's domain, preserving the 2026-06-26 state.

**F9 — six `*Delta` enums are documented as "a resolved edit"** (`delta.rs:57`, `181`, `289`,
`395`, `523`, `651`; confirmed; FB1 and FE5, found independently by both areas). The stereo deltas
(`:751`, `:878`) and `ConstraintDelta` (`:1310`) say "a resolved change". Claim: nomenclature
guide §"Edit and undo" ("**Not:** *delta*") and §"Delta and Update"; the guide is normative for
explanatory documentation. Defense: the module's founding term ("Add Delta (resolved Edit)
type"); lowercase "edit" as plain English; the entry post-dates the prose. Verdict: each sentence
is individually wrong under the settled entries and the module carries two words for one
doc-comment position;
the conforming word is "change".

**F10 — `Normalize for EntitySpan<T>` states no contract for its `Modified → Unchanged`
collapse** (`delta.rs:1710-1727`; variant docs `1654-1657`; confirmed at low severity; FE6). The
impl landed 2026-08-29, after the `# Semantic properties` convention; `impl Normalize for
ReactionSpan` (`reaction_span.rs:1908`) is also undocumented, and `reaction_span.rs:145` documents
that construction does not collapse. Claim: property-tests guide §"Semantic properties belong with
the public API" lines 96–105 (state the relation the assertion uses); data-types guide lines
297–301 and 308–310. Defense: a common law may be documented once at a discoverable shared
location. Verdict: no such location exists.

**F11 — the `Delta` type doc describes a hypothetical SqPO application policy**
(`delta.rs:1361-1364`; confirmed as narrowed; FE7). The sentence beginning "An SqPO application
policy could instead cascade" restates doc 184 §"Deltas under SqPO application". Claim: AGENTS.md
("Code states current behavior; discussion documents preserve reasoning and history");
code-reviews guide §Documentation. Defense: it explains why completeness is not applicability.
Verdict: the preceding sentence already states the behavioral point; the SqPO sentence goes.

**F12 — `Normalize for Deltas` rustdoc says "inconsistent" for an unsatisfiable set**
(`delta.rs:3312`; reduced; FB7). Claim: nomenclature guide §Normalize ("returns
`Err(Contradiction)` for an unsatisfiable represented value"); the §Contradiction and
§Inconsistency entries name each other as confusable neighbours. Defense: ordinary English; the
whitepaper and `ApplyPreconditionError::InconsistentReaction` use the adjective. Verdict: reduced
to a wording change ("unsatisfiable" or "contradictory"); the collision class did not survive;
the consumer variant is out of scope and noted under open items.

### Naming

**F13 — four unit-test names encode behavior** (`delta.rs:4666`, `4672`, `4690`, `5021`;
confirmed; FB2). `_field_noop_dropped`, `_created_absorbs_field`, `_created_then_removed_cancels`,
`_remove_subsumes_field`. Claim: test-writing skill §"Names and assertions" and §Prohibitions;
the rule existed from 2026-05-09, the tests date from 2026-06-25. Defense: single-scenario tests
with no case table, so the suffix is the only scenario statement. Verdict: confirmed; the skill
supplies the pattern (scenario noun, or rows of `test_deltas_normalize`).

**F14 — "canonical" names normalized equivalence in two test names and a local**
(`delta.rs:3631`, `4612`, `4244`; confirmed with corrected dating; FB3). Both tests decide by
`normalized_eq` (`diff_field` at `1577`, `superimpose` at `1736`). The names date from 2026-07-04,
when entity-level `canonical_eq` was the value comparison; the module-wide rename to
"normalize" on 2026-08-10 migrated the sibling names and left these. Claim: nomenclature guide
§"Canonical and canonicalize" ("**Not:** *normalize*") and §Normalize. Defense: the names
followed the vocabulary of their day; other `_canonical` test names exist in the crate (they
compare aggregates across order, a different concept). Verdict: residues of an executed
migration are individually wrong; the conforming term is "normalized".

### Test witnesses

The confirmed and reduced witness findings rest on the code-reviews guide's two-direction rule
(every stated law has a test) and the data-type-contracts skill §"Derive verification from the
contract" (exact unit cases for each meaningful success and failure boundary). Coverage counts
are from the full lib and property runs and were re-measured by the refutation pass.

**F15 — the created-entity-on-removed-atom guard has a witness for bonds only**
(`delta.rs:3371`, `3382`, `3393`, `3404`, `3420`, `3431` count 0; witness `5086-5099`; confirmed;
FD5). Defense: the ordinary arms copy the bond arm; the DPO validator owns dangling for lhs
entities. Verdict: the two stereo arms differ structurally (site chained for stereo atoms,
omitted for stereo bonds), so the copy argument does not cover them.

**F16 — cross-entity order independence is stated three times but has only a trivial witness**
(`delta.rs:3252-3258`, `3308-3312`; witness `5064-5075`; confirmed; FD7). The unit witness permutes
two singleton groups of different variants, which the final sort orders identically regardless of
the fold; the property asserts idempotence only. Interleaved two-entity chains do normalize
equal. Defense: independence follows structurally from per-id grouping and sorting. Verdict: a
stated law whose witness cannot fail. Owner: a shuffle of interleaved multi-element groups.

**F17 — created-entity contradiction paths have no witness** (`delta.rs:1916`, `2929`, `2953`,
`2963`, `3108`, `3118` count 0; confirmed; FD9): an operation after a created entity's `Remove`,
a second stereo `Add`, and a non-connecting stereo configuration chain. Defense: the stereo folds
mirror the generic fold; doc 229 S2b's constructor `DuplicateReference` rejects double `Add` for
`Reaction` inputs. Verdict: each is a failure boundary of a public operation. The
operation-after-`Remove` rule is stated nowhere in rustdoc.

**F18 — the stereo preserved-entity fold is largely unwitnessed; the stereo-bond preserved
fold has none** (`delta.rs:3014`, `3034`, `3037-3041`, `3064`, `3070`, `3169`, `3189`,
`3192-3196`, `3217-3221`, `3225` count 0; confirmed; FD10): constraint-chain fusion, the
`(None, None)` skip, constraint reversion on `Remove`, configuration mismatch on `Remove`, and
operation after `Remove`. The only stereo preserved witnesses are the three stereo-atom rows at
`4187-4251`, without constraints. Defense: documented twin of the atom fold; S2a scoped the
created path only. Verdict: S2a's scope explains but does not defer the absence.

**F19 — the involution property's generator never yields `Constraint` deltas and yields
`ModifyConstraint` for one kind only** (`tests/property/delta.rs:271-276` over
`comprehensive_reaction_strategy`; confirmed at low severity; FD12). Statistics in the table
above. Defense: the shared generator meets its documented aim; `inverse` is a field swap; unit
witnesses exist for `ConstraintDelta::inverse`. Verdict: a variant-coverage finding against the
law stated over all deltas. Owner: a delta-level strategy over the nine `Delta` variants.

**F20 — stereo delta `FrameTransport` `None` paths and ordinary-kind `ModifyConstraint`
transport are never executed** (`delta.rs:1195`, `1211`, `1233`, `1259`, `1275`, `1297`;
`1039-1043`, `1081-1092`, `1130-1141`, `1173-1177`; confirmed; FD14). The stereo-atom delta has no
`None` witness at all. Claim: property-tests guide lines 423–427 and 437–438 name deltas among the
carriers whose `None` behavior must be tested; `traits.rs:178-181`. Defense: doc 214 S0v placed the
systematic `None` properties in `tests/property/frame.rs` at form level; delta transport is thin
dispatch. Verdict: the guide names deltas explicitly.

**F21 — the stereo-atom created-path old-value mismatch is never executed** (`delta.rs:2981`
count 0 while its stereo-bond twin `3136` is executed; reduced; FD8). Did not survive: per-kind
"matrix holes" for dative, multicenter, and noncovalent old-value mismatches and the combined
reordered-plus-intervening-change path, because doc 229's requirement is category-level and the
generic boundary (`:1936`) is executed.

**F22 — `Remove` after an `Electrons` `ModifyField` on an aromatic or multicenter entity has no
witness** (`delta.rs:2447-2449`, `2568-2570` count 0; reduced; FD11). A vector-valued,
frame-sensitive revert through `apply_field`'s `normalized_eq` check, not covered by the atom
witness. Did not survive: dative and noncovalent, whose `field_inverse` are one-line delegations
of an executed generic path.

**F23 — the `ModifyConstraint` arm of each `apply_*_change` adapter has no witness**
(`delta.rs:2802-2803` and the seven twins, count 0; reduced; FD15). The arm with logic
(`compare_and_set`) is executed by no test in the crate. Did not survive: the `Add`/`Remove`
no-op arms, which are trivially total.

**F24 — `EntitySpan::superimpose`'s `Added` and `(None, None)` outcomes and
`ConstraintSpan::normalize`'s `Removed` arm have no direct witness** (`delta.rs:1739`, `1740`,
`1799` count 0; reduced; FD16). The five-outcome contract at `1730-1734` has a one-row test. Did
not survive: `ConstraintSpan::lhs`/`rhs`, trivial accessors with no stated law.

### Generators

**F25 — no generated input reaches the S2a cancellation path** (`tests/property/delta.rs:172-258`;
reduced; FD2). All 574 sampled `Deltas` with a same-id `Add`+`Remove` pair normalize to
`Contradiction` because `Add` and `Remove` draw independent forms whose pairwise `normalized_eq`
rate is zero; the reaction generators produce no such pair. Before S2a every pair cancelled.
Survives as a degeneracy finding: S2a changed what the property exercises without a generator
change. Did not survive: the reading of doc 229's S2a result as a coverage claim; it lists the
property among passing gates.

**F26 — `deltas_strategy` emits only atom, bond, and constraint deltas**
(`tests/property/delta.rs:249-258`;
reduced; FD3). The only bare-`Deltas` normal-form property samples none of the six overlay and
stereo kinds and states no domain. Did not survive: the implication that the overlay and stereo
folds lack generative coverage altogether; consumer-level properties
(`test_reaction_compose_normalized_deltas`,
the nine-law suite over `Reaction`) reach them.

**F27 — `test_deltas_normalize_idempotent` discards 28.5 % of samples and 81 % of the samples
that fold anything** (`tests/property/delta.rs:264-269`; reduced; FD4). Only 6.3 % of samples reach
a multi-operation fold on a consistent input. Survives as a degeneracy observation under the
property-tests guide's preference for deriving valid aggregates over rejection. Did not survive:
the named-contradictory-domain requirement, which post-dates the generator.

**F28 — six `*_delta_diff_apply` properties assert the patch law under `==` without stating
the relation or the normalized-domain precondition** (`tests/property/delta.rs:278-283`, `309-314`,
`331-338`, `356-363`, `381-388`, `406-413`; reduced; FD1). The two stereo siblings in the same file
assert `normalized_eq`. Reproduction: charge `Lit(1)` against `lit_set([1])` yields an empty diff,
`applied == rhs` false, `applied.normalized_eq(&rhs)` true; the property never sees it because
every sampled form is normalized (0 non-normalized in 10,000 for all eight form strategies).
Claim: nomenclature guide §"Patch algebra" ("read under `normalized_eq`, not `==`");
property-tests guide lines 58–60 (a precondition the assertion requires must be reflected in its
documentation). Defense: on a normalized domain `==` is a valid stronger check, and the guide's
own example at `property-tests.md:135` uses `==`. Did not survive: the wrong-semantics class.

## Open questions from the refutation

- **The `id` surface of the eight `*Delta` enums** (FC3, merged with structure open item O3).
  `StereoAtomDelta::id` and `StereoBondDelta::id` are `pub` inherent; the six other kinds expose
  `id` only through the `pub(crate)` `EntityFold` trait whose doc restricts its use to the fold
  and span lowering. `Delta::collect_frame_action_domain` and `Delta::reframe_by_actions`
  (`delta.rs:1402-1440`, `1467-1517`) therefore hand-match four variants per overlay kind, and
  `reaction.rs:2522-2659` repeats the matches. Neither shape is required by any document; whether
  every `*Delta` exposes `pub fn id`, and which `EntityFold` members are supported operations, is
  a public-API decision.
- **Whether exact transport cases must pin action direction** (FD6). All seven S2a reordered
  rows use single transpositions, so swapping the arguments of `between` in any `remove_matches`
  or stereo fold would leave all 27 rows green; the code is correct (3-cycle probe). The
  property-tests guide's sentence on nonidentity actions (lines 425–427) is written for generated
  values. Settled by extending it to exact transport cases (then a 3-cycle aromatic or multicenter
  row and a stereo row on a kind with more than two cosets are required) or by recording that
  exact cases need not be mutation-resistant.

## Refuted findings

Recorded with their dismissing citations so they are not re-flagged.

- **Asymmetric reference guard** (FA5): the guard rejects a created stereo atom on a
  net-removed atom but admits a created stereo bond on a net-removed bond and a
  `ConstraintDelta::Add` on a removed atom. Dismissed by doc 229 §"Proposed changes" lines
  1573–1576 (wider consistency checks stay with their existing first consumers) and the integrity
  guide's lazy-property list; `to_reaction_span` reports `Contradiction` for the bond case without
  panic and no integrity contract is breached. The finder flagged the tie; ties favor the original.
  F1 is not covered by that deferral.
- **`EntityFold::split` collides with the guide's *Split* entry** (FB4): the entry (2026-08-06)
  post-dates the method (2026-06-25); reclassified as a migration (below).
- **`StereoConfigFold` clips "configuration" onto the `*Config` suffix** (FB5): both the suffix
  row and the clipped-abbreviation rule post-date the enum (2026-07-01); folded into the
  clipped-identifier migration.
- **`remove_matches` says "matches" for a `normalized_eq` check** (FB6): doc 214 line 216 and
  doc 215 lines 180–181 name the operation "removal matching" with the value relation left to the
  caller, and the nomenclature guide's Matching entry asks only to qualify the term, which
  `remove_` does.
- **"set" for the `Modify` vocabulary in case labels and a private variant** (FB8): no guide
  entry or retired row governs "set"; the only basis is a rank-4 record of variant names, which
  the identifiers honour.
- **"before/after" and "previous" for the `old`/`new` sides in prose** (FB9): no guide entry
  reserves the words; the identifiers carry the distinction; the change is aesthetic.
- **`EntitySpan`/`ConstraintSpan` housed in the delta module** (FC1): the grouping rule
  (2026-08-31) post-dates the placement (2026-06-26) and the split is deliberately tracked (docs
  168 and 227, doc 229 line 1200); folded into the split proposal, which moves them to a span
  module.
- **Focus types `Delta` and `Deltas` defined last** (FC2): the rule post-dates the layout; the
  crate precedent invoked is unwritten and not a normative source; the sibling `edit.rs` has the
  same shape; folded into the split proposal, whose owning module begins with them.
- **No exact `inverse` witness for four overlay kinds** (FD13): no normative source requires a
  witness per method; the involution law is property-tested for every kind the generator emits;
  `inverse` has no failure boundary.

## Migrations

Conventions adopted after the reviewed code was written; recorded once each and addressed in the
proposed design rather than the findings list.

1. **Contract sections.** No `# Errors` or `# Semantic properties` section exists in `delta.rs`.
   The laws are stated in prose and tested: `Delta::inverse` involution (`:1358-1360`, `:1382`),
   the `Deltas` order contract (`:3252-3258`), the `Normalize for Deltas` normal form
   (`:3308-3312`, to which S2a appended its cancellation sentence rather than opening a section),
   `EntitySpan::superimpose` (`:1730-1734`), the `None` conditions of the eight `FrameTransport`
   impls, and the `Contradiction` conditions of `EntitySpan::normalize`, `ConstraintSpan::normalize`,
   and the macro-generated `apply_field`/`apply_constraint`.
2. **Module size.** 5157 lines (3462 non-test) against the guide's roughly 1000. The file was 922
   lines on 2026-06-25, 3438 by 2026-07-12, and 5361 when the guide landed; the split is tracked
   in doc 168 §"Module and file organization", doc 227, and doc 229 line 1200. The proposed split
   and its visibility map follow.
3. **Inline test module.** 1693 lines inline against "long test suites in separate test
   submodules"; crate-wide 80 inline suites against 3 separate `tests.rs` files; S2a added 322
   inline lines after the guide landed. Target: a `tests.rs` child, per the crate precedents
   `ir/molecule/tests.rs`, `dsl/molecule/tests.rs`, `ir/canonicalize/tests.rs`.
4. **`pub(crate)` markers.** 22 in `delta.rs`, all predating the AGENTS.md visibility rule
   (2026-09-03); 6 are forced by `pub(crate)` parameter types (`ConstraintFrameActionDomain`,
   `ConstraintFrameActions`); the remaining 16 are the six `uses_participant_frame` methods, the
   `EntityOp`/`EntityFold` pair, and the eight `apply_*_change` functions. The regroupings that
   would remove them are public-API decisions (below).
5. **Clipped identifiers.** `EntityOp`/`op`/`ops`, the lowercase `config` locals, and
   `StereoConfigFold`/`fold_stereo_config` predate the complete-words rule (2026-08-08); the same
   clips are crate-wide (`RelOp`, `MemOp`, `CosetOp`, `config` in `stereo.rs` and `reaction.rs`).
   Whether the rule reaches crate-private identifiers is open. Conforming spelling if adopted:
   `StereoConfigurationFold`, `fold_stereo_configuration`.
6. **`EntityFold::split`/`rebuild` under the Split entry.** The guide reserves `split` for
   decomposition into connected components (2026-08-06); the method (2026-06-25) deconstructs a
   delta into its `EntityOp`. The guide supplies no term for delta deconstruction; consult before
   naming.
7. **Strategy-domain documentation.** None of the delta-suite strategies states that it emits
   normalized values (`property-tests.md:211`, 2026-08-06); all sampled 100 % normalized. F26 and
   F28 record the consequences.

### Proposed split and visibility map

Layout per crate precedent (`reaction.rs` with `reaction/{dpo,integrity}.rs`; `constraint.rs`
with per-kind children and `pub use`; `molecule.rs` with `molecule/tests.rs`); no `mod.rs`.
Module names are open (nomenclature guide §"Module names"; consult before naming).

| Proposed file | Content (current lines) | Non-test lines |
| --- | --- | --- |
| `ir/delta.rs` (owning) | header; `mod` declarations; `pub use` of children; `Deltas` with inherent methods, `FromIterator`, `IntoIterator` (3250–3305); `Normalize for Deltas` (3307–3461); `Delta` with impl and `FrameTransport` (1352–1554); `ConstraintDelta` with impls (1310–1350); `options_normalized_eq` (1617–1623), `transport_optional` (1002–1010); `diff_field_ops!`, `fold_field_ops!` (1558–1613, 1627–1648), placed before the `mod` declarations for textual scope | ~640 |
| per-kind child, atom | `AtomDelta` (57–179); `EntityPatch`/`EntityFold` impls (2017–2118); `apply_atom_change` (2796–2807) | ~240 |
| per-kind child, bond | 181–287; 2120–2227; 2809–2820 | ~230 |
| per-kind child, dative | 289–393; `FrameTransport` 1012–1047; 2229–2354; 2822–2835 | ~290 |
| per-kind child, aromatic | 395–521; 1049–1096; 2356–2471; 2837–2850 | ~310 |
| per-kind child, multicenter | 523–649; 1098–1145; 2473–2592; 2852–2865 | ~315 |
| per-kind child, noncovalent | 651–749; 1147–1181; 2594–2707; 2867–2880 | ~270 |
| per-kind child, stereo | `StereoAtomDelta`, `StereoBondDelta` (751–1000); `FrameTransport` (1183–1308); `EntityPatch` impls (2709–2790); `apply_stereo_*_change` (2882–2910); `StereoConfigFold`, `fold_stereo_config`, both stereo group folds (2912–3248) | ~830 |
| fold child | `EntityOp`, `EntityFold`, `fold_group`, `fold_created`, `fold_preserved` (1812–2015) | ~215 |
| span sibling under `ir/` | `EntitySpan`, `ConstraintSpan` with impls (1650–1810) | ~175 |
| `ir/delta/tests.rs` | the inline `mod tests` (3463–5157) | ~1690 |

Visibility map (current → proposed). "pub in private child" is a `pub` item inside a private
`mod` that only the owning module and its descendants can name; it introduces no scoped marker.

| Symbol | Current | Proposed |
| --- | --- | --- |
| `Deltas` and its methods; `Delta`, `Delta::inverse`; `ConstraintDelta`, `::inverse`, `FrameTransport` | pub | pub in the owning module; `ir.rs` re-exports unchanged |
| `Delta::collect_frame_action_domain`, `Delta::reframe_by_actions`, `ConstraintDelta::collect_frame_action_domain`, `ConstraintDelta::reframe_by_actions`, `ConstraintSpan::collect_frame_action_domain`, `ConstraintSpan::reframe_by_actions` | pub(crate) | pub(crate), forced by `pub(crate)` parameter types (`constraint.rs:39-41`, `constraint/molecule.rs:34,140`); consumers `reaction.rs:559,561,2796`, `reaction_span.rs:1963,2037` |
| eight `*Delta` enums, their `inverse` and `for_update` | pub | pub in per-kind child, re-exported by the owning module; `ir.rs` unchanged; `for_update` consumers `dsl/reaction.rs:728-1043` |
| `StereoAtomDelta::id`, `StereoBondDelta::id` | pub | pub; extension to the six other kinds is the open `id`-surface question |
| six `*Delta::uses_participant_frame` | pub(crate) | pub(crate) retained (consumers: the owning module and `reaction.rs:2532-2747`); removable if `reaction.rs::reframe_application_deltas` routes through `Delta::reframe_by_actions` |
| `FrameTransport`, `EntityPatch`, `EntityFold` impls | trait impls | in per-kind children |
| eight `apply_*_change` | pub(crate) fn | pub in per-kind child plus `pub(crate) use` re-export (consumers `reaction_span.rs:2321-2561`); removable by an inherent apply method per kind, a public-API decision tied to the nomenclature guide's `EntityPatch::apply` open item |
| `EntityOp<F>` | pub(crate) | pub in the fold child, re-exported `pub(crate)` only while `EntityFold` is crate-visible |
| `EntityFold` | pub(crate) trait | pub in the fold child plus `pub(crate) use` (consumer `reaction_span.rs:32,1455-1534`); alternatives: make it `pub` beside `EntityPatch`, or move `id`/`into_delta` onto `EntityPatch` and keep fold-only members private (public-API decision) |
| `fold_group`; `fold_stereo_atom_group`, `fold_stereo_bond_group` | private | pub in private child (called by the owning module) |
| `fold_created`, `fold_preserved`, `StereoConfigFold`, `fold_stereo_config` | private | private |
| `options_normalized_eq`, `transport_optional`, the two macros | private / module-local | private in the owning module; descendants access ancestors' private items |
| `EntitySpan<T>` and `ConstraintSpan` with their inherent methods and impls | pub | pub in the span sibling; `ir.rs` re-exports from it; the six overlay modules, `incidence.rs`, `canonicalize.rs`, `reaction_span.rs`, and `dsl/reaction_span.rs` change their import path; no cycle (the span module depends only on `traits`, `constraint`, `entity`, `frame`, `error`) |
| `mod tests` | inline | `#[cfg(test)] mod tests;` in `delta/tests.rs`; `use super::*` keeps working |

Cross-module dependencies the split creates resolve by `pub` inside private children and
ancestor-private access; the 22 existing `pub(crate)` markers map to 22, or to the 6 forced ones
if the three public-API regroupings above (inherent apply per kind, `EntityFold` publicity,
routing reaction transport through `Delta::reframe_by_actions`) are approved. No `pub(super)`, no
`pub(in ...)`, no new public import path.

## Open items

1. **Owner of the cancelled-created-id check (F1).** Either `Normalize for Deltas` treats a
   cancelled created id as absent for every surviving `Add` (or returns `Contradiction` when a
   surviving delta references it), keeping bare `Deltas::normalize` self-consistent; or
   `Reaction::normalize` and `Reaction::reframe` gate their published value, leaving bare
   `Deltas::normalize` as it is. The first keeps the check with the operation that already owns
   the sibling guard; the second keeps `Deltas` an open carrier and makes the reaction
   transformations invoke the gate the data-types guide names.
2. **Equivalence preserved by `Deltas::normalize`.** The data-types guide (lines 308–310) asks an
   operation deriving a delta normal form to state its roundtrip under the relevant semantic
   equality; no relation between a `Deltas` value and its normal form is stated anywhere. F2
   narrows the answer space.
3. **`kind` stamp selection under mixed input (F4)**: derive the stamp deterministically (from the
   entity's configuration chain or each constraint's own stamp) or reject mixed stamps as
   `Contradiction`. Related: the bespoke stereo folds outlive their recorded reason (doc 134 I6a
   recorded relative operations with no `EntityOp` image; those were removed 2026-08-11); the
   remaining obstacle to the generic fold is the `kind` field on stereo `ModifyConstraint`, whose
   place in the IR delta is a design question.
4. **The `id` surface and `EntityFold` membership** (open question above).
5. **Direction-detecting exact transport cases** (open question above).
6. **Unit `Contradiction` as the only failure of bare `Deltas::normalize`.** Representation-domain
   failures (oversized stereo frames, repeated participants, electron-count length mismatch) are
   reported under the semantic label; deferred by the data-types guide (lines 1105–1106,
   "pending the repository-wide error review") and fixed by the trait signature.
7. **Guide discrepancies.** The nomenclature guide's Patch algebra entry names `EntityPatch::apply`,
   which has never existed (the operation is the eight `apply_*_change` functions); tied to the
   guide's open issue 3. The property-tests guide writes the update law with `==` (line 135)
   while the nomenclature guide reads both spellings under `normalized_eq`;
   `test_unpaired_electrons_form_difference_to` asserts `==` while its siblings assert
   `normalized_eq`.
8. **Nomenclature items without a guide term**: the generic incidence carrier named `atoms`
   (`EntityFold::Atoms`, `EntityOp::{Add, Remove}::atoms`, the `remove_matches` parameters) where
   the guide's generic noun is participant; entity-kind truncation in crate-private names
   (`apply_dative_change` against `apply_stereo_atom_change`); "fuse" carrying three senses
   (`fuse_field`, fused reframe, fused rings); "a dynamic entity" in the `EntitySpan::Modified`
   doc (the guide says "relabeled"); "id space" against "entity-id frame" inside the guide itself;
   names for `split`/`rebuild` (migration 6); module names for the split.
9. **Test-module import style**: the inline suite mixes `use super::super::x` with
   `use crate::ir::x`; precedents are mixed and no guide states a rule.
10. **Unreachable or uncalled code**: `EntityPatch::diff_field`/`diff_constraints` have no caller
    because every impl overrides `diff` with `for_update`; the `fold_preserved` `Add` arm
    (`:1963`) and its stereo twins (`:3017`, `:3172`) cannot execute; the stereo
    `modify_field`/`modify_constraint` fallbacks are documented as unused. A `ModifyConstraint`
    with different `old`/`new` keys passes the preserved path unchecked (lazy continuity;
    `:1980` has coverage count 0).
11. **Generator design**: `stereo_coset_for_kind` emits only `Lit` and `Undetermined` cosets, so
    the stereo diff/apply properties also run on a fully normalized domain; whether
    `deltas_strategy` should sample stereo kinds with more than two cosets so that transport
    direction becomes observable.
12. **Consumer-side twins, out of scope**: `reaction_span.rs:2366-2371` indexes a map inside a
    `Result`-returning operation (unreachable after the F1 fix; whether to index or
    `ok_or(Contradiction)` is for the reaction-span review); `apply` reports
    `CorrespondenceMismatch` for the F1 input, a label naming the correspondence rather than the
    deltas; `ApplyPreconditionError::InconsistentReaction` is the consumer twin of F12.

Observations recorded by the refutation pass that no agent pushed, for triage: the `charge_set`
test constructor helper (`delta.rs:4623-4631`) against the test-writing skill's no-helper rule in
force since 2026-05-09; "the rhs-hand value" at `delta.rs:2794`; coverage count 0 at `delta.rs:1980`
and at `Deltas::len` (`:3275-3277`); partially unexecuted arms of
`Delta::collect_frame_action_domain`
and `Delta::reframe_by_actions`.

## Proposed design

### `Normalize for Deltas`

- Reduce every carried payload in the rebuild step of `fold_created`, `fold_preserved`,
  `fold_stereo_atom_group`, and `fold_stereo_bond_group`: `attributes`, `old`, and `new` pass
  through their own `Normalize::normalize`, matching `Molecule::normalize` and the existing
  `ConstraintDelta` arm. This resolves F2 and F3 together; the stereo folds then compare with
  `normalized_eq` like the generic path. Document the resulting normal form and the equivalence
  it preserves (open item 2) in a `# Semantic properties` section, and its `Contradiction`
  conditions under `# Errors`, closing migration 1 for this operation.
- Close F1 at the owner chosen under open item 1.
- Settle the `kind` stamp under open item 3 and implement the chosen rule in both stereo folds.

### Documentation corrections

- `StereoAtomDelta::ModifyConstraint::kind` (`:776-779`, `:902`): state that normalization
  threads and derives the stamp and that frame transport gates on it (F5).
- `test_deltas_normalize_idempotent` doc: state the idempotence law it tests (F6).
- `Normalize for Deltas` doc: replace "Different entities are independent" with the acceptance
  rule the F1 owner settles (F7); "unsatisfiable" or "contradictory" for the `Err` condition
  (F12).
- `EntityOp` doc: describe the per-kind `Atoms` carrier over all implementors (F8).
- Six `*Delta` docs: "A resolved change to a single …" (F9); "modification" at `:1352` rides
  along.
- `Normalize for EntitySpan<T>`: a `# Semantic properties` section stating the `Modified →
  Unchanged` collapse and its relation (F10).
- `Delta` type doc: drop the SqPO sentence (F11).
- Contract-section adoption for the remaining prose laws (migration 1) once the normal-form
  relation is settled.

### Tests and generators

- Exact rows in `test_deltas_normalize_error` and `test_deltas_normalize`: one created-entity
  guard row per guarded kind including the stereo-atom site and stereo-bond ligand variants
  (F15); operation after a created entity's `Remove`, stereo double `Add`, non-connecting stereo
  configuration chain (F17); a stereo-bond preserved-fold table mirroring the stereo-atom one,
  plus constraint-chain, constraint-revert, configuration-mismatch, and `(None, None)` rows for
  both stereo kinds (F18); `stereo_atom_old_value` (F21); `Remove` after an `Electrons` change on
  aromatic and multicenter entities (F22); an interleaved multi-element cross-entity shuffle
  (F16); and, if open item 5 settles for detectability, a 3-cycle aromatic or multicenter row and
  a stereo row on a kind with more than two cosets.
- `FrameTransport` `None` rows for the stereo-atom delta (degree mismatch, inadmissible kind),
  the stereo-bond `Add`/`Remove` degree checks, and `ModifyConstraint` transport for the four
  ordinary overlay kinds (F20).
- A table over the eight `apply_*_change` functions with `ModifyConstraint` rows (F23);
  `EntitySpan::superimpose` rows for all five outcomes and a `ConstraintSpan::normalize` `Removed`
  row (F24).
- Generator revision for `deltas_strategy`: derive `Remove.attributes` from the `Add` state in a
  fraction of same-id samples, including a transported frame (F25); add overlay and stereo delta
  arms (F26); chain `ModifyField` values within a group and separate a named contradictory domain
  with its own exact-error property (F27); document the emitted domain (migration 7).
- A delta-level strategy over all nine `Delta` variants for the involution property (F19).
- `*_delta_diff_apply` docs and assertions: state the relation and the normalized-domain
  precondition, or assert `normalized_eq` as the stereo siblings do (F28; interacts with open
  item 7).
- Rename the four behavior-named tests or fold them into `test_deltas_normalize` rows (F13);
  rename the two `_canonical` tests and the `canon` local to "normalized" (F14).

### Structure

The split in migration 2 with its visibility map, pending the module names (open item 8) and the
`id`-surface decision (open item 4); the test relocation of migration 3 travels with it.

### Records

- Doc 229's S2a result lists "the delta normal-form property" among passing gates; recording the
  property's evidence scope there (no generated cancellation input after S2a) prevents the
  reading that the property covers the change.
- Guide clarifications under open item 7 and open item 5 belong to the living guides, not to
  code.

## Process notes

- Dating a convention by `git log --follow` on the current skill path misses its earlier
  `.claude/skills` home; date rule text with `git log -S` across both catalogs. This changed one
  migration record into a confirmed finding (F13).
- The area briefs introduced finding classes the guide does not define (`duplication`,
  "public method with no test"). Where no guide or skill rule exists the finding cannot survive
  (FD13, the duplication half of FC3); where the data-type-contracts skill's boundary-case rule
  applies, witness findings survive on that basis. The next brief should name that skill section
  as the basis for witness findings and should not steer an area toward pushing under a guide
  that post-dates the code (FC1, FC2).
- Two deferral checks missed governing text that changed the verdict: doc 229 §"Proposed
  changes" lines 1573–1576 (FA5) and doc 134 §6 I6e (FA4). Both are reachable from the status
  index; the brief's governing-document list should name completed implementation records at
  stage granularity, not only by document.
