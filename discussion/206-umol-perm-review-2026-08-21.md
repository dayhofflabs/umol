# 206 — umol-perm review cycle

Status: Completed
Date: 2026-08-21
Relates: [119](119-umol-perm-review-2026-06-21.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[157](157-umol-perm-fallibility-and-arguments-2026-07-20.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[code reviews guide](../docs/development/code-reviews.md),
[data type contracts guide](../docs/development/data-types.md)

## Purpose

First review-cycle record for the umol-perm crate, produced under the
`review-cycle` skill and the code-reviews guide. Five area reviews
(construction/integrity/fallibility, nomenclature, module structure and
visibility, tests and generators, documentation) ran against pinned commit
`836f7402`, followed by a refutation pass over the pooled findings. This
document records the surviving findings with both arguments, the refuted
findings, the consolidated open items, and the proposed design. It contains no
staged implementation plan; implementation is scheduled separately after
triage.

Verdict tally: 9 confirmed, 7 reduced, 2 refuted, before consolidation of
cross-area duplicates. Headline: check placement is fully compliant — of
roughly 30 fallible surfaces, every check classifies as tier-1 integrity or
first-requiring-operation validation, with zero defensive, premature, or
inherited-fallibility occurrences, and checked/asserted constructor pairs
share one implementation. Doc 157's settled decisions are implemented exactly.
The surviving defects are concentrated in contract documentation, missing test
witnesses, and a small set of semantically inaccurate doc comments.

Provenance: at the pinned commit, `docs/development/code-reviews.md` and the
"Contract sections in doc comments" section of
`docs/development/data-types.md` were uncommitted working-tree state; the
review applied the working-tree versions and the refutation pass verified that
no confirmed fallibility or documentation finding rests solely on the
uncommitted texts (the committed `data-type-contracts` skill independently
requires the contract headings). Test evidence: 163 unit cases and 26
properties, all green; measured line coverage 94.8 % overall (permutation,
group, error 100 %; coset 97.5 %; oriented 97.1 %; class 80.0 %).

## Confirmed and reduced findings

Each finding records the claim, the strongest case for the original (the
defense), and the refutation verdict. Line references are at commit
`836f7402`.

### Semantic accuracy

**F1 — `CosetSpace::count` formula is wrong for cis/trans and axial spaces**
(`coset.rs:89-92`; confirmed). The doc "The number of cosets, `n! / |R|`"
is false for the CT/AX spaces, whose parent is the order-8 wreath group
`S2 ≀ S2`, not `Sn`: `count()` is 8/4 = 2 while the formula gives 24/4 = 6.
The constructor asserts `numbering.len() == parent.order() / group.order()`
and the unit tests pin CisTrans → 2, Axial → 2. Defense: exact for every
`Sn`-parent class; the formula is inherited from doc 104's era when all
parents were `Sn`. Verdict: the struct's own doc describes non-`Sn` parents
two paragraphs above; correct statement is `|P| / |R|`.

**F2 — struct doc says left cosets `P/R`; the algebra is right cosets `R\P`**
(`coset.rs:25`; confirmed; found independently by the nomenclature and
documentation reviews). The module header, the nomenclature guide's Coset
entry ("the prose establishes right cosets `Rσ`, so `R\P` is correct and the
struct comment is wrong"), and the implementation (`coset_rep` minimizes over
`r.compose(permutation)`, i.e. over `Rσ`) all fix `R\P`. Defense: none
available — the only counter-reading (the module header is the wrong one) is
adjudicated against by the guide itself.

**F3 — `dihedral` doc claims "(order 2n)"; false at degrees 1 and 2**
(`group.rs:63-64`; confirmed). `reflection(2)` is the identity, so
`dihedral(2)` is {e, (0 1)} of order 2; `dihedral(1)` is trivial. Both inputs
are reachable without panicking; tests cover degrees 3–5 only. Defense: "of
the n-gon" conventionally implies n ≥ 3; no current stereo class uses lower
degrees; 157 settled only degree-0 documentation. Verdict: 157's rationale —
edge behavior documented rather than discovered — applies equally, and the
doc already carries the degree-0 edge.

**F4 — `from_image` doc omits the maximum-degree panic and is affirmatively
misleading** (`permutation.rs:41-45`; confirmed). "Panics unless `image` is a
bijection of `0..image.len()`" is incomplete: a valid length-7 bijection
panics via `ImageTooLong`, because `try_from` rejects the length before
bijectivity. All four siblings (`between`, `between_all`, `unrank`,
`from_cycles`) state the fixed-maximum panic. Defense: the type doc directly
above states the domain. Verdict: a misleading affirmative condition on the
designated construction path for external fixed tables beats implicit
context.

**F5 — `ClassKey::space` is documented as total but panics on parse-reachable
keys, permanently poisoning the registry** (`class.rs:138-148`; reduced from
the fallibility review's F1, absorbing the documentation review's F4).
`"Cyc0".parse::<ClassKey>()` succeeds (`from_str` checks only
`degree > MAX_DEGREE`), and `.space()` then asserts inside the held
`REGISTRY` lock; every later `space()` call on any key panics
"coset-space registry poisoned" — an undocumented process-wide consequence.
The same undocumented panic reaches any directly constructed family key
outside the fixed domain. Defense: the asserted contracts and the poisoning
assertion are settled behavior (157 S4b; poisoning recorded in doc 119 as out
of scope). Verdict: the defense covers the behavior, not the silence; the
crate's own convention is to name panic conditions on the invoked method.

**F6 — the oriented wrapper states no inherited contracts**
(`oriented.rs:61-92`, `oriented.rs:110-118`; reduced from the fallibility
review's F1). `OrientedPermutation::{identity, apply, compose}` and both
`generate` operations inherit point-domain, degree-mismatch, and
maximum-degree panics; 157 settled "the same contracts as their underlying
permutation operations", but that statement appears in no rustdoc. Defense:
the module docs state the crate-wide bound. Verdict: the inherited-contract
statement exists nowhere a caller of the oriented API can see it.

### Contract documentation

**F7 — no `# Errors`/`# Panics` sections anywhere in the crate; contracts are
prose-only or absent** (crate-wide; confirmed). Zero contract headings exist
in umol-perm, against 14 `# Panics` and 12 `# Errors` in umol-graph-ir at the
same commit; the `TryFrom<&[usize]>` impl is entirely undocumented and
`FromStr` lacks `# Errors`. The committed data-type-contracts skill requires
the headings "exactly when applicable", so the finding does not rest on the
uncommitted guides. The fallibility review's remaining per-operation items
(the maximum-degree panics of `identity`, `symmetric`, `alternating`,
`cyclic`, `dihedral`, `generate`) reduce into this finding: the domain is
stated in module and type prose but the per-operation entries are missing.
Defense: the prose is accurate for most operations, and conversion is
mechanical. Verdict: confirmed as the compliance criterion; ordinary `None`
conditions correctly remain in prose.

**F8 — `# Semantic properties` sections are absent crate-wide, with three
concrete unstated laws** (reduced). The crate-wide absence is the known
migration state the property-tests guide acknowledges. Three tested laws are
stated nowhere in any form: the `enantiomer` involution (which requires
improper² ∈ R); the `TryFrom<&[usize]>` contract (five pinned outcomes, no
doc at all); the `ClassKey` `Display`/`FromStr` roundtrip (both impls
undocumented). Claimed orphans that were dropped: the `index`/`unindex` and
`reindex` identity laws are derivable from the operations' own prose, and
`order == elements().len()` is definitional.

**F9 — the degree-bound prose contradicts its own single-named-home promise**
(`permutation.rs:1`, `group.rs:1` vs `permutation.rs:15-18`; reduced).
`MAX_DEGREE`'s doc promises raising it is "a one-line change here" while two
module docs hardcode "6"; doc 119 §F's "fix prose" item is not in its
Implemented record. Defense: both statements are currently true, and a
concrete 6 reads faster. Verdict: the narrow inconsistency survives; resolve
by rewording one side consistently.

### Test witnesses

**F10 — the five settled `Option` failure contracts of the coset lookup
family have no witnesses, and doc 119's record of them is inaccurate**
(`coset.rs` `coset_rep`/`index`/`unindex`/`reindex`/`enantiomer`; confirmed,
strengthened). No test in the crate asserts a `None` result for any of the
five; `reindex`'s parent-membership guard (`coset.rs:119`) is executed by no
test. Doc 119 records "umol-perm unit cases for the new `None` paths" as
landed; git history shows they never existed at any commit — the
consolidation defense is disproven, and doc 119's Implemented record is
inaccurate on this point. Defense (surviving context): generated `None`
coverage was consciously placed on `observable_coset`/`orbit_reps`, which
have exact `None` tables.

**F11 — the improper-query-on-improper-free-group branch of `contains` is
untested** (`oriented.rs:168`; reduced). The `None => false` arm — a group
with no mirror component contains no improper operation — has zero
executions; the unit table's group has an improper part and the property
groups are only queried with their own elements. The S2c-violation framing
was dropped (S2c's text does not name this case); the surviving basis is the
skill's failure-boundary rule. A one-row table addition closes it.

**F12 — family `ClassKey` spaces and `Coset` equality have no executable
witnesses** (`class.rs:38-57`, `class.rs:130-133`, `class.rs:262-264`;
reduced). The four family arms of `build()` are executed by no test anywhere;
no `.space()` call on a family key exists in the workspace; nothing defers
family-space tests. The hand-written pointer-identity `PartialEq for Coset`
(cross-class inequality at equal index) has no witness, and `umol_perm::Coset`
has zero importers outside the crate, so the transitively-exercised defense
is factually wrong. Dropped: the 80 %-baseline framing (the baseline is met)
and the `new_unchecked` item (testing a documented bypass's non-check favors
the original).

**F13 — `coset_rep`'s stated min-rank canonicality has no non-circular
witness** (`coset.rs:94-97`; reduced from the documentation review's F7).
The numbering is keyed by `coset_rep` itself, so a consistently wrong
canonical choice still passes every roundtrip; nothing asserts the
minimum-Lehmer-rank property. The companion claims were dismissed: `elements`
sortedness is load-bearing for `contains`' binary search (indirect but
genuine validation), and `improper_rep` minimality has a direct exact
witness.

### Naming

Naming replacements below are candidates recorded from the reviews; the final
names require consultation and are not settled by this document.

**F14 — property test `test_coset_space` states no law and its stem collides
with the `CosetSpace` unit tests** (`tests/property.rs:339`; confirmed,
found independently by two reviews). The test asserts the S4b interning law
(`ptr::eq(coset.space(), key.space())`), which its sibling
`test_class_key_space_interning` names for the identical law; seven
`test_coset_space_*` unit names denote `CosetSpace` methods, making the bare
stem read as a `CosetSpace` test. One law currently carries two naming shapes
in one file.

**F15 — the class-key strategies misstate or omit their emitted domains**
(`tests/property.rs:84-108`; confirmed). `class_key_text()` emits `ClassKey`
values across all ten families — "text" names its consumer, not its domain —
and `class_key()` silently narrows every space/coset property to the six
geometry classes; neither strategy carries a doc comment. The settled S7c
scope defends the narrowing, not its invisibility. Recorded candidates:
`class_key_text` → `class_key` (full domain) and `class_key` →
`geometry_class_key`.

**F16 — the property suite's definition-level references are unexplained**
(`tests/property.rs`; reduced). The full-`Sn` enumeration reference in
`test_permutation_between_all` and the BFS orbit-traversal reference in
`test_coset_space_orbit_reps` carry no role comments, and the
`permutation_group` degree cap is uncommented while the identical
`oriented_group` cap is commented. The broader claim against the one-line
module doc did not survive (flat suite, per-strategy docs).

## Refuted findings

Recorded so they are not re-flagged.

- **Ordering inside `class.rs`** (public `space` after the private `build`
  tables): the guide's "focus types and functions first" governs module-level
  grouping, not method order inside one impl block; the focus types do lead
  their groups, and reviews are functional, not aesthetic. Ties favor the
  original.
- **"slots beyond `degree` are identity"** (`permutation.rs:4`): the sentence
  describes the backing `[u32; MAX_DEGREE]` storage cells, which `try_from`
  identity-fills and `inverse`/`compose` rely on; `apply` panics for points
  beyond degree, so the proposed rewrite ("points beyond degree are fixed")
  would state a false public contract. Storage cells and domain points are
  distinct concepts; no synonym violation.

## Open items

1. **`Cyc0`/`Dih0` boundary.** `ClassKey::from_str` accepts degree 0 for
   cyclic/dihedral, but the resulting key's `space()` panics and poisons the
   registry (F5). Reject at parse (extending `ParseClassKeyError`), constrain
   construction, or document only — the doc contract cannot be written
   cleanly until this is decided.
2. **`MAX_DEGREE` visibility.** Doc 119 settled `pub(crate)`; evidence
   postdating that decision cuts against it: the property suite hardcodes 6
   in seven strategies, public rustdoc names the constant, and the value
   leaks through the error variants' `maximum` fields.
3. **`Coset`'s public surface.** Zero importers outside umol-perm; whether an
   unconsumed producer-facing type keeps two public construction routes.
4. **`Coset` placement.** The type lives in `class.rs` while the `coset`
   module holds `CosetSpace`; the layering defense (the registry dependency)
   is recorded, the reader-navigation asymmetry remains.
5. **OpenSMILES verification cluster.** The TH/@AL index-to-arrangement
   correspondence is pinned by no umol-perm test; the TB/OH enantiomer
   @↔@@ pairing remains the known open item from doc 119; the
   fixed-point-freeness of TB/OH axial swaps — design data for doc 103's `~`
   operator — is unpinned (an exact scan during review confirmed it holds:
   moved cosets TH 2/2, AX 2/2, TB 20/20, OH 30/30; CT 0/2, SP 0/3). Decide
   whether to pin now or together with the pairing verification.
6. **`reindex` right-action compatibility.**
   `reindex(i, g∘h) == reindex(reindex(i, g), h)` — the substrate of doc
   103's `^`-composition rule — is untested in umol-perm; decide whether
   downstream coverage suffices.
7. **Per-assert panic pinning.** `apply`/`compose`/`identity`/`unrank`,
   `from_cycles`' degree assert, and `Coset::new`'s range assert have no
   `#[should_panic]` cases while the 119/157-enumerated set does; decide
   whether the current selection is the settled scope.
8. **Law-stating property names.** Six property names state only the
   operation (`test_permutation_between`, `test_coset_space_orbit_reps`, …);
   no guide entry adjudicates whether the operation name suffices.
9. **Strategy stem spelling.** `perm_of`/`perm_pair`/`perm_triple` vs
   `permutation()` in one file; the complete-words rule is scoped to public
   identifiers, so unification direction is open.
10. **Module-naming rule's normative home.** Doc 118's noun/verb module rule
    was never migrated into a living guide; decide whether to record it in
    the nomenclature guide, and whether the private `oriented` module becomes
    a noun under it.
11. **Doc 119 record correction.** F10 shows doc 119's Implemented record is
    inaccurate for the `None`-path unit cases; a dated addendum to doc 119 is
    the appropriate correction vehicle.
12. **Governing-text provenance.** The code-reviews guide, the data-types
    contract section, and the review-cycle skill were uncommitted
    working-tree state at the pinned commit; they should be committed so
    future cycles review against tree-recorded standards.

## Proposed design

No staged plan here; stages follow after triage.

### Contract sections for umol-perm

One documentation pass over the crate, applying the data-types contract
sections:

- Convert every prose panic/error contract to per-operation `# Panics` /
  `# Errors` entries naming the property (F7), including the currently
  undocumented `TryFrom<&[usize]>` and the `# Errors`-less `FromStr`.
- Correct the inaccurate texts while converting: `|P| / |R|` (F1), `R\P`
  (F2), dihedral order claim scoped to degree ≥ 3 (F3), `from_image`
  maximum-length condition (F4), `ClassKey::space` panic domain and the
  registry-poisoning consequence (F5), and the oriented wrapper's inherited
  contracts (F6).
- Add `# Semantic properties` for the three unstated laws (F8): enantiomer
  involution, `TryFrom` contract, `ClassKey` text roundtrip.
- Add the crate's two genuine `# Assumes` cases: `Coset::new_unchecked`
  (`index < count`, established by the caller) and `CosetSpace::new`
  (`improper ∈ parent`, unchecked in the body, established only by
  `ClassKey::build`'s fixed tables — `enantiomer`, `is_chiral`, and
  `observable_coset` correctness silently depend on it).
- Add `# Establishes` entries per the review inventory; the non-obvious ones:
  `identity`/`from_image`/`TryFrom` establish the identity-filled
  representation tail that derived `Eq`/`Ord`/`Hash` rely on; `generate`
  establishes Lehmer-sorted closure; `cycles` establishes the canonical
  decomposition and its `from_cycles` roundtrip; `coset_rep` establishes
  minimum-Lehmer-rank representatives; `unindex`/`reindex`/`enantiomer`
  establish their roundtrip/identity/involution laws; interning establishes
  pointer-stable `&'static` spaces; `OrientedPermutation::inverse`
  establishes orientation preservation; `orbit_reps` establishes the
  minimum-index representative map (idempotent, `result[i] <= i`).
- Reword the degree-bound prose consistently with `MAX_DEGREE`'s
  single-named-home promise (F9).

### Test additions

- Unit `None` cases for `coset_rep`, `index`, `unindex`, `reindex` (both
  `None` causes, including the currently unexecuted parent-membership
  guard), and `enantiomer` (F10).
- One `contains` table row: improper query against an improper-free group
  (F11).
- Family-space witnesses: construct `Symmetric`/`Alternating`/`Cyclic`/
  `Dihedral` spaces once each and pin count/degree basics (F12), unless
  triage records an explicit deferral instead.
- A `Coset` equality case: cross-class inequality at equal index (F12).
- A non-circular `coset_rep` minimality witness: independent minimum over
  `Rσ` by Lehmer rank (F13).
- Role comments for the two definition-level references and the domain
  comment for the uncommented degree cap (F16).

### Naming changes (pending consultation)

- `test_coset_space` → a law-stating name (F14).
- `class_key_text` / `class_key` → domain-stating names (F15); strategy doc
  comments stating emitted domains either way.

### Records

- Dated addendum to doc 119 correcting the `None`-path unit-case record
  (open item 11).
- Commit the governing normative texts (open item 12).

## Second cycle: delta at commit a52651ee (2026-08-31)

This cycle reviews only the umol-perm surfaces changed since the first
cycle's pin (`836f7402..a52651ee`): the new `dynamic` module
(`DynPermutation`), the `CosetSpace` additions (`is_partitioned`, `allows`,
`normalizer`, the `observable_coset` parameter rename), the `MAX_DEGREE`
publication, the `between_all` removal, and the associated tests. Two
merged-area review agents and one refutation agent ran under the same
stance. Verdict tally: six confirmed, one reduced, one open question, one
refuted. Check placement in the delta is again fully compliant — two
class-1 checked constructors, four class-2 `Option` guards, one internal
`debug_assert` with an established producer, zero defensive or premature
checks.

### First-cycle conclusions rechecked

Every first-cycle finding on unchanged code stands; oriented.rs and
group.rs were not touched. Deltas:

- **F9 resolved.** Doc 214 S0a published `MAX_DEGREE` with a rewritten doc
  and removed the inaccurate "one-line change here" promise — the surviving
  core of F9. The module docs' concrete "6" remains as plain statement.
- **Open item 2 resolved** by the same settled subitem (214 S0a), with a
  dedicated integration test using the symbol.
- **F16 reduced.** The `between_all` removal (settled by 214 S0p) took the
  full-`Sn` enumeration reference and `repeated_orderings` with it; only
  the BFS orbit reference remains unexplained. Open item 8's list shrinks
  by `test_permutation_between_all`.
- **Inventory maintenance.** The `between_all` entries in the proposed
  contract sections are moot; `CosetSpace::new`'s `# Assumes` set gains a
  second member (D4 below); the delta's new operations join the F7/F8
  documentation pass (D3).
- **Open item 10 recontoured.** Adjective names for data modules are
  settled acceptable when they are the most concise description
  (`oriented`, `dynamic`); the remaining question is only recording the
  operation-module rule (doc 118) in the nomenclature guide.
- **Open item 12 closed.** The governing texts are committed on this
  branch ("Update code review guidance").

### Second-cycle findings

**D1 — `observable_coset` doc asserts a one-sided generated supergroup;
the implementation computes double cosets** (`coset.rs:208-215`;
confirmed). The delta replaced "under a fluxional supergroup" with "the
supergroup generated by the proper-rotation group and
`supergroup_generators`" — a checkable group-theoretic claim that is
false: cosets are left-defined (`Rσ`) while the generators act on cosets
by right multiplication, so the classes are the right-`⟨generators⟩`
orbits on cosets (double cosets `R\P/⟨generators⟩`). Reproduction
(independently confirmed by the refutation agent): trigonal-bipyramidal
with one transposition yields 10 observable classes; the claimed generated
group has order 120 and would admit exactly one. Defense: `orbit_reps`
directly above states the right action correctly and the sole consumer's
semantics are unaffected — the code is right, the sentence is wrong.
Fix: reword on `orbit_reps`' model.

**D2 — `DynPermutation::degree` doc says "the number of positions moved"**
(`dynamic.rs:41-44`; confirmed). `degree()` returns the frame's position
count; `identity(8).degree() == 8` with zero moved points, and "moves" is
the crate's own term for non-fixed points. One-word class fix.

**D3 — both `TryFrom` impls for `DynPermutation` carry no rustdoc and no
`# Errors` entries** (`dynamic.rs:105-123`; confirmed). New post-adoption
code — the migration shelter does not apply; the two named failure
properties are tested but stated nowhere. The remedy may land inside the
F7 pass, but the defect is the delta's.

**D4 — the `is_partitioned` producer contract is stated nowhere**
(`coset.rs:32,40-46,160-172`; producer `class.rs:37`; confirmed). The flag
⟺ degree-4 two-block parent coupling is established only by
`ClassKey::build`'s tables, is load-bearing for `normalizer`'s published
`allows`-membership law, and is guarded only by a degree-4 `debug_assert`.
Joins `improper ∈ parent` as the constructor's second genuine `# Assumes`
entry.

**D5 — a private helper precedes the focus type in `dynamic.rs`**
(`dynamic.rs:5-23`; confirmed). `validate_image` opens the module;
the first cycle's refuted ordering finding located the rule at module
level, which is exactly this case. Move the helper below the impl, beside
its two callers. (The `Decomposition` lead-in in coset.rs is
distinguishable: a constructor-parameter type on the construction
surface.)

**D6 — the three new coset tests were prepended out of declaration order**
(`coset.rs:376-482`; confirmed). The suite parallels method declaration
order; `allows`/`normalizer` tests belong between `is_chiral`'s and
`orbit_reps`' tests per the test-writing skill's ordering rule.

**D7 — the slice `TryFrom` has no in-crate witness** (`dynamic.rs:104-113`;
reduced). All umol-perm tests and the strategy dispatch to the `Vec` impl;
the slice impl's only executions are downstream success-path literals.
Reduced to a minor completeness item (both impls share `validate_image`),
closed by one table row; not a settled-scope violation.

### Refuted (recorded)

- **`normalizer` property-witness demand**: doc 209 S1a — the governing
  subitem for `allows`/`normalizer`, settled and implemented line by line
  — prescribes exactly the delivered exact-table evidence, and the
  generative selection law is deliberately placed downstream in the
  graph-IR frame-selection properties (211/214; property-tests.md).
  Deliberate settled scope is not a flaw. Both review agents missed 209
  because the index marks it Superseded while its S0–S1 remain
  settled-implemented.

### New open items

13. **`normalizer`/`allows` guide registration.** The names are settled
    (doc 211), but three tensions are unadjudicated: the Normalize
    glossary scopes *normalize* to frame-preserving operations while this
    operation selects a frame; the agent-noun reservation; and the
    group-theoretic meaning of "normalizer" (`N_P(R)`) adjacent to this
    exact type. Register, qualify, or migrate — consultation needed.
14. **`supergroup_generators` naming.** The rename from `fluxional`
    matches the actual consumer but has no recorded consultation, and if
    D1's rewording lands, the parameter name should be revisited with it —
    no settled repo term exists for "generators of the right-acting merge
    group over cosets."
15. **`DynPermutation` cross-degree `Ord`.** The derived ordering compares
    across degrees by lexicographic prefix; no consumer orders these
    values and no law states the semantics. Stated contract or documented
    unspecified — undecided.
16. **Strategy literals vs `MAX_DEGREE`.** Eight bound-shaped literals in
    tests/property.rs could now spell the public constant (completing
    214 S0a's anti-duplication rationale locally) or remain deliberate;
    the `2..=4` tractability caps stay literal either way.
17. **Stated home for `DynPermutation`'s tested laws.** Whether the F8
    derivability precedent covers the dyn carrier's group laws or the
    documentation pass adds a `# Semantic properties` section.

### Process notes

- Two findings cited normative text that does not exist (a nomenclature
  "degree entry"; a literal "every trait impl tested" rule) and survived
  only on other grounds — citation verification by the refutation pass is
  earning its cost.
- A Superseded status can hide settled scope from reviewers who follow the
  index: doc 209's S0–S1 govern part of this delta. Review agents should
  read a superseded document's status note before excluding it.

### Proposed design (delta additions)

- Reword `observable_coset` on the right-action model (D1) and fix
  `degree`'s doc (D2) during the contract pass; add `# Errors` to both
  `TryFrom` impls (D3) and the `is_partitioned` producer contract to
  `CosetSpace::new`'s `# Assumes` (D4).
- Move `validate_image` below the impl (D5); reorder the three coset tests
  (D6); add one slice-`TryFrom` row (D7).
- Naming consultations per open items 13–14.

## Triage record (2026-08-31)

- F1–F4: doc corrections applied.
- F5 and open item 1 settled and implemented: parsing rejects degrees outside
  the family's supported domain (`InvalidDegree`, widened to cover
  out-of-domain degrees). The domain is single-sourced in
  `ClassKey::degree_domain`, uniformly `1..=MAX_DEGREE` for the four families
  — the zero row was dropped entirely: `Cyc0`/`Dih0` name no group, and
  `Sym0`/`Alt0` are unconsumed induction roots. Direct construction retains
  its panic, now documented in `# Panics` on `space` together with the
  registry-poisoning consequence. An exhaustive enumeration test over every
  representable key asserts parse ⟺ domain ⟺ buildable agreement — a total
  check of the finite key space, not a sample — and incidentally supplies
  F12's family-arm witnesses. `build` deliberately gains no entry assertion:
  the callee asserts already enforce the domain, and a duplicate would be a
  class-3 check.
- F6–F8 implemented as one contract-section pass, absorbing D3 and D4:
  `# Panics`/`# Errors` headings on every panicking or `Result`-returning
  operation in the crate, including the previously undocumented `TryFrom`
  impls of both permutation carriers; the oriented wrapper's inherited
  contracts stated per operation (F6); `# Semantic properties` for the three
  unstated laws (F8: the `enantiomer` involution, the `Permutation` one-line
  construction inverse, the `ClassKey` text roundtrip); and the two genuine
  `# Assumes` cases (`Coset::new_unchecked`; `CosetSpace::new`'s
  `improper ∈ parent` together with the `is_partitioned` producer contract,
  closing D4). The `# Establishes` entries remain the outstanding piece of
  the contract pass.
- F9 closed as resolved: the "one-line change" sentence was self-talk, not a
  statement of fact, and doc 214's rewrite removed it; nothing further.
- F10–F13 implemented: `None` witnesses for all five coset lookups —
  `coset_rep` (wrong degree), `index` (wrong degree; outside the restricted
  parent), `unindex` (out of range), `reindex` (both causes, including the
  previously unexecuted parent-membership guard), `enantiomer` (out of
  range); the improper-query-on-proper-group `contains` row (F11); a `Coset`
  equality table with cross-class inequality at equal index (F12's residue —
  the family-arm witnesses were already supplied by the degree-domain
  enumeration); and a non-circular `coset_rep` minimality witness (F13): a
  definition-level reference computing the minimum over `Rσ` exhaustively
  across every permutation of three spaces' degrees.
- `# Establishes` pass settled and implemented under the aggressive discharge
  reading. The disambiguation is recorded in data-types.md: `# Semantic
  properties` states laws of the operation algebra; `# Assumes`/`# Establishes`
  state data properties downstream code relies on; a property established by
  every public producer of a type is stated once on the type. All three
  type-wide invariants were already stated in place (the identity-filled tail,
  the sorted element store, `Coset`'s pointer-plus-index identity). Entries
  added: `index` (a returned index is `< count`), `space` (pointer-canonical
  interning — the property `Coset` equality relies on), and `normalizer`'s
  `allows`-membership moved from `# Semantic properties` to `# Establishes`
  under the disambiguation. The code-reviews two-direction rule was relaxed to
  traceability: every tested law must be traceable to a stated source so a
  failing property identifies its semantic reason; heading ceremony is not
  required for laws derivable from the operations' prose.
- The `class_key_text` strategy and its roundtrip property were deleted as
  superseded: the exhaustive degree-domain enumeration asserts the same law
  over the complete key space rather than a sample. F15's `class_key_text`
  half dissolves with it; the residue is `class_key`'s undocumented
  six-geometry narrowing.
- Open item 15 resolved: `DynPermutation` ordering is degree-first, then
  lexicographic on the image — a manual `Ord`/`PartialOrd` replacing the
  derive (which interleaved degrees), stated on the type doc and pinned by a
  comparison table including the cross-degree case. No workspace consumer
  orders these values, so the semantics change is unobserved.
- F14–F16 and D1, D2, D5–D7 implemented, all decided by the recorded
  conventions without further consultation: the interning property renamed
  `test_coset_space_interning` (law-stating, per its sibling); the `class_key`
  strategy renamed `geometry_class_key` (name states the emitted domain);
  role comment on the BFS orbit reference and the missing degree-cap doc on
  `permutation_group` (F16); `observable_coset` reworded on the right-action
  double-coset model with the parameter name left to open item 14 (D1);
  `degree`'s doc corrected (D2); `validate_image` moved beside its `TryFrom`
  callers (D5); the three coset tests moved into declaration order (D6); the
  slice `TryFrom` witnessed per table row (D7).
- Open items 10 and 11 executed: the module-naming rule (result-noun/verb for
  operation modules; nouns or concise adjectives for data modules) is recorded
  in the nomenclature guide, and doc 119 carries the dated step-5 correction.
- Items 3, 4, 6, 7, and 9 dispositioned and executed: `Coset::new_unchecked`
  removed (no call sites; trivially re-added if a producer appears); `Coset`
  moved to `coset.rs` with the public path unchanged; the `reindex`
  composition law stated as a semantic property and pinned by a generated
  property over parent-element pairs; per-assert panic pinning declined —
  the 119/157 selection is the settled scope, the panic origins being
  obvious; strategy stems unified to `permutation_*`, and the `perm`/`perms`
  locals in umol-geometric renamed. Item 8 is settled by the traceability
  rule.
- Item 5 stays open as external-evidence work: verify the arrangement
  numbering against the OpenSMILES document in `materials/formats/opensmiles`
  and cross-validate with RDKit; natural trigger is the next umol-io
  conformance work.
- Item 16 executed: the domain-shaped strategy bounds spell `MAX_DEGREE`; the
  tractability caps remain deliberate literals.
- Item 13 executed: the nomenclature guide registers `Normalizer` — the
  returned action, reframe-side despite the stem, not the group-theoretic
  `N_P(R)`.
- Item 14 executed: `observable_coset`'s parameter is plain `generators`,
  matching `orbit_reps`; the corrected doc sentence carries the right-action
  semantics, so the wrong qualifier is deleted rather than replaced.
- Item 5 moved to doc 153 as task T9. With that, every finding and open item
  of both cycles is closed or moved, and this document is complete.
