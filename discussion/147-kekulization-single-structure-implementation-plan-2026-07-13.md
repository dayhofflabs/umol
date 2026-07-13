# Experiment B: single-structure kekulization implementation plan

## Objective and boundary

Implement deterministic single-structure kekulization for the closed-shell cases settled in
discussion 145:

- ordinary systems whose contribution vector is all `1`: perfect matching;
- one prescribed exposed `2e` donor or `0e` acceptor: delete that atom from the matching graph
  and require a perfect matching of the remainder;
- one delocalized charge on an otherwise all-`1` monoelement system: compute a general maximum
  matching with deficiency one and localize the system charge on its exposed atom;
- non-bipartite systems: use Edmonds, including azulene and odd charged rings;
- map every selected bond and exposed atom through the existing molecule correspondence rather than
  rebuilding host/subgraph maps or introducing another matching result carrier.

The first release explicitly rejects undetermined or malformed contribution vectors, contributions
outside `{0,1,2}`, open-shell aromatic systems, more than one prescribed hole, system-charge
magnitude greater than one, and mixed prescribed/mobile demands. It implements one deterministic
structure only. Enumeration, symmetry quotienting, weighted/HMO selection, arbitrary
degree-constrained matching, and multiple mobile holes remain outside Experiment B.

The four open questions at the end of discussion 145 are therefore resolved for this implementation
boundary as follows:

1. charge magnitude greater than one is unsupported rather than interpreted as a hole count;
2. the exposed site and matching follow the caller-supplied canonical atom order, with no chemical
   score;
3. Experiment B has no enumeration API;
4. symmetry reduction is irrelevant to a single deterministic result and remains separate.

The localized output preserves atom-local lone-pair and charge information already present for
heterogeneous prescribed holes. Charge equalization has moved only monoelement delocalized charge
onto the aromatic-system record; kekulization reverses that move by adding the system charge to the
chosen exposed atom. It does not synthesize a lone pair merely from contribution `2`.

One correction to the shorthand Experiment B list in discussion 145 is required by the existing
aromaticity contract: neutral borepin boron contributes `0` and is a prescribed hole, whereas
boratabenzene anion retains atom-local `B-` and contribution `1` because heterogeneous systems
are not charge-equalized. Boratabenzene is therefore a charged, no-hole preservation case rather
than a second 0e-hole case.

All data-backed chemical cases live in a dedicated `umol-graph` integration-test directory and use
the existing molecule EDN syntax. Every edited test follows the test-writing conventions; touched
test modules are normalized rather than extending local deviations.

## S0 — Demand vocabulary and trusted inputs

### S0a — Add the internal matching-demand classifier

**Module:** `umol-graph/src/ops/transform/kekulizer.rs`

**Kind:** additive (green)

**Dependencies:** [dep: none]

Add a private `MatchingInput` carrying required-covered atoms, required-exposed atoms, and the exact
exposed count. Add a private demand mode/classification used only to dispatch the two supported
solutions:

- prescribed mode, where the required-exposed set is the complete exposed set;
- one-mobile-hole mode, where no site is prescribed and the exact exposed count is one.

Derive it from each aromatic system's positional `ElectronCountsAst`, literal charge, and spin.
In prescribed mode, contributions `1` become required-covered; `0` and `2` identify a prescribed
exposed atom while retaining whether its local π state is acceptor or donor. In the all-`1`,
charge-`±1` mobile case, all atoms remain flexible and the exact exposed count is one. An all-`1`,
zero-charge system is prescribed mode with every atom required-covered and exposed count zero.

Extend `KekulizerError` with specific, system-identified variants for:

- undetermined electrons, charge, or spin;
- electron-vector/member-count mismatch;
- unsupported contribution value;
- unsupported open shell;
- unsupported multiple prescribed holes;
- unsupported system-charge magnitude;
- mixed prescribed and mobile demand.

Classification is read-only and must not partially mutate the molecule on error. Unit table tests
cover every supported class and every error variant, including atom positions in positional vectors.

**Verification:** focused kekulizer unit tests, graph crate check, and Clippy.

### S0b — Establish the Experiment B integration corpus

**Module:** `umol-graph/tests/kekulization.rs`,
`umol-graph/tests/kekulization/mod.rs`, and
`umol-graph/tests/kekulization/data/*.edn`

**Kind:** additive test foundation (green)

**Dependencies:** [dep: S0a]

Add locally owned, unambiguous EDN inputs for:

- benzene and pyridine;
- pyrrole, furan, and thiophene;
- borepin and boratabenzene;
- cyclopentadienyl anion and tropylium;
- azulene;
- one fused heterocycle containing a prescribed donor and ordinary covered atoms.

At this stage the integration harness verifies that each fixture parses as a ground molecule and
that its aromatic-system participants, positional contributions, charge, and spin encode the
intended demand. It does not add failing expectations for unimplemented kekulization. Keep expected
localized outputs adjacent and clearly suffixed only when the corresponding implementation subitem
lands.

**Verification:** the new integration target compiles and its fixture-structure table passes.

**Stage exit:** every supported and rejected demand has a precise internal representation, and the
chemical corpus independently confirms the intended inputs.

## S1 — Correspondence-native general matching

### S1a — Replace manual extraction maps with molecule correspondences

**Module:** `umol-graph/src/ops/transform/kekulizer.rs`

**Kind:** behavior-preserving refactor (green)

**Dependencies:** [dep: S0a]

Refactor `plan_systems` to use the `MoleculeCorrespondence` returned by
`MoleculeAst::induced_subgraph`:

- order the system atoms by filtering the caller-supplied canonical `node_order`;
- materialize the extracted molecule through that correspondence;
- map selected subgraph bonds with `correspondence.bonds().right_of(...)`;
- map subgraph atoms with `correspondence.atoms().right_of(...)`.

Remove the manual `HashMap`, sorted-host reconstruction, and positional assumptions about
subgraph bond IDs. Validate that the filtered canonical order contains every system atom exactly
once; report a specific `InvalidNodeOrder` error for missing or duplicate participants while
allowing unrelated atoms from other systems in the global order.

Do not add `Matching → Correspondence` conversion or a public matching-transport API. The caller
needs selected edges and exposed vertices in the host molecule, both already supplied by the
existing correspondence.

Tests prove behavior preservation on benzene and on multiple disjoint aromatic systems, plus
non-identity host/subgraph atom and bond ID mappings.

**Verification:** kekulizer unit tests and graph correspondence tests remain green.

### S1b — Migrate the kekulizer model from perfect DFS to general Edmonds

**Module:** `umol-graph/src/ops/transform/kekulizer.rs`,
`umol-graph/src/ops/transform.rs`, and all `KekulizationModel` callers

**Kind:** breaking public model-field migration (red→green within S1)

**Dependencies:** [dep: S1a]

Change `KekulizationModel::algorithm` from `PerfectMatchingAlgorithm` to
`MaxMatchingAlgorithm`, with `Edmonds` as the default. Migrate the constructor, exports, tests,
and all workspace callers in the same subitem.

For the zero-hole path, compute a maximum matching and accept it only when
`matching.is_perfect(graph.node_count())`; otherwise preserve the system-specific no-matching
error. The caller-supplied canonical atom order controls extracted node IDs and therefore the
deterministic Edmonds result. Keep `PerfectMatchingAlgorithm` and `Graph::perfect_matching` in
graph-core for unrelated callers; Experiment B does not retire that API.

Add cases for ordinary benzene/pyridine, atom-locally charged boratabenzene, and non-bipartite
azulene. Assert exact host bond IDs for a fixed canonical order, matching validity, deterministic
repeatability, and preservation of the boron-localized charge.

**Verification:** graph-core tests, kekulizer tests, workspace check, and Clippy.

**Stage exit:** the existing no-hole functionality is correspondence-native, deterministic, and
uses the general-graph solver required by the remaining Experiment B cases.

## S2 — Prescribed-hole kekulization

### S2a — Plan a perfect matching after prescribed-hole deletion

**Module:** `umol-graph/src/ops/transform/kekulizer.rs`

**Kind:** additive (green)

**Dependencies:** [dep: S0a, S1b]

For prescribed mode:

1. remove the single required-exposed atom from the canonically ordered system atom list;
2. obtain a residual `MoleculeCorrespondence` directly from the host AST;
3. extract the residual molecule and run Edmonds;
4. require a perfect matching of every residual atom;
5. map residual matched bonds back to host `BondId` values through the correspondence;
6. record the prescribed host atom as exposed in `SystemPlan`.

Verify that every required-covered atom is matched and the exposed set is exact. Distinguish an
unsupported demand from a supported demand with no feasible matching.

Unit tests cover five- and seven-membered residual paths, a fused heterocycle, nontrivial host IDs,
and an impossible prescribed hole.

**Verification:** focused unit tests plus graph-core matching tests.

### S2b — Apply prescribed donor/acceptor plans

**Module:** `umol-graph/src/ops/transform/kekulizer.rs` and
`umol-graph/tests/kekulization/`

**Kind:** additive behavior (green)

**Dependencies:** [dep: S2a, S0b]

Extend the application pass so prescribed holes receive no double bond while all residual matched
bonds become double and all other aromatic-system bonds become single. Remove aromatic bond and atom
constraints and the aromatic-system record as before.

Preserve the exposed atom's existing local charge, lone pairs, hydrogens, and spin. Contribution
`2` is evidence that the existing local donor state supplies the π pair; contribution `0` is an
existing acceptor state. Do not infer or overwrite those fields from element identity.

Add exact localized-output integration cases for pyrrole, furan, thiophene, borepin, and the fused
heterocycle. Assert bond orders, exposed atom identity, preserved atom fields, removed aromatic
constraints/system records, total charge, and valence electron accounting.

**Verification:** focused integration target, graph crate tests, and Clippy.

**Stage exit:** fixed 0e/2e heteroatom holes kekulize without weakening matching correctness or
discarding their pre-existing local chemistry.

## S3 — Mobile charged holes and localized validation

### S3a — Plan a deterministic one-mobile-hole maximum matching

**Module:** `umol-graph/src/ops/transform/kekulizer.rs`

**Kind:** additive (green)

**Dependencies:** [dep: S0a, S1b]

For one-mobile-hole mode, run Edmonds on the full canonically ordered system graph. Require
`node_count - 2 * matching.size() == 1`; a different deficiency is a specific
`MatchingDeficiency` error. Determine the unique exposed subgraph atom with
`Matching::is_matched`, map it to its host `AtomId`, and include it in `SystemPlan`.

Tests use multiple canonical atom orders for odd cycles to prove that:

- each result is maximum and has exactly one exposed atom;
- repeating an order is deterministic;
- changing the canonical order may choose a different, but correspondingly canonical, localization
  site;
- the planner does not enumerate alternative sites.

**Verification:** focused planner tests and graph-core Edmonds tests.

### S3b — Localize system charge transactionally

**Module:** `umol-graph/src/ops/transform/kekulizer.rs` and
`umol-graph/tests/kekulization/`

**Kind:** additive behavior (green)

**Dependencies:** [dep: S3a, S0b]

For a mobile charged plan, add the literal system charge (`-1` or `+1`) to the exposed atom's
literal charge, then remove the system record. Preserve every other atom field. Apply all plans to
a candidate clone rather than the caller's AST, run model-independent entity-structure, valence,
and spin invariant validation on the localized candidate, and replace the caller's AST only after
all systems and validations succeed.

Add `KekulizerError` variants that preserve the specific post-localization structural or invariant
failure. Any planning, localization, or validation error must leave the original AST unchanged.

Integration cases cover cyclopentadienyl anion and tropylium. Assert the canonical exposed site,
localized charge sign, exact bond orders, total charge conservation, closed-shell state, removal of
the system, and input immutability on failure.

**Verification:** focused integration target, validator tests, and workspace check.

### S3c — Complete inverse π-electron accounting for localized holes

**Module:** `umol-graph/src/ops/transform/aromatizer.rs` and its tests

**Kind:** additive correction (green)

**Dependencies:** [dep: S2b, S3b]

Extend `electrons_from_kekule` for the localized states introduced or exposed by Experiment B:

- exposed neutral aromatic boron: `0`;
- exposed `C+`: `0`;
- exposed `C-`: `2`;
- the existing N/O/S/Se/P/As donor cases remain `2`;
- exactly one incident double bond remains `1`.

Use literal element, charge, bond-order, and local-state preconditions; return `None` for ambiguous
or chemically unsupported states. Normalize the touched aromatizer tests to the test-writing
conventions, then add table cases for all new and existing branches.

Add integration round-trip accounting checks over every Experiment B fixture: compare total charge
and the reconstructed per-atom π contributions with the pre-kekulization aromatic-system data. Do
not require the aromatizer to restore the same aromatic-system identity or the same arbitrary charge
localization.

**Verification:** aromatizer unit tests, Experiment B integration tests, and aromaticity tests.

**Stage exit:** both prescribed and mobile one-hole systems produce transactionally validated,
charge-conserving localized molecules with a tested inverse electron-accounting path.

## S4 — Error contract and acceptance

### S4a — Complete rejected-case integration coverage

**Module:** `umol-graph/tests/kekulization/`

**Kind:** additive verification (green)

**Dependencies:** [dep: S2b, S3b]

Add dedicated EDN fixtures and a table of exact errors for:

- undetermined, length-mismatched, and out-of-domain electron vectors;
- open-shell aromatic systems;
- two prescribed holes;
- system charge magnitude greater than one;
- mixed prescribed and mobile demand;
- a supported prescribed demand whose residual graph has no perfect matching;
- a mobile demand whose maximum-matching deficiency is not one;
- incomplete or duplicate canonical node order.

Every error case asserts both the precise `KekulizerError` value and that the input molecule remains
structurally unchanged.

**Verification:** full kekulization integration target.

### S4b — Run the complete Experiment B conformance matrix

**Module:** `umol-graph/tests/kekulization/` and discussion 145 or an adjacent result note

**Kind:** additive verification (green)

**Dependencies:** [dep: S3c, S4a]

Run the supported corpus as one acceptance table covering:

- no holes: benzene, pyridine, and azulene;
- prescribed 2e holes: pyrrole, furan, thiophene, and the fused heterocycle;
- prescribed 0e holes: borepin;
- atom-local heterogeneous charge with no hole: boratabenzene;
- mobile charged holes: cyclopentadienyl anion and tropylium.

For every case assert matching validity, exact exposed count, required coverage, deterministic
repeatability, charge conservation, localized valence/spin invariants, absence of aromatic
constraints/system entries, and inverse π-electron accounting. Record the accepted boundary and any
fixture-specific findings without committing generated snapshots as the sole oracle.

Run `cargo fmt --all`, graph-core tests with `proptest`, focused graph tests, workspace tests with
`umol-py/.venv` activated, workspace Clippy over all targets with warnings denied, and
`git diff --check`.

**Stage exit:** Experiment B supports the agreed closed-shell single-hole chemistry and rejects every
out-of-scope demand explicitly and without partial mutation.

## S5 — Optional dispatch optimization

### S5a — Add a verified bipartite fast path

**Module:** `umol-graph-core/src/algorithms/matching.rs`,
`umol-graph/src/ops/transform/kekulizer.rs`, and matching/kekulization benchmarks

**Kind:** deferrable additive optimization (green)

**Dependencies:** [dep: S4b]

Only after the general Edmonds path is accepted, resolve the current `HopcroftKarp` name/complexity
contract and optionally dispatch bipartite residual graphs to it. Cross-check the chosen matching's
validity and chemical output against Edmonds on the full Experiment B corpus and benchmark enough
repeated small systems to justify the extra branch.

The core Experiment B deliverable is complete without S5.

## Critical path

`S0a → S1a → S1b → S2a → S2b → S3c → S4b` establishes prescribed-hole support.
`S0a → S1a → S1b → S3a → S3b → S3c → S4b` establishes mobile charged-hole support.
`S0a → S0b → S2b/S3b → S4a → S4b` is the data-backed conformance path.

S2 and S3a may proceed independently after S1b. S3c is the join point for the common inverse
electron-accounting contract. S5 is explicitly deferrable.

## Deferred work

- more than one prescribed or mobile hole;
- mixed fixed/flexible degree-constrained matching;
- charge magnitude greater than one;
- open-shell and excited-state kekulization;
- weighted/HMO bond and localization scoring;
- enumeration of hole placements or Kekulé structures;
- symmetry-inequivalent output;
- automatic FKT counting during single-structure kekulization;
- a public matching-demand or matching-mapping API without another concrete non-chemical caller.
