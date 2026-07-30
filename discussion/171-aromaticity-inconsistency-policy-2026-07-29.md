# 171 — Aromaticity resolution and conformance

Status: In Progress
Date: 2026-07-29
Relates: [166](166-molecule-ops-2026-07-27.md)

Aromatic input that the chemistry model declines currently leaves projections behind with no relation
to project from. This document states the problem, the resolver policy that follows from
kekulization being an explicit operation, and the related charge-delocalization and conformance work
that must land as one coherent unit.

## Symptom

`MoleculeAst.from_smiles("c1ccoc1", chemistry_model=<mdl>)` returns a molecule with **no** aromatic
system whose atoms still carry `#a` (`#a2` on the oxygen) and whose bonds still carry `#a`.
`reset_aromatic_valence = true` does not change the output. Verified 2026-07-29, Python bindings:

| SMILES | `daylight()` | `mdl()` |
| --- | --- | --- |
| furan `c1ccoc1` | 1 aromatic system | 0 |
| thiophene `c1ccsc1` | 1 | 0 |
| pyrrole `c1cc[nH]c1` | 1 | 0 |
| pyridine `c1ccncc1` | 1 | 1 |
| benzene `c1ccccc1` | 1 | 1 |

`mdl()` is `HueckelRule` over {C, N} with `min_ring_size: 6`, so furan fails on both the element scope
and the ring size.

Rendered output for furan under `mdl()`:

```clojure
{:atoms ["C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" … "O#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" …]
 :bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] …]}
```

## Why it is wrong

Per-atom `#a` (aromatic valence) and per-bond `#a` (aromatic incidence) are projections of an
aromatic-system relation onto its participants. With no system present they project from nothing. This
contradicts the entity model: a derivable constraint whose relation is absent has no derived side.

## Why the obvious repair does not work

Aromaticity is the only overlay carrying information that the localized bonds would otherwise have to
hold. In the output above every ring bond is order 1 and every ring carbon is `#v2 #h1 #a1`. Strip the
`#a` and each ring carbon has three bonds' worth of valence and needs four, so **reset alone yields an
invalid molecule**. Discarding aromaticity requires *supplying* the alternating bond orders it stood
in for.

Stereo does not have the same bond-order dependency. Removing an unrealizable stereo assertion leaves
the constitution and localized bond orders intact, so `Strip` is a sound explicit stereo-resolution
policy. A shared policy enum would therefore either expose an unsound aromatic option or remove a
sound stereo option.

## Constraint: perception is flag-driven

Aromaticity perception in the ingest path is driven by the input's per-atom aromatic flags, not by ring
analysis of a Kekulé structure. Verified: `C1=CC=CC=C1` yields **zero** aromatic systems under every
model, benzene included, while `c1ccccc1` yields one under both `daylight()` and `mdl()`.

Consequence: once the flags are discarded there is no flag-driven path that rediscovers the system.
This rules out `Strip`. It does not rule out `Keep`: retaining the assertions preserves the information
needed by a later operation or a different chemistry model.

## Policy

**Kekulization and aromatization are explicit operations, not part of any sanitization.** A resolver
must therefore not kekulize on the caller's behalf, which removes the only variant that could have
repaired the structure in place.

- **`Error`** — resolution returns `Contradictory` when the input asserts aromaticity the model
  declines. This is the default.
- **`Keep`** — retain the unmatched assertion and continue without materializing the declined
  relation. This is an explicit opt-in for callers that need the staging information.

`AromaticityInconsistencyPolicy` therefore has `Keep` and `Error`.
`StereoInconsistencyPolicy` retains `Keep`, `Strip`, and `Error`. Both are operational resolver
configuration and live beside `AromaticityResolveConfig` and `StereoResolveConfig`, respectively.
The current generic stereo `InconsistencyPolicy` is renamed rather than generalized.

Validators have no inconsistency policy. They validate the AST against the selected model and return
`Contradictory` for an asserted projection that is absent, unrealizable, or inconsistent with its
materialized relation.

`Kekulize` is an operation, not a policy, because it changes localized bond orders.

Raw parse/raise has no aromatic-system relation for `Kekulizer` to consume. A caller who wants furan
localized before applying `mdl()` must first resolve it under a model that accepts the source aromatic
form, then run the explicit transformation:

```text
parse/raise → resolve(source model) → kekulize → resolve/validate(target model)
```

The current Python surface does not expose this full chain.

## Derivation and verification

The consistency decision belongs with the derivation that has enough information to make it. Policy
is applied only by resolvers:

- `AromaticityPerception::derive` returns an `AromaticityDerivation` containing independently
  accepted systems and exact unmatched atom, bond, and existing-system IDs.
- `StereoPerception::derive` returns a `StereoDerivation` containing realizable atom and bond
  materializations and exact unrealizable or mismatched sites.
- The `*Perception` types carry the model and perform the operation. The `*Derivation` values are
  policy-free results consumed by resolvers and validators.
- Existing aromatic and stereo relations are compared with independently perceived results; their
  presence is not itself evidence of conformance.
- `AromaticityResolver` applies `Keep` or `Error` to unmatched aromatic projections before planning
  additions.
- `StereoResolver` applies `Keep`, `Strip`, or `Error` to unrealizable stereo projections before
  planning additions or removals.
- `AromaticityConformanceValidator` and `StereoConformanceValidator` consume the same policy-free
  derivation results but always report mismatches as contradictions.

Under aromatic `Keep`, an existing aromatic system rejected by the selected model remains unchanged;
resolution does not delete a materialized relation. `Error` reports it as a contradiction.

`Undetermined` constraints are vacuous and are not missing projections. Validators ignore an absent
or `Undetermined` projection. A non-vacuous assertion without a matching relation, or an explicit
projection that contradicts its relation, is `Contradictory`. A non-ground but non-vacuous assertion
that cannot yet be decided is `Underdetermined`.

`reset_aromatic_valence` clears the aromatic-valence constraints on atoms of newly materialized
systems by setting them to `Undetermined`. It is not the lever for a system the model declines.

### Public perception API

`AromaticityPerception` retains `find_systems` as the low-level algorithmic operation whose caller
supplies the π-electron source. `derive` is the standard AST-facing operation: it reads aromatic
assertions, calls `find_systems`, and compares the accepted systems with existing relations and
non-vacuous projections.

```rust
pub struct AromaticityDerivation {
    pub systems: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    pub mismatches: Vec<AromaticityMismatch>,
}

pub enum AromaticityMismatch {
    AtomProjection { atom: AtomId },
    BondProjection { bond: BondId },
    ExistingSystem { system: AromaticSystemId },
    ElectronContribution {
        system: AromaticSystemId,
        atom: AtomId,
    },
}

impl AromaticityPerception {
    pub fn new(model: &AromaticityModel) -> Self;

    pub fn find_systems<F>(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemAst)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(&AtomView<'_>) -> Option<u8>;

    pub fn derive(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
    ) -> Result<
        Solution<AromaticityDerivation, AromaticityContradiction>,
        AromaticityError,
    >;
}
```

`AromaticityResolver` and `AromaticityConformanceValidator` consume `derive`.
`Aromatizer` calls `find_systems` directly with its Kekulé electron source. The current mutating
`AromaticityPerception::add_systems` is removed; resolvers and transformers own materialization.

`StereoPerception::derive_stereo_atom` and `derive_stereo_bond` are the public per-entity operations.
The molecule-wide `derive` calls them for non-vacuous `#T` and `#C` assertions and compares their
results with existing stereo relations.

```rust
pub struct StereoPerception {
    model: StereoModel,
}

pub struct StereoDerivation {
    pub atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>,
    pub bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>,
    pub mismatches: Vec<StereoMismatch>,
}

pub enum StereoMismatch {
    UnrealizableAtom { atom: AtomId },
    UnrealizableBond { bond: BondId },
    AtomRelation { stereo_atom: StereoAtomId },
    BondRelation { stereo_bond: StereoBondId },
}

impl StereoPerception {
    pub fn new(model: &StereoModel) -> Self;

    pub fn derive(&self, ast: &MoleculeAst) -> StereoDerivation;

    pub fn derive_stereo_atom(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoAtomAst)>;

    pub fn derive_stereo_bond(
        &self,
        ast: &MoleculeAst,
        bond: BondId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoBondAst)>;
}
```

`StereoResolver` applies its policy to `StereoDerivation`; `StereoConformanceValidator` treats its
mismatches as contradictions before running the existing relation-shape and graph-symmetry checks.
The helper that derives the two ligands at one end of a cis-trans bond remains private.

## Resolver boundary

Charge delocalization moves out of `AromaticityResolver` into the explicit transformer specified
below, and `AromaticityResolveConfig::delocalize_charge` is retired in the same public-config
migration that adds the aromatic inconsistency policy.

The lasting resolver invariant is that resolution does not change a localized bond order, add or
remove an atom or localized bond, or otherwise alter the constitution. Kekulization and aromatization
change bond orders, which is why they remain caller-invoked transformations.

## Charge delocalization

Move charge delocalization out of `AromaticityResolver` and into a
model-independent `DelocalizeCharge` transformer under `umol-graph/src/ops/transform`. Retire
`AromaticityResolveConfig::delocalize_charge`.

Delocalization rewrites one resolved representation into another; it does not fill undetermined
state. For an aromatic system whose participants carry literal formal charges, the transformer moves
the summed charge onto the system, sets the contributing atoms to literal zero, and adjusts their
π-electron contributions so the total is preserved. It leaves undetermined or non-literal charges
unchanged and is idempotent.

No `LocalizeCharge` inverse is added: the choice of which atom receives the charge is not canonical,
so localization is not a function of the delocalized structure.

The localized and delocalized representations are observably distinct. A `"C#c-1"` pattern matches
the localized form and not the delocalized form; `"C#a2"` and `"C#a1"` likewise distinguish their
stored contributions. The choice cannot be made silently.

Consequently, `ingest::smiles` returns the localized form for inputs such as `[cH-]1cccc1` unless the
caller explicitly applies `DelocalizeCharge`. Resolution may still add aromatic and stereo
relations. With the explicit projection-reset options disabled, it changes no already determined
field or constraint.

## Aromaticity conformance

`AromaticityConformanceValidator` validates every stored representation of aromaticity against
independent perception:

- stored aromatic-system participant sets match the perceived sets;
- stored per-atom electron contributions match the perceived contributions;
- every non-vacuous localized-bond aromatic assertion agrees with the bonds induced by the stored
  systems;
- every non-vacuous aromatic-valence assertion agrees with the corresponding system contribution;
- an asserted aromatic atom or bond without a matching system, an extra stored system, or an explicit
  negative assertion on a system participant is a contradiction;
- absent and `Undetermined` projections are ignored;
- a non-ground, non-vacuous contribution that cannot be decided yields `Underdetermined`, not a false
  match.

This work was previously listed in doc 166 and belongs here because the resolver and validator must
use the same policy-free projection comparison.

## Open

- Whether the contradiction should name the offending atoms and the reason (element out of scope
  against ring size), beyond the exact unmatched atom and bond IDs required initially.
- Whether kekulization is reachable inside a reaction at all; it is not currently a primitive delta.
- Which high-level Python transformation surface should expose the source-resolve → kekulize → target
  resolve/validate chain.

## Staged implementation plan

Every Rust subitem carries focused `#[rstest]` coverage and leaves its affected crate green; Python
subitems carry focused pytest coverage. Tests use exact `Solution`, contradiction, edit-plan, error,
or molecule equality rather than summary-only assertions.

### S0 — Policy-free foundations

#### S0a — Charge-delocalization transformer **Done**

**Module:** `umol-graph/src/ops/transform/delocalize_charge.rs`,
`umol-graph/src/ops/transform.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add public `DelocalizeCharge` implementing `Transformer`. Move the model-independent charge and
π-contribution calculation behind this transformer without changing the existing resolver path yet.
Use `Infallible` as its error type: systems with non-literal data are left unchanged. Derive every
selected system update from the immutable input before applying the complete transformation.

Focused tables cover cyclopentadienyl anion, tropylium cation, a heterogeneous aromatic system,
non-literal input, and multiple systems. They assert charge conservation, π-electron conservation,
field preservation, and idempotence.

#### S0b — Aromatic perception and derivation

**Module:** `umol-graph/src/ops/aromaticity.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add `AromaticityPerception::derive`, returning public `AromaticityDerivation`. The result contains
the systems independently accepted by the selected model together with sorted unmatched atom, bond,
and existing-system IDs. Compare non-vacuous atom aromatic-valence assertions, non-vacuous
localized-bond aromatic assertions, stored aromatic-system participant sets, and stored electron
contributions. Existing systems must be re-perceived; their presence does not count as acceptance.
Absent and `Undetermined` projections are ignored.

Focused tables cover MDL furan rejection, Daylight furan acceptance, missing and extra atom/bond
assertions, vacuous projections, contribution mismatch, an already conformant system, and an
existing system rejected by the selected model. They assert the complete `AromaticityDerivation`.

#### S0c — Stereo perception and derivation

**Module:** `umol-graph/src/ops/stereo.rs`, `umol-graph/src/ops.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add public model-carrying `StereoPerception` and policy-free result `StereoDerivation`. Move the
current atom ligand-frame, bond ligand-frame, and side-ligand calculations into
public `StereoPerception::derive_stereo_atom` and `derive_stereo_bond`; only the side-ligand helper
remains private. The molecule-wide `derive` calls the per-entity methods, leaving the resolver wired
to its existing methods until the new implementation has parity. The result contains realizable
materializations and exact unrealizable or mismatched non-vacuous `#T` and `#C` sites without
choosing `Keep`, `Strip`, or `Error`. It also compares an asserted site already covered by a stereo
relation instead of skipping it.

Ordered table tests cover successful tetrahedral and cis-trans derivation, disabled kinds, element
scope, ligand arity, aromatic exclusion, existing elements, and exact atom/bond contradictions.

**Stage exit:** charge delocalization is available as an explicit transformer, and aromaticity and
stereo can derive and compare their relations without an operational policy.
`cargo test -p umol-graph` passes.

### S1 — Resolver contract migration

#### S1a — Aromatic resolver and Python configuration

**Module:** `umol-graph/src/ops/aromaticity.rs`,
`umol-graph/src/ops/resolve/aromaticity.rs`, `umol-graph/src/ops/resolve.rs`,
`umol-graph/src/parse.rs`, `umol-graph/src/ingest.rs`, `umol-py/src/resolve.rs`,
`umol-py/src/lib.rs`, and the affected Rust and Python configuration tests

**Kind:** breaking public config migration (red → green)

**Dependencies:** `[dep: S0a, S0b]`

Add `AromaticityInconsistencyPolicy::{Keep, Error}` beside `AromaticityResolveConfig`, defaulting to
`Error`, and add it as the config's `inconsistency` field. In the same public-config migration, remove
`delocalize_charge`; `AromaticityPerception::add_systems` and `AromaticityResolver` must no longer
equalize charge implicitly.

`AromaticityResolver::plan` consumes the S0b derivation before constructing edits. `Error` returns the
exact contradiction without mutation. `Keep` retains unmatched projections and adds only independently
accepted systems that are not already present; it also retains an existing system rejected by the
selected model. Repeated resolution of a conformant molecule is an identity. Add a
composite-resolver case in `ops/resolve.rs` proving that an aromatic contradiction rolls back edits
already applied by valence resolution.

Migrate all Rust config literals and the Python wrapper in the same subitem. Python exposes
`AromaticityInconsistencyPolicy`, makes `inconsistency` keyword-only, removes
`delocalize_charge`, and updates constructor signatures, getters, conversions, equality, repr, import
coverage, and all affected workflow/molecule/reaction tests. The localized post-ingest charge
representation becomes the new exact expected value.

#### S1b — Stereo policy rename

**Module:** `umol-graph/src/ops/resolve/stereo.rs`, `umol-graph/src/ops/resolve.rs`,
`umol-graph/src/ingest.rs`, `umol-py/src/resolve.rs`, `umol-py/src/lib.rs`, and the affected Rust and
Python configuration tests

**Kind:** breaking public type rename (red → green)

**Dependencies:** `[dep: S0c]`

Rename the existing public `InconsistencyPolicy` to `StereoInconsistencyPolicy` in Rust and Python.
Retain `Keep`, `Strip`, and `Error`, retain `Error` as the default, and change no stereo-resolution
behavior. Migrate every Rust import/config literal and every Python import, annotation, conversion,
repr, and exact enum/config test in the same subitem. Rewire `StereoResolver` to the S0c derivation
through `StereoPerception` without changing its planned edits.

**Stage exit:** each resolver owns a policy whose variants are sound for that operation; aromatic
resolution changes no determined value, while stereo retains its explicit stripping option.
`cargo test -p umol-graph`, `cargo test -p umol-py`, `maturin develop`, and
`pytest -q umol-py/tests` pass with Python 3.13 active.

### S2 — Conformance completion

#### S2a — Aromatic projection conformance

**Module:** `umol-graph/src/ops/validate/aromaticity.rs`, `umol-graph/src/ops/validate.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S0b, S1a]`

Keep `AromaticityConfig` as the validator's algorithm configuration; do not add a validation
inconsistency policy. Replace the current count-and-atom-set-only comparison with the S0b projection
comparison. An asserted projection without a matching relation, a stored relation rejected by
perception, or an explicit projection/relation mismatch is `Solution::Contradictory`. Absent and
`Undetermined` projections are ignored; a non-ground, non-vacuous contribution that cannot be
decided is `Solution::Underdetermined`. Extend `AromaticityValidatorContradiction` with exact
deterministic payloads and preserve the existing setup-error boundary.

Focused tables cover participant-set mismatch, per-atom contribution mismatch, localized-bond flag
mismatch, aromatic-valence mismatch, vacuous projections, model rejection of a stored system,
non-ground non-vacuous contributions, and a fully conformant existing system.

#### S2b — Stereo projection conformance

**Module:** `umol-graph/src/ops/validate/stereo.rs`, `umol-graph/src/ops/validate.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S0c, S1b]`

Keep `StereoValidateConfig` limited to graph-symmetry algorithms and iteration limits; do not add an
inconsistency field. Before the existing relation and symmetry checks, use `StereoDerivation` to
from `StereoPerception::derive` to require every non-vacuous asserted `#T` and `#C` site to have a
realizable, matching stereo relation. Absent and `Undetermined` projections are ignored. Add exact
atom/bond contradiction variants for an unrealizable, absent, or mismatched relation.

Focused tables cover absent and unrealizable tetrahedral/cis-trans relations, relation/assertion
mismatch, conformant existing relations, and preservation of every existing graph-symmetry
validation result.

#### S2c — Composite conformance outcomes

**Module:** `umol-graph/src/ops/validate.rs`

**Kind:** additive tests (green)

**Dependencies:** `[dep: S2a, S2b]`

Extend the composite validator tables with exact aromatic and stereo projection contradictions and
an underdetermined aromatic contribution. Assert that the contradiction is wrapped in the correct
`ValidatorContradiction` variant and that validation does not mutate the input.

**Stage exit:** standalone and composite validators reject projection/relation inconsistencies
without resolver policy or mutation. `cargo test -p umol-graph` passes.

### S3 — Public ingestion acceptance

#### S3a — Rust ingestion propagation

**Module:** `umol-graph/src/ingest.rs`, `umol-graph/src/parse.rs`,
`umol-graph/tests/resolution/`

**Kind:** additive (green)

**Dependencies:** `[dep: S1a]`

Add molecule and reaction SMILES table cases proving that MDL furan, thiophene, and pyrrole return
the exact aromatic contradiction by default, while explicit `Keep` preserves the unmatched
projections without adding a system. Retain positive MDL pyridine/benzene and Daylight
furan/thiophene/pyrrole references. Update charge-sensitive resolution fixtures and snapshots to the
localized representation; do not normalize them through `DelocalizeCharge` in expected-value
construction.

#### S3b — Python ingestion propagation

**Module:** `umol-py/tests/test_molecule.py`, `umol-py/tests/test_reaction.py`,
`umol-py/tests/test_workflow.py`

**Kind:** additive (green)

**Dependencies:** `[dep: S1a, S3a]`

Verify that configured molecule and reaction SMILES ingestion maps default aromatic rejection to the
existing `ContradictionError`, that explicit `Keep` returns the exact preserved representation, and
that the default output retains localized charge because Python has no implicit delocalization
option. No change to `umol-py/src/error.rs` is planned because the existing resolver-contradiction
mapping already has the correct public category. This stage does not add transformers or validators
to the Python surface.

**Stage exit:** Rust and Python ingestion demonstrate default rejection, explicit retention, and the
new localized charge representation. The feature-gated resolution conformance suite and installed
Python suite pass.

## Critical path and deferral

The aromatic resolver path is `S0a + S0b → S1a → S3a → S3b`; aromatic conformance is
`S0b → S1a → S2a → S2c`. The stereo cleanup is `S0c → S1b → S2b → S2c`.

No stage in this work unit is deferrable. Rich model-rejection reasons, `LocalizeCharge`,
reaction-level kekulization, and Python transformation/validation methods remain separate proposed
work.

Final verification is:

1. `cargo fmt --all`
2. `cargo test -p umol-graph`
3. `cargo test -p umol-graph --features conformance --test resolution`
4. `cargo test --workspace`
5. `cargo clippy --workspace --all-targets -- -D warnings`
6. With `umol-py/.venv` active and `python` confirmed as Python 3.13,
   `maturin develop` and `pytest -q umol-py/tests`
7. `git diff --check`
