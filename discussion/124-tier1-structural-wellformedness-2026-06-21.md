# 124 — Tier-1 structural well-formedness validation

Status: Active
Date: 2026-06-21

## Motivation

`test_substructure_cross_validation_planted` fails on a persisted seed: the
`Incidence` strategy returns a duplicate occurrence (`[[0,3],[0,3]]`) where
`GraphAndOverlays` returns one (`[[0,3]]`). The host in that seed carries **two
parallel dative bonds** (acceptor `3` ← donor `{0}`, twice). In the Levi
(incidence) graph the two interchangeable dative pseudonodes admit two distinct
embeddings that collapse to the same atom correspondence; the strategy does not
deduplicate by atom correspondence, so it double-counts. `GraphAndOverlays` runs
the subiso on the plain bond graph and verifies overlays once per atom map, so it
does not.

Deduplicating the `Incidence` results would patch the symptom. The base cause is
that the molecule is structurally invalid: two identical dative relations. The
generator produces such molecules because nothing checks structural
well-formedness. This is a tier-1 (data-integrity) check that does not exist.

The fix is to add the missing structural checks to the tier-1 validator and make
the cross-validation property filter its inputs through it. The match path then
assumes well-formed input as a precondition; no `Incidence` dedup is added.

## Home

`EntityStructureValidator` in `umol-ast/src/ast/validate.rs` — already the tier-1
integrity validator. It currently checks only that `electrons.len()` matches the
participant count for aromatic systems and multicenter bonds. The new structural
rules extend it with new `EntityStructureContradiction` variants.

Construction stays infallible; the validator is a standalone pass. (Note: the
module doc comment claims it runs "at AST construction/raise", but no construction
or raise path currently calls it. Reconciling that — wiring it into raise vs.
correcting the comment — is a separate decision, deferred below.)

## Rules

Cross-*type* parallelism is always allowed (borazine: a localized bond and a
dative bond on the same atom pair coexist). The no-parallel rules are per type.

### Within a single relation

| Entity | Rule |
|---|---|
| Localized bond | endpoints distinct (no self-loop) |
| Dative bond | donor list has no duplicates; acceptor ∉ donors |
| Noncovalent bond | the two endpoints distinct |
| Aromatic system | participants distinct |
| Multicenter bond | participants distinct |

### Between relations of the same type

| Entity | Rule | Allowed examples |
|---|---|---|
| Localized bond | ≤1 per unordered atom pair | — |
| Dative bond | for a shared acceptor, donor sets must be vertex-disjoint; violation iff two datives share an acceptor **and** their donor sets share a vertex | Cp→Fe←Cp (same acceptor, disjoint donors); inverse sandwich (same donors, distinct acceptors) |
| Noncovalent bond | no two of the same `kind` on the same unordered pair | vdW + H-bond on one pair |
| Aromatic system | pairwise vertex-disjoint (no shared atom) | — |
| Multicenter bond | no two with identical participant set; partial overlap allowed | B₂H₆ B–H–B bridges overlap on both B |
| Stereo atom | site atoms pairwise distinct | — |
| Stereo bond | site bonds pairwise distinct | — |

The dative rule is deliberately strict as a starting point (it forbids two
distinct dative bonds from overlapping donor groups to the same acceptor); revisit
if a real case needs it relaxed.

## Proposed contradiction variants

Following the existing `<Entity><Problem>` naming (`AromaticSystemElectronsLengthMismatch`).
Names tentative — to confirm.

- `BondSelfLoop { atom }`
- `BondsParallel { atoms: [AtomId; 2] }`
- `DativeBondDonorDuplicate { acceptor, donor }`
- `DativeBondAcceptorIsDonor { atom }`
- `DativeBondsParallel { acceptor, shared_donor }`
- `NoncovalentBondSelfLoop { atom }`
- `NoncovalentBondsParallel { atoms: [AtomId; 2], kind }`
- `AromaticSystemDuplicateParticipant { atom }`
- `AromaticSystemsOverlap { atom }`
- `MulticenterBondDuplicateParticipant { atom }`
- `MulticenterBondsIdentical { atoms }`
- `StereoAtomSitesDuplicate { atom }`
- `StereoBondSitesDuplicate { bond }`

`validate` returns the first contradiction found, matching the current
short-circuit style. Within/between checks run per entity type.

## Test consumption

`molecule_strategy()` may emit structurally-invalid molecules. The
cross-validation invariant (all strategies + subiso algorithms agree) holds only
for well-formed hosts/patterns. Both `test_substructure_cross_validation` and
`test_substructure_cross_validation_planted` gain a `prop_filter` that drops
molecules for which `EntityStructureValidator::validate` returns `Contradictory`.
The persisted regression seed is then a filtered (skipped) input.

No change to `substructure_matches_incidence`: with well-formed input there is no
parallel-pseudonode multiplicity, so the two strategies already agree.

## Implementation

Done in `EntityStructureValidator` with the variants above. Notes:

- The bond check does not allocate. Each atom's CSR adjacency is sorted by atom
  id, so a self-loop is a neighbor equal to the atom and a parallel bond is two
  adjacent equal neighbors; a single `O(2E)` scan of `ast.neighbors(_)` suffices.
- The other checks use a small set/map per entity type (relations are sparse).
- Validator tests build inputs from the molecule DSL via `mol!`. This surfaced a
  roundtrip gap: the dative writer already emits `:donor` as a scalar (one donor)
  or a vector (many), but `read_dative_bond_entry` / `parse_dative_bond_entry`
  only read a scalar, so a multi-donor dative did not round-trip. Both readers now
  accept scalar-or-vector; `DativeBondEntryInput.donor` became `donors: Vec`. With
  that, every contradiction case is DSL-expressible (including donor-duplicate via
  `:donor [1 1]`).
- `test_substructure_cross_validation{,_planted}` now filter through the validator;
  both pass, and the planted regression seed is filtered.

## Deferred

- Whether raise (`umol-io`) and/or construction should run integrity validation
  and reject invalid input, and reconciling the `validate.rs` doc comment with
  actual call sites.
- Whether the generator itself should be constrained to emit only well-formed
  molecules (benefits all consumers, more invasive) instead of filtering per test.
- Stereo ligand-distinctness and other intra-relation stereo shape rules beyond
  site uniqueness.
