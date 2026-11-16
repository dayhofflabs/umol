# 248 — Explicit repair during import

Status: Proposed
Date: 2026-11-16
Relates: [247](247-aromatic-system-bonds-2026-11-16.md),
[249](249-release-0-11-preparation-2026-11-16.md),
[data type contracts](../docs/development/data-types.md),
[operation contracts](../docs/development/operations.md)

## Scope and release boundary

Preserve the discussion of explicit repair infrastructure for external molecular
inputs whose descriptions do not fit umol's semantic model. The motivating
examples are the azulene and quinoid inputs in
[issue 45](https://github.com/dayhofflabs/umol/issues/45).

The general direction is accepted for further design: repair is an explicitly
requested operation, separate from ordinary raising and resolution, and returns
a record of the fixes it applies. The detailed semantics, names and interfaces
remain provisional and require further discussion.

This work is deferred beyond the 0.11 development cycle. No implementation is
authorized or started, and this record contains no implementation plan.

## Why a separate operation

Ordinary raising preserves representable source information. Resolution refines
that information and creates compatible entities; it does not silently discard
an incompatible description to make the input acceptable. Importing with fixes
may deliberately change that description and therefore has a different
preservation contract.

Resolver::resolve and Resolver::resolve_with_report currently produce the same
molecular result. The optional report describes unresolved completions and
tie-break uses. Adding fixes to that optional reporting arrangement would let
the reportless operation suppress the evidence of a changed source description.

A separately invoked repair operation makes permission to change the description
explicit before information is lost. Its result would always include the fix
record. A caller can still discard a returned record; the interface cannot
guarantee inspection. The explicit invocation and the returned record serve
different purposes: permission and inspection, respectively.

The separate operation may reuse resolution internally. The public distinction
does not require a second implementation of ordinary resolution.

## Motivating corrections

### Azulene with an explicit single fusion bond

The input `c1ccc2cccc-2cc1` describes the fusion bond as single while the perimeter
is aromatic. A repair could relax the corresponding negative bond-aromaticity
assertion so that umol can represent the ten-atom aromatic system. The resulting
system includes all induced edges. The input's distinction between the single
fusion edge and the aromatic perimeter would not be preserved on ordinary
aromatic export.

This must be distinguished from a single bond between separate systems, as in
biphenyl. The exact eligibility rule for relaxing an assertion is not settled;
blanket removal of negative bond-aromaticity assertions is not an agreed repair.

### Quinoid structures labeled aromatic

For p-xylylene, `C=c1ccc(=C)cc1`, repair could replace the source's aromatic
description with a localized Kekulé structure such as `C=C1C=CC(=C)C=C1`, with no
native aromatic system. This requires assigning ring bond orders consistently
with the explicit exocyclic double bonds and atom valences. Merely clearing
aromatic labels would not provide the required localized structure.

The two examples belong to the same explicit import-repair infrastructure but
need different corrections. One permits a native aromatic system while giving
up an edge distinction; the other replaces the aromatic description with
localized bonds. Neither correction promises to roundtrip the original aromatic
description. The applicability of quinoid repair beyond this example remains
open.

## Provisional interface

The proposed consuming interface accepts a faithfully raised Molecule. The names
and signatures below are a discussion sketch, not a completed SIN decision.

```rust
impl<'a> Repairer<'a> {
    pub fn new(
        model: &'a ChemistryModel,
        config: RepairConfig,
    ) -> Self;

    pub fn repair(
        &self,
        molecule: Molecule,
    ) -> Result<
        Solution<
            (Molecule, RepairReport),
            RepairContradiction,
            (),
        >,
        RepairError,
    >;
}
```

RepairConfig would specify which kinds of fixes are permitted and contain the
ordinary resolution configuration. The caller could permit the azulene
correction, quinoid kekulization, or both. Enabling a fix permits its application
where applicable; it does not require changing every input.

Determined would return the resolved molecule and its fix record together. An
already acceptable input would return the ordinary resolution result with an
empty fix record. Underdetermination, contradiction and execution failure would
remain distinct, following Resolver's consuming contract; those outcomes would
not publish a repaired molecule.

RepairReport would record the concrete semantic fixes actually applied. Ordinary
completion and justified assertion discharge remain resolution work, rather
than repairs. The proposal has one successful return form containing both the
molecule and report, with no reportless repair counterpart.

The individual fix record has not been designed. It still needs a precise
account of the affected sites, information discarded or replaced, and how those
references relate to the input and resulting molecule. No report fields,
fix-kind enum, configuration layout or defaults have been settled.

## Preliminary evidence and its limits

Scratch trials on 2026-11-16 showed that current resolution can find azulene's
ten-atom aromatic system after the generated negative bond-aromaticity assertion
is omitted. The trial used the SMILES valence model with the default chemistry
model's other components. The analogous guaiazulene trial also produced a
ten-atom system; fluoranthene produced systems of sizes six and ten, and biphenyl
two six-atom systems. These trials isolate the effect of the assertion; they do
not establish a general rule for selecting eligible bonds.

RDKit 2026.9.2 reperceives localized p-xylylene as aromatic and emits
`C=c1ccc(=C)cc1`. Its donor-count rule treats an exocyclic carbon differently
from an exocyclic oxygen. In umol, the two exocyclic-double-bearing ring carbons
have zero internal pi-electron contribution in the localized description, while
the other four contribute one each. The resulting four-electron ring is not
accepted by the Hückel rule. Reading the localized input succeeds without an
aromatic system. Thus accepting the foreign aromatic description does not
require adopting RDKit's aromaticity classification.

These observations motivate the proposed corrections, not a general repair
algorithm. They are recorded here so the design does not depend on preserving
temporary trial files.

## Remaining design questions

- Define each permitted fix's applicability and preservation guarantees,
  including how internal chord bonds are distinguished from bonds between
  systems and which quinoid descriptions can be localized.
- Verify that the faithfully raised Molecule retains the information needed to
  identify and report the intended corrections. Settle operation ownership and
  integration with external-format entry points.
- Define what happens when several corrections or localized assignments are
  possible, including selection, interaction between fixes and search limits.
- Settle the configuration layout, defaults, individual fix records, reference
  namespaces, and diagnostics for underdetermined or unsuccessful attempts.
- Complete SIN review of the public and internal interfaces, including Rust and
  Python exposure, before implementation planning.

The aromatic-system definition and HMO treatment recorded in discussion 247
remain unchanged. No bond exclusions are proposed. Fractional double-bond
character computed over Kekulé structures remains an independent calculation
and does not decide this import-repair design.
