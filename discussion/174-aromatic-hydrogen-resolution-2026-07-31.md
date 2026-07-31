# 174 — Implicit hydrogen on bare aromatic heteroatoms

Status: Proposed
Date: 2026-07-31
Relates: [053](053-molecule-validation-scheme-2026-02-17.md),
[054](054-molecule-validation-phases-2026-02-17.md),
[171](171-aromaticity-inconsistency-policy-2026-07-29.md)

A lowercase aromatic heteroatom written without brackets has an implicit-hydrogen count that cannot be
computed before aromaticity is perceived, and an aromaticity that cannot be perceived before the
hydrogen count is known. This document records what the OpenSMILES specification actually says, what
the circularity costs today, and why the fix is a change of intermediate type rather than an iterative
resolver.

## Motivating input

`CHEBI:16848`, coenzyme γ-F420-2, SMILES copied verbatim from the ChEBI entry:

```
C[C@H](OP(=O)(O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2nc(=O)nc(=O)c-2cc2ccc(O)cc21)C(=O)N[C@@H](CCC(=O)N[C@@H](CCC(=O)O)C(=O)O)C(=O)O
```

Both toolkits reject it:

```
umol   ContradictionError: no atom-typing match for AtomId(18) (element C, charge Some(0))
RDKit  Can't kekulize mol.  Unkekulized atoms: 16 17 20 23 24 25 26 27 28 30 31
```

The pyrimidine-2,4-dione ring is written fully aromatic — `c2nc(=O)nc(=O)c-2` — with both ring
nitrogens bare. C2 and C4 spend their double bonds exocyclically and the C4a–C10a bond is pinned
single by the `-`, so no alternating assignment remains. The ring has no Kekulé structure as written.

**The entry is internally inconsistent.** Adding the lactam hydrogen at N3 — `nc(=O)[nH]c(=O)` —
kekulizes and yields C29H36N5O18P, average 773.60, monoisotopic 773.1793. ChEBI's own entry states
**C29H36N5O18P, average 773.598, monoisotopic 773.17930**. The published SMILES therefore cannot
correspond to the formula printed beside it. \umol\ reads the Kekulé form of the corrected structure
and agrees independently: 53 heavy atoms, 55 bonds, 6 stereocentres, C29H36N5O18P.

The H position is **unique** for this molecule: H on N3 kekulizes, H on N1 does not.

## What OpenSMILES actually says

Checked against `materials/formats/opensmiles/OpenSMILES specification.pdf` and the upstream
asciidoc source; they agree.

§3.1.5, Organic Subset, states one rule for the whole subset with **no aromatic exception**:

> The implicit hydrogen count is determined by summing the bond orders of the bonds connected to the
> atom. If that sum is equal to a known valence for the element or is greater than any known valence
> then the implicit hydrogen count is 0. Otherwise the implicit hydrogen count is the difference
> between that sum and the next highest known valence.

§3.6 lists the hydrogens that *must* be written explicitly — charged, H–H, bridging, deuterium and
tritium. Aromatic heteroatoms are not among them. The specification never mentions `[nH]` as a
requirement and never assigns zero hydrogens to a bare lowercase heteroatom.

**So "absent a bracket H, assume `#h0`" is not in the specification.** That reading is ours.

**But the specification's rule is inapplicable as written.** "Summing the bond orders" presupposes
integer bond orders, which a lowercase ring does not have until it is kekulized, and §3.1.5 never says
to kekulize first. This is a defect in OpenSMILES, not in our reading of it. Kekulize-then-count is the
only interpretation that makes the rule executable.

Observed behaviour on bare aromatic heteroatoms (2026-07-31):

| input | RDKit | umol |
| --- | --- | --- |
| `c1cccn1` pyrrole, bare n | rejected | accepted → C4H4N, **0 aromatic systems** |
| `c1cncn1` imidazole, bare n | rejected | accepted → C3H3N2, **0 aromatic systems** |
| bare-n indole | rejected | accepted, 0 aromatic systems |
| `O=c1ccccn1` 2-pyridone, bare n | rejected | rejected |
| `c1cc[nH]c1`, `c1c[nH]cn1`, `c1ccoc1` | accepted | accepted |

RDKit refuses rather than guessing. We accept and emit a structure that is not a molecule: `c1cccn1`
returns **C4H4N**, a neutral species with 25 valence electrons and no radical marker, carrying `#a`
on every atom with no aromatic system to project from. That is the dangling-projection state of doc
171 reached by a second route — ingest-time kekulization failure rather than model refusal.

**Accepting is worse than rejecting here.** Doc 171's policy should be widened to cover any path that
leaves aromatic projections without a materialized relation, not only a model declining a system.

## The circularity

Implicit hydrogen needs bond orders; bond orders need a Kekulé assignment; the assignment needs the
per-atom π contributions, which depend on the hydrogen count. Pyridine-type nitrogen contributes one
π electron and carries no hydrogen; pyrrole-type contributes two and carries one.

Underdetermination is genuine. `c1cccn1` with the hydrogen free admits at least three completions,
all C4H5N and all real compounds: pyrrole (N–H), 2H-pyrrole, 3H-pyrrole. The aromatic form does not
determine which. This is why the `[nH]` convention exists in practice, and why pinning `#h0` was a
defensible engineering choice even though the specification does not license it.

Uniqueness is a property of the instance, not of the problem: the F420 ring above has exactly one
completion. An operation here should therefore return a **solution set with its cardinality**, not a
solution.

## What the resolver does today

The registry already carries both candidate states:

```toml
"N #n #v2 #a",   # C5H5N        pyridine-type: one lone pair, one pi electron, no H
"N #v2 #a2 #h",  # C4H4NH       pyrrole-type: no lone pair, two pi electrons, one H
```

Resolution was run directly on aromatic-flagged rings with charge pinned to zero and `#h` left
undetermined (`MoleculeDefaults::new()`, default `ChemistryModel`):

| input | verdict | result |
| --- | --- | --- |
| pyrrole, `#h` open | `Determined` | **correct** — `[1,1,1,1,2]#c0` materialized, N is `#h1 #a2` |
| pyrrole, `#h0` pinned | `Determined` | no aromatic system; N `#h0 #n1 #a1`; `#a` flags dangling |
| pyridine, `#h` open | `Determined` | **wrong** — N is `#h1 #a2` in a six-ring; 7 π electrons; no system |
| imidazole, `#h` open | `Determined` | **wrong** — *both* nitrogens `#a2`; no system |

Leaving the hydrogen open fixes the case that was broken and breaks the two that worked. All four
report `Determined`, including the three that are chemically wrong.

**The cause is where the candidate set collapses, not the direction of the pipeline.** The valence
phase narrows to a set and then selects a single winner with `compare_valence_preference`, a local
ordering that cannot know which choice lets the ring satisfy 4n+2. Aromaticity never votes.

## Proposal

Doc 053 already specified the shape: `Phase: Set<CandidateState> → Set<CandidateState>`, with the
empty set as failure and no backtracking because commitment is deferred. The implementation collapses
one phase too early.

- The valence phase should **emit the candidate set** when narrowing leaves more than one member,
  rather than choosing by preference.
- The aromaticity phase selects from it, using the criterion it alone has: whether the resulting π
  system satisfies the model.
- An empty set at any point is a contradiction. A set with more than one member surviving the final
  phase is `Underdetermined`, not a silent pick.
- `compare_valence_preference` remains a tie-break of last resort, applied only where the ambiguity
  survives every phase, and its use should be visible in the verdict.

The pipeline stays one-directional — valence → aromaticity → stereo. Only the type of the value
passed between phases changes. No fixpoint iteration is required.

Consequence for ingest: the SMILES reader should stop committing `#h0` on bare aromatic heteroatoms
and leave the field undetermined, letting the phases decide. That is a prerequisite, not a separate
change — the resolver cannot narrow a field the reader has already pinned.

## Regression triple

Any implementation that gets all three right has deferred the collapse rather than retuned the
preference order:

| input | expected |
| --- | --- |
| `c1cccn1` | N `#h1 #a2`, one aromatic system `[1,1,1,1,2]` |
| `c1ccncc1` | N `#h0 #n1 #a1`, one aromatic system `[1,1,1,1,1,1]` |
| `c1cncn1` | one N `#a1`, one N `#h1 #a2`, one aromatic system |

Add the corrected F420 structure as a fourth case once completion is implemented: unique solution,
C29H36N5O18P.

## Open

- Whether hydrogen completion is offered at all, or whether unkekulizable input is simply rejected
  with a diagnostic naming the ring. Rejection is defensible and is what RDKit does; completion is
  more useful for database input, which is where such structures come from.
- If completion is offered, whether it belongs in ingest configuration or as an explicit
  transformation the caller invokes, by the doc 166 boundary. Completion changes localized bond
  orders, which argues for a transformation.
- Whether the same premature collapse affects other phases — stereo in particular, where a similar
  local preference may exist.
- Separately, and not part of this circularity: umol rejects an aromatic carbon bearing an exocyclic
  double bond even in correctly written input. 2-pyridone, uracil and 4-pyranone all fail as `[nH]`
  forms that RDKit accepts. The registry has no row for a carbon at localized valence 4 with aromatic
  valence 0. This is coverage, and it also blocks the F420 structure independently of the hydrogen
  question.
