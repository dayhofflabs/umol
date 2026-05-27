# OpenSMILES Spec Validation Strategy

Date: 2026-03-06

## Context

The UMOL OpenSMILES formalization (`spec/opensmiles-spec.md`) and its diagnostics registry
(`spec/opensmiles-errors.md`) define the normative behavior of the SMILES parser. This document
records strategies for validating spec completeness and generating conformance test inputs.

## Validation approaches

### 1. Error code coverage (low effort, high signal)

For every code in `opensmiles-errors.md`, verify that at least one test case triggers it. This is
mechanical: enumerate the error codes and check for corresponding test expectations.

Status: not yet implemented. The `ParseError` enum in `error.rs` has variants that are not yet
registered in the errors doc (see below). Conversely, some error codes in the doc correspond to
semantic validation that is not yet implemented (Part 3 stubs).

### 2. Graph-first test generation (medium effort, high coverage)

Generate valid SMILES by construction rather than by filtering. This sidesteps the
context-sensitivity of SMILES ring closures (which are mildly context-sensitive due to crossing
ring indices; see notes below).

#### Algorithm

1. **Stoichiometry.** Pick numbers of atoms of each element from the organic subset (excluding H).
2. **Valence and degree.** Pick a valid valence for each atom. Pick degree (number of bonds) ≤
   valence. The degree sequence must be graphical (Erdős–Gallai condition). Valence and degree are
   coupled — elements with valence 1 (F, Cl, Br, I) force degree 1.
3. **Connected graph construction.** Build a random spanning tree to guarantee connectedness:
   start with a connected set containing atom 0; repeatedly pick a random unconnected atom and
   attach it to a random atom in the connected set. Then add extra edges to reach target degrees:
   pick pairs where both atoms have residual degree > 0, no edge exists, and add the edge. Repeat
   until all residuals are zero or no valid pair remains. If the residual sequence is not
   realizable on the remaining non-edges, retry from step 2.
   A simpler variant that avoids retries: pick a target edge count M, build the spanning tree
   (N−1 edges), add M−(N−1) random edges (avoiding duplicates/self-loops), and derive degrees
   from the result rather than prescribing them.
4. **Atom assignment.** Assign atoms within each degree class to nodes with corresponding degrees.
   Compute valence defect = valence − degree for each atom.
5. **Triple bonds.** Iterate over edges. If defect ≥ 2 for both endpoints, replace the single bond
   with a triple bond with some probability. Adjust defects by −2.
6. **Double bonds.** Iterate over edges. If defect ≥ 1 for both endpoints, replace with a double
   bond with some probability. Adjust defects by −1. Steps 5–6 are greedy and may fail to find a
   valid assignment even when one exists; this is acceptable for test generation.
7. **Charge.** For parser testing: assign random charge values in a valid range, independent of
   valence bookkeeping. For chemically plausible molecules: the charge–valence relationship is
   element-dependent (e.g., CH3+ has valence 3 with defect 1; NH4+ has valence 4 with defect −1),
   and charge should be decided before or jointly with valence.
8. **Unpaired electrons.** Same as charge: random for parser testing, element-dependent for
   chemical plausibility.
9. **Implicit hydrogens.** Assign implicit H = remaining defect.
10. **Serialize to SMILES.** DFS-linearize the graph. Back edges become ring closures.
11. **Decorations (optional).** Layer on bracket fields (isotope, chirality, atom class),
    aromatic lowercasing, stereochemistry markers as independent random passes.

#### Notes on SMILES context-sensitivity

SMILES ring closures are mildly context-sensitive: ring indices can cross (e.g., `C1C2CC1CC2`),
which is a cross-serial dependency that no CFG can express. Asymmetric delimiters (`{1`/`}1`)
would not help because crossing matched brackets still require more than a pushdown automaton.
Since the ring index space is finite (0–99), the tracking state is finite, so SMILES is
CFG ∩ Regular in theory — but the grammar would need ~2^100 nonterminals, which is useless in
practice. Graph-first generation avoids the issue entirely: ring closures arise naturally from
DFS back edges.

#### Notes on sampling distributions

Uniform sampling over all graphs with a given degree sequence is studied (configuration model,
Boltzmann sampling). Uniform sampling over chemically valid molecular graphs with a given
stoichiometry is much harder due to the irregular constraint surface. For parser conformance
testing, uniformity does not matter — structural diversity (varying ring counts, branching
patterns, bond orders, charges) is more useful than statistical representativeness.

#### Integration

For `proptest`, implement a custom `Strategy` that builds a `petgraph` graph per the algorithm
above and serializes it. The current fuzzing suite (`umol-fuzz`) only verifies absence of panics.
Graph-first generation would enable property-based testing of round-trip invariants
(parse → serialize → parse = same graph) and semantic properties (atom count, bond count).

### 3. Normative clause tagging (high effort, highest assurance)

Extract every MUST/SHALL/MUST NOT statement from the spec, assign each a stable clause ID, and tag
test cases with the clauses they exercise. Any clause without both a positive and a negative test
is a coverage gap.

This is labor-intensive to set up but provides the strongest completeness guarantee. Could be
automated with a clause extractor that parses the markdown for RFC 2119 keywords.

### 4. Mutation testing (medium effort, complementary)

Take valid inputs, apply small mutations (delete a character, swap tokens, duplicate a bracket
field), and verify the parser rejects each. Tests the negative constraints of the spec. Can be
combined with graph-first generation: generate valid → mutate → assert rejection.

## Current gaps

### ParseError variants not in errors doc

The following `ParseError` variants in `error.rs` are not yet registered in `opensmiles-errors.md`:

- `MismatchedRingBondDonations` — dative bond ring closure conflict (CXSMILES extension only)
- `InvalidCxTag` — CXSMILES tag parsing error (extended parser only)
- `MissingReactionArrow` — reaction SMILES (extended parser only)
- `AtomIndexOutOfBounds` — CX extension atom index reference error
- `BondIndexOutOfBounds` — CX extension bond index reference error
- `MismatchedAtomBondIndices` — CX extension consistency error
- `SgroupIndexOutOfBounds` — CX extension S-group index reference error

These are all extended-parser-specific. The basic parser's error surface is fully covered.

### Part 3 semantic stubs

The following sections in `opensmiles-spec.md` are marked as "to be specified":

- Valence and Hydrogen Validation
- Stereochemistry Validation
- Aromaticity Validation

The corresponding error codes exist in `opensmiles-errors.md` (VAL\_\*, AROM\_\*, STEREO\_\*) as
forward-looking definitions.

## Recommendation

Start with error code coverage (approach 1) as a quick win. Follow with graph-first generation
(approach 2) for systematic positive-case coverage. Clause tagging (approach 3) can be deferred
until Part 3 semantics are specified.
