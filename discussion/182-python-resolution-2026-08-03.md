# 182 — Expose resolution on the Python surface

Status: Proposed
Date: 2026-08-03
Relates: [178](178-python-lattice-ops-2026-08-01.md),
[179](179-python-editing-and-transactions-2026-08-02.md)

Third instance of the same gap. Doc 178 exposed the lattice operations, doc 179 the edit vocabulary,
and each time the argument was that an operation central to the model was unreachable from the
interface most users have. Resolution is the remaining one, and it is the largest.

A molecule built by editing cannot be resolved from Python. `MoleculeAst.from_smiles` accepts
`chemistry_model` and `resolve_config` because ingest invokes resolution; `MoleculeAst.parse` accepts
only `defaults` and does not resolve; there is no `resolve` method. The three `*ResolveConfig` classes
are exported and nothing outside reaction application consumes them.

## Justification

**Resolution is never automatic** (author, 2026-08-03), and that is the design. It follows that the
explicit operation has to be callable, otherwise the only way to reach it is to route a structure
back through SMILES.

**Two sections of the whitepaper are about it.** \Cref{sec:validity} is resolution and the chemistry
model; \Cref{sec:lattice} is the order that narrowing moves along. A reader can execute the lattice
operations after 178 and the edits after 179, and cannot execute the operation those two exist to
support.

**Section 9 needs it to state its own division of labour.** Mutation performs what was written;
resolution fills what was left open, when asked. Building methylamine has two routes — state the
hydrogen count explicitly, or leave `#h*` open and resolve — and only the first can be shown in a
listing today.

## Scope

**In:**

- Resolution on `MoleculeAst`, taking `chemistry_model` and `resolve_config`, both already exported.
- The verdict. See the open question below.
- Whatever is needed for the same operation on a reaction, if that is not already reachable through
  reaction application.

**Out:**

- The phase resolvers (`ValenceResolver`, `AromaticityResolver`, `StereoResolver`, `BondsResolver`,
  `MulticenterBondsResolver`). The composite operation is the interface; the phases are internal, and
  their ordering is a live design question (doc 174).
- `ResolverError` as a value. Rollback failure indicates a defect, not an outcome a caller plans for;
  an exception is right.

## Open question: verdict or exception

The Rust signature returns `Solution<(), ResolverContradiction>` — `Determined`, `Underdetermined`,
`Contradictory` — and Python currently has `ContradictionError` and `UnderdeterminedError` but no
`Solution`.

Two precedents point opposite ways.

- **178 argued for values.** `meet` and `join` return `None` when no bound exists, because "failure to
  relate two otherwise valid values is an ordinary result, not a Python exception." By that reasoning
  `Underdetermined` is an ordinary outcome of resolution and belongs in a returned verdict.
- **`from_smiles` raises today**, and consistency with it argues for raising.

The distinction that may resolve it: for ingest, a structure that will not resolve is a failure of the
call, so raising fits. For a molecule under construction, `Underdetermined` is frequently the expected
state and the caller wants to inspect it. If so, the answer is that ingest keeps raising and the
explicit operation returns a verdict, which is a difference in kind rather than an inconsistency.

Recommend deciding this before implementing; it determines the shape of every listing in
\Cref{sec:validity}.

## Settled semantics

- Rust mutates the AST in place. Python should follow 178's precedent and **return the resolved
  molecule without mutating the receiver**, consistent with `canonicalize`.
- `chemistry_model` and `resolve_config` are optional and default as they do for `from_smiles`.
- Resolution is implemented as an edit transaction with a rollback journal, so a contradiction leaves
  the input unchanged. That property should hold at the Python boundary and be tested.

## Verification

Follow 178 and 179: algebraic properties stay in Rust, and the Python tests check availability and
representative cross-boundary results. At minimum, the three verdicts each reached from a constructed
molecule, and the rollback property — a contradictory resolution leaves the input untouched.

The reader-facing check is that \Cref{sec:validity} and \Cref{sec:mutation} listings execute against
the built module before those sections ship.

## Note

Naming is unaffected by doc [176](176-ast-naming-2026-07-31.md); only the class this hangs off would
move.
