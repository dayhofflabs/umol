# Code reviews

## Purpose

This is a normative guide for structured review cycles over existing code. It defines what a
review examines, what qualifies as a finding, the evidence a finding requires, the review and
refutation stance, and the form of the recorded result. The `review-cycle` repository skill
defines the operational procedure for running a cycle.

Reviews are functional, not aesthetic. Rewriting for its own sake is not a goal; a finding must
show that the code contradicts a settled design, violates repository policy, or forfeits a
property the design promises.

## Priority order

Logical correctness and adherence to the layered architecture precede efficiency. Efficiency lost
to a sound design can usually be regained by improving the design; a design constrained by
premature optimization is much harder to recover. Efficiency findings therefore classify two
ways:

- a design that forfeits efficiency — for example, a boundary that forces re-derivation or
  copying that a better-shaped API would avoid — is a design finding at full priority;
- locally suboptimal code under a sound design is recorded and ranked below correctness and
  layering.

Premature optimization that already constrains a design is the highest-value finding of all.

## Normative basis

A finding must cite the source it violates. The sources rank as follows:

1. the living guides in this directory;
2. the repository skills;
3. the `AGENTS.md` architectural policies and public API contract gates;
4. settled discussion documents for the reviewed scope — `discussion/000-status.md` is
   authoritative for status, and superseded reasoning must not be cited.

Current code and tests are authoritative for implemented behavior. Deliberate deferral is not a
flaw: before flagging an absence, check the status index and the governing document. Open scope
is described as deferred or open, never as a gap or debt.

## Review areas

**Construction, integrity, and fallibility.** Audit check placement against the contract
sections of [data-types.md](data-types.md). Classify every fallible constructor,
`Result`-returning operation, panic, and assertion as one of:

1. tier-1 representation integrity at a settled boundary;
2. a property this operation is the first to require, validated here;
3. re-validation of a property established upstream;
4. construction-time enforcement of a property no operation requires;
5. fallibility inherited from a callee that need not be fallible.

Classes 1 and 2 are compliant; a compliant check carries an `# Errors` or `# Panics` entry naming
its property. Classes 3–5 are findings, and each finding names the operation that should own the
check.

**Nomenclature.** Every symbol — types, functions, enum variants, struct fields, modules, and
test names, private as much as public — is checked against [nomenclature.md](nomenclature.md).
The two failure classes are synonyms (two terms for one concept, including across crates) and
collisions (one term for two concepts).

**Module structure and visibility.** Modules are organized by logical grouping, not
implementation history: focus types and functions first, roughly 1000 lines or fewer, long test
suites in separate test submodules. Visibility prefers private and `pub`; cross-dependencies
created by a split are resolved first by regrouping and only then with `pub(crate)`;
`pub(super)` is for test submodules and rare otherwise; `pub(in ...)` is not used. A proposed
split includes its visibility map.

**Tests and generators.** Property tests are reviewed as specifications: a suite can encode
wrong semantics correctly, so each stated property is traced to the design decision that
justifies it. Generator instance distributions are examined statistically — variant coverage,
size distribution, degeneracy rates — by sampling, not by reading the strategy code alone.
Coverage of at least 80 % per module is a baseline, not a goal; where a testable property
exists, property tests are preferred over unit tests.

**Documentation.** Doc comments are verified against actual behavior. Contract sections
(`# Assumes`, `# Establishes`, `# Errors`, `# Panics`, `# Semantic properties`) must be present
where required and accurate, and the semantic-properties sections must agree with the property
suites in both directions: every stated law has a test, and every tested law is traceable to a
stated source — the operation's prose or contract sections, or the property's own name and
documentation — so that a failing property identifies its semantic reason.
Shortening verbose prose counts as much as writing missing documentation.

## Findings and evidence

A finding is a contradiction with settled design or a violation of repository policy. A
genuinely unsettled question is an open item, not a finding.

A convention adopted after the reviewed code was written is a migration, not a defect: the
crate-wide absence of the convention is recorded once, and its adoption belongs to the proposed
design rather than the findings list. An individually wrong statement, or a specific law stated
nowhere in any form, remains a finding. Missing specification that is historical rather than
structural is filled, not litigated.

Evidence requirements:

- every finding cites its normative source;
- a semantic contradiction includes a reproduction — a concrete input with divergent or
  incorrect results, or a violated stated property;
- a generator finding includes the sampled statistics.

Findings are actionable only.

## Review and refutation stance

Both roles argue both sides and self-check before pushing; they differ in emphasis, weighted
three to one.

The review agent weighs three parts against the original to one part for it. Every finding it
pushes records two arguments: the violation claim with its citation, and the strongest case that
the original is correct or intentional, which requires the deferral check. A finding without its
defense argument is invalid. An unresolved tie is pushed, flagged uncertain.

The refutation agent inverts the weighting. Before dismissing a finding it first states the
finding's strongest form — narrowing an overbroad claim to the part that actually violates a
standard — and then argues the original's case. An unresolved tie favors the original. Verdicts
are graded:

- **Confirmed** — the finding stands as pushed;
- **Reduced** — the narrowed form survives; the original claim did not;
- **Open question** — neither side settles it; recorded as an open item;
- **Refuted** — dismissed, with the citation that dismisses it.

Findings that survive keep both recorded arguments attached.

## Output

A review cycle produces a dated discussion document containing scope, objective, the surviving
findings with verdicts and both arguments, open items, and a proposed design. It contains no
staged implementation plan, and the review changes no code. The document records the reviewed
commit, and accepted findings enter the ordinary design and planning lifecycle.
