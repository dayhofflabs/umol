# Property tests as executable specifications

## Scope

Property tests have two distinct roles in umol:

1. they detect defects by exercising implementations over many inputs;
2. they state semantic properties of the library in an executable form.

The second role should be visible to library users. A property hidden only in
test source still contributes validation evidence, but it does not communicate
the public contract at the point where users need it.

This document defines the terminology and documentation policy for exposing
those properties without introducing a second hand-maintained inventory of the
property suite. It uses `umol-graph-core` as the clearer case, where the
properties are predominantly topological, and `umol-ast` as the more realistic
case, where one public operation is commonly exercised over several distinct
operational domains.

This guide does not claim that every existing property is already documented publicly.

## Specification and evidence

The semantic assertion and the evidence supporting it are different things:

| Concept | Meaning |
| --- | --- |
| Semantic property | A quantified assertion about public behavior. It is part of the API contract and is independent of any particular generator. |
| Operational domain | The class of inputs or states over which one execution of the property is meaningful: raw values, canonical values, valid DSL values, malformed reactions, simple graphs, multigraphs, and so on. |
| Validation method | The means used to check the assertion: direct expected result, definition-level reference implementation, comparison between implementations, transformation relation, bounded exhaustive collection, or literature result. |
| Evidence scope | The actual generated distribution, size bound, checked-in collection, or set of named examples over which the validation ran. |
| Property test | Executable code combining an operational domain, a validation method, and an assertion. |
| Regression case | A retained input known to have violated a property. It preserves a discovered example but does not define the property. |

The semantic property is normative. The operational domain and evidence scope
describe how the implementation has been tested. A generator must not silently
become the definition of the public behavior.

For example:

```text
Semantic property:
    update(x, difference_to(x, y)) = y

Operational domains:
    canonical entity ASTs
    entity ASTs containing independent undetermined fields
    each of the eight entity families

Validation method:
    construct the update, apply it, and compare with y

Evidence scope:
    the strategies and case count used by the current test run
```

Changing the generated distribution does not change the semantic property.
Discovering that the assertion requires a precondition does change the public
property and must be reflected in its documentation.

## Relation to verification

Property-based testing is specification-based testing: it makes general
statements executable and searches for counterexamples. It is not, in general,
a proof of those statements.

The appropriate claims are:

| Validation | Public claim |
| --- | --- |
| Generated inputs | “Generatively validated over …” |
| Every member of a finite domain | “Exhaustively checked for every … through …” |
| Simple reference implementation | “Checked against a definition-level implementation over …” |
| Independent implementation | “Cross-validated against … over …” |
| Related executions | “Validated under relabeling / reversal / roundtrip / composition …” |
| Published examples | “Checked against the results reported for …” |
| Arbitrary bytes or syntax | “Fuzz-tested for robustness” or “parser parity checked over arbitrary strings” |

“Verified” or “proved” should not be used without a precise qualification. A
bounded exhaustive test may establish the assertion for the stated finite
domain, assuming the completeness of the input collection and correctness of
the assertion machinery. It says nothing by itself about inputs outside the
bound. A comparison between two implementations establishes agreement, not
independent correctness; shared defects remain possible.

The useful analogy to verification is therefore methodological rather than
logical. Property tests require the library to state its invariants,
equivalences, preservation laws, and failure behavior precisely. Their
execution supplies repeatable evidence for those statements, but sampled
execution does not acquire the guarantee of a proof merely by increasing the
number of cases.

## Public documentation policy

### Semantic properties belong with the public API

Public types, traits, and operations with non-trivial algebraic or semantic
contracts should have a `# Semantic properties` rustdoc section. The section
should state:

- the quantified values and required preconditions;
- the equality or equivalence relation used by the assertion;
- the result or error behavior;
- any preservation, completeness, uniqueness, or roundtrip requirement.

Use a more specific term inside the section where helpful:

- **law** for associativity, commutativity, absorption, idempotence, or
  involution;
- **invariant** for a condition preserved by an operation;
- **roundtrip property** for paired representations;
- **soundness** and **completeness** for matching or enumeration;
- **transaction property** for application and rollback;
- **transformation property** for relabeling, reversal, permutation, or
  canonicalization.

“Metamorphic property” is the established testing term for the last category,
but the public documentation should prefer the concrete transformation when
one is available.

### Validation evidence remains separate from the contract

A stable description of the evidence supporting a property may follow the contract in crate-,
module-, trait-, or operation-family documentation. Do not require a separate `# Validation`
section on every method. Exact strategy weights and run counts belong in test configuration and CI
output rather than public rustdoc.

For example:

```rust
/// # Semantic properties
///
/// For every valid `x` and `y`,
/// `x.update(&x.difference_to(y)) == *y`.
///
/// An update derived from `x` to itself is empty.
```

The operation-family documentation may then state that the properties are checked generatively for
every entity AST family, including values with independently undetermined spin fields. The example
is not a requirement to repeat the same text on eight implementations. A shared trait or
operation-family location is preferable when it accurately states the common contract.

### The property target remains the authoritative inventory

The executable test target and `cargo test --test property -- --list` remain
the authoritative inventory of property tests. The existing
`umol-ast/tests/property.rs` policy is retained: a property-suite README would
duplicate the source and drift.

A crate-level “Correctness and validation” page may:

- define the evidence vocabulary above;
- identify the major semantic property families;
- link to the public APIs where the properties are stated;
- explain how to run the property suite.

It must not reproduce a test-by-test coverage table. Public rustdoc is
authoritative for the semantic contract; the test target is authoritative for
the executable validation inventory.

### Test source should point back to the property

Property modules should explain:

- which public properties they validate;
- the operational domains generated in that module;
- the role of any definition-level or independent reference;
- why apparently overlapping properties are distinct.

The test implementation need not duplicate the complete public statement.
Where the connection is not obvious from the operation under test, a short
comment should name the documented property. Test names continue to follow the
workspace test-writing conventions; they are not a substitute for the
specification.

## Property-suite organization

Organize a property suite along three independent axes:

| Concern | Location |
| --- | --- |
| Semantic property | Subject and operation modules |
| Operational domain | Domain-specific strategy modules |
| Validation method | Operation-local reference support, exhaustive support, captured results, or literature fixtures |

The property modules are the executable inventory. Their hierarchy should make the public subject
and operation visible without reproducing the production source tree mechanically. Uniform law
families such as lattice laws or entity-update laws may remain together even when they cover many
types; molecule, reaction, stereo, or algorithm families should split by operation once one flat
module obscures the asserted properties.

### Strategy modules

- Split strategies by the domain generated: values, entity ASTs, molecules, reactions, edits,
  metadata, and other independently meaningful domains.
- Import production symbols directly in the property module that tests them. A strategy module must
  not act as a wildcard prelude that re-exports the production API.
- Import only the strategies a property module uses, so its operational domain is visible at the
  call site.
- Keep an operation-specific scenario beside that operation's properties unless more than one
  property family genuinely shares it.
- Generate simple independent parameters and derive structurally valid aggregates deterministically
  where practical. Prefer this to broad rejection with `prop_filter`, which hides the admitted
  domain and weakens shrinking.
- Generate invalid inputs by applying one named defect to a valid base when the property concerns an
  exact failure boundary. Carry independently derived expected information in a scenario value when
  recomputing it through the production operation would make the assertion circular.
- Keep raw, canonical, ground, structurally valid, semantically valid, and deliberately malformed
  domains distinct. A strategy's name and documentation must say which domain it emits.

### Reference and fixture support

- Keep definition-level reference implementations independent of the optimized production path and
  separate from data generation.
- Locate reference support with the operation it validates. Do not accumulate unrelated references
  in a global helper module.
- Preserve the provenance of captured and literature results. Checked-in captured output must not
  make the external producer a runtime test dependency.
- Keep regression files at the nearest stable subject root so reorganizing operation modules does
  not orphan minimized failures.

### Structural constraints

- Follow the repository module convention of `<module>.rs` plus `<module>/<child>.rs`; do not use
  `<module>/mod.rs`.
- Use the narrowest test-crate visibility supported by the hierarchy without widening production
  visibility.
- Avoid `common`, `shared`, or `utils` dumping grounds. Name support by its domain or validation
  role.
- Module documentation states the properties, operational domains, and validation methods represented
  below it, including why deliberate overlaps are distinct.
- Do not create a property-suite README or a manually maintained property table. The test target and
  `cargo test --test property -- --list` remain authoritative.

## `umol-graph-core`: clear topological properties

The cycle operations illustrate a comparatively clean specification. The
objects and properties are topological, the size bounds are explicit, and
small graphs admit direct definition-level enumeration.

### Simple-cycle enumeration

A public semantic statement for simple-cycle enumeration can be organized as
follows:

```text
For graph G and maximum source-edge size b:

- every emitted Cycle is a cycle of G with at most b source edges;
- every such cycle is emitted exactly once;
- node relabeling changes identifiers but not the represented cycle set;
- enumeration and complete visitation produce the same cycle sequence;
- direct and fallback operations agree on simple graphs;
- the direct operation rejects a non-simple graph before emitting a cycle;
- the combined operation agrees with the fallback on non-simple graphs.
```

The current property suite supplies several kinds of evidence:

| Assertion | Operational domain | Validation |
| --- | --- | --- |
| Soundness, completeness, size bound, and uniqueness | Generated multigraphs with bounded nodes and edges | Compare with enumeration of every edge subset satisfying the definition of a cycle. |
| Direct/fallback/combined behavior | Generated simple and non-simple graphs | Compare all operation paths and their exact errors or emitted cycles. |
| Visitor/enumerator agreement | Generated multigraphs | Collect visitor emissions and compare with enumeration. |
| Relabeling behavior | Generated multigraphs under a fixed node reversal | Compare normalized source-edge sets. |
| Small simple-graph completeness | Every checked-in simple graph through order six | Compare with definition-level exhaustive enumeration. |
| Wider independent agreement | Captured simple-graph results through order eight and captured bounded multigraph results | Compare exact cycle collections. |
| Chemically and graph-theoretically significant cases | Named literature fixtures | Compare source-stated results and preserve provenance. |

The table communicates why the overlaps matter. The generated comparison finds
small counterexamples and exercises multigraph behavior; the bounded
simple-graph collection removes sampling within its bound; captured results
provide an independently produced comparison; literature fixtures retain
meaningful named structures. None is merely “more cases” for the same purpose.

### Minimum bases, relevant cycles, and ring families

The same organization applies to the higher cycle operations:

- a minimum cycle basis has the cycle-space dimension, consists of linearly
  independent cycles, and has minimum total length;
- a relevant cycle belongs to at least one minimum cycle basis;
- Unique Ring Families partition relevant cycles according to the documented
  relation;
- stored family nodes, edges, weights, and counts agree with the represented
  cycles;
- relabeling preserves the corresponding mathematical result.

For small generated graphs, the suite constructs all cycle bases directly,
selects those with minimum weight, derives relevant cycles from their union,
and derives ring families from the definition. The production Horton,
Vismara, and Kolodzik implementations are then compared with those results.
For the larger checked-in collection, captured independent results extend the
evidence beyond the range where direct basis enumeration is practical.

This is a strong example of property tests as communication: the exhaustive
test support describes the mathematical definitions more directly than the
optimized implementations do.

### Algorithm cross-validation is relative evidence

The subgraph-isomorphism suite compares the named implementations with VF2 on
generated labeled graphs and on planted subgraphs. This establishes agreement
over the generated domain and ensures that the planted cases have at least one
match. It should be described as cross-validation, not as an unconditional
proof that every implementation is correct.

Where a separate theorem justifies a reduction, as recorded for ArcMatch, that
argument and the cross-validation evidence have different roles:

- the theorem justifies the transformation;
- the property test checks the implementation and its agreement with the
  reference path.

## `umol-ast`: properties over multiple operational domains

`umol-ast` cannot be described adequately by saying that it generates
arbitrary ASTs and checks laws. The same surface participates in several
different semantics, and the operational domain is part of understanding the
evidence.

### Lattice values

The lattice suite distinguishes at least two domains:

| Domain | Purpose |
| --- | --- |
| Canonical, satisfiable values | Check the full lattice laws and laws that depend on canonical representation. |
| Raw, possibly non-canonical values | Check input-canonicality-independent behavior, canonical folding, `matches`/`meet` consistency, and canonicality of results. |

The public `Lattice` contract should state the algebraic laws once. The
validation description should then say that they are exercised both on
canonical values and on raw expression trees. Omitting the second domain would
leave the canonicalization path largely untested even if every generated
canonical value satisfied the lattice equations.

### Serialization and parsing

Serialization contains several related but non-identical properties:

| Property | Operational domain | Comparison |
| --- | --- | --- |
| `parse(display(x)) = x` | Generated entity DSL values whose surface is lossless | Exact representation equality |
| EDN tree roundtrip | Generated AST or DSL values | Exact representation equality |
| Streaming/tree parser parity | Valid rendered values and arbitrary strings | Equal parsed value or equal rejection |
| Stable rendering | Successfully parsed canonical surface | Equality of the first and second rendering |
| Defaults projection roundtrip | Values interpreted under a specified defaults configuration | Equality after applying the same configuration semantics |
| Vacuous constraint elision | Explicitly constructed vacuous constraints | Equality with the corresponding bare surface |

These cannot be collapsed into a single “serialization roundtrips” claim.
Each property uses a different domain and comparison relation. In particular,
surface normalization may require canonical equality rather than stored
representation equality, while a lossless AST EDN roundtrip requires `==`.

### Updates, deltas, and edits

The entity update properties state a direct public law:

```text
x.update(x.difference_to(y)) = y
```

The delta properties are related but distinct:

- deriving field and constraint deltas and applying them reaches the target;
- diffing an entity with itself emits no deltas;
- delta inversion is an involution;
- canonicalizing a consistent `Deltas` value is idempotent;
- edit application followed by rollback restores the original molecule;
- appending transaction batches preserves the sequential application and
  rollback semantics.

The direct update law is about entity ASTs. Delta properties are about a
reaction-side change representation. Edit and transaction properties are about
host-molecule mutation. The similar equations do not make the three
operational domains interchangeable.

### Molecule comparison

The comparison suite already records distinct relations:

- `==` compares stored representation;
- `equiv` compares molecule semantics in the current ID and participant frame;
- `equiv_under` compares after an explicit correspondence;
- entity `canonical_eq` compares canonical entity semantics.

The properties exercise:

- reflexivity and symmetry of `equiv`;
- agreement of `equiv` with `==` on canonical ASTs;
- reduction of `equiv_under` to `equiv` under the identity correspondence;
- symmetry under reversing a correspondence;
- composition of correspondence-aware equivalence on generated atom
  reorderings.

This example shows why the comparison relation must appear in every public
property. Writing only “the result equals the input” is ambiguous for an AST
with representation equality, canonical equality, same-frame equivalence, and
correspondence-aware equivalence.

It also exposes a useful documentation question: if `equiv` is publicly
presented as an equivalence relation, transitivity is part of that claim. The
current suite directly checks reflexivity and symmetry, while its transitivity
property is expressed through `equiv_under` on a restricted atom-reordering
domain. A later documentation and coverage pass should either state and test
general `equiv` transitivity or document a narrower contract.

### Reaction operations

Reaction properties span several operational domains:

| Domain | Representative properties |
| --- | --- |
| Generated well-formed reactions | Canonicalization idempotence, reaction/span reconstruction, derivation reversal, composition, and serialization |
| Comprehensive entity reactions | Roundtrips and transformations across all eight entity families |
| Host-relative refinements | Pattern-relative updates lower against the matched host value rather than replacing it with the pattern value |
| Explicit correspondences | `apply_at` agrees with a matching-derived application for the same match |
| Malformed reactions | Invalid references, incidence mismatches, discontinuous updates, and invalid stereo configurations return exact typed errors without panics |
| Transaction failures | A fatal application error is emitted once and terminates the iterator |

The error properties are part of the executable specification just as much as
successful roundtrips are. Generating only valid reactions would validate the
happy path while leaving the public failure contract unstated.

Algorithm choices in these properties must remain explicit. A reaction
application property using `GraphAndOverlays` and VF2 validates the operation
under those selected algorithms; it must not imply that an unmentioned hidden
default is part of the semantic contract.

## Public documentation shape

The preferred public shape is:

```rust
/// Performs the operation.
///
/// # Semantic properties
///
/// - State the property and its preconditions.
/// - Name the equality or equivalence relation.
/// - State exact failure behavior where it is part of the contract.
```

Stable evidence categories and operational domains belong in the applicable crate-, module-,
trait-, or operation-family documentation, with links to definitions or published sources where
they govern the assertion. Not every method needs a semantic-properties section. Trivial accessors
need none, and a
family of operations should document a common law once when the shared
location is discoverable. The documentation is warranted when the property
clarifies semantics that are not apparent from the type signature.

## Change policy

- A new semantic operation should identify its principal properties before its
  property tests are treated as complete.
- A changed semantic property requires corresponding changes to public
  documentation and executable validation.
- A changed generator or case count changes evidence, not the contract, unless
  it reveals that the documented operational domain was inaccurate.
- A minimized regression remains attached to the property that it violated.
- Deliberate overlap should be justified by different operational domains,
  validation methods, or evidence scopes.
- A passing release test suite means that the documented properties were
  validated by that release's checked-in evidence. It does not turn sampled
  properties into proofs.

## Precedents

- QuickCheck describes the programmer-provided properties as a specification
  of the program and tests them over generated cases:
  <https://hackage.haskell.org/package/QuickCheck>.
- John Hughes' *How to Specify It!* treats property-based tests as tests
  against a specification rather than a set of examples:
  <https://research.chalmers.se/en/publication/517894>.
- ScalaCheck calls collections of properties specifications directly:
  <https://scalacheck.org/>.
- Rust's `Ord` documentation is a conventional API example of publishing
  algebraic consistency and transitivity requirements independently of their
  tests: <https://doc.rust-lang.org/std/cmp/trait.Ord.html>.
- Chen, Cheung, and Yiu introduced metamorphic testing for predictable
  relations between related executions:
  <https://arxiv.org/abs/2002.12543>.
- Sullivan et al. describe bounded exhaustive testing as exhaustive coverage
  of an input space through a stated size bound:
  <https://www.coppit.org/papers/issta_2004_bounded_exhaustive_testing.pdf>.
- The Hypothesis discussion of complete specifications states the essential
  boundary: a property set may characterize the intended function completely,
  while execution still checks only a finite set of examples:
  <https://hypothesis.works/articles/tests-as-complete-specifications/>.
