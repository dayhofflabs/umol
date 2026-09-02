# 219 — Path constraints in substructure matching

Status: Proposed
Date: 2026-09-01
Relates: [060](060-molecule-builder-dsl-2026-03-21.md),
[079](079-pattern-language-design-2026-04-10.md),
[080](080-unified-constraint-ast-2026-04-10.md),
[087](087-constraint-taxonomy-2026-04-17.md),
[165](165-ast-api-worklist-2026-07-27.md),
[193](193-subpattern-constraints-2026-08-09.md),
[195](195-molecule-constraint-matching-2026-08-12.md),
[197](197-deferred-dsl-features-2026-08-16.md),
[207](207-reaction-network-spike-2026-08-24.md),
[DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Design one molecule-level constraint that relates two pattern atoms through the host topology
between them: not bonded, same or different connected component, a path length in a range, and
later a path whose interior bonds and atoms satisfy patterns. The document settles the
graph-theoretic model, the placement under the homoiconicity requirement, the leaf and expression
shape that grows into regular path queries without a second leaf, and the two-stage evaluation that
makes the leaf usable in matching before the graph-core search learns about it.

Names introduced below are proposals. Doc 079 used `Path` and `PathExpr`; the constructor names
follow regular-expression vocabulary.

## Problem

Substructure matching is a monomorphism: every pattern bond maps to a host bond, and a pattern atom
pair with no pattern bond is unconstrained. A two-atom pattern with no bond therefore matches every
atom pair of the host, bonded or not, in one component or in two.

The reaction-network rules in doc 207 meet this directly. A bond-forming rule such as the carbanion
plus carbocation rule has a bondless left-hand side. It matches already-bonded pairs, and those
applications are rejected only at checked product publication, when the parallel bond fails
integrity and `apply` reports a structural conflict. A rule that closes rings only, or associates
two components only, cannot be written.

The constraint vocabulary has one pair-topology leaf, the connectedness leaf with an atom subset,
and negation around it expresses different components. Nothing evaluates it: the matcher rejects
every pattern with a non-empty molecule-scope list under the doc 195 gate. The intramolecular
predicate on noncovalent bonds derives same-component membership by reachability over localized
bonds and is the one existing decision that "connected" means the localized-bond topology. No leaf
expresses non-adjacency, a distance, or a path. Doc 079 sketched all three; doc 080 declined a
dormant path variant.

Graph-core provides connected components, biconnected blocks, a bounded breadth-first neighborhood,
and bounded simple-path enumeration. The subgraph-isomorphism search accepts one predicate per node
and one per edge and nothing indexed by a non-adjacent query pair.

## Model

### Matching preserves positive facts

A pattern and a host are relational structures over one signature: the atom relation, the
localized-bond relation, and one relation per overlay kind. A match is an injective homomorphism.
Homomorphisms preserve positive facts: a pattern bond becomes a host bond, a pattern aromatic system
becomes a host aromatic system. A non-bond is not a fact but the absence of one, and no homomorphism
preserves absence.

There are two ways to constrain absence. Strengthen the morphism globally, which is the induced
embedding; or add negated and derived atoms to the query, which is conjunctive queries with negation
over relations definable from the bond relation. This document selects the second.

### Monomorphism and induced embeddings

Under a monomorphism every pattern bond maps onto a host bond and pattern non-bonds say nothing.
Under an induced embedding every pattern non-bond additionally maps onto a host non-bond, so the host
bonds among the image atoms are exactly the pattern bonds. The propane pattern matches cyclopropane
as a monomorphism and not as an induced embedding.

Induced reads the pattern's bond relation as complete over the pattern's atoms. That is the
closed-world reading, and doc 195 fixed the open-world reading for matching: a pattern constrains
only what it mentions. Induced is also a property of the embedding rather than of the term and has
no truth value on a ground molecule, which is why it could only ever live as an operation parameter.
An induced match is a monomorphism with non-adjacency asserted on every pattern non-edge, so per-pair
predicates are its open-world form, and no embedding-kind parameter is adopted for substructure
matching. The nomenclature guide's embedding-kind entry describes the common-subgraph operations,
which are the only operations that take it.

### The shortest-path metric

Let d be the shortest-path distance over localized bonds, with unreachable pairs at infinity. The
metric layer of pair statements:

| Statement | In terms of d |
|---|---|
| bonded | d = 1 |
| not bonded | d ≥ 2 |
| same component | d < ∞ |
| different component | d = ∞ |
| separated by k to l bonds | k ≤ d ≤ l |

A monomorphism maps paths to paths, so it is non-expanding: the host distance between two images is
at most the pattern distance between the atoms. The pattern's own bonds therefore supply upper
bounds on host distance, and upper bounds need no predicate. Lower bounds, which "not bonded" and
"different component" are, can never come from the pattern's structure and must be written. Induced
embeddings reflect adjacency, and an isometric embedding would preserve d exactly; both sit on the
same axis, and a per-pair predicate is a two-sided bound on d for one pair.

### Regular path queries

The second layer relates two atoms by the label sequence of some path between them: a regular
expression over steps, where a step is a localized bond satisfying a bond pattern, with optional
tests on the atoms passed. Conjunctive queries whose atoms are regular path queries are the standard
graph-database formalism, from Mendelzon and Wood (1995, *Finding regular simple paths in graph
databases*, SIAM J. Comput. 24) through Barceló (2013, *Querying graph databases*, PODS) to Cypher,
GQL (ISO/IEC 39075:2024), and SPARQL property paths.

The metric layer is the unlabeled case with negation: "some path of length k to l" is any bond
repeated k to l times, and a lower bound on shortest distance is the negation of an upper-bounded
existential. Regular path queries are positive statements, so absence comes from the boolean
combinators already in the constraint tree, never from the expression.

Paths are simple paths. Under walk semantics a lower bound on length is meaningless on an undirected
graph, since a walk may retrace one bond. Simple-path evaluation is NP-hard on general graphs and
irrelevant at molecular size with bounded lengths.

### Relations outside path queries

Two atoms lie on a common cycle exactly when they lie in one biconnected block of at least three
atoms (Whitney 1932). That relation is not a regular path query: a common cycle is two internally
disjoint paths, and a path query speaks about one path's labels and cannot make a path return
through a named atom once interior atoms are unbound. Biconnectivity is definable in monadic
second-order logic, which is why graph databases do not offer it. It is a derived pair relation with
its own evaluator over the block decomposition graph-core already has.

Same block and same ring differ chemically. In naphthalene every atom shares a block with every
other, since all lie on the ten-membered perimeter cycle, while atoms of the two six-rings share no
ring. Same ring means some relevant cycle contains both, which depends on the ring model and is
what the ring-count predicates already use; a size-qualified form is a lookup in the same ring set.
These are further pair leaves with their own evaluators, sharing the pair index and matcher
admission designed here and nothing else. They are outside this document's scope.

## Homoiconicity

Doc 060 and the specification state the requirement: a ground term is a degenerate pattern, one
grammar, differing by evaluation context. The test that decides a placement is whether the predicate
has a truth value on a ground term by itself.

The existing derived predicates show the scheme. The degree predicate written on an atom is a claim
about that term's own topology. On a ground molecule it has a truth value discharged by derivation.
On a pattern it is evaluated at the image in the host, and the pattern atom's own degree is not
consulted. Same syntax, same meaning, evaluated on whichever structure the term is.

Pair predicates are that scheme at arity two. The index is a pair of the term's own atoms, never
host atoms and never a match parameter. On a ground molecule "atoms 3 and 5 are joined by a path of
two to four bonds" is a derived fact, redundant or false, exactly as a degree predicate is. On a
pattern it constrains the images. The derived relations are definable from the bond relation and
add nothing to the signature, so ground and pattern remain structures over one signature, one
complete and one partial.

Doc 079 settled the home: the constraints list is the general-arity location, and the inline atom
and bond forms are its arity-one sugar. The connectedness leaf is this predicate at arity n. No
different shape from the molecule and pattern pair is needed.

Three placements fail the test and are rejected. A pattern-only overlay entity has no ground
meaning. A parameter of the match operation is not on the term. An induced flag has no ground truth
value.

A path predicate whose degenerate case is one step of any bond states that a bond exists between
the images. That has the same truth as a bond entry, produces no bond correspondence, and cannot be
edited. The distinction between entity and predicate is the existing one between mapped and
filtered, and unbound interior atoms are what keep a path a predicate.

The intramolecular predicate on noncovalent bonds is this pair predicate placed on an overlay
instead of on the atom pair. It is subsumed rather than kept alongside; its removal is a separate
decision.

## Design

### Leaf and expression

The leaf is the path leaf from the start, and its expression type is the recursive algebra from the
start, populated with the constructors the metric layer needs.

```rust
pub enum MoleculeConstraint {
    // existing variants unchanged
    Path {
        from: AtomId,
        to: AtomId,
        expr: PathExpr,
    },
}

pub enum PathExpr {
    /// One localized bond satisfying the pattern.
    Bond(BondForm),
    /// The inner expression repeated; the repetition count satisfies the form.
    Repeat(Box<PathExpr>, NumForm),
}
```

The metric layer in this vocabulary, with `Not` from the constraint tree:

| Statement | Constraint |
|---|---|
| bonded | `Path(Bond(any))` |
| not bonded | `Not(Path(Bond(any)))` |
| same component | `Path(Repeat(Bond(any), 1..))` |
| different component | `Not(Path(Repeat(Bond(any), 1..)))` |
| some path of length k to l | `Path(Repeat(Bond(any), k..=l))` |
| shortest distance at least k | `Not(Path(Repeat(Bond(any), 1..=k-1)))` |

The later constructors are additions to the same enum: `Atom(AtomForm)` as a zero-length test on
the atom reached, `Concat(Vec<PathExpr>)`, and `Alt(Vec<PathExpr>)`. That is the nested regular
expression with node tests. Nothing written for the initial fragment moves, and no second leaf
appears: a distance leaf, a non-adjacency leaf, or a pair form of the connectedness leaf would each
duplicate an expression above.

Three properties carry the extension:

- **Recursive from day one.** A flat step-with-bounds struct is smaller today, but adding
  concatenation turns it into the `Repeat(Bond(..))` case of a recursive enum and rewrites the
  payload type, its serialization, and its Python class. The two-constructor enum is already the
  final type.
- **Steps reuse the pattern bond type.** A step's `BondForm` is the value a pattern bond entry
  carries, so a path step and a bond entry have one matching semantics, including derived
  predicates such as ring membership on the step bond, and no second bond-pattern grammar exists.
- **The count slot is a `NumForm`.** It reads, serializes, and matches like every other numeric
  slot and inherits the same deferred treatment of variable and expression forms. The alternative
  is an explicit minimum and optional maximum in regular-expression style, which simplifies the
  evaluator at the price of a second way to write a numeric constraint. Open below.

### Semantics

- **Endpoints are ordered.** An expression with concatenation reads from `from` to `to`.
  Normalization may swap the endpoints only by reversing the expression, which is the identity for
  the initial fragment. `from` and `to` must be distinct atoms of the term; distinctness and
  reference validity are tier-1 integrity.
- **Paths are simple paths in the host, blind to the embedding.** Interior atoms are distinct from
  each other and from both endpoints, are not part of the correspondence, and may coincide with
  images of other pattern atoms. This keeps the predicate binary in the two endpoints, which is
  what stage two below depends on.
- **Scope is the host.** Evaluation reads the host's localized bonds at the images. The pattern's
  own topology is not consulted, exactly as for the degree predicate.
- **Topology is localized bonds.** Overlays do not carry paths. This follows the topology
  definition in the nomenclature guide and the existing intramolecular derivation.
- **Negation is tree-level.** The expression algebra has no complement; `Not` wraps the leaf.
- **Lower bounds alone.** `Repeat(Bond(any), k..)` with k above one asks for some simple path of
  length at least k, a long-path question rather than reachability. It is well defined and
  computable at molecular size, and the shortest-distance readings above avoid it. The operative
  shapes are bounded counts and the open range from one.

### DSL surface

One production joins `molecule-constraint` in specification §7.12. A structured form with one map
per constructor and the count as a value-expression string is the proposal:

```
molecule-constraint ::= ...
  | { :path { :from atom-ref :to atom-ref :expr path-expr } }

path-expr ::= { :bond bond-string }
            | { :repeat [ path-expr value-expr ] }
```

The bond-string is the §7.4 form, so `{:bond "*"}` is any bond and `{:bond "1#R+"}` is a single
ring bond. Later constructors add `{:atom atom-string}`, `{:concat [path-expr+]}`, and
`{:alt [path-expr+]}`. The keyword names and the shape of the repeat form are open. Reading and
rendering follow the DSL serialization rules for molecule-scope leaves; structural atom refs are
accepted as everywhere.

### Evaluation

The evaluator is one function of the host, two host atoms, and the expression:

```rust
fn path_satisfied(
    host: &Molecule,
    from: AtomId,
    to: AtomId,
    expr: &PathExpr,
    bond_satisfies: &impl Fn(BondId, &BondForm) -> bool,
) -> bool
```

`bond_satisfies` is the per-key satisfies closure the matcher already builds for pattern bonds from
its constraint tables and ring context, so step semantics are borrowed. For the initial fragment the
body is a bounded depth-first search over satisfying bonds, with plain reachability when the count
has no upper bound and an adjacency lookup for a single step. Host component labels are built once
per match run for the reachability case, the way the ring set is built once when a ring key occurs.
When discharge of molecule-scope constraints on ground molecules is designed under doc 195, this
function is its evaluator for the path leaf.

**Stage one, graph-IR only.** The blanket gate in `visit_substructure_matches` becomes an admission
walk. A molecule-scope list whose entries are path leaves, alone or under `Not` and `And`, is
admitted; any other entry returns the existing rejection. Admitted trees are evaluated in the
visitor after overlay verification by calling the function above on the images. This is a local
change to the two strategy functions, and the incidence strategy is covered without further work
because the function reads localized bonds on the host regardless of the graph searched.

**Stage two, graph-core pair check.** The subgraph-isomorphism search gains a pair check beside the
node and edge predicates: a set of constrained query pairs and a callback over their two images,
invoked at the extension step when the second endpoint of a pair becomes mapped.

```rust
pub fn visit_subgraph_isomorphisms<B, F>(
    &self,
    query: &Graph,
    constrained_pairs: &[(NodeId, NodeId)],
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    pair_match: &mut impl FnMut(NodeId, NodeId, NodeId, NodeId) -> bool,
    alg: SubgraphIsomorphismAlgorithm,
    visitor: F,
) -> ControlFlow<B>
```

That is the loop each algorithm already runs over edges to mapped neighbors. It touches neither
vertex ordering nor candidate generation, so every selector takes it uniformly, and domain reduction
in ArcMatch may ignore the pairs soundly. The callback is the evaluator with the host bound.
Graph-IR compiles the admitted tree into per-pair predicates, conjunctions and negations of leaves
on one pair, and anything that does not factor into pairs, such as a disjunction across two pairs,
keeps using stage-one post-verification. Nothing at the IR, DSL, or Python level changes between
the stages, and the post-verification path remains as the permanent fallback. Whether the pairs are
passed as a slice or as a second graph over the query nodes is open.

**Later, domain restriction.** A positive pair predicate restricts the second endpoint's candidates
once the first is mapped: the first image's component for reachability, a ball of bounded radius
for a bounded count. This is the move the RDKit VF2 variant makes for bonds, it is where the
disconnected rule left-hand sides in doc 207 stop seeding the second atom from every host atom, and
it lives inside candidate generation per algorithm. It is separable from both stages and not
designed here.

### Reactions

A path constraint on a reaction left-hand side is a precondition. It is evaluated by matching and is
not transported to the product; the product's constraints are the host's own. `apply` matches
through `substructure_matches`, so stage-one admission is sufficient for the doc 207 catalogs. The
six tests ignored under the doc 195 gate stay ignored until the leaf kinds they carry are admitted.

The pair leaf completes a vocabulary the bond side already has. The inverse of a bond formation with
a different-component precondition is a bond removal whose bond carries the not-in-ring predicate,
since a bond whose removal disconnects is exactly a non-ring bond. The inverse of a ring closure
with a same-component precondition is a bond removal on a ring bond. Doc 207's reversibility
invariant can state both.

## Consequences

- **Specification.** §6.2 changes from rejecting every non-empty molecule-scope list to rejecting
  entries outside the admitted set; §7.12 gains the production above; §6.1 lists the path leaf
  among derived predicates. Admission under the doc 197 policy, normative with conformance coverage
  or experimental, is open.
- **Graph IR.** One variant, one expression type, one evaluator, the admission walk, and a sibling
  arm in every match over `MoleculeConstraint` that the compiler enumerates: reference validation,
  distinct-endpoint integrity, compaction, remapping, normalization, the atom-set cascade used by
  edits, and the reaction span.
- **Graph-core.** Stage two only; nothing in stage one.
- **Python.** The molecule constraint class gains the variant and the expression type with
  construction parity.
- **Nomenclature guide.** The embedding-kind entry is corrected to the common-subgraph operations;
  an entry for the path constraint and its simple-path semantics is added.
- **Doc 197.** The molecule-scope matching entry records that the path leaf is admitted and the
  rest remains deferred to doc 195.
- **Doc 195.** The evaluation surface it designs receives its first leaf; scope and evaluation
  point are settled here for this leaf only.
- **Doc 165.** The optional atom-subset scope item is unaffected; the fixed pair index does not use
  that shape.
- **Tests as specifications.** The constrained match set equals the unconstrained match set minus
  the embeddings the host relation excludes; `Not(Path)` selects exactly the complement within the
  unconstrained set; a single any-bond step holds exactly when the images are bonded; reachability
  is symmetric and never lost under bond addition; a bounded count never gains under bond removal.

## Rejected alternatives

- **A length form over an extended domain with an unreachable value.** Duplicates
  `Repeat(Bond(any), _)` under negation and adds a lattice type.
- **Boolean leaves for not bonded and same component.** Smallest today, no seam to counts or
  paths, and the third leaf arrives anyway.
- **A path kind in the bond relation.** Reads like Cypher, but places a non-bond in the topology
  that every consumer of bonds would skip, against the topology definition and the doc 079 rule
  that constraints are the escape hatch.
- **A ninth entity kind.** Has no ground meaning and costs deltas, spans, remapping, and bindings.
- **An embedding-kind parameter on substructure matching.** Closed-world, an embedding property
  rather than a term property, and subsumed by per-pair predicates.

## Open questions

- Leaf keyword, expression constructor names, and the EDN shape of the repeat form.
- Count slot as `NumForm` or as an explicit minimum and optional maximum.
- Admission of the leaf to the normative specification with conformance coverage, or experimental
  status under doc 197.
- Whether the existing connectedness leaf is evaluated in stage one with the same component labels,
  retained as the n-ary convenience, or folded into the path leaf.
- Whether the rejection error for non-admitted entries keeps its count or names the construct, as
  doc 195 asks.
- The stage-two carrier for constrained pairs: a slice of pairs or a second graph over the query
  nodes.
