# 223 — Formalizing SMARTS and SMIRKS

Status: Informational
Date: 2026-09-05
Relates: [079](079-pattern-language-design-2026-04-10.md),
[115](115-variable-facility-2026-06-16.md),
[132](132-reaction-ast-implementation-plan-2026-06-25.md),
[185](185-python-reaction-span-2026-08-04.md),
[193](193-subpattern-constraints-2026-08-09.md),
[195](195-molecule-constraint-matching-2026-08-12.md),
[197](197-deferred-dsl-features-2026-08-16.md),
[219](219-path-constraints-2026-09-01.md),
[DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md),
[OpenSMILES formalization](../umol-io/spec/opensmiles-spec.md),
[integrity guide](../docs/development/integrity.md)

## Aim

Explore a formalization of SMARTS and SMIRKS analogous to the repository's OpenSMILES
formalization: precise lexical and syntactic rules, explicit semantics for each construct,
and declared decisions where existing definitions are ambiguous or unsuitable. The ambition
is to give the pattern languages a coherent meaning in umol's existing molecular and reaction
semantics, rather than define their meaning by a particular toolkit's execution.

**This looks feasible, and umol already supplies the foundation.** Molecules and molecular
patterns share a representation; so do reactions and reaction patterns. The native model need
not adopt SMILES defaults or SMIRKS restrictions, but the boundary formalization should preserve
existing language definitions wherever they can be interpreted faithfully in that model.
The question is how to formalize the languages against that foundation, not whether umol
needs a new pattern model.

The investigation may identify an explicitly documented variant rather than a faithful account
of every historical dialect. That is acceptable; the smaller fallback is a SMARTS/SMIRKS
boundary representation with a well-defined translation into existing Molecule and ReactionSpan
objects, with Reaction supplying the operational form. This record is exploration, not a normative specification, API design, or implementation
plan. It does not impose a restriction to ground products or chemically conservative reactions.

## Governing rule for compatibility

Hew as closely as possible to existing SMARTS/SMIRKS definitions. A surface convention is not
an invitation to choose a preferred alternative: where its definition is clear and fits the
semantic model, inherit it and specify its interpretation. This applies to operator precedence,
shorthands, omissions, and other source-language distinctions, even when another spelling would
be simpler.

Deviate only when an existing definition prevents clear parsing or cannot be accommodated by
umol's current semantics or extensions without breaking its structural concepts. Missing parser
or evaluator support alone is not such a conflict. Extensions may supply needed observations
or constraints while preserving the existing distinction between entities, asserted properties,
correspondences, and reaction spans.

For each proposed deviation, identify the source rule, the concrete parsing or structural
conflict, why faithful interpretation does not resolve it, and the replacement meaning. Earlier
deviation candidates must meet this criterion; preference alone does not settle them. Where
sources leave a question open or dialects disagree, state that uncertainty and the adopted
interpretation explicitly.

For every question, first check whether known umol semantics already answers it. Then establish
the existing source-language definition and its translation, and isolate any actual incompatibility. It should not ask the user to redesign
already defined syntax one construct at a time.

## Background and established foundation

Read the whitepaper source at `/Users/dr/Documents/paper/umol/main.tex`, especially Basic
Principles, Identity and Constraints, Attribute Lattice, Mutation, and Reactions. Its central
construction is already the answer to much of the representation question:

- Atoms, localized bonds, and overlays carry lattice-valued attributes. A concrete value,
  a set, and a wildcard belong to the same attribute domain.
- A pattern attribute p matches t when p ∧ t = t. Compatibility is a different question:
  p ∧ t ≠ bottom. Topological occurrence is supplied by subgraph matching.
- Derived properties are projections of the graph and overlays. A pattern can assert a
  projection without reproducing the structure from which a host derives it. Aromatic
  membership does not require the pattern to contain an entire aromatic system.
- Reaction is an lhs plus deltas, with lattice-valued definitions on both. Changes can
  create, delete, and modify entities; their numbers need not be conserved. Chemistry
  validation and resolution are separate operations. Reaction application does not
  automatically repair hydrogen counts or impose a valence model.

These are starting points, not proposed additions. The whitepaper's pyridine-type pattern,
`N#a1#h0#R(6)1`, directly illustrates the asserted/derived relationship relevant to SMARTS.
Its distinction between aromatic nonmembership and zero electron contribution also matters:
`#a!` and `#a0` are different values, so aromaticity cannot be reduced to positive contribution.

The [DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md) and current
[substructure](../umol-graph-ir/src/ir/substructure.rs),
[constraints](../umol-graph-ir/src/ir/constraint.rs), and
[reaction](../umol-graph-ir/src/ir/reaction.rs) code were inspected at umol commit
`cd902742b20eec4cf4aac65c0d0e5afa8e93fee9`. They substantiate the shared representation.
Some evaluator coverage remains separate: molecule-scope constraint matching rejects those
constraints today, and recursive matching remains the work of docs 193 and 195. Those limits
do not make patterns impossible or determine what a language formalization may describe.

The OpenSMILES document supplies the useful organizational precedent: lexical rules, grammar,
semantic interpretation, deviations, diagnostics, and conformance examples. Its unfinished
semantic sections and some historical implementation notes mean it is a precedent for the
attempt, not a claim that all its contents are complete or current. This exploration does
not revise that document.

## What the external languages actually distinguish

Daylight defines molecular SMARTS, reaction SMARTS for querying reaction records, and SMIRKS for
transforms. It explicitly notes that SMILES and SMARTS defaults differ. Unspecified SMARTS
properties impose no restriction. Recursive SMARTS anchors an inner query at its first atom;
component grouping adds restrictions beyond disconnected query topology.
[Daylight SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)

Daylight SMIRKS restricts the query language to make changes interpretable: pairwise mapped atoms,
concrete bond expressions, and restrictions on atoms whose bonding changes. Its explicit-hydrogen
rules address the distinction between hydrogen atoms and hydrogen-count predicates. Thus the
query/construction tension is present in the original design.
[Daylight SMIRKS](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)

RDKit documents a different transform dialect under “reaction SMARTS,” including product mapped
dummies that copy reactant atoms and product any-bonds that copy corresponding reactant bonds.
These are constructive operations, not ordinary wildcard predicates. There is documentation;
there is no reason to treat it as the universal semantics of reaction patterns.
[RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)

Local source inspection reinforces the separation of syntax, predicates, and execution:

| Checkout and inspected files | Evidence relevant to this spike |
| --- | --- |
| RDKit `28fbe8a69b6debf9bc7d12fc0506757932be17b1`; [smarts.yy](../materials/codes/rdkit/Code/GraphMol/SmilesParse/smarts.yy), [ReactionRunner.cpp](../materials/codes/rdkit/Code/GraphMol/ChemReactions/ReactionRunner.cpp) | Parser actions build Boolean and recursive query objects; the reaction runner separately handles product properties, missing values, and query bonds. Parsing does not supply the entire transform semantics. |
| CDK `b2996d233d21aea840513ba0f77cc74ae2e30e9a`; [Smarts.java](../materials/codes/cdk/tool/smarts/src/main/java/org/openscience/cdk/smarts/Smarts.java), [Expr.java](../materials/codes/cdk/base/isomorphism/src/main/java/org/openscience/cdk/isomorphism/matchers/Expr.java), [Smirks.java](../materials/codes/cdk/tool/smarts/src/main/java/org/openscience/cdk/smirks/Smirks.java) | Explicit flavors change interpretation, including degree and ring-count conventions. Expressions distinguish total/implicit H, degree, ring count, and ring size. Recursive matching is rooted; transforms are compiled separately. |
| Open Babel `889c350feb179b43aa43985799910149d4eaa2bc`; [parsmart.cpp](../materials/codes/openbabel/src/parsmart.cpp) | Atom and bond expression evaluators call host observations such as aromaticity, explicit/total degree, hydrogen counts, and ring membership. The host model is part of predicate meaning. |

These are source observations, not executed differential results. The three implementations
have not been established equivalent. Their versions and host preparation would need pinning
for any future comparison.

## What would actually be formalized

The semantic account should be independent of the textual traversal, while specifying exactly
how that traversal denotes a pattern. A grammar alone is insufficient, but a new general
logic framework is unnecessary as a starting point.

| Layer | SMARTS | SMIRKS |
| --- | --- | --- |
| Lexical | Tokens, numeric domains, element names, operator tokens, bracket context. | Shared tokens plus reaction separators and mapping syntax. |
| Syntactic | Atom/bond expressions, precedence, branches, closures, grouping, recursive queries. | Participant sections, side-local graphs, map-label scope, permitted side expressions. |
| Structural interpretation | Query atoms and bonds, closure references, component conditions, stereo frames. | Shared versus side-only entities, changed incidence, and left/right frame correspondence. |
| Attribute interpretation | Inherent-field predicates and asserted projections, including Boolean expressions. | The same attribute meanings on each side, together with explicitly determined changes. |
| Operational meaning | What constitutes an occurrence and what each predicate observes. | What a rule application means under the existing delta and preservation semantics. |
| Conformance | Positive/negative examples and explicit dialect differences. | Examples distinguishing unchanged, modified, created, and removed entities and attributes. |

A semantic clause should say what property of a molecule or reaction makes an expression true
or what change it denotes. It should not say merely “call the toolkit's ring routine” or
“use its product sanitization.”

### Molecular clauses: a small illustrative fragment

For a chosen occurrence e of the query topology in a host M, interpret each atom or bond
expression against the corresponding host entity. A literal denotes a singleton in its value
lattice; a wildcard denotes top. An asserted derived property is compared against its host
projection using the same refinement rule as an inherent property.

Here M is the host molecular structure, e is a candidate structure-preserving correspondence
from query entities to host entities, and φ is a condition on those entities. The notation
`M, e ⊨ φ` means “φ holds in M under e”; `⊨` is the satisfaction symbol. For a query atom x,
`e(x)` is its host image. A, B, and φ below denote formulas, not atoms.

Boolean composition has the ordinary satisfaction clauses:

```text
M, e ⊨ A ∧ B    iff M, e ⊨ A and M, e ⊨ B
M, e ⊨ A ∨ B    iff M, e ⊨ A or  M, e ⊨ B
M, e ⊨ ¬A       iff not (M, e ⊨ A)
```

This fragment assumes that the observations involved have defined values on the host.
An unavailable observation is not false. For partially specified hosts, the specification
must distinguish refinement from compatibility rather than silently switch between them.

The following summarize interpretation clauses. Decisions recorded below settle part of this
surface; the remaining entries identify semantics still to specify.

| Construct | Semantic reading to specify |
| --- | --- |
| `[#6]` | Carbon element, independent of aromatic membership. |
| `[C]`, `[c]` | Carbon plus, respectively, nonmembership or membership in an aromatic system. The membership predicate includes zero-contribution members. |
| An omitted atom property | No assertion about that property; no obligation to import a SMILES completion default. |
| `D`, `X`, `h`, `v` | Precisely named degree, total-degree, hydrogen, and valence observations. Uppercase H counts total hydrogens; lowercase h counts implicit hydrogens, as detailed below. |
| `R`, `r`, `x`, bond `@` | Precisely defined ring observations. A ring-closure digit introduces connectivity, not a ring-membership assertion. |
| `~`, `-`, `:`, omitted bond | Explicit predicates on a host bond and its projections; any legacy shorthand expansion is a language decision. |
| `C.C` | Two distinct query atoms with no query bond requirement; host images may be adjacent, farther apart, or disconnected. |
| Component grouping | Explicit restrictions on which host components contain the selected entities. |
| `@`, `@@`, directional bonds | Configuration constraints interpreted in the traversal-induced participant frame and transported to the host frame. |

Boolean scope matters even though ordinary lattice-valued patterns already exist. For example,
`(element = C ∧ charge = 0) ∨ (element = N ∧ charge = +1)` must retain that correlation.
Replacing it by independent element and charge sets admits extra cases. Formalizing arbitrary
SMARTS Boolean expressions therefore requires specifying their composition, not claiming that
patterns themselves need to be invented. For the same reason, operator precedence and the
special lexical treatment of hydrogen deserve explicit grammar productions.

For recursive SMARTS, a clause can be stated without defining general recursion:

```text
M, e ⊨ $(Q) at x
    iff there is an occurrence f of Q in M with f(root(Q)) = e(x).
```

Here Q is the nested query, root(Q) its distinguished first atom, and f a correspondence
for that nested query. Thus f is a mapping, not a formula.

Each nested query has its own entity scope and is internally injective. Both roots map to
the outer atom being tested. Non-root witnesses may overlap each other and images of other
outer query entities; no additional disjointness constraint is imposed. This is the same
independence from the outer correspondence used for simple-path interiors in doc 219.
Negation means no such rooted occurrence exists. Finite nesting does not require a fixed-point
semantics. This identifies what doc 193 would eventually implement; it does not settle that
feature's public representation.

### Reaction clauses: ReactionSpan is the semantic counterpart of SMIRKS

ReactionSpan already represents the two-sided meaning of SMIRKS input. It materializes the
superimposed left and right structures with shared entity identity and per-entity before/after
forms. Reaction is the operational representation of that meaning as an lhs plus deltas.
The distinction is between two existing representations, not a missing reaction-pattern facility.

This role was anticipated. [Doc 132](132-reaction-ast-implementation-plan-2026-06-25.md)
explicitly identifies the span as the landing point for mapped reaction SMILES, SMIRKS, and GML,
followed by conversion to lhs-plus-deltas. [Doc 185](185-python-reaction-span-2026-08-04.md)
later motivates exposing the span and its construction paths as bridges from existing rule
representations. The whitepaper's two-sided presentation describes the same relationship.

The semantic route to formalize is therefore:

```text
SMIRKS boundary expression → ReactionSpan → Reaction
                            two-sided       lhs + deltas
```

This is a semantic decomposition, not a requirement for a parser to allocate each intermediate.
[ReactionSpan](../umol-graph-ir/src/ir/reaction_span.rs) represents each entity as Unchanged,
Modified, Added, or Removed. Unchanged and Modified entities persist across sides; Removed
entities occur only on the left and Added entities only on the right. Its lhs and rhs
projections are Molecule values and may carry nonliteral forms. Its to_reaction conversion
provides the operational representation. It is a semantic object, not a carrier for every
source-format annotation or spelling.

The formalization needs to specify how the input denotes this span: which entities are shared,
which side forms they carry, and how participant frames correspond. The resulting changes
then follow through the existing span-to-reaction conversion, rather than a second,
SMIRKS-specific reaction mechanism. Historical restrictions that simplified a particular
transform interpreter need not be imported.

At minimum, it must distinguish:

- Identity shared between sides, as declared by mapping, from a mere equality of attribute values.
- An entity present only on one side from a preserved entity with changed attributes.
- A property used to select the lhs from an instruction changing that property.
- A source RHS wildcard with preservation semantics from a native undetermined form to assign.
  Mapped RHS * and corresponding-bond ~ use the preservation reading recorded below; omissions
  require their own source-language clauses.
- A changed stereo glyph from a changed configuration after accounting for participant order.

The deltas may themselves contain nonliteral forms. There is no requirement here that every
rule produce a uniquely resolved concrete product. A well-defined transformation of a partial
structure remains a meaningful reaction; resolving or validating the result is a separate request.
What must be unambiguous is the transformation denoted by the expression, not every attribute
of the resulting structure.

A nonliteral RHS describes the admissible product values. A set-valued product and its concrete
alternatives have the same denotation, with the same surrounding structure, preservation conditions,
and correlations. Assigning that symbolic form and requiring the corresponding product condition
are not different denotational semantics. Symbolic storage, resolution, and enumeration concern
how the result is represented and computed. ReactionSpan supplies the two-sided meaning; its
conversion to Reaction supplies the operational representation. The mechanics of carrying and
materializing alternatives through application remain to be worked out.

Questions such as unmapped atoms, repeated mapping classes, agents, omitted fields, and any-bonds
on the RHS need their own clauses. Preserve, set, and unconstrained must not be conflated by
implementation convention. Daylight reaction-query map classes and RDKit product-copy behavior
are evidence for these distinctions, not mandatory limitations of the umol interpretation.
Reaction SMARTS as a query on existing reactions should be distinguished explicitly from SMIRKS
as a transformation, even if some syntax is shared.

The graph-and-overlay preservation rules are already part of umol. A SMIRKS translation has to
express its changes in those terms, including incident overlays. It need not silently perceive,
repair, or chemically validate the product. If legacy notation omits information needed to
identify the intended overlay change, that is a specific translation ambiguity to document.

## Variables and paths as related semantic work

### Value bindings and entity correspondence

[Doc 115](115-variable-facility-2026-06-16.md) identifies three separate properties of a
variable: its value domain, identity, and binding scope. It places correlations at the lowest
scope containing their variables and distinguishes anonymous field-local bounds from deliberate
reuse of a named variable. Its proposed facility is not implemented by virtue of being described;
the index marks it Proposed, and its historical implementation details should not be read as
current behavior. For example, current NumForm has explicit range forms, so the document's old
fixed-name encoding of anonymous bounds is not the present representation of every bound.

This distinction sharpens the satisfaction notation. If named value variables are admitted,
write `M, e, β ⊨ φ`, where β assigns values to typed variables in the current scope. The
correspondence e maps entities; β binds attribute values. A reaction map label declaring that
an entity persists across sides is not automatically an attribute variable. Nor does a repeated
anonymous wildcard introduce a shared binding.

For example, two fields each constrained to {1,2} are independent. Binding both to the same
variable n restricts them to (1,1) or (2,2). A reaction expression referring to an lhs-bound n
also needs a rule for carrying that binding into its delta expressions. These are questions
about the meaning of shared names, not about choosing a solver or a variable-environment API.

Nested queries make scope consequential: which names are local, whether outer variables may be
referenced, and whether inner bindings can escape. Existentially hidden inner variables differ
from values exported as part of a match. Negation must quantify over the intended local bindings;
it cannot create a binding by failing to find one. The agreed
[scope rules in doc 115](115-variable-facility-2026-06-16.md#nested-query-scope-decision-2026-09-05)
use lexical ownership, read-only enclosing bindings, and existential locals with no export.
They initially reject shadowing and restrict evaluation to values supplied by parameters or
positive structural matches. Doc 115 records the SQL, SPARQL, and Datalog precedents and
distinguishes that evaluation restriction from the denotation of symbolic patterns.

Named value variables need not become new SMARTS syntax. They are relevant both as a precise
account of binding during interpretation and as a possible extension if the desired variant
admits shared-value constraints. Any such extension must be distinguished from the legacy
constructs being formalized.

### Path predicates and their ground meaning

[Doc 219](219-path-constraints-2026-09-01.md) makes the homoiconicity criterion explicit:
a predicate has a truth value on a ground term's own structure, and a pattern asserts the
same predicate at its images in the host. Its proposed path constraint relates two of the
term's atoms through localized-bond topology. It is not a new bond or overlay entity.

This directly informs component grouping and negative connectivity conditions. A missing query
bond says nothing; non-adjacency and different-component requirements are explicit predicates.
An asserted path filters an occurrence without adding its interior atoms or bonds to the
correspondence. Those interior atoms may coincide with images of other query atoms. A drawn
chain, by contrast, introduces entities that are mapped and may be addressed by reaction deltas.

Doc 219 selects simple host paths and distinguishes an existential path length from shortest-path
distance. A triangle has a two-edge simple path between adjacent atoms, but their shortest distance
is one. Consequently “some path has length at least two” does not express “not bonded.” This
is exactly the sort of distinction a language semantics must make independently of its evaluator.

Its proposed regular path expressions reuse ordinary bond patterns and later atom tests, with
negation outside the path expression. This offers a related extension of the same pattern
vocabulary, rather than a reason to reinterpret recursive SMARTS as arbitrary-length paths:
a rooted subquery describes a finite graph; a repeated path describes a family of paths.
Component conditions connect to reachability immediately, while general regular-path syntax
would be an explicitly identified extension.

A path condition on a reaction lhs is a selection condition, not an instruction to copy its
witness into the product. If the formalism later permits binding a path length or referring to
an interior entity, that crosses into doc 115's binding and scope questions and needs an explicit
meaning. Neither feature should be inferred from an existential path predicate alone.

## Where clean semantics and compatibility can diverge

The relevant question for each construct is whether its meaning follows from the existing
model, requires a coherent extension, or conflicts with clear parsing or existing structural
concepts. A convention is not declined merely because an alternative seems preferable.

Ring predicates are a useful example. The whitepaper defines the native reference as Relevant
rings through size 22; legacy dialects may use other counts or minimum-ring conventions. The
chosen formalization adopts the native reference. Daylight names SSSR in its primitive table,
so this is an explicit choice of ring semantics, not a claim that the manual contains no ring
requirement. It avoids dependence on an unspecified choice among minimum cycle bases. Exact
primitive meanings and behavior outside the reference cutoff still need spelling out.

Aromatic queries are another. The paper already explains how asserted atom/bond projections
avoid specifying complete overlays. The remaining language question is exactly which projection
an aromatic token denotes, including how it behaves on structures beyond conventional organic
aromaticity. Similarly, degree and hydrogen-count predicates can intentionally be sensitive to
explicit versus implicit hydrogen representation; invariance should be claimed only where the
chosen observation warrants it.

Default values, restrictions to SMILES-like reaction-center atoms, and mapping conventions
belong to the boundary language rather than axioms of the native model. Their boundary meanings
should be retained unless the governing deviation criterion is met. Once deviations are substantive, the specification should name itself as a defined
SMARTS/SMIRKS variant rather than imply universal toolkit compatibility.

## Decisions from discussion — 2026-09-05

These settle the stated cases for this formalization; they do not claim compatibility with
every dialect or fill in unaddressed syntax by inference.

| Construct | Decision |
| --- | --- |
| `[#6]`, `[C]`, `[c]` | Carbon alone, carbon with aromatic nonmembership, and carbon with aromatic membership, respectively. Other fields remain unconstrained. |
| Boolean precedence | Retain Daylight: !, then & or juxtaposition, then comma, then semicolon. `[C,N+]` is C or (N and +1); `[C,N;+]` is (C or N) and +1. |
| `[H]`, `[#1]` | Explicit hydrogen atom. Hydrogen cannot belong to an aromatic system in this formalization. |
| `[XH]`, `[#zH]`, `[XHk]`, `[#zHk]` | Total hydrogen count on the identified element: k hydrogens, or one when the count is omitted. X here is a metavariable for an element symbol, z its atomic number, and k the count; X is not the SMARTS connectivity primitive in these schematic examples. |
| `C.C` | No query bond requirement; 0 < d ≤ infinity for the two distinct host images. |
| `(C.C)` | 0 < d < infinity: distinct atoms in the same component, possibly adjacent. |
| `(C).(C)` | d = infinity: disconnected. |
| `[$(*C);$(*CC)]` | Both predicates have the outer atom as their root. Their remaining witnesses may overlap each other and other matched entities; each occurrence remains internally injective. |
| `[!$(C=O)]` | No rooted C=O occurrence: the root fails C, or has no double-bonded neighbor satisfying O. This negates existence of a matching neighbor, not each neighbor independently. C and O retain their aliphatic meanings. |
| Ring reference | Use umol's native reference: Relevant rings through size 22. |
| `Rn` | The atom belongs to exactly n reference rings. |
| `rn` | The smallest reference ring containing the atom has size n. |
| `xn` | Exactly n incident localized bonds belong to reference rings. |
| `c:c` | The colon asserts aromatic membership on the bond, through the existing bond #a projection; the endpoint c expressions also assert carbon and aromatic membership. No whole aromatic overlay is constructed by the query. |

Here d is shortest-path distance over host localized bonds between the two selected atoms;
infinity denotes no path. Distinctness follows from injective matching. A dot introduces no
query bond and asserts no host non-adjacency, including between larger fragments. Grouping
adds component restrictions. Non-adjacency requires an explicit constraint, as in doc 219.

### Hydrogen semantics and syntactic disambiguation

Total-H semantics is accepted: it is the existing assertable and checkable #H constraint,
counting implicit hydrogens plus explicit hydrogen neighbors. Lowercase h denotes the separate
implicit-H count (#h). The earlier proposed change of uppercase H to implicit-only counting is
withdrawn. The concern was the overloaded spelling, not the derived total-H observation.

The inspected RDKit and CDK parsers resolve that overloading by recognizing a special complete
bracket form for an explicit hydrogen atom: optional isotope, H, optional charge, optional map,
then the closing bracket. Outside that form, H in an atom expression is a total-H-count
predicate, with an optional numeral defaulting to one. Atomic number #1 remains an unambiguous
way to write an explicit hydrogen atom inside a compound Boolean expression.

| Expression | Reading in the inspected RDKit/CDK grammar |
| --- | --- |
| `[H]`, `[2H]`, `[H+]`, `[H:1]` | Hydrogen atom, optionally isotope-qualified, charged, or mapped. |
| `[H1]`, `[*H1]` | Any atom with exactly one total attached hydrogen. |
| `[CH]` | Aliphatic carbon with exactly one total attached hydrogen. |
| `[H;+]` | Any atom with one total attached hydrogen and charge +1. The semicolon prevents the special hydrogen-atom production. |
| `[#1;+]` | Hydrogen atom with charge +1. |
| `[h1]`, `[h]` | Exactly one implicit hydrogen, or at least one implicit hydrogen, respectively. Lowercase h has no aromatic meaning. |

RDKit enumerates the special alternatives in hydrogen_atom in
[smarts.yy](../materials/codes/rdkit/Code/GraphMol/SmilesParse/smarts.yy); ordinary H_TOKEN
productions construct H-count queries. CDK's parseExplicitHydrogen in
[Smarts.java](../materials/codes/cdk/tool/smarts/src/main/java/org/openscience/cdk/smarts/Smarts.java)
tries the complete special form first and resets the input position on failure before ordinary
expression parsing. No host chemistry is consulted to choose between the meanings. These are
source-derived readings, not newly executed parser conformance cases.

Decision: retain the existing hydrogen parsing exception in the formalization, including the
overloaded spellings. Adding an operator can unexpectedly change H from an element into a count
predicate, but the grammar is disambiguatable and does not break the semantic model. State the
whole-bracket exception explicitly and interpret its results as the distinct hydrogen element,
total-H constraint (#H), and implicit-H field (#h). Total and implicit counts are not aliases.
This settles the language direction; parser implementation and conformance verification remain
future work.

### Mapped RHS wildcards

Decision: incorporate the broader dialect's preservation reading for mapped RHS * and for RHS ~
when a corresponding reactant bond exists. These source constructs retain the matched atom's
element and the matched bond's type, respectively; they do not assign native undetermined values.
For example, `[C:1]=[O:2]>>[C:1][*:2]` preserves the oxygen element while changing the bond.
The corresponding fields have the same before/after forms in ReactionSpan, producing no edit
for those fields. Other explicitly specified changes remain independent.

RDKit documents these behaviors in its
[reaction SMARTS section](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts).
They belong to the broader dialect: Daylight's strict SMIRKS restricts bond queries and atom
expressions at changing centers. This preservation reading fits the existing span model;
it is not a general rule that every native wildcard means copy. Unmatched RHS wildcards
instead have the open-attribute interpretation below.

### Unmatched RHS wildcards

Decision: an unmapped RHS * adds an atom with undetermined element; an RHS ~ with no
corresponding lhs bond adds a bond with undetermined type. The entity's existence is asserted;
the open attribute describes its admissible alternatives. Thus `[C:1]>>[C:1]-*` adds an atom
attached by a single bond, while `[C:1].[O:2]>>[C:1]~[O:2]` adds a bond between the preserved
carbon and oxygen. Other stated constraints and structural integrity still apply.

This follows the existing Added semantics and the decision on nonliteral RHS values:
matching checks the product constraints, and application denotes the admissible products.
It does not require choosing a concrete element or bond type immediately. Preservation
applies only when there is a corresponding lhs entity supplying the attribute to retain.

This is an explicit boundary interpretation for the chosen variant. RDKit's
[reaction examples](../materials/codes/rdkit/Docs/Book/RDKit_Book.rst) leave unmapped RHS
dummy atoms as dummies; here * means an undetermined element, not a distinguished dummy
element. RHS bond queries also extend beyond strict Daylight SMIRKS. No new transformation
semantics or enumeration mechanism is introduced by admitting these cases.

### Omitted RHS charge

Decision: preserve omitted RHS charge on a mapped atom whose element is unchanged; an explicit
RHS charge sets that value. RDKit's reaction-SMARTS path copies
the host charge when the product has no charge query, and uses an explicitly supplied charge
otherwise. The parser enables this behavior with the implicit-properties flag;
[Reaction.h](../materials/codes/rdkit/Code/GraphMol/ChemReactions/Reaction.h) documents the
omission convention, and updateImplicitAtomProperties in
[ReactionRunner.cpp](../materials/codes/rdkit/Code/GraphMol/ChemReactions/ReactionRunner.cpp)
implements it. This is a source-derived account, not an executed comparison or a claim about
every SMIRKS dialect.

Thus `[O-:1]>>[O:1]` retains charge -1 in this path, whereas `[O-:1]>>[O+0:1]` explicitly sets
charge zero. The first has equal charge forms across the span; the second has a charge change.
When the lhs does not constrain charge either, preserving its open charge form on both sides
produces no charge edit and leaves the actual host value intact. Independently parsing the
omitted RHS charge as zero or as a different form would lose this distinction.

This behavior fits ReactionSpan and is accepted for this case. It must not be generalized
to element changes, newly created atoms, or omitted hydrogen/stereo fields without checking
those cases. Daylight's mixed SMILES/SMARTS restrictions also require their own source reading;
the strict and broader dialects are not assumed identical.

### Model-dependent hydrogen completion

Decision: do not infer or perform hydrogen-count adjustments in situ. The implicit-H count
appropriate after a reaction is chemistry-model-dependent. Compiling another toolkit's assumed
valence completion into an explicit arithmetic delta would still impose that model; moving the
adjustment into the input translation does not make it model-independent. Named variables or
count expressions do not resolve this problem.

Where the input leaves the product implicit-H count to completion, it may remain undetermined.
Resolution can subsequently determine it under an explicitly chosen chemistry model. Reaction
application does not perform that resolution, repair hydrogens, or assume a unique valid product
count. Explicitly specified H fields and constraints retain their declared meanings.

Independent resolution of the two sides under a chosen chemistry model can subsume valence-based
hydrogen completion: each side determines its own count from its topology, fields, and assertions.
There is no additional chemistry in a count adjustment that cannot instead be expressed as the
difference between those resolved counts. The count inferred for the lhs must not accidentally
become an asserted rhs count when the rhs leaves it open.

This applies to complete side structures, including preserved host context when applying a
pattern. Resolving an isolated query fragment as though it were a complete molecule would infer
hydrogens at its open boundary and change its matching meaning. Explicit H atoms, total-H
constraints, and other supplied facts remain inputs to resolution rather than being replaced
by unconstrained valence saturation. A selected model can leave a side underdetermined or reject
its constraints; the formalism does not promise every side resolves uniquely.

Source inspection clarifies the earlier toolkit comparison. CDK's atomTypeOps in
[Smirks.java](../materials/codes/cdk/tool/smarts/src/main/java/org/openscience/cdk/smirks/Smirks.java)
emits AdjustH for changes in explicitly supplied total-H counts and for suppressed-H bookkeeping;
its calcImplH also performs built-in valence completion for unbracketed product atoms. These
are not all instances of a generic bond-order-driven H adjustment. RDKit marks differences in
mapped template degree in
[Reaction.cpp](../materials/codes/rdkit/Code/GraphMol/ChemReactions/Reaction.cpp).
Its updateImplicitAtomProperties copies the stored H count and no-implicit flag when degree is
unchanged and the RHS has no H-count query; its
[Atom.cpp](../materials/codes/rdkit/Code/GraphMol/Atom.cpp) separately calculates implicit valence
using its allowed-valence rules. Degree here counts neighbors; it is not bond-order sum.

Those toolkits therefore do have specialized completion machinery, though not umol's general
explicit resolution operation. The difference is where the model is selected and when completion
runs, not a fundamental need for in-situ arithmetic. With compatible models and the same explicit
facts, independent side resolution is the intended replacement for that completion. Exact toolkit
parity has not been established by executing comparisons. Model-independent reaction application
still performs only its specified changes. The exact source cases that preserve a count versus leave it undetermined,
and their representation through span conversion and application, remain to specify; this decision
does not equate an unchanged wildcard field with an instruction to clear a concrete host value.

### Entity creation and deletion

Decision: shared atom-map labels identify preserved entities, with independently specified
side attributes. Unmapped lhs atoms are Removed; unmapped rhs atoms are Added. Thus
`[C:1]Br>>[C:1]O` preserves carbon, removes bromine and its bond, and adds oxygen and its bond.
Unspecified product attributes may remain open for resolution.

Even `[C:1]O>>[C:1]O` does not declare preserved oxygen identity: the unmapped oxygen on each
side denotes removal and addition. Equal descriptions do not imply correspondence. Bonds and
overlays are interpreted through their participants and side presence in the existing span model.

The existing dangling condition remains in force. Deleting an atom must account for its incident
bonds and overlays; an application cannot silently delete additional host context merely because
the source omitted it. No new structural entity or special deletion mechanism is introduced.

Decision: follow [Daylight's pairwise mapping rule](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)
at the SMIRKS boundary. Each atom-map label used on the reactant or product side must occur
exactly once on each side. A one-sided label or multiple occurrences on either side is invalid;
additions and deletions use unmapped atoms. A valid pair declares one preserved entity whose
side attributes are interpreted as described above.

This is a restriction on SMIRKS interpretation, not on ReactionSpan or the existing
reaction-SMILES parser. The latter retains one-sided and repeated mapping classes; those
classes do not acquire an inferred pairing through this decision. No new correspondence
machinery is needed for the accepted pairwise fragment.

### Agent sections

Decision: reject nonempty agent sections for now, consistent with current reaction-SMILES
ingestion into the graph model. The initial supported fragment uses an empty middle section
(`L>>R`). Do not silently discard agent patterns or translate them into unchanged span entities.

Agent support is deferred until a concrete need arises; a later boundary extension can admit
the middle section. No agent-context representation or evaluation machinery is proposed here.
The existing reaction-SMILES boundary parser can retain agents, but graph ingestion rejects
them with AgentsUnsupported; this decision does not change that parser.

### Stereo omission and configuration changes

Decision: interpret the source stereo cases as follows, with configurations compared after
transporting their participant frames rather than comparing @/@@ glyphs directly:

| Stereo specified | Interpretation |
| --- | --- |
| Neither side | Preserve existing host stereo where its frame can be transported. |
| Both sides | Retain or change configuration according to the two side specifications and frames. |
| Lhs only | Erase the configuration specification. |
| Rhs only | Set the specified configuration. |

RDKit documents the corresponding source behaviors in its
[reaction chirality discussion](../materials/codes/rdkit/Docs/Book/RDKit_Book.rst).
For umol, erasure means an undetermined rhs configuration when the stereo entity and its frame
remain applicable. It does not assert that the site is non-stereogenic. Removal of the stereo
site or loss of an applicable participant frame is a separate structural matter, not a synonym
for erasing the configuration value.

Decision: reuse the existing SMILES treatment of stereo context and its translation into umol
participant frames. Incomplete ligand context is not a new pattern-specific design question.
The [OpenSMILES specification](../umol-io/spec/opensmiles-spec.md#double-bond-stereochemistry)
and existing SMILES interpretation own source ordering and incomplete directional information;
the [DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md) owns the ordered ligand frame,
virtual ligands, and configuration relative to that frame. This exploration introduces no
additional blanket rejection rule for incomplete stereo context. The reaction-specific work
is the correspondence and frame transport already described above. Any further dialect-specific
matching clause must first be checked against that existing treatment.

### Nonliteral RHS values and product alternatives

Decision: a product whose charge is {0,+1} and the corresponding pair of concrete products
have equivalent semantics. The same applies to an undetermined RHS value and its admissible
concrete alternatives. These statements retain all other structural facts, context preservation,
and correlations; they do not replace correlated alternatives by independent field sets.

“Assign the partial value” and “require the product to satisfy that value” describe the same
product denotation, not competing semantic choices. Enumerating concrete results is a possible
materialization of that denotation. A wildcard need not determine a single concrete outcome for
the transformation to have a clear meaning. Matching checks membership in the described reaction
relation; application obtains the admissible products for a given input.

The mechanics of symbolic products, refinement, and concrete enumeration are not settled here.
That does not make the set-valued semantics ambiguous. In particular, the representation's ability
to store a set is not a claim that all enumeration or application mechanics are implemented.

### Correspondence with the existing DSL ring constraints

The [DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md), §§5.5 and 7.3, already provides
these observations through #R (total or size-scoped ring count) and #x (incident ring-bond count).
The numbered forms have the following translations; the atom subject is omitted here:

| Pattern primitive | Existing DSL constraints |
| --- | --- |
| `Rn` | `#Rn` (with `#R!` also spelling count zero). |
| `xn` | `#xn`. This counts bonds, unlike #y, which sums their orders. |
| `rn`, n ≥ 3 | `#R(n)+`, together with `#R(k)!` for every ring size 3 ≤ k < n. |

For example, `r6` requires `#R(6)+#R(3)!#R(4)!#R(5)!`. These are existing size-scoped constraints;
no new smallest-ring primitive is necessary. All counts use the same bounded reference, so no
positive size membership exists above 22.

There is a discrepancy in the DSL spec's current SMARTS-parity annotation: it labels `rn` as
`#R(n)+` alone. That expresses membership in a ring of size n, not that n is the smallest size.
The formalization retains the settled smallest-ring reading and records the extra exclusions.
The existing DSL constraint semantics are unchanged; its parity annotation needs a focused
correction when this correspondence is incorporated into the specification.

Omitted numerals and special cases still require surface clauses. The colon decision identifies
its aromatic projection; the full bond-expression clause should also state how it treats
localized order.

## Review of the remaining probes

Check existing umol semantics before treating a probe as a new design question. The following
cases already have answers in this record or its owning discussions:

| Probe | Existing answer |
| --- | --- |
| Mapped RHS fields and wildcards | The decisions above distinguish preservation, explicit changes, and independently resolved H counts. Other field-specific omissions still need their source clauses checked. |
| Added/deleted atoms | Shared mapped identity and side-only entities become the existing span states, under the existing dangling condition. |
| One-sided or repeated atom-map labels | Follow Daylight SMIRKS: each reactant/product map label occurs exactly once on each side; additions and deletions use unmapped atoms. |
| A nonliteral RHS attribute | A symbolic product and its concrete alternatives have the same denotation; materialization mechanics remain separate. |
| Reordered mapped stereo neighbors | Transport the participant frame before comparing configurations; glyph changes alone do not establish inversion. |
| Independent ranges versus one shared variable | Independent ranges permit independent values; reuse of one variable correlates its occurrences, within the scope where its identity is shared (doc 115). |
| A drawn chain versus an existential path | A chain contributes mapped entities; path interiors are existential witnesses and do not enter the correspondence (doc 219). |
| A two-edge path between adjacent triangle vertices | Such a simple path exists, while shortest distance is one; existential path length and shortest distance are different observations (doc 219). |

Named-value scope across nested queries is settled by the
[decision in doc 115](115-variable-facility-2026-06-16.md#nested-query-scope-decision-2026-09-05):
outer bindings are readable, inner variables remain existential and local, and negation exports
no bindings. This answers the nested visibility question left open in docs 193 and 197;
binding syntax and implementation remain deferred. It concerns the variable-facility extension,
not a missing definition of ordinary SMARTS atom-map syntax.

Agent sections are unsupported in the initial fragment, as in current reaction-SMILES ingestion.
Unmatched RHS wildcards use Added entities with open attributes, as specified above.
Stereo context follows the existing SMILES treatment; it is not an unresolved pattern-specific
probe. The
probe inventory is exploratory, not an implementation sequence or a newly written test suite.

For any chosen fragment, semantic preservation means its language interpretation agrees with
its translation into umol, including occurrences and reaction changes, not merely that it parses
or prints again. Relabeling and traversal changes must preserve meaning after frame transport.
Toolkit comparison can expose disagreements but cannot decide the intended semantics by majority
vote. No differential executions or new conformance tests were performed for this exploration.

## Assessment

A formalization is a credible goal because the molecular and reaction semantics already exist.
The work is to make SMARTS/SMIRKS constructs precise against them, and to expose the places where
format fidelity conflicts with that meaning. A useful result could follow the OpenSMILES spec's
structure while being more explicit about asserted projections, Boolean scope, nested-query scope,
and the interpretation of two-sided reaction notation as ReactionSpan, with its existing
conversion to deltas.

The ambitious result is a coherent language semantics, independent of a parser and a toolkit.
The smaller result is an explicitly bounded variant with a faithful boundary translation. The
same construct-by-construct analysis serves both. Existing implementation follow-ups retain their
owners in docs 115, 193, 195, 197, and 219; none is made a prerequisite for exploring the formalization, and
this record tracks no implementation work.

## What the exploration establishes

For the agreed fragment, the existing molecule, constraint, correspondence, ReactionSpan,
and reaction models appear semantically sufficient. No new foundational reaction object has
emerged from the probes. This is evidence for sufficiency, not a completed coverage proof,
a conformance result, or a claim that every evaluator is implemented. Much of the conceptual
work preceded this spike; the boundary analysis tests whether those abstractions accommodate
the motivating legacy languages.

Multiple embeddings and multiple admissible products remain meaningful. The pairwise SMIRKS
restriction excludes ambiguous or one-sided mapping labels at that boundary, not multiple
matches. Non-pairwise map-class semantics also occur in
[Daylight reaction SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html);
they are not an RDKit invention.

Reaction recognition and rule application still need explicit operation contracts. A query
can recognize a specified change while leaving other changes unconstrained. Applying a rule
preserves context it does not change. Likewise, an unmapped atom in a general reaction query
does not necessarily assert deletion, whereas the accepted SMIRKS interpretation does.
These are distinct uses of the semantic objects, not interchangeable interpretations of every
two-sided expression. Formalizing their relationship remains part of the scope.

### Coupling and the original motivation

The difficulty exposed by legacy reaction processing is not attributable solely to bugs,
hydrogen heuristics, or accumulated syntax cases. Molecular representation, chemical perception,
identity tracking, query interpretation, and product construction can jointly determine what
an operation means. Their coupling makes users manage interactions through operational protocols.

The local Colibri and Colibri2 implementations provide concrete evidence. In
`/Users/dr/Dropbox/Source/python/colibri/colibri/task/rdkitlibgen.py`, the exact-mapping
reaction path expands hydrogens, Kekulizes while clearing aromatic flags, applies the reaction,
sanitizes products, disables implicit hydrogens, and removes hydrogens again. The corresponding
Colibri2 generator additionally recovers identity through RDKit's old_mapno and react_atom_idx
properties. These observations are not a finding that every step works around a bug; some
accommodate deliberate representation conventions.

The first four discussion records, inspected in git at f3445b7ed, connect the project to those
frustrations. In particular, discussion/02-design-summary-2025-02-16.md already calls for
multiple explicit molecular models, explicit lossy conversion, and DPO reaction transformations.
The apparent simplicity of the present semantic account is therefore consistent with the
foundations having absorbed the problem that motivated them.

This supports the whitepaper's thesis that separating concepts gives each more freedom to be
manipulated independently. It does not prove that thesis or establish comparative superiority
over a toolkit.

## Questions beyond SMARTS and SMIRKS

The following are exploratory questions, not an implementation plan. Begin each by checking
what existing umol semantics already determines. The opportunity is to describe families of
changes and relationships between changes, with SMARTS/SMIRKS as one input language.

### Reaction specialization and overlap

These are the principal questions for further exploration:

- **Specialization:** when does one reaction description specialize another? A candidate
  denotational reading is inclusion of the concrete reactions described, with correspondence
  and preservation conditions retained. What follows already from attribute refinement and
  constraints, and what requires reasoning about structural occurrences and context?
- **Overlap:** when do two descriptions admit a common concrete reaction? Can that overlap
  itself be represented as a reaction description or a family of descriptions? What evidence
  witnesses overlap or establishes disjointness?
- **Equivalence:** can two differently written rules denote the same reactions? How should
  this differ from merely producing equal products for a particular input, or from equality
  after discarding correspondence information?
- **Composition and net change:** composition already belongs to the model. Can we ask whether
  a sequence of reactions satisfies a net-change constraint, or whether its composite
  specializes a given reaction description, before selecting concrete substrates?
- **Scope of comparison:** do we compare partial reaction observations, context-preserving
  rule applications, or complete concrete reaction instances? Fixing this is necessary before
  treating inclusion or intersection as an algorithmic problem. How are multiple occurrences,
  alternative correspondences, and model-dependent resolution accounted for?

Clear denotations do not by themselves provide decision procedures. Which useful fragments
support effective specialization, overlap, and equivalence checks without general constraint
solving? Which results can remain symbolic? No algorithm, complexity bound, or public API is
settled here.

### Other expressive directions

- **Relations between before/after values:** can a rule state that an atom loses one electron,
  or that two charge changes sum to zero, without enumerating all starting states? The variable
  facility supplies a route to correlated changes; its scope and evaluation restrictions
  remain relevant.
- **Overlay changes:** how can reaction descriptions directly express changes to coordination,
  multicenter bonding, aromatic systems, and stereo frames using their existing entity
  semantics, rather than encoding every change through atom and localized-bond fields?
- **Uncertain correspondence:** how can alternative correspondences express incomplete
  knowledge about an observed reaction? The initial pairwise SMIRKS boundary restriction
  need not limit the eventual language of reaction knowledge. Any representation must retain
  the correlations between alternatives rather than merge them into independent choices.
