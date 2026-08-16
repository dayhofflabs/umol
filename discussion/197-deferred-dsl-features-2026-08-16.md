# 197 — Deferred DSL features

Status: Proposed
Date: 2026-08-16
Relates: [115](115-variable-facility-2026-06-16.md),
[163](163-release-preparation-2026-07-26.md),
[164](164-dsl-edn-worklist-2026-07-27.md),
[195](195-molecule-constraint-matching-2026-08-12.md),
[DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md)

## Purpose

This document collects language and operation ideas deliberately omitted from
the normative DSL specification because their semantics or operational support
are not complete. It is an index of possible future work, not a roadmap, a
compatibility commitment, or a promise that any item will retain its current
internal representation or surface syntax.

The code may already parse, render, store, or partially operate on some forms
listed here. That implementation is experimental unless and until a complete
contract returns to the normative specification with conformance coverage.
Removing a form from the specification does not require removing it from the
code.

## Admission policy

A feature returns to the normative specification only when its public syntax,
semantics, supported operations, failure behavior, and conformance evidence are
settled. A partial implementation may remain undocumented, change shape, or be
removed without first resolving this document.

Detailed designs should continue in focused discussion documents. This file
records the common inventory and links to an existing owner where one exists.

## Deferred inventory

### Variables and expression evaluation

The deferred variable facility includes:

- numeric, element, isotope, and stereo-coset variables;
- variable domains, binding identity, and scope;
- boolean and arithmetic expressions whose evaluation depends on variables;
- reuse across attributes, atoms, molecule patterns, reaction sides, or
  networks;
- match-result bindings and reaction-product substitution.

The existing syntax and AST variants do not settle the eventual design. In
particular, a later design may distinguish anonymous field-local expressions
from named variables and may choose different scoping or surface notation.
Doc 115 contains the existing design prospecting.

### Relative and nonliteral stereo expressions

Coset variables and sets, mirror or involution operators, and permutation
actions are deferred. Their current parser and storage forms do not commit the
eventual matching, binding, normalization, or reaction semantics.

The higher stereo classes axial, square-planar, trigonal-bipyramidal, and
octahedral are likewise outside the normative DSL surface until their complete
matching and model behavior are selected and covered. Tetrahedral atom stereo
and cis/trans bond stereo remain the normative classes.

### Coordination and haptic bonding

The normative dative form is currently one donor directed to one acceptor.
Multi-donor coordination and haptic interactions require a settled entity split
and electron-counting model rather than being promised as a generalized dative
entry. This includes deciding how atom-side donated- and accepted-pair
constraints such as `#d` and `#t` project over the eventual entities.

Existing multi-donor storage is not a commitment to preserve that
representation.

### Overlay-entity constraints in matching

The DSL can represent assertions on dative bonds, aromatic systems,
multicenter bonds, noncovalent bonds, stereo atoms, and stereo bonds. Complete
substructure evaluation against derived host readings remains deferred,
including:

- aromatic- and multicenter-system total-electron assertions;
- dative aromatic and ring assertions;
- noncovalent intramolecular assertions;
- stereo ligand-symmetry, fluxionality, topicity, and stereogenicity
  assertions.

Representation, parsing, rendering, validation, resolution, and matching are
separate operation contracts. Support in one does not imply support in the
others.

### Molecule-scope constraints in matching

Evaluation of entity leaves lifted to molecule scope, relational constraints,
molecule-wide constraints, and boolean combinators during substructure matching
is deferred to doc 195. The current matcher rejects patterns carrying such
constraints rather than weakening them silently.

### Dative ring projection

Deriving dative-bond ring membership requires a ring model whose topology
includes dative overlays. The model and its interaction with localized-ring
queries remain unsettled. A stored dative ring assertion alone does not select
that future model.

### Compound rule evaluation

Higher-order or compound rules whose products are themselves rule-bearing
terms are outside the current reaction contract. The reaction DSL remains an
LHS molecule plus ordered edits; no rule-composition or higher-order evaluation
surface is promised.

## Removed specification text

The following snapshots preserve the specific grammar and operative language
removed from the DSL specification on 2026-08-16. Former-placement notes record
where each passage came from; return-placement notes record where a future,
settled contract would belong. The excerpts are historical design inputs only;
their syntax and semantics are not selected by this document.

### Variable contexts and expression grammar

Former placement: §1, the non-normative term-algebra sketch:

> - **Query / Pattern**: adds wildcards (`*`), element/numeric sets,
>   `bool-expr`, and `?id` references. Sufficient for substructure queries.
> - **Rule**: adds element variables and cross-atom `id` scope (§6).
>   Sufficient for transformation rules.

Former placement: §3, Evaluation context:

```text
| Query | Wildcards and constraints per §5.1 and §7. |
| Rule  | Binds, sets, boolean expressions, and arithmetic per §5.1, §6, and §7. |
```

Former placement: §5.1, `value-expr` and its infix grammar:

```ebnf
value-expr ::= '*'
             | nat-set
             | nat
             | range
             | '?' id
             | '?' id '::' nat-set
             | bool-expr

bool-expr  ::= or-expr
or-expr    ::= and-expr ( '|' and-expr )*
and-expr   ::= not-expr ( '&' not-expr )*
not-expr   ::= '!' not-expr
             | '(' bool-expr ')'
             | rel-expr
rel-expr   ::= mem-expr ( rel-op mem-expr )?
mem-expr   ::= add-expr ( ( '::' | '!:' ) nat-set )?
add-expr   ::= mult-expr ( add-op mult-expr )*
mult-expr  ::= unary-expr ( mult-op unary-expr )*
unary-expr ::= sign* base-expr
base-expr  ::= nat | '?' id | '(' add-expr ')'
rel-op     ::= '<=' | '>=' | '==' | '<' | '>'
add-op     ::= '+' | '-'
mult-op    ::= '*' | '/' | '%'
id         ::= [a-zA-Z][a-zA-Z0-9_]*
```

Former placement: §5.1, top-level variables and same-string reuse:

> A `value-expr` may be only a numeric variable `?id`, optionally with a
> finite in-domain `?id :: nat-set`. A bare `?id` lowers to `ArithExpr::Var`;
> the domain form lowers to a membership predicate.

> For one `atom-string` at match time, numeric `?id` values in `bool-expr` are
> fixed from the matched target atom so that all constraints on that `id` hold
> together. The same `id` may appear in several payloads on one `atom-string`;
> one value satisfies every use. Cross-atom reuse of the same `id` is not fixed
> here.

Former placement: §5.2, `element`:

```ebnf
element ::= element-literal
          | '*'
          | element-set
          | '!' element-literal
          | '!' element-set
          | element-var

element-var    ::= '?' id [ ( '::' | '!:' ) element-domain ]
element-domain ::= element-set | element-literal
```

The removed prose made `element-var` a Query/Rule form, allowed an optional
membership or exclusion domain, and required a bare variable to have an
existing rule-scope binding.

Former placement: §5.3, `isotope-payload`:

```ebnf
isotope-payload ::= '=' | '*' | signed-int | nat-set
                  | '!' signed-int | '!' nat-set
                  | '?' id [ '::' isotope-domain ]
isotope-domain  ::= nat-set | '!' signed-int | '!' nat-set
```

Former placement: §6.3, Bindings:

> For a fixed embedding of the LHS pattern into the target, each `id`
> introduced by an element variable or by `?id` in `bool-expr` / `value-expr`
> has exactly one value in that match. The engine first chooses an embedding,
> then the match binding is fixed: numeric for `?id`, element symbol for
> nominal variables.

> **Multiple results from one ground target.** Ambiguity does **not** require
> an indeterminate target. The **same** ground molecule can admit **several**
> distinct **embeddings** of the **same** LHS (e.g. two equivalent
> substituents). Each embedding yields its own **match binding**. Whether the
> rule **fires once**, **once per embedding**, or **aggregates** products is
> **policy** for the rule evaluator, not fixed by this specification.

> An element variable carries element-symbol values. `?id` in `bool-expr`
> carries numeric bind/use for that attribute. Arithmetic applies only to
> numeric `id` values. Nominal `id` may be reused on the RHS via an element
> variable with the same name.

> On one `atom-string`, the same numeric `id` may appear in multiple predicate
> payloads and all uses are jointly satisfied. Whether `id` may also be shared
> across atom strings on a rule LHS or RHS was not fixed.

Former placement: §6.3, Identifier scope:

> On **one** **atom-string**, the same numeric **`id`** may appear in
> **multiple** predicate payloads; all uses denote **one** value and are
> **jointly** satisfied (**§5.1.2**). Whether **`id`** may also be shared
> **across** atom-strings on a rule LHS (or RHS) is **not** fixed here;
> implementations **SHOULD** document cross-atom **`id`** rules. **Order** of
> **predicates** (**§7.3**) is arbitrary; evaluation **MUST** treat all
> constraints on **`id`** as **simultaneous**, not sequential by textual order.

Former placement: §9.3, Substructure query example:

```clojure
{:atoms [[:C "C#h(?h >= 2)"]
         [:N "N"]]
 :bonds [[:C :N :single]]}
```

The normative example now uses the anonymous range `C#h(2..)`, which expresses
the supported predicate without advertising a named variable.

Return placement: context distinctions in §1 and §3; numeric expressions in
§5.1; element and isotope variables in §5.2 and §5.3; binding and substitution
semantics in §6 and §8; examples only after those contracts agree.

### Nonliteral cosets

Former placement: §5.8, `coset-expr`:

```ebnf
config     ::= '*' | '!' | '+' | nat | coset-expr
coset      ::= '*' | nat | coset-expr
coset-expr ::= '~' coset-expr
             | '\'' coset-expr
             | coset-expr '^' cycles
             | nat
             | '?' id [ '::' coset-set ]
             | coset-set
coset-set  ::= '{' nat (',' nat)* '}'
```

> **Coset operators (reserved).** The **`~`** (involution) and
> **`^`*cycles*** (group action by a permutation in 0-indexed disjoint-cycle
> notation, **§7.11**) operators, and the coset variable / set forms
> (**`?id`**, **`?id :: {…}`**, **`{…}`**), **parse** as **`coset-expr`** and
> **round-trip**, but their **matching** semantics are **staged** —
> relative-stereo binding and non-tetrahedral coset domains land incrementally.
> Only **ground literal cosets** (and the **`*`** / **`!`** / **`+`**
> sentinels) are presently matched; a conforming matcher **MAY** reject an
> operator / variable **`coset-expr`** until the corresponding stage lands.

The normative grammar now contains only sentinels and literal coset indices.

Return placement: surface grammar in §5.8, structured equivalents in §7.12,
and matching or reaction semantics in §6 and §8.

### Higher stereo classes

Former placement: §7.11, `class` and Class:

```ebnf
class ::= 'Th' | 'Ct' | 'Ax' | 'Sp' | 'Tb' | 'Oh'
```

> `Th` tetrahedral, `Ct` cis/trans, `Ax` axial (allene-type), `Sp`
> square-planar, `Tb` trigonal-bipyramidal, `Oh` octahedral. Matching presently
> realizes `Th` and `Ct`; `Ax` / `Sp` / `Tb` / `Oh` parse and round-trip but
> their matching is staged.

Former placement: §7.12, structured constraint grammar:

```ebnf
stereo-kind ::= :tetrahedral | :cis-trans | :axial | :square-planar
              | :trigonal-bipyramidal | :octahedral
coset-form  ::= int | :undetermined | [ int+ ] | "coset-string"
```

Return placement: class syntax and semantics in §7.11, structured kinds in
§7.12, and grounding and matching behavior in §6.

### Multi-donor dative entries

Former placement: §4, `dative-bond-entry`:

```ebnf
dative-bond-entry ::= { [:id keyword] :donors [ atom-ref+ ]
                        :acceptor atom-ref :attrs dative-bond-spec }
```

> A dative bond entry binds a single acceptor to one or more donors (a
> coordination center): `:acceptor` names the atom accepting the electron
> pair(s); `:donors` is a vector of one or more donating atoms. The leading
> `order` token records the number of donated pairs; one shared `:attrs` covers
> every donor-to-acceptor edge of the entry.

Former placements: §7.12 structural refs, §8 additions/removals, and §8.1 span
modifications used the corresponding plural productions:

```ebnf
dative-bond-ref      ::= int | keyword | { :donors [atom-ref+] :acceptor atom-ref }
dative-bond-addition ::= { :donors [ atom-handle* ] :acceptor atom-handle ... }
dative-bond-removal  ::= { :id dative-bond-handle :donors [ atom-handle* ] ... }
dative-bond-modify   ::= { ... :donors [ atom-ref+ ] :acceptor atom-ref ... }
```

These productions are narrowed to a one-element `:donors` vector in the
normative specification. The field remains plural to match the existing data
shape without promising multi-donor semantics.

Return placement: entity shape and validity in §4 and §4.1, constraint refs in
§7.12, and operational and span forms in §8 and §8.1. A future coordination or
haptic entity split may instead require new sections rather than restoring
these productions verbatim.

### Matching promises for deferred constraints

Former placement: §6.1, Derived predicates:

> Every predicate admitted in the DSL that is not an inherent field is a
> derived predicate: a topological query evaluated against the target graph
> once an embedding is proposed. This includes per-aromatic, per-multicenter,
> per-dative ring-membership predicates and the molecule-wide entries of
> §7.12. Derived predicates filter matches.

Former placement: §6.1, Symmetry-derived stereo predicates:

> The stereo entity predicates `#p` ligand symmetry, `#o` topicity, and `#g`
> stereogenicity are derived from the resolved molecule's graph automorphisms.
> As derived predicates they filter matches. `#f` fluxionality is matched as a
> stored predicate.

Former placement: §6.2, Molecule-level match:

> A molecule-map pattern matches a target iff every atom string and bond string
> matches field-wise, each structural relation matches per its own inherent
> fields, and every derived predicate holds on the resulting embedding.

The replacement contract promises derived evaluation only for atom and
localized-bond constraints, preserves inherent-field and participant matching
for overlays, and rejects a non-empty molecule-scope constraints list.

Return placement: supported derived predicates in §6.1 and the complete
molecule-level match contract in §6.2, after each constraint family has
conformance coverage.

### Dative ring projection

Former placement: §5.5 and §7.7:

> Dative-bond `#R` remains an asserted `RingMembershipForm` value: deriving it
> requires a ring model that includes dative overlays, whose semantics are not
> defined by this specification.

The §7.7 table described its storage as "asserted constraint; topology
derivation deferred." The normative text now defines only the stored assertion.

Return placement: the selected ring projection in §5.5, the dative predicate
meaning in §7.7, and its matching behavior in §6.1.

### Compound rule evaluation

Former placement: §1, EDN and rules:

> **EDN and rules.** EDN carries the relational structure. **Rule** evaluation
> (pattern **LHS** → product **RHS**, **§6**) is a separate computation layer:
> it **MAY** consume and produce molecule maps that use the same surface
> notation. The **reaction map** (**§8**) is the operational encoding of such a
> rule — the **LHS** plus an ordered **`:deltas`** edit list whose application
> yields the **RHS**.

Former placement: §1, the non-normative term-algebra sketch:

> - **Compound**: rules whose RHS produces molecule maps that are themselves
>   other terms; higher-order composition.

Return placement: the term-algebra overview in §1 and an operational contract
adjacent to §8, if higher-order rule composition is selected.

## Open questions

- Whether related items should return together or as smaller independently
  conforming subsets.
- Which existing parser extensions are useful enough to retain while their
  semantics remain experimental.
- Whether a deferred item should receive a focused design document before work
  begins or can be removed from this inventory without replacement.
