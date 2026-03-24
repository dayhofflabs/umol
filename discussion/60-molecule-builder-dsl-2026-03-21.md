# Molecule DSL

## Preface

The DSL rests on three connected ideas.

**Homoiconicity.** Data and patterns share notation. A ground molecule and a substructure
query differ only in whether their attribute slots contain concrete values or constraint
expressions. An atom spec `Cv4H2` and an atom query `C?v4H*` share the same grammar;
the spec is a degenerate query that matches exactly one atom type. This collapses the
usual asymmetry between data formats and query languages, and is the reason that a single
representation covers both `Molecule`/`MoleculeBuilder` objects and query objects.

**Relational structure.** A molecule is a set of named relations — an atom relation, a
covalent-bond relation, a dative-bond relation, an aromatic-system relation, and so on —
matching `MoleculeBuilder`'s fields directly. Each relation is a collection of typed
tuples. Fragment composition is relation union. Substructure matching is relational
containment.

**EDN + Datalog.** EDN (Extensible Data Notation) represents the relational structure
natively: a molecule is an EDN map whose keys are relation names and whose values are
vectors of tuples. The semantic operations — matching, transformation, resolution — are
instances of Datalog-with-arithmetic evaluation: rules of the form
`LHS-pattern → RHS-pattern` with bound variables and arithmetic guards. EDN is the data
layer; the rule evaluator is the computation layer. No JVM is required: a Rust
EDN/Clojure reader (`clojure-reader`) can serve as the host reader while the Datalog
rule evaluator remains part of umol itself.


## Term algebra

Four levels, defined by mathematical structure.

### Level 1 — Ground terms

Fully resolved molecules. All attribute slots are concrete values. Extension = singleton.
Values of type `Molecule`.

```edn
{:atoms   {:c1 #atom "Cv4H3"  :o1 #atom "Ov2H1"}
 :bonds   [[:c1 :o1 :single]]
 :charge  0
 :spin    #spin "^0x1"}
```

### Level 2 — Constraint terms

Predicates over the countably infinite discrete set of all chemical graphs. Some attribute
slots contain wildcards, variable bindings, or guards. Extension = arbitrary (possibly
infinite) subset of chemical space.

```edn
{:atoms   {:c1 #atom "C?v4H*"  :o1 #atom "O?v*H*"}
 :bonds   [[:c1 :o1 :single]]
 :charge  nil
 :spin    nil}
```

**Partial structures and queries are the same object.** Both are constraint terms. The
distinction is purely in the evaluation context:

- *Builder context*: the term is expected to narrow to a singleton under resolution.
  `MoleculeBuilder` fails if the constraint is underdetermined after all resolution passes.
- *Query context*: the term selects an arbitrary-sized subset. The same constraint term
  that fails as a builder spec may be a valid substructure query.

Within level-2 semantics, the atomic properties provided in the #atom token are treated
as exact match predicates, wildcards define constraints. Moreover, variables, expressions,
and guards are allowed.

### Level 3 — Rules

Pairs `{:lhs <constraint-term> :rhs <constraint-term>}`. The LHS is matched against a
molecule (or multiset of molecules for n-ary rules); the RHS is constructed from the
match, with LHS-bound variables in scope.

```
unimolecular:  Molecule → Option<Set<Molecule>>
n-ary:         MultiSet<Molecule> → Option<Set<MultiSet<Molecule>>>
```

`None` = LHS did not match. The output is a set because a match may produce multiple
valid products (e.g. multiple Kekulé structures, multiple tautomers).

**Manipulations and reactions are the same type.** Kekulization, tautomerization,
aromatization, and named reactions all have this signature. The distinction is purely
about LHS specificity:

| LHS | Example |
|---|---|
| Structural motif | Kekulization (any aromatic system) |
| Reaction template | Alcohol oxidation (any primary or secondary alcohol) |
| Fully determined molecule | A specific named substrate |

Variables bound on the LHS may appear in arithmetic expressions on the RHS (see
[Grammar — arithmetic in rules](#arithmetic-in-rules)).

For n-ary rules, `:lhs` is a vector of constraint terms; atom labels are namespaced to
distinguish reactant molecules (`:mol1/c1`, `:mol2/c3`).

### Level 4 — Compound terms

EDN lists representing function application. Arguments may be any level-1–3 terms or
further compound terms.

```clojure
(apply-rule oxidize-alcohol methanol)
(mix (enantiomer mol) (invert-all (enantiomer mol)))
(fuse benzene pyrrole {:benzene/c3 :pyrrole/c0})
(select-tautomer lactam criterion/lowest-energy)
```

Return types are ground terms, sets of ground terms, or further compound terms. Level-4
terms are EDN lists, not molecule maps; they are computation above the data layer and do
not alter the molecule map schema.


## Grammar

### Molecule map schema

```
molecule-map ::=
  { :atoms    atom-collection
    :bonds    bond-list
    [:dative   bond-list]
    [:aromatic aromatic-list]
    [:mc       mc-list]
    [:nc       bond-list]
    [:charge   int | nil]
    [:spin     spin-literal | nil]
    [:expect   { :charge int :multiplicity keyword }]   ; test checksum only
    [:guards   [ logic-expr* ]]                          ; rule context only
  }

atom-collection ::= { keyword → atom-literal }+    ; named form  (authoring surface)
                  | [ atom-literal* ]               ; indexed form (canonical)

bond-list    ::= [ bond-entry* ]
bond-entry   ::= [ keyword keyword bond-spec ]     ; full form
               | [ keyword keyword bond-keyword ]  ; shorthand (see below)

aromatic-list ::= [ { :atoms    [keyword+]
                      :electrons nat
                      [:charge  int]
                      [:spin    spin-literal] }* ]

mc-list ::= [ { :atoms    [keyword+]
                :electrons (nat | nil)
                [:charge  int] }* ]
```

Absent optional keys = "not applicable" (e.g. no aromatic systems in this molecule).
`nil` values = "unknown" in builder context, "unconstrained" in query context. The
distinction matters: absent `:aromatic` means no aromatic systems; `:aromatic nil` would
mean the aromatic systems are not yet known.

**Named vs indexed atoms.** The named form (`:atoms` as a map, keys are EDN keywords) is
the authoring surface. The indexed form (`:atoms` as a vector) is the serialized canonical
form. Round-trip: named → indexed by assigning positions in insertion order; indexed →
named by restoring stored labels from `MoleculeBuilder.labels` or generating synthetic
labels (`:a0`, `:a1`, …). Fragment-scoped names use EDN keyword namespaces:
`:benzene/c1`, `:pyrrole/n1`. Map `merge` on `:atoms` plus `concat` on all other sections
is well-defined fragment composition with no index arithmetic.

**`:expect`** is a test-only checksum: `{:charge 0 :multiplicity :singlet}`. Validated at
parse time; not stored in the model. See [Design Notes — charge and spin](#charge-and-spin).

**`:guards`** appears only in rules (Level 3). It holds predicate expressions over
variables bound in the LHS atom specs that cannot be expressed inline (e.g. a guard
spanning two atom attributes). Guards inline in atom specs (e.g. `H?h>=1`) are
preferred for locality.

### Atom spec

The atom spec is a compact positional token string encoding all per-atom attributes.
Depending on the context, the spec is submitted to the atom or the query parser.

```
atom-spec    ::= element spec*
spec         ::= charge | valence | H-spec | lone-pairs | spin-inline | av | mv | ap | dp

element      ::= [A-Z][a-z]*
charge       ::= [+-] value-expr
valence      ::= 'v' value-expr
H-spec       ::= 'H' value-expr
lone-pairs   ::= '/' value-expr
spin-inline  ::= '^' value-expr 'x' value-expr   ; (unpaired, 2S+1) — shared with #spin
av           ::= 'a' value-expr                   ; aromatic valence
mv           ::= 'm' value-expr                   ; multicenter valence
ap           ::= '<' value-expr                   ; accepted pairs
dp           ::= '>' value-expr                   ; donated pairs


value-expr   ::= nat                              ; L1, L2, L3: concrete value
               | '*'                              ; L2, L3:     wildcard (unconstrained)
               | '?' id                           ; L3 LHS:     bind variable
               | '?' id cmp nat                   ; L3 LHS:     bind with guard
               | '(' id op nat ')'                ; L3 RHS:     arithmetic over bound var

cmp          ::= '>=' | '<=' | '>' | '<' | '='
op           ::= '+' | '-' | '*' | '/'
```

Tag:
- `#atom "Cv4H2"` — ground spec; accepted in ground and constraint contexts
- `#atom "C?v4H*"` — constraint spec; accepted in constraint/rule contexts

The binding forms (`?id`, `?id>=n`) and arithmetic expressions (`(?h-1)`) are only
well-formed inside a rule (Level 3). They are parse errors in standalone molecule maps.

Examples:

| Payload | Meaning |
|---|---|
| `Cv4H4` | Carbon, valence 4, 4 implicit H (methane) |
| `Nv3H2/1` | Nitrogen, valence 3, 2H, 1 lone pair (amine) |
| `C?v4H*` | Carbon, valence 4 required, H count unconstrained |
| `C?v4H?h>=1` | Carbon, valence 4, bind H count to `h`, guard h ≥ 1 |
| `C?v4H(?h-1)` | Carbon, valence 4, H count = h−1 (rule RHS) |

### Bond spec

The bond spec encodes bond *properties*. It does not encode atom indices — those are
structural and live in the `:bonds` section of the molecule map. The `#bond` tag
dispatches to the bond spec parser; no outer delimiters in the payload.

```
bond-spec     ::= order charge? spin-inline? aromatic-hint?

order         ::= 'o' value-expr
charge        ::= [+-] value-expr
spin-inline   ::= '^' value-expr 'x' value-expr
aromatic-hint ::= 'a' ('0' | '1')

value-expr    ::= (same grammar as atom spec value-expr)
```

The canonical bond entry is `[:atom-a :atom-b #bond "..."]`. Shorthands:

| Shorthand | Expands to |
|---|---|
| `[:a :b :single]` | `[:a :b #bond "o1"]` |
| `[:a :b :double]` | `[:a :b #bond "o2"]` |
| `[:a :b :triple]` | `[:a :b #bond "o3"]` |
| `[:a :b :aromatic]` | `[:a :b #bond "o1a1"]` |

Full `#bond` form is required for charged bonds, radical bonds, or any bond where the
shorthand is ambiguous. Examples: `#bond "o1+1^1x2"` (radical cation single bond),
`#bond "o2^2x3"` (triplet double bond).

**Aromatic hint lifecycle.** The `a1` flag in the bond spec is a SMILES/MOL parse
artifact. It marks a bond that originated from an aromatic bond token (lowercase atom
or `:` bond type in SMILES). It lives in `BondBuilder.aromatic_hint` during
construction; `Bond` has no aromatic flag.

- *During parsing*: SMILES parser emits `a1` bonds for aromatic-context bonds.
- *After aromatic perception*: perception consumes `a1`, constructs the `:aromatic`
  section, and sets all formerly-`a1` bonds to `:single`. The hint is gone.
- *Canonical ground term*: no `a1` present. Bonds between aromatic atoms are `:single`;
  the π system is encoded in `:aromatic` exclusively. Canonical serializers (SMILES→DSL,
  SDF→DSL, `MoleculeBuilder`→DSL) never emit `a1`.
- *After kekulization*: `:aromatic` section is removed; bonds are promoted to
  alternating single/double. `a1` is not involved — perception already ran.
- *Mixed Kekulé/aromatic* (e.g. biphenyl with one ring aromatic, one Kekulé): the two
  representations coexist natively. One ring appears in `:aromatic` with `:single`
  bonds; the other has explicit double bonds in `:bonds`. No `a1`, no implicit
  sanitization.

The presence of `a1` in a DSL document signals intermediate parser state — aromatic
perception has not yet run. This is a valid transient representation but not a
canonical ground term.

Bond pattern variables follow the same grammar as atom spec value expressions:
`#bond "o?k"` binds bond order to `k`; `#bond "o?k>=2"` adds a guard; `#bond "o?k"`
on the rule RHS preserves the matched order.

### Spin spec

```
spin-spec ::= '^' value-expr 'x' value-expr
```

The spin sub-grammar is shared across all three spec types. The `#spin` tag dispatches
to the spin spec parser for standalone annotations (aromatic system spin, molecular spin
override). The same `^nxm` token pair appears inline within atom specs and bond specs.

Tag: `#spin "^2x3"` — triplet, 2 unpaired electrons.

The three parsers (`#atom`, `#bond`, `#spin`) are independent but share the
spin sub-grammar. Atom and bond specs embed spin inline via `^`/`x` tokens; `#spin`
provides standalone spin literals at the molecule map level.

### Arithmetic in rules

SMARTS/SMIRKS handles atom mapping but cannot express arithmetic on attributes. A single
rule for alcohol oxidation must handle primary (`CH₂OH → CHO`, H: 2→1) and secondary
(`CHOH → CO`, H: 1→0) cases with one template.

```edn
{:lhs {:atoms {:ca #atom "C?v4H?h>=1"  :o #atom "O?v2H1"}
       :bonds [[:ca :o :single]]}
 :rhs {:atoms {:ca #atom "C?v4H(?h-1)" :o #atom "O?v2H0"}
       :bonds [[:ca :o :double]]}}
```

- `?h` binds the matched H-count on the LHS.
- `?h>=1` is an inline guard (prevents H going negative).
- `(?h-1)` is a computed attribute on the RHS.
- The same rule applies to both cases because `?h` is universally quantified.

The grammar extension is additive: `value-expr` in the atom spec string gains binding and
arithmetic forms that are parse errors outside rule context.

### Document identity

Every molecule map carries two namespaced keys identifying the model kind and schema version:

```edn
{:umol/kind    :molecular-graph
 :umol/version 1
 :atoms {...}
 :bonds  [...]}
```

`:umol/kind :molecular-graph` distinguishes graph DSL documents from future geometric DSL
documents (`umol-models-geometric`). The two models are not in a hierarchy and share no
canonical conversion — a constitutional molecule is a graph-model object; a 3D
conformation is a geometric-model object. The tag makes the distinction explicit at the
document boundary. `:umol/version` enables schema evolution. Both are regular map keys,
not EDN metadata (`^{...}`), for portability across all EDN readers.

### Topology notation (sugar over named atoms)

Topology notation separates graph structure from vertex/edge coloring. It compiles to a
flat named-atom map before `MoleculeBuilder` ingestion — it is not a separate data model.

Three additional keys, all optional and all expanded away before processing:

```
:types    { keyword → atom-literal }       ; atom type aliases
:topology [ [nat nat]* ]                   ; abstract adjacency list
:nodes    [ keyword* ]                     ; positional vertex coloring (index → type)
:names    { nat → keyword }                ; promote positions to named anchors
```

Expansion rules:
1. Expand `:types` aliases everywhere they appear.
2. For each position `i` in `:nodes`: add atom `(get :names i | :node-i)` with the
   specified type to `:atoms`.
3. For each edge `[i j]` in `:topology`: add bond between the corresponding atom labels
   to `:bonds` (`:single` by default; override via `:edge-types {[i j] bond-spec}`).
4. Merge any explicit `:atoms` and `:bonds` entries.

`:names` is the bridge between the positional and named namespaces. Promote only atoms
that need to be cross-referenced (in `:bonds`, `:aromatic`, queries, or rules).
Everything else stays as `:node-N`.

```edn
; 2-chlorobutane — only the substituted carbon needs a name
{:umol/kind    :molecular-graph
 :umol/version 1
 :types    {:ct   #atom "Cv1H3"
            :cmcl #atom "Cv2H1"
            :cm   #atom "Cv2H2"
            :cl   #atom "Clv1H0"}
 :topology [[0 1] [1 2] [2 3] [1 4]]
 :nodes    [:ct :cmcl :cm :ct :cl]
 :names    {1 :c-alpha}}

; Expands to
{:umol/kind    :molecular-graph
 :umol/version 1
 :atoms {:node-0  #atom "Cv1H3"
         :c-alpha #atom "Cv2H1"
         :node-2  #atom "Cv2H2"
         :node-3  #atom "Cv1H3"
         :node-4  #atom "Clv1H0"}
 :bonds [[:node-0  :c-alpha :single]
         [:c-alpha :node-2  :single]
         [:node-2  :node-3  :single]
         [:c-alpha :node-4  :single]]}
```

Topology notation is most useful for congeneric series (same connectivity, different
decoration) and symmetric structures where the type table collapses repetition. It is
**unsuitable for rule patterns** — see Atom mapping below. Rules use named atoms directly.

The `:path` and `:ring` shorthands cover the common cases of linear chains and rings:

```edn
; hexane — open chain with named termini
{:path [:c0 #atom "Cv1H3"  #chain ["Cv2H2" 4]  :c1 #atom "Cv1H3"]}

; cyclohexane — ring with one named anchor
{:ring [:c0 #atom "Cv2H2"  #chain ["Cv2H2" 5]]}
```

A `:path` or `:ring` vector interleaves named-anchor / spec pairs with `#chain [spec n]`
fixed-length segments. Named anchors are referenceable across all sections; `#chain`
segments expand to `:node-N` atoms with sequential single bonds. `:ring` additionally
emits a closing bond from the last atom back to the first (`:single` by default;
override with `#chain [spec n bond-spec]`). Both keys compile to `:atoms` + `:bonds`
entries by the same expansion rules above.

The naming is symmetric with the query notation: `:path`/`:ring` are ground term sugar
(fixed length, no Kleene star); `:paths`/`:rings` are the query-context counterparts
(variable length, Kleene star permitted).

### Atom mapping in rules

Named atoms give implicit atom mapping for free. A label appearing in both LHS and RHS
names the same atom — it persists through the transformation with its attributes updated
by the RHS spec and any computed expressions. Labels present only in the LHS are
consumed; labels present only in the RHS are created.

Bond persistence is implicit via label pairs: `[:ca :o ...]` in both LHS and RHS means
the bond persists, with order/charge updated to the RHS spec.

**Environment atoms** — atoms in the matched molecule not in the LHS pattern — are
preserved implicitly. The rule specifies only what changes; environment bonds reconnect
automatically during graph rewriting.

**Variable scoping.** Each label binds independently, so two LHS atoms matching the same
type carry independent variable bindings:

```edn
; Diol oxidation — :h1 and :h2 do not interfere
{:lhs {:atoms {:c1 #atom "C?v*H?h1>=1" :o1 #atom "O?v2H1"
               :c2 #atom "C?v*H?h2>=1" :o2 #atom "O?v2H1"}
       :bonds [[:c1 :o1 :single] [:c2 :o2 :single] [:c1 :c2 :single]]}
 :rhs {:atoms {:c1 #atom "C?v*H(?h1-1)" :o1 #atom "O?v2H0"
               :c2 #atom "C?v*H(?h2-1)" :o2 #atom "O?v2H0"}
       :bonds [[:c1 :o1 :double] [:c2 :o2 :double] [:c1 :c2 :single]]}}
```

**Bimolecular rules.** `:lhs` is a vector of constraint terms; labels are scoped to their
containing pattern. Mapping is established by which labels appear in the RHS. Atoms can
migrate between product molecules:

```edn
; SN2: amine + alkyl chloride → secondary amine + HCl
{:lhs [{:atoms {:n #atom "N?v3H2"}}
       {:atoms {:c #atom "C?v4H?m" :lg #atom "Clv1"}
        :bonds [[:c :lg :single]]}]
 :rhs [{:atoms {:n #atom "N?v4H1" :c #atom "C?v4H?m"}
        :bonds [[:n :c :single]]}
       {:atoms {:lg #atom "Clv1H1"}}]}
```

`:n`, `:c`, `:lg` are mapped; `:lg` migrates to the second product. No additional mapping
syntax is required.

**Multiple matches.** When the LHS matches at multiple sites, the application strategy
is a parameter of the Level 4 call, not of the rule itself:

```clojure
(apply-rule r mol)        ; apply to one match (deterministic by graph ordering)
(apply-rule-all r mol)    ; apply to each match independently → set of products
(apply-rule-once r mol)   ; apply to all non-conflicting matches simultaneously
```

**Constraint on topology notation in rules.** Positional mapping (`node-0` in LHS →
`node-0` in RHS) is ambiguous when topology changes between LHS and RHS (bonds added,
broken, reordered). Rules must use named atoms. If the rule's structure is regular enough
to tempt topology notation, the atoms at chemically active sites still need `:names`
promotion — at which point topology sugar adds no value over writing named atoms directly.

### Molecule map examples

**Methanol (ground term):**

```edn
{:umol/kind    :molecular-graph
 :umol/version 1
 :atoms   {:c #atom "Cv4H3"  :o #atom "Ov2H1"}
 :bonds   [[:c :o :single]]
 :charge  0
 :spin    #spin "^0x1"}
```

**Indole (ground term, named atoms):**

```edn
{:umol/kind    :molecular-graph
 :umol/version 1
 :atoms   {:n   #atom "Nv2a2H1"
           :c2  #atom "Cv2a1H1"  :c3  #atom "Cv2a1H1"
           :c3a #atom "Cv3a2"
           :c4  #atom "Cv2a1H1"  :c5  #atom "Cv2a1H1"
           :c6  #atom "Cv2a1H1"  :c7  #atom "Cv2a1H1"
           :c7a #atom "Cv3a2"}
 :bonds   [[:n :c2 :single] [:c2 :c3 :single] [:c3 :c3a :single]
           [:c3a :c4 :single] [:c4 :c5 :single] [:c5 :c6 :single]
           [:c6 :c7 :single] [:c7 :c7a :single] [:c7a :n :single]
           [:c3a :c7a :single]]
 :aromatic [{:atoms [:n :c2 :c3 :c3a :c4 :c5 :c6 :c7 :c7a]
             :electrons 10}]
 :charge  0
 :spin    #spin "^0x1"}
```

**Substructure query (primary or secondary alcohol carbon):**

```edn
{:atoms {:ca #atom "C?v4H?h>=1"  :o #atom "O?v2H1"}
 :bonds [[:ca :o :single]]}
```


## Design notes

### Charge and spin

Molecular charge is not a top-level stored field. It is computed:

```
charge = Σ atom.charge + Σ bond.charge + Σ aromatic.charge + Σ mc_set.charge
```

Motivating cases: `HO₃⁺` (charge on atom), `(Br₂)⁺` (charge on bond, `Bond.charge = 1`),
`Cp⁻` (charge on aromatic system, `AromaticSystem.charge = −1`), `(I₃)⁻` (charge on
multicenter set). All are already representable in the struct fields. `Molecule::charge()`
is arithmetic over existing fields. The `:expect {:charge n}` DSL key is a test checksum,
not a primary field; it imposes no storage obligation.

Molecular spin is not a sum — angular momentum coupling requires CG coefficients.
`SpinState::is_constructible_from` implements sequential coupling. Molecular spin is
either provided explicitly (`:spin` key on the molecule map) or validated against the
space of achievable couplings from constituent spins. It is never silently inferred. See
`discussion/61-spin-state-builder-2026-03-22.md` for the uniform `SpinStateBuilder`
design across atom, bond, aromatic, and multicenter features.

Charge and spin annotations in the DSL belong on the features that carry them, not on a
synthetic top-level field. The top-level `:charge` and `:spin` keys exist only as
convenience summaries / override points.

### Implementation patterns

The DSL and resolution pipeline map to named patterns from Fowler's *Domain Specific
Languages* and Parr's *Language Implementation Patterns*:

- **Semantic Model** (Fowler): `Molecule`/`MoleculeBuilder` are the semantic model.
  Parsing the DSL produces a `MoleculeBuilder`, not an interpreted AST. The Rust structs
  are the only authoritative representation.

- **Production Rule System** (Fowler): rules `{:lhs ... :rhs ...}` are production rules.
  One engine covers kekulization, tautomers, and named reactions; there is no separate
  "manipulation" infrastructure.

- **Attribute Grammar** (Parr): resolution phases are attribute grammar computations.
  Valence candidate filtering is synthesized (bottom-up from bond connectivity);
  aromaticity assignment may be inherited (top-down from ring detection).

- **Symbol Table / Scoping** (Parr): named atoms require a symbol table (keyword → index).
  EDN keyword namespaces provide lexical scoping for fragment composition; `fuse` merges
  two scopes with an explicit interface.

- **Computed Attributes**: `H(?h-1)` on a rule RHS is a computed attribute in the
  attribute grammar. The implementation threads a binding environment through LHS
  matching and evaluates arithmetic expressions during RHS construction — standard in
  attributed Datalog systems such as Soufflé.


## Open questions

No open questions currently.


## Resolved questions

**O1 — `#atom` tag unification and type-directed lowering.** Use `#atom` as the single
tag for atom literals in all contexts. Parsing is split into two phases:
(1) parse to a shared `AtomExpr` intermediate representation (same surface grammar in all
contexts), then (2) lower to target type based on context/purpose:

- ground molecule context -> `AtomTypeSpec`
- query/rule context -> `AtomTypeQuery`
- conformance input context -> `TableIR::Atom`
- builder context -> `AtomBuilder`
- resolved molecule (`GraphIR::Atom`) is never parsed directly; it is produced by resolution

This keeps the domain types separate while unifying syntax and frontend parsing. The same
pattern applies to bonds and aromatic systems: parse once to neutral expression types,
then lower to context-specific targets (`TableBond`/bond constraints; aromatic query terms /
builder aromatic systems).

**O2 — Lone-pair token.** Canonical token is `/` (`lone-pairs ::= '/' value-expr`).

**O3 — Donated/accepted pairs.** Keep both in atom grammar:
`dp ::= '>' value-expr`, `ap ::= '<' value-expr`.

**O4 — Label storage location.** Do not add `MoleculeBuilder.labels` now. Keep labels in
DSL parse/serialization layer (symbol table); core builder remains label-agnostic.

**O5 — Aromatic system decomposition.** Keep per-atom aromatic valence on atoms (`av`) and
system membership/metadata in `:aromatic`. Keep `:electrons` authorable but validate against
member contributions and charge. `rings` stays derived from topology/perception.

**O6 — Dative directionality.** In `:dative` entries, tuple order is semantic:
`[:donor :acceptor ...]`.

**O7 — `#chain` grammar.** Make grammar explicit and prefer explicit nested atom literals:
`#chain [#atom "Cv2H2" 4 bond-keyword?]` with default `:single`.

**O8 — Sugar keys in schema.** Include `:types`, `:topology`, `:nodes`, `:names`, `:path`,
`:ring`, and `:edge-types` in formal schema as preprocessing sugar expanded before builder
ingestion.

**O9 — `:aromatic` bond shorthand.** Keep `[:a :b :aromatic] -> [:a :b #bond "o1a1"]` as
input-only transient sugar; never emit in canonical ground-term serialization.

**O10 — Top-level `:charge` in canonical terms.** `:charge` remains optional summary/guard
field. If present, it must equal the charge computed from feature-local charges; it is not
authoritative over atom/bond/system charges.

**O11 — Noncovalent structure.** `:nc` uses a dedicated entry structure (interaction type +
optional constraints) rather than reusing covalent `bond-spec`.

**O12 — Charge type range.** Moot in this document (molecular charge type updated in code).

**O13 — Host reader and parser responsibility split.** The host reader should be
`clojure-reader` (replacing earlier `edn-rs` mention). EDN/Clojure syntax parsing stays in
the host reader crate; umol implements domain semantics (`#atom`, `#bond`, `#spin`,
context checks, lowering) in its own code. No custom EDN parser should be written.

**O14 — Parse mode naming.** Mode name is `Conformance` (not `ConformanceAtom`) to keep
scope open for conformance-level parsing of atoms, bonds, and aromatic-system sections.

**Q1 — Repeat notation.** Superseded by topology notation and Level 4 composition.
Topology notation (`:topology` + `:nodes` + `:types`) handles congeneric series and
symmetric structures more expressively than a raw repeat count. The `:path` /
`#chain` shorthand handles linear chains with named anchors. Level 4 operations
(`(chain n type)`, `(ring n type)`) handle programmatic construction. No repeat syntax
is needed in the static DSL.

**Q4 — N-ary rule atom scoping.** Labels in a bimolecular rule are scoped to their
containing pattern (each element of the `:lhs` vector). No namespace prefixes are
required — labels are unambiguous because they are local to their pattern. EDN keyword
namespaces (`:mol1/c1`) are available for documentation clarity but carry no semantic
weight. The mapping is determined solely by which labels appear in the `:rhs`.

**Q6 — Canonical indexed form.** Canonical form is insertion order: hand-authored
named-atom maps serialize in the order atoms were written; Level 4-produced maps use
construction order. Graph-isomorphism canonicalization is deferred to query-time
deduplication (a separate concern from authoring). Two hand-authored representations of
the same molecule may serialize differently — this is acceptable because the authoring
surface is named atoms, not integer indices, and identity is checked semantically by the
rule engine, not textually.

**Q2 — Bond pattern variables.** Confirmed. Bond specs use identical `value-expr` grammar
to atom specs: `#bond "o?k"` binds bond order; `#bond "o?k>=2"` adds an inline guard;
`#bond "o?k"` on a rule RHS passes the matched order through unchanged. Maximum
consistency between atom and bond query grammar; implementation work only.

**Q3 — Feature matching semantics.** Resolved by a uniform compositional principle:
each feature type is matched by checking predicates at its own level plus the predicates
of its constituent lower-level features.

- *Atom*: match if all atom-spec predicates hold on the target atom.
- *Bond*: match if atom predicates hold on both endpoint atoms AND bond predicates hold
  on the bond itself.
- *Aromatic system*: match if atom predicates hold on each member atom, bond predicates
  hold on each bond internal to the system, and system-level predicates (`:electrons`,
  `:charge`, `:spin`) hold on the system as a whole. A query specifying `k` atoms matches
  any aromatic system containing at least those `k` atoms with those properties — the
  query system need not specify all atoms in the target system (subgraph semantics within
  the aromatic relation).
- *Multicenter set / noncovalent*: same principle — member atom predicates + system
  predicates.

This is subgraph matching applied uniformly across all five `MoleculeBuilder` collections.
The aromatic `av` token on individual atom specs is an atom-level predicate (this atom
participates in *some* aromatic system); the `:aromatic` section in the query specifies
system-level predicates. The two are orthogonal and can be combined.

**Q5 — `:expect` as guards.** `:expect` is syntactic sugar over guards on computed
molecule properties. `{:charge 0 :multiplicity :singlet}` is shorthand for the guard
expressions `(= molecule/charge 0)` and `(= molecule/multiplicity :singlet)`, where
`molecule/charge` is the sum-over-features charge and `molecule/multiplicity` is the
resolved spin multiplicity. Any guard expression that can be written in the `:guards`
list can be used in `:expect`. The `{:charge n :multiplicity k}` map form is retained as
readable shorthand for the common case. Error behavior: parse-time assertion failure —
`:expect` is a test checksum, not a structural field, and violation is always an error
at molecule-map ingestion time regardless of context.

**Q7 — DSL as canonical parser output.** Confirmed. The DSL ground-term schema is the
canonical serialization for resolved `Molecule` / `GraphIR` objects. All umol parsers
(SMILES, SDF, CIF, and any future format) emit DSL molecule maps as their primary output.
SMILES → DSL round-trip is a correctness criterion: the round-trip must produce a
semantically equivalent molecule map. `QueryMolecule` and `Reaction` as DSL objects are
planned but the schema details are deferred until those types are designed.

**Q8 — `:path`/`:ring`: `:path` (open chain) and `:ring` (closed cycle) are syntactic
sugar for homogenous molecule fragments. Syntax is identical — interleaved named-anchor/spec
pairs with `#chain [spec n]` fixed-length segments — with one addition: `:ring` emits
a closing bond from last to first atom. These entries are chosen for symmetry with the
query notation where `:paths`/`:rings` are the variable-length query counterparts.

**Q9 — Aromatic hint serialization.** The `a1` bond spec flag is valid as parser input
but is never emitted by canonical serializers. The invariant: a canonical ground term
either has an `:aromatic` section (with `:single` bonds between aromatic atoms, no `a1`)
or has explicit Kekulé bond orders (no `:aromatic` section, no `a1`). `a1` presence in
a document signals intermediate state. After kekulization, `:aromatic` is removed and
bonds gain explicit orders; after aromatization, `:aromatic` is populated and bonds
revert to `:single`. Mixed Kekulé/aromatic molecules (e.g. biphenyl with one ring each)
are represented natively with no special syntax — the two regions coexist without
normalization.
