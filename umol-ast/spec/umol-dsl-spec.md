# umol DSL specification

Normative definition of the molecule **EDN** surface, **contexts**, **molecule map** shape, **value expressions**, **bindings**, and **atom-string** / **bond-string** subgrammars.

The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, **RECOMMENDED**, **MAY**, and **OPTIONAL** in this document are to be interpreted as described in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

**Encoding.** Atom-string and bond-string payloads are Unicode text; **UTF-8** is the usual encoding on the wire and on disk. The subgrammars in **§7** define meaning only for **ASCII**: graphic characters U+0020–U+007E and the ASCII whitespace tokens named in **§7.1**. Any other code point in an atom-string or bond-string is **invalid**.

---

## 1. General

**Homoiconicity.** The same string grammars denote **ground** data and **constraint** patterns. A ground atom-string is a degenerate case of a query atom-string that matches exactly one interpretation.

**Relational molecule.** A molecule map (**§4**) is a set of named relations (atoms, localized bonds, optional dative / aromatic / multicenter sections, and global fields). Fragment composition is relation-wise merge where defined.

**EDN and rules.** EDN carries the relational structure. **Rule** evaluation (pattern **LHS** → product **RHS**, **§6**) is a separate computation layer: it **MAY** consume and produce molecule maps that use the same surface notation.

**Term algebra levels (non-normative sketch).** Four levels of expressiveness exist in this algebra, each a strict superset of the previous:

- **L1 — Ground**: fully instantiated molecules; no wildcards, binds, or logic. Every numeric slot is a concrete integer; element is a single symbol.
- **L2 — Constraint / Query**: adds wildcards (**`*`**), element/numeric sets, **`bool-expr`**, and **`?id`** references. Sufficient for substructure queries.
- **L3 — Rule**: adds **`element-bind`** / **`element-ref`**, cross-atom **`id`** scope (**§6**), and **`:guards`** on molecule maps. Sufficient for transformation rules.
- **L4 — Compound**: rules whose RHS produces molecule maps that are themselves L2/L3 terms; higher-order composition.

The subgrammars in **§5**, **§7** define forms that are syntactically valid across all levels; which forms are *semantically* allowed depends on which level is in force (**§3**).

**Case sensitivity.** **`atom-string`**, **`bond-string`**, and **`value-expr`** lexing (**§5**, **§7**) is **case-sensitive** throughout: e.g. **`#a`** and **`#A`** are distinct predicate tags; **`?x`** and **`?X`** are distinct **`id`**s; **`element-literal`** (**§7.4**) follows **IUPAC** element casing (**`Cl`**, **`Br`**, not arbitrary case folding). Implementations **MUST NOT** treat these fragments as case-insensitive (unlike **Fortran**-style languages).

---

## 2. EDN representation

Molecule data **SHOULD** use **EDN** maps and **Clojure-style** **keywords** for atom site labels.

An **atom literal** is an **EDN string** whose contents are an *atom-string* (**§7.3** / **§7.4**).

A **bond literal** in **full** form is an **EDN string** whose contents are a *bond-string* (**§7.5**). **§7.7** defines **keyword** shorthands (e.g. **`:single`**) that expand to an equivalent bond-string payload.

**`#` only appears inside the string.** Within an **atom-string** or **bond-string**, **`#`** starts a **predicate tag** (**§7.3**, **§7.5**); no reader-dispatch tag is used at the EDN layer.

**`:atoms`** is a **vector** of atom entries. Each entry is either:

- an atom literal (an EDN string carrying an **atom-string** payload),
- an alias reference (keyword matching a key in **`:atom-aliases`**), or
- an inline-id entry **`[`** *keyword* *atom-spec* **`]`** — a 2-element vector where the first element is a keyword id and the second is an atom literal or alias reference.

**`:bonds`** is a **vector** of bond entries (**§4**). Other optional keys are listed in **§4**.

---

## 3. Evaluation context

The same token grammars are used in every context; which **`value-expr`** shapes are legal depends on context.

| Context   | Constraints on strings |
|-----------|-------------------------|
| **Ground** | Wildcards, element **`*`**, element **brace-set**, **`?`**, **`bool-expr`**, **`element-bind`**, and **`element-ref`** are **invalid** unless this specification explicitly allows them for that slot. |
| **Query**  | Wildcards and constraints per **§5** and **§7**. |
| **Rule**   | Binds, guards, sets, boolean guards, and arithmetic per **§5**, **§6**, and **§7**. |

**Builder-oriented** use (expecting a unique ground resolution) and **query-oriented** use (selecting a set of matches) differ only in evaluation policy, not in the grammars.

**Ground / Query / Rule split.** **Ground** input uses the **same** grammars as **Query** / **Rule**; only allowed token shapes differ. Implementations **MAY** provide a **Ground-only** parser (or fast path) and a **full** parser. That split is **not** a second language: any string valid as **Ground** **MUST** be interpreted the same whether handled by the restricted implementation or the full one.

---

## 4. Molecule map

A **molecule map** has the following **normative** keys. Optional keys **MAY** be absent; absent means “not applicable” for structural sections. **`nil`** on optional value-typed keys means “unknown” in builder-oriented evaluation and “unconstrained” in query-oriented evaluation, when implementations distinguish those modes.

```
molecule-map ::=
  { :atoms              atom-list
    :bonds              bond-list
    [:dative-bonds      dative-bond-list]
    [:aromatic-systems  aromatic-system-list]
    [:multicenter-bonds multicenter-bond-list]
    [:noncovalent-bonds noncovalent-bond-list]
    [:stereo-atoms      stereo-atom-list]
    [:stereo-bonds      stereo-bond-list]
    [:atom-aliases      atom-alias-list]
    [:constraints       constraint-list]
    [:guards            [ logic-expr* ]]
  }

atom-alias-list      ::= [ keyword "atom-string" ]*

atom-list ::= [ atom-entry* ]

atom-entry ::= atom-spec
             | [ keyword atom-spec ]

atom-spec  ::= "atom-string" | keyword

bond-list ::= [ bond-entry* ]
dative-bond-list   ::= [ dative-bond-entry* ]

atom-ref ::= int | keyword

bond-entry ::= { [:id keyword] :a atom-ref :b atom-ref :type bond-spec }
             | [ atom-ref atom-ref bond-spec ]
dative-bond-entry   ::= { [:id keyword] :donor ( atom-ref | [ atom-ref+ ] )
                          :acceptor atom-ref :type dative-bond-spec }

bond-spec ::= "bond-string" | bond-keyword
dative-bond-spec ::= "dative-string" | dative-keyword

aromatic-system-list    ::= [ aromatic-system-entry* ]
multicenter-bond-list ::= [ multicenter-bond-entry* ]
noncovalent-bond-list ::= [ noncovalent-bond-entry* ]

aromatic-system-entry    ::= { [:id keyword] :atoms [ atom-ref+ ] [:electrons [ value-expr+ ]] :type "aromatic-string" }
multicenter-bond-entry ::= { [:id keyword] :atoms [ atom-ref+ ] [:electrons [ value-expr+ ]] :type "multicenter-string" }
noncovalent-bond-entry ::= { [:id keyword] :a atom-ref :b atom-ref :type noncovalent-bond-spec }

noncovalent-bond-spec ::= "noncovalent-string" | noncovalent-keyword
noncovalent-keyword ::= :h-bond | :halogen-bond | :chalcogen-bond | :ionic | :van-der-waals

stereo-atom-list ::= [ stereo-atom-entry* ]
stereo-bond-list ::= [ stereo-bond-entry* ]

stereo-atom-entry ::= { [:id keyword] :site atom-ref :ligands [ ligand-ref+ ] :type stereo-spec }
stereo-bond-entry ::= { [:id keyword] :site bond-ref :ligands [ ligand-ref+ ] :type stereo-spec }

ligand-ref  ::= atom-ref | [ :h atom-ref ] | [ :lp atom-ref ]
stereo-spec ::= "stereo-string" | stereo-keyword
stereo-keyword ::= :ccw | :cw | :z | :e
```

**`:type` is mandatory.** Every structural entry that has a DSL surface — **`bond-entry`**, **`dative-bond-entry`**, **`aromatic-system-entry`**, **`multicenter-bond-entry`**, **`noncovalent-bond-entry`**, **`stereo-atom-entry`**, **`stereo-bond-entry`** — **MUST** carry a **`:type`** key. The payload is a subgrammar string (or its EDN-keyword shorthand where defined); an entry without **`:type`** is a parse error.

Whether an empty string **`:type ""`** parses depends on the subgrammar:

- **`aromatic-string`** (**§7.9**) and **`multicenter-string`** (**§7.10**) admit the empty payload — the grammar is zero-or-more **`#…`** predicates, so **`:type ""`** is the canonical "no inline state" form.
- **`bond-string`** (**§7.5**), **`dative-string`** (**§7.11**), and **`noncovalent-string`** (**§7.12**) require a leading inherent-field token (bond order, dative order, noncovalent kind), so an empty payload is a parse error. Use the appropriate keyword shorthand (e.g. **`:single`**) or the literal token (e.g. **`"1"`**, **`"Hbd"`**).
- **`stereo-string`** (**§7.14**) requires a leading **`class`** token (**`Th`** / **`Ct`** / **`Sp`** / **`Tb`** / **`Oh`**) followed by a **`coset`**, so an empty payload is a parse error. Use a literal token (e.g. **`"Th1"`**) or a **`stereo-keyword`** (**`:ccw`** / **`:cw`** / **`:z`** / **`:e`**).

**Dative bond entry.** A dative bond entry binds a **single acceptor** to **one or more donors** (a coordination center): **`:acceptor`** names the atom accepting the electron pair(s); **`:donor`** names the donating atom, or a **vector** of donating atoms when several donate to the same acceptor. The leading **`order`** token of the **`dative-string`** payload (**§7.12**) records the number of donated pairs; one shared **`:type`** covers every donor→acceptor edge of the entry. The **`:acceptor`** and every **`:donor`** **MUST** reference distinct atom sites. The mandatory **`:type`** slot carries a **`dative-string`** (**§7.12**) — order plus optional aromatic flag (**`#a`**) and the ring-membership predicate (**`#R`**); its leading order parallels the bond-string's (**§7.5** / **§7.6**). The dative-string has **no** direction token — direction is expressed entirely by the **`:donor`** / **`:acceptor`** assignment.

**Multicenter entry.** The mandatory **`:type`** slot carries a **`multicenter-string`** payload (**§7.11**) encoding per-system charge, spin, and the optional asserted total electron count (**`#e<n>`**). The **`multicenter-string`** subgrammar is independent from **`aromatic-string`** even though they share the same predicate shape. A vacuous string (**`""`**) is admissible when the multicenter bond carries no inline state.

**`:electrons` vector (aromatic and multicenter entries).** Both **`aromatic-system-entry`** and **`multicenter-bond-entry`** **MAY** carry an optional **`:electrons`** key holding the **per-atom** electron contributions as a **single atomic value**: **either** **`:undetermined`** (the whole vector unknown) **or** a vector of **concrete integers**, one per member atom. A concrete vector **MUST** have the same length as the entry's **`:atoms`** vector — position **`i`** is the contribution of the atom at position **`i`** of **`:atoms`**. When **`:electrons`** is omitted the contributions are **`:undetermined`**. Its content is independent of the optional **`#e`** total on the entry's **`:type`** — when both are present, downstream validation **MAY** require **`sum(:electrons) == #e`** on ground inputs.

**Noncovalent kind.** A **`noncovalent-bond-entry`** **MUST** carry **`:type`**. The value is either a **`noncovalent-keyword`** (ground shorthand) or a **`noncovalent-string`** (**§7.13**) carrying the kind as an expression. The five **`noncovalent-keyword`** values expand to the corresponding ground literal in **§7.13** and are accepted wherever a **`noncovalent-bond-spec`** is expected.

**Stereo atom / stereo bond entry.** A **`stereo-atom-entry`** overlays a coordination-stereo configuration (tetrahedral and the higher geometries) on an atom site; a **`stereo-bond-entry`** overlays a cis/trans configuration on a bond site.

- **`:site`** names the bearing entity — an **`atom-ref`** for a **`stereo-atom-entry`**, a **`bond-ref`** for a **`stereo-bond-entry`**. A site **MUST** carry at most one stereo element (**§4.1**).
- **`:ligands`** is the **ordered** local reference frame against which the **`:type`** **`coset`** index is numbered; **order is significant**. Each **`ligand-ref`** is either a plain **`atom-ref`** (a neighbor atom) or a **virtual ligand** — **`[:h atom-ref]`** for an implicit hydrogen borne by the named atom, or **`[:lp atom-ref]`** for a lone pair on the named atom. For a **`stereo-atom-entry`** the bearing atom of every virtual ligand is the site atom; for a **`stereo-bond-entry`** it is the relevant double-bond terminus.
- **`:type`** carries a **`stereo-string`** (**§7.14**) — the **`class`** (**`Th`** / **`Ct`** / **`Sp`** / **`Tb`** / **`Oh`**) plus the **`coset`** index — or a **`stereo-keyword`** shorthand (**§7.14**). The coset index is a dense per-class arrangement number relative to the **`:ligands`** order, not a permutation rank.

**`:id`**. Each structural entry **MAY** include **`:id`** with an EDN **keyword** value. When present, **`:id`** values **MUST** be **pairwise distinct** across **all** entries in the **same** **molecule map** (every list combined).

**Keyword namespace disjointness.** All keyword-shaped identifiers within a single molecule definition — atom ids, atom alias names, structural entry **`:id`** values, and future keyword namespaces (bond alias names) — **MUST** be drawn from **mutually disjoint** namespaces. No two identifier kinds **MAY** share a keyword name within the same molecule map. Alias names **MUST NOT** be valid element symbols (**§7.4**).

**`:atom-aliases`**. The **`atom-alias-list`** defines named atom shorthands scoped to the enclosing molecule map. It is a flat vector of alternating keyword/atom-spec pairs. Each value **MUST** be an **EDN string** carrying an **atom-string** payload. An **`atom-entry`** that is a bare **keyword** (not a string and not in a **`[id entry]`** position) is an alias reference and **MUST** resolve to a key in **`:atom-aliases`**. Aliases are resolved at parse time; the resolved **`atom-string`** is substituted as if written inline. A reference to an undefined alias is an error. Alias definitions **MUST** be bijective: no two alias names **MAY** map to the same atom definition.

**`:constraints`**. Molecule-wide and per-entity constraints, cross-entity relational predicates, sub-pattern anchors, and boolean combinators live here. The canonical grammar appears in **§7.9**. Whole-molecule charge and spin assertions are written as **`{:charge-sum {:sum n}}`** and **`{:spin-sum {:spin spin-literal}}`** entries (omit `:atoms` to range over the whole molecule); a subset is selected by adding `:atoms [...]`. There is no top-level **`:charge`** or **`:spin`** key on the molecule map.

**Inline ids.** An **`atom-entry`** of the form **`[`** *keyword* *atom-spec* **`]`** assigns the keyword as an **id** to the atom at that position. Ids enable symbolic reference from bond endpoints (instead of positional index). Entries with and without ids **MAY** be freely mixed within the same **`:atoms`** vector.

**Endpoints.** Every atom site referenced from a structural relation **MUST** exist under **`:atoms`**, either by positional index (integer) or by id keyword:

- **`bond-entry`** **`:a`** and **`:b`**
- **`dative-bond-entry`** **`:donor`** and **`:acceptor`**
- **`noncovalent-bond-entry`** **`:a`** and **`:b`**
- every reference in an **`aromatic-system-entry`** or **`multicenter-bond-entry`** **`:atoms`** vector

**Positional index endpoints.** Let **`n`** be the length of **`:atoms`**. An integer endpoint **`i`** with **0 ≤ i < n** denotes the atom at position **`i`**. A keyword endpoint denotes the atom whose inline id is that keyword.

**Empty molecule.** The **`atom-list`** **MAY** have length **0** (**`[]`**). If **`:atoms`** is empty, **`:bonds`** **MUST** be **`[]`**, and **`:dative-bonds`**, **`:aromatic-systems`**, **`:multicenter-bonds`**, **`:noncovalent-bonds`**, **`:stereo-atoms`**, and **`:stereo-bonds`** **MUST** be absent or **empty** lists — no bond-like or stereo entry **MAY** name a site that is not in **`:atoms`**.

**`bond-keyword`** and **`dative-keyword`** (shorthands) are defined in **§7.7**.

### 4.1 Structural validity (within one map)

These rules apply **within** a single **molecule map**. **Constraints across** relation kinds (e.g. the same atom pair in **`:bonds`** and **`:dative-bonds`**) are **not** specified here.

**`:bonds` (localized).** The list **MUST NOT** contain two **`bond-entry`** values with the same **unordered** pair of atom sites **{`:a`, `:b`}** (endpoints as a set).

**`:dative-bonds`.** No **donor→acceptor** edge may repeat — counting each donor in an entry's **`:donor`** list against that entry's **`:acceptor`**, across all entries. A donor→acceptor edge and the reverse acceptor→donor edge between the **same** two atoms also violate this rule.

**`:aromatic-systems`.** For every two distinct **`aromatic-system-entry`** values, the sets of keywords in their **`:atoms`** vectors **MUST** be disjoint. Aromatic systems **MUST NOT** share an atom.

**`:noncovalent-bonds`.** For every two distinct **`noncovalent-bond-entry`** values, the sets **{`:a`, `:b`}** **MUST** be disjoint. Noncovalent bonds **MUST NOT** share an atom with another noncovalent bond in the same map.

**`:stereo-atoms` / `:stereo-bonds`.** The lists **MUST NOT** contain two **`stereo-atom-entry`** values with the same **`:site`** atom, nor two **`stereo-bond-entry`** values with the same **`:site`** bond. Each atom (resp. bond) **MUST** be the site of at most one stereo element.

**`:guards`** **MAY** appear only in **Rule** context; it holds predicates over variables bound on the **LHS** that are not expressed inline in atom-strings.

---

## 5. Value expressions

**`value-expr`** appears **only** inside a **predicate payload** (**§7.3**, **§7.5**) after **`#` *tag***. The character **`#`** **MUST NOT** appear inside a **`value-expr`** (it is reserved for starting the next predicate).

```
digit      ::= '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
nat        ::= digit+
decimal-tail ::= digit*

nat-set    ::= '{' nat (',' nat)* '}'

An empty **`nat-set`** **`{ }`** is **invalid**.

value-expr ::= '*'
             | nat-set
             | nat
             | '?' id                  (* top-level Ref — bind reference     *)
             | '?' id '::' nat-set     (* top-level Bind — named domain     *)
             | bool-expr

bool-expr  ::= or-expr

or-expr    ::= and-expr ( '|' and-expr )*
and-expr   ::= not-expr ( '&' not-expr )*
not-expr   ::= '!' not-expr
             | '(' bool-expr ')'
             | rel-expr

rel-expr   ::= mem-expr ( rel-op mem-expr )?

mem-expr   ::= add-expr ( '::' nat-set )?

add-expr ::= mult-expr ( add-op mult-expr )*
mult-expr ::= unary-expr ( mult-op unary-expr )*

unary-expr ::= sign* base-expr

sign ::= '+' | '-'

base-expr ::= nat
               | '?' id
               | '(' add-expr ')'

rel-op ::= '<=' | '>=' | '==' | '<' | '>'
add-op ::= '+' | '-'
mult-op ::= '*' | '/' | '%'

id  ::= [a-zA-Z][a-zA-Z0-9_]*
```

**`(`** **`bool-expr`** **`)`** appears under **`not-expr`** so logic can be **grouped** without being mistaken for **arithmetic** parentheses.

**Top-level `nat`.** A **`nat`** forms a complete top-level **`value-expr`** only when it is **cut** by **end of the substring being tokenized** or by the **next predicate** (**`#`** *tag* on the atom-string / bond-string), after optional whitespace — e.g. **`#h1#v3`** yields **`1`** then **`3`**. If the next non-whitespace character is anything else (e.g. **`+`** in **`1+2`**), parsing **MUST NOT** treat the **`nat`** as this alternative; it falls through to **`bool-expr`**. Implementations **MAY** represent this form as **`Lit`** distinct from a trivial relational **`bool-expr`**. The common **Ground** case (**`#h3`**, **`#v0`**, …) is typically this shape.

**Top-level `nat-set`.** A **`value-expr`** may be **only** a **`nat-set`** (after the usual ignored whitespace between **`value-expr`** tokens, **§7.1**). It denotes a **finite numeric disjunction** for the **one** quantity fixed by the enclosing predicate tag (**§7.3**, **§7.5**): that quantity **MUST** equal one of the listed **`nat`** values. This is the same constraint **shape** as a top-level **`nat-set`** in **bond-string** **`order`** (**§7.5**) and **`element-set`** for the **element** prefix, applied at the **predicate payload** level (e.g. **`#h{1,2,3}`** with payload **`{1,2,3}`**). It **MUST NOT** introduce a numeric **`?id`**; implementations **MAY** lower it to **`bool-expr`** internally. The form ***arith* `::` *nat-set*** on **`mem-expr`** is unchanged: it constrains the **arithmetic** value on the left of **`::`**, not an implicit slot quantity by bare **`{…}`** alone.

**Top-level `?` *id* and `?` *id* `::` *nat-set*.** A **`value-expr`** may be **only** **`?` *id*** (a numeric **bind reference**) or **`?` *id* `::` *nat-set*** (a **named-bind** with a finite admissible domain). These shapes parse at the **value-expr** level before falling through to **`bool-expr`**, and produce the AST variants **`ValueAst::Ref`** and **`ValueAst::Bind`** respectively. Inside a compound expression (e.g. **`?h + 1`**, **`?h == 0`**), the same **`?` *id*** appears as **`Expr::Var`** inside **`bool-expr`** — the discriminator is whether the surrounding context is the whole value or an operand of a larger operator.

**Paren-transparency for top-level bind/ref.** Outer parentheses around a top-level **`?` *id*** or **`?` *id* `::` *nat-set*** are **optional** and **semantically transparent**: implementations **MUST** accept the bare forms and any nesting depth of outer parens (**`(?h)`**, **`((?h))`**, **`(?h :: {1,2})`**, **`((?h :: {1,2}))`**) as identical AST. The **canonical** rendered form is **bare** (no outer parens). Disambiguation against larger expressions like **`(?h + 1)`** or **`(?a :: {0}) & 0 <= 0`** is handled by requiring a **terminator** (end-of-payload or next **`#`** predicate) after the parenthesized bind/ref before the arm fires; otherwise the parens are interpreted as **`bool-expr`** grouping (**§5.1**).

**`unary-expr`** is **`sign`*** **`base-expr`**: zero or more leading **`+`** / **`-`**, then **`nat`**, **`?id`**, or parenthesized **`add-expr`**. Examples: **`#c+1`**, **`#c-2`**, **`#c--1`**. A **`sign`** with **no** following **`base-expr`** is **invalid** in the general grammar; **`#c`** additionally accepts a payload consisting **only** of **`+`** or **`-`** (after trimming whitespace) as **+1** or **−1** (**§7.3**).

**Equality** is **`==`** (not **`=`**). **Finite numeric membership** uses the **`::`** token and a **`nat-set`**: **`?h + 1 :: {2,3}`** parses as **`(?h + 1) :: {2,3}`** — the full **additive** form is built before **`::`**, which sits **below** **relations** and **logic** only (same layering as former **`in`**).

**Meaning of `::`.** The same token **`::`** appears in two shapes:

- In **`element-bind`** (**§7.4**), **`?` *id* `::` *element-set*** means: the nominal variable is constrained to **membership in a set of element symbols** (chemical **`element-literal`** values).
- In **`mem-expr`** (inside **`value-expr`**), ***arith* `::` *nat-set*** means: the left-hand **arithmetic** value **MUST** be a member of the **numeric** set. After matching, concrete values **MUST** fit the slot’s type: **`u8`** for most atom/bond numeric predicates (**§7.2**), **`i8`** for formal charge (**`#c`**), **`u32`** for isotope mass (**`#i`**).

**Ground** **`value-expr`** (predicate payloads where allowed) are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (non-empty, entries valid for the slot), with optional leading **`sign`** sequence on a bare **`nat`** / **`decimal-tail`** only: no **`?`**, **`::`**, relations, logic, **`*`**, or **parentheses** (**§5.4**), except **`#c`** (**§7.3**) also allows a payload consisting **solely** of **`+`** or **`-`** (**+1** / **−1**). Implementations **MAY** use a restricted parser for **Ground** (**§3**).

**Numeric evaluation.** Where **`add-expr`** / **`mult-expr`** are evaluated to **concrete** counts (e.g. **Rule** RHS), intermediate and final **numeric** results **SHOULD** be computed in a range consistent with the target slot (**§7.2**, typically **`u8`** where that slot is **`u8`**). Values **outside** the slot’s range **MUST** be rejected at validation.

**`bool-expr`** is the **`value-expr`** form for **constraints** in **Query** / **Rule** (only in slots implementations allow). A **plain numeral** payload is an **`add-expr`** that is only a **`nat`**. **Parentheses** **`(`** **`add-expr`** **`)`** group **arithmetic** only.

**Additional form (same precedence as other `add-expr` operands):** **`nat` `add-op` `?` `id`** (e.g. **`4-?h`**) is covered by **`mult-expr`** / **`add-expr`**; a **parenthesized** spelling **`(` `nat` `add-op` `?` `id` `)`** is equivalent and **MAY** be used for clarity. **Extending** to arbitrary **binary** nesting of **`arith-expr`** is listed in **§7.8**.

Inside **`nat-set`**, optional whitespace is allowed only adjacent to commas, like **element-set** (**§7.4**).

**`?id`** in **`bool-expr`** introduces or uses a **numeric** bind (**§6**). **`element-ref`** (**§7.4**) is **not** **`bool-expr`**. **Variables are not surface-typed**; illegal combinations are rejected when lowering / validating, not by this grammar’s token shapes.

### 5.1 Precedence and parentheses (infix)

Binding strength **decreases** down the table (rows tie left-to-right within the row where conventional):

| Precedence | Constructs |
|------------|------------|
| tightest | **`(`** **`add-expr`** **`)`**, **`(`** **`bool-expr`** **`)`**, **`base-expr`**, **`sign`*** prefix on **`base-expr`** |
| multiplicative | **`*` `/` `%`** |
| additive | **`+` `-`** |
| membership | **`::`** **`nat-set`** as suffix on **`add-expr`** (parse **`?h + 1 :: {2}`** as **`(?h + 1) :: {2}`**) |
| (entire **`value-expr`**) | top-level **`nat-set`** — disjunction of numeric literals for the predicate slot; not an infix operator |
| relational | **`<` `>` `<=` `>=` `==`** between **`mem-expr`** operands |
| unary | prefix **`!`** (boolean) |
| conjunction | **`&`** |
| loosest | **`|`** |

**`!`** binds **tighter** than **`&`** and **`|`** (e.g. **`! ?h == 1 & ?v == 2`** parses as **`(! (?h == 1)) & (?v == 2)`**).

Use **`(`** **`bool-expr`** **`)`** to override precedence for logic.

**Parsing a predicate payload:** the **payload** substring (**§7.1**) is parsed as **`value-expr`** (or a **Ground** subset). It **MUST NOT** contain **`#`**.

### 5.2 Truth and repeated `id` on one atom-string

For one **atom-string** at match time, **numeric** **`?id`** values in **`bool-expr`** are fixed from the **matched** target atom so that **all** constraints on that **`id`** **hold together** (**§6**). Prefix **`!`** is classical boolean negation on its operand **`bool-expr`**. Example: **`#h(!?h==1)`** or **`#h( ! ?h == 1 )`**: **`h`** is the implicit H count; the guard holds iff **`h ≠ 1`**.

The same **`id`** may appear in **several** payloads on **one** **atom-string**; one value satisfies **every** use. Example: **`#h(?h <= 4)#v(4-?h)`** — **`h ≤ 4`** and **`v = 4 − h`**.

**Cross-atom** reuse of the same **`id`** is **not** fixed here; implementations **SHOULD** document whether and how it is allowed.

### 5.3 `decimal-tail` and omitted numeral = 1

**`decimal-tail`** is **`digit`*** (**§5**). Its numeric meaning is **1** when there are **no** digits; otherwise the usual base‑10 value of the digit sequence (no separate **`nat`** required for the all-zero-digits case).

When a **predicate** (**§7.3**, **§7.5**) allows a **decimal-only** payload and the payload is **only** a **`decimal-tail`** (not **`*`**, not **`(`**, not **`bool-expr`**, not a **`sign`**-only **`#c`** form), the **omitted** form (zero digits after the tag) means **1**. Lexing **MUST** take the **longest** contiguous run of **`digit`**s as that numeral (**greedy**). **Special predicate** payloads (**§7.3**) are **not** **`decimal-tail`**-only forms.

| Atom tag | Omitted numeral = 1 (decimal-only payloads) |
|----------|-----------------------------------------------|
| **`#c`** | **no** — charge **MUST** be explicit (**`#c0`**, **`#c+`**, **`#c-`**, **`#c+2`**, **`#c-2`**, …); empty **`#c`** is **invalid** |
| **`#h` `#n` `#u` `#s` `#v` `#d` `#t`** | yes, when the payload is decimal-only |
| **`#a`** | yes when decimal-only; **`#a*`**, **`#a+`**, **`#a!`** are **special** (**§7.3**) |
| **`#m`** | yes when decimal-only; **`#m*`**, **`#m+`**, **`#m!`** are **special** (**§7.3**) |
| **`#i`** | yes, when the payload is decimal-only; bare **`#i`** denotes isotope mass **1** |
| **`#R`** | yes when decimal-only (incl. the sized **`#R(<size>)`** form); **`#R*`**, **`#R+`**, **`#R!`** are **special** (**§7.3**) |

Bond predicates that use **decimal-only** payloads follow the same **`decimal-tail`** rule where applicable (**§7.5**).

In **Query** and **Rule**, any predicate slot that allows a full **`value-expr`** may use **`bool-expr`**, **`*`**, top-level **`nat-set`**, **`decimal-tail`**, etc., as allowed for that tag.

### 5.4 Wildcards, sets, logic, arithmetic

- The **`*`** **wildcard** is allowed in **`value-expr`**, **`element`**, and **`order`**
- **`bool-expr`**: **infix** **`&` `|` `!`**, **relations**, **`::`**, **`+ - * / %`**, unary **`-`**, **`?id`**, **`nat`**, **`(`** **`add-expr`** **`)`**.

**Ground:** no **`bool-expr`** (no **`?`**, **`::`**, relations, logic), no **`element-bind`**, no **`element-ref`**, no top-level negation (**`!`** *literal* or **`!`** *set*); predicate payloads are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (and tag-specific literals such as **`#i=`** for natural isotope) only where allowed.

**Query:** **`bool-expr`** where allowed; **`decimal-tail`**; **element** / **order** extensions as allowed.

**Rule:** full **`value-expr`**; **element** may use **`element-bind`** / **`element-ref`** (**§6**).

**`<` `>` `<=` `>=` `==`** appear **only** inside **`value-expr`** (predicate payloads). **Dative** donated / accepted pair counts use predicates **`#d`** / **`#t`** (**§7.3**), not bare **`<` `>`** at the top level of the atom-string.

---

## 6. Match semantics and bindings

**Ground molecule, pattern LHS.** The target is **ground** (fully instantiated). The **LHS** of a rule (or query) may still contain **wildcards**, **sets**, **binds**, and **guards**: that is **pattern** data, not an indeterminate molecule.

### 6.1 Inherent fields and derived predicates

**Inherent fields.** Each AST form — atom, localized bond, aromatic system, multicenter bond, dative bond, noncovalent bond, stereo atom, stereo bond — carries a fixed set of **inherent fields**. An inherent field's value **identifies** the entity at that slot. An entity is **ground** iff every inherent field holds a single concrete value (a literal; not a wildcard, set, bind, ref, or unresolved symbolic state). Nothing else affects grounding.

| Form | Inherent fields |
|------|-----------------|
| atom | element, isotope mass (**`#i`**), charge (**`#c`**), implicit hydrogens (**`#h`**), lone pairs (**`#n`**), spin (unpaired **`#u`**, multiplicity **`#s`**) |
| localized bond | order, charge (**`#c`**), spin (**`#u`**, **`#s`**) |
| aromatic system | charge (**`#c`**), spin (**`#u`**, **`#s`**), π-electron count (**`#e`**) |
| multicenter bond | charge (**`#c`**), spin (**`#u`**, **`#s`**), electron count (**`#e`**) |
| dative bond | a single **`:acceptor`** and its one-or-more **`:donor`** atoms — the assignment on the map entry (**§4**) — plus the leading **`order`** token of the dative-string (number of donated electron pairs; **§7.12**). |
| noncovalent bond | interaction kind (**`:h-bond`**, **`:halogen-bond`**, **`:chalcogen-bond`**, **`:ionic`**, **`:van-der-waals`**) |
| stereo atom | coordination **`class`** (geometry) and **`coset`** configuration index (the **`:type`** payload, **§7.14**). The bearing **`:site`** atom and ordered **`:ligands`** frame (**§4**) are the relation's participants, not payload fields. |
| stereo bond | cis/trans **`class`** and **`coset`** configuration (the **`:type`** payload, **§7.14**). The bearing **`:site`** bond and ordered **`:ligands`** frame (**§4**) are the relation's participants. |

**Derived predicates.** Every predicate admitted in the DSL that is not an inherent field is a **derived predicate** — a topological query evaluated against the target graph once an embedding is proposed. This includes per-atom **`#D`**, **`#X`**, **`#x`**, **`#H`**, **`#R`** (**§7.3**); the bond-namespace **`#R`**; per-aromatic, per-multicenter, per-dative ring-membership predicates; and the molecule-wide entries of **§7.9**. Derived predicates **filter** matches; they do **not** carry identity and **do not** affect grounding. Adding a derived predicate — even a wildcard-valued one — to a pattern never makes a ground target stop being ground.

**Symmetry-derived stereo predicates.** The stereo entity predicates — **`#p`** ligand symmetry, **`#o`** topicity, **`#g`** stereogenicity (**§7.14** / **§7.9**) — are **derived** from the resolved molecule's **graph automorphisms** (the ligand-frame symmetry group of the stereo element), not from the local string. As derived predicates they **filter** matches and **do not** affect grounding: a stereo element is ground iff its **`class`** + **`coset`** are concrete (**§6.1** table), regardless of which **`#p`**/**`#o`**/**`#g`** assertions it carries. The validator computes the molecule-wide symmetry once on the **resolved** AST and **cross-checks** the derived value against each stored constraint — when both are ground and inconsistent (including a kind/degree mismatch, or a **`'`** value on an achiral class), the molecule is rejected — exactly as the topology-derived fields are cross-checked against the stored inherent fields. **`#f`** fluxionality is a stored dynamical assertion (not derivable from a static graph); it is matched as a stored predicate, not cross-checked.

### 6.2 Pattern–target match

**Match as solution-set inclusion.** Each attribute slot has a **solution set** — the set of ground values the slot admits. A **literal** (e.g. **`C`**, **`3`**, **`+1`**) admits exactly itself; a **set** (**`{C,N}`**, top-level **`nat-set`**) admits its members; a **negation** (**`!H`**, **`!12`**) admits everything in the slot's value domain *except* the named literal; a **negative set** (**`!{F,Cl}`**, **`!{12,13}`**) admits the complement of the listed entries; a **wildcard** (**`*`**) admits everything in the slot's value domain; a **`bool-expr`** admits every value for which the expression holds (**§5**); a **special-symbolic** payload (**`#i=`**, **`#a*`**, **`#a+`**, **`#a!`**, **`#m*`**, **`#m+`**, **`#m!`**, **`#R*`**, **`#R+`**, **`#R!`**) admits only its named symbolic state (**§7.3**). For a given slot, the **pattern** matches the **target** iff `solution-set(pattern)` ⊇ `solution-set(target)` — the pattern admits every value the target admits. Match is **not** symmetric.

| pattern kind | target kind | matches iff |
|--------------|-------------|-------------|
| wildcard (**`*`**) | any | always |
| non-wildcard | wildcard | never (target's set is strictly larger) |
| literal | literal | values equal |
| literal | set | set is exactly that singleton |
| set | literal | literal is a set member |
| set **P** | set **T** | **T ⊆ P** |
| **`bool-expr`** | literal | expression holds on the literal (**§5**) |
| **`bool-expr`** | set | expression holds on **every** set member |
| **`bool-expr`** | wildcard / **`bool-expr`** | **undefined** in general; implementations **MAY** reject |
| special-symbolic **s** | special-symbolic **t** | **s == t** |
| special-symbolic | literal / set | never (disjoint domains) |

**Element matching.** Parallels the above: **`element-literal`** against **`element-literal`** by equality; **`element-set`** against a target iff the target's admissible symbols are a subset; **`element-bind`** behaves as its inner **`element-set`** for the match (the nominal binding is a side effect, not a match filter, **§6.3**); **`element-ref`** outside a resolved rule-scope binding context matches nothing.

**Noncovalent-kind matching.** Same shape as element matching, over the five-literal domain **`{:h-bond, :halogen-bond, :chalcogen-bond, :ionic, :van-der-waals}`**.

**Molecule-level match.** A molecule-map pattern matches a target iff (a) every **`atom-string`** matches its corresponding target atom-string field-wise — element and every inherent-field predicate payload, per the rules above; (b) every **`bond-string`** matches field-wise; (c) each structural relation (**`:aromatic-systems`**, **`:multicenter-bonds`**, **`:dative-bonds`**, **`:noncovalent-bonds`**, **`:stereo-atoms`**, **`:stereo-bonds`**) matches per its own inherent fields; (d) every derived predicate holds on the resulting embedding (**§6.1**). Any failure rejects the embedding.

### 6.3 Bindings

**One binding per match.** For a **fixed** embedding of the LHS pattern into the target (one way of mapping pattern sites to concrete atoms/bonds that satisfies all constraints), each **`id`** introduced by **`element-bind`** or by **`?id`** in **`bool-expr`** / **`value-expr`** has **exactly one** value in that match. There is no separate “CSP over the whole molecule without choosing an embedding”: the engine first chooses an embedding (or enumerates them — see below), then **the match binding** is fixed — i.e. the mapping from each such **`id`** to its concrete value for that embedding (numeric for **`?id`**, element symbol for nominal binds).

**Multiple results from one ground target.** Ambiguity does **not** require an indeterminate target. The **same** ground molecule can admit **several** distinct **embeddings** of the **same** LHS (e.g. two equivalent substituents). Each embedding yields its own **match binding**. Whether the rule **fires once**, **once per embedding**, or **aggregates** products is **policy** for the rule evaluator, not fixed by this specification.

**Nominal vs numeric.** **`element-bind`** and **`element-ref`** carry **element**-symbol values. **`?id`** in **`bool-expr`** carries **numeric** bind / use for that attribute. Arithmetic applies only to **numeric** **`id`** values. **Nominal** **`id`** may be **re-used** on the RHS via **`element-ref`** with the same name; no arithmetic on those.

**Identifier scope.** On **one** **atom-string**, the same numeric **`id`** may appear in **multiple** predicate payloads; all uses denote **one** value and are **jointly** satisfied (**§5.2**). Whether **`id`** may also be shared **across** atom-strings on a rule LHS (or RHS) is **not** fixed here; implementations **SHOULD** document cross-atom **`id`** rules. **Order** of **predicates** (**§7.3**) is arbitrary; evaluation **MUST** treat all constraints on **`id`** as **simultaneous**, not sequential by textual order.

---

## 7. Subgrammars

### 7.1 Whitespace and `#`

- **ASCII whitespace** (space, tab, CR, LF) is **not** significant: it **MAY** appear between **tokens** and is **ignored**, except where this section **forbids** it.
- **Leading and trailing** whitespace on the whole **atom-string** or **bond-string** is ignored.

**Whitespace is forbidden:**

- Between **`#`** and the **tag letter** of a **predicate** (**§7.3**, **§7.5**): **`#h`** is valid; **`# h`** is **invalid**.
- Inside multi-character operators: **`<=`**, **`>=`**, **`==`**, **`::`** (**§5**).
- Between **`?`** and the first character of **`id`** in **`?id`** (numeric bind) and in **`element-ref`** / **`element-bind`**.

**`#` (U+0023).**

- **Inside** an **atom-string** or **bond-string**, **`#`** **MUST** appear **only** as the first character of a **predicate** (**`#` *tag***). **`#`** **MUST NOT** appear inside **`element`**, **`order`**, or inside any **`value-expr`** / **payload** substring.

**Payload extraction.** A **predicate** is **`#`**, one **tag** character **`[A-Za-z_]`**, and a **payload** consisting of all following characters up to (but not including) the **next** **`#`** or **end of string**, after **whitespace normalization** for the purpose of **tokenizing** the payload as **`value-expr`**: the payload text **MAY** contain ignored whitespace between **`value-expr`** tokens as in **§5**. The **payload** **MUST NOT** contain **`#`**.

**Examples (atom):** **`C`**, **`C#h3`**, **`C#h*`**, **`!H`**, **`!{F,Cl}`**, **`?e`**, **`?e :: {Cl,Br}`**, **`?e :: !{F,Cl}`**, **`C#a*`**, **`C#a !`**, **`C#c+`**, **`C#c-`**, **`C#c +`**.

- A **`nat`** and an **`id`** contain **no** internal whitespace.
- A **relational** token is **`<=`**, **`>=`**, **`==`**, or a **single** **`<`** or **`>`** that is **not** part of **`<=` `>=`**. **Multi-character** tokens are one lexical unit.
- An **arithmetic** token is **`+`**, **`-`**, **`*`**, **`/`**, **`%`**. Leading **`+`** / **`-`** on a **`base-expr`** are **`sign`** tokens (**§5**); binary **`+`** / **`-`** appear between **`mult-expr`** operands.
- **Inside** an **element** or **order** **brace set** `{…}`, optional whitespace is allowed **only** immediately before or after a comma separating entries. No whitespace inside an **`element-literal`** or **`order-entry`** (**`nat`**).

**`=` (U+003D).** Plain **`=`** is **not** an operator in **`value-expr`**; equality is **`==`**. **`=`** **MUST NOT** appear in a **Ground** string in a **`value-expr`** position (no **`bool-expr`**). It **MAY** appear inside payloads only as part of **`==`** or other defined tokens.

**Vacuous-payload elision (canonical rendering).** Implementations **MAY** elide vacuous payloads — predicates whose value is **`Undetermined`** (`*` in surface form), and inherent fields with prefixed tags (**`#c`**, **`#u`**, **`#s`**, **`#e`**, …) whose value is **`Undetermined`** — from the canonical rendered form. Both forms remain admissible **on parse**, so a renderer that elides them still accepts a string in which they appear explicitly. **Leading unprefixed** inherent fields (bond **order**, atom **element**, noncovalent bond **type**) are **exempt** from this elision because they fix the entity-string's start position; for these, **`Undetermined`** **MUST** render as **`*`**. Round-trip identity at the AST level therefore holds only for ASTs whose constraint and inherent-field payloads are non-vacuous (or whose vacuous payloads sit on a leading unprefixed field).

### 7.2 Numerical limits

**Chemical elements.** Any **`element-literal`**, any entry in an **`element-set`**, and any **nominal** binding or reference (**`element-bind`**, **`element-ref`**) **MUST** refer only to elements from **hydrogen** (**H**) through **oganesson** (**Og**). Implementations **MUST** reject symbols outside that range in **Ground**; **Query** / **Rule** **SHOULD** use the same restriction unless explicitly documented otherwise.

**Charges.** **Formal charge** on atoms (**`#c`**), **formal bond charge** (**`#c`** on **bond-string**), **aromatic-system charge** (**`#c`** on **aromatic-string**, **§7.10**), and atom-subset charge sums **`{:charge-sum {:atoms [...] :sum n}}`** (**§7.9**) where integral **MUST** fit a **signed 8-bit** integer (**−128…127**). The **`#c`** payload is a **`value-expr`** (or **Ground** subset) that evaluates to the signed charge, including the **special** forms **`+`** / **`-`** for **±1** (**§7.3**), e.g. **`#c2`**, **`#c-2`**, **`#c+`**, **`#c-`**.

**Isotope mass number.** The numeric value carried by **`#i`** in **Ground** **MUST** fit an **unsigned 32-bit** integer.

**Unsigned 8-bit numeric slots (0…255).** After parsing, **concrete Ground** values for the following **MUST** fit **`u8`** (atom tags, **§7.3**):

| Slot | Atom predicate |
|------|------------------|
| Implicit H count | **`#h`** |
| Lone pairs (nonbonding) | **`#n`** |
| Unpaired electrons | **`#u`** |
| Spin multiplicity (2S+1) | **`#s`** |
| σ localized valence | **`#v`** |
| Dative donated pairs | **`#d`** |
| Dative accepted pairs | **`#t`** (“accepted”) |
| Ring membership (total / per size) | **`#R`** / **`#R(<size>)`** |
| Aromatic π contribution (numeric) | **`#a`** |
| Multicenter valence | **`#m`** |

**Bond-string** **`#u`** / **`#s`** / **`#c`** use the **bond** namespace (**§7.5**); meanings parallel **unpaired electrons**, **multiplicity**, and **bond formal charge**.

**Aromatic-string** **`#c`** / **`#u`** / **`#s`** / **`#e`** use the **aromatic-system** namespace (**§7.10**); **`#e`** denotes the total π-electron count (**`u8`**), other tags parallel the bond namespace.

**Lexical** **`nat`** in the grammar is unbounded; **Ground** validation **MUST** reject values outside the **u8** (or **i8** / **u32** as above) range for the corresponding slot.

**Bond order** (**§7.6**) uses a **discrete** model; **fractional** bond orders **MUST NOT** appear in the **bond-string**. **Aromatic** connectivity **MUST NOT** be encoded as a bond **order**; use the molecule map’s **`:aromatic-systems`** section (**§4**) and ordinary **`:bonds`** entries.

### 7.3 Atom subgrammar

```
atom-string ::= element atom-predicate*

atom-predicate ::= '#' tag payload

tag ::= [A-Za-z_]
```

- **`element`** is first (**§7.4**).
- **Zero or more** **`atom-predicate`** units follow. **Optional** ASCII whitespace **MAY** appear **between** **`element`** and the first **`#`**, and **between** successive predicates.
- **At most one** predicate per **tag letter** per **`atom-string`** (each row of the table below is a **kind**).
- **Canonical order** of predicates after **`element`** (stable serialized form): **`#i`**, **`#c`**, **`#h`**, **`#n`**, **`#u`**, **`#s`**, **`#v`**, **`#d`**, **`#t`**, **`#a`**, **`#m`**, then the uppercase derived predicates **`#D`**, **`#X`**, **`#x`**, **`#H`**, **`#R`**, **`#T`**. Implementations **MAY** specify further ordering for fields not listed here.

**`payload` parsing.** After trimming leading / trailing whitespace on the **payload** substring, parse as follows:

1. **`#c`**, **Ground** or **Query** / **Rule**: if the trimmed payload is **exactly** **`+`** or **`-`**, the formal charge is **+1** or **−1** (same meaning as **`#c+1`** / **`#c-1`**). Otherwise parse as **`value-expr`** (**§5**) (or the **Ground** subset in **§5.4**).
2. **`#i`**: parsed by a **dedicated isotope subgrammar** (see below), not as **`value-expr`**.
3. **Any other tag**: parse the payload as **`value-expr`** (or **Ground** subset) unless the payload matches a **special** form below.

**`#i` isotope subgrammar.** The isotope-mass slot uses its own subgrammar, not **`value-expr`**, because isotope mass numbers are tagged enum-like and have no arithmetic-on-numerics use. Empty payload (bare **`#i`**) denotes mass **1** (per §5.3 decimal-tail).

```
isotope-payload ::= '='                           (* Natural — naturally most abundant *)
                  | '*'                            (* Undetermined — wildcard           *)
                  | signed-int                     (* Lit                                *)
                  | nat-set                        (* Set — finite mass disjunction      *)
                  | '!' signed-int                 (* Not — cofinite singleton          *)
                  | '!' nat-set                    (* NotSet — cofinite multi           *)
                  | '?' id                         (* Ref — bind reference              *)
                  | '?' id '::' isotope-domain     (* Bind — named domain                *)

isotope-domain  ::= nat-set                        (* MemOp::In                          *)
                  | '!' signed-int                 (* MemOp::NotIn (singleton)           *)
                  | '!' nat-set                    (* MemOp::NotIn                       *)
```

**Paren-transparency.** Outer parentheses around **`?` *id*** or **`?` *id* `::` *isotope-domain*** are **optional** and **semantically transparent** (same rule as element §7.4 and value-expr §5). Canonical render is bare.

**Natural is its own channel.** **`=`** (Natural) **does not unify** with numeric variants in the lattice: **`Natural ∧ Lit(n) = ⊥`**, **`Natural ∨ Lit(n) = Undetermined`** for any **`n`**. Natural is the "no specific mass committed" state and is **disjoint** from any explicit mass number. **Ground** isotope is either **`Natural`** or a single **`Lit(n)`**.

**Special predicate payloads** (trimmed; **not** parsed as boolean **`!`** — these are **opaque** lexemes for the given tag):

| Form | Tag | Meaning |
|------|-----|---------|
| **`=`** | **`#i`** | **Natural isotope**: mass number of the **naturally most abundant** isotope of **`element`**. This is the default / expected isotope for each element. |
| **`*`** | **`#h`** | **Wildcard** implicit H count (**Query** / **Rule**). |
| **`*`** | **`#a`** | **No constraint** on aromatic π contribution — equivalent to omitting **`#a`** entirely. |
| **`+`** | **`#a`** | **Sugar** for the constraint **`?a >= 0`**: atom is a member of some aromatic system with an unspecified π contribution (**Query** / **Rule**). |
| **`!`** | **`#a`** | Atom is **not** a member of any aromatic system. Distinct from **`#a0`**: a **`#a0`** atom *is* in an aromatic system and contributes **zero** π electrons (e.g. a carbocation with an empty p orbital participating in a ring current); a **`#a!`** atom has no aromatic membership at all. Cross-checked against **`:aromatic-systems`** membership during validation (inconsistency is a validator error, not a parse / ground error). |
| **`*`** | **`#m`** | **No constraint** on multicenter valence — equivalent to omitting **`#m`** entirely. |
| **`+`** | **`#m`** | **Sugar** for the constraint **`?m >= 0`**: atom is a member of some multicenter bond with an unspecified multicenter-valence count (**Query** / **Rule**). |
| **`!`** | **`#m`** | Atom is **not** a member of any multicenter bond. Parallels **`#a!`**; cross-checked against **`:multicenter-bonds`** membership during validation. |
| **`*`** | **`#R`** | **No constraint** on ring count — equivalent to omitting **`#R`** entirely. |
| **`+`** | **`#R`** | **Sugar** for the constraint **`?r >= 1`**: atom is in **at least one** ring (**Query** / **Rule**). With a size, **`#R(<size>)+`** means at least one ring of that size (SMARTS **`r<size>`**). |
| **`!`** | **`#R`** | Ring count **0**: atom is in **no** ring (acyclic). With a size, **`#R(<size>)!`** means no ring of that size. |
| **`*`** | **`#T`** | **No constraint** on tetrahedral stereo — equivalent to omitting **`#T`** entirely (**`Undetermined`**). |
| **`!`** | **`#T`** | Atom is **not** a tetrahedral stereocenter (**`NotStereo`**). |
| **`+`** | **`#T`** | Atom **is** a tetrahedral stereocenter with an **unspecified** coset (**`Stereo(Undetermined)`**). |
| **`+`** / **`-`** (alone) | **`#c`** | **+1** / **−1** formal charge (**§7.3** above). |

Other **`#h`** / **`#a`** / **`#m`** payloads use the usual **`value-expr`** / **`decimal-tail`** rules (**§5**, **§5.3**).

| Tag | Meaning |
|-----|---------|
| **`#i`** | Isotope mass; **special** **`#i=`** (natural isotope, **§7.3**) |
| **`#c`** | Formal charge |
| **`#h`** | Implicit H count |
| **`#n`** | Lone pairs (nonbonding pair count) |
| **`#u`** | Unpaired electron count |
| **`#s`** | Spin multiplicity (2S+1) |
| **`#v`** | **Localized valence**: sum of **bond orders** of **localized** **`:bonds`** edges to **non-hydrogen** neighbors in the **molecular graph** (multiple bonds count by full order). **Excludes** implicit H (**`#h`**), **dative**, **aromatic**-section bonding, **multicenter**, and **noncovalent** contributions — those are separate fields. |
| **`#d`** | Dative **donated** pair count (electrons donated **by** this atom) |
| **`#t`** | Dative **accepted** pair count (“accepted”; electrons accepted **by** this atom) |
| **`#a`** | Aromatic π contribution; **special** **`#a*`**, **`#a+`**, **`#a!`** (**§7.3**) |
| **`#m`** | Multicenter valence; **special** **`#m*`**, **`#m+`**, **`#m!`** (**§7.3**) |
| **`#D`** | **Degree**: number of neighbors in the molecular graph (SMARTS `D`). Derived predicate evaluated against the target; **not** a ground atom field. |
| **`#X`** | **Connectivity**: degree plus implicit-H count (SMARTS `X`). Derived. |
| **`#x`** | **Ring connectivity**: number of ring bonds at the atom (SMARTS `x`). Derived. |
| **`#H`** | **Total hydrogens**: implicit H count plus explicit H neighbors (SMARTS `H`). Derived. |
| **`#R`** | **Ring membership**. **`#R<count>`** bounds the **total** ring count; **`#R(<size>)<count>`** bounds the count of rings of that **size**. Each count follows the **§5.3** omitted-numeral convention: bare **`#R`** / **`#R(<size>)`** means count **1**; **`#R<n>`** / **`#R(<size>)<n>`** means exactly **n**. **Special** **`#R*`** (no constraint), **`#R+`** (sugar for **`?r >= 1`**, "in at least one ring"), **`#R!`** (count **0**). SMARTS parity: **`R`** → **`#R+`**, **`Rn`** → **`#Rn`**, **`R0`** → **`#R!`**, **`rn`** → **`#R(n)+`**. Derived. |
| **`#T`** | **Tetrahedral stereo** configuration at the atom (SMARTS-style stereo query). Payload is a **`config`** (**§7.14**): **special** **`#T*`** / **`#T!`** / **`#T+`**, a coset literal **`#T<n>`** (e.g. **`#T1`**, **`#T2`**), or a coset operator-expression. Canonical constraint form **`{:atom [i {:tetrahedral-stereo …}]}`** (**§7.9**). |

**Case convention.** **Lowercase** tag letters denote the atom's own state fields (isotope, charge, spin, implicit H, localized valence, π contribution, …) plus the SMARTS-lowercase ring-connectivity predicate **`#x`**. **Uppercase** tag letters denote **derived predicates** (topology queries over the surrounding graph) in the SMARTS-parity set. The namespaces are **disjoint**: **`#h`** (implicit H slot) and **`#H`** (total H count) coexist, and **`#x`** (ring connectivity) and **`#X`** (connectivity) coexist, without collision. Ring membership is the single uppercase **`#R`** — total count, or **`#R(<size>)`** for one size; there is **no** lowercase **`#r`** (the former ring-size tag folds into **`#R(<size>)`**).

### 7.4 Element and bond **`order`** (via **`value-expr`**)

The **`element`** nonterminal (**atom-string** prefix) is **literal** | **wildcard `*`** | **brace set** | **negation** | **bind** | **ref** (**§7.4** grammar below). The **bond-string** **`order`** prefix (**§7.5**) is a single **`value-expr`** (**§5**), which **subsumes** literal **`nat`**, **`*`**, brace **`nat-set`**, **`?` *id***, **`?` *id* `::` *nat-set***, and **arithmetic** / logic (e.g. **`1+1`**, **`?o+1`**) where allowed by context.

```
element ::= element-literal
          | '*'
          | element-set
          | '!' element-literal
          | '!' element-set
          | element-bind
          | element-ref

element-set ::= '{' element-literal (',' element-literal)* '}'

element-bind   ::= '?' id '::' element-domain
element-domain ::= element-set
                 | '!' element-literal
                 | '!' element-set
element-ref    ::= '?' id

element-literal ::= [A-Z][a-z]*
```

- **`element-literal`**: one chemical symbol; **§7.2** (H–Og).
- **`*`**: any element; **invalid** in **Ground** unless narrowed by a containing rule outside this specification.
- **`element-set`**: finite non-empty disjunction of **one or more** **`element-literal`** entries; **§7.2**. **Query** / **Rule** when **Ground** disallows wildcards.
- **`!` *element-literal*** / **`!` *element-set***: cofinite **negation** — admits everything in the element domain **except** the named literal / set members. **§7.2** range applies to the excluded entries. **Invalid** in **Ground**.
- **`element-bind`**: **Query** / **Rule** only. Introduces a **nominal** variable **`id`** constrained to **membership in** (**`MemOp::In`**) or **exclusion from** (**`MemOp::NotIn`**) an **`element-domain`** (**§6**). **`::`** here means **set membership in a set of element symbols** (**§5**). The `!` prefix flips the operator to **`NotIn`**. **Invalid** in **Ground**.
- **`element-ref`**: **Query** / **Rule** only. **Nominal reference**: **`id`** must already be bound as a nominal in rule scope (**§6**). Appears only in the **element** position at the start of the atom-string. No arithmetic on nominal variables.

**Paren-transparency.** Outer parentheses around an **`element-bind`** or **`element-ref`** are **optional** and **semantically transparent**: implementations **MUST** accept the bare forms (**`?e`**, **`?e :: {C,N}`**) and any nesting depth of outer parens (**`(?e)`**, **`((?e :: {C,N}))`**, …) as identical AST. The **canonical** rendered form is **bare** (no outer parens).

Optional ASCII whitespace inside **`element-bind`** around **`::`** and around commas in the inner **`element-set`**, per **§7.1**.

### 7.5 Bond subgrammar

**Bond-string** uses a **separate** namespace from **atom-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on bonds (**§7.2**). The **`order`** prefix is a **`value-expr`** (**§5**); see **§7.4** for the parallel with **`element`**.

```
bond-string ::= order bond-predicate*

bond-predicate ::= '#' tag payload

order ::= value-expr
```

**Segmentation.** Let **`order-text`** be the substring of the **`bond-string`** from the first character after any leading whitespace up to (but not including) the first **`#`** that starts a **`bond-predicate`** (**`#`** immediately followed by a **tag** letter, with **no** whitespace between **`#`** and the tag), or the whole **trimmed** string if there is no such **`#`**. **`order-text`** **MUST** be parsed as **`value-expr`** (**§5**).

**Bond predicates.** **Zero or more** **`bond-predicate`** units follow **`order`**. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`a`**, **`C`**; **`#R`** **MAY** appear **multiple** times (one per ring scope — total and/or per-size). **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#a`**, **`#R`** (total first, then by size ascending), **`#C`**.

**Whitespace** between **`#`** and the **tag** letter is **invalid** (**§7.1**).

**`#c` (bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.)

**`#u`** / **`#s`.** After **`#u`** or **`#s`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **No** extra lookahead is required beyond **`value-expr`** termination and the next predicate or end of string.

**`#R` (bond ring membership).** Same forms as atom-level **`#R`** (**§7.3**): **`#R<count>`** (total ring count) or **`#R(<size>)<count>`** (count of rings of that size); bare means **1**; **`#R*`** means no constraint; **`#R+`** is sugar for **`?r >= 1`** ("bond lies in at least one ring"); **`#R!`** means count **0**.

| Tag | Meaning (bond namespace) |
|-----|---------------------------|
| **`#c`** | Bond formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (bond centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (bond centered); **`u8`** |
| **`#a`** | Bond is a member of some aromatic system declared under **`:aromatic-systems`**. Derived predicate; no payload. Canonical form is **`{:bond [i :aromatic]}`** under **`:constraints`**. |
| **`#R`** | **Ring membership**: **`#R<count>`** bounds the bond's **total** ring count, **`#R(<size>)<count>`** the count of rings of that **size**. Omitted-numeral convention (**§5.3**); **special** **`#R*`** / **`#R+`** / **`#R!`** (as for atoms, **§7.3**). Derived. |
| **`#C`** | **Cis/trans stereo** configuration at the bond (SMARTS-style stereo query). Payload is a **`config`** (**§7.14**): **special** **`#C*`** / **`#C!`** / **`#C+`**, a coset literal **`#C<n>`** (e.g. **`#C1`**, **`#C2`**), or a coset operator-expression. Canonical constraint form **`{:bond [i {:cis-trans-stereo …}]}`** (**§7.9**). |

**Bond order values** in the **bond-string** **MUST NOT** be **fractional** after evaluation (**§7.6**). **Aromatic** bond **order** as a distinct category **MUST NOT** be used in **`order`**; use **§4** instead.

### 7.6 Bond order

**Semantic model** for **localized** bond order in the **bond-string** (the **`order`** nonterminal is **`value-expr`**, **§7.5**):

- **Discrete** orders **1**, **2**, **3**, and **4** after any **arithmetic** and binding.
- **Any** order: **`*`** as **`value-expr`** (**Query** / **Rule**).
- **Finite set**: top-level **`nat-set`** in **`value-expr`** (e.g. **`{1,2,3}`** or **`{2}`**).
- **Arithmetic and constraints**: full **`value-expr`** on **`order`**, including **`add-expr`**, **`::`** **`nat-set`**, **`bool-expr`**, and **`?id`** binds, subject to **Ground** restrictions below.

In **Ground**, **`order-text`** **MUST** denote a single definite order in **{1,2,3,4}**: **`*`**, a top-level **`nat-set`** whose entries are **only** **1**–**4**, **`(?` *id* `::` *set* `)`** / **`(?` *id* `)`** only where the implementation resolves them to one value, or **`value-expr`** that is **only** **`sign`*** **`nat`** with value **1**–**4** (no **`?`**, **`::`**, relations, logic, or **`(`** … **`)`**). **Query** / **Rule** **MAY** use the full **`value-expr`** grammar on **`order`**.

This section does **not** define **`bond-keyword`** shorthands; see **§7.7**.

### 7.7 Bond and atom literals

**Bond entry shorthands.** A **`bond-keyword`** as the **`:type`** value of a **`bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent bond-string payload. Normative expansion table:

| Keyword | Expands to | Bond order |
|---------|-----------|------------|
| **`:single`** | **`"1"`** | 1 |
| **`:double`** | **`"2"`** | 2 |
| **`:triple`** | **`"3"`** | 3 |
| **`:quadruple`** | **`"4"`** | 4 |

Implementations **MUST** accept these four keywords wherever **`bond-spec`** is expected. No other **`bond-keyword`** values are defined; unrecognized keywords **MUST** be rejected.

**Dative entry shorthands.** A **`dative-keyword`** as the **`:type`** value of a **`dative-bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent dative-string payload. Normative expansion table:

| Keyword | Expands to | Pairs donated | Example |
|---------|-----------|---------------|---------|
| **`:single`** | **`"1"`** | 1 | NH₃→BF₃ |
| **`:double`** | **`"2"`** | 2 | Ni(C₄H₄)₂ |
| **`:triple`** | **`"3"`** | 3 | |
| **`:quadruple`** | **`"4"`** | 4 | uranocene U(C₈H₈)₂ |

Implementations **MUST** accept these four keywords wherever **`dative-bond-spec`** is expected. Higher pair counts and any non-trivial dative payload **MUST** use the **`dative-string`** form (**§7.12**); unrecognized **`dative-keyword`** values **MUST** be rejected.

**Atom literals.** Atom literals are **EDN strings** whose contents are **atom-string** payloads (**§7.3** / **§7.4**). Keyword-shaped atom shorthands (via **`:atom-aliases`**) are defined in **§4**.

### 7.9 Constraint grammar

Molecule-wide constraints live under the **`:constraints`** key on a **molecule-map** (**§4**). Each entry is a **single-key map** whose key names the constraint kind. Constraints fall into four categories:

- **Entity** — a value-only predicate over one entity. Same payload as the inline string form on that entity; lift/inline (below) move them between scopes.
- **Relational** — a predicate that ties one DAMN entity (dative bond, aromatic system, multicenter bond, noncovalent bond) or stereo element (stereo atom, stereo bond) to atoms, bonds, or atom-predicates. Cannot appear inline.
- **Molecule-scope** — predicates over the molecule as a whole or an arbitrary atom/bond subset.
- **Combinator** — **`:and`** / **`:or`** / **`:not`** over any constraint-entry.

```
constraint-list ::= [ constraint-entry* ]

constraint-entry ::=
    entity-constraint
  | relational-constraint
  | molecule-constraint
  | combinator-constraint

entity-constraint ::=
    { :atom              [atom-ref            atom-constraint-form] }
  | { :bond              [bond-ref            bond-constraint-form] }
  | { :dative-bond       [dative-bond-ref     dative-bond-constraint-form] }
  | { :aromatic-system   [aromatic-system-ref aromatic-system-constraint-form] }
  | { :multicenter-bond  [multicenter-bond-ref multicenter-bond-constraint-form] }
  | { :noncovalent-bond  [noncovalent-bond-ref noncovalent-bond-constraint-form] }
  | { :stereo-atom       [stereo-atom-ref      stereo-atom-constraint-form] }
  | { :stereo-bond       [stereo-bond-ref      stereo-bond-constraint-form] }

relational-constraint ::=
    { :dative-bond-donor              [dative-bond-ref atom-ref]              }
  | { :dative-bond-acceptor           [dative-bond-ref atom-ref]              }
  | { :dative-bond-parallels          [dative-bond-ref bond-ref]              }
  | { :dative-bond-donor-satisfies    [dative-bond-ref atom-constraint-form]  }
  | { :dative-bond-acceptor-satisfies [dative-bond-ref atom-constraint-form]  }
  | { :aromatic-system-atoms          [aromatic-system-ref [atom-ref+]]       }
  | { :aromatic-system-contains       [aromatic-system-ref atom-ref]          }
  | { :aromatic-system-contains-all   [aromatic-system-ref [atom-ref+]]       }
  | { :aromatic-system-all-atoms      [aromatic-system-ref atom-constraint-form] }
  | { :aromatic-system-any-atom       [aromatic-system-ref atom-constraint-form] }
  | { :multicenter-bond-atoms         [multicenter-bond-ref [atom-ref+]]      }
  | { :multicenter-bond-contains      [multicenter-bond-ref atom-ref]         }
  | { :multicenter-bond-contains-all  [multicenter-bond-ref [atom-ref+]]      }
  | { :multicenter-bond-all-atoms     [multicenter-bond-ref atom-constraint-form] }
  | { :multicenter-bond-any-atom      [multicenter-bond-ref atom-constraint-form] }
  | { :noncovalent-bond-ends          [noncovalent-bond-ref [atom-ref atom-ref]] }
  | { :noncovalent-bond-contains      [noncovalent-bond-ref atom-ref]         }
  | { :noncovalent-bond-ends-satisfy  [noncovalent-bond-ref
                                       [atom-constraint-form atom-constraint-form]] }
  | { :stereo-atom-site               [stereo-atom-ref atom-ref]              }
  | { :stereo-atom-contains           [stereo-atom-ref atom-ref]              }
  | { :stereo-atom-ligands            [stereo-atom-ref [atom-ref+]]           }
  | { :stereo-atom-all-ligands        [stereo-atom-ref atom-constraint-form]  }
  | { :stereo-atom-any-ligand         [stereo-atom-ref atom-constraint-form]  }
  | { :stereo-bond-site               [stereo-bond-ref bond-ref]              }
  | { :stereo-bond-contains           [stereo-bond-ref atom-ref]              }
  | { :stereo-bond-ligands            [stereo-bond-ref [atom-ref+]]           }
  | { :stereo-bond-all-ligands        [stereo-bond-ref atom-constraint-form]  }
  | { :stereo-bond-any-ligand         [stereo-bond-ref atom-constraint-form]  }

molecule-constraint ::=
    { :charge-sum     { [:atoms [atom-ref+]]? :sum value-expr } }
  | { :spin-sum       { [:atoms [atom-ref+]]? :spin spin-literal } }
  | { :bond-order-sum { [:bonds [bond-ref+]]? :sum value-expr } }
  | { :connected      { [:atoms [atom-ref+]]? } }
  | { :sub-pattern    { :anchor anchor-spec :pattern molecule-map } }

combinator-constraint ::=
    { :and [constraint-entry+] }
  | { :or  [constraint-entry+] }
  | { :not constraint-entry }

atom-constraint-form ::=
    { :valence             value-expr }
  | { :aromatic-valence    ( value-expr | :not-aromatic | :undetermined ) }
  | { :multicenter-valence ( value-expr | :not-multicenter | :undetermined ) }
  | { :donated-pairs       value-expr }
  | { :accepted-pairs      value-expr }
  | { :degree              value-expr }
  | { :connectivity        value-expr }
  | { :ring-connectivity   value-expr }
  | { :total-hydrogens     value-expr }
  | { :ring-membership     ring-membership-form }
  | { :tetrahedral-stereo  stereo-config-form }

bond-constraint-form ::=
    :aromatic
  | { :ring-membership ring-membership-form }
  | { :cis-trans-stereo stereo-config-form }

dative-bond-constraint-form ::=
    :aromatic
  | { :ring-membership ring-membership-form }

(* One ring-membership fact. With :size, a count of rings of that size;     *)
(* without :size, the total ring count. :count is required.                 *)
ring-membership-form ::= { :size nat :count value-expr } | { :count value-expr }

aromatic-system-constraint-form  ::= { :electron-count value-expr }
multicenter-bond-constraint-form ::= { :electron-count value-expr }
noncovalent-bond-constraint-form ::= (* uninhabited — no value-only variants yet *)

(* A stereo entity-constraint form carries the element's :kind (its stereo *)
(* subtype) plus exactly one predicate key. :kind is mandatory — a detached *)
(* molecule-scope constraint is not well-formed without the subtype, since *)
(* a permutation cannot recover its degree (the kind is many-to-one on degree). *)
stereo-atom-constraint-form ::= { :kind stereo-kind stereo-predicate-pair }
stereo-bond-constraint-form ::= { :kind stereo-kind stereo-predicate-pair }

stereo-predicate-pair ::=
    :ligand-symmetry ligand-symmetry-form
  | :fluxionality    permutation-form
  | :topicity        { :pair ligand-pair :relation topicity-relation }
  | :stereogenicity  stereogenicity-relation

stereo-kind          ::= :tetrahedral | :cis-trans | :axial | :square-planar
                       | :trigonal-bipyramidal | :octahedral
permutation-form     ::= [ cycle* ]                 (* vector of disjoint cycles; identity [] *)
cycle                ::= [ nat+ ]                    (* p0→p1→…→p0, 0-indexed positions *)
ligand-symmetry-form ::= { :perm permutation-form [:orientation (:proper | :improper)]
                                                   [:member (:in | :not-in)] }
ligand-pair          ::= [ nat nat ]                (* two ligand-frame positions *)
topicity-relation       ::= :homotopic | :enantiotopic | :diastereotopic | #{ keyword+ } | :undetermined
stereogenicity-relation ::= :symmetric | :prochiral | :stereogenic       | #{ keyword+ } | :undetermined

stereo-config-form ::= :undetermined | :not-stereo | { :stereo coset-form }
coset-form ::= int | :undetermined | [ int+ ] | "coset-string"

anchor-spec ::=
    { [:atoms             [[atom-ref atom-ref]+]]?
      [:bonds             [[bond-ref bond-ref]+]]?
      [:dative-bonds      [[dative-bond-ref dative-bond-ref]+]]?
      [:aromatic-systems  [[aromatic-system-ref aromatic-system-ref]+]]?
      [:multicenter-bonds [[multicenter-bond-ref multicenter-bond-ref]+]]?
      [:noncovalent-bonds [[noncovalent-bond-ref noncovalent-bond-ref]+]]?
      [:stereo-atoms      [[stereo-atom-ref stereo-atom-ref]+]]?
      [:stereo-bonds      [[stereo-bond-ref stereo-bond-ref]+]]? }

atom-ref             ::= int | keyword
bond-ref             ::= int | keyword
dative-bond-ref      ::= int | keyword
aromatic-system-ref  ::= int | keyword
multicenter-bond-ref ::= int | keyword
noncovalent-bond-ref ::= int | keyword
stereo-atom-ref      ::= int | keyword
stereo-bond-ref      ::= int | keyword
```

**Ref resolution.** An integer ref is the **positional** index into the corresponding entity vector on the molecule map: **`atom-ref`** → **`:atoms`**, **`bond-ref`** → **`:bonds`**, **`dative-bond-ref`** → **`:dative-bonds`**, **`aromatic-system-ref`** → **`:aromatic-systems`**, **`multicenter-bond-ref`** → **`:multicenter-bonds`**, **`noncovalent-bond-ref`** → **`:noncovalent-bonds`**, **`stereo-atom-ref`** → **`:stereo-atoms`**, **`stereo-bond-ref`** → **`:stereo-bonds`**. A keyword ref resolves against the **`:id`** declared on the corresponding entry (**§4**). On serialization, implementations **MUST** emit the **`:id`** keyword when one is declared on the referenced entry, falling back to the positional integer otherwise.

**Narrow inner forms for DAMN entities.** **`:aromatic-system`** and **`:multicenter-bond`** narrow leaves carry only the **`:electron-count`** value-only variant; every other predicate on those entities is a relational leaf instead. **`:noncovalent-bond`** narrow leaves have no inhabited inner form yet; every noncovalent predicate is a relational leaf.

**Stereo entity constraints carry the kind.** The **`:stereo-atom`** / **`:stereo-bond`** entity-constraint forms (**`#p`** / **`#f`** / **`#o`** / **`#g`**) carry a mandatory **`:kind`** naming the element's stereo subtype, plus exactly one predicate key. The **`:kind`** is redundant with the referenced element at the **entity** level (the inline form omits it — the **`:type`** **`class`** supplies it, **§7.14**) but is **required** at molecule scope, where the constraint is detached from its element: a permutation payload cannot recover its degree, and **`stereo-kind`** is many-to-one on degree (**`:tetrahedral`** and **`:square-planar`** are both degree 4). The kind/degree (and the chiral-class restriction on **`'`** values) is cross-checked against the resolved element by the validator (**§6.1**); **`inline_constraints`** drops the carried kind back into the element. These constraints are **distinct** from the atom/bond **`:tetrahedral-stereo`** (**`#T`**) / **`:cis-trans-stereo`** (**`#C`**) inline configurations (which assert the local **coset** at the bearing atom/bond) and from the stereo **relational leaves** (**`:stereo-atom-…`** / **`:stereo-bond-…`**, no inline form).

**Anchor cardinality.** Each keyed slot in **`anchor-spec`** is optional and may appear at most once; if present, it is a vector of **`(target-side-ref, pattern-side-ref)`** pairs of the same entity kind. An empty **`anchor-spec`** denotes an unanchored sub-pattern (the pattern can embed anywhere). Target-side refs resolve against the outer molecule's metadata; pattern-side refs against the pattern molecule's metadata.

**Molecule-scope subset selectors.** `:charge-sum`, `:spin-sum`, `:bond-order-sum`, and `:connected` accept an **optional** `:atoms` (or `:bonds`) vector. When **omitted**, the predicate ranges over **every** atom (or bond) in the molecule, including atoms added by future structural growth. When **present**, the predicate ranges over the listed entities only. An empty vector `[]` is **distinct** from omission: it selects no entities.

**Sub-pattern materialization.** A **`:sub-pattern`** **`:pattern`** is a full **molecule-map**; its inner constraints are evaluated independently from the outer constraint tree. The pattern carries **no defaults** — values pass through verbatim — so a pattern's atom **`charge: undetermined`** stays **`undetermined`** at match time and behaves as a wildcard (**§5.4**); zero-defaulting that would apply to a ground input does **not** apply inside a pattern.

**Sugar (inline string equivalents).** Narrow leaves whose entity has a string subgrammar admit two interchangeable serializations:

- the inline **`#tag`** payload on the entity's **`:type`** string (or, for atoms, on the atom literal directly);
- the **`{:<entity> [ref form]}`** entry in **`:constraints`**.

Parsers **MUST** accept both. Bare per-entity predicates (not nested under **`:and`** / **`:or`** / **`:not`** / **`:sub-pattern`**) **MAY** be emitted in the sugared inline form; nested predicates **MUST** be emitted as **`:constraints`** entries since the inline form has no logical context.

**Inline-form coverage by entity:**

- **Atom** (**§7.3**): all `atom-constraint-form` variants except the derived ones (`#D`, `#X`, `#x`, `#H`, `#R`) lift to inline atom predicates, including `:tetrahedral-stereo` → `#T`; the derived predicates also have inline tags but are pattern-only.
- **Bond** (**§7.5**): all `bond-constraint-form` variants (`:aromatic`, `:ring-membership`, `:cis-trans-stereo`) have inline forms (`#a`, `#R`, `#C`).
- **Dative bond** (**§7.12**): all `dative-bond-constraint-form` variants (`:aromatic`, `:ring-membership`) have inline forms (`#a`, `#R`).
- **Aromatic system** (**§7.10**): the single `aromatic-system-constraint-form` variant `:electron-count` has the inline form `#e<n>`.
- **Multicenter bond** (**§7.11**): the single `multicenter-bond-constraint-form` variant `:electron-count` has the inline form `#e<n>`.
- **Noncovalent bond** (**§7.13**): `noncovalent-bond-constraint-form` is uninhabited; no inline-form question arises.
- **Stereo atom / stereo bond** (**§7.14**): all four `stereo-atom-constraint-form` / `stereo-bond-constraint-form` predicates (`:ligand-symmetry`, `:fluxionality`, `:topicity`, `:stereogenicity`) have inline forms (`#p`, `#f`, `#o`, `#g` on the `:type` string). On **inline** the kind is omitted (the `:type` `class` supplies it); on **lift** the element's kind is written into the molecule-scope form's `:kind`. The atom/bond `:tetrahedral-stereo` / `:cis-trans-stereo` predicates (`#T` / `#C`) are separate atom/bond inline constraints; the `:stereo-atom-…` / `:stereo-bond-…` predicates are relational leaves with no inline form.

**Relational leaves** (**§7.9** `relational-constraint`) and **molecule-scope leaves** (`molecule-constraint`) have **no** inline form regardless of which entity they reference.

**Lift / inline.** The two storage scopes — inline on the entity (`AtomAst::constraints` etc.) and at molecule scope (`MoleculeAst::constraints` as `{:atom [ref form]}` peers) — are interchangeable for the inline-capable narrow leaves. Implementations **SHOULD** expose:

- **`lift_constraints`** (entity → molecule): drains every inline store into the molecule list as `{:<entity> [ref form]}` peers.
- **`inline_constraints`** (molecule → entity): drains top-level inline-capable narrow leaves from the molecule list into the targeted entity's inline store.

Combinator subtrees, relational leaves, and molecule-scope leaves are never moved by either operation. With multiple top-level entries targeting the same (entity, kind), `inline_constraints` resolves the collision via the entity store's per-kind insert policy (last-wins for unique-kind variants).

**Multiple constraints per entity.** Each per-entity constraint serializes as its **own** entity-constraint entry; implementations **MUST NOT** bundle multiple constraints on the same entity into a single map.

### 7.10 Aromatic system subgrammar

**Aromatic-string** uses a **separate** namespace from **atom-string** and **bond-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on aromatic systems (**§7.2**). It carries **per-aromatic-system** state — overall **charge** (**`#c`**) and **spin** (**`#u`**, **`#s`**) as inherent fields, and an **optional asserted total π-electron count** (**`#e<n>`**) as an inline constraint — as the **`:type`** value of an **`aromatic-system-entry`** (**§4**). The **per-atom** π contributions live in the entry's **`:electrons`** vector (**§4**), not in this string.

```
aromatic-string ::= aromatic-predicate*

aromatic-predicate ::= '#' tag payload
```

**Aromatic predicates.** **Zero or more** **`aromatic-predicate`** units. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (aromatic-system formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.) Same convention as atom (**§7.3**) and bond (**§7.5**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **`#e`** omitted means **1** π-electron.

| Tag | Meaning (aromatic-system namespace) | Storage |
|-----|---------------------------------------|----------|
| **`#c`** | Aromatic-system formal charge (**`i8`**, **§7.2**) | inherent field |
| **`#u`** | Unpaired electrons (system-centered); **`u8`** | inherent field |
| **`#s`** | Spin multiplicity (2S+1) (system-centered); **`u8`** | inherent field |
| **`#e`** | Asserted total π-electron count; **`u8`** | inline constraint (`AromaticSystemConstraint::ElectronCount`) |

**`#e<n>` semantics.** **`#e<n>`** asserts the system's total π-electron count and parses to an inline aromatic-system constraint (`AromaticSystemConstraint::ElectronCount(n)`) on the entry's constraint store, **not** to a direct field. The per-atom contributions in the entry's **`:electrons`** vector (**§4**) are the canonical data; **`#e<n>`** is the optional total assertion that downstream validation cross-checks against **`sum(:electrons)`** on ground inputs. **`#e`** is omitted from the canonical entity-string form when no `ElectronCount` constraint is present.

**No canonical-constraint equivalent for charge / spin.** Aromatic-system charge (**`#c`**), unpaired electrons (**`#u`**), and spin multiplicity (**`#s`**) live as direct fields on the aromatic-system entity (set by the aromatic-string predicates above) and have **no** canonical **`:constraints`** form.

### 7.11 Multicenter-bond subgrammar

**Multicenter-string** uses a **separate** namespace from **atom-string**, **bond-string**, and **aromatic-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on multicenter bonds (**§7.2**). It carries **per-multicenter-bond** state — overall **charge** (**`#c`**) and **spin** (**`#u`**, **`#s`**) as inherent fields, and an **optional asserted total electron count** (**`#e<n>`**) as an inline constraint — as the **`:type`** value of a **`multicenter-bond-entry`** (**§4**). The **per-atom** electron contributions live in the entry's **`:electrons`** vector (**§4**), not in this string.

```
multicenter-string ::= multicenter-predicate*

multicenter-predicate ::= '#' tag payload
```

**Multicenter predicates.** **Zero or more** **`multicenter-predicate`** units. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (multicenter-bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. Same convention as atom (**§7.3**), bond (**§7.5**), and aromatic (**§7.10**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **`#e`** omitted means **1** bonded electron.

| Tag | Meaning (multicenter-bond namespace) | Storage |
|-----|----------------------------------------|----------|
| **`#c`** | Multicenter-bond formal charge (**`i8`**, **§7.2**) | inherent field |
| **`#u`** | Unpaired electrons (bond-centered); **`u8`** | inherent field |
| **`#s`** | Spin multiplicity (2S+1) (bond-centered); **`u8`** | inherent field |
| **`#e`** | Asserted total bonded electron count; **`u8`** | inline constraint (`MulticenterBondConstraint::ElectronCount`) |

**`#e<n>` semantics.** **`#e<n>`** asserts the multicenter bond's total electron count and parses to an inline multicenter-bond constraint (`MulticenterBondConstraint::ElectronCount(n)`), parallel to the aromatic-system case (**§7.10**). Per-atom contributions in the entry's **`:electrons`** vector (**§4**) are the canonical data; **`#e<n>`** is the optional total assertion that downstream validation cross-checks against **`sum(:electrons)`** on ground inputs.

**Per-atom participation.** The atom-side **`#m`** predicate (**§7.3**) is a per-atom multicenter-membership marker; the per-atom electron share for a given multicenter bond lives in that bond's **`:electrons`** vector (**§4**), not inside the multicenter-string. Endpoint references (which atoms the bond spans) live in the **`:atoms`** vector of the **`multicenter-bond-entry`** (**§4**); they **MUST NOT** be encoded inside the multicenter-string.

### 7.12 Dative-bond subgrammar

**Dative-string** carries the **bond order** (number of donated electron pairs) and optional **aromatic** and **ring-membership** constraints on a single **`dative-bond-entry`** (**§4**). The grammar parallels **bond-string** (**§7.5**): a leading **`order`** token followed by zero or more **`#…`** predicates. The dative-string has **no** inherent-field tags beyond order and **no** direction token; direction is expressed entirely by the **`:donor`** / **`:acceptor`** assignment on the containing entry.

```
dative-string ::= order dative-predicate*

order            ::= value-expr | '*'
dative-predicate ::= '#' tag payload
```

**Order.** The leading **`order`** token is a **`value-expr`** (**§5**) — typically a positive integer literal — that records how many electron pairs are donated. **`*`** means **`Undetermined`**. The **`dative-keyword`** shorthands (**§7.7**) — **`:single`**, **`:double`**, **`:triple`**, **`:quadruple`** — expand to the literal forms **`"1"`**, **`"2"`**, **`"3"`**, **`"4"`**.

**Dative predicates.** **Zero or more** **`dative-predicate`** units after the order token. **Optional** ASCII whitespace **MAY** appear between the order and the first **`#`**, and between successive predicates. **At most one** **`#a`**; **`#R`** **MAY** appear **multiple** times (one per ring scope — total and/or per-size). **Canonical** predicate order (stable serialization): order, then **`#a`**, then **`#R`** (total first, then by size ascending).

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#a` (dative-bond aromatic flag).** A bare **`#a`** with no payload marks the dative bond as participating in an aromatic system. Examples: the N→B π-donation of borazine, O→B of boroxine, or a C→M coordination spanning a metallaaromatic ring. The flag carries no value (parallel to the bond-namespace **`#a`** of **§7.5**); aromatic-ring perception cross-checks the flag against actual ring membership.

**`#R` (dative-bond ring membership).** Same forms as the atom-level and bond-level **`#R`** (**§7.3**, **§7.5**): **`#R<count>`** (total ring count) or **`#R(<size>)<count>`** (count of rings of that size); bare means **1**; **`#R*`** means no constraint; **`#R+`** is sugar for **`?r >= 1`** ("dative bond lies in at least one ring"); **`#R!`** means count **0**.

| Tag | Meaning (dative-bond namespace) | Storage |
|-----|-----------------------------------|----------|
| (leading) | **Order**: number of donated electron pairs (**`u8`**, **§7.2**) | inherent field |
| **`#a`** | **Aromatic**: the dative bond is part of an aromatic system. | derived |
| **`#R`** | **Ring membership**: **`#R<count>`** (total) / **`#R(<size>)<count>`** (per size); **special** **`#R*`** / **`#R+`** / **`#R!`** (**§7.3**). | derived |

**Direction.** Dative bonds are intrinsically directional. Direction is carried entirely by the ordered **`:donor`** / **`:acceptor`** assignment on the containing **`dative-bond-entry`** (**§4**); the dative-string itself has **no** direction token. Under pattern matching (**§6**), the embedding MUST map pattern **`:donor`** atoms to target donors and the pattern **`:acceptor`** to the target acceptor — a donor/acceptor swap across the embedding rejects the match.

**Donor / acceptor / cross-bond references.** Donor-side and acceptor-side constraints on the endpoint atoms (equivalent to the **`:donated-pairs`** / **`:accepted-pairs`** atom-constraint forms of **§7.9**, or atom-string **`#d`** / **`#t`** of **§7.3** pinned to one endpoint) attach via the molecule-wide **`:constraints`** section; they **MUST NOT** be encoded inside the dative-string. The same holds for the "parallels another bond" relation and any reference to other molecule-level entities.

### 7.13 Noncovalent-bond subgrammar

**Noncovalent-string** encodes the **interaction kind** of a single **`noncovalent-bond-entry`** (**§4**) as an expression. It has **no** predicates; the kind is the whole payload.

```
noncovalent-string ::= noncovalent-kind-expr

noncovalent-kind-expr ::= noncovalent-kind-literal
                        | '*'
                        | noncovalent-kind-set
                        | noncovalent-kind-bind
                        | noncovalent-kind-ref

noncovalent-kind-set  ::= '{' noncovalent-kind-literal
                              (',' noncovalent-kind-literal)* '}'
noncovalent-kind-bind ::= '(' '?' id '::' noncovalent-kind-set ')'
noncovalent-kind-ref  ::= '(' '?' id ')'

noncovalent-kind-literal ::= 'Hbd' | 'Xbd' | 'Ybd' | 'Ion' | 'Vdw'
```

**Literal meanings.**

| Literal | Interaction kind | Keyword equivalent (**§4**) |
|---------|------------------|-------------------------------|
| **`Hbd`** | hydrogen bond | **`:h-bond`** |
| **`Xbd`** | halogen bond | **`:halogen-bond`** |
| **`Ybd`** | chalcogen bond | **`:chalcogen-bond`** |
| **`Ion`** | ionic interaction | **`:ionic`** |
| **`Vdw`** | van der Waals interaction | **`:van-der-waals`** |

Each **`noncovalent-kind-literal`** is exactly three ASCII characters: one leading uppercase letter followed by two lowercase letters. The parser consumes the full three-character token; partial prefixes (**`H`**, **`Hb`**, …) **MUST** be rejected. Leading / trailing whitespace on the whole **noncovalent-string** is ignored (**§7.1**).

**Wildcard and sets.** **`*`** admits any kind. An **`noncovalent-kind-set`** **`{Hbd,Ion}`** admits its members. These forms are **invalid** in **Ground**; **Query** / **Rule** **MAY** use them.

**Bind and ref.** **`noncovalent-kind-bind`** introduces a **nominal** variable **`id`** constrained to membership in the given set (**§6**); **`noncovalent-kind-ref`** references a nominal binding established elsewhere in the rule scope. Both are **invalid** in **Ground**. **`::`** here means **set membership in a set of noncovalent-kind symbols**, parallel to its use in **`element-bind`** (**§7.4**, **§5**).

**Relation to `noncovalent-keyword` (`§4`).** The five **`noncovalent-keyword`** values (**`:h-bond`** etc.) are a ground-only EDN shorthand on the **`noncovalent-bond-entry`**'s **`:type`**; they expand to the corresponding **`noncovalent-kind-literal`** and are semantically identical. The keyword shorthand **MUST NOT** be used inside a **`noncovalent-kind-set`** or a **`noncovalent-kind-bind`**; those take kind literals only.

### 7.14 Stereo subgrammar

**Stereo-string** uses a **separate** namespace from the atom / bond / aromatic / multicenter strings (**§7.2**). It is the **`:type`** payload of a **`stereo-atom-entry`** / **`stereo-bond-entry`** (**§4**) and names the coordination **`class`** plus the realized **`coset`** index over the entry's ordered **`:ligands`** frame. The **same** **`config`** grammar — **`coset`** preceded by two extra sentinels and **without** the leading **`class`** — is the payload of the atom **`#T`** / bond **`#C`** inline constraints (**§7.3** / **§7.5**) and the source form behind the **`{:stereo coset-form}`** EDN constraint (**§7.9**).

```
stereo-string ::= class coset stereo-predicate*
config        ::= '*' | '!' | '+' | nat | coset-expr

class ::= 'Th' | 'Ct' | 'Sp' | 'Tb' | 'Oh'
coset ::= '*' | nat | coset-expr

coset-expr ::= '~' coset-expr             (* involution operator             *)
             | '\'' coset-expr            (* mirror operator                 *)
             | coset-expr '^' cycles      (* group action by a permutation   *)
             | nat                        (* literal coset index             *)
             | '?' id [ '::' coset-set ]  (* coset variable / domain         *)
             | coset-set                  (* literal coset set               *)

coset-set ::= '{' nat (',' nat)* '}'

stereo-predicate ::=
    '#p' ['!'] ( '~' | ['\''] cycles )    (* ligand symmetry: ! not-in, ' improper, ~ kind involution *)
  | '#f' ( '~' | cycles )                 (* fluxionality (proper move); ~ kind involution            *)
  | '#o' relation ligand-pair             (* topicity of a ligand pair                                *)
  | '#g' relation                         (* stereogenicity classification                            *)

cycles      ::= '()' | ( '(' nat (',' nat)* ')' )+   (* disjoint cycles, 0-indexed, identity ()       *)
ligand-pair ::= '(' nat ',' nat ')'                  (* two 0-indexed ligand-frame positions          *)
relation    ::= '*' | ['!'] glyph                    (* * undetermined; ! set-complement; glyph singleton *)
glyph       ::= '=' | '\'' | '/'
```

**Class.** **`Th`** tetrahedral, **`Ct`** cis/trans, **`Sp`** square-planar, **`Tb`** trigonal-bipyramidal, **`Oh`** octahedral. A **`stereo-atom-entry`** carries an atom-centered class (**`Th`** / **`Sp`** / **`Tb`** / **`Oh`**); a **`stereo-bond-entry`** carries **`Ct`**. Matching presently realizes **`Th`** and **`Ct`**; **`Sp`** / **`Tb`** / **`Oh`** parse and round-trip but their matching is **staged**.

**Coset.** The **`coset`** is a **dense per-class arrangement index** over the entry's ordered **`:ligands`** frame (**§4**) — the OpenSMILES arrangement number for the class, **not** a Lehmer / permutation rank. For **`Th`**: **`1`** = anticlockwise (**`@`**), **`2`** = clockwise (**`@@`**). For **`Ct`**: **`1`** = **Z** (cis), **`2`** = **E** (trans). **`Sp`** / **`Tb`** / **`Oh`** follow the OpenSMILES numbering for their class. **`*`** is an **undetermined** (open) coset.

**`config`** (atom **`#T`** / bond **`#C`** / **`{:stereo …}`**).** The constraint payload extends **`coset`** with two leading sentinels: **`*`** = **`Undetermined`** (no stereo constraint — equivalent to omitting the predicate), **`!`** = **`NotStereo`** (the site is **not** a stereocenter), **`+`** = **`Stereo`** with an **undetermined** coset (the site **is** a stereocenter, coset unspecified). A bare **`nat`** / **`coset-expr`** is **`Stereo`** with that coset. The EDN equivalents are **`:undetermined`** / **`:not-stereo`** / **`{:stereo :undetermined}`** / **`{:stereo coset-form}`** (**§7.9**); a **`coset-set`** serializes to the EDN vector form **`[ int+ ]`**, every other **`coset-expr`** to a **`"coset-string"`**.

**Coset operators (reserved).** The **`~`** (involution) and **`^`*cycles*** (group action by a permutation in 0-indexed disjoint-cycle notation, **§7.14**) operators, and the coset variable / set forms (**`?id`**, **`?id :: {…}`**, **`{…}`**), **parse** as **`coset-expr`** and **round-trip**, but their **matching** semantics are **staged** — relative-stereo binding and non-tetrahedral coset domains land incrementally. Only **ground literal cosets** (and the **`*`** / **`!`** / **`+`** sentinels) are presently matched; a conforming matcher **MAY** reject an operator / variable **`coset-expr`** until the corresponding stage lands.

**Inline predicates (`#p` / `#f` / `#o` / `#g`).** After the leading **`class coset`**, a **`stereo-string`** carries **zero or more** stereo predicates — the inline form of the per-element stereo constraints (the molecule-scope structured peers are **§7.9**). Each predicate's permutation degree is the **`class`** degree (number of ligand positions). The four predicates are:

- **`#p`** — **ligand symmetry**: asserts a ligand permutation is (or, with **`!`**, is not) a symmetry of the element. Payload **`['!'] (['`''] cycles | '~')`**: **`!`** asserts **non**-membership (default: membership); **`'`** marks the permutation **improper** (orientation-reversing; default proper); **`cycles`** is the permutation in **disjoint-cycle notation** (below). The **`~`** sugar denotes the **class involution** (the orientation-reversing generator for chiral classes, the configuration-swapping ligand permutation for achiral classes) and already carries the class-appropriate orientation, so it is not combined with **`'`**.
- **`#f`** — **fluxionality**: a proper ligand permutation realized by dynamics. Payload **`(cycles | '~')`**; no **`!`**/**`'`** (it is a bare proper move). **`~`** is the class involution.
- **`#o`** — **topicity** of a ligand pair. Payload **`relation ligand-pair`** where **`ligand-pair`** is **`(i,j)`** (two 0-indexed positions in the **`:ligands`** frame). The glyphs are **`=`** homotopic, **`'`** enantiotopic, **`/`** diastereotopic.
- **`#g`** — **stereogenicity** classification. Payload **`relation`**: **`=`** symmetric, **`'`** prochiral, **`/`** stereogenic.

**Disjoint-cycle notation.** **`cycles`** is a product of disjoint cycles **`(p0,p1,…)(q0,q1,…)`**, **0-indexed** over the ligand frame, each cycle mapping **`p0→p1→…→p0`**; the identity is **`()`**. Cycle points **MUST** be in range (**`< class degree`**) and disjoint. This is the same permutation the **`Permutation`** `Display` emits and the structured **`permutation-form`** (**§7.9**) encodes as a vector of cycles.

**Relation completeness.** **`relation`** is **`*`** (**`Undetermined`** — the full domain), a single glyph (a **singleton**), or **`!`** + a glyph (the 2-element **complement** of that glyph's value). Over the 3-element topicity / stereogenicity domain these cover every non-empty subset. On serialization, a singleton renders as its glyph and a 2-element set as **`!`** + the complement glyph. A full-domain (**`Undetermined`**) relation is a **vacuous** predicate: like the atom **`#a*`** / **`#T*`** special forms (**§7.3**) it is **admissible on parse** but **elided** from the canonical rendered string (**§7.1**) — **`#o*`** / **`#g*`** are dropped on render, equivalent to omitting the predicate.

**`~` rendering.** A **`#p`** / **`#f`** permutation equal to the class involution (and, for **`#p`**, matching the involution's orientation) renders as **`~`**; otherwise as explicit **`cycles`**.

**Chiral-class restriction.** The **`'`** value — **`#p`** improper, **`#o`** enantiotopic, **`#g`** prochiral — is meaningful only on a **chiral class** (**`Th`** atom; **`Ct`** / **`Sp`** are achiral). It **parses** on any class; an inconsistent class/value pairing is rejected by the **validator** (the resolved-symmetry cross-check, **§6.1**), not at parse.

**`stereo-keyword` shorthand (`§4`).** The four **`stereo-keyword`** values expand to canonical **`class`**+**`coset`** literals: **`:ccw`** → **`Th0`**, **`:cw`** → **`Th1`**, **`:z`** → **`Ct0`**, **`:e`** → **`Ct1`**. They are a ground EDN shorthand on the **`stereo-spec`**'s **`:type`** and are semantically identical to the expanded string. On serialization, implementations **MUST** emit the **`stereo-keyword`** for these four canonical shapes **only when the element carries no inline predicates**, falling back to the **`stereo-string`** otherwise.

---

## 8. Molecule map examples (non-normative)

Examples use the vector **`:atoms`** form with inline ids. Bond entries show **`:id`** where useful.

### 8.1 Methanol (CH₃OH) — Ground, L1

```clojure
{:atoms [[:C  "C#h3"]
         [:O  "O#h1"]
         [:H  "H"]]
 :bonds [[:C :O :single]
         [:O :H :single]]}
```

The **`H`** atom here represents an **explicit** hydrogen (e.g. a hydroxyl H one wishes to name). Implicit H counts on **`C`** (**`#h3`**) and **`O`** (**`#h1`**) already account for the remaining hydrogens.

### 8.2 Indole — Ground, L1, aromatic ring

```clojure
{:atoms [[:N   "N#h1"]
         [:C2  "C"]
         [:C3  "C"]
         [:C3a "C"]
         [:C4  "C#h1"]
         [:C5  "C#h1"]
         [:C6  "C#h1"]
         [:C7  "C#h1"]
         [:C7a "C"]]
 :bonds [[:N :C2 :single] [:C2 :C3 :double] [:C3 :C3a :single]
         [:C3a :C7a :single] [:C7a :N :single] [:C3a :C4 :single]
         [:C4 :C5 :double] [:C5 :C6 :single] [:C6 :C7 :double]
         [:C7 :C7a :single]]
 :aromatic-systems [{:id :ar1 :atoms [:N :C2 :C3 :C3a :C7a] :type ""}
            {:id :ar2 :atoms [:C3a :C4 :C5 :C6 :C7 :C7a] :type ""}]}
```

Localized bonds carry the σ-skeleton orders; the aromatic π system is expressed in **`:aromatic-systems`**.

### 8.3 Substructure query — L2

Match any carbon with at least two implicit hydrogens that is directly bonded to a nitrogen:

```clojure
{:atoms [[:C "C#h(?h >= 2)"]
         [:N "N"]]
 :bonds [[:C :N :single]]}
```

**`(?h >= 2)`** is a **`bool-expr`** payload on **`#h`**; **`?h`** is bound to the matched atom's implicit H count.

### 8.4 Transformation rule — L3

Replace a primary amine carbon (C with three H and bonded to NH₂) with a quaternary carbon (no H, same bond to nitrogen):

```clojure
{:lhs {:atoms [[:C "C#h3"]
               [:N "N#h2"]]
       :bonds [[:C :N :single]]}
 :rhs {:atoms [[:C "C#h0"]
               [:N "N#h2"]]
       :bonds [[:C :N :single]]}}
```

(The **`:lhs`** / **`:rhs`** wrapping is a rule-level convention, not a molecule map key — not normative here.)
