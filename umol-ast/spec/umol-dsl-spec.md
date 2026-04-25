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

- **L1 — Ground**: fully instantiated molecules; no wildcards, binds, or logic. Every numeric slot is a concrete integer; element is a single symbol. The atom-string and bond-string parsers in the current reference implementation target this level.
- **L2 — Constraint / Query**: adds wildcards (**`*`**), element/numeric sets, **`bool-expr`**, and **`?id`** references. Sufficient for substructure queries.
- **L3 — Rule**: adds **`element-bind`** / **`element-ref`**, cross-atom **`id`** scope (**§6**), and **`:guards`** on molecule maps. Sufficient for transformation rules.
- **L4 — Compound**: rules whose RHS produces molecule maps that are themselves L2/L3 terms; higher-order composition. Not further specified in this revision.

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

**Implementation note (non-normative).** **Ground** input uses the **same** grammars as **Query** / **Rule**; only allowed token shapes differ. For speed, implementations **MAY** provide a **Ground-only** parser (or fast path) and a **full** parser, analogous to **basic** vs **extended** entry points in **`umol-models-graph`** I/O (e.g. **`umol-models-graph/src/io/smiles/parser.rs`**, **`umol-models-graph/src/io/ctfile/parser.rs`**). That split is **not** a second language: any string valid as **Ground** **MUST** be interpreted the same whether handled by the restricted implementation or the full one.

**Open issue (blocks full Query / Rule lowering).** In the current **reference model**, a **Ground** atom is **fully specified** in every numeric slot (**isotope mass** on the resolved atom remains **optional** as today). **Query** and **Rule** introduce **partial** or **wildcard** attribute states. Allowing **indeterminacy** in the **electron budget** (e.g. lone pairs, unpaired electrons, and related fields not all fixed independently) has **not** been specified: how many **valid completions** exist for a given pattern, and how **`#n`** (lone pairs) vs **`#u`** / **`#s`** (unpaired electrons vs spin multiplicity) interact under partial data, are **TBD**. Until that is defined, implementations **MAY** reject or restrict constraint atom-strings beyond **literal Ground**-shaped payloads.

---

## 4. Molecule map

A **molecule map** has the following **normative** keys. Optional keys **MAY** be absent; absent means “not applicable” for structural sections. **`nil`** on optional value-typed keys means “unknown” in builder-oriented evaluation and “unconstrained” in query-oriented evaluation, when implementations distinguish those modes.

```
molecule-map ::=
  { :atoms         atom-collection
    :bonds         localized-bond-list
    [:dative       dative-bond-list]
    [:aromatic     aromatic-list]
    [:multicenter  multicenter-list]
    [:noncovalent  noncovalent-list]
    [:atom-aliases alias-list]
    [:constraints  constraint-list]
    [:guards   [ logic-expr* ]]
  }

alias-list      ::= [ keyword "atom-string" ]*

atom-collection ::= [ atom-entry* ]

atom-entry ::= atom-spec
             | [ keyword atom-spec ]

atom-spec  ::= "atom-string" | keyword

localized-bond-list ::= [ localized-bond-entry* ]
dative-bond-list   ::= [ dative-bond-entry* ]

atom-ref ::= int | keyword

localized-bond-entry ::= { [:id keyword] :a atom-ref :b atom-ref :type bond-spec }
                      | [ atom-ref atom-ref bond-spec ]
dative-bond-entry   ::= { [:id keyword] :donor atom-ref :acceptor atom-ref
                          [:type "dative-string"] }

bond-spec ::= "bond-string" | bond-keyword

aromatic-list    ::= [ aromatic-entry* ]
multicenter-list ::= [ multicenter-entry* ]
noncovalent-list ::= [ noncovalent-entry* ]

aromatic-entry    ::= { [:id keyword] :atoms [ atom-ref+ ] [:type "aromatic-string"] }
multicenter-entry ::= { [:id keyword] :atoms [ atom-ref+ ] [:type "multicenter-string"] }
noncovalent-entry ::= { [:id keyword] :a atom-ref :b atom-ref :type noncovalent-spec }

noncovalent-spec ::= "noncovalent-string" | noncovalent-keyword
noncovalent-keyword ::= :h-bond | :halogen-bond | :chalcogen-bond | :ionic | :van-der-waals
```

**Dative bond entry.** A dative bond's only identity-bearing content at the molecule-map level is the ordered endpoint pair: **`:donor`** names the atom donating the electron pair, **`:acceptor`** names the atom accepting it. The bond order is fixed at two electrons by definition. **`:donor`** and **`:acceptor`** **MUST** reference distinct atom sites. The optional **`:type`** slot carries a **`dative-string`** payload (**§7.12**) encoding ring-membership constraints (**`#R`**, **`#r`**); the dative-string itself has **no** inherent-field tags and **no** direction token — direction is expressed entirely by the **`:donor`** / **`:acceptor`** assignment. When **`:type`** is absent, the dative bond has no inline constraints.

**Multicenter entry.** The optional **`:type`** slot carries a **`multicenter-string`** payload (**§7.11**) encoding per-system charge, spin, and total electron count. The **`multicenter-string`** subgrammar is independent from **`aromatic-string`** even though they share the same predicate shape.

**Noncovalent kind.** A **`noncovalent-entry`** **MUST** carry **`:type`**. The value is either a **`noncovalent-keyword`** (ground shorthand) or a **`noncovalent-string`** (**§7.13**) carrying the kind as an expression. The five **`noncovalent-keyword`** values expand to the corresponding ground literal in **§7.13** and are accepted wherever a **`noncovalent-spec`** is expected.

**`:id`**. Each structural entry **MAY** include **`:id`** with an EDN **keyword** value. When present, **`:id`** values **MUST** be **pairwise distinct** across **all** entries in the **same** **molecule map** (every list combined).

**Keyword namespace disjointness.** All keyword-shaped identifiers within a single molecule definition — atom ids, atom alias names, structural entry **`:id`** values, and future keyword namespaces (bond alias names) — **MUST** be drawn from **mutually disjoint** namespaces. No two identifier kinds **MAY** share a keyword name within the same molecule map. Alias names **MUST NOT** be valid element symbols (**§7.4**).

**`:atom-aliases`**. The **`alias-list`** defines named atom shorthands scoped to the enclosing molecule map. It is a flat vector of alternating keyword/atom-spec pairs. Each value **MUST** be an **EDN string** carrying an **atom-string** payload. An **`atom-entry`** that is a bare **keyword** (not a string and not in a **`[id entry]`** position) is an alias reference and **MUST** resolve to a key in **`:atom-aliases`**. Aliases are resolved at parse time; the resolved **`atom-string`** is substituted as if written inline. A reference to an undefined alias is an error. Alias definitions **MUST** be bijective: no two alias names **MAY** map to the same atom definition.

**`:constraints`**. Molecule-wide and per-entity constraints, cross-entity relational predicates, sub-pattern anchors, and boolean combinators live here. The canonical grammar appears in **§7.9**. Whole-molecule charge and spin assertions are written as **`{:charge-sum {:sum n}}`** and **`{:spin-sum {:spin spin-literal}}`** entries (omit `:atoms` to range over the whole molecule); a subset is selected by adding `:atoms [...]`. There is no top-level **`:charge`** or **`:spin`** key on the molecule map.

**Inline ids.** An **`atom-entry`** of the form **`[`** *keyword* *atom-spec* **`]`** assigns the keyword as an **id** to the atom at that position. Ids enable symbolic reference from bond endpoints (instead of positional index). Entries with and without ids **MAY** be freely mixed within the same **`:atoms`** vector.

**Endpoints.** Every atom site referenced from a structural relation **MUST** exist under **`:atoms`**, either by positional index (integer) or by id keyword:

- **`localized-bond-entry`** **`:a`** and **`:b`**
- **`dative-bond-entry`** **`:donor`** and **`:acceptor`**
- **`noncovalent-entry`** **`:a`** and **`:b`**
- every reference in an **`aromatic-entry`** or **`multicenter-entry`** **`:atoms`** vector

**Positional index endpoints.** Let **`n`** be the length of **`:atoms`**. An integer endpoint **`i`** with **0 ≤ i < n** denotes the atom at position **`i`**. A keyword endpoint denotes the atom whose inline id is that keyword.

**Empty molecule.** The **`atom-collection`** **MAY** have length **0** (**`[]`**). If **`:atoms`** is empty, **`:bonds`** **MUST** be **`[]`**, and **`:dative`**, **`:aromatic`**, **`:multicenter`**, and **`:noncovalent`** **MUST** be absent or **empty** lists — no bond-like entry **MAY** name a site that is not in **`:atoms`**.

**`aromatic-entry`** and **`multicenter-entry`** **MAY** carry a **`:type`** key whose value is an **aromatic-string** (**§7.10**) encoding per-system charge, spin, and π-electron count.

**`bond-keyword`** (shorthand) is defined in **§7.7**. **`logic-expr`** and **`spin-literal`** are not fully specified in this document (TODO).

### 4.1 Structural validity (within one map)

These rules apply **within** a single **molecule map**. **Constraints across** relation kinds (e.g. the same atom pair in **`:bonds`** and **`:dative`**) are **not** specified here.

**`:bonds` (localized).** The list **MUST NOT** contain two **`localized-bond-entry`** values with the same **unordered** pair of atom sites **{`:a`, `:b`}** (endpoints as a set).

**`:dative`.** The list **MUST NOT** contain two **`dative-bond-entry`** values with the same **unordered** pair **{`:donor`, `:acceptor`}**. A donor→acceptor bond and the reverse acceptor→donor bond between the **same** two atoms violate this rule.

**`:aromatic`.** For every two distinct **`aromatic-entry`** values, the sets of keywords in their **`:atoms`** vectors **MUST** be disjoint. Aromatic systems **MUST NOT** share an atom.

**`:noncovalent`.** For every two distinct **`noncovalent-entry`** values, the sets **{`:a`, `:b`}** **MUST** be disjoint. Noncovalent bonds **MUST NOT** share an atom with another noncovalent bond in the same map.

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
| **`#h`** | yes, when the payload is decimal-only; **`#h=`**, **`#h*`**, etc. are **special** (**§7.3**) |
| **`#n` `#u` `#s` `#v` `#d` `#t` `#r` `#m`** | yes, when the payload is decimal-only |
| **`#a`** | yes when decimal-only; **`#a*`**, **`#a+`**, **`#a!`** are **special** (**§7.3**) |
| **`#i`** | yes, when the payload is decimal-only; bare **`#i`** denotes isotope mass **1** (unusual but permitted) |
| **`#R`** | yes when decimal-only; **`#R*`**, **`#R+`** are **special** (**§7.3**) |

Bond predicates that use **decimal-only** payloads follow the same **`decimal-tail`** rule where applicable (**§7.5**).

In **Query** and **Rule**, any predicate slot that allows a full **`value-expr`** may use **`bool-expr`**, **`*`**, top-level **`nat-set`**, **`decimal-tail`**, etc., as allowed for that tag.

### 5.4 Wildcards, sets, logic, arithmetic

- The **`*`** **wildcard** is allowed in **`value-expr`**, **`element`**, and **`order`**
- **`bool-expr`**: **infix** **`&` `|` `!`**, **relations**, **`::`**, **`+ - * / %`**, unary **`-`**, **`?id`**, **`nat`**, **`(`** **`add-expr`** **`)`**.

**Ground:** no **`bool-expr`** (no **`?`**, **`::`**, relations, logic), no **`element-bind`**, no **`element-ref`**; predicate payloads are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (and tag-specific literals such as **`#h=`**) only where allowed.

**Query:** **`bool-expr`** where allowed; **`decimal-tail`**; **element** / **order** extensions as allowed.

**Rule:** full **`value-expr`**; **element** may use **`element-bind`** / **`element-ref`** (**§6**).

**`<` `>` `<=` `>=` `==`** appear **only** inside **`value-expr`** (predicate payloads). **Dative** donated / accepted pair counts use predicates **`#d`** / **`#t`** (**§7.3**), not bare **`<` `>`** at the top level of the atom-string.

---

## 6. Match semantics and bindings

**Ground molecule, pattern LHS.** The target is **ground** (fully instantiated). The **LHS** of a rule (or query) may still contain **wildcards**, **sets**, **binds**, and **guards**: that is **pattern** data, not an indeterminate molecule.

### 6.1 Inherent fields and derived predicates

**Inherent fields.** Each AST form — atom, localized bond, aromatic system, multicenter bond, dative bond, noncovalent bond — carries a fixed set of **inherent fields**. An inherent field's value **identifies** the entity at that slot. An entity is **ground** iff every inherent field holds a single concrete value (a literal; not a wildcard, set, bind, ref, or unresolved symbolic state). Nothing else affects grounding.

| Form | Inherent fields |
|------|-----------------|
| atom | element, isotope mass (**`#i=`**), charge (**`#c`**), implicit hydrogens (**`#h=`**), lone pairs, spin (unpaired **`#u`**, multiplicity **`#s`**) |
| localized bond | order, charge (**`#c`**), spin (**`#u`**, **`#s`**) |
| aromatic system | charge (**`#c`**), spin (**`#u`**, **`#s`**), π-electron count (**`#e`**) |
| multicenter bond | charge (**`#c`**), spin (**`#u`**, **`#s`**), electron count (**`#e`**) |
| dative bond | ordered endpoint pair — the **`:donor`** / **`:acceptor`** assignment on the map entry (**§4**). The dative-string payload has no inherent-field tags; order is two electrons by definition. |
| noncovalent bond | interaction kind (**`:h-bond`**, **`:halogen-bond`**, **`:chalcogen-bond`**, **`:ionic`**, **`:van-der-waals`**) |

**Derived predicates.** Every predicate admitted in the DSL that is not an inherent field is a **derived predicate** — a topological query evaluated against the target graph once an embedding is proposed. This includes per-atom **`#D`**, **`#X`**, **`#x`**, **`#H`**, **`#R`**, **`#r`** (**§7.3**); the bond-namespace **`#R`**, **`#r`**; per-aromatic, per-multicenter, per-dative ring-membership and ring-size predicates; and the molecule-wide entries of **§7.9**. Derived predicates **filter** matches; they do **not** carry identity and **do not** affect grounding. Adding a derived predicate — even a wildcard-valued one — to a pattern never makes a ground target stop being ground.

### 6.2 Pattern–target match

**Match as solution-set inclusion.** Each attribute slot has a **solution set** — the set of ground values the slot admits. A **literal** (e.g. **`C`**, **`3`**, **`+1`**) admits exactly itself; a **set** (**`{C,N}`**, top-level **`nat-set`**) admits its members; a **wildcard** (**`*`**) admits everything in the slot's value domain; a **`bool-expr`** admits every value for which the expression holds (**§5**); a **special-symbolic** payload (**`#h=`**, **`#i=`**, **`#a*`**, **`#a+`**, **`#a!`**, **`#R*`**, **`#R+`**) admits only its named symbolic state (**§7.3**). For a given slot, the **pattern** matches the **target** iff `solution-set(pattern)` ⊇ `solution-set(target)` — the pattern admits every value the target admits. Match is **not** symmetric.

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

**Molecule-level match.** A molecule-map pattern matches a target iff (a) every **`atom-string`** matches its corresponding target atom-string field-wise — element and every inherent-field predicate payload, per the rules above; (b) every **`bond-string`** matches field-wise; (c) each structural relation (**`:aromatic`**, **`:multicenter`**, **`:dative`**, **`:noncovalent`**) matches per its own inherent fields; (d) every derived predicate holds on the resulting embedding (**§6.1**). Any failure rejects the embedding.

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

**Examples (atom):** **`C`**, **`C#h3`**, **`C#h*`**, **`C#h=`**, **`C#a*`**, **`C#a !`**, **`C#c+`**, **`C#c-`**, **`C#c +`**, **`(?e :: {Cl,Br})#v(?q == 2)`**.

- A **`nat`** and an **`id`** contain **no** internal whitespace.
- A **relational** token is **`<=`**, **`>=`**, **`==`**, or a **single** **`<`** or **`>`** that is **not** part of **`<=` `>=`**. **Multi-character** tokens are one lexical unit.
- An **arithmetic** token is **`+`**, **`-`**, **`*`**, **`/`**, **`%`**. Leading **`+`** / **`-`** on a **`base-expr`** are **`sign`** tokens (**§5**); binary **`+`** / **`-`** appear between **`mult-expr`** operands.
- **Inside** an **element** or **order** **brace set** `{…}`, optional whitespace is allowed **only** immediately before or after a comma separating entries. No whitespace inside an **`element-literal`** or **`order-entry`** (**`nat`**).

**`=` (U+003D).** Plain **`=`** is **not** an operator in **`value-expr`**; equality is **`==`**. **`=`** **MUST NOT** appear in a **Ground** string in a **`value-expr`** position (no **`bool-expr`**). It **MAY** appear inside payloads only as part of **`==`** or other defined tokens.

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
| Ring size | **`#r`** |
| Aromatic π contribution (numeric) | **`#a`** |
| Multicenter valence | **`#m`** |

**Bond-string** **`#u`** / **`#s`** / **`#c`** use the **bond** namespace (**§7.5**); meanings parallel **unpaired electrons**, **multiplicity**, and **bond formal charge**.

**Aromatic-string** **`#c`** / **`#u`** / **`#s`** / **`#e`** use the **aromatic-system** namespace (**§7.10**); **`#e`** denotes the total π-electron count (**`u8`**), other tags parallel the bond namespace.

**`#n`** and **`#s`** together constrain spin-related fields; when both appear, consistency rules are **TBD** (**§3**, open issue).

**Lexical** **`nat`** in the grammar is unbounded; **Ground** validation **MUST** reject values outside the **u8** (or **i8** / **u32** as above) range for the corresponding slot.

**Bond order** (**§7.6**) uses a **discrete** model; **fractional** bond orders **MUST NOT** appear in the **bond-string**. **Aromatic** connectivity **MUST NOT** be encoded as a bond **order**; use the molecule map’s **`:aromatic`** section (**§4**) and ordinary **`:bonds`** entries.

### 7.3 Atom subgrammar

```
atom-string ::= element atom-predicate*

atom-predicate ::= '#' tag payload

tag ::= [A-Za-z_]
```

- **`element`** is first (**§7.4**).
- **Zero or more** **`atom-predicate`** units follow. **Optional** ASCII whitespace **MAY** appear **between** **`element`** and the first **`#`**, and **between** successive predicates.
- **At most one** predicate per **tag letter** per **`atom-string`** (each row of the table below is a **kind**).
- **Canonical order** of predicates after **`element`** (stable serialized form): **`#i`**, **`#c`**, **`#h`**, **`#n`**, **`#u`**, **`#s`**, **`#v`**, **`#d`**, **`#t`**, **`#a`**, **`#m`**, then the uppercase derived predicates **`#D`**, **`#X`**, **`#x`**, **`#H`**, **`#R`**, followed by **`#r`** (ring size). Implementations **MAY** specify further ordering for fields not listed here.

**`payload` parsing.** After trimming leading / trailing whitespace on the **payload** substring, parse as follows:

1. **`#c`**, **Ground** or **Query** / **Rule**: if the trimmed payload is **exactly** **`+`** or **`-`**, the formal charge is **+1** or **−1** (same meaning as **`#c+1`** / **`#c-1`**). Otherwise parse as **`value-expr`** (**§5**) (or the **Ground** subset in **§5.4**).
2. **Any other tag**: parse the payload as **`value-expr`** (or **Ground** subset) unless the payload matches a **special** form below.

**`#i`** follows the usual **§5.3** decimal-only rule: bare **`#i`** denotes mass **1**. This is chemically unusual and implementations **SHOULD** warn, but the form is **valid**.

**Special predicate payloads** (trimmed; **not** parsed as boolean **`!`** — these are **opaque** lexemes for the given tag):

| Form | Tag | Meaning |
|------|-----|---------|
| **`=`** (payload is a single equals sign) | **`#h`** | **Normal / valence-model implicit hydrogen**: implicit H count is whatever the valence model assigns for this **`element`** and the rest of the atom’s fields (**Query** / **Rule** indeterminacy **MAY** apply). |
| **`=`** | **`#i`** | **Natural isotope**: mass number of the **naturally most abundant** isotope of **`element`**. This is the default / expected isotope for each element. |
| **`*`** | **`#h`** | **Wildcard** implicit H count (**Query** / **Rule**). |
| **`*`** | **`#a`** | **No constraint** on aromatic π contribution — equivalent to omitting **`#a`** entirely. |
| **`+`** | **`#a`** | **Sugar** for the constraint **`?a >= 0`**: atom is a member of some aromatic system with an unspecified π contribution (**Query** / **Rule**). |
| **`!`** | **`#a`** | Atom is **not** a member of any aromatic system. Distinct from **`#a0`**: a **`#a0`** atom *is* in an aromatic system and contributes **zero** π electrons (e.g. a carbocation with an empty p orbital participating in a ring current); a **`#a!`** atom has no aromatic membership at all. In **Ground**, a **`#a!`** atom **MUST NOT** appear in any **`:aromatic`** entry. |
| **`*`** | **`#R`** | **No constraint** on ring count — equivalent to omitting **`#R`** entirely. |
| **`+`** | **`#R`** | **Sugar** for the constraint **`?r >= 1`**: atom is in **at least one** ring (**Query** / **Rule**). |
| **`+`** / **`-`** (alone) | **`#c`** | **+1** / **−1** formal charge (**§7.3** above). |

Other **`#h`** / **`#a`** payloads use the usual **`value-expr`** / **`decimal-tail`** rules (**§5**, **§5.3**).

| Tag | Meaning |
|-----|---------|
| **`#i`** | Isotope mass; **special** **`#i=`** (natural isotope, **§7.3**) |
| **`#c`** | Formal charge |
| **`#h`** | Implicit H count; **special** **`#h=`**, **`#h*`** (**§7.3**) |
| **`#n`** | Lone pairs (nonbonding pair count) |
| **`#u`** | Unpaired electron count |
| **`#s`** | Spin multiplicity (2S+1) |
| **`#v`** | **Localized valence**: sum of **bond orders** of **localized** **`:bonds`** edges to **non-hydrogen** neighbors in the **molecular graph** (multiple bonds count by full order). **Excludes** implicit H (**`#h`**), **dative**, **aromatic**-section bonding, **multicenter**, and **noncovalent** contributions — those are separate fields. |
| **`#d`** | Dative **donated** pair count (electrons donated **by** this atom) |
| **`#t`** | Dative **accepted** pair count (“accepted”; electrons accepted **by** this atom) |
| **`#a`** | Aromatic π contribution; **special** **`#a*`**, **`#a+`**, **`#a!`** (**§7.3**) |
| **`#m`** | Multicenter valence |
| **`#D`** | **Degree**: number of neighbors in the molecular graph (SMARTS `D`). Derived predicate evaluated against the target; **not** a ground atom field. |
| **`#X`** | **Connectivity**: degree plus implicit-H count (SMARTS `X`). Derived. |
| **`#x`** | **Ring connectivity**: number of ring bonds at the atom (SMARTS `x`). Derived. |
| **`#H`** | **Total hydrogens**: implicit H count plus explicit H neighbors (SMARTS `H`). Derived. |
| **`#R`** | **Ring count**. Follows the **§5.3** omitted-numeral convention: bare **`#R`** means **1** ring; **`#R<n>`** means exactly **n** rings. **Special** **`#R*`** (no constraint) and **`#R+`** (sugar for **`?r >= 1`**, "in at least one ring"). Derived. |
| **`#r`** | **Ring size**: the atom belongs to a ring of the given size. Derived. |

**Case convention.** **Lowercase** tag letters denote the atom's own state fields (isotope, charge, spin, implicit H, localized valence, π contribution, …) plus the SMARTS-lowercase ring predicates **`#r`** (ring size) and **`#x`** (ring connectivity). **Uppercase** tag letters denote other **derived predicates** (topology queries over the surrounding graph) in the SMARTS-parity set. The two namespaces are **disjoint**: **`#h`** (implicit H slot) and **`#H`** (total H count) coexist, **`#r`** (ring size) and **`#R`** (ring count) coexist, and **`#x`** (ring connectivity) and **`#X`** (connectivity) coexist, without collision.

### 7.4 Element and bond **`order`** (via **`value-expr`**)

The **`element`** nonterminal (**atom-string** prefix) is **literal** | **wildcard `*`** | **brace set** | **`(?` *id* `::` *set* `)`** | **`(?` *id* `)`** (**§7.4** grammar below). The **bond-string** **`order`** prefix (**§7.5**) is a single **`value-expr`** (**§5**), which **subsumes** literal **`nat`**, **`*`**, brace **`nat-set`**, **`(?` *id* `::` *set* `)`**, **`(?` *id* `)`**, and **arithmetic** / logic (e.g. **`1+1`**, **`?o+1`**) where allowed by context.

```
element ::= element-literal | '*' | element-set | element-bind | element-ref
element-set ::= '{' element-literal (',' element-literal)* '}'
element-bind ::= '(' '?' id '::' element-set ')'
element-ref ::= '(' '?' id ')'
element-literal ::= [A-Z][a-z]*
```

- **`element-literal`**: one chemical symbol; **§7.2** (H–Og).
- **`*`**: any element; **invalid** in **Ground** unless narrowed by a containing rule outside this specification.
- **`element-set`**: finite non-empty disjunction of **one or more** **`element-literal`** entries; **§7.2**. **Query** / **Rule** when **Ground** disallows wildcards.
- **`element-bind`**: **Query** / **Rule** only. Introduces a **nominal** variable **`id`** constrained to **membership in** the **`element-set`** (**§6**). **`::`** here means **set membership in a set of element symbols** (**§5**). **Invalid** in **Ground**.
- **`element-ref`**: **Query** / **Rule** only. **Nominal reference**: **`id`** must already be bound as a nominal in rule scope (**§6**). Appears only in the **element** position at the start of the atom-string. No arithmetic on nominal variables.

Optional ASCII whitespace inside **`element-bind`** after **`(`**, before **`)`**, around **`::`**, and around commas in the inner **`element-set`**, per **§7.1**.

### 7.5 Bond subgrammar

**Bond-string** uses a **separate** namespace from **atom-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on bonds (**§7.2**). The **`order`** prefix is a **`value-expr`** (**§5**); see **§7.4** for the parallel with **`element`**.

```
bond-string ::= order bond-predicate*

bond-predicate ::= '#' tag payload

order ::= value-expr
```

**Segmentation.** Let **`order-text`** be the substring of the **`bond-string`** from the first character after any leading whitespace up to (but not including) the first **`#`** that starts a **`bond-predicate`** (**`#`** immediately followed by a **tag** letter, with **no** whitespace between **`#`** and the tag), or the whole **trimmed** string if there is no such **`#`**. **`order-text`** **MUST** be parsed as **`value-expr`** (**§5**).

**Bond predicates.** **Zero or more** **`bond-predicate`** units follow **`order`**. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`a`**, **`R`**, **`r`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#a`**, **`#R`**, **`#r`**.

**Whitespace** between **`#`** and the **tag** letter is **invalid** (**§7.1**).

**`#c` (bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.)

**`#u`** / **`#s`** / **`#r`.** After **`#u`**, **`#s`**, or **`#r`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **No** extra lookahead is required beyond **`value-expr`** termination and the next predicate or end of string.

**`#R` (bond ring count).** Same **special** payloads as atom-level **`#R`** (**§7.3**): bare **`#R`** means **1**; **`#R*`** means no constraint; **`#R+`** is sugar for **`?r >= 1`** ("bond lies in at least one ring").

| Tag | Meaning (bond namespace) |
|-----|---------------------------|
| **`#c`** | Bond formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (bond centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (bond centered); **`u8`** |
| **`#a`** | Bond is a member of some aromatic system declared under **`:aromatic`**. Derived predicate; no payload. Canonical form is **`{:bond [i :aromatic]}`** under **`:constraints`**. |
| **`#R`** | **Ring count**: the bond lies in this many rings. Follows the **§5.3** omitted-numeral convention; **special** **`#R*`**, **`#R+`** (as for atoms, **§7.3**). Derived. |
| **`#r`** | **Ring size**: the bond lies in a ring of this size. Derived. |

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

**Bond entry shorthands.** A **`bond-keyword`** as the **`:bond`** value of a **`localized-bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent bond-string payload. Normative expansion table:

| Keyword | Expands to | Bond order |
|---------|-----------|------------|
| **`:single`** | **`"1"`** | 1 |
| **`:double`** | **`"2"`** | 2 |
| **`:triple`** | **`"3"`** | 3 |
| **`:quadruple`** | **`"4"`** | 4 |

Implementations **MUST** accept these four keywords wherever **`bond-spec`** is expected. No other **`bond-keyword`** values are defined in this revision; unrecognized keywords **MUST** be rejected.

**Atom literals.** Atom literals are **EDN strings** whose contents are **atom-string** payloads (**§7.3** / **§7.4**). Keyword-shaped atom shorthands (via **`:atom-aliases`**) are defined in **§4**.

### 7.8 Future extensions

**Full binary arithmetic.** Generalize **`add-expr`** / **`mult-expr`** to a single **`arith-expr`** with arbitrary **binary** nesting (not only **`nat`** and **`?id`** as **leaves** and one level of **`(` arith `)`**):

```
arith-expr ::= arith-expr mult-op arith-expr
             | arith-expr add-op arith-expr
             | base-arith

base-arith ::= nat | '?' id | '(' arith-expr ')'
```

Precedence of **`mult-op`** over **`add-op`** unchanged. **Membership** **`::`**, **relations** (**`==`**, **`<`**, …), and **logic** (**§5.1**) would still apply **outside** this **`arith-expr`** layer as today. **`sign`*** prefixes (**§5**) would extend to **`arith-expr`** leaves as needed.

**Other** (non-normative placeholders): functions (**`min`**, **`abs`**, …), typed variables, cross-atom **`id`** scope, **chained** relations **`a < b < c`**, etc.

### 7.9 Constraint grammar

Molecule-wide constraints live under the **`:constraints`** key on a **molecule-map** (**§4**). Each entry is a **single-key map** whose key names the constraint kind. Constraints fall into four categories:

- **Entity** — a value-only predicate over one entity. Same payload as the inline string form on that entity; lift/inline (below) move them between scopes.
- **Relational** — a predicate that ties one DAMN entity (dative bond, aromatic system, multicenter bond, noncovalent bond) to atoms, bonds, or atom-predicates. Cannot appear inline.
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
  | { :ring-count          value-expr }
  | { :ring-size           value-expr }

bond-constraint-form ::=
    :aromatic
  | { :ring-count value-expr }
  | { :ring-size value-expr }

dative-bond-constraint-form ::=
    { :ring-count value-expr }
  | { :ring-size value-expr }

aromatic-system-constraint-form  ::= (* uninhabited — no value-only variants yet *)
multicenter-bond-constraint-form ::= (* uninhabited — no value-only variants yet *)
noncovalent-bond-constraint-form ::= (* uninhabited — no value-only variants yet *)

anchor-spec ::=
    { [:atoms             [[atom-ref atom-ref]+]]?
      [:bonds             [[bond-ref bond-ref]+]]?
      [:dative-bonds      [[dative-bond-ref dative-bond-ref]+]]?
      [:aromatic-systems  [[aromatic-system-ref aromatic-system-ref]+]]?
      [:multicenter-bonds [[multicenter-bond-ref multicenter-bond-ref]+]]?
      [:noncovalent-bonds [[noncovalent-bond-ref noncovalent-bond-ref]+]]? }

atom-ref             ::= int | keyword
bond-ref             ::= int | keyword
dative-bond-ref      ::= int | keyword
aromatic-system-ref  ::= int | keyword
multicenter-bond-ref ::= int | keyword
noncovalent-bond-ref ::= int | keyword
```

**Ref resolution.** An integer ref is the **positional** index into the corresponding entity vector on the molecule map: **`atom-ref`** → **`:atoms`**, **`bond-ref`** → **`:bonds`**, **`dative-bond-ref`** → **`:dative`**, **`aromatic-system-ref`** → **`:aromatic`**, **`multicenter-bond-ref`** → **`:multicenter`**, **`noncovalent-bond-ref`** → **`:noncovalent`**. A keyword ref resolves against the **`:id`** declared on the corresponding entry (**§4**). On serialization, implementations **MUST** emit the **`:id`** keyword when one is declared on the referenced entry, falling back to the positional integer otherwise.

**Uninhabited narrow inner forms.** **`:aromatic-system`**, **`:multicenter-bond`**, and **`:noncovalent-bond`** narrow leaves are reserved keys whose inner **`*-constraint-form`** has no inhabited variants in this revision; every cross-entity predicate on those entities is a relational leaf instead. Future revisions **MAY** add value-only variants here without reshaping the surrounding grammar.

**Anchor cardinality.** Each keyed slot in **`anchor-spec`** is optional and may appear at most once; if present, it is a vector of **`(target-side-ref, pattern-side-ref)`** pairs of the same entity kind. An empty **`anchor-spec`** denotes an unanchored sub-pattern (the pattern can embed anywhere). Target-side refs resolve against the outer molecule's metadata; pattern-side refs against the pattern molecule's metadata.

**Molecule-scope subset selectors.** `:charge-sum`, `:spin-sum`, `:bond-order-sum`, and `:connected` accept an **optional** `:atoms` (or `:bonds`) vector. When **omitted**, the predicate ranges over **every** atom (or bond) in the molecule, including atoms added by future structural growth. When **present**, the predicate ranges over the listed entities only. An empty vector `[]` is **distinct** from omission: it selects no entities.

**Sub-pattern materialization.** A **`:sub-pattern`** **`:pattern`** is a full **molecule-map**; its inner constraints are evaluated independently from the outer constraint tree. The pattern carries **no defaults** — values pass through verbatim — so a pattern's atom **`charge: undetermined`** stays **`undetermined`** at match time and behaves as a wildcard (**§5.4**); zero-defaulting that would apply to a ground input does **not** apply inside a pattern.

**Sugar (inline string equivalents).** Narrow leaves whose entity also has a string subgrammar **`packed-string`** form (atom **§7.3**, bond **§7.5**, dative **§7.12**) admit two interchangeable serializations:

- the inline **`#tag`** payload on the entity's **`:type`** string (or, for atoms, on the atom literal directly);
- the **`{:<entity> [ref form]}`** entry in **`:constraints`**.

Parsers **MUST** accept both. Bare per-entity predicates (not nested under **`:and`** / **`:or`** / **`:not`** / **`:sub-pattern`**) **MAY** be emitted in the sugared inline form; nested predicates **MUST** be emitted as **`:constraints`** entries since the inline form has no logical context. Aromatic-system, multicenter, noncovalent, and *any* relational leaf has **no** inline form.

**Lift / inline.** The two storage scopes — inline on the entity (`AtomAst::constraints` etc.) and at molecule scope (`MoleculeAst::constraints` as `{:atom [ref form]}` peers) — are interchangeable for the inline-capable narrow leaves. Implementations **SHOULD** expose:

- **`lift_constraints`** (entity → molecule): drains every inline store into the molecule list as `{:<entity> [ref form]}` peers.
- **`inline_constraints`** (molecule → entity): drains top-level inline-capable narrow leaves from the molecule list into the targeted entity's inline store.

Combinator subtrees, relational leaves, and molecule-scope leaves are never moved by either operation. With multiple top-level entries targeting the same (entity, kind), `inline_constraints` resolves the collision via the entity store's per-kind insert policy (last-wins for unique-kind variants).

**Multiple constraints per entity.** Each per-entity constraint serializes as its **own** entity-constraint entry; implementations **MUST NOT** bundle multiple constraints on the same entity into a single map.

### 7.10 Aromatic system subgrammar

**Aromatic-string** uses a **separate** namespace from **atom-string** and **bond-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on aromatic systems (**§7.2**). It carries **per-aromatic-system** state — overall **charge**, **spin**, and total **π-electron count** — as the **`:type`** value of an **`aromatic-entry`** or **`multicenter-entry`** (**§4**).

```
aromatic-string ::= aromatic-predicate*

aromatic-predicate ::= '#' tag payload
```

**Aromatic predicates.** **Zero or more** **`aromatic-predicate`** units. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (aromatic-system formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.) Same convention as atom (**§7.3**) and bond (**§7.5**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **`#e`** omitted means **1** π-electron; this is chemically unusual but grammatically valid.

| Tag | Meaning (aromatic-system namespace) |
|-----|---------------------------------------|
| **`#c`** | Aromatic-system formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (system-centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (system-centered); **`u8`** |
| **`#e`** | Total π-electron count; **`u8`** |

**No canonical-constraint equivalent.** Aromatic-system charge (**`#c`**), unpaired electrons (**`#u`**), spin multiplicity (**`#s`**), and π-electron count (**`#e`**) live as direct fields on the aromatic-system entity (set by the aromatic-string predicates above) and have **no** canonical **`:constraints`** form in this revision.

### 7.11 Multicenter-bond subgrammar

**Multicenter-string** uses a **separate** namespace from **atom-string**, **bond-string**, and **aromatic-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on multicenter bonds (**§7.2**). It carries **per-multicenter-bond** state — overall **charge**, **spin**, and total bonded **electron count** — as the **`:type`** value of a **`multicenter-entry`** (**§4**).

```
multicenter-string ::= multicenter-predicate*

multicenter-predicate ::= '#' tag payload
```

**Multicenter predicates.** **Zero or more** **`multicenter-predicate`** units. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (multicenter-bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. Same convention as atom (**§7.3**), bond (**§7.5**), and aromatic (**§7.10**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **`#e`** omitted means **1** bonded electron; chemically unusual but grammatically valid.

| Tag | Meaning (multicenter-bond namespace) |
|-----|----------------------------------------|
| **`#c`** | Multicenter-bond formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (bond-centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (bond-centered); **`u8`** |
| **`#e`** | Total bonded electron count; **`u8`** |

**Per-atom participation.** The **`#e`** payload names the total electrons in the multicenter bond; the per-atom share of that count is expressed on each participating atom-string via the **`#m`** predicate (**§7.3**), not inside the multicenter-string. Endpoint references (which atoms the bond spans) live in the **`:atoms`** vector of the **`multicenter-entry`** (**§4**); they **MUST NOT** be encoded inside the multicenter-string.

### 7.12 Dative-bond subgrammar

**Dative-string** carries **ring-membership** constraints on a single **`dative-bond-entry`** (**§4**). It has **no** leading **`order`** token (the order of a dative bond is fixed at two electrons) and **no** inherent-field predicates.

```
dative-string ::= dative-predicate*

dative-predicate ::= '#' tag payload
```

**Dative predicates.** **Zero or more** **`dative-predicate`** units. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`R`**, **`r`**. **Canonical** predicate order (stable serialization): **`#R`**, **`#r`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#R` (dative-bond ring count).** Same **special** payloads as the atom-level and bond-level **`#R`** (**§7.3**, **§7.5**): bare **`#R`** means **1**; **`#R*`** means no constraint; **`#R+`** is sugar for **`?r >= 1`** ("dative bond lies in at least one ring").

**`#r` (dative-bond ring size).** After **`#r`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means ring size **1** (same convention as **§5.3** for decimal-only slots).

| Tag | Meaning (dative-bond namespace) |
|-----|-----------------------------------|
| **`#R`** | **Ring count**: the dative bond lies in this many rings. Derived predicate. |
| **`#r`** | **Ring size**: the dative bond lies in a ring of this size. Derived. |

**Direction.** Dative bonds are intrinsically directional. Direction is carried entirely by the ordered **`:donor`** / **`:acceptor`** assignment on the containing **`dative-bond-entry`** (**§4**); the dative-string itself has **no** direction token. Under pattern matching (**§6**), the embedding MUST map a pattern **`:donor`** to a target **`:donor`** and a pattern **`:acceptor`** to a target **`:acceptor`** — a donor/acceptor swap across the embedding rejects the match.

**Donor / acceptor / cross-bond references.** Donor-side and acceptor-side constraints on the endpoint atoms (equivalent to the **`:donated-pairs`** / **`:accepted-pairs`** atom-constraint forms of **§7.9**, or atom-string **`#d`** / **`#t`** of **§7.3** pinned to one endpoint) attach via the molecule-wide **`:constraints`** section; they **MUST NOT** be encoded inside the dative-string. The same holds for the "parallels another bond" relation and any reference to other molecule-level entities.

### 7.13 Noncovalent-bond subgrammar

**Noncovalent-string** encodes the **interaction kind** of a single **`noncovalent-entry`** (**§4**) as an expression. It has **no** predicates; the kind is the whole payload.

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

**Relation to `noncovalent-keyword` (`§4`).** The five **`noncovalent-keyword`** values (**`:h-bond`** etc.) are a ground-only EDN shorthand on the **`noncovalent-entry`**'s **`:type`**; they expand to the corresponding **`noncovalent-kind-literal`** and are semantically identical. The keyword shorthand **MUST NOT** be used inside a **`noncovalent-kind-set`** or a **`noncovalent-kind-bind`**; those take kind literals only.

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
 :aromatic [{:id :ar1 :atoms [:N :C2 :C3 :C3a :C7a]}
            {:id :ar2 :atoms [:C3a :C4 :C5 :C6 :C7 :C7a]}]}
```

Localized bonds carry the σ-skeleton orders; the aromatic π system is expressed in **`:aromatic`**.

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
