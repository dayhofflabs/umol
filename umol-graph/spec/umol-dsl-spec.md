# umol DSL specification

Normative definition of the molecule **EDN** surface, **contexts**, **molecule map** shape, **value expressions**, **bindings**, and **atom-string** / **bond-string** subgrammars.

The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, **RECOMMENDED**, **MAY**, and **OPTIONAL** in this document are to be interpreted as described in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

**Encoding.** Atom-string and bond-string payloads are Unicode text; **UTF-8** is the usual encoding on the wire and on disk. The subgrammars in **§7** define meaning only for **ASCII**: graphic characters U+0020–U+007E and the ASCII whitespace tokens named in **§7.1**. Any other code point in an atom-string or bond-string is **invalid**.

---

## 1. General

**Homoiconicity.** The same string grammars denote **ground** data and **constraint** patterns. A ground atom-string is a degenerate case of a query atom-string that matches exactly one interpretation.

**Relational molecule.** A molecule map (**§4**) is a set of named relations (atoms, covalent bonds, optional dative / aromatic / multicenter sections, and global fields). Fragment composition is relation-wise merge where defined.

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

An **atom literal** is written **`#atom "`** *atom-string* **`"`** (the **atom-string** grammar is **§7.3** / **§7.4**).

A **bond literal** in **full** form is **`#bond "`** *bond-string* **`"`** (**§7.5**). **§7.7** defines **keyword** shorthands (e.g. **`:single`**) that expand to an equivalent **`#bond`** payload.

**EDN `#atom` / `#bond` vs `#` inside the string.** The **reader dispatch** tokens **`#atom`** and **`#bond`** apply to the whole literal. **Inside** the quoted **atom-string** or **bond-string**, **`#`** starts a **predicate tag** (**§7.3**, **§7.5**). The two uses are **not** the same syntactic class (analogy: section headings **`#`** vs **`##`** in Markdown).

**`:atoms`** is a **vector** of atom entries. Each entry is either:

- an atom literal (string),
- an alias reference (keyword matching a key in **`:atom-aliases`**), or
- an inline-id entry **`[`** *keyword* *atom-entry* **`]`** — a 2-element vector where the first element is a keyword id and the second is an atom literal or alias reference.

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
    :bonds         covalent-bond-list
    [:dative       dative-bond-list]
    [:aromatic     aromatic-list]
    [:multicenter  multicenter-list]
    [:noncovalent  noncovalent-list]
    [:atom-aliases alias-list]
    [:charge       int | nil]
    [:spin         spin-literal | nil]
    [:expect   { :charge int :multiplicity keyword }]
    [:guards   [ logic-expr* ]]
  }

alias-list      ::= [ keyword #atom "atom-string" ]*

atom-collection ::= [ atom-entry* ]

atom-entry ::= #atom "atom-string"
             | keyword
             | [ keyword atom-entry ]

covalent-bond-list ::= [ covalent-bond-entry* ]
dative-bond-list   ::= [ dative-bond-entry* ]

atom-ref ::= int | keyword

covalent-bond-entry ::= { [:id keyword] :a atom-ref :b atom-ref :bond bond-spec }
                      | [ atom-ref atom-ref bond-spec ]
dative-bond-entry   ::= { [:id keyword] :donor atom-ref :acceptor atom-ref :bond bond-spec }

bond-spec ::= #bond "bond-string" | bond-keyword

aromatic-list    ::= [ aromatic-entry* ]
multicenter-list ::= [ multicenter-entry* ]
noncovalent-list ::= [ noncovalent-entry* ]

aromatic-entry    ::= { [:id keyword] :atoms [ keyword+ ] }
multicenter-entry ::= { [:id keyword] :atoms [ keyword+ ] }
noncovalent-entry ::= { [:id keyword] :a keyword :b keyword :bond noncovalent-spec }

noncovalent-spec ::= bond-spec
```

**`:id`**. Each structural entry **MAY** include **`:id`** with an EDN **keyword** value. When present, **`:id`** values **MUST** be **pairwise distinct** across **all** entries in the **same** **molecule map** (every list combined).

**Keyword namespace disjointness.** All keyword-shaped identifiers within a single molecule definition — atom ids, atom alias names, structural entry **`:id`** values, and future keyword namespaces (bond alias names) — **MUST** be drawn from **mutually disjoint** namespaces. No two identifier kinds **MAY** share a keyword name within the same molecule map. Alias names **MUST NOT** be valid element symbols (**§7.4**).

**`:atom-aliases`**. The **`alias-list`** defines named atom shorthands scoped to the enclosing molecule map. It is a flat vector of alternating keyword/atom-spec pairs. Each value **MUST** be a **`#atom`** tagged literal. An **`atom-entry`** that is a bare **keyword** (not a **`#atom`** tagged literal and not in a **`[id entry]`** position) is an alias reference and **MUST** resolve to a key in **`:atom-aliases`**. Aliases are resolved at parse time; the resolved **`atom-string`** is substituted as if written inline. A reference to an undefined alias is an error. Alias definitions **MUST** be bijective: no two alias names **MAY** map to the same atom definition.

**Inline ids.** An **`atom-entry`** of the form **`[`** *keyword* *atom-entry* **`]`** assigns the keyword as an **id** to the atom at that position. Ids enable symbolic reference from bond endpoints (instead of positional index). Entries with and without ids **MAY** be freely mixed within the same **`:atoms`** vector.

**Endpoints.** Every atom site referenced from a structural relation **MUST** exist under **`:atoms`**, either by positional index (integer) or by id keyword:

- **`covalent-bond-entry`** **`:a`** and **`:b`**
- **`dative-bond-entry`** **`:donor`** and **`:acceptor`**
- **`noncovalent-entry`** **`:a`** and **`:b`**
- every reference in an **`aromatic-entry`** or **`multicenter-entry`** **`:atoms`** vector

**Positional index endpoints.** Let **`n`** be the length of **`:atoms`**. An integer endpoint **`i`** with **0 ≤ i < n** denotes the atom at position **`i`**. A keyword endpoint denotes the atom whose inline id is that keyword.

**Empty molecule.** The **`atom-collection`** **MAY** have length **0** (**`[]`**). If **`:atoms`** is empty, **`:bonds`** **MUST** be **`[]`**, and **`:dative`**, **`:aromatic`**, **`:multicenter`**, and **`:noncovalent`** **MUST** be absent or **empty** lists — no bond-like entry **MAY** name a site that is not in **`:atoms`**.

**`aromatic-entry`** and **`multicenter-entry`** **MAY** include additional keys (**`:electrons`**, **`:charge`**, **`:spin`**, …) when specified here or by an implementation; those fields are **not** fully normative in this revision beyond **`:id`** and **`:atoms`**.

**`bond-keyword`** (shorthand) is defined in **§7.7**. **`logic-expr`** and **`spin-literal`** are not fully specified in this document (TODO).

### 4.1 Structural validity (within one map)

These rules apply **within** a single **molecule map**. **Constraints across** relation kinds (e.g. the same atom pair in **`:bonds`** and **`:dative`**) are **not** specified here.

**`:bonds` (covalent).** The list **MUST NOT** contain two **`covalent-bond-entry`** values with the same **unordered** pair of atom sites **{`:a`, `:b`}** (endpoints as a set).

**`:dative`.** The list **MUST NOT** contain two **`dative-bond-entry`** values with the same **unordered** pair **{`:donor`, `:acceptor`}**. A donor→acceptor bond and the reverse acceptor→donor bond between the **same** two atoms violate this rule.

**`:aromatic`.** For every two distinct **`aromatic-entry`** values, the sets of keywords in their **`:atoms`** vectors **MUST** be disjoint. Aromatic systems **MUST NOT** share an atom.

**`:noncovalent`.** For every two distinct **`noncovalent-entry`** values, the sets **{`:a`, `:b`}** **MUST** be disjoint. Noncovalent bonds **MUST NOT** share an atom with another noncovalent bond in the same map.

**`:expect`** is a **checksum** for tests: implementations **MAY** validate it at parse time and **MUST NOT** treat it as part of the persistent chemical model.

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
| **`#n` `#u` `#s` `#v` `#d` `#r` `#m`** | yes, when the payload is decimal-only |
| **`#a`** | yes when decimal-only; **`#a*`**, **`#a!`** are **special** (**§7.3**) |
| **`#i`** | **no** — isotope mass requires a non-empty **`nat`** (or non-**Ground** **`value-expr`** where allowed) |

Bond predicates that use **decimal-only** payloads follow the same **`decimal-tail`** rule where applicable (**§7.5**).

In **Query** and **Rule**, any predicate slot that allows a full **`value-expr`** may use **`bool-expr`**, **`*`**, top-level **`nat-set`**, **`decimal-tail`**, etc., as allowed for that tag.

### 5.4 Wildcards, sets, logic, arithmetic

- The **`*`** **wildcard** is allowed in **`value-expr`**, **`element`**, and **`order`**
- **`bool-expr`**: **infix** **`&` `|` `!`**, **relations**, **`::`**, **`+ - * / %`**, unary **`-`**, **`?id`**, **`nat`**, **`(`** **`add-expr`** **`)`**.

**Ground:** no **`bool-expr`** (no **`?`**, **`::`**, relations, logic), no **`element-bind`**, no **`element-ref`**; predicate payloads are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (and tag-specific literals such as **`#h=`**) only where allowed.

**Query:** **`bool-expr`** where allowed; **`decimal-tail`**; **element** / **order** extensions as allowed.

**Rule:** full **`value-expr`**; **element** may use **`element-bind`** / **`element-ref`** (**§6**).

**`<` `>` `<=` `>=` `==`** appear **only** inside **`value-expr`** (predicate payloads). **Dative** donated / accepted pair counts use predicates **`#d`** / **`#r`** (**§7.3**), not bare **`<` `>`** at the top level of the atom-string.

---

## 6. Bindings and rule application

**Ground molecule, pattern LHS.** The target is **ground** (fully instantiated). The **LHS** of a rule (or query) may still contain **wildcards**, **sets**, **binds**, and **guards**: that is **pattern** data, not an indeterminate molecule.

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
- **Outside** the quoted string, **EDN** **`#atom`** / **`#bond`** are unrelated (**§2**).

**Payload extraction.** A **predicate** is **`#`**, one **tag** character **`[A-Za-z_]`**, and a **payload** consisting of all following characters up to (but not including) the **next** **`#`** or **end of string**, after **whitespace normalization** for the purpose of **tokenizing** the payload as **`value-expr`**: the payload text **MAY** contain ignored whitespace between **`value-expr`** tokens as in **§5**. The **payload** **MUST NOT** contain **`#`**.

**Examples (atom):** **`C`**, **`C#h3`**, **`C#h*`**, **`C#h=`**, **`C#a*`**, **`C#a !`**, **`C#c+`**, **`C#c-`**, **`C#c +`**, **`(?e :: {Cl,Br})#v(?q == 2)`**.

- A **`nat`** and an **`id`** contain **no** internal whitespace.
- A **relational** token is **`<=`**, **`>=`**, **`==`**, or a **single** **`<`** or **`>`** that is **not** part of **`<=` `>=`**. **Multi-character** tokens are one lexical unit.
- An **arithmetic** token is **`+`**, **`-`**, **`*`**, **`/`**, **`%`**. Leading **`+`** / **`-`** on a **`base-expr`** are **`sign`** tokens (**§5**); binary **`+`** / **`-`** appear between **`mult-expr`** operands.
- **Inside** an **element** or **order** **brace set** `{…}`, optional whitespace is allowed **only** immediately before or after a comma separating entries. No whitespace inside an **`element-literal`** or **`order-entry`** (**`nat`**).

**`=` (U+003D).** Plain **`=`** is **not** an operator in **`value-expr`**; equality is **`==`**. **`=`** **MUST NOT** appear in a **Ground** string in a **`value-expr`** position (no **`bool-expr`**). It **MAY** appear inside payloads only as part of **`==`** or other defined tokens.

### 7.2 Numerical limits

**Chemical elements.** Any **`element-literal`**, any entry in an **`element-set`**, and any **nominal** binding or reference (**`element-bind`**, **`element-ref`**) **MUST** refer only to elements from **hydrogen** (**H**) through **oganesson** (**Og**). Implementations **MUST** reject symbols outside that range in **Ground**; **Query** / **Rule** **SHOULD** use the same restriction unless explicitly documented otherwise.

**Charges.** **Formal charge** on atoms (**`#c`**), **formal bond charge** (**`#c`** on **bond-string**), and molecule-map **`:charge`** where integral **MUST** fit a **signed 8-bit** integer (**−128…127**). The **`#c`** payload is a **`value-expr`** (or **Ground** subset) that evaluates to the signed charge, including the **special** forms **`+`** / **`-`** for **±1** (**§7.3**), e.g. **`#c2`**, **`#c-2`**, **`#c+`**, **`#c-`**.

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
| Dative accepted pairs | **`#r`** (“received”) |
| Aromatic π contribution (numeric) | **`#a`** |
| Multicenter valence | **`#m`** |

**Bond-string** **`#u`** / **`#s`** / **`#c`** use the **bond** namespace (**§7.5**); meanings parallel **unpaired electrons**, **multiplicity**, and **bond formal charge**.

**`#n`** and **`#s`** together constrain spin-related fields; when both appear, consistency rules are **TBD** (**§3**, open issue).

**Lexical** **`nat`** in the grammar is unbounded; **Ground** validation **MUST** reject values outside the **u8** (or **i8** / **u32** as above) range for the corresponding slot.

**Bond order** (**§7.6**) uses a **discrete** model; **fractional** bond orders **MUST NOT** appear in the **bond-string**. **Aromatic** connectivity **MUST NOT** be encoded as a bond **order** in **`#bond`**; use the molecule map’s **`:aromatic`** section (**§4**) and ordinary **`:bonds`** entries.

### 7.3 Atom subgrammar

```
atom-string ::= element atom-predicate*

atom-predicate ::= '#' tag payload

tag ::= [A-Za-z_]
```

- **`element`** is first (**§7.4**).
- **Zero or more** **`atom-predicate`** units follow. **Optional** ASCII whitespace **MAY** appear **between** **`element`** and the first **`#`**, and **between** successive predicates.
- **At most one** predicate per **tag letter** per **`atom-string`** (each row of the table below is a **kind**).
- **Canonical order** of predicates after **`element`** (stable serialized form): **`#i`**, **`#c`**, **`#h`**, **`#n`**, **`#u`**, **`#s`**, **`#v`**, **`#d`**, **`#r`**, **`#a`**, **`#m`**. Implementations **MAY** specify further ordering for fields not listed here.

**`payload` parsing.** After trimming leading / trailing whitespace on the **payload** substring, parse as follows:

1. **`#c`**, **Ground** or **Query** / **Rule**: if the trimmed payload is **exactly** **`+`** or **`-`**, the formal charge is **+1** or **−1** (same meaning as **`#c+1`** / **`#c-1`**). Otherwise parse as **`value-expr`** (**§5**) (or the **Ground** subset in **§5.4**).
2. **Any other tag**: parse the payload as **`value-expr`** (or **Ground** subset) unless the payload matches a **special** form below.

**`#i`** **MUST** have a non-empty payload in **Ground** (at least one **`digit`**, or full **`value-expr`** where allowed).

**Special predicate payloads** (trimmed; **not** parsed as boolean **`!`** — these are **opaque** lexemes for the given tag):

| Form | Tag | Meaning |
|------|-----|---------|
| **`=`** (payload is a single equals sign) | **`#h`** | **Normal / valence-model implicit hydrogen**: implicit H count is whatever the valence model assigns for this **`element`** and the rest of the atom’s fields (**Query** / **Rule** indeterminacy **MAY** apply). |
| **`*`** | **`#h`** | **Wildcard** implicit H count (**Query** / **Rule**). |
| **`*`** | **`#a`** | **Wildcard** aromatic π contribution (**Query** / **Rule**). |
| **`!`** | **`#a`** | Atom is **not** a member of any aromatic system. Distinct from **`#a0`**: a **`#a0`** atom *is* in an aromatic system and contributes **zero** π electrons (e.g. a carbocation with an empty p orbital participating in a ring current); a **`#a!`** atom has no aromatic membership at all. In **Ground**, a **`#a!`** atom **MUST NOT** appear in any **`:aromatic`** entry. |
| **`+`** / **`-`** (alone) | **`#c`** | **+1** / **−1** formal charge (**§7.3** above). |

Other **`#h`** / **`#a`** payloads use the usual **`value-expr`** / **`decimal-tail`** rules (**§5**, **§5.3**).

| Tag | Meaning |
|-----|---------|
| **`#i`** | Isotope mass |
| **`#c`** | Formal charge |
| **`#h`** | Implicit H count; **special** **`#h=`**, **`#h*`** (**§7.3**) |
| **`#n`** | Lone pairs (nonbonding pair count) |
| **`#u`** | Unpaired electron count |
| **`#s`** | Spin multiplicity (2S+1) |
| **`#v`** | **Localized valence**: sum of **bond orders** of **covalent** **`:bonds`** edges to **non-hydrogen** neighbors in the **molecular graph** (multiple bonds count by full order). **Excludes** implicit H (**`#h`**), **dative**, **aromatic**-section bonding, **multicenter**, and **noncovalent** contributions — those are separate fields. |
| **`#d`** | Dative **donated** pair count (electrons donated **by** this atom) |
| **`#r`** | Dative **accepted** pair count (“received”; electrons accepted **by** this atom) |
| **`#a`** | Aromatic π contribution; **special** **`#a*`**, **`#a!`** (**§7.3**) |
| **`#m`** | Multicenter valence |

**Convention:** use **lowercase** ASCII for **`tag`** letters in authoring.

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

**Bond predicates.** **Zero or more** **`bond-predicate`** units follow **`order`**. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**. No other **`#tag`** letters are defined in **bond-string** in this revision.

**Whitespace** between **`#`** and the **tag** letter is **invalid** (**§7.1**).

**`#c` (bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.)

**`#u`** / **`#s`.** After **`#u`** or **`#s`**, parse a **`value-expr`** (**§5**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.3** for decimal-only slots). **No** extra lookahead is required beyond **`value-expr`** termination and the next predicate or end of string.

| Tag | Meaning (bond namespace) |
|-----|---------------------------|
| **`#c`** | Bond formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (bond centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (bond centered); **`u8`** |

**Bond order values** in **`#bond`** **MUST NOT** be **fractional** after evaluation (**§7.6**). **Aromatic** bond **order** as a distinct category **MUST NOT** be used in **`order`**; use **§4** instead.

### 7.6 Bond order

**Semantic model** for **covalent** bond order in **`#bond`** (the **`order`** nonterminal is **`value-expr`**, **§7.5**):

- **Discrete** orders **1**, **2**, **3**, and **4** after any **arithmetic** and binding.
- **Any** order: **`*`** as **`value-expr`** (**Query** / **Rule**).
- **Finite set**: top-level **`nat-set`** in **`value-expr`** (e.g. **`{1,2,3}`** or **`{2}`**).
- **Arithmetic and constraints**: full **`value-expr`** on **`order`**, including **`add-expr`**, **`::`** **`nat-set`**, **`bool-expr`**, and **`?id`** binds, subject to **Ground** restrictions below.

In **Ground**, **`order-text`** **MUST** denote a single definite order in **{1,2,3,4}**: **`*`**, a top-level **`nat-set`** whose entries are **only** **1**–**4**, **`(?` *id* `::` *set* `)`** / **`(?` *id* `)`** only where the implementation resolves them to one value, or **`value-expr`** that is **only** **`sign`*** **`nat`** with value **1**–**4** (no **`?`**, **`::`**, relations, logic, or **`(`** … **`)`**). **Query** / **Rule** **MAY** use the full **`value-expr`** grammar on **`order`**.

This section does **not** define **`bond-keyword`** shorthands; see **§7.7**.

### 7.7 Bond and atom literals

**Bond entry shorthands.** A **`bond-keyword`** as the **`:bond`** value of a **`covalent-bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent **`#bond`** payload. Normative expansion table:

| Keyword | Expands to | Bond order |
|---------|-----------|------------|
| **`:single`** | **`#bond "1"`** | 1 |
| **`:double`** | **`#bond "2"`** | 2 |
| **`:triple`** | **`#bond "3"`** | 3 |
| **`:quadruple`** | **`#bond "4"`** | 4 |

Implementations **MUST** accept these four keywords wherever **`bond-spec`** is expected. No other **`bond-keyword`** values are defined in this revision; unrecognized keywords **MUST** be rejected.

**Atom literals.** The **EDN** **`#atom`** tag and **atom-string** payload are fully defined by **§7.3** / **§7.4**. **Additional** tagged atom forms or keywords **MAY** be listed here when introduced.

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

---

## 8. Molecule map examples (non-normative)

Examples use the vector **`:atoms`** form with inline ids. Bond entries show **`:id`** where useful.

### 8.1 Methanol (CH₃OH) — Ground, L1

```clojure
{:atoms [[:C  #atom "C#h3"]
         [:O  #atom "O#h1"]
         [:H  #atom "H"]]
 :bonds [[:C :O :single]
         [:O :H :single]]}
```

The **`H`** atom here represents an **explicit** hydrogen (e.g. a hydroxyl H one wishes to name). Implicit H counts on **`C`** (**`#h3`**) and **`O`** (**`#h1`**) already account for the remaining hydrogens.

### 8.2 Indole — Ground, L1, aromatic ring

```clojure
{:atoms [[:N   #atom "N#h1"]
         [:C2  #atom "C"]
         [:C3  #atom "C"]
         [:C3a #atom "C"]
         [:C4  #atom "C#h1"]
         [:C5  #atom "C#h1"]
         [:C6  #atom "C#h1"]
         [:C7  #atom "C#h1"]
         [:C7a #atom "C"]]
 :bonds [[:N :C2 :single] [:C2 :C3 :double] [:C3 :C3a :single]
         [:C3a :C7a :single] [:C7a :N :single] [:C3a :C4 :single]
         [:C4 :C5 :double] [:C5 :C6 :single] [:C6 :C7 :double]
         [:C7 :C7a :single]]
 :aromatic [{:id :ar1 :atoms [:N :C2 :C3 :C3a :C7a]}
            {:id :ar2 :atoms [:C3a :C4 :C5 :C6 :C7 :C7a]}]}
```

Covalent bonds carry the σ-skeleton orders; the aromatic π system is expressed in **`:aromatic`**.

### 8.3 Substructure query — L2

Match any carbon with at least two implicit hydrogens that is directly bonded to a nitrogen:

```clojure
{:atoms [[:C #atom "C#h(?h >= 2)"]
         [:N #atom "N"]]
 :bonds [[:C :N :single]]}
```

**`(?h >= 2)`** is a **`bool-expr`** payload on **`#h`**; **`?h`** is bound to the matched atom's implicit H count.

### 8.4 Transformation rule — L3

Replace a primary amine carbon (C with three H and bonded to NH₂) with a quaternary carbon (no H, same bond to nitrogen):

```clojure
{:lhs {:atoms [[:C #atom "C#h3"]
               [:N #atom "N#h2"]]
       :bonds [[:C :N :single]]}
 :rhs {:atoms [[:C #atom "C#h0"]
               [:N #atom "N#h2"]]
       :bonds [[:C :N :single]]}}
```

(The **`:lhs`** / **`:rhs`** wrapping is a rule-level convention, not a molecule map key — not normative here.)
