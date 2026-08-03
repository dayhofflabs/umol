# umol DSL specification

Normative definition of the molecule **EDN** surface, **contexts**, **molecule map** shape, **value expressions**, **bindings**, and **atom-string** / **bond-string** subgrammars.

The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, **RECOMMENDED**, **MAY**, and **OPTIONAL** in this document are to be interpreted as described in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

**Encoding.** Atom-string and bond-string payloads are Unicode text; **UTF-8** is the usual encoding on the wire and on disk. The subgrammars in **§7** define meaning only for **ASCII**: graphic characters U+0020–U+007E and the ASCII whitespace tokens named in **§7.1**. Any other code point in an atom-string or bond-string is **invalid**.

---

## 1. General

**Homoiconicity.** The same string grammars denote **ground** data and **constraint** patterns. A ground atom-string is a degenerate case of a query atom-string that matches exactly one interpretation.

**Relational molecule.** A molecule map (**§4**) is a set of named relations (atoms, localized bonds, optional dative / aromatic / multicenter sections, and global fields). Fragment composition is relation-wise merge where defined.

**EDN and rules.** EDN carries the relational structure. **Rule** evaluation (pattern **LHS** → product **RHS**, **§6**) is a separate computation layer: it **MAY** consume and produce molecule maps that use the same surface notation. The **reaction map** (**§8**) is the operational encoding of such a rule — the **LHS** plus an ordered **`:deltas`** edit list whose application yields the **RHS**.

**Term algebra levels (non-normative sketch).** Four levels of expressiveness exist in this algebra, each a strict superset of the previous:

- **Ground**: fully instantiated molecules; no wildcards, binds, or logic. Every numeric slot is a concrete integer; element is a single symbol.
- **Query / Pattern**: adds wildcards (**`*`**), element/numeric sets, **`bool-expr`**, and **`?id`** references. Sufficient for substructure queries.
- **Rule**: adds **element variables** and cross-atom **`id`** scope (**§6**). Sufficient for transformation rules.
- **Compound**: rules whose RHS produces molecule maps that are themselves other terms; higher-order composition.

The subgrammars in **§5.1**, **§7** define forms that are syntactically valid across all levels; which forms are *semantically* allowed depends on which level is in force (**§3**).

**Case sensitivity.** **`atom-string`**, **`bond-string`**, and **`value-expr`** lexing (**§5.1**, **§7**) is **case-sensitive** throughout: e.g. **`#a`** and **`#A`** are distinct predicate tags; **`?x`** and **`?X`** are distinct **`id`**s; **`element-literal`** (**§5.2**) follows **IUPAC** element casing (**`Cl`**, **`Br`**, not arbitrary case folding). Implementations **MUST NOT** treat these fragments as case-insensitive (unlike **Fortran**-style languages).

---

## 2. EDN representation

Molecule data **SHOULD** use **EDN** maps and **Clojure-style** **keywords** for atom site labels.

An **atom literal** is an **EDN string** whose contents are an *atom-string* (**§7.3** / **§5.2**).

A **bond literal** in **full** form is an **EDN string** whose contents are a *bond-string* (**§7.4**). **§7.6** defines **keyword** shorthands (e.g. **`:single`**) that expand to an equivalent bond-string payload.

**`#` only appears inside the string.** Within an **atom-string** or **bond-string**, **`#`** starts a **predicate tag** (**§7.3**, **§7.4**); no reader-dispatch tag is used at the EDN layer.

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
| **Ground** | Wildcards, element **`*`**, element **brace-set**, **`?`**, **`bool-expr`**, and **element variables** are **invalid** unless this specification explicitly allows them for that slot. |
| **Query**  | Wildcards and constraints per **§5.1** and **§7**. |
| **Rule**   | Binds, sets, boolean expressions, and arithmetic per **§5.1**, **§6**, and **§7**. |

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
  }

atom-alias-list      ::= [ keyword "atom-string" ]*

atom-list ::= [ atom-entry* ]

atom-entry ::= atom-spec
             | [ keyword atom-spec ]

atom-spec  ::= "atom-string" | keyword

bond-list ::= [ bond-entry* ]
dative-bond-list   ::= [ dative-bond-entry* ]

atom-ref ::= int | keyword

bond-entry ::= { [:id keyword] :atoms [ atom-ref atom-ref ] :type bond-spec }
             | [ atom-ref atom-ref bond-spec ]
dative-bond-entry   ::= { [:id keyword] :donors [ atom-ref+ ]
                          :acceptor atom-ref :type dative-bond-spec }

bond-spec ::= "bond-string" | bond-keyword
dative-bond-spec ::= "dative-string" | dative-keyword

aromatic-system-list    ::= [ aromatic-system-entry* ]
multicenter-bond-list ::= [ multicenter-bond-entry* ]
noncovalent-bond-list ::= [ noncovalent-bond-entry* ]

aromatic-system-entry    ::= { [:id keyword] :atoms [ atom-ref+ ] :type "aromatic-string" }
multicenter-bond-entry ::= { [:id keyword] :atoms [ atom-ref+ ] :type "multicenter-string" }
noncovalent-bond-entry ::= { [:id keyword] :atoms [ atom-ref atom-ref ] :type noncovalent-bond-spec }

noncovalent-bond-spec ::= "noncovalent-string"

stereo-atom-list ::= [ stereo-atom-entry* ]
stereo-bond-list ::= [ stereo-bond-entry* ]

stereo-atom-entry ::= { [:id keyword] :site atom-ref :ligands [ ligand-ref+ ] :type stereo-spec }
stereo-bond-entry ::= { [:id keyword] :site bond-ref :ligands [ ligand-ref+ ] :type stereo-spec }

ligand-ref  ::= atom-ref | [ :h atom-ref ] | [ :lp atom-ref ]
stereo-spec ::= "stereo-string" | stereo-keyword
stereo-keyword ::= :ccw | :cw | :z | :e
```

**`:type` is mandatory.** Every structural entry that has a DSL surface — **`bond-entry`**, **`dative-bond-entry`**, **`aromatic-system-entry`**, **`multicenter-bond-entry`**, **`noncovalent-bond-entry`**, **`stereo-atom-entry`**, **`stereo-bond-entry`** — **MUST** carry a **`:type`** key. The payload is a subgrammar string (or its EDN-keyword shorthand where defined); an entry without **`:type`** is a parse error.

An empty string **`:type ""`** is a **parse error** for **every** structural subgrammar — each **MUST** begin with a leading inherent-field token:

- **`aromatic-string`** (**§7.8**) and **`multicenter-string`** (**§7.9**) lead with an **electron counts** specification — **`*`** (undetermined) or a **`[n,n,…]`** vector (the counterpart of the bond-string's order or the atom-string's element). Empty is invalid; use **`"*"`** for undetermined counts.
- **`bond-string`** (**§7.4**), **`dative-string`** (**§7.7**), and **`noncovalent-string`** (**§7.10**) **MUST** begin with a leading inherent-field token (bond order, dative order, noncovalent kind). Use the appropriate keyword shorthand (e.g. **`:single`**) or the literal token (e.g. **`"1"`**, **`"Hbd"`**).
- **`stereo-string`** (**§7.11**) **MUST** begin with a leading **`class`** token (**`Th`** / **`Ct`** / **`Ax`** / **`Sp`** / **`Tb`** / **`Oh`**) followed by a **`coset`**. Use a literal token (e.g. **`"Th1"`**) or a **`stereo-keyword`** (**`:ccw`** / **`:cw`** / **`:z`** / **`:e`**).

**Dative bond entry.** A dative bond entry binds a **single acceptor** to **one or more donors** (a coordination center): **`:acceptor`** names the atom accepting the electron pair(s); **`:donors`** is a **vector** of one or more donating atoms. The leading **`order`** token of the **`dative-string`** payload (**§7.7**) records the number of donated pairs; one shared **`:type`** covers every donor→acceptor edge of the entry. The **`:acceptor`** and every donor **MUST** reference distinct atom sites. The mandatory **`:type`** slot carries a **`dative-string`** (**§7.7**) — order plus optional aromatic constraint (**`#a`**) and the ring-membership predicate (**`#R`**); its leading order parallels the bond-string's (**§7.4** / **§7.5**). The dative-string has **no** direction token — direction is expressed entirely by the **`:donors`** / **`:acceptor`** assignment.

**Multicenter entry.** The mandatory **`:type`** slot carries a **`multicenter-string`** payload (**§7.9**) — a leading **electron counts** specification then per-system charge, unpaired-electron count and multiplicity, and the optional asserted total electron count (**`#e<n>`**). The **`multicenter-string`** subgrammar is independent from **`aromatic-string`** even though they share the same predicate shape.

**Per-atom electron counts (aromatic and multicenter entries).** The **per-atom** electron contributions are the **mandatory leading** specification of the **`aromatic-string`** / **`multicenter-string`** (**§7.8** / **§7.9**), not a map key: a leading **`*`** (the whole vector undetermined) **or** a **`[n,n,…]`** vector of concrete integers, one per member atom. A concrete vector **MUST** have the same length as the entry's **`:atoms`** vector — position **`i`** is the contribution of the atom at position **`i`** of **`:atoms`**. The electron counts are independent of the optional **`#e`** total — when both are present, downstream validation **MAY** require their **sum** to equal **`#e`** on ground inputs.

**Noncovalent kind.** A **`noncovalent-bond-entry`** **MUST** carry **`:type`**. The value is a **`noncovalent-string`** (**§7.10**) carrying the interaction kind (e.g. **`"Hbd"`**).

**Stereo atom / stereo bond entry.** A **`stereo-atom-entry`** overlays a coordination-stereo configuration (tetrahedral and the higher geometries) on an atom site; a **`stereo-bond-entry`** overlays a cis/trans configuration on a bond site.

- **`:site`** names the bearing entity — an **`atom-ref`** for a **`stereo-atom-entry`**, a **`bond-ref`** for a **`stereo-bond-entry`**. A site **MUST** carry at most one stereo element (**§4.1**).
- **`:ligands`** is the **ordered** local reference frame against which the **`:type`** **`coset`** index is numbered; **order is significant**. Each **`ligand-ref`** is either a plain **`atom-ref`** (a neighbor atom) or a **virtual ligand** — **`[:h atom-ref]`** for an implicit hydrogen borne by the named atom, or **`[:lp atom-ref]`** for a lone pair on the named atom. For a **`stereo-atom-entry`** the bearing atom of every virtual ligand is the site atom; for a **`stereo-bond-entry`** it is the relevant double-bond terminus.
- **`:type`** carries a **`stereo-string`** (**§7.11**) — the **`class`** (**`Th`** / **`Ct`** / **`Sp`** / **`Tb`** / **`Oh`**) plus the **`coset`** index — or a **`stereo-keyword`** shorthand (**§7.11**). The coset index is a dense per-class arrangement number relative to the **`:ligands`** order, not a permutation rank.

**`:id`**. Each structural entry **MAY** include **`:id`** with an EDN **keyword** value. When present, **`:id`** values **MUST** be **pairwise distinct** across **all** entries in the **same** **molecule map** (every list combined).

**Keyword namespace disjointness.** All keyword-shaped identifiers within a single molecule definition — atom ids, atom alias names, structural entry **`:id`** values, and future keyword namespaces (bond alias names) — **MUST** be drawn from **mutually disjoint** namespaces. No two identifier kinds **MAY** share a keyword name within the same molecule map. Alias names **MUST NOT** be valid element symbols (**§5.2**).

**Resolved representation and metadata.** The literal **`:id`** field and inline atom-id form assign DSL **keywords**; they do not replace the numerical identifiers used by the resolved AST. A metadata-preserving parser **MUST** retain a bidirectional association between each keyword and the resolved entity, including its entity kind, together with the atom-alias definitions. A metadata-free parse **MAY** discard this surface information after resolving all references. Metadata-preserving rendering **MUST** reject metadata that names an entity absent from the paired molecule or reaction span, or a reaction entity not introduced in the appropriate lhs or delta scope.

**`:atom-aliases`**. The **`atom-alias-list`** defines named atom shorthands scoped to the enclosing molecule map. It is a flat vector of alternating keyword/atom-spec pairs. Each value **MUST** be an **EDN string** carrying an **atom-string** payload. An **`atom-entry`** that is a bare **keyword** (not a string and not in a **`[id entry]`** position) is an alias reference and **MUST** resolve to a key in **`:atom-aliases`**. Aliases are resolved at parse time; the resolved **`atom-string`** is substituted as if written inline. A reference to an undefined alias is an error. Alias definitions **MUST** be bijective: no two alias names **MAY** map to the same atom definition.

**`:constraints`**. Molecule-wide and per-entity constraints, cross-entity relational predicates, sub-pattern anchors, and boolean combinators live here. The canonical grammar appears in **§7.12**. Whole-molecule charge and unpaired-electron coupling assertions are written as **`{:charge-sum {:sum n}}`** and **`{:unpaired-electron-coupling {:unpaired-electrons {:count n :multiplicity m}}}`** entries (omit `:atoms` to range over the whole molecule); a subset is selected by adding `:atoms [...]`. There is no top-level **`:charge`** or **`:spin`** key on the molecule map.

**Inline ids.** An **`atom-entry`** of the form **`[`** *keyword* *atom-spec* **`]`** assigns the keyword as an **id** to the atom at that position. Ids enable symbolic reference from bond endpoints (instead of positional index). Entries with and without ids **MAY** be freely mixed within the same **`:atoms`** vector.

**Endpoints.** Every atom site referenced from a structural relation **MUST** exist under **`:atoms`**, either by positional index (integer) or by id keyword:

- **`bond-entry`** **`:atoms`** (exactly two)
- **`dative-bond-entry`** **`:donors`** and **`:acceptor`**
- **`noncovalent-bond-entry`** **`:atoms`** (exactly two)
- every reference in an **`aromatic-system-entry`** or **`multicenter-bond-entry`** **`:atoms`** vector

**Positional index endpoints.** Let **`n`** be the length of **`:atoms`**. An integer endpoint **`i`** with **0 ≤ i < n** denotes the atom at position **`i`**. A keyword endpoint denotes the atom whose inline id is that keyword.

**Empty molecule.** The **`atom-list`** **MAY** have length **0** (**`[]`**). If **`:atoms`** is empty, **`:bonds`** **MUST** be **`[]`**, and **`:dative-bonds`**, **`:aromatic-systems`**, **`:multicenter-bonds`**, **`:noncovalent-bonds`**, **`:stereo-atoms`**, and **`:stereo-bonds`** **MUST** be absent or **empty** lists — no bond-like or stereo entry **MAY** name a site that is not in **`:atoms`**.

**`bond-keyword`** and **`dative-keyword`** (shorthands) are defined in **§7.6**.

### 4.1 Structural validity (within one map)

These rules apply **within** a single **molecule map**. **Constraints across** relation kinds (e.g. the same atom pair in **`:bonds`** and **`:dative-bonds`**) are **not** specified here.

**`:bonds` (localized).** The list **MUST NOT** contain two **`bond-entry`** values with the same **unordered** pair of atom sites (their two **`:atoms`**, as a set).

**`:dative-bonds`.** A **donor→acceptor** edge **MUST NOT** repeat — counting each donor in an entry's **`:donors`** list against that entry's **`:acceptor`**, across all entries. A donor→acceptor edge and the reverse acceptor→donor edge between the **same** two atoms also violate this rule.

**`:aromatic-systems`.** For every two distinct **`aromatic-system-entry`** values, the sets of keywords in their **`:atoms`** vectors **MUST** be disjoint. Aromatic systems **MUST NOT** share an atom.

**`:multicenter-bonds`.** The list **MUST NOT** contain two **`multicenter-bond-entry`** values with the **same** set of **`:atoms`**. Entries **MAY** share atoms — an atom **MAY** participate in several multicenter bonds (e.g. a bridging boron) — but no two **MAY** have **identical** atom sets.

**`:noncovalent-bonds`.** A **`noncovalent-bond-entry`**'s two **`:atoms`** **MUST** be distinct (no self-loop). The list **MUST NOT** contain two entries with the **same** unordered **`:atoms`** pair: **at most one** noncovalent interaction per atom pair, **regardless of kind**.

**`:stereo-atoms` / `:stereo-bonds`.** The lists **MUST NOT** contain two **`stereo-atom-entry`** values with the same **`:site`** atom, nor two **`stereo-bond-entry`** values with the same **`:site`** bond. Each atom (resp. bond) **MUST** be the site of at most one stereo element.


---

## 5. Leaf types

The **leaf grammars** in this section — the numeric **`value-expr`** and the **element**, **isotope**, **boolean**, **ring-membership**, **electron-counts**, **noncovalent-kind**, **coset**, and **relation** leaves — are the atomic value-types that the entity-string subgrammars (**§7**) compose. Each entity string references these leaves rather than redefining them.

### 5.1 Value expression

**`value-expr`** appears **only** inside a **predicate payload** (**§7.3**, **§7.4**) after **`#` *tag***. The character **`#`** **MUST NOT** appear inside a **`value-expr`** (it is reserved for starting the next predicate).

```
digit      ::= '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
nat        ::= digit+
decimal-tail ::= digit*

nat-set    ::= '{' nat (',' nat)* '}'

An empty **`nat-set`** **`{ }`** is **invalid**.

value-expr ::= '*'
             | nat-set
             | nat
             | range                   (* half-open numeric range            *)
             | '?' id                  (* variable                          *)
             | '?' id '::' nat-set     (* variable with in-domain           *)
             | bool-expr

range      ::= '(' signed-int '..' ')'   (* RangeFrom: bound <= value       *)
             | '(' '..' signed-int ')'    (* RangeTo:   value < bound        *)

bool-expr  ::= or-expr

or-expr    ::= and-expr ( '|' and-expr )*
and-expr   ::= not-expr ( '&' not-expr )*
not-expr   ::= '!' not-expr
             | '(' bool-expr ')'
             | rel-expr

rel-expr   ::= mem-expr ( rel-op mem-expr )?

mem-expr   ::= add-expr ( ( '::' | '!:' ) nat-set )?

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

**Top-level `nat-set`.** A **`value-expr`** may be **only** a **`nat-set`** (after the usual ignored whitespace between **`value-expr`** tokens, **§7.1**). It denotes a **finite numeric disjunction** for the **one** quantity fixed by the enclosing predicate tag (**§7.3**, **§7.4**): that quantity **MUST** equal one of the listed **`nat`** values. This is the same constraint **shape** as a top-level **`nat-set`** in **bond-string** **`order`** (**§7.4**) and **`element-set`** for the **element** prefix, applied at the **predicate payload** level (e.g. **`#h{1,2,3}`** with payload **`{1,2,3}`**). It **MUST NOT** introduce a numeric **`?id`**; implementations **MAY** lower it to **`bool-expr`** internally. The form ***arith* `::` *nat-set*** on **`mem-expr`** is unchanged: it constrains the **arithmetic** value on the left of **`::`**, not an implicit slot quantity by bare **`{…}`** alone.

**Top-level `range`.** A **`value-expr`** may be a half-open **`range`**: **`(i..)`** is **`RangeFrom(i)`** (admits every value **≥ `i`**), **`(..j)`** is **`RangeTo(j)`** (admits every value **< `j`**). Bounds are **`signed-int`** (so charge ranges admit negatives). The **both-bounded** form **`(i..j)`** is **not** a range — it is the finite set **`{i, …, j−1}`** and **MUST** be written as a **`nat-set`** (**`(i..i)`** is the empty set, a contradiction); restricting **`range`** to half-open keeps it always non-empty and canonical. Ranges are **anonymous** (no **`id`**), so distinct occurrences never couple — unlike a **`?id >= i`** predicate, whose variable would be shared. **`#R+`** lowers to **`RangeFrom(1)`**; **`(i..)`**/**`(..j)`** are the general forms (e.g. **`#R(1..)`**, **`#R(6)(1..)`**). A **`range`** participates in matching by **solution-set inclusion** (**§6.2**): **`RangeFrom(i)`** matches a target **`Lit(n)`** iff **`n ≥ i`**, a target set iff all members are **≥ `i`**, and a target **`RangeFrom(j)`** iff **`j ≥ i`**.

**Top-level `?` *id* and `?` *id* `::` *nat-set*.** A **`value-expr`** may be **only** a numeric **variable** **`?` *id*** — optionally with a finite in-**domain** **`?` *id* `::` *nat-set***. A bare **`?` *id*** lowers to **`ValueTerm::Var`**; the domain form lowers to a membership predicate (**`ValuePredicate::Mem`** with **`MemOp::In`**). Inside a compound expression (e.g. **`?h + 1`**, **`?h == 0`**), the same **`?` *id*** appears as **`ValueTerm::Var`** — the discriminator is whether the surrounding context is the whole value or an operand of a larger operator.

**Paren-transparency for top-level variables.** Outer parentheses around a top-level **`?` *id*** or **`?` *id* `::` *nat-set*** are **optional** and **semantically transparent**: implementations **MUST** accept the bare forms and any nesting depth of outer parens (**`(?h)`**, **`((?h))`**, **`(?h :: {1,2})`**, **`((?h :: {1,2}))`**) as identical AST. The **canonical** rendered form is **bare** (no outer parens). Disambiguation against larger expressions like **`(?h + 1)`** or **`(?a :: {0}) & 0 <= 0`** is handled by requiring a **terminator** (end-of-payload or next **`#`** predicate) after the parenthesized variable before the arm fires; otherwise the parens are interpreted as **`bool-expr`** grouping (**§5.1.1**).

**`unary-expr`** is **`sign`*** **`base-expr`**: zero or more leading **`+`** / **`-`**, then **`nat`**, **`?id`**, or parenthesized **`add-expr`**. Examples: **`#c+1`**, **`#c-2`**, **`#c--1`**. A **`sign`** with **no** following **`base-expr`** is **invalid** in the general grammar; **`#c`** additionally accepts a payload consisting **only** of **`+`** or **`-`** (after trimming whitespace) as **+1** or **−1** (**§7.3**).

**Equality** is **`==`** (not **`=`**). **Finite numeric membership** uses the **`::`** token and a **`nat-set`**: **`?h + 1 :: {2,3}`** parses as **`(?h + 1) :: {2,3}`** — the full **additive** form is built before **`::`**, which sits **below** **relations** and **logic** only (same layering as former **`in`**).

**Meaning of `::` / `!:`.** The membership operators — **`::`** (in) and **`!:`** (not-in) — appear in two shapes:

- In an **element variable** (**§5.2**), the domain form **`?` *id* `::` *element-domain*** / **`?` *id* `!:` *element-domain*** constrains the variable to **membership in** / **exclusion from** a set (or single) of element symbols (chemical **`element-literal`** values).
- In **`mem-expr`** (inside **`value-expr`**), ***arith* `::` *nat-set*** / ***arith* `!:` *nat-set*** asserts the left-hand **arithmetic** value **is** / **is not** a member of the **numeric** set. After matching, concrete values **MUST** fit the slot’s type: **`u8`** for most atom/bond numeric predicates (**§7.2**), **`i8`** for formal charge (**`#c`**), **`u32`** for isotope mass (**`#i`**).

**Ground** **`value-expr`** (predicate payloads where allowed) are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (non-empty, entries valid for the slot), with optional leading **`sign`** sequence on a bare **`nat`** / **`decimal-tail`** only: no **`?`**, **`::`**, relations, logic, **`*`**, or **parentheses** (**§5.1.4**), except **`#c`** (**§7.3**) also allows a payload consisting **solely** of **`+`** or **`-`** (**+1** / **−1**). Implementations **MAY** use a restricted parser for **Ground** (**§3**).

**Numeric evaluation.** Where **`add-expr`** / **`mult-expr`** are evaluated to **concrete** counts (e.g. **Rule** RHS), intermediate and final **numeric** results **SHOULD** be computed in a range consistent with the target slot (**§7.2**, typically **`u8`** where that slot is **`u8`**). Values **outside** the slot’s range **MUST** be rejected at validation.

**`bool-expr`** is the **`value-expr`** form for **constraints** in **Query** / **Rule** (only in slots implementations allow). A **plain numeral** payload is an **`add-expr`** that is only a **`nat`**. **Parentheses** **`(`** **`add-expr`** **`)`** group **arithmetic** only.

**Additional form (same precedence as other `add-expr` operands):** **`nat` `add-op` `?` `id`** (e.g. **`4-?h`**) is covered by **`mult-expr`** / **`add-expr`**; a **parenthesized** spelling **`(` `nat` `add-op` `?` `id` `)`** is equivalent and **MAY** be used for clarity.

Inside **`nat-set`**, optional whitespace **MAY** appear only adjacent to commas, like **element-set** (**§5.2**).

**`?id`** in **`bool-expr`** introduces or uses a **numeric** variable (**§6**). An **element variable** (**§5.2**) is **not** **`bool-expr`**. **Variables are not surface-typed**; illegal combinations are rejected when lowering / validating, not by this grammar’s token shapes.

#### 5.1.1 Precedence and parentheses (infix)

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

#### 5.1.2 Truth and repeated `id` on one atom-string

For one **atom-string** at match time, **numeric** **`?id`** values in **`bool-expr`** are fixed from the **matched** target atom so that **all** constraints on that **`id`** **hold together** (**§6**). Prefix **`!`** is classical boolean negation on its operand **`bool-expr`**. Example: **`#h(!?h==1)`** or **`#h( ! ?h == 1 )`**: **`h`** is the implicit H count; the guard holds iff **`h ≠ 1`**.

The same **`id`** may appear in **several** payloads on **one** **atom-string**; one value satisfies **every** use. Example: **`#h(?h <= 4)#v(4-?h)`** — **`h ≤ 4`** and **`v = 4 − h`**.

**Cross-atom** reuse of the same **`id`** is **not** fixed here; implementations **SHOULD** document whether and how it is allowed.

#### 5.1.3 `decimal-tail` and omitted numeral = 1

**`decimal-tail`** is **`digit`*** (**§5.1**). Its numeric meaning is **1** when there are **no** digits; otherwise the usual base‑10 value of the digit sequence (no separate **`nat`** required for the all-zero-digits case).

When a **predicate** (**§7.3**, **§7.4**) allows a **decimal-only** payload and the payload is **only** a **`decimal-tail`** (not **`*`**, not **`(`**, not **`bool-expr`**, not a **`sign`**-only **`#c`** form), the **omitted** form (zero digits after the tag) means **1**. Lexing **MUST** take the **longest** contiguous run of **`digit`**s as that numeral (**greedy**). **Special predicate** payloads (**§7.3**) are **not** **`decimal-tail`**-only forms.

| Atom tag | Omitted numeral = 1 (decimal-only payloads) |
|----------|-----------------------------------------------|
| **`#c`** | **no** — charge **MUST** be explicit (**`#c0`**, **`#c+`**, **`#c-`**, **`#c+2`**, **`#c-2`**, …); empty **`#c`** is **invalid** |
| **`#h` `#n` `#u` `#s` `#v` `#d` `#t`** | yes, when the payload is decimal-only |
| **`#a`** | yes when decimal-only; **`#a*`**, **`#a+`**, **`#a!`** are **special** (**§7.3**) |
| **`#m`** | yes when decimal-only; **`#m*`**, **`#m+`**, **`#m!`** are **special** (**§7.3**) |
| **`#i`** | yes, when the payload is decimal-only; bare **`#i`** denotes isotope mass **1** |
| **`#R`** | yes when decimal-only (incl. the sized **`#R(<size>)`** form); **`#R*`**, **`#R+`**, **`#R!`** are **special** (**§7.3**) |

Bond predicates that use **decimal-only** payloads follow the same **`decimal-tail`** rule where applicable (**§7.4**).

In **Query** and **Rule**, any predicate slot that allows a full **`value-expr`** may use **`bool-expr`**, **`*`**, top-level **`nat-set`**, **`decimal-tail`**, etc., as allowed for that tag.

#### 5.1.4 Wildcards, sets, logic, arithmetic

- The **`*`** **wildcard** **MAY** appear in **`value-expr`**, **`element`**, and **`order`**
- **`bool-expr`**: **infix** **`&` `|` `!`**, **relations**, **`::`**, **`+ - * / %`**, unary **`-`**, **`?id`**, **`nat`**, **`(`** **`add-expr`** **`)`**.

**Ground:** no **`bool-expr`** (no **`?`**, **`::`**, relations, logic), no **element variable**, no top-level negation (**`!`** *literal* or **`!`** *set*); predicate payloads are **`decimal-tail`** / **`nat`** / top-level **`nat-set`** (and tag-specific literals such as **`#i=`** for natural isotope) only where allowed.

**Query:** **`bool-expr`** where allowed; **`decimal-tail`**; **element** / **order** extensions as allowed.

**Rule:** full **`value-expr`**; **element** may use an **element variable** (**§6**).

**`<` `>` `<=` `>=` `==`** appear **only** inside **`value-expr`** (predicate payloads). **Dative** donated / accepted pair counts use predicates **`#d`** / **`#t`** (**§7.3**), not bare **`<` `>`** at the top level of the atom-string.

### 5.2 Element

The **`element`** nonterminal (**atom-string** prefix) is **literal** | **wildcard `*`** | **brace set** | **negation** | **variable** (**§5.2** grammar below). The **bond-string** **`order`** prefix (**§7.4**) is a single **`value-expr`** (**§5.1**), which **subsumes** literal **`nat`**, **`*`**, brace **`nat-set`**, **`?` *id***, **`?` *id* `::` *nat-set***, and **arithmetic** / logic (e.g. **`1+1`**, **`?o+1`**) where allowed by context.

```
element ::= element-literal
          | '*'
          | element-set
          | '!' element-literal
          | '!' element-set
          | element-var

element-set ::= '{' element-literal (',' element-literal)* '}'

element-var    ::= '?' id [ ( '::' | '!:' ) element-domain ]
element-domain ::= element-set
                 | element-literal

element-literal ::= [A-Z][a-z]*
```

- **`element-literal`**: one chemical symbol; **§7.2** (H–Og).
- **`*`**: any element; **invalid** in **Ground** unless narrowed by a containing rule outside this specification.
- **`element-set`**: finite non-empty disjunction of **one or more** **`element-literal`** entries; **§7.2**. **Query** / **Rule** when **Ground** disallows wildcards.
- **`!` *element-literal*** / **`!` *element-set***: cofinite **negation** — admits everything in the element domain **except** the named literal / set members. **§7.2** range applies to the excluded entries. **Invalid** in **Ground**.
- **`element-var`**: **Query** / **Rule** only. A **nominal** variable **`id`**, **optionally** carrying a **domain** — **membership in** (**`::`**, **`MemOp::In`**) or **exclusion from** (**`!:`**, **`MemOp::NotIn`**) an **`element-domain`** (a set or single element symbol; **§5.1**, **§6**). The operator carries the polarity; the domain itself is unnegated (unlike the literal complement **`!{…}`**, which is a field value, not a variable). With **no** domain the variable is a bare reference whose **`id`** **MUST** already be bound in rule scope (**§6**); no arithmetic on nominal variables. **Invalid** in **Ground**.

**Paren-transparency.** Outer parentheses around an **`element-var`** are **optional** and **semantically transparent**: implementations **MUST** accept the bare forms (**`?e`**, **`?e :: {C,N}`**) and any nesting depth of outer parens (**`(?e)`**, **`((?e :: {C,N}))`**, …) as identical AST. The **canonical** rendered form is **bare** (no outer parens).

Optional ASCII whitespace inside **`element-var`** around **`::`** and around commas in the inner **`element-set`**, per **§7.1**.

### 5.3 Isotope

**`#i` isotope subgrammar.** The isotope-mass slot uses its own subgrammar, not **`value-expr`**, because isotope mass numbers are tagged enum-like and have no arithmetic operations. Empty payload (bare **`#i`**) denotes mass **1** (per §5.1.3 decimal-tail).

```
isotope-payload ::= '='                           (* Natural — naturally-occurring ratios *)
                  | '*'                            (* Undetermined — wildcard           *)
                  | signed-int                     (* Lit                                *)
                  | nat-set                        (* Set — finite mass disjunction      *)
                  | '!' signed-int                 (* Not — cofinite singleton          *)
                  | '!' nat-set                    (* NotSet — cofinite multi           *)
                  | '?' id [ '::' isotope-domain ] (* Var — variable, optional in-domain *)

isotope-domain  ::= nat-set                        (* in-only; isotope variable has no not-in *)
                  | '!' signed-int                 (* MemOp::NotIn (singleton)           *)
                  | '!' nat-set                    (* MemOp::NotIn                       *)
```

**Paren-transparency.** Outer parentheses around **`?` *id*** or **`?` *id* `::` *isotope-domain*** are **optional** and **semantically transparent** (same rule as the element variable above and value-expr §5.1). Canonical render is bare.

**Undetermined** value: **`\*`** describes an undetermined value. For constraints, **`\*`** means **no constraint** and MAY be elided.

**Natural is its own channel.** **`=`** (Natural) **does not unify** with numeric variants in the lattice: **`Natural ∧ Lit(n) = ⊥`**, **`Natural ∨ Lit(n) = Undetermined`** for any **`n`**. Natural is the "no specific mass committed" state and is **disjoint** from any explicit mass number. **Ground** isotope is either **`Natural`** or a single **`Lit(n)`**.

### 5.4 Boolean

A **`boolean`** describes a truth value (*`true`* or *`false`*) with an additional undetermined variant. It is used as the payload of the aromatic constraint **`#a`** (bond, dative), the intramolecular constraint **`#I`** (noncovalent), and the trailing value of the stereo predicates **`#p`** / **`#f`**. Its compact grammar is a single trailing sentinel:

```
boolean ::= '' | '+' | '!' | '*'
```

- **`''`** (omitted) or **`'+'`** — **`true`**.
- **`'!'`** — **`false`**.
- **`'*'`** — **`Undetermined`** (no constraint); **vacuous** and **elided** from the canonical rendered string (**§7.1**), equivalent to omitting the predicate.

It lowers to **`BooleanAst`** (**`Undetermined`** | **`Lit(bool)`**); its structured EDN form is **`bool ::= true | false | :undetermined`**.

### 5.5 Ring membership

The **`ring-membership`** leaf is the payload of the **`#R`** predicate on an atom (**§7.3**), bond (**§7.4**), or dative bond (**§7.7**): a ring **`count`**, optionally scoped to rings of one **`size`**.

```
ring-membership ::= [ '(' size ')' ] count
size            ::= nat
count           ::= '*' | '+' | '!' | value-expr
```

- **`#R<count>`** bounds the **total** ring count; **`#R(<size>)<count>`** bounds the count of rings of that **`size`**. Bare (**`#R`** / **`#R(<size>)`**) means count **1** (**§5.1.3**).
- **special counts**: **`*`** = **`Undetermined`** (no constraint, **elided** on render, **§7.1**); **`+`** = **`RangeFrom(1)`** ("in at least one ring"); **`!`** = **`Lit(0)`** (acyclic, or no ring of that size).
- SMARTS parity: **`R`** → **`#R+`**, **`Rn`** → **`#Rn`**, **`R0`** → **`#R!`**, **`rn`** → **`#R(n)+`**.

It lowers to **`RingMembershipAst`** (a **`count`** value plus a **`RingScope`** of **`All`** or **`Size(n)`**); its structured EDN form is **`ring-membership-form`** (**§7.12**). The **`#R`** predicate **MAY** appear **multiple** times on one entity — one per ring scope (total and/or per-size).

**Ring enumeration parameters.** Derived ring predicates use one fixed projection; their syntax does not carry a ring-set selector or an enumeration configuration.

| Parameter | Value for derived DSL predicates |
|-----------|----------------------------------|
| **`RingModel.kind`** | **`RingSetKind::Relevant`** — the union of all minimum cycle bases |
| **`RingModel.max_ring_size`** | **22** bonds |
| **`RingConfig.relevant_cycle_algorithm`** | **`RelevantCycleEnumerationAlgorithm::Vismara`** by default |
| **`RingConfig.simple_cycle_algorithm`** | unused for this projection |

The model fields define the observable projection: **`All`** counts relevant rings of at most 22 bonds, and **`Size(n)`** counts rings of exactly **n** bonds within that projection (therefore zero when **n > 22**). The algorithm field is operational and MUST NOT change the resulting ring set. The general **`MoleculeAst::rings`** API accepts an explicit **`RingModel`** and **`RingConfig`**, but those parameters are not part of molecule DSL syntax and do not alter the meaning of atom **`#R`**, atom **`#x`**, atom **`#y`**, or localized-bond **`#R`**.

This projection is defined over the localized atom-bond graph. Dative-bond **`#R`** remains an asserted **`RingMembershipAst`** value: deriving it requires a ring model that includes dative overlays, whose semantics are not defined by this specification.

### 5.6 Electron counts

The **`electron-counts`** leaf is the mandatory leading specification of the aromatic-string (**§7.8**) and multicenter-string (**§7.9**) — the **per-atom** electron contributions, the counterpart of the bond-string's leading order.

```
electron-counts ::= '*' | '[' int (',' int)* ']'
```

- **`*`** — the whole vector **`Undetermined`** (per-atom electrons unspecified).
- **`[n,n,…]`** — a vector of concrete integers, **one per member atom** (whitespace ignored); position **`i`** is the contribution of the atom at position **`i`** of the entry's **`:atoms`** vector (**§4**). A concrete vector **MUST** have the same length as **`:atoms`**.

It lowers to **`ElectronCountsAst`**. The optional **`#e<n>`** total (aromatic / multicenter) is cross-checked against the **sum** of the electron-counts on ground inputs.

### 5.7 Noncovalent kind

The **`noncovalent-kind`** leaf is the interaction-kind field of the noncovalent-string (**§7.10**).

```
noncovalent-kind-expr    ::= noncovalent-kind-literal | '*'
noncovalent-kind-literal ::= 'Hbd' | 'Xbd' | 'Ybd' | 'Ion' | 'Vdw'
```

**Literal meanings.**

| Literal | Interaction kind |
|---------|------------------|
| **`Hbd`** | hydrogen bond |
| **`Xbd`** | halogen bond |
| **`Ybd`** | chalcogen bond |
| **`Ion`** | ionic interaction |
| **`Vdw`** | van der Waals interaction |

Each **`noncovalent-kind-literal`** is exactly three ASCII characters: one leading uppercase letter followed by two lowercase letters. The parser consumes the full three-character token; partial prefixes (**`H`**, **`Hb`**, …) **MUST** be rejected.

**`*`** admits any kind and **MUST NOT** appear in **Ground** (it is a **Query** / **Rule** form only). The kind **MUST** be either a single **`noncovalent-kind-literal`** or the wildcard **`*`**; it has no set, variable, or domain forms. It lowers to **`NoncovalentBondKindAst`** (**`Undetermined`** | **`Lit`**).

### 5.8 Coset

The **`coset`** leaf is the realized configuration index of a stereo element (**§7.11**); extended with two leading sentinels as **`config`**, it is the payload of the atom **`#T`** / bond **`#C`** inline constraints (**§7.3** / **§7.4**) and the source form behind the **`{:stereo coset-form}`** EDN constraint (**§7.12**).

```
config     ::= '*' | '!' | '+' | nat | coset-expr
coset      ::= '*' | nat | coset-expr

coset-expr ::= '~' coset-expr             (* involution operator             *)
             | '\'' coset-expr            (* mirror operator                 *)
             | coset-expr '^' cycles      (* group action by a permutation (cycles, §7.11) *)
             | nat                        (* literal coset index             *)
             | '?' id [ '::' coset-set ]  (* coset variable / domain         *)
             | coset-set                  (* literal coset set               *)

coset-set  ::= '{' nat (',' nat)* '}'
```

**Coset.** The **`coset`** is a **dense, 0-based per-class arrangement index** over the entry's ordered **`:ligands`** frame (**§4**) — the OpenSMILES arrangement order for the class, renumbered from **`0`**, **not** a Lehmer / permutation rank. For **`Th`**: **`0`** = anticlockwise (**`@`**), **`1`** = clockwise (**`@@`**). For **`Ct`**: **`0`** = **Z** (cis), **`1`** = **E** (trans). **`Ax`** / **`Sp`** / **`Tb`** / **`Oh`** follow the OpenSMILES arrangement order (0-based) for their class. **`*`** is an **undetermined** (open) coset.

**`config`** (atom **`#T`** / bond **`#C`** / **`{:stereo …}`**).** The constraint payload extends **`coset`** with two leading sentinels: **`*`** = **`Undetermined`** (no stereo constraint — equivalent to omitting the predicate), **`!`** = **`NotStereo`** (the site is **not** a stereocenter), **`+`** = **`Stereo`** with an **undetermined** coset (the site **is** a stereocenter, coset unspecified). A bare **`nat`** / **`coset-expr`** is **`Stereo`** with that coset. The EDN equivalents are **`:undetermined`** / **`:not-stereo`** / **`{:stereo :undetermined}`** / **`{:stereo coset-form}`** (**§7.12**); a **`coset-set`** serializes to the EDN vector form **`[ int+ ]`**, every other **`coset-expr`** to a **`"coset-string"`**.

**Inline ligand frame (`#T` / `#C` without `:ligands`).** An atom **`#T`** or bond **`#C`** inline coset (and the **`{:atom [i {:tetrahedral-stereo …}]}`** / **`{:bond [i {:cis-trans-stereo …}]}`** EDN forms) carries **no** **`:ligands`** vector, so its index is numbered against an **implicit frame derived from the molecular graph**:

- **`#T`** (tetrahedral, atom site): the atom's neighbor atoms in **ascending atom-index** order, then — when there are **exactly three** explicit neighbors — **one** virtual ligand (an implicit hydrogen or a lone pair) appended **last**. The site **MUST** present **3 or 4** ligands this way (three explicit + one virtual, or four explicit); any other count is invalid.
- **`#C`** (cis/trans, double-bond site): for **each** double-bond terminus in turn, that terminus's substituents (its neighbors other than the far terminus) in **ascending atom-index** order, then — when a side has a **single** explicit substituent — **one** virtual ligand appended **last within that side's group**. Each side **MUST** present **at least one** explicit substituent.

A **`stereo-atom-entry`** / **`stereo-bond-entry`** **`:ligands`** vector **overrides** this implicit frame; the two coincide when **`:ligands`** lists the same neighbors in the same order.

**Coset operators (reserved).** The **`~`** (involution) and **`^`*cycles*** (group action by a permutation in 0-indexed disjoint-cycle notation, **§7.11**) operators, and the coset variable / set forms (**`?id`**, **`?id :: {…}`**, **`{…}`**), **parse** as **`coset-expr`** and **round-trip**, but their **matching** semantics are **staged** — relative-stereo binding and non-tetrahedral coset domains land incrementally. Only **ground literal cosets** (and the **`*`** / **`!`** / **`+`** sentinels) are presently matched; a conforming matcher **MAY** reject an operator / variable **`coset-expr`** until the corresponding stage lands.

### 5.9 Relation

The **`relation`** leaf is the value of the stereo topicity **`#o`** and stereogenicity **`#g`** predicates (**§7.11**), over the 3-glyph topicity / stereogenicity domain.

```
relation  ::= '*' | ['!'] ( glyph | glyph-set )    (* * undetermined; bare glyph = singleton; a set = members; ! = complement *)
glyph-set ::= '{' glyph (',' glyph)* '}'
glyph     ::= '=' | '\'' | '/'
```

**Relation forms.** A **`relation`** has four surface forms, each **faithful to its stored variant** (representation, not canonicalization): **`*`** (**`Undetermined`** — the full domain); a bare **`glyph`** (a **`Lit`** singleton, e.g. **`=`**); a **`glyph-set`** **`{a,b,…}`** (an explicit **`LitSet`** of members, e.g. **`{=,'}`**); or a leading **`!`** on a glyph or set (a **`NotSet`** — the **complement** of the named member(s): **`!/`** = not diastereotopic, **`!{=,'}`** = neither homotopic nor enantiotopic). This mirrors the EDN, which distinguishes the member vector **`[a b]`** (**`LitSet`**) from the complement **`{:not-in [x]}`** (**`NotSet`**, **§7.12**). Over the 3-element topicity / stereogenicity domain every non-empty subset is expressible several ways; **canonicalization** (a **separate** pass) reduces a set to the smaller of its positive / complement side — a 2-set to **`!x`**, a 1-set to a bare glyph — but the surface **round-trips whichever variant the AST holds**. A full-domain (**`Undetermined`**) relation is a **vacuous** predicate: like the atom **`#a*`** / **`#T*`** special forms (**§7.3**) it is **admissible on parse** but **elided** from the canonical rendered string (**§7.1**) — **`#o*`** / **`#g*`** are dropped on render, equivalent to omitting the predicate.

---

## 6. Match semantics and bindings

**Ground molecule, pattern LHS.** The target is **ground** (fully instantiated). The **LHS** of a rule (or query) may still contain **wildcards**, **sets**, and **binds**: that is **pattern** data, not an indeterminate molecule.

### 6.1 Inherent fields and derived predicates

**Inherent fields.** Each AST form — atom, localized bond, aromatic system, multicenter bond, dative bond, noncovalent bond, stereo atom, stereo bond — carries a fixed set of **inherent fields**. An inherent field's value **identifies** the entity at that slot. An entity is **ground** iff every inherent field holds a single concrete value (a literal; not a wildcard, set, bind, ref, or unresolved symbolic state). Nothing else affects grounding.

| Form | Inherent fields |
|------|-----------------|
| atom | element, isotope mass (**`#i`**), charge (**`#c`**), implicit hydrogens (**`#h`**), lone pairs (**`#n`**), unpaired-electron count (**`#u`**) and multiplicity (**`#s`**) |
| localized bond | order, charge (**`#c`**), unpaired-electron count (**`#u`**) and multiplicity (**`#s`**) |
| aromatic system | charge (**`#c`**), unpaired-electron count (**`#u`**) and multiplicity (**`#s`**), π-electron count (**`#e`**) |
| multicenter bond | charge (**`#c`**), unpaired-electron count (**`#u`**) and multiplicity (**`#s`**), electron count (**`#e`**) |
| dative bond | a single **`:acceptor`** and its one-or-more **`:donors`** atoms — the assignment on the map entry (**§4**) — plus the leading **`order`** token of the dative-string (number of donated electron pairs; **§7.7**). |
| noncovalent bond | interaction kind (**`Hbd`**, **`Xbd`**, **`Ybd`**, **`Ion`**, **`Vdw`**) |
| stereo atom | coordination **`class`** (geometry) and **`coset`** configuration index (the **`:type`** payload, **§7.11**). The bearing **`:site`** atom and ordered **`:ligands`** frame (**§4**) are the relation's participants, not payload fields. |
| stereo bond | cis/trans **`class`** and **`coset`** configuration (the **`:type`** payload, **§7.11**). The bearing **`:site`** bond and ordered **`:ligands`** frame (**§4**) are the relation's participants. |

**Derived predicates.** Every predicate admitted in the DSL that is not an inherent field is a **derived predicate** — a topological query evaluated against the target graph once an embedding is proposed. This includes per-atom **`#D`**, **`#X`**, **`#V`**, **`#x`**, **`#y`**, **`#H`**, **`#R`** (**§7.3**); the bond-namespace **`#R`**; per-aromatic, per-multicenter, per-dative ring-membership predicates; and the molecule-wide entries of **§7.12**. Derived predicates **filter** matches; they do **not** carry identity and **do not** affect grounding. Adding a derived predicate — even a wildcard-valued one — to a pattern never makes a ground target stop being ground.

**Symmetry-derived stereo predicates.** The stereo entity predicates — **`#p`** ligand symmetry, **`#o`** topicity, **`#g`** stereogenicity (**§7.11** / **§7.12**) — are **derived** from the resolved molecule's **graph automorphisms** (the ligand-frame symmetry group of the stereo element), not from the local string. As derived predicates they **filter** matches and **do not** affect grounding: a stereo element is ground iff its **`class`** + **`coset`** are concrete (**§6.1** table), regardless of which **`#p`**/**`#o`**/**`#g`** assertions it carries. The validator computes the molecule-wide symmetry once on the **resolved** AST and **cross-checks** the derived value against each stored constraint — when both are ground and inconsistent (including a kind/degree mismatch, or a **`'`** value on an achiral class), the molecule is rejected — exactly as the topology-derived fields are cross-checked against the stored inherent fields. **`#f`** fluxionality is a stored dynamical assertion (not derivable from a static graph); it is matched as a stored predicate, not cross-checked.

### 6.2 Pattern–target match

**Match as solution-set inclusion.** Each attribute slot has a **solution set** — the set of ground values the slot admits. A **literal** (e.g. **`C`**, **`3`**, **`+1`**) admits exactly itself; a **set** (**`{C,N}`**, top-level **`nat-set`**) admits its members; a **negation** (**`!H`**, **`!12`**) admits everything in the slot's value domain *except* the named literal; a **negative set** (**`!{F,Cl}`**, **`!{12,13}`**) admits the complement of the listed entries; a **wildcard** (**`*`**) admits everything in the slot's value domain; a **`bool-expr`** admits every value for which the expression holds (**§5.1**); a **`range`** admits its half-line (**`(i..)`** every value **≥ `i`**, **`(..j)`** every value **< `j`**, **§5.1**); a **special-symbolic** payload (**`#i=`**, **`#a*`**, **`#a+`**, **`#a!`**, **`#m*`**, **`#m+`**, **`#m!`**) admits only its named symbolic state (**§7.3**). The **`#R`** family is ordinary numeric — **`#R*`** = wildcard, **`#R+`** = **`RangeFrom(1)`**, **`#R!`** = **`Lit(0)`**, **`#Rn`** = **`Lit(n)`** — matched by the wildcard / range / literal rules. For a given slot, the **pattern** matches the **target** iff `solution-set(pattern)` ⊇ `solution-set(target)` — the pattern admits every value the target admits. Match is **not** symmetric.

| pattern kind | target kind | matches iff |
|--------------|-------------|-------------|
| wildcard (**`*`**) | any | always |
| non-wildcard | wildcard | never (target's set is strictly larger) |
| literal | literal | values equal |
| literal | set | set is exactly that singleton |
| set | literal | literal is a set member |
| set **P** | set **T** | **T ⊆ P** |
| **`bool-expr`** | literal | expression holds on the literal (**§5.1**) |
| **`bool-expr`** | set | expression holds on **every** set member |
| **`bool-expr`** | wildcard / **`bool-expr`** | **undefined** in general; implementations **MAY** reject |
| special-symbolic **s** | special-symbolic **t** | **s == t** |
| special-symbolic | literal / set | never (disjoint domains) |

**Element matching.** Parallels the above: **`element-literal`** against **`element-literal`** by equality; **`element-set`** against a target iff the target's admissible symbols are a subset; an **element variable** with a domain behaves as its inner **`element-set`** for the match (the nominal binding is a side effect, not a match filter, **§6.3**); a bare **element variable** outside a resolved rule-scope binding context matches nothing.

**Noncovalent-kind matching.** A wildcard (**`*`**) matches any kind; a literal matches by equality over the five-literal domain **`{Hbd, Xbd, Ybd, Ion, Vdw}`**. There are **no** set / bind / ref forms.

**Molecule-level match.** A molecule-map pattern matches a target iff (a) every **`atom-string`** matches its corresponding target atom-string field-wise — element and every inherent-field predicate payload, per the rules above; (b) every **`bond-string`** matches field-wise; (c) each structural relation (**`:aromatic-systems`**, **`:multicenter-bonds`**, **`:dative-bonds`**, **`:noncovalent-bonds`**, **`:stereo-atoms`**, **`:stereo-bonds`**) matches per its own inherent fields; (d) every derived predicate holds on the resulting embedding (**§6.1**). Any failure rejects the embedding.

### 6.3 Bindings

**One binding per match.** For a **fixed** embedding of the LHS pattern into the target (one way of mapping pattern sites to concrete atoms/bonds that satisfies all constraints), each **`id`** introduced by an **element variable** or by **`?id`** in **`bool-expr`** / **`value-expr`** has **exactly one** value in that match. There is no separate “CSP over the whole molecule without choosing an embedding”: the engine first chooses an embedding (or enumerates them — see below), then **the match binding** is fixed — i.e. the mapping from each such **`id`** to its concrete value for that embedding (numeric for **`?id`**, element symbol for nominal variables).

**Multiple results from one ground target.** Ambiguity does **not** require an indeterminate target. The **same** ground molecule can admit **several** distinct **embeddings** of the **same** LHS (e.g. two equivalent substituents). Each embedding yields its own **match binding**. Whether the rule **fires once**, **once per embedding**, or **aggregates** products is **policy** for the rule evaluator, not fixed by this specification.

**Nominal vs numeric.** An **element variable** carries **element**-symbol values. **`?id`** in **`bool-expr`** carries **numeric** bind / use for that attribute. Arithmetic applies only to **numeric** **`id`** values. **Nominal** **`id`** may be **re-used** on the RHS via an **element variable** (no domain) with the same name; no arithmetic on those.

**Identifier scope.** On **one** **atom-string**, the same numeric **`id`** may appear in **multiple** predicate payloads; all uses denote **one** value and are **jointly** satisfied (**§5.1.2**). Whether **`id`** may also be shared **across** atom-strings on a rule LHS (or RHS) is **not** fixed here; implementations **SHOULD** document cross-atom **`id`** rules. **Order** of **predicates** (**§7.3**) is arbitrary; evaluation **MUST** treat all constraints on **`id`** as **simultaneous**, not sequential by textual order.

---

## 7. Subgrammars

### 7.1 Whitespace and `#`

- **ASCII whitespace** (space, tab, CR, LF) is **not** significant: it **MAY** appear between **tokens** and is **ignored**, except where this section **forbids** it.
- **Leading and trailing** whitespace on the whole **atom-string** or **bond-string** is ignored.

**Whitespace is forbidden:**

- Between **`#`** and the **tag letter** of a **predicate** (**§7.3**, **§7.4**): **`#h`** is valid; **`# h`** is **invalid**.
- Inside multi-character operators: **`<=`**, **`>=`**, **`==`**, **`::`** (**§5.1**).
- Between **`?`** and the first character of **`id`** in **`?id`** (numeric variable) and in an **element variable**.

**`#` (U+0023).**

- **Inside** an **atom-string** or **bond-string**, **`#`** **MUST** appear **only** as the first character of a **predicate** (**`#` *tag***). **`#`** **MUST NOT** appear inside **`element`**, **`order`**, or inside any **`value-expr`** / **payload** substring.

**Payload extraction.** A **predicate** is **`#`**, one **tag** character **`[A-Za-z_]`**, and a **payload** consisting of all following characters up to (but not including) the **next** **`#`** or **end of string**, after **whitespace normalization** for the purpose of **tokenizing** the payload as **`value-expr`**: the payload text **MAY** contain ignored whitespace between **`value-expr`** tokens as in **§5.1**. The **payload** **MUST NOT** contain **`#`**.

**Examples (atom):** **`C`**, **`C#h3`**, **`C#h*`**, **`!H`**, **`!{F,Cl}`**, **`?e`**, **`?e :: {Cl,Br}`**, **`?e !: {F,Cl}`**, **`C#a*`**, **`C#a !`**, **`C#c+`**, **`C#c-`**, **`C#c +`**.

- A **`nat`** and an **`id`** contain **no** internal whitespace.
- A **relational** token is **`<=`**, **`>=`**, **`==`**, or a **single** **`<`** or **`>`** that is **not** part of **`<=` `>=`**. **Multi-character** tokens are one lexical unit. Plain **`=`** is **not** a relational operator.
- An **arithmetic** token is **`+`**, **`-`**, **`*`**, **`/`**, **`%`**. Leading **`+`** / **`-`** on a **`base-expr`** are **`sign`** tokens (**§5.1**); binary **`+`** / **`-`** appear between **`mult-expr`** operands.
- **Inside** an **element** or **order** **brace set** `{…}`, optional whitespace **MAY** appear **only** immediately before or after a comma separating entries. No whitespace inside an **`element-literal`** or **`order-entry`** (**`nat`**).

**Vacuous-payload elision (canonical rendering).** Implementations **MAY** elide vacuous payloads — predicates whose value is **`Undetermined`** (`*` in surface form), and inherent fields with prefixed tags (**`#c`**, **`#u`**, **`#s`**, **`#e`**, …) whose value is **`Undetermined`** — from the canonical rendered form. Both forms remain admissible **on parse**, so a renderer that elides them still accepts a string in which they appear explicitly. **Leading unprefixed** inherent fields (bond **order**, atom **element**, noncovalent bond **type**) are **exempt** from this elision because they fix the entity-string's start position; for these, **`Undetermined`** **MUST** render as **`*`**. Round-trip identity at the AST level therefore holds only for ASTs whose constraint and inherent-field payloads are non-vacuous (or whose vacuous payloads sit on a leading unprefixed field).

### 7.2 Numerical limits

**Chemical elements.** Any **`element-literal`**, any entry in an **`element-set`**, and any **nominal** **element variable** (**`element-var`**) **MUST** refer only to elements from **hydrogen** (**H**) through **oganesson** (**Og**). Implementations **MUST** reject symbols outside that range.

**Charges.** **Formal charge** on atoms (**`#c`**), **formal bond charge** (**`#c`** on **bond-string**), **aromatic-system charge** (**`#c`** on **aromatic-string**, **§7.8**), and atom-subset charge sums **`{:charge-sum {:atoms [...] :sum n}}`** (**§7.12**) where integral **MUST** fit a **signed 8-bit** integer (**−128…127**). The **`#c`** payload is a **`value-expr`** that evaluates to the signed charge, including the **special** forms **`+`** / **`-`** for **±1** (**§7.3**), e.g. **`#c2`**, **`#c-2`**, **`#c+`**, **`#c-`**.

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

**Bond-string** **`#u`** / **`#s`** / **`#c`** use the **bond** namespace (**§7.4**); meanings parallel **unpaired electrons**, **multiplicity**, and **bond formal charge**.

**Aromatic-string** **`#c`** / **`#u`** / **`#s`** / **`#e`** use the **aromatic-system** namespace (**§7.8**); **`#e`** denotes the total π-electron count (**`u8`**), other tags parallel the bond namespace.

**Lexical** **`nat`** in the grammar is unbounded; **Ground** validation **MUST** reject values outside the **u8** (or **i8** / **u32** as above) range for the corresponding slot.

**Bond order** (**§7.5**) uses a **discrete** model; **fractional** bond orders **MUST NOT** appear in the **bond-string**. **Aromatic** connectivity **MUST NOT** be encoded as a bond **order**; use the molecule map’s **`:aromatic-systems`** section (**§4**) and ordinary **`:bonds`** entries.

### 7.3 Atom subgrammar

```
atom-string ::= element atom-predicate*

atom-predicate ::= '#' tag payload

tag ::= [A-Za-z_]
```

- **`element`** is first (**§5.2**).
- **Zero or more** **`atom-predicate`** units follow. **Optional** ASCII whitespace **MAY** appear **between** **`element`** and the first **`#`**, and **between** successive predicates.
- **At most one** predicate per **tag letter** per **`atom-string`** (each row of the table below is a **kind**).
- **Canonical order** of predicates after **`element`** (stable serialized form): **`#i`**, **`#c`**, **`#h`**, **`#n`**, **`#u`**, **`#s`**, **`#v`**, **`#d`**, **`#t`**, **`#a`**, **`#m`**, **`#T`**, then the derived predicates **`#D`**, **`#X`**, **`#V`**, **`#x`**, **`#y`**, **`#H`**, **`#R`**. Implementations **MAY** specify further ordering for fields not listed here.

**`payload` parsing.** After trimming leading / trailing whitespace on the **payload** substring, parse as follows:

1. **`#c`**: if the trimmed payload is **exactly** **`+`** or **`-`**, the formal charge is **+1** or **−1** (same meaning as **`#c+1`** / **`#c-1`**). Otherwise parse as **`value-expr`** (**§5.1**) (or the **Ground** subset in **§5.1.4**).
2. **`#i`**: parsed by a **dedicated isotope subgrammar** (**§5.3**, after **`element`**), not as **`value-expr`**.
3. **Any other tag**: parse the payload as **`value-expr`** (or **Ground** subset) unless the payload matches a **special** form below.

**Special predicate payloads**:

| Form | Tag | Meaning |
|------|-----|---------|
| **`=`** | **`#i`** | **Natural isotope**: the **naturally-occurring isotopic ratios** of **`element`**, following OpenSMILES. This is the default / expected state for each element. It is **not** a mass number: `Natural` is its own channel and is disjoint from every `Lit(n)` (**§5.3**). |
| **`+`** | **`#a`** | Aromatic, **unspecified** aromatic valence. Note that **`#a0`** is a valid aromatic valence (aromatic boron, tropylium). |
| **`!`** | **`#a`** | Not aromatic is **not** a member of any aromatic system. |
| **`+`** | **`#m`** | Multicenter, **unspecified** multicenter valence. Note that **`#m0`** is valid multicenter valence (B2H6). |
| **`!`** | **`#m`** | Not in multicenter bond. |
| **`+`** | **`#R`** | Atom is in **at least one** ring. With size, **`#R(<size>)+`** means at least one ring of that size (SMARTS **`r<size>`**). |
| **`!`** | **`#R`** | Atom is in not in a ring (acyclic). With a size, **`#R(<size>)!`** means no ring of that size. |
| **`!`** | **`#T`** | Atom is **not** a tetrahedral stereocenter (**`NotStereo`**). |
| **`+`** | **`#T`** | Atom **is** a tetrahedral stereocenter with an **unspecified** coset. |
| **`+`** / **`-`** | **`#c`** | **+1** / **−1** formal charge (**§7.3** above). |

Other **`#h`** / **`#a`** / **`#m`** payloads use the usual **`value-expr`** / **`decimal-tail`** rules (**§5.1**, **§5.1.3**).

| Tag | Meaning |
|-----|---------|
| **`#i`** | Isotope mass; **special** **`#i=`** (natural isotope, **§5.3**) |
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
| **`#T`** | **Tetrahedral stereo** configuration at the atom (SMARTS-style stereo query). Payload is a **`config`** (**§5.8**): **special** **`#T*`** / **`#T!`** / **`#T+`**, a coset literal **`#T<n>`** (e.g. **`#T1`**, **`#T2`**), or a coset operator-expression. Canonical constraint form **`{:atom [i {:tetrahedral-stereo …}]}`** (**§7.12**). |
| **`#D`** | **Degree**: number of neighbors in the molecular graph (SMARTS `D`). Derived predicate evaluated against the target; **not** a ground atom field. |
| **`#X`** | **Total degree** (connectivity): degree plus implicit-H count (SMARTS `X`). Derived. |
| **`#V`** | **Total valence**: localized valence plus implicit hydrogens, aromatic valence, and multicenter valence. Derived. |
| **`#x`** | **Ring degree** (ring connectivity): number of ring bonds at the atom (SMARTS `x`). Derived. |
| **`#y`** | **Ring valence**: sum of **bond orders** of the atom's **ring** bonds. Derived. |
| **`#H`** | **Total hydrogens**: implicit H count plus explicit H neighbors (SMARTS `H`). Derived. |
| **`#R`** | **Ring membership**. **`#R<count>`** bounds the **total** ring count; **`#R(<size>)<count>`** bounds the count of rings of that **size**. Each count follows the **§5.1.3** omitted-numeral convention: bare **`#R`** / **`#R(<size>)`** means count **1**; **`#R<n>`** / **`#R(<size>)<n>`** means exactly **n**. **Special** **`#R*`** (no constraint), **`#R+`** (the range **`RangeFrom(1)`**, "in at least one ring"), **`#R!`** (count **0**). SMARTS parity: **`R`** → **`#R+`**, **`Rn`** → **`#Rn`**, **`R0`** → **`#R!`**, **`rn`** → **`#R(n)+`**. Derived. |

### 7.4 Bond subgrammar

**Bond-string** uses a **separate** namespace from **atom-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on bonds (**§7.2**). The **`order`** prefix is a **`value-expr`** (**§5.1**); see **§5.2** for the parallel with **`element`**.

```
bond-string ::= order bond-predicate*

bond-predicate ::= '#' tag payload

order ::= value-expr
```

**Segmentation.** Let **`order-text`** be the substring of the **`bond-string`** from the first character after any leading whitespace up to (but not including) the first **`#`** that starts a **`bond-predicate`** (**`#`** immediately followed by a **tag** letter, with **no** whitespace between **`#`** and the tag), or the whole **trimmed** string if there is no such **`#`**. **`order-text`** **MUST** be parsed as **`value-expr`** (**§5.1**).

**Bond predicates.** **Zero or more** **`bond-predicate`** units follow **`order`**. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`a`**, **`C`**; **`#R`** **MAY** appear **multiple** times (one per ring scope — total and/or per-size). **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#a`**, **`#R`** (total first, then by size ascending), **`#C`**.

**Whitespace** between **`#`** and the **tag** letter is **invalid** (**§7.1**).

**`#c` (bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5.1**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.)

**`#u`** / **`#s`.** After **`#u`** or **`#s`**, parse a **`value-expr`** (**§5.1**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.1.3** for decimal-only slots). **No** extra lookahead is required beyond **`value-expr`** termination and the next predicate or end of string.

**`#R` (bond ring membership).** Same forms as atom-level **`#R`** (**§7.3**): **`#R<count>`** (total ring count) or **`#R(<size>)<count>`** (count of rings of that size); bare means **1**; **`#R*`** means no constraint; **`#R+`** is the range **`RangeFrom(1)`** ("bond lies in at least one ring"); **`#R!`** means count **0**.

| Tag | Meaning (bond namespace) |
|-----|---------------------------|
| **`#c`** | Bond formal charge (**`i8`**, **§7.2**) |
| **`#u`** | Unpaired electrons (bond centered); **`u8`** |
| **`#s`** | Spin multiplicity (2S+1) (bond centered); **`u8`** |
| **`#a`** | **Aromatic** membership; **`bool`** |
| **`#R`** | **Ring membership**: **`#R<count>`** gives the **total** ring count, **`#R(<size>)<count>`** the count of rings of that **size**. Omitted-numeral convention (**§5.1.3**); Derived. |
| **`#C`** | **Cis/trans stereo** configuration at the bond (SMARTS-style stereo query). Payload is a **`config`** (**§5.8**): **special** **`#C*`** / **`#C!`** / **`#C+`**, a coset literal **`#C<n>`** (e.g. **`#C1`**, **`#C2`**), or a coset operator-expression. Canonical constraint form **`{:bond [i {:cis-trans-stereo …}]}`** (**§7.12**). |

**Special predicate payloads**:

| Form | Tag | Meaning |
|------|-----|---------|
| **`!`** | **`#a`** | Bond **not** aromatic. |
| **`+`** | **`#a`** | Bond aromatic. Equivalent to **`#a`**. |
| **`!`** | **`#R`** | Bond **not** in a ring. When used with a ring size specification, **`#R(6)!`**, bond not in a ring of **that** size. |
| **`+`** | **`#R`** | Bond in **at least one** ring. When used with a ring size specification, **`#R(6)+`**, bond is at least one ring of **that** size. |
| **`!`** | **`#C`** | Bond **not** cis/trans stereogenic. |
| **`+`** | **`#C`** | Bond **is** cis/trans stereogenic with an **unspecified** coset. |
| **`*`** | **`#C`** | Cis/trans coset **undetermined** — no constraint; **elided** on render. |

**Bond order values** in the **bond-string** **MUST NOT** be **fractional** after evaluation (**§7.5**). **Aromatic** bond **order** as a distinct category **MUST NOT** be used in **`order`**; use **§4** instead.

### 7.5 Bond order

**Semantic model** for **localized** bond order in the **bond-string** (the **`order`** nonterminal is **`value-expr`**, **§7.4**):

- **Discrete** orders **1**, **2**, **3**, and **4** after any **arithmetic** and binding.
- **Any** order: **`*`** as **`value-expr`** (**Query** / **Rule**).
- **Finite set**: top-level **`nat-set`** in **`value-expr`** (e.g. **`{1,2,3}`** or **`{2}`**).
- **Arithmetic and constraints**: full **`value-expr`** on **`order`**, including **`add-expr`**, **`::`** **`nat-set`**, **`bool-expr`**, and **`?id`** binds, subject to **Ground** restrictions below.

In **Ground**, **`order-text`** **MUST** denote a single definite order in **{1,2,3,4}**: **`*`**, a top-level **`nat-set`** whose entries are **only** **1**–**4**, **`(?` *id* `::` *set* `)`** / **`(?` *id* `)`** only where the implementation resolves them to one value, or **`value-expr`** that is **only** **`sign`*** **`nat`** with value **1**–**4** (no **`?`**, **`::`**, relations, logic, or **`(`** … **`)`**). **Query** / **Rule** **MAY** use the full **`value-expr`** grammar on **`order`**.

This section does **not** define **`bond-keyword`** shorthands; see **§7.6**.

### 7.6 Bond and atom literals

**Bond entry shorthands.** A **`bond-keyword`** as the **`:type`** value of a **`bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent bond-string payload. Normative expansion table:

| Keyword | Expands to | Bond order |
|---------|-----------|------------|
| **`:single`** | **`"1"`** | 1 |
| **`:double`** | **`"2"`** | 2 |
| **`:triple`** | **`"3"`** | 3 |
| **`:quadruple`** | **`"4"`** | 4 |
| **`:aromatic`** | **`"1#a"`** | 1 with an aromatic-participation constraint |

Implementations **MUST** accept these five keywords wherever **`bond-spec`** is expected. No other **`bond-keyword`** values are defined; unrecognized keywords **MUST** be rejected.

**Dative entry shorthands.** A **`dative-keyword`** as the **`:type`** value of a **`dative-bond-entry`** (**§4**) is a fixed **EDN keyword** that expands to an equivalent dative-string payload. Normative expansion table:

| Keyword | Expands to | Pairs donated | Example |
|---------|-----------|---------------|---------|
| **`:single`** | **`"1"`** | 1 | NH₃→BF₃ |
| **`:double`** | **`"2"`** | 2 | Ni(C₄H₄)₂ |
| **`:triple`** | **`"3"`** | 3 | |
| **`:quadruple`** | **`"4"`** | 4 | uranocene U(C₈H₈)₂ |

Implementations **MUST** accept these four keywords wherever **`dative-bond-spec`** is expected. Higher pair counts and any non-trivial dative payload **MUST** use the **`dative-string`** form (**§7.7**); unrecognized **`dative-keyword`** values **MUST** be rejected.

**Atom literals.** Atom literals are **EDN strings** whose contents are **atom-string** payloads (**§7.3** / **§5.2**). Keyword-shaped atom shorthands (via **`:atom-aliases`**) are defined in **§4**.

### 7.7 Dative-bond subgrammar

**Dative-string** carries the **bond order** (number of donated electron pairs) and optional **aromatic** and **ring-membership** constraints on a single **`dative-bond-entry`** (**§4**). The grammar parallels **bond-string** (**§7.4**): a leading **`order`** token followed by zero or more **`#…`** predicates. The dative-string has **no** inherent-field tags beyond order and **no** direction token; direction is expressed entirely by the **`:donors`** / **`:acceptor`** assignment on the containing entry.

```
dative-string ::= order dative-predicate*

order            ::= value-expr | '*'
dative-predicate ::= '#' tag payload
```

**Order.** The leading **`order`** token is a **`value-expr`** (**§5.1**) — typically a positive integer literal — that records how many electron pairs are donated. **`*`** means **`Undetermined`**. The **`dative-keyword`** shorthands (**§7.6**) — **`:single`**, **`:double`**, **`:triple`**, **`:quadruple`** — expand to the literal forms **`"1"`**, **`"2"`**, **`"3"`**, **`"4"`**.

**Dative predicates.** **Zero or more** **`dative-predicate`** units after the order token. **Optional** ASCII whitespace **MAY** appear between the order and the first **`#`**, and between successive predicates. **At most one** **`#a`**; **`#R`** **MAY** appear **multiple** times (one per ring scope — total and/or per-size). **Canonical** predicate order (stable serialization): order, then **`#a`**, then **`#R`** (total first, then by size ascending).

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#a` (dative-bond aromatic constraint).** A **boolean** constraint asserting whether the dative bond participates in an aromatic system. **`#a`** / **`#a+`** assert it **does** (**`true`**); **`#a!`** asserts it does **not** (**`false`**); **`#a*`** is **`undetermined`** — no constraint, **elided** on render. Examples of the **`true`** case: the N→B π-donation of borazine, O→B of boroxine, or a C→M coordination spanning a metallaaromatic ring. The semantics parallel the bond-namespace **`#a`** of **§7.4**; aromatic-ring perception cross-checks the constraint against actual ring membership.

**`#R` (dative-bond ring membership).** Same forms as the atom-level and bond-level **`#R`** (**§7.3**, **§7.4**): **`#R<count>`** (total ring count) or **`#R(<size>)<count>`** (count of rings of that size); bare means **1**; **`#R*`** means no constraint; **`#R+`** is the range **`RangeFrom(1)`** ("dative bond lies in at least one ring"); **`#R!`** means count **0**.

| Tag | Meaning (dative-bond namespace) | Storage |
|-----|-----------------------------------|----------|
| (leading) | **Order**: number of donated electron pairs (**`u8`**, **§7.2**) | inherent field |
| **`#a`** | **Aromatic**: boolean constraint — the dative bond **is** (**`#a`** / **`#a+`**) / **is not** (**`#a!`**) part of an aromatic system; **`#a*`** = **`undetermined`**. | boolean constraint |
| **`#R`** | **Ring membership**: **`#R<count>`** (total) / **`#R(<size>)<count>`** (per size); **special** **`#R*`** / **`#R+`** / **`#R!`** (**§7.3**). | asserted constraint; topology derivation deferred |

**Direction.** Dative bonds are intrinsically directional. Direction is carried entirely by the ordered **`:donors`** / **`:acceptor`** assignment on the containing **`dative-bond-entry`** (**§4**); the dative-string itself has **no** direction token. Under pattern matching (**§6**), the embedding MUST map pattern **`:donors`** atoms to target donors and the pattern **`:acceptor`** to the target acceptor — a donor/acceptor swap across the embedding rejects the match.

**Donor / acceptor / cross-bond references.** Donor-side and acceptor-side constraints on the endpoint atoms (equivalent to the **`:donated-pairs`** / **`:accepted-pairs`** atom-constraint forms of **§7.12**, or atom-string **`#d`** / **`#t`** of **§7.3** pinned to one endpoint) attach via the molecule-wide **`:constraints`** section; they **MUST NOT** be encoded inside the dative-string. The same holds for the "parallels another bond" relation and any reference to other molecule-level entities.

### 7.8 Aromatic system subgrammar

**Aromatic-string** uses a **separate** namespace from **atom-string** and **bond-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on aromatic systems (**§7.2**). It carries **per-aromatic-system** state — overall **charge** (**`#c`**), unpaired-electron count (**`#u`**), and multiplicity (**`#s`**) as inherent fields, and an **optional asserted total π-electron count** (**`#e<n>`**) as an inline constraint — as the **`:type`** value of an **`aromatic-system-entry`** (**§4**). The **per-atom** π contributions are the **mandatory leading `electron-counts`** of this string (below), not a separate map key.

```
aromatic-string ::= electron-counts aromatic-predicate*

aromatic-predicate ::= '#' tag payload
```

**Electron counts.** The string **MUST** begin with the **`electron-counts`** leaf (**§5.6**) — the per-atom π contributions.

**Aromatic predicates.** **Zero or more** **`aromatic-predicate`** units following the **`electron-counts`**. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (aromatic-system formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5.1**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. (So e.g. **`#c+2`** is charge **+2** via **`value-expr`**, not **`#c+`** followed by junk.) Same convention as atom (**§7.3**) and bond (**§7.4**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5.1**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.1.3** for decimal-only slots). **`#e`** omitted means **1** π-electron.

| Tag | Meaning (aromatic-system namespace) | Storage |
|-----|---------------------------------------|----------|
| **`#c`** | Aromatic-system formal charge (**`i8`**, **§7.2**) | inherent field |
| **`#u`** | Unpaired electrons (system-centered); **`u8`** | inherent field |
| **`#s`** | Spin multiplicity (2S+1) (system-centered); **`u8`** | inherent field |
| **`#e`** | Asserted total π-electron count; **`u8`** | inline constraint (`AromaticSystemConstraint::ElectronCount`) |

**`#e<n>` semantics.** **`#e<n>`** asserts the system's total π-electron count and parses to an inline aromatic-system constraint (`AromaticSystemConstraint::ElectronCount(n)`) on the entry's constraint store, **not** to a direct field. The per-atom contributions in the string's leading **`electron-counts`** are the canonical data; **`#e<n>`** is the optional total assertion that downstream validation cross-checks against their **sum** on ground inputs. **`#e`** is omitted from the canonical entity-string form when no `ElectronCount` constraint is present.

**No canonical-constraint equivalent for charge / unpaired-electron state.** Aromatic-system charge (**`#c`**), unpaired electrons (**`#u`**), and spin multiplicity (**`#s`**) live as direct fields on the aromatic-system entity (set by the aromatic-string predicates above) and have **no** canonical **`:constraints`** form.

### 7.9 Multicenter-bond subgrammar

**Multicenter-string** uses a **separate** namespace from **atom-string**, **bond-string**, and **aromatic-string**: the **same** **`tag`** letter **MAY** denote a **different** meaning on multicenter bonds (**§7.2**). It carries **per-multicenter-bond** state — overall **charge** (**`#c`**), unpaired-electron count (**`#u`**), and multiplicity (**`#s`**) as inherent fields, and an **optional asserted total electron count** (**`#e<n>`**) as an inline constraint — as the **`:type`** value of a **`multicenter-bond-entry`** (**§4**). The **per-atom** electron contributions are the **mandatory leading `electron-counts`** of this string (below), not a separate map key.

```
multicenter-string ::= electron-counts multicenter-predicate*

multicenter-predicate ::= '#' tag payload
```

**Electron counts.** The string **MUST** begin with the **`electron-counts`** leaf (**§5.6**) — the per-atom electron contributions.

**Multicenter predicates.** **Zero or more** **`multicenter-predicate`** units following the **`electron-counts`**. **Optional** ASCII whitespace **MAY** appear before the first **`#`** and between successive predicates. **At most one** predicate per **tag** letter among **`c`**, **`u`**, **`s`**, **`e`**. **Canonical** predicate order (stable serialization): **`#c`**, **`#u`**, **`#s`**, **`#e`**.

**Whitespace** between **`#`** and the tag letter is **invalid** (**§7.1**).

**`#c` (multicenter-bond formal charge).** After **`#c`**, parse **either** a full **`value-expr`** (**§5.1**) **first**, **or** if that fails, a payload consisting **solely** of **`+`** (meaning **+1**) or **solely** of **`-`** (meaning **−1**), with **no** space between **`c`** and **`+`** / **`-`**. Same convention as atom (**§7.3**), bond (**§7.4**), and aromatic (**§7.8**) **`#c`**.

**`#u` / `#s` / `#e`.** After **`#u`**, **`#s`**, or **`#e`**, parse a **`value-expr`** (**§5.1**) **first**; if that fails, the **omitted** payload means numeric slot **1** (same convention as **§5.1.3** for decimal-only slots). **`#e`** omitted means **1** bonded electron.

| Tag | Meaning (multicenter-bond namespace) | Storage |
|-----|----------------------------------------|----------|
| **`#c`** | Multicenter-bond formal charge (**`i8`**, **§7.2**) | inherent field |
| **`#u`** | Unpaired electrons (bond-centered); **`u8`** | inherent field |
| **`#s`** | Spin multiplicity (2S+1) (bond-centered); **`u8`** | inherent field |
| **`#e`** | Asserted total bonded electron count; **`u8`** | inline constraint (`MulticenterBondConstraint::ElectronCount`) |

**`#e<n>` semantics.** **`#e<n>`** asserts the multicenter bond's total electron count and parses to an inline multicenter-bond constraint (`MulticenterBondConstraint::ElectronCount(n)`), parallel to the aromatic-system case (**§7.8**). Per-atom contributions in the string's leading **`electron-counts`** are the canonical data; **`#e<n>`** is the optional total assertion that downstream validation cross-checks against their **sum** on ground inputs.

**Per-atom participation.** The atom-side **`#m`** predicate (**§7.3**) is a per-atom multicenter-membership marker; the per-atom electron share for a given multicenter bond is the leading **`electron-counts`** of that bond's multicenter-string. Endpoint references (which atoms the bond spans) live in the **`:atoms`** vector of the **`multicenter-bond-entry`** (**§4**); they **MUST NOT** be encoded inside the multicenter-string.

### 7.10 Noncovalent-bond subgrammar

**Noncovalent-string** encodes the **interaction kind** of a single **`noncovalent-bond-entry`** (**§4**), optionally followed by the **`#I`** intramolecular predicate. The leading kind is the inherent field; **`#I`** is its one inline constraint.

```
noncovalent-string ::= noncovalent-kind-expr intramolecular?

intramolecular ::= '#I' ( '' | '+' | '!' | '*' )
                   (* '' / '+' intramolecular (true); '!' intermolecular (false); '*' undetermined *)
```

The **`noncovalent-kind-expr`** leading field is the **noncovalent-kind** leaf (**§5.7**). Leading / trailing whitespace on the whole **noncovalent-string** is ignored (**§7.1**).

**Intramolecular predicate (`#I`).** A trailing **`#I`** asserts whether the interaction is **intramolecular** — its two atoms lie in the **same** covalent connected component. The trailing polarity sets the truth value: **`#I`** / **`#I+`** intramolecular (true), **`#I!`** intermolecular (false), **`#I*`** undetermined. It is the noncovalent bond's **only** inline constraint (**`NoncovalentBondConstraint::Intramolecular`**); its structured EDN form is **`{:intramolecular bool}`** (**§7.12**). A **`#I*`** (undetermined) predicate is **vacuous** and **elided** from the canonical rendered string (**§7.1**), equivalent to omitting it.

### 7.11 Stereo subgrammar

**Stereo-string** uses a **separate** namespace from the atom / bond / aromatic / multicenter strings (**§7.2**). It is the **`:type`** payload of a **`stereo-atom-entry`** / **`stereo-bond-entry`** (**§4**) and names the coordination **`class`** plus the realized **`coset`** index over the entry's ordered **`:ligands`** frame. The **same** **`config`** grammar — **`coset`** preceded by two extra sentinels and **without** the leading **`class`** — is the payload of the atom **`#T`** / bond **`#C`** inline constraints (**§7.3** / **§7.4**) and the source form behind the **`{:stereo coset-form}`** EDN constraint (**§7.12**).

```
stereo-string ::= class coset stereo-predicate*

class ::= 'Th' | 'Ct' | 'Ax' | 'Sp' | 'Tb' | 'Oh'

stereo-predicate ::=
    '#p' ( '~' | ['\''] cycles ) boolean  (* ligand symmetry: ' improper, ~ kind involution; boolean trailing *)
  | '#f' ( '~' | cycles ) boolean         (* fluxionality (proper move); ~ kind involution; boolean trailing  *)
  | '#o' ligand-pair relation             (* topicity: ligand pair, then relation                             *)
  | '#g' relation                         (* stereogenicity classification                                    *)

cycles      ::= '()' | ( '(' nat (',' nat)* ')' )+   (* disjoint cycles, 0-indexed, identity ()       *)
ligand-pair ::= '(' nat ',' nat ')'                  (* two 0-indexed ligand-frame positions          *)
```

**Class.** **`Th`** tetrahedral, **`Ct`** cis/trans, **`Ax`** axial (allene-type), **`Sp`** square-planar, **`Tb`** trigonal-bipyramidal, **`Oh`** octahedral. A **`stereo-atom-entry`** carries an atom-centered class (**`Th`** / **`Ax`** / **`Sp`** / **`Tb`** / **`Oh`**); a **`stereo-bond-entry`** carries **`Ct`**. Matching presently realizes **`Th`** and **`Ct`**; **`Ax`** / **`Sp`** / **`Tb`** / **`Oh`** parse and round-trip but their matching is **staged**.

**Inline predicates (`#p` / `#f` / `#o` / `#g`).** After the leading **`class coset`**, a **`stereo-string`** carries **zero or more** stereo predicates — the inline form of the per-element stereo constraints (the molecule-scope structured peers are **§7.12**). Each predicate's permutation degree is the **`class`** degree (number of ligand positions). The four predicates are:

- **`#p`** — **ligand symmetry**: asserts whether a ligand permutation **is** (`boolean` **`+`** or omitted), **is not** (**`!`**), or is **undetermined** (**`*`**) a symmetry of the element. Payload **`(['`''] cycles | '~') boolean`**: **`'`** marks the permutation **improper** (orientation-reversing; default proper); **`cycles`** is the permutation in **disjoint-cycle notation** (below); the **trailing** **`boolean`** is the assertion's truth value — a scalar polarity written **after** the permutation (**`+`** / omitted true, **`!`** false, **`*`** undetermined). The **`~`** sugar denotes the **class involution** (the orientation-reversing generator for chiral classes, the configuration-swapping ligand permutation for achiral classes) and already carries the class-appropriate orientation, so it is not combined with **`'`**.
- **`#f`** — **fluxionality**: whether a proper ligand permutation is realized by dynamics. Payload **`(cycles | '~') boolean`**; the permutation carries no **`'`** (it is a bare proper move), the **trailing** **`boolean`** its truth value. **`~`** is the class involution.
- **`#o`** — **topicity** of a ligand pair. Payload **`ligand-pair relation`** — the **`(i,j)`** pair (two 0-indexed positions in the **`:ligands`** frame) **first**, then the **`relation`**. Its negation is **not** a trailing polarity like **`#p`** / **`#f`**: the relation's set-complement is a **leading** **`!`** on the glyph (**`!/`** = not diastereotopic; **Relation completeness** below). The glyphs are **`=`** homotopic, **`'`** enantiotopic, **`/`** diastereotopic.
- **`#g`** — **stereogenicity** classification. Payload **`relation`**: **`=`** symmetric, **`'`** prochiral, **`/`** stereogenic.

Each predicate places its **parameter** (the permutation for **`#p`** / **`#f`**, the ligand pair for **`#o`**, none for **`#g`**) directly after the tag, then its **value**. The value forms differ: **`#p`** / **`#f`** take a **trailing** scalar **`boolean`** — its negation marker **`!`** follows the permutation (**`#p(0,1)!`**); **`#o`** / **`#g`** take a **`relation`** whose set-complement is a **leading** **`!`** on the glyph (**`#o(0,1)!/`**). The two are different value types — a three-valued boolean vs. a set with complement — so the **`!`** trails on **`#p`** / **`#f`** but leads on **`#o`** / **`#g`**.

**Disjoint-cycle notation.** **`cycles`** is a product of disjoint cycles **`(p0,p1,…)(q0,q1,…)`**, **0-indexed** over the ligand frame, each cycle mapping **`p0→p1→…→p0`**; the identity is **`()`**. Cycle points **MUST** be in range (**`< class degree`**) and disjoint. This is the same permutation the **`Permutation`** `Display` emits and the structured **`permutation-form`** (**§7.12**) encodes as a vector of cycles.

**`~` rendering.** A **`#p`** / **`#f`** permutation equal to the class involution (and, for **`#p`**, matching the involution's orientation) renders as **`~`**; otherwise as explicit **`cycles`**.

**Chiral-class restriction.** The **`'`** value — **`#p`** improper, **`#o`** enantiotopic, **`#g`** prochiral — is meaningful only on a **chiral class** (**`Th`** atom; **`Ct`** / **`Sp`** are achiral). It **parses** on any class; an inconsistent class/value pairing is rejected by the **validator** (the resolved-symmetry cross-check, **§6.1**), not at parse.

**`stereo-keyword` shorthand (`§4`).** The four **`stereo-keyword`** values expand to canonical **`class`**+**`coset`** literals: **`:ccw`** → **`Th0`**, **`:cw`** → **`Th1`**, **`:z`** → **`Ct0`**, **`:e`** → **`Ct1`**. They are a ground EDN shorthand on the **`stereo-spec`**'s **`:type`** and are semantically identical to the expanded string. On serialization, implementations **MUST** emit the **`stereo-keyword`** for these four canonical shapes **only when the element carries no inline predicates**, falling back to the **`stereo-string`** otherwise.

### 7.12 Constraint grammar

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
    { :dative-bond-donors              [dative-bond-ref [atom-ref+]]           }
  | { :dative-bond-donor               [dative-bond-ref atom-ref]              }
  | { :dative-bond-contains-all-donors [dative-bond-ref [atom-ref+]]           }
  | { :dative-bond-all-donors          [dative-bond-ref atom-constraint-form]  }
  | { :dative-bond-any-donor           [dative-bond-ref atom-constraint-form]  }
  | { :dative-bond-acceptor            [dative-bond-ref atom-ref]              }
  | { :dative-bond-acceptor-satisfies  [dative-bond-ref atom-constraint-form]  }
  | { :dative-bond-parallels           [dative-bond-ref bond-ref]              }
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
    { :charge-sum                 { [:atoms [atom-ref+]]? :sum value-expr } }
  | { :unpaired-electron-coupling { [:atoms [atom-ref+]]?
                                     :unpaired-electrons unpaired-electrons-form } }
  | { :bond-order-sum             { [:bonds [bond-ref+]]? :sum value-expr } }
  | { :connected                  { [:atoms [atom-ref+]]? } }
  | { :sub-pattern                { :anchor anchor-spec :pattern molecule-map } }

combinator-constraint ::=
    { :and [constraint-entry+] }
  | { :or  [constraint-entry+] }
  | { :not constraint-entry }

atom-constraint-form ::=
    { :valence             value-expr }
  | { :donated-pairs       value-expr }
  | { :accepted-pairs      value-expr }
  | { :aromatic-valence    ( value-expr | :not-aromatic | :undetermined ) }
  | { :multicenter-valence ( value-expr | :not-multicenter | :undetermined ) }
  | { :tetrahedral-stereo  stereo-config-form }
  | { :degree              value-expr }
  | { :total-degree        value-expr }
  | { :total-valence       value-expr }
  | { :ring-degree         value-expr }
  | { :ring-valence        value-expr }
  | { :total-hydrogens     value-expr }
  | { :ring-membership     ring-membership-form }

bond-constraint-form ::=
    { :aromatic boolean }
  | { :ring-membership ring-membership-form }
  | { :cis-trans-stereo stereo-config-form }

dative-bond-constraint-form ::=
    { :aromatic boolean }
  | { :ring-membership ring-membership-form }

(* A boolean lattice value: the two literals plus the top (no constraint).  *)
boolean ::= true | false | :undetermined

(* One ring-membership fact. With :size, a count of rings of that size;     *)
(* without :size, the total ring count. :count is required.                 *)
ring-membership-form ::= { :size nat :count value-expr } | { :count value-expr }

unpaired-electrons-form ::= { :count value-expr :multiplicity value-expr }

aromatic-system-constraint-form  ::= { :electron-count value-expr }
multicenter-bond-constraint-form ::= { :electron-count value-expr }
noncovalent-bond-constraint-form ::= { :intramolecular bool }

(* A stereo entity-constraint form is a positional 2-vector: the element's *)
(* :kind (its stereo subtype) first, then a single-key predicate map. Kind is *)
(* first so the permutation degree is known before the predicate value is read *)
(* — a detached molecule-scope constraint cannot recover the degree otherwise *)
(* (the kind is many-to-one on degree). The position is fixed by the vector, *)
(* not by map-key order, so a streaming reader sees the kind before the value. *)
stereo-atom-constraint-form ::= [ stereo-kind stereo-predicate-map ]
stereo-bond-constraint-form ::= [ stereo-kind stereo-predicate-map ]

stereo-predicate-map ::=
    { :ligand-symmetry ligand-symmetry-form }
  | { :fluxionality    fluxionality-form }
  | { :topicity        topicity-form }
  | { :stereogenicity  stereogenicity-form }

stereo-kind          ::= :tetrahedral | :cis-trans | :axial | :square-planar
                       | :trigonal-bipyramidal | :octahedral
permutation-form     ::= [ cycle* ]                 (* vector of disjoint cycles; identity [] *)
cycle                ::= [ nat+ ]                    (* p0→p1→…→p0, 0-indexed positions *)
bool                 ::= true | false | :undetermined   (* the #p / #f trailing boolean, §7.11 *)
ligand-symmetry-form ::= { :permutation permutation-form [:orientation (:proper | :improper)]
                                                         [:invariant bool] }
fluxionality-form    ::= { :permutation permutation-form [:active bool] }
ligand-pair          ::= [ nat nat ]                (* two ligand-frame positions *)
topicity-form        ::= { :pair ligand-pair :relation topicity-relation }
stereogenicity-form  ::= { :relation stereogenicity-relation }
topicity-relation       ::= :homotopic | :enantiotopic | :diastereotopic | [ keyword+ ] | { :not-in [ keyword+ ] } | :undetermined
stereogenicity-relation ::= :symmetric | :prochiral | :stereogenic       | [ keyword+ ] | { :not-in [ keyword+ ] } | :undetermined

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
bond-ref             ::= int | keyword | { :atoms [atom-ref atom-ref] }
dative-bond-ref      ::= int | keyword | { :donors [atom-ref+] :acceptor atom-ref }
aromatic-system-ref  ::= int | keyword | { :atoms [atom-ref+] }
multicenter-bond-ref ::= int | keyword | { :atoms [atom-ref+] }
noncovalent-bond-ref ::= int | keyword | { :atoms [atom-ref atom-ref] }
stereo-atom-ref      ::= int | keyword | { :site atom-ref :ligands [ligand-ref+] }
stereo-bond-ref      ::= int | keyword | { :site bond-ref :ligands [ligand-ref+] }
```

**Ref resolution.** An integer ref is the **positional** index into the corresponding entity vector on the molecule map: **`atom-ref`** → **`:atoms`**, **`bond-ref`** → **`:bonds`**, **`dative-bond-ref`** → **`:dative-bonds`**, **`aromatic-system-ref`** → **`:aromatic-systems`**, **`multicenter-bond-ref`** → **`:multicenter-bonds`**, **`noncovalent-bond-ref`** → **`:noncovalent-bonds`**, **`stereo-atom-ref`** → **`:stereo-atoms`**, **`stereo-bond-ref`** → **`:stereo-bonds`**. A keyword ref resolves against the **`:id`** declared on the corresponding entry (**§4**). On serialization, implementations **MUST** emit the **`:id`** keyword when one is declared on the referenced entry, falling back to the positional integer otherwise.

**Structural refs.** The **map** form of a **non-atom** ref names the entity by its **participants** instead of by position or id: **`:atoms`** for a bond / noncovalent bond (a 2-vector) or an aromatic system / multicenter bond (the atom set); **`:donors`** + **`:acceptor`** for a dative bond; **`:site`** + **`:ligands`** for a stereo atom / stereo bond (the bearing site plus the ordered ligand frame — the same **`(site, ligand-multiset)`** key that identifies the element). It resolves by looking the participant key up among the entities of that kind; because each kind's participants are unique on a molecule (**§4.1**), at most one entity matches, and an **unmatched** key is a **parse error**. A structural ref map **MUST NOT** carry **`:type`** or **`:id`** (those mark an entity *definition*, not a ref). **`atom-ref`** has **no** structural form. Structural refs are **accepted wherever a ref is** — entity entries (a **`stereo-*-entry`** **`:site`**), entity / relational / molecule-scope constraints, sub-pattern anchor pairs, **`:bond-order-sum`** **`:bonds`**, reaction deltas, and reaction-span refs — and by **both** the tree and streaming parsers. They are **input-only**: the emission priority is **keyword > positional**, and a structural ref is **never** re-emitted (serialization produces only the **`:id`** keyword or the positional integer).

**Narrow inner forms for DAMN entities.** **`:aromatic-system`** and **`:multicenter-bond`** narrow leaves carry only the **`:electron-count`** value-only variant; every other predicate on those entities is a relational leaf instead. **`:noncovalent-bond`** narrow leaves carry only the **`:intramolecular`** value-only variant (**`#I`**, **§7.10**); every other noncovalent predicate is a relational leaf.

**Stereo entity constraints carry the kind.** The **`:stereo-atom`** / **`:stereo-bond`** entity-constraint forms (**`#p`** / **`#f`** / **`#o`** / **`#g`**) are a positional **2-vector** **`[stereo-kind stereo-predicate-map]`** — the element's stereo subtype first, then a single-key predicate map (so the leaf is **`{:stereo-atom [<ref> [<kind> {<predicate>}]]}`**). The kind is redundant with the referenced element at the **entity** level (the inline form omits it — the **`:type`** **`class`** supplies it, **§7.11**) but is **REQUIRED** at molecule scope, where the constraint is detached from its element: a permutation payload cannot recover its degree, and **`stereo-kind`** is many-to-one on degree (**`:tetrahedral`** and **`:square-planar`** are both degree 4). It is **first** (positional, container-fixed — not a map key) so the degree is known before the predicate value is read. The kind/degree (and the chiral-class restriction on **`'`** values) is cross-checked against the resolved element by the validator (**§6.1**); **`inline_constraints`** drops the carried kind back into the element. These constraints are **distinct** from the atom/bond **`:tetrahedral-stereo`** (**`#T`**) / **`:cis-trans-stereo`** (**`#C`**) inline configurations (which assert the local **coset** at the bearing atom/bond) and from the stereo **relational leaves** (**`:stereo-atom-…`** / **`:stereo-bond-…`**, no inline form).

**Anchor cardinality.** Each keyed slot in **`anchor-spec`** is optional and may appear at most once; if present, it is a vector of **`(target-side-ref, pattern-side-ref)`** pairs of the same entity kind. An empty **`anchor-spec`** denotes an unanchored sub-pattern (the pattern can embed anywhere). Target-side refs resolve against the outer molecule's metadata; pattern-side refs against the pattern molecule's metadata.

**Molecule-scope subset selectors.** `:charge-sum`, `:unpaired-electron-coupling`, `:bond-order-sum`, and `:connected` accept an **optional** `:atoms` (or `:bonds`) vector. When **omitted**, the predicate ranges over **every** atom (or bond) in the molecule, including atoms added by future structural growth. When **present**, the predicate ranges over the listed entities only. An empty vector `[]` is **distinct** from omission: it selects no entities.

**Sub-pattern materialization.** A **`:sub-pattern`** **`:pattern`** is a full **molecule-map**; its inner constraints are evaluated independently from the outer constraint tree. The pattern carries **no defaults** — values pass through verbatim — so a pattern's atom **`charge: undetermined`** stays **`undetermined`** at match time and behaves as a wildcard (**§5.1.4**); zero-defaulting that would apply to a ground input does **not** apply inside a pattern.

**Anonymous patterns.** A sub-pattern **`:pattern`** **MUST** be **anonymous**: an entry **MUST NOT** carry **`:id`** and the map **MUST NOT** declare **`:atom-aliases`**. Pattern entities are named positionally (or structurally) by the anchor and the pattern's inner constraints, so a symbolic namespace inside the pattern is meaningless; a pattern with either is a **parse error**.

**Sugar (inline string equivalents).** Narrow leaves whose entity has a string subgrammar admit two interchangeable serializations:

- the inline **`#tag`** payload on the entity's **`:type`** string (or, for atoms, on the atom literal directly);
- the **`{:<entity> [ref form]}`** entry in **`:constraints`**.

Parsers **MUST** accept both. Bare per-entity predicates (not nested under **`:and`** / **`:or`** / **`:not`** / **`:sub-pattern`**) **MAY** be emitted in the sugared inline form; nested predicates **MUST** be emitted as **`:constraints`** entries since the inline form has no logical context.

**Inline-form coverage by entity:**

- **Atom** (**§7.3**): all `atom-constraint-form` variants except the derived ones (`#D`, `#X`, `#V`, `#x`, `#y`, `#H`, `#R`) lift to inline atom predicates, including `:tetrahedral-stereo` → `#T`; the derived predicates also have inline tags but are pattern-only.
- **Bond** (**§7.4**): all `bond-constraint-form` variants (`:aromatic`, `:ring-membership`, `:cis-trans-stereo`) have inline forms (`#a`, `#R`, `#C`).
- **Dative bond** (**§7.7**): all `dative-bond-constraint-form` variants (`:aromatic`, `:ring-membership`) have inline forms (`#a`, `#R`).
- **Aromatic system** (**§7.8**): the single `aromatic-system-constraint-form` variant `:electron-count` has the inline form `#e<n>`.
- **Multicenter bond** (**§7.9**): the single `multicenter-bond-constraint-form` variant `:electron-count` has the inline form `#e<n>`.
- **Noncovalent bond** (**§7.10**): the single `noncovalent-bond-constraint-form` predicate (`:intramolecular`) has an inline form (`#I` on the `:type` string).
- **Stereo atom / stereo bond** (**§7.11**): all four `stereo-atom-constraint-form` / `stereo-bond-constraint-form` predicates (`:ligand-symmetry`, `:fluxionality`, `:topicity`, `:stereogenicity`) have inline forms (`#p`, `#f`, `#o`, `#g` on the `:type` string). On **inline** the kind is omitted (the `:type` `class` supplies it); on **lift** the element's kind is written as the **first element** of the molecule-scope form's 2-vector (**§7.12**). The atom/bond `:tetrahedral-stereo` / `:cis-trans-stereo` predicates (`#T` / `#C`) are separate atom/bond inline constraints; the `:stereo-atom-…` / `:stereo-bond-…` predicates are relational leaves with no inline form.

**Relational leaves** (**§7.12** `relational-constraint`) and **molecule-scope leaves** (`molecule-constraint`) have **no** inline form regardless of which entity they reference.

**Combining the inline and `:constraints` forms.** An entity **MAY** carry per-entity constraints in **both** serializations at once. They apply **conjunctively** — an entity's effective constraints are its inline predicates **together with** every molecule-scope per-entity entry that references it; neither serialization overrides the other. A same-kind clash with **conflicting** values (e.g. inline **`#v4`** and **`{:atom [i {:valence 3}]}`** on the same atom) is an unsatisfiable conjunction — a **contradiction** — and **MUST** be rejected as an error.

**Lift / inline.** The two storage scopes — inline on the entity (`AtomAst::constraints` etc.) and at molecule scope (`MoleculeAst::constraints` as `{:atom [ref form]}` peers) — are interchangeable for the inline-capable narrow leaves. Implementations **SHOULD** expose:

- **`lift_constraints`** (entity → molecule): drains every inline store into the molecule list as `{:<entity> [ref form]}` peers.
- **`inline_constraints`** (molecule → entity): drains top-level inline-capable narrow leaves from the molecule list into the targeted entity's inline store.

Combinator subtrees, relational leaves, and molecule-scope leaves are never moved by either operation. With multiple top-level entries targeting the same (entity, kind), `inline_constraints` resolves the collision via the entity store's per-kind insert policy (last-wins for unique-kind variants).

**Multiple constraints per entity.** Each per-entity constraint serializes as its **own** entity-constraint entry; implementations **MUST NOT** bundle multiple constraints on the same entity into a single map.


---

## 8. Edit and reaction documents

### 8.1 Standalone edit document

A standalone edit document is a bare ordered vector of host-specific edits. It carries no molecule
or molecule metadata: positional handles identify entities in the initial host, while **`{:new n}`**
identifies the **`n`**th entity of the corresponding kind created by an earlier edit in the same
vector. Creation ordinals are independent for all eight entity kinds.

```
edits-document ::= [ edit* ]

edit ::=
    { :atom                atom-edit }
  | { :bond                bond-edit }
  | { :dative-bond         dative-bond-edit }
  | { :dative-bonds        dative-bonds-edit }
  | { :aromatic-system     aromatic-system-edit }
  | { :aromatic-systems    aromatic-systems-edit }
  | { :multicenter-bond    multicenter-bond-edit }
  | { :multicenter-bonds   multicenter-bonds-edit }
  | { :noncovalent-bond    noncovalent-bond-edit }
  | { :noncovalent-bonds   noncovalent-bonds-edit }
  | { :stereo-atom         stereo-atom-edit }
  | { :stereo-atoms        stereo-atoms-edit }
  | { :stereo-bond         stereo-bond-edit }
  | { :stereo-bonds        stereo-bonds-edit }
  | { :topology            topology-edit }
  | { :constraint          edit-constraint }

atom-handle                ::= nat | { :new nat }
bond-handle                ::= nat | { :new nat }
dative-bond-handle         ::= nat | { :new nat }
aromatic-system-handle     ::= nat | { :new nat }
multicenter-bond-handle    ::= nat | { :new nat }
noncovalent-bond-handle    ::= nat | { :new nat }
stereo-atom-handle         ::= nat | { :new nat }
stereo-bond-handle         ::= nat | { :new nat }

atom-edit ::=
    { :add    atom-spec }
  | { :remove atom-handle }
  | { :modify [ atom-handle checked-atom-update ] }

bond-edit ::=
    { :add    [ atom-handle atom-handle bond-spec ] }
  | { :remove bond-handle }
  | { :modify [ bond-handle checked-bond-update ] }

checked-atom-update ::= { :expect partial-atom-string :update partial-atom-string }
checked-bond-update ::= { :expect partial-bond-string :update partial-bond-string }

dative-bond-edit ::=
    { :add    dative-bond-addition }
  | { :modify [ dative-bond-handle checked-dative-update ] }
dative-bonds-edit ::= { :remove [ dative-bond-removal* ] }

dative-bond-addition ::= { :donors [ atom-handle* ] :acceptor atom-handle
                            :type dative-bond-spec }
dative-bond-removal  ::= { :id dative-bond-handle :donors [ atom-handle* ]
                            :acceptor atom-handle :type dative-bond-spec }
checked-dative-update ::= { :expect partial-dative-string :update partial-dative-string }

aromatic-system-edit ::=
    { :add    aromatic-system-addition }
  | { :modify [ aromatic-system-handle checked-aromatic-update ] }
aromatic-systems-edit ::= { :remove [ aromatic-system-removal* ] }

aromatic-system-addition ::= { :atoms [ atom-handle* ] :type "aromatic-string" }
aromatic-system-removal  ::= { :id aromatic-system-handle :atoms [ atom-handle* ]
                               :type "aromatic-string" }
checked-aromatic-update ::= { :expect partial-aromatic-string
                              :update partial-aromatic-string }

multicenter-bond-edit ::=
    { :add    multicenter-bond-addition }
  | { :modify [ multicenter-bond-handle checked-multicenter-update ] }
multicenter-bonds-edit ::= { :remove [ multicenter-bond-removal* ] }

multicenter-bond-addition ::= { :atoms [ atom-handle* ] :type "multicenter-string" }
multicenter-bond-removal  ::= { :id multicenter-bond-handle :atoms [ atom-handle* ]
                                :type "multicenter-string" }
checked-multicenter-update ::= { :expect partial-multicenter-string
                                 :update partial-multicenter-string }

noncovalent-bond-edit ::=
    { :add    noncovalent-bond-addition }
  | { :modify [ noncovalent-bond-handle checked-noncovalent-update ] }
noncovalent-bonds-edit ::= { :remove [ noncovalent-bond-removal* ] }

noncovalent-bond-addition ::= { :atoms [ atom-handle atom-handle ]
                                 :type noncovalent-bond-spec }
noncovalent-bond-removal  ::= { :id noncovalent-bond-handle
                                 :atoms [ atom-handle atom-handle ]
                                 :type noncovalent-bond-spec }
checked-noncovalent-update ::= { :expect partial-noncovalent-string
                                 :update partial-noncovalent-string }

edit-ligand ::= atom-handle | [ :h atom-handle ] | [ :lp atom-handle ]

stereo-atom-edit ::=
    { :add    stereo-atom-addition }
  | { :modify [ stereo-atom-handle checked-stereo-update ] }
stereo-atoms-edit ::= { :remove [ stereo-atom-removal* ] }
stereo-atom-addition ::= { :site atom-handle :ligands [ edit-ligand* ] :type stereo-spec }
stereo-atom-removal  ::= { :id stereo-atom-handle :site atom-handle
                            :ligands [ edit-ligand* ] :type stereo-spec }

stereo-bond-edit ::=
    { :add    stereo-bond-addition }
  | { :modify [ stereo-bond-handle checked-stereo-update ] }
stereo-bonds-edit ::= { :remove [ stereo-bond-removal* ] }
stereo-bond-addition ::= { :site bond-handle :ligands [ edit-ligand* ] :type stereo-spec }
stereo-bond-removal  ::= { :id stereo-bond-handle :site bond-handle
                            :ligands [ edit-ligand* ] :type stereo-spec }
checked-stereo-update ::= { :expect partial-stereo-string :update partial-stereo-string }

topology-edit ::= { :remove { :atoms [ atom-handle* ] :bonds [ bond-handle* ] } }

edit-constraint ::=
    { :add    edit-constraint-entry }
  | { :remove edit-constraint-entry }
```

**Constraint handles.** An **`edit-constraint-entry`** has the complete **`constraint-entry`** shape
from **§7.12**, except that every reference into the target molecule is the corresponding typed
handle above. This includes references nested under **`:and`**, **`:or`**, and **`:not`**, relational
constraints, subset-valued molecule constraints, and the target side of a **`:sub-pattern`** anchor.
The nested pattern and the pattern side of each anchor pair retain the ordinary pattern-local refs of
**§7.12**. Keyword ids and structural refs are not accepted in standalone edits because the document
does not carry the host's metadata or structure. Normalized internal constraint slots never appear in
the surface form.

**Ordering and handle allocation.** Edit order and duplicate entries are semantic and **MUST** be
preserved. A **`{:new n}`** handle may refer only to an earlier same-kind creation; forward,
out-of-range, and removed handles are rejected during application. Parsing and collection rebuild
the next creation ordinal for every kind in one pass over the ordered entries. Atom and bond
additions are singleton surface entries; rendering a batched internal addition emits one entry per
created entity without changing its position relative to surrounding edits.

**Checked updates.** Every **`:modify`** carries both the expected partial value and the replacement
partial value. The two partials **MUST** address the same field or constraint keys. An undetermined
constraint denotes absence, so **`:expect "#v*"`** followed by **`:update "#v4"`** adds a valence
constraint, and the reverse removes it. Parsing reconstructs checked field and entity-constraint
edits; it does not validate chemical semantics against a host molecule.

**Removal preconditions and batching.** Atom-only and bond-only singleton removals use their
singular entity forms. **`:topology :remove`** is the inseparable operation for simultaneous atom and
bond removal, allowing one pre-removal snapshot and one compaction after cascading removal.
Overlay removals use the plural family key and retain a vector of complete removal records: the
entity handle, its participants, and its full AST value. The recorded values are checked during
application and retained for exact rollback. These semantic batches **MUST NOT** be lowered to a
sequence of independent removals.

**Defaults.** The same **`MoleculeDefaults`** value governs parsing and rendering. Defaults apply to
full definitions in additions and recorded overlay removals. They never apply to the partial
**`:expect`** or **`:update`** values. A semantic render/parse round trip therefore uses the same
defaults in both directions.

**Serialization.** Canonical serialization is the bare vector in stored order. It preserves
duplicates and semantic removal batches. Each **`Id(n)`** renders as **`n`** and each **`New(n)`** as
**`{:new n}`**, including handles nested inside constraints.

### 8.2 Reaction map

A reaction has two interchangeable surface forms, denoting the same transformation: the **operational** form (a left-hand side plus an edit list) and the **span** form (the superimposed `L ∪_K R` graph). The operational form is defined first; the span form follows.

**Operational form.** A **reaction map** describes a graph transformation as a left-hand-side molecule (**`:lhs`**) together with an **ordered** list of **`:deltas`** that edit it. The transformed (right-hand-side) molecule is the result of applying the deltas in order to **`:lhs`**; it is **not** written out. This is the **operational** surface — the lhs plus the edits — not a superimposed L∪K∪R graph.

```
reaction-map ::=
    { :lhs    molecule-map
      :deltas [ delta* ]
      [:atom-aliases atom-alias-list]? }

delta ::=
    { :atom             atom-delta }
  | { :bond             bond-delta }
  | { :dative-bond      dative-bond-delta }
  | { :aromatic-system  aromatic-system-delta }
  | { :multicenter-bond multicenter-bond-delta }
  | { :noncovalent-bond noncovalent-bond-delta }
  | { :stereo-atom      stereo-atom-delta }
  | { :stereo-bond      stereo-bond-delta }
  | { :constraint       constraint-delta }

atom-delta ::=
    { :add    atom-entry }
  | { :remove atom-ref }
  | { :modify [ atom-ref partial-atom-string ] }

bond-delta ::=
    { :add    bond-entry }
  | { :remove bond-ref }
  | { :modify [ bond-ref partial-bond-string ] }

(* The four DAMN overlays share one delta shape (§7.12 refs, §4 entries, §7.7–7.10 partials). *)
dative-bond-delta ::=
    { :add    dative-bond-entry }
  | { :remove dative-bond-ref }
  | { :modify [ dative-bond-ref partial-dative-string ] }

aromatic-system-delta ::=
    { :add    aromatic-system-entry }
  | { :remove aromatic-system-ref }
  | { :modify [ aromatic-system-ref partial-aromatic-string ] }

multicenter-bond-delta ::=
    { :add    multicenter-bond-entry }
  | { :remove multicenter-bond-ref }
  | { :modify [ multicenter-bond-ref partial-multicenter-string ] }

noncovalent-bond-delta ::=
    { :add    noncovalent-bond-entry }
  | { :remove noncovalent-bond-ref }
  | { :modify [ noncovalent-bond-ref partial-noncovalent-string ] }

(* Stereo adds three relative-op verbs; each carries an explicit stereo-kind (§7.12). *)
stereo-atom-delta ::=
    { :add    stereo-atom-entry }
  | { :remove stereo-atom-ref }
  | { :modify [ stereo-atom-ref partial-stereo-string ] }
  | { :swap   [ stereo-atom-ref stereo-kind ] }
  | { :mirror [ stereo-atom-ref stereo-kind ] }
  | { :apply  [ stereo-atom-ref stereo-kind "cycles" ] }

stereo-bond-delta ::=
    { :add    stereo-bond-entry }
  | { :remove stereo-bond-ref }
  | { :modify [ stereo-bond-ref partial-stereo-string ] }
  | { :swap   [ stereo-bond-ref stereo-kind ] }
  | { :mirror [ stereo-bond-ref stereo-kind ] }
  | { :apply  [ stereo-bond-ref stereo-kind "cycles" ] }

constraint-delta ::=
    { :add    constraint-entry }
  | { :remove constraint-entry }
```

**`:lhs`** is a **`molecule-map`** (**§4**); **`atom-entry`** / **`bond-entry`** are the **§4** entry forms (a bare spec, an **`[id spec]`** / **`[a b spec]`** vector, or the **`{:id … :atoms … :type …}`** map); **`constraint-entry`** is a single **§7.12** constraint (a per-entity narrow leaf, a relational leaf, a molecule-scope leaf, or a combinator). **`atom-ref`** / **`bond-ref`** are **§7.12** refs. The **overlay** deltas reuse the same pieces per family: their **entries** are the **§4** overlay entry maps, their **`:modify`** **partials** are the **§7.7–7.11** compact strings, and their **refs** (**`dative-bond-ref`** … **`stereo-bond-ref`**) are **§7.12** refs. **`:lhs`** and **`:deltas`** are **REQUIRED**; **`:atom-aliases`** is **OPTIONAL**.

**Reference id spaces.** A delta **`:remove`** / **`:modify`** target names an **existing lhs** entity, resolved in the **lhs id space** (positional index into **`:lhs`**'s **`:atoms`** / **`:bonds`**, or its declared **`:id`**). A created atom (**`:atom :add`**) **extends** the namespace; **bond endpoints** (in a **`:bond :add`**) and every **ref inside a `:constraint`** delta resolve against the **union** of lhs entities and reaction-created entities (lhs ∪ created). The same integer index that addresses an lhs entity in the lhs id space addresses a created entity once allocated: created atoms take indices continuing past the lhs atom count, in delta order.

**No forward references.** A delta **MAY** reference only entities present in **`:lhs`** or created by an **earlier** delta in **`:deltas`**. A reference to an entity created **later** is an error.

**Create vs. edit.** **`:remove`** and **`:modify`** **MUST** target an **lhs** entity; removing or modifying an entity **created in the same reaction** is an error (collapse the creation into its final state instead). **`:add`** introduces a new entity.

**`:modify` payload.** The **`partial-atom-string`** (**`partial-bond-string`**) is a compact **atom-string** (**bond-string**, **§7.3** / **§7.4**) carrying **only** the changes: a field left **`undetermined`** (e.g. an omitted element) keeps the lhs value; a field with a definite value **overwrites** it; a constraint predicate **sets** that constraint; an **undetermined** predicate written as **`#tag*`** **removes** it (**§7.1** — the same vacuous form that is elided on a full render is, on a **`:modify`** partial, the explicit **removal marker**). Consecutive **`:modify`** edits to the **same** entity (of any family) **coalesce** on serialization into a **single** **`:modify`** with one merged partial. The overlay partials work the same way over their own strings (**§7.7–7.11**).

**Overlay deltas.** The six overlay families — **`:dative-bond`**, **`:aromatic-system`**, **`:multicenter-bond`**, **`:noncovalent-bond`**, **`:stereo-atom`**, **`:stereo-bond`** — take **singular** delta keys, matching **`:atom`** / **`:bond`** (the **plural** **`:dative-bonds`** … keys name the **`:lhs`** molecule-map **collections**, **§4**, not deltas). Each shares the atom/bond delta shape: **`:add`** an entry, **`:remove`** a ref, **`:modify`** an **`[ref partial]`** pair. A **`:remove`** / **`:modify`** target resolves in the **lhs id space** of that family; **`:add`** allocates the next id of the family and (like a created atom) participants resolve against the lhs ∪ created union.

**Stereo relative ops.** **`:stereo-atom`** / **`:stereo-bond`** add three verbs that transform the coset in place: **`:swap`** (the class involution), **`:mirror`** (the enantiomer), and **`:apply`** (a ligand-frame permutation in disjoint-cycle notation, **§7.11**). Each carries an **explicit `stereo-kind`** — **`[ref stereo-kind]`**, or **`[ref stereo-kind "cycles"]`** for **`:apply`**. The kind is **REQUIRED** because the coset algebra is parametrized by it (a relative op is uninterpretable without it) and carrying it makes the delta **self-contained** — independent of the lhs entity, so it is well-formed even when the lhs coset is open. The **`"cycles"`** permutation's degree is the **`stereo-kind`** degree (**§7.11**).

**Stereo `:modify` partial.** The **`partial-stereo-string`** is the modify-variant of the stereo-string (**§7.11**): the **`coset`** is **optional** (omitted = unchanged — it keeps the lhs coset), but the **`class`** **MUST** be present once a coset or predicate appears, since the predicates render and parse against it. So **`"*"`** alone (undetermined, no predicates), or **`"Th"`** / **`"Th1"`** / **`"Th#o(0,1)="`** — but **`"*#o…"`** (a predicate with no class) is a parse error.

**`:constraint` deltas.** **`{:constraint {:add …}}`** / **`{:constraint {:remove …}}`** add or remove one molecule-scope or per-entity constraint (**§7.12**); refs inside the **`constraint-entry`** resolve against the lhs ∪ created union.

**`:atom-aliases`.** As in **§4**, with the alias namespace spanning lhs ∪ reaction. Aliases are resolved **after** the entire map is read, independent of tree vs. streaming parse, so their position in the top-level map is **not** significant; canonical serialization emits **`:atom-aliases`** **last**.

**Serialization.** A reaction map serializes its keys in the order **`:lhs`**, **`:deltas`**, then **`:atom-aliases`** (only when aliases are present). **`:lhs`** renders per **§4**; deltas render in **stored order** (the canonical AST order, not source order); each ref renders as its **`:id`** keyword when one is declared on the referenced entry, falling back to the positional integer (**§7.12**). Serializing a reaction that carries **no** surface metadata emits the **positional** form throughout (no **`:id`** keywords, no aliases); **`:id`** / alias output requires retaining the declared ids and aliases alongside the structural graph.

**Span form.** A reaction may instead be written as its **superimposed span** — the `L ∪_K R` graph overlaying the before and after states — which shares the **molecule-map shape** (**`:atoms`** / **`:bonds`** / **`:constraints`** / **`:atom-aliases`**). Each entry is either a **bare** molecule entry (**unchanged** — present and identical on both sides) or that entry wrapped in a single-key **verb** map, **`:add`** / **`:modify`** / **`:remove`** (the **same** verbs as the operational deltas).

```
span-map ::=
    { :atoms             [ atom-span* ]
      [:bonds            [ bond-span* ]]?
      [:dative-bonds     [ dative-bond-span* ]]?
      [:aromatic-systems [ aromatic-system-span* ]]?
      [:multicenter-bonds [ multicenter-bond-span* ]]?
      [:noncovalent-bonds [ noncovalent-bond-span* ]]?
      [:stereo-atoms     [ stereo-atom-span* ]]?
      [:stereo-bonds     [ stereo-bond-span* ]]?
      [:constraints      [ constraint-span* ]]?
      [:atom-aliases     atom-alias-list]? }

atom-span      ::= atom-span-body | [ keyword-id atom-span-body ]
atom-span-body ::=
    atom-value                              (* Unchanged *)
  | { :add    atom-value }
  | { :remove atom-value }
  | { :modify [ atom-value atom-value ] }   (* [left right] *)

bond-span ::=
    bond-entry                              (* Unchanged *)
  | { :add    bond-entry }
  | { :remove bond-entry }
  | { :modify ( [ atom-ref atom-ref [ bond-value bond-value ] ]
              | { [:id keyword]? :atoms [ atom-ref atom-ref ] :type [ bond-value bond-value ] } ) }

(* The six overlay spans mirror bond-span: a bare §4 entry (Unchanged), or an :add / :remove /   *)
(* :modify wrapper. :modify restates the entry map with a two-element [left right] :type pair    *)
(* (participants once). <x>-value is the family's :type string (§7.7–7.11).                     *)
dative-bond-span      ::= dative-bond-entry      | { :add dative-bond-entry }      | { :remove dative-bond-entry }      | { :modify dative-bond-modify }
aromatic-system-span  ::= aromatic-system-entry  | { :add aromatic-system-entry }  | { :remove aromatic-system-entry }  | { :modify aromatic-system-modify }
multicenter-bond-span ::= multicenter-bond-entry | { :add multicenter-bond-entry } | { :remove multicenter-bond-entry } | { :modify multicenter-bond-modify }
noncovalent-bond-span ::= noncovalent-bond-entry | { :add noncovalent-bond-entry } | { :remove noncovalent-bond-entry } | { :modify noncovalent-bond-modify }
stereo-atom-span      ::= stereo-atom-entry      | { :add stereo-atom-entry }      | { :remove stereo-atom-entry }      | { :modify stereo-atom-modify }
stereo-bond-span      ::= stereo-bond-entry      | { :add stereo-bond-entry }      | { :remove stereo-bond-entry }      | { :modify stereo-bond-modify }

dative-bond-modify      ::= { [:id keyword]? :donors [ atom-ref+ ] :acceptor atom-ref :type [ dative-value dative-value ] }
aromatic-system-modify  ::= { [:id keyword]? :atoms  [ atom-ref+ ]                   :type [ aromatic-value aromatic-value ] }
multicenter-bond-modify ::= { [:id keyword]? :atoms  [ atom-ref+ ]                   :type [ multicenter-value multicenter-value ] }
noncovalent-bond-modify ::= { [:id keyword]? :atoms  [ atom-ref atom-ref ]           :type [ noncovalent-value noncovalent-value ] }
stereo-atom-modify      ::= { [:id keyword]? :site atom-ref :ligands [ ligand+ ]     :type [ stereo-value stereo-value ] }
stereo-bond-modify      ::= { [:id keyword]? :site bond-ref :ligands [ ligand+ ]     :type [ stereo-value stereo-value ] }

constraint-span ::=
    constraint-entry                        (* Unchanged *)
  | { :add    constraint-entry }
  | { :remove constraint-entry }            (* no :modify — constraints are a by-value multiset *)
```

**`atom-value`** is an atom literal or an **`:atom-aliases`** keyword (the value position of a **§4** atom entry); **`bond-value`** is a **`bond-string`** / **`bond-keyword`** (a bond entry's **`:type`**); each **`<x>-value`** is the corresponding overlay's **`:type`** payload (**`dative-value`** a **`dative-string`**, **`aromatic-value`** an **`aromatic-string`**, and so on, **§7.7–7.11**); **`atom-entry`** / **`bond-entry`** / the overlay entries / **`constraint-entry`** are the **§4** / **§7.12** forms; **`ligand`** is a **§4** stereo ligand.

**Union id space.** Span entries are in the **union id space** (`L ∪ R`): every entity — unchanged, added, removed, or modified — occupies a slot, and positions are the ids (no allocation, unlike the operational form). An atom's optional **`[:id …]`** and a bond's **`:id`** name it; refs (bond endpoints, constraint refs) resolve against this single id space.

**`:modify` carries both sides** — `[left right]`, **complete** values (atoms `[left right]`; bonds carry endpoints once, `[a b [left right]]`). This is the one place the span differs from the operational **`:modify`**, which carries only the new value: the span is **self-contained** (it has no **`:lhs`** to recover the old value from).

**Constraints** are a by-value multiset — bare = unchanged, **`{:add c}`** / **`{:remove c}`**; there is **no** **`:modify`**.

**Per-side consistency.** The **left** projection (unchanged ∪ removed ∪ modified-left) and the **right** projection (unchanged ∪ added ∪ modified-right) **MUST** each be a valid molecule — every ref a side uses **MUST** resolve within that side (e.g. a bond present on the left **MUST** have both endpoints present on the left; an overlay present on a side **MUST** have all its participants — donors/acceptor, atoms, endpoints, or stereo site + ligands — present on that side; a stereo bond's site bond **MUST** be present too).

**Homoiconicity.** A plain molecule map (every entry bare) is a valid span — the identity reaction; the molecule *is* the span.

**Equivalence.** The span and operational forms denote the same reaction: the span's **left** projection is the operational **`:lhs`**, and applying the **`:deltas`** yields the **right** projection. Both reuse the whole molecule entity / constraint grammar unchanged; the span adds only the **`{:add|:modify|:remove}`** wrapper.

A substitution (replace an O by an N on a carbon) as a span:

```clojure
{:atoms ["C"                     ; Unchanged
         {:remove "O"}           ; Removed — left only
         {:add "N"}]             ; Added — right only
 :bonds [{:remove [0 1 :single]} ; Removed — C–O
         {:add [0 2 :single]}]}  ; Added — C–N
```

A bond-order change, with both atoms unchanged:

```clojure
{:atoms ["C" "C"]
 :bonds [{:modify [0 1 [:single :double]]}]}  ; order 1 → 2, endpoints once
```

## 9. Molecule map examples (non-normative)

Examples use the vector **`:atoms`** form with inline ids. Bond entries show **`:id`** where useful.

### 9.1 Methanol (CH₃OH) — Ground

```clojure
{:atoms [[:C  "C#h3"]
         [:O  "O#h1"]
         [:H  "H"]]
 :bonds [[:C :O :single]
         [:O :H :single]]}
```

The **`H`** atom here represents an **explicit** hydrogen (e.g. a hydroxyl H one wishes to name). Implicit H counts on **`C`** (**`#h3`**) and **`O`** (**`#h1`**) already account for the remaining hydrogens.

### 9.2 Indole — Ground, aromatic ring

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
 :aromatic-systems [{:id :ar1 :atoms [:N :C2 :C3 :C3a :C7a] :type "*"}
            {:id :ar2 :atoms [:C3a :C4 :C5 :C6 :C7 :C7a] :type "*"}]}
```

Localized bonds carry the σ-skeleton orders; the aromatic π system is expressed in **`:aromatic-systems`**.

### 9.3 Substructure query

Match any carbon with at least two implicit hydrogens that is directly bonded to a nitrogen:

```clojure
{:atoms [[:C "C#h(?h >= 2)"]
         [:N "N"]]
 :bonds [[:C :N :single]]}
```

**`(?h >= 2)`** is a **`bool-expr`** payload on **`#h`**; **`?h`** is bound to the matched atom's implicit H count.

### 9.4 Reaction (transformation rule)

A reaction is **`:lhs`** plus an ordered **`:deltas`** edit list (**§8**); the right-hand side is the result of applying the deltas, not written out.

Strip the three implicit H from a primary-amine carbon (C bonded to NH₂), leaving a quaternary carbon — a single field edit on the lhs **`:C`**:

```clojure
{:lhs {:atoms [[:C "C#h3"]
               [:N "N#h2"]]
       :bonds [[:C :N :single]]}
 :deltas [{:atom {:modify [:C "#h0"]}}]}
```

A reaction that grows the graph and asserts a constraint — add a hydroxyl O to the carbon, then require the result be connected:

```clojure
{:lhs {:atoms [[:C "C#h3"]]}
 :deltas [{:atom {:add [:O "O#h1"]}}
          {:bond {:add [:C :O :single]}}
          {:constraint {:add {:connected {}}}}]}
```

The added **`:O`** atom is visible to the later **`:bond :add`** (no forward references, **§8**); the **`:connected`** molecule constraint (**§7.12**) ranges over the post-edit graph.
