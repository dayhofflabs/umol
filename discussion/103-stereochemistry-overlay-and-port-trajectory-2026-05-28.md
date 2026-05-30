# Stereochemistry: overlay deliverable and port-model trajectory — 2026-05-28

Status: **Active / decision record.** Synthesis of the design conversation that ran over
[101-stereochemistry-framework-2026-05-28.md](101-stereochemistry-framework-2026-05-28.md)
(framework + reference review) and
[102-stereochemistry-framework-research-2026-05-28.md](102-stereochemistry-framework-research-2026-05-28.md)
(independent research + maximal proposal). This doc records the chosen path, the overlay
model to build now, the maximal port model to keep in view, the contract steps to take
now, the forward-compat requirements, and the `umol-port` spike.

## Decision summary

- **Ship a stereochemistry overlay on the existing atom/bond graph now.** Stereo lives at
  the **molecule level** as a separate element collection that *references* structural
  sites.
- **Spin up an experimental `umol-port` crate off the critical path**, timeboxed, to
  validate the maximal port/entity model. **Migrate to ports only if the spike succeeds.**
- The overlay and the eventual port model share **one contract** (DSL surface + public
  semantic type). Three decisions (D1–D3, below) make that contract invariant across the
  migration, so the overlay is the first backend of a durable contract, not a throwaway.
- Rationale for sequencing: the overlay has bounded, known scope (no graph-core change, no
  parser change, additive only); the port model carries genuine research risk (no
  reference implementation commits to a port-factored electron model as storage). The
  deadline path must hold the bounded item; the unbounded item stays a spike.

## Scope of the current deliverable

In:

- **Organic central chirality** (tetrahedral) and **stereobonds** (cis/trans, i.e. E/Z
  geometry). Standard organic stereo only.
- **Parse and preserve** SMILES `@`/`@@`, `/`/`\`, and the stereo tokens of SMIRKS
  reaction patterns; round-trip through the DSL `:stereo` key (class in the config head).
- **Assert / match** stereo properties in patterns (query side), including **relative
  stereo** via binds and the involution `~` ("the other config" — TH meso/d,l, CT E/Z relation).
- **Convert to 3D**: a stereo element drives the local geometry during embedding
  (umol-geometric).
- The operator notation (`^` generic apply, `~` = `^1` involution) is **defined in the
  contract**; only `^1` / `~` for the binary classes (TH/CT) is implemented. A distinct mirror
  operator is future-additive (no breaking change).

Out (deferred; reachable additively — see forward-looking section):

- **CIP** R/S labels — explicitly **not needed** for this deliverable.
- Non-tetrahedral *perception* (SP / TB / OH / allene): parse and preserve the arrangement
  index opaquely for round-trip, do not perceive.
- Ports, the entity/electron model, lone-pair and haptic stereo.
- Full reaction **stereo-transfer** semantics (retain / invert / erase on rule
  application). Parsing and matching stereo in SMIRKS is in scope; generating products
  with transferred stereo is forward-looking (D3 keeps it additive).

## Decisions and rationale (conversation record)

- **Molecule-level, not per-atom/per-bond.** Per-atom chirality / per-bond E/Z is
  context-dependent under renumbering, substructure extraction, and reactions, and
  excludes axial/planar/helical chirality. Every mature general-chemistry toolkit except
  RDKit already stores stereo as a molecule-level element collection (CDK, Indigo,
  LillyMol, OpenBabel, StereoMolGraph); CDK ships it as the primary store.
- **A stereo element is `(focus, class, ordered carriers, configuration)`** with carriers
  **always explicit** — including virtual non-bond carriers (implicit-H / lone-pair) as
  distinct sentinel entries. (Superseded: an earlier draft derived carriers from
  focus+topology+counts; that is withdrawn — see "Superseded directions.")
- **Class is a structural discriminator, not an attribute — no algebra over it.** The class
  fixes *what the focus is* and *how many carriers of what kind*, so a class set/wildcard is
  ill-typed (a set of orbit groups is not a config space) and `"*…"` over class is meaningless.
  Carried as a **fixed enum in the config head** (`"Th1"`, `"Ct2"`), one flat `:stereo` key —
  no algebra (a value-expr in the class position is a parse error), and it **parallels the
  atom/bond constraint** `#T1`/`#C1` (class with config at both levels). Chosen over per-class
  top-level keys because it makes the config head self-describing — `:cw`/`:e` sugars carry
  their class — at the cost of a softer analogy to `:bonds`/`:dative`.
- **Configuration is definitional → it is the head of the per-entry definition, not a
  `#`-predicate.** The head is the inherent definition (without it the entity is empty);
  `#`-predicates are assertions layered on top (`#e6` asserts about an aromatic system, it
  does not define it). Config is definitional, so it is the head. (Superseded: the `#o` /
  `#g` two-tag scheme is withdrawn.)
- **Storage is equivariant, not canonical (refines D2).** The configuration is the `Sₙ/R`
  coset relative to the **explicit `:carriers` list order**. Equivariant: a structural edit is
  a clean **relabel of the carrier refs with the config untouched** (config is over list
  positions); attribute mutation never moves carriers. The atom/bond-numbering trade extended
  to stereo — not canonical in storage; canonicalization is external. Canonical labeling / CIP
  / the `A_s` quotient are **derived on-demand** views (equality-up-to-isomorphism,
  serialization, meso). The user-facing "index-independent notation" is this derived canonical
  view; storage is equivariant.
- **Stereo is irreducibly a rotation datum; the atom/bond constraint is a scalar parity, not
  the carrier-ordered config.** Stereo = a combinatorial-map / rotation-system datum (ordered
  incident edges; the order is irreducible — see "Adjacent-field precedents"). The full
  carrier-ordered configuration lives only in the molecule-level element. What projects to the
  atom/bond DSL **without a variadic entity-ref list** is a **scalar parity**; a dichotomy
  (canonical index-independent-but-global vs incidence-order equivariant-but-local) is **resolved
  to the incidence frame** — canonical rank is unstable under topology edits *and* attribute
  mutation, breaking the "edits are remaps" invariant. The constraint is an **uppercase per-class
  tag** (derived-predicate namespace, like `#R`/`#D`): `#T` tetrahedral, `#C` cis-trans (reserved
  `#A`/`#Q`/`#B`/`#O`), a **derived projection** of the element onto the atom's local incidence
  order — distinct notation from the element's `:type` head, so the frame is unmistakable.
- **CIP is a boundary codec, not the internal basis** — even with unlimited time.
  Graph-canonical is the better internal basis: it is the same machinery canonicalization
  needs anyway, is uniform across all stereo classes, and avoids CIP Rule 5's
  self-reference (priority depending on remote descriptors). Para-stereocenters are handled
  by a refinement **fixpoint** (fold parity into the color, re-refine), the same fixpoint
  InChI uses. CIP, when wanted, is an I/O dialect on both read and write.
- **Topicity and prochirality fall out of the automorphism action.** Homotopic /
  enantiotopic / diastereotopic = orbit relations under `Â` (proper) vs `Â*` (with
  parity-inversion); a site is prochiral iff it has an enantiotopic carrier pair. The
  existing `AtomAutomorphism` (orbits, canonical labeling, group order) is the constitutional
  `Â`; stereo adds the parity-inversion extension.
- **Ports endorsed.** Undirected half-edges suffice for storage: donation direction is read
  from per-port electron contribution, so no directed darts are needed at the storage layer
  (direction matters only in the stereo configuration's view, downstream). The general model
  is the Ugi–Dugundji bond-electron matrix **factored through ports** (an entity–port
  incidence matrix with port-ownership) — equivalently a localized-MO coefficient matrix
  over a minimal valence basis, block-structured by owner. The appeal is uniformity; the
  cost is novelty.
- **Ship atoms/bonds first; migrate only if the spike succeeds.**

## The overlay model (build now)

### Semantic type

A molecule-level `StereoElement` collection on `MoleculeAst`. The storage container (a
dedicated overlay vs. one of the existing relation-set abstractions) is an internal choice,
leaning a dedicated overlay; it is insulated from the contract by D1/D2 and confirmed at
implementation kickoff. Conceptually each element carries:

- **class** — `Tetrahedral` / `CisTrans` (deliverable); a **structural discriminator** fixing
  the focus type, carrier count, and orbit group `R`. No algebra over it. Carried as a **fixed
  enum in the config head** (`"Th…"`, `"Ct…"`), one flat `:stereo` key.
- **focus** — a ref whose type is fixed by the class: an **atom** (tetrahedral) or a **bond**
  (stereobond). Polymorphic across classes (CDK's `IStereoElement<Focus,Carrier>`), determined
  per entry by the head class.
- **ordered carriers** — an **explicit** list, length fixed by the class. Entries are
  **neighbor atom-refs** for real substituents and **virtual sentinels** `[:implicit-hydrogen ref]` /
  `[:lone-pair ref]` for non-bond carriers (implicit-H / lone-pair), each a distinct entry. Always
  listed; real neighbors in neighbor order, non-bond carriers trailing (H then LP).
- **configuration** — the orbit in `Sₙ/R` relative to the **explicit carrier-list order**.
  Three-valued plus achiral, per StereoMolGraph: `None` (undefined), `0` (defined, achiral),
  `±1` (chiral handedness). This is the **definition**, written in the per-entry head (not a
  `#`-predicate).
- **specified-ness** — `Specified` / `Undefined` / `Unknown`, carried by the config head's
  `*`/`+` channels (distinct from achiral).
- **group** — enhanced stereo I/O sugar: `Abs` / `And(n)` / `Or(n)`.

**Carriers are always explicit, including non-bond sentinels.** Required, not optional —
deriving them consults the graph (elision), does not generalize to allene/axial/planar
(whose carriers are not the focus's incident bonds), and an optional-here-required-there rule
is incoherent. Non-bond carriers are listed as distinct sentinel entries (not implied by
count), which preserves prostereogenic expressivity — e.g. 1,1-dichloroethene's `=CH₂` needs
its two H as **distinct** carrier slots, or substitution to different isomers can't be
expressed — and makes structural edits a clean ref-relabel (the murky "coset update, not a
remapping" only arose from the implicit `[bond, bond]` count scheme).

**What half-edges / ports buy for stereo specifically.** Stereo does not order "atoms" in the
abstract; it orders **attachment sites around a focus**. In the deliverable, a neighbor atom is
usually a readable name for that site. The half-edge view becomes useful only when that shorthand
is no longer one-to-one:

- if the same atom participates through more than one local attachment site, atom refs collapse
  sites that stereo must distinguish;
- if the carrier is non-atomic (implicit H, lone pair, vacancy, haptic ligand group), an atom ref
  is at best a sentinel rather than the site itself;
- if a rewrite preserves the neighboring atom but changes which attachment site is retained,
  a half-edge/port names the stereochemical slot directly.

This is the stereo-specific reason to keep `:bond` / `:port` / `:fragment` carrier arms reserved.
It is not an argument that the current molecule model must become port-based before the stereo
overlay ships.

No `chirality` field is added to `AtomAst`; no `stereo` field is added to `BondAst`. This
firewall is what keeps the port migration open. Stereo lives **only** at the molecule level;
atom/bond-level stereo is a scalar-parity constraint / derived predicate (see "Ordering, the
atom/bond constraint, and perception"), never a stored carrier-ordered config on the atom.

### DSL: one `:stereo` key, class in the config head

The head-vs-predicate distinction sharpens the structural/scalar principle:

> **Structural references → map keys. The inherent definition → the head of the type-string.
> Assertions → `#`-predicates on top.** (Head = the minimal definition without which the entity
> is empty: `C` for an atom, `1` for a bond, `Hbd` for a noncovalent bond, `""` for an aromatic
> system; `#e6` is a predicate because it *asserts about* an aromatic system, not defines it.)
> For stereo: `:focus` / `:carriers` are structural refs → map keys; the **class + configuration
> together are the inherent definition** → the **head** of the per-entry config string.

The **class is a fixed enum in the config head** (`Th`, `Ct`, …), *not* a value-expr — so there
is no algebra over class (`*Th1` / `{Th,Sp}1` are parse errors), which is what keeps the
focus/carrier interpretation determined (the parser reads the head class). This **parallels the
atom/bond constraint** `#T1` / `#C1` (class in the tag, config in the payload — see the
constraint section): class travels with config at both levels, distinguished by frame and
position. One flat `:stereo` key holds class-tagged entries — chosen over per-class top-level
keys because it parallels the constraint and makes the config head self-describing (so the
`:cw`/`:e` sugars carry their class, and `"Th1"` is interpretable on its own); the cost is a
softer analogy to `:bonds`/`:dative`.

```
:stereo       ::= [ entry* ]
entry         ::= { [:id keyword]
                    :focus    ref                   ; type fixed by the head class (atom | bond)
                    :carriers [ carrier+ ]          ; EXPLICIT, length fixed by class (order below)
                    :type     config-string         ; class HEAD + config — or keyword sugar
                    [:group   stereo-group] }

carrier       ::= atom-ref                          ; sugar: real substituent by neighbor atom
                | carrier-vector                    ; typed non-atom / future extension carrier

carrier-vector ::= [ :implicit-hydrogen atom-ref ]   ; virtual: one implicit H on that atom
                 | [ :lone-pair atom-ref ]           ; virtual: one lone pair on that atom (deferred)
                 | [ :bond bond-ref atom-ref ]       ; reserved: half-edge of a bond at this focus
                 | [ :port port-ref ]                ; reserved: first-class port carrier
                 | [ :fragment fragment-ref ]        ; reserved: fragment / haptic group carrier
                 | [ :opaque keyword edn-value* ]    ; reserved extension form, preserved not interpreted
                ; Duplicates distinct by position: two [:implicit hydrogen :C2] = two distinct H slots.
stereo-group  ::= :abs | {:and int} | {:or int}       ; ABS/AND/OR I/O sugar

config-string ::= class config stereo-predicate*     ; class head + config + future #-assertions
class         ::= 'Th' | 'Ct'                         ; fixed enum (NO algebra); reserved Al Sp Tb Oh
config        ::= '*'                                 ; Undetermined — no constraint (pattern)
                | '+'                                 ; stereogenic, configuration unspecified (wavy)
                | config-value
config-value  ::= coset-term
                | '~' config-value                    ; "the other" involution (= ^1); involutive
                | config-value '^' nat                 ; generic apply (GAP x^g); only ^1 now
coset-term    ::= nat | nat-set | '?' id | '?' id '::' nat-set
                ; nat = OpenSMILES-style arrangement index relative to the :carriers order
```

**Surface carrier contract.** The compact form is still the ordinary one: a real substituent is
written as its neighbor atom ref. Typed vector forms exist for carriers that are not a neighbor
atom, or for a future model where the stereo site must be named more precisely. This keeps common
organic stereo readable while avoiding an untyped extension hole.

Examples:

```clojure
; ordinary atom carriers
:carriers [:F :Cl :Br [:implicit-hydrogen :C]]

; prostereogenic duplicate implicit-H slots: distinct by vector position
:carriers [:Cl1 :Cl2 [:implicit-hydrogen :C2] [:implicit-hydrogen :C2]]

; later: a lone pair as a stereochemical carrier
:carriers [:O :Me :Ph [:lone-pair :S]]

; later: an explicitly named half-edge when atom-ref sugar is not precise enough
:carriers [[:bond :b1 :C] [:bond :b2 :C] [:implicit-hydrogen :C] [:lone-pair :C]]

; later: direct port carrier after a port model exists
:carriers [[:port :p17] [:port :p18] [:implicit-hydrogen :C] [:lone-pair :C]]

; later: haptic/fragment carrier, e.g. an eta-bound ligand group as one site
:carriers [[:fragment :cp-ring] :CO :PPh3 :Cl]
```

Unknown carrier vectors are **not** silently interpreted. A strict implementation rejects
unknown vector tags; an extension-preserving implementation may round-trip only the explicit
`:opaque` form. This follows the rest of the DSL: open maps can grow, but recognized grammar
arms have fixed meanings.

**Why atom refs first.** For the deliverable's atom/bond graph, a real carrier is the
**neighbor atom**, not the connecting bond. The bonds-as-carriers idea (an earlier draft) is
withdrawn as the default: its main motivation was multi-edge disambiguation, which does not arise
in umol's `:bonds` graph (it uses **order-n edges, not parallel edges** — a double bond is one
edge of order 2; dative/multicenter interactions live in separate relation sets), so a neighbor
atom names the direction unambiguously in the common TH/CT cases. Atoms also win on readability,
on not forcing optional bond ids, and on aligning with the neighbor-order convention. The typed
`:bond` arm is reserved for the uncommon case where the stereochemical **site** is an incidence
rather than just the neighboring atom. (Ports later: an atom-carrier resolves to "the focus's
port toward that atom" when that port is unique, so atom sugar can lower to port refs without
changing the written contract.)

**Carrier order.** Real carriers in **neighbor order** — the atom's neighbor ordering
(atom-focus / `#T`), or atom-1's neighbors then atom-2's (bond-focus / `#C`) — then **implicit-H
slots (× count), then lone-pair slots (× count)** trailing, since the non-bond count is variable.
This is *neighbor* order (keyed on neighbors), neither bond-id nor `:bonds`-list order; it is
identical for the atom view and the bond view. A `stereo_carriers()` accessor (atom view / bond
view) returns this grounded ordered list, making the count and the constraint's local frame
explicit to the user. SMILES's implicit-H-first ordering is converted to H-last at `raise`.

(`!` "not stereogenic" is an atom/bond *assertion*, not an element config — see the constraint
section. SMILES `@`/`@@`/`/`/`\` are parser tokens mapped to a coset at `raise`, never a DSL
config value.)

Keyword sugar — **self-contained** (carries the class): `:ccw` = `:type "Th1"`, `:cw` =
`"Th2"`, `:e` = `"Ct1"`, `:z` = `"Ct2"`.

```clojure
; bromochlorofluoromethane — tetrahedral; 4 explicit carriers incl. the implicit-H sentinel
{:atoms [[:C "C#h1"] [:F "F"] [:Cl "Cl"] [:Br "Br"]]
 :bonds [[:C :F :single] [:C :Cl :single] [:C :Br :single]]
 :stereo [{:id :s1 :focus :C
           :carriers [:F :Cl :Br [:implicit-hydrogen :C]]
           :type "Th1"}]}                          ; or :type :ccw

; 1,1-dichloroethene =CH₂ — prostereogenic: TWO distinct H carriers on C2 (the case that
;   forces explicit non-bond carriers; a count-based scheme cannot express it)
{:atoms [[:C1 "C"] [:C2 "C#h2"] [:Cl1 "Cl"] [:Cl2 "Cl"]]
 :bonds [[:C1 :Cl1 :single] [:C1 :Cl2 :single]
         {:id :bd :a :C1 :b :C2 :type :double}]
 :stereo [{:focus :bd
           :carriers [:Cl1 :Cl2 [:implicit-hydrogen :C2] [:implicit-hydrogen :C2]]
           :type "Ct+"}]}                           ; stereogenic, config unspecified

; relative stereo via binds: d/l share ?r; meso is ~?r ("the other")
{:stereo [{:focus :C2 :carriers [...] :type "Th?r :: {1,2}"}
          {:focus :C3 :carriers [...] :type "Th(~?r)"}]}   ; meso

; the SAME centre asserted as an atom CONSTRAINT (local frame, user-writable, distinct notation):
;   atom-string  "C#h#T1"  — carbon, one implicit H, tetrahedral config 1
```

The config index is the OpenSMILES arrangement index **relative to the explicit `:carriers`
list order**. This is **equivariant storage**, not canonical (D2, refined): a structural edit
relabels the carrier refs with the config untouched (the config is over list *positions*),
attribute mutation never moves carriers. SMILES `@`/`@@`/`/`/`\` are parser tokens mapped to a
coset at `raise`; the index-independent canonical form is **derived on demand**
(equality-up-to-isomorphism, canonical serialization, meso) — never the storage frame. `~`
("the other") is frame-independent for the binary classes (the unique non-identity element,
regardless of carrier order), which is why `~?r` expresses meso under any listing. The same
index form (`5` etc.) round-trips higher-coordination centres (Sp/Tb/Oh) without perceiving
them.

**Undetermined (`*`) vs Unknown (`+`).** Deliberately distinct, mirroring `#a*` vs `#a+`:
`*` is the lattice top (no assertion; matches anything; the site may not be a stereocenter at
all), `+` is a positive assertion that the site **is** stereogenic with configuration
deliberately either/unspecified (the wavy bond, molfile wedge-4, CX `w`). The distinction is
required for round-trip fidelity — a wavy bond must serialize back as wavy, not as "no
stereo" (information lost) nor as the explicit racemic set `{1,2}` (meaning changed).

### Configuration operations (the config algebra)

Configurations are cosets in `R\Sₙ` (`Sₙ`/proper-rotation group). The operations that move
between them are **group actions on the coset space, not integer arithmetic on the index** —
`3 - ?r` is rejected because the coset index is an arbitrary labeling and arithmetic on it is
not invariant. There is one generic operation and one distinguished shorthand; both are
**defined in the contract now** (so the notation is forward-stable), with the per-class
*coding* filled for the binary classes (TH/CT) and reserved for the rest:

- **`^` (generic, infix):** `config ^ k` applies effective-group operation `k` (GAP's `x^g`
  action convention). The operation space is the class's **effective operation group** — `Sₙ`
  modulo the core of `R` (largest normal subgroup of `Sₙ` inside `R`); raw `Sₙ` Lehmer codes
  over-parametrize (many permutations act identically on cosets). This is the small
  class-specific group; only `^1` is exercised by the deliverable:

  | Class | R | effective operation group | size | codes implemented now |
  |---|---|---|--:|---|
  | TH | A₄ | S₄/A₄ ≅ Z₂ | 2 | `{0: id, 1: the other}` |
  | CT | — | Z₂ | 2 | `{0: id, 1: E/Z swap}` |
  | SP | D₄ (core V₄) | S₄/V₄ ≅ S₃ | 6 | reserved |
  | TB | D₃ | Sₙ on 20 cosets | — | reserved |
  | OH | O | Sₙ on 30 cosets | — | reserved |

- **`~` (the predefined involution, unary):** `= ^1`, where code 1 is the class's **canonical
  generating involution** — *"the other configuration."* Uniform across classes (the operation
  is "go to the other coset," not "reflect"), visible, involutive (`~~ = id` by construction —
  the coding assigns an involution to code 1; for binary classes it is the unique non-identity
  element). This is the **only operator the deliverable exercises**:
  - TH: `~?r` = the other config = the enantiomer → meso (`?r` homochiral, `~?r` meso).
  - CT: `~?r` = the other config = the E/Z swap.

  Choosing `~` = "the other" (rather than `~` = reflection) is deliberate: a *reflection*
  glyph would be the **identity on CT** (E/Z are achiral — their mirror is themselves), a
  footgun where `~` silently does nothing. "The other" is non-trivial and uniform for both
  binary classes.

`~~?r == ?r` everywhere; for binary classes `~` is the unique non-identity element, so it is
frame-independent (Z₂ has one "other" regardless of labeling). `~?r == ?r ^ 1` by definition.
All operators distribute over the channels (`~*`/`* ^ k` = `*`, same for `+`); composition
`?r ^ a ^ b` applies `b ∘ a`.

A **distinct mirror / reflection operator** (the improper point-group op, `= ^ reflection_code`)
is *not* needed now — for TH it coincides with `~`, for CT it is the identity, and it only
diverges from `~` for multi-config classes (octahedral / TB coordination complexes). It can be
added later as another shorthand over `^` **without any breaking change**.

Precedence (unary `~` binds tighter than infix `^`) is conventional; glyph details are
settleable without blocking, since only `~` is exercised by the deliverable.

**Operations stay scalar — they are not structural.** An operation is named by its code,
*not* by carrier references (`swap :F :Cl` would pull structural refs into the scalar config
string, violating the separation principle). The code references *positions* in the explicit
`:carriers` order, the same coupling the config value already has — no new structural
dependence. The only genuinely structural permutation is a carrier-*naming* exchange used as a
**transformation/edit** (a reaction rule that says "exchange these two specific substituents");
that is a different operation from expressing a config value and lives in the transformation
layer, never in the config string.

**`:group`** (sample-layer ABS/AND/OR) stays a map key for I/O round-trip, but it is not the
main assertion mechanism. The native assertion mechanism is ordinary variable sharing and config
operations across elements (`Th?o`, `Th~?o`, etc.); that is strictly more compositional inside
umol's existing query/rule term algebra.

### Ordering, the atom/bond constraint, and perception

**Storage frame (equivariant).** The configuration is the orbit relative to the **explicit
`:carriers` list order**. A structural edit is a clean relabel of the carrier refs with the
config (over list *positions*) untouched; attribute mutation never moves carriers. The
atom/bond-numbering trade extended to stereo — equivariant in storage, canonicalization
external.

**Adjacent-field precedents — chemistry is not unique.** Stereo is a *combinatorial map /
rotation system* (a.k.a. ribbon graph, fat graph): a graph plus, at each vertex, the cyclic
order of incident half-edges (darts). The order is irreducible — there is no order-free
encoding of the full configuration. The same structure recurs:
- **AFMS** (Andersen–Flamm–Merkle–Stadler, *Chemical Graph Transformation with
  Stereo-Information*, 2017): the "ordered list method" — an order on each vertex's incident
  edges, permutation groups giving equivalence classes; supports *partial* stereo. The direct
  precedent for stereo in graph rules.
- **Knot diagrams**: PD (planar-diagram) notation is a per-crossing *fixed-length 4-tuple* of
  arcs in a canonical local rotation; Gauss code is a *signed* sequence (over/under). The
  knot's over/under is our cis/trans / R/S.
- **Rotation systems** (topological graph theory): a pair (σ, ρ) of dart permutations — the
  canonical math of graph + local orientation.
All of them either localize the rotation (AFMS list, PD tuple, SMARTS string) or hoist it to a
molecule-level carrier list (CDK). The order is always present.

**The atom/bond constraint — a scalar parity, and the dichotomy.** The full carrier-ordered
configuration lives only in the molecule-level element (it is a rotation datum — listing it on
the atom would be the variadic entity-ref list we reject). What *can* project onto the
atom/bond DSL as a **scalar, ref-list-free** field is a **parity**, and there is a genuine
dichotomy in which frame it is relative to:

| form | order is… | scalar / no ref-list | index-independent | locality |
|---|---|---|---|---|
| **incidence-order** (SMARTS `@`/`@@`) | the incident bonds' order in `:bonds` | yes | **no** (equivariant) | **local** — matchable on a fragment |
| **canonical** (R/S-style; graph-canonical parity now, CIP later) | a canonical rank | yes | **yes** | **global** — rank needs the whole molecule |

Fundamental, not a gap: index-independence needs a *total* canonical order on the neighbors,
whose tie-breaking needs global information — so a scalar local descriptor is either
index-equivariant (SMARTS) or index-independent-but-global (R/S), never both. Either way the
field is a **scalar parity** with **no entity-ref list** (the order comes from `:bonds`, not a
list inside the atom field), so it satisfies the hard constraint and is **symmetric** with a
derived predicate (computed parity ↔ asserted parity), exactly like `#R`/`#D`.

**Resolution — incidence frame, uppercase per-class tags.** Canonical-rank parity is
**rejected**: a canonical index is unstable under *both* topology edits and attribute mutation
(any recolor that permutes the rank changes it), so it would break the "edits are remaps"
invariant every other field holds. The constraint uses the **incidence-order** frame
(equivariant, stable, SMARTS-faithful) — the same trade as equivariant element storage.

- **Distinct notation, by the existing case convention.** The constraint is an **assertion about
  a projection** of the element onto the atom's local incidence order — and uppercase `#`-tags
  are already the *derived-predicate* namespace (`#D` degree, `#X` connectivity, `#H` total-H,
  `#R` ring count).
  So stereo constraints are **uppercase per-class tags**, one per class, **mirroring the head
  class token** (`Th`→`#T`, `Ct`→`#C`): `#T` tetrahedral, `#C` cis-trans stereobond; reserved
  `#A` allene, `#Q` square-planar, `#B` trigonal-bipyramidal, `#O` octahedral (deferred letters
  settleable; guideline: mirror the class token, disambiguate collisions). The uppercase tag and
  position make the **frame** unmistakable while paralleling the element: `#T1` (atom-string,
  local-incidence frame) ↔ `:type "Th1"` (element, explicit-carrier frame) — same config
  algebra (`1` / `{1,2}` / `?r` / `~?r` / `*` / `+`), different frame. (Letters coexist with
  lowercase by case: `#T`/`#t`, `#C`/`#c`, `#O`/(no `#o`); none collide with the taken
  `#D`/`#X`/`#H`/`#R`, and joint-domain `#E` is untouched.)
- **User-writable transport and assertion.** `#T` / `#C` are constraints, not storage fields.
  They have two jobs: (1) a compact transport layer for SMILES/SMARTS/MOL local stereo notation,
  and (2) a way to assert a local projection of a molecule-level stereo element. `C#h#T1` is a
  valid atom-string (carbon, one implicit H, tetrahedral config 1). For a ground molecule, a
  validator can derive the corresponding local projection from `:stereo` and compare it against
  the asserted `#T1`; for a pattern, `#T1` is matched SMARTS-style by mapping the pattern's
  local incidence order into the target's local incidence order. Thus `#T/#C` are parallel to
  derived atom predicates: computed from elements when available, assertable in patterns, and
  cross-checkable on ground data.
- **Survives ports.** The projection is onto the atom's incident *bonds* now, incident *ports*
  later — the scalar tag+coset never lists carriers, so nothing in `#T<coset>` changes.

`!` "not stereogenic" lives here too (a coarse, automorphism-derived atom assertion). **Remaining
sub-decision:** the precise local order — incident bonds in **bond-id** vs **`:bonds`-appearance**
order — and the placement of special substituents (implicit-H / lone-pair) within it.

**Perception — `raise` builds elements directly from TableIR.** There is *no* intermediate
carrier-ordered atom constraint. TableIR stays the faithful atom-centric input
(`atom.chirality` = `@` relative to the input neighbor order, which it preserves); `raise`
reads that order and emits the molecule-level element directly: focus = the atom, **carriers =
the input-order neighbor atoms + explicit virtual sentinels**, config = the coset for `@`
relative to that listed order. So the atom-centric → molecule-centric bridge is `raise` itself,
with no side channel and nothing stored on the atom. This is the disanalogy with aromaticity:
`#a` works because π-contribution is a *genuine coarse projection*; a stereo arrangement has no
coarse projection (it is the irreducible rotation), so there is no atom-level config — only the
molecule-level element, plus the scalar-parity constraint above for matching and the symmetric
framework. Stereogenicity / achiral / prochiral classification runs during/after `raise` on the
elements (+ automorphism).

**Derived views.** Canonical labeling (`AtomAutomorphism::canonical_labeling`), CIP, meso, and
stereogenicity are computed **on demand**, never the storage frame. Para-stereocenters use a
refinement **fixpoint** (fold parity into the colour, re-refine to convergence). Equality
up-to-isomorphism canonicalizes both sides; structural `PartialEq` works directly on the
equivariant stored form.

### Achiral, nonstereogenic, and prochiral sites

Three distinct cases, handled at different layers (meso is a fourth, molecule-level case,
covered under decisions above):

- **Achiral configuration** (defined, its own mirror) — stored as an element with achiral
  parity (`0`); the isomer distinction lives in the **carrier arrangement** (the orbit),
  not a handedness sign. In scope this is exactly the **cis/trans stereobond**: an entry with
  `:type "Ct1"`/`"Ct2"` is parity-0, with E-vs-Z carried by the arrangement, and equality
  compares modulo the permutation group with **no inversion orbit**. `~` on a CT element is the
  E/Z swap (not a no-op — `~` is "the other config", not a reflection). Tetrahedral centres are
  the chiral (`±1`) case.
- **Nonstereogenic atom** (defined shape, one distinguishable arrangement — `CH3`,
  `CH2XY`, an `AABC` centre) — **no element emitted.** The automorphism pass detects it:
  the site stabilizer `A_s` identifies the candidate configurations because swapping the
  redundant ligands is an automorphism. This is the `is_actually_chiral` test; we take the
  lean route (emit nothing) rather than instantiating a one-assignment shape object.
- **Prochiral atom** (nonstereogenic now, enantiotopic ligands — ethanol `CH2`, glycerol
  C2) — **not stored.** Prochirality is a derived topicity classification from the
  `Â`/`Â*` action over `AtomAutomorphism`: ligands swapped by a proper automorphism are
  homotopic, by an improper one only are enantiotopic (⇒ prochiral), by neither are
  diastereotopic. It is a query capability — the same machinery as stereogenicity — not
  part of the `:stereo` layer.

Forward-looking: one stronger model is to instantiate candidate shape objects at every plausible
site, including nonstereogenic/prochiral ones, so that a prochiral→chiral desymmetrization is a
state change rather than the creation of a new element. That is not necessary for this overlay.
With **explicit H atoms**, the prochiral ligand sites are already structural atoms and topicity is
derivable from the automorphism action. The gap is mostly **implicit** sites: two implicit H
slots on a `CH₂` atom are not structural atoms, so a query that wants to talk about "replace one
of the enantiotopic hydrogens" needs either (a) temporary virtual carrier slots during matching
or (b) a stored candidate-site layer. The chosen overlay uses (a): no element is stored for a
nonstereogenic/prochiral site, but perception/matching may materialize virtual H/lone-pair slots
when evaluating topicity or stereo-forming patterns.

### I/O and 3D

- **Parsing → molecule-level elements** (detail in "Ordering, the atom/bond constraint, and
  perception"). At `table_ir::raise`, `atom.chirality` (`@`/`@@`), `bond.wedge` (`/`/`\`), and
  MOL parity are read **directly into molecule-level elements** — focus, explicit carriers (the
  input-order neighbor atoms + virtual sentinels), config from the input order — with no
  intermediate atom-level carrier-ordered constraint. TableIR `StereoInterpretation` +
  `cx_data.stereo_groups` populate specified-ness (`*`/`+`) and `:group`. The scalar-parity
  atom/bond constraint (dichotomy above) is for matching / the symmetric framework, derived
  from the elements, not the transport vehicle.
- 3D conversion: the configuration drives local geometry during embedding (signed tetrahedron
  volume for centres; substituent side for stereobonds).

## The maximal port model (keep in view; can cut scope)

The "do-it-right" target the spike probes. One entity type — an electron-weighted set of
ports — with port-ownership:

- **Ports** are undirected half-edges. Per-port electron contribution encodes donation
  direction (dative = contributions `(2,0)`), so no darts.
- **Entities** are variadic electron-weighted port sets; localized bond / dative /
  aromatic / multicenter / lone pair / radical are *classifications* by arity + electron
  pattern, not separate structures. The current four relation kinds collapse into one table.
- **Lone pairs** are monovalent entities (one port, 2 e⁻, no partner); **unpaired
  electrons** are singly-occupied ports (1 e⁻). Both leave `AtomAst`. A lone pair's port
  re-pairing with an acceptor port *is* dative-bond formation.
- **Ports are polymorphic** — a port carries a localized-valence order; a double bond is
  one order-2 entity over one port-pair. No σ/π decomposition.
- **Port ownership is `port → atom` or `port → opaque fragment`.** Fragment-owned group
  ports are ligand group orbitals (SALCs) — the principled treatment of haptic bonds
  (η⁶-Cp⁻), tying into the umol-msym SALC machinery.
- **Valence ports** (bounded minimal basis, electron-conserving) vs **interaction ports**
  (non-covalent, unbounded, no valence budget) — the existing valence/non-covalent split
  re-expressed.
- **The atom-level graph is a derived, cached projection** (contract ports per atom).
  Classical polynomial algorithms run on the projection; only electron-counting, stereo,
  and reactions consult the full port structure. Ports are the normal form; the atom graph
  is a functor of it. This is the guardrail against orbital-centering every algorithm.
- Lineage: Ugi–Dugundji BE-matrices (reaction algebra), combinatorial maps / rotation
  systems (the half-edge structure), SALCs (delocalized/haptic), Andersen–Flamm–Merkle–
  Stadler (stereo-aware graph rewriting).

## Contract steps to do now

- **Add the `:stereo` top-level key** (class in the config head; grammar above) and the
  molecule-level `StereoElement` semantic type, plus the uppercase per-class atom/bond
  constraints (`#T`/`#C`, user-writable). Nothing carrier-ordered is *stored* on `AtomAst` /
  `BondAst` (the firewall) — the constraint is a scalar-parity projection (dichotomy above).
- **D1 — carrier slot is an extensible typed `carrier` sum** (`atom-ref | [:implicit-hydrogen atom-ref] |
  [:lone-pair atom-ref]`; reserved typed arms: `[:bond bond-ref atom-ref]`, `[:port port-ref]`,
  `[:fragment fragment-ref]`, `[:opaque ...]`). Implement the atom + virtual-H arms now.
  Unknown carrier tags are rejected unless they use the explicit `:opaque` preservation form.
  This is what lets ports slot in additively later without making today's grammar ambiguous.
- **D2 — store configuration as the equivariant `Sₙ/R` coset relative to the
  explicit `:carriers` list order**, not the graph-canonical labeling. Equivariant: structural
  edits relabel carrier refs with the config untouched; attribute mutation never moves carriers.
  Not canonical-in-storage — the atom/bond-numbering trade extended to stereo. Canonical
  labeling / CIP / the `A_s` quotient are **derived on-demand** views (equality-up-to-iso,
  serialization, meso), never the storage frame. Not a double-coset, not CIP-in-storage.
- **D3 — when the rule (SMIRKS) layer is specified, give it a total, explicit lhs↔rhs
  correspondence** (the shared-`:id` mechanism already does this). Stereo transfer semantics are
  deliberately not designed here; the correspondence only keeps that layer additive.

## Forward-looking design (smooth port transition)

How each maximal-model aspect reaches the contract, given D1–D3:

| Aspect | Contract impact |
|---|---|
| Carriers as ports | **Additive** — port-ref arm of the `carrier-ref` sum (D1). |
| Stereo-aware canonicalization (fixpoint) | **Invisible** — a derivation; storage unchanged. Version only persisted canonical keys. |
| Reaction stereo transfer | **Additive** — `:stereo` on lhs/rhs + explicit D3 correspondence leave room for later transfer semantics; no semantics fixed here. |
| Stability axis (configurational/conformational/fluxional) | **Additive** — new optional element field. |
| Assertion layer | **Mostly already present** — variable sharing and config operations express relative assertions; `:group` is only ABS/AND/OR I/O sugar. |
| Port/entity molecule model (lone pairs/radicals off atoms; bonds re-based on ports; opaque fragments) | **Breaking to the molecule model only.** Stage via new optional `:ports`/`:entities` keys + dual-encoding + a `:model` version gate. **The stereo contract is invariant across it** because carriers are refs (D1) and config is a model-agnostic orbit (D2). |

The linchpin is the `carrier-ref` sum: an atom-ref carrier becomes a port-ref under the
same slot, and the orbit configuration does not care whether the connectivity beneath it is
atom-based or port-based. So the only genuinely model-breaking item (the entity/port model)
never reaches the stereo layer.

## The `umol-port` spike

- **Form**: a separate experimental crate, off the deliverable critical path. Hard timebox
  (≈2 days of real effort; treat optimistic estimates as ~10× under).
- **Goal**: validate the entity–port incidence storage + the atom-graph projection, the
  electron model (lone pairs / radicals as entities; BE balance), and two hard cases —
  haptic bonding via fragment-owned group ports, and a lone-pair centre.
- **Success criteria**: (1) the atom projection is cheap enough that classical algorithms
  are not materially slowed; (2) electron bookkeeping balances on representative
  organics + an organometallic; (3) a stereo element's carriers map to ports cleanly,
  confirming the D1 ref-sum is sufficient.
- **Risk (carried, not mitigated)**: pioneering — no open-source toolkit commits to a
  port-factored electron model as storage with a full algorithm suite. The lineage exists
  in pieces (BE-matrices, combinatorial maps, SALCs, stereo-aware rewriting); the integrated
  performant implementation does not.
- **Gating**: migration of the molecule model to ports happens **only if the spike meets
  its success criteria**. Otherwise the overlay remains the model and the contract (built
  for both) loses nothing.

## Superseded directions (record of what was tried and rejected)

The design converged through several drafts; these are kept so the reasoning is not relitigated.

- **`#o` config predicate + `#g` geometry predicate (headless stereo-string).** Rejected: the
  configuration is *definitional*, so it belongs in the **head**, not a `#`-predicate (heads are
  the inherent definition; predicates assert about an already-defined entity). Class travels with
  config in the head (`"Th1"`).
- **Class as the top-level key (`:tetrahedral` / `:cis-trans`, per-class keys).** Considered and
  superseded by **class in the config head** (`"Th1"`, one flat `:stereo` key): it parallels the
  atom/bond constraint `#T1`/`#C1` and makes the config head self-describing (the `:cw`/`:e`
  sugars carry their class; `"Th1"` is interpretable alone). Per-class keys would have been more
  like `:bonds`/`:dative`, but lose self-containment. (No-class-algebra holds either way — the
  head class is a fixed enum, not a value-expr.)
- **Derived carriers (`:carriers` optional, inferred from focus + incident bonds + counts).**
  Rejected: requires consulting the graph (elision); does not generalize to allene/axial/planar
  (carriers ≠ the focus's incident bonds); an optional-here-required-there rule is incoherent;
  and the count-based scheme cannot express prostereogenic cases (1,1-dichloroethene's two `=CH₂`
  hydrogens need two **distinct** carrier slots). Carriers are always explicit, with non-bond
  carriers as distinct typed vectors (`[:implicit-hydrogen ref]` / `[:lone-pair ref]`).
  A prior vector-with-type-string form (`[ref "#h"]`) was replaced because typed keyword heads
  make future carrier arms unambiguous.
- **Bonds (bond-refs) as the real carriers.** Argued early for multi-edge disambiguation +
  port-readiness; withdrawn. Multi-edge does not arise in umol's `:bonds` graph (order-n edges,
  not parallel edges), so a neighbor atom names the direction unambiguously; atoms also win on
  readability, optional-bond-ids, neighbor-order alignment, and sentinel uniformity; and atoms
  migrate to port-refs cleanly. Real carriers are **neighbor atom-refs** in the compact surface;
  the reserved `[:bond bond-ref atom-ref]` arm is for the narrower case where a stereochemical
  carrier is an incidence/half-edge rather than just the neighboring atom.
- **Bond-id-order storage frame.** Superseded by the **explicit `:carriers` list order** (now
  that carriers are always explicit), which makes structural edits a clean ref-relabel with the
  config untouched.
- **"The atom/bond constraint = the element, one frame" (constraint as the carrier-ordered
  config on the atom).** Rejected: that is either a variadic entity-ref list on the atom
  (forbidden) or "twice the same info written differently." Replaced by the **scalar-parity
  constraint** in the **incidence frame** as an **uppercase per-class tag** (`#T`/`#C`; a derived
  projection of the element), with the carrier-ordered config living only in the molecule-level
  element, and `raise` building elements directly from TableIR.
- **Canonical-rank frame for the constraint.** Rejected: a canonical index is unstable under
  topology edits *and* attribute mutation (recolor permutes the rank), breaking the
  "edits are remaps" invariant. The incidence frame (equivariant) is used instead.

## Sources (this revision)

- [SMARTS chirality semantics — SMARTS101](https://smarts101.dev/how-to-smarts) · [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html) (chirality = neighbor *order*, not R/S; unspecified matches both)
- [Andersen, Flamm, Merkle, Stadler — *Chemical Graph Transformation with Stereo-Information* (2017)](https://link.springer.com/chapter/10.1007/978-3-319-61470-0_4) (the "ordered list method" + permutation groups; partial stereo)
- [Rotation system — Wolfram MathWorld](https://mathworld.wolfram.com/RotationSystem.html) · [Rotation system / combinatorial map — Wikipedia](https://en.wikipedia.org/wiki/Rotation_system)
- [Planar diagram (PD) notation — Knot Atlas](https://katlas.org/wiki/Planar_Diagrams) · [Knot notation — KnotInfo](https://knotinfo.math.indiana.edu/descriptions/notation.html)
