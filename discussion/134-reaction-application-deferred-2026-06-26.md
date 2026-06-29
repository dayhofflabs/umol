# 134 — Reaction-application: deferred items

Five items were left to work out while landing the `Edit`/`Delta`/`Undo` vocabulary refactor
(`Add`/`Modify`/`Remove` ↔ `Added`/`Modified`/`Removed`, key-based `ModifyConstraint`, by-value
molecule `Add`/`Remove`, `CascadedConstraints`) and designing the reaction EDN surface (doc 133) —
two reaction-application features, plus structural entity refs (item 3) and a constraint-store
cleanup (item 4), which surfaced while implementing 133, plus a serde streaming-parser / naming
review (item 5) that surfaced while wiring the R12 tests. The vocabulary/transact groundwork is
landed and green; all are increment-2, captured here so they are not lost, to be tackled after the
reaction DSL (133) is implemented.

## 1 — molecule-level constraints don't apply in reactions

**State.** `ReactionAst::apply_at` has `Delta::Constraint(_) => {}` — a molecule-level
`ConstraintDelta` (add/remove on the flat `Vec<Constraint>`) is **silently dropped** on apply.
(The span side is resolved separately: `ReactionSpanAst` gains a `constraints: Vec<ConstraintSpan>`
field and `to_reaction_span` populates it — doc 133, decision 14. This item is the *operational*
`apply_at` path only.)

**Mechanism is ready.** The transact layer now supports general by-value
`Edit::Add/RemoveMoleculeConstraint` (true multiset; `Remove` is last-match-by-value with the
position captured for rollback via `CascadedConstraints`). So this is *closable* — only the
lowering is missing.

**Blocker.** Applying a rule-frame constraint to the host needs its atom/bond refs remapped through
the match `m` (lhs atom → host atom, created atoms → appended). `IdRemapping` / `Constraint::remap`
only express *removal-compaction* (old→new after deleting ids); the match is an arbitrary
injection, so there is no existing remap that fits.

**Proposed solution.**
- Add a total match-based ref-map: `Constraint::map_topology_refs(self, atom: impl Fn(AtomId) ->
  AtomId, bond: impl Fn(BondId) -> BondId) -> Self`, mirroring `Constraint::remap`'s tree traversal
  (entity-leaves, `Relational`, `Molecule`-scope `:atoms`/`:bonds` lists and sub-pattern anchors,
  `And`/`Or`/`Not`). It is total (no drops). Overlay entity-leaf refs (dative/aromatic/multicenter/
  noncovalent/stereo) are out of increment-1 — the embedding maps only atoms/bonds — so those error
  (and are subsumed by item 2).
- Wire `apply_at`: for `Delta::Constraint(ConstraintDelta::Add/Remove(c))`, emit
  `Edit::Add/RemoveMoleculeConstraint { c.map_topology_refs(lhs→host via m, created→appended) }`,
  ordered **before** `RemoveTopology` so transact's own renumbering (`remap_with_update`) finishes
  the job (and drops the constraint if it referenced a removed atom).
- Test: a reaction with a molecule-level constraint delta applies and round-trips through rollback.

Smaller of the two: one remap method + the `apply_at` arm.

## 2 — overlay entity `Delta` infrastructure

**State.** The AST `Delta` is `Atom(AtomDelta) | Bond(BondDelta) | Constraint(ConstraintDelta)`
only. The six overlay entities — dative bond, aromatic system, multicenter bond, noncovalent bond,
stereo atom, stereo bond — have **no `*Delta`**. So a reaction cannot add / modify / remove an
overlay entity at all.

**Mismatch.** The reaction EDN surface (doc 133) and the `apply` / `to_reaction_span` design are
written for all eight entities, but the AST `Delta` covers three (atom / bond / constraint). This is
the increment-1 / increment-2 split in the localized-topology scoping, left to work out.

**Proposed solution.**
- Add the six overlay `*Delta` enums on the `EntityDelta` pattern: `Add` / `ModifyField` /
  `ModifyConstraint` / `Remove`, extend `Delta`, and the `apply_at` / `to_reaction_span` lowering.
- The `EntityDelta` generic path (the `diff` / `deltas_from_states` default methods + the
  `field_ops!` macro) means each new entity is mostly mechanical — per entity only: its `*Delta`
  enum, the `EntityDelta` impl (the `(variant ⇒ field)` map), the `Atoms` type, `into_delta`, the
  `MoleculeAst` field, the `Edit` variant, and its `EntitySpan` span slot.
- **Span generalization (tied to 133 work item 2).** Every `MoleculeAst` entity collection is
  parameterized by its AST payload type `R`: atoms/bonds as `Vec<R>`, the six overlays as `*Relation`
  / `*Birelation` sets that split a topological part (participants + incidence over `NodeId` /
  `EdgeId`) from a `data: Vec<R>` payload column. `ReactionSpanAst` is the identical structure with
  that payload lifted `R → EntitySpan<R>` in every collection — topology shared, only the per-entry
  value becomes a span; `EntitySpan<T>` itself is unchanged. It's done for the two vecs; lifting the
  six overlay relation sets' `data` columns to `EntitySpan<…>` is the remaining work.

Larger of the two: six overlay deltas + the span generalization.

## 3 — structural entity refs

**State.** Every `<entity>-ref` is `int | keyword` (position or id) — reaction surface (doc 133) and
spec (§7.9). A bond/overlay with no id can only be named by position; to reference it otherwise you
must give it an `:id`.

**Want.** Name an entity by its constituents instead — a bond by its endpoints, an aromatic /
multicenter system by its members, a dative bond by donors+acceptor, a stereo element by
site+ligands. (Atoms are the base; no structural form.)

**Form — explicit map, uniform.** `<entity>-ref ::= int | keyword | <structural-map>`, the §4 entry
form minus `:type`/`:id`:

| entity | structural ref |
|---|---|
| bond, noncovalent-bond | `{:atoms [atom-ref atom-ref]}` |
| aromatic-system, multicenter-bond | `{:atoms [atom-ref+]}` |
| dative-bond | `{:donors [atom-ref+] :acceptor atom-ref}` |
| stereo-atom | `{:site atom-ref [:ligands [atom-ref+]]}` |
| stereo-bond | `{:site bond-ref [:ligands [...]]}` |

Map form (not the bare `[atom-ref …]` vector) so a structural ref is self-delimiting and never
collides with the vectors it nests inside (anchor pairs `[[ref ref]+]`, relational `[ref target]`).
Checked: the forms paired with a ref never key on `:atoms`/`:donors`/`:site`, and stereo forms are
vectors, so `[<ref> <form>]` stays unambiguous by shape/key.

**Sites.** One production, so it reaches every non-atom ref at once: reaction `:remove`/`:modify`;
entity constraints; relational leading refs plus the embedded `bond-ref` in `:dative-bond-parallels`
/ `:stereo-bond-site`; `:bond-order-sum :bonds`; anchor pairs; stereo-bond entry `:site`. In an
entity constraint a structural ref is the *query* form (point at an id-less entity by its parts); the
relational leaves still assert membership of a named entity.

Extending a production propagates everywhere automatically:

  - A. Reaction deltas (doc 133, not yet in spec): :remove <ref>, :modify <ref> …. ← the proposed site.
  - B. Entity constraints (§7.9 677–683): {:bond [bond-ref form]} + the 6 other non-atom entities.
  - C. Relational constraints (§7.9 686–717): the leading <entity>-ref of every leaf, plus two embedded cross-refs to a bond-ref: :dative-bond-parallels [dative-bond-ref bond-ref] and
  :stereo-bond-site [stereo-bond-ref bond-ref].
  - D. Molecule-scope (§7.9 722): :bond-order-sum {:bonds [bond-ref+] …}.
  - E. Anchor-spec (§7.9 794–800): pair-lists [[<entity>-ref <entity>-ref]+] for bonds, dative-bonds, aromatic-systems, multicenter-bonds, noncovalent-bonds, stereo-atoms, stereo-bonds.
  - F. Stereo-bond entry :site (§4 line 148): a bond-ref to an existing bond.

**Cost.** Not a `define_ref` tweak. The macro emits a uniform `enum { Index, Id }` resolving over an
`id → index` map. The structural variant has a per-entity payload (`[AtomRef;2]` / `Vec<AtomRef>` /
donors+acceptor / site+ligands) and resolves by matching the host's entity *collections* by
constituents — a different payload and resolution per entity. §4.1 uniqueness (no two bonds on a
pair, aromatic systems disjoint, …) makes the match ≤1 hit, so the semantics are clean; the code
shape is the work.

## 4 — entity constraints: uniqueness-by-key explicit

**State.** Each entity constraint store (`AtomConstraints`, plus the bond / overlay analogs) is
addressed two ways: by `kind` (`AtomConstraintKind`), which is **multi-valued** — `RingMembership`
has one entry per size — and by `key` (`AtomConstraintKey`), which is **always unique** (a kind's
subkeys distinguish its entries). `add` calls this distinction `is_unique` ("unique kinds replace,
ring appends") — misleading: by key *everything* is unique; "non-unique" only means a kind has
subkeys.

**Want.** Redesign the constraint-store API so uniqueness-by-key is the clear primary invariant —
every entry is addressable by a unique key, a kind is just a group of keys — and retire the
kind-centric `is_unique`/"unique" nomenclature. Applies to all entity constraint stores.

**Now.** `AtomAst::update` (reaction `:modify` resolution) already works per key via `remove_by_key`;
this item is the broader API/nomenclature cleanup, not blocking.

## 5 — FromEdn streaming-parser audit + parse-direction naming review

Two related reviews surfaced while wiring the R12 tests (doc 133).

**5a — streaming `from_edn_str` audit.** Per dsl-serialization, every EDN-shaped `*Dsl` must override
`FromEdn::from_edn_str` to drive its hand-written streaming `read_*` — never the trait default, which
delegates to the tree (`read_string` + `from_edn`), violating "streaming must not delegate to tree."
Today only the top-level / full-entity types do (`MoleculeDsl`, `ReactionDsl`, `AtomDsl`, `BondDsl`,
`DativeBondDsl`, … via `read_subgrammar_all`); the constituent boundary types inherit the
tree-delegating default. Audit every `FromEdn` impl and add the override where missing. The full
type list and per-type wiring notes are in doc 133's R12 deferred note; in summary:

- Most delegate to an existing single-arg `read_<type>_dsl` (mechanical); a few readers
  (`read_topicity`, `read_stereogenicity`, the stereo constraint readers) are private and need
  `pub(super)`; the ref and stereo constraint `*Dsl` macros need the override in the macro body.
- `read_molecule_constraint_dsl(de, key)` and `read_relational_constraint_dsl(de, key)` take a
  pre-consumed dispatch key (the `{`+key is consumed by the umbrella `read_constraint_dsl`), so they
  have **no** standalone reader — a thin standalone reader or a small restructure is needed, a design
  decision rather than a mechanical change.
- `fuzz_constraints` calls `ConstraintDsl`/`ConstraintsDsl::from_edn_str` directly, so its "streaming"
  arm is currently the tree path — fixing the override makes that target exercise the real streaming
  reader.

**5b — parse-direction naming.** Review the naming across the two parse directions and decide a
consistent scheme: string → AST (`FromStr`/`parse_<type>` → AST via defaults) versus the EDN-DSL
boundary (`FromEdn::from_edn` tree / `from_edn_str` streaming, `read_<type>` / `read_edn_<type>`,
`*Dsl` vs `*Ast` trait targets). Confirm the `parse_`/`read_`/`read_edn_` split reads unambiguously
for "parse a string into an AST" vs "parse into the EDN/DSL form," and settle the prospective
`read_edn_` → `pick_` rename (dsl-serialization, not yet adopted).

## Sequencing

All five are increment-2, after the reaction DSL (133) lands. Item 1 is independent and small. Item
2 is the larger one and subsumes item 1's overlay-ref limitation. Item 3 rides with item 2 — six of
its seven structural forms are overlay entities and share the resolution rework; bond is the only
localized entity gaining one, and it is unblocked (name an existing bond by `:id` or index until
then). Item 4 is independent API/naming cleanup. Item 5 is independent serde cleanup (5a mechanical
plus the two dispatch-key readers; 5b a naming pass). None blocks the 133 surface design, which is
atom/bond/constraint with index|id refs and complete on its own.
