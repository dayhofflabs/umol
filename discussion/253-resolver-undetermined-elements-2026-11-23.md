# 253 — Undetermined elements in resolution

Status: In Progress
Date: 2026-11-23
Relates: [operations guide](../docs/development/operations.md),
[data types guide](../docs/development/data-types.md)

## Purpose

A molecule can hold an atom whose element is undetermined. SMILES writes one as `*`, and generic
compound records use it for a group that is left open. Resolution of such a molecule returns
underdetermined and no molecule, so the atoms whose elements are known get no hydrogen counts,
lone pairs, or bond defaults either. `*C` is the smallest case.

The [partial resolution discussion](#partial-resolution) defines an opt-in policy for wildcard
entities. The default still requires a concrete result. The [staged implementation plan](#staged-implementation-plan)
follows the earlier element-only design, which is retained for reference.

## Current behavior

`Resolver::resolve` runs isotope completion, valence admission, aromatic selection, selection among
completions, the constitution edits, stereo, and the defaults of localized and multicenter bonds.
The isotope plan is applied whatever its outcome. Each later phase that is underdetermined or
contradictory ends resolution with that outcome.

| Phase | With an atom of undetermined element |
| --- | --- |
| Isotope | Reads only the isotope field. `Natural` completes an undetermined isotope whatever the element. |
| Valence | `AtomTypingValence::admit` and `CountsValence::admit` each begin by returning underdetermined, with no completions, when any atom has a non-literal element. Resolution of `*C(O)C` ends here. |
| Aromaticity | Not reached in `Resolver`. Called alone, see below. |
| Stereo | Not reached in `Resolver`. Called alone, see below. |
| Bonds, multicenter bonds | Neither reads elements. |

SMILES ingestion gives a bare `*` a natural isotope, charge zero, and no unpaired electrons, and
no aromatic assertion. It reads `[1*]` as an undetermined element with isotope mass 1.

### Final test of `Resolver`

After the last phase, `resolve`, `resolve_with_report`, `resolve_into`, and
`resolve_into_with_report` test `Molecule::is_concrete`. A molecule that fails it is not returned,
and the outcome is underdetermined. The opening comment of `ops/resolve.rs` and the paragraph on
`Resolver::resolve` in the operations guide state the test as the contract.

With the test switched off in a trial, the outcome changes for six kinds of input, in which every
phase finishes and something is left non-concrete:

| Left non-concrete | Example |
| --- | --- |
| An isotope: undetermined under `Strict`, or given as a set or a variable | `C#c0#h4` under `Strict` |
| A stereo entity with an open configuration, created from `#T+` or `#C+` or stored | `C#h#T+` with three halogens; a double bond marked "either" in a MOL file |
| A noncovalent bond with undetermined attributes | two methane molecules with a noncovalent bond `"*"` |
| An atom field that the selected completion leaves undetermined | `C` with a registry whose carbon entry is `C#c0#v0` |
| A bond of undetermined order | `C~C` read with the lenient SMILES syntax |
| A stored aromatic system with undetermined attributes | benzene with `#a` atoms and a system `"*"` |

The trial covered the unit tests of `umol-graph`, its resolution and property suites, and the
tests of `umol-py`. Six unit tests, two property tests, and four Python tests depend on the test;
the resolution suite does not.

### Aromaticity resolver called alone

`AromaticityResolver::plan_resolve` on a six-membered ring of single bonds with five atoms
`C#h#a` and a sixth atom as given:

| Sixth atom | Rule and scope | Outcome |
| --- | --- | --- |
| `*#a` | Hückel, any element | Determined: one system of the six atoms |
| `*#a` | Hückel, Daylight allow-list | Contradictory: `AromaticValenceFailure` at the `*` |
| `*#a+` | Hückel, any element | Underdetermined |
| `*` | Hückel, any element | Contradictory: `AromaticValenceFailure` at a carbon |
| `*#a` | HMO, any element | Contradictory: `HmoInvalidInput`, odd electron count |
| `*#a` | Clar | Contradictory: `ClarNonBenzenoid` |

A `*` without an aromatic assertion that is bonded to a ring atom of benzene written with `#a`
gives a determined plan with the system of the six carbons, under both scopes.

The Hückel rule reads the scope without requiring a literal element: under `ElementScope::Any`
every atom is in scope, and under an allow-list an atom is in scope only with a literal element of
the list. HMO needs parameters of the element and Clar accepts carbon only.

### Stereo resolver called alone

Before S2b, `StereoResolver::plan_resolve` gave the following outcomes with the default
failure policies. The [stereo scope table](#stereo) gives the implemented frame rules.

| Input | Scope of the tetrahedral kind | Outcome |
| --- | --- | --- |
| Stored stereo atom at `*` with four atom ligands | any element | Underdetermined |
| The same | allow-list of carbon | Underdetermined |
| Stored stereo atom at a carbon with `*` among its ligands | any element | Determined |
| `#T+` on `*` with literal hydrogen and lone-pair counts and four neighbors | any element | Underdetermined |
| The same | allow-list of carbon | Underdetermined |

At that point, `derive_stereo_atom`, `derive_stereo_bond`, and `identify_virtual_ligands`
each required a literal element before checking the scope. A non-literal element left the
row unresolved under either scope.

### Constructors

| Resolver | Settings | Constructors |
| --- | --- | --- |
| `Resolver`, `AromaticityResolver`, `StereoResolver` | a config struct | `new(model)`, `with_config(model, config)` |
| `ValenceResolver` | none | `new(model)` |
| `BondsResolver`, `MulticenterBondsResolver` | none | `new()` |
| `IsotopeResolver` | `IsotopePolicy` | `new(policy)` |

`IsotopeResolver` has a setting and no `with_config`; its default comes from the derived `Default`.

## Partial resolution

This section records the agreed design for partial resolution.

### Question

`Resolver` returns determined only for a concrete molecule; the test of `Molecule::is_concrete` at
its end ensures that. This semantics is coherent and too narrow, since it has no partial
resolution. Removing the test, or adding `Admit`, without a definition of partial resolution
leaves open what happens to a molecule with undetermined atoms or bonds and what an
underdetermined outcome means. A plan is of no use before that is defined.

Indeterminacy has two sources. An input can state it, as `*` and `~` do. Or a format leaves
something out: external formats omit hydrogen counts and have no notation for lone pairs. Filling
what is left out is the main purpose of resolution. The DSL and raising do not record the source,
and the internal representation takes no SMILES or MOL semantics, so the two have to be told apart
without it.

### What the phases do today

| Undetermined in the input | Today |
| --- | --- |
| Element of any atom | `AtomTypingValence::admit` and `CountsValence::admit` return underdetermined for the whole molecule, and resolution ends. |
| Order of a bond | `AtomTypingValence::admit` returns underdetermined for the whole molecule. `CountsValence::admit` gives the two atoms no completions and the rest is resolved. |
| Charge of an atom | `CountsValence::admit` gives the atom no completions. |
| Aromatic contribution of a ring atom | `AromaticityResolver::plan_resolve` returns an empty underdetermined plan. |
| Ligand evidence of a requested stereo entity | `StereoResolver::plan_resolve` returns an empty underdetermined plan. |
| Stereo configuration, kind of a noncovalent bond, isotope set or variable, electrons of a stored aromatic system | No phase fills them. Only the final test reports them. |

A phase either ends resolution or passes over the entity. Removing only the final test would
make the outcome depend on which path the phase takes.

### Principal fields

Each entity form has principal fields, a subset of its inherent fields. They are what the DSL
string of the entity begins with. The other fields have no name.

| Entity | Principal field | Other inherent fields |
| --- | --- | --- |
| Atom | element | isotope mass, charge, implicit hydrogens, lone pairs, unpaired electrons |
| Localized bond | order | charge, unpaired electrons |
| Dative bond | order | none |
| Aromatic system | electrons | charge, unpaired electrons |
| Multicenter bond | electrons | charge, unpaired electrons |
| Noncovalent bond | kind | none |
| Stereo atom, stereo bond | configuration, which holds the kind and the coset | none |

Resolution fills the other fields. It never fills a principal field of a stored entity. That gives
the two sources of indeterminacy without reference to a format:

- An undetermined principal field is stated indeterminacy. A format needs a notation for it,
  because nothing restores an element, an order, or a configuration.
- An undetermined value in another field is an omission as far as resolution is concerned. Calling
  resolution asks for it to be filled.

A **wildcard entity** has an undetermined principal field or an undetermined component of that
field. For stereo, both `*` (undetermined kind) and `Th*` (tetrahedral kind, undetermined coset)
qualify. A set or variable in place of a literal expresses a different form of indeterminacy.
The wildcard predicate can be evaluated in one pass over the molecule's entities without a model.

Three facts of the present code bear on the rule:

- Resolution creates aromatic systems and stereo entities from `#a`, `#T+`, and `#C+` assertions
  and so writes their principal fields. A created stereo entity gets its kind and an undetermined
  coset. Resolution can therefore produce a wildcard entity even when the input contains none.
- No place was found that fills the charge or the unpaired electrons of a stored aromatic system,
  as `BondsResolver` and `MulticenterBondsResolver` do for their entities.
- The valence sources currently require a literal element. A set, a variable, and `*` all stop
  admission today. The wildcard policy admits only `*` in the element field.

### Outcomes

`Solution` follows constraint satisfaction and keeps its three variants. Three distinctions of
that field carry over; they are stated from memory of the literature and were not looked up.

- A model has parameters and decision variables. The parameters come with the data, and solving is
  defined for a model whose parameters are given. A principal field is a parameter.
- A solution is a total assignment. A variable may stay unassigned only when nothing else depends
  on it.
- An answer may leave variables free and then holds for every value of them. A goal that needs a
  free variable is delayed; when delayed goals remain, the answer is the unknown case with the
  residue attached.

Under `Admit`, `Determined` means the resolver completed the requested work while leaving the
admitted principal values open. For every assignment to those values that full resolution
accepts under the same model, the fully resolved result must refine the partial one. Resolution
may construct an entity and discharge an assertion that the entity entails. When a required
completion depends on choosing an admitted value, the outcome is `Underdetermined`.

| Input | What resolution can establish | Outcome under `Admit` |
| --- | --- | --- |
| `*C` | The carbon and the bond are concrete; the `*` is unchanged. | Determined |
| `C~C` | The hydrogen counts of the carbons need the order and are not filled. | Underdetermined |
| Six-membered Hückel ring with `*#a1`, scope any element | The system records one electron from the wildcard atom and discharges its `#a1`. The element stays open. | Determined |
| The same ring with a bare `*` | The missing contribution prevents resolution of the requested system. | Underdetermined |

For the bare `*` in the ring, the standalone aromaticity resolver currently reports
`AromaticValenceFailure` at a carbon. Returning `Underdetermined` for this missing contribution
is a change from the [current behavior](#aromaticity-resolver-called-alone).

### Standalone and aggregate outcomes

A standalone resolver returns `Determined` when it completes its own work. For example, take an
otherwise concrete molecule with a `#T+` assertion, an established tetrahedral frame, and a model
that permits creation of the requested entity.

| Call | Outcome |
| --- | --- |
| Standalone `StereoResolver::resolve` | `Determined`. It creates `Th*` and discharges `#T+`. The known kind entails the assertion. |
| Aggregate resolution under `Strict` | `Underdetermined`. The coset is still open. |
| Aggregate resolution under `Admit` | `Determined`. The molecule contains `Th*`. |

The same outcomes apply to a stored entity with that frame and open coset, and to the
corresponding `#C+` and `Ct*` cases. Discovery of unasserted sites belongs to `StereoAdder`.

`resolve`, `resolve_with_report`, `resolve_into`, and `resolve_into_with_report` derive their
outcome from the phases. A phase accounts for fields left open on the entities it handles,
including entities it creates. This removes the final `is_concrete()` test from all four entry
points. `has_wildcards()` and `is_concrete()` are caller queries. Consuming calls return a molecule
only on `Determined`, and in-place calls commit only on `Determined`. The report records unresolved
valence candidates and tie-breaks.

### Phase outcomes under `Admit`

- **Isotope.** A wildcard atom keeps its supplied isotope, including an undetermined one.
  Other atoms follow the isotope policy. An unresolved isotope on one of those atoms counts
  toward the aggregate outcome while later phases can still establish a contradiction.
  With `IsotopePolicy::Natural`, projection retains natural composition on a wildcard atom.
  Projecting a molecule with ground isotopes and resolving its isotopes under the same policies
  then recovers the original molecule exactly.
- **Valence.** Both the counts table and atom typing skip candidate generation for wildcard
  atoms and complete other atoms from the available evidence. The checks of `#v`, `#d`, and
  `#t` still apply to every atom. An open localized or dative order yields `Underdetermined`
  when a required completion or assertion needs the order. If a selected candidate leaves a
  required inherent field open on an atom of known element, the aggregate is underdetermined.
- **Aromaticity.** A wildcard outside the requested system can be immaterial. For a member,
  the selected rule and element scope determine which evidence is needed. Under Hückel with
  any element allowed, literal member contributions can establish a stored system whose
  electron vector is open. Its vector stays open and each member assertion that still supplies
  its contribution is retained. A missing required contribution yields `Underdetermined`, as
  does a required system charge or spin that remains unresolved. An element allow-list
  excludes a wildcard element.
- **Stereo.** A stored kindless entity can be retained when its explicit frame is established
  and the requested work leaves the kind open. If satisfying `#T+` would require choosing the
  kind of that stored entity, the outcome is `Underdetermined`. Insufficient evidence for a
  requested frame likewise gives `Underdetermined`.
- **Bond defaults and noncovalent bonds.** Localized and multicenter bonds receive the usual
  charge and spin defaults independently of an open order or electron vector. The multicenter
  phase still checks `#m`. No phase needs to choose a noncovalent bond's kind. The aggregate
  retains a bond with open kind.

At each phase, an established conflict between an assertion and an independent reading remains
contradictory. The existing failure policies apply to unrealizable assertions and entities.

### Policy

One field of `ResolveConfig`, `wildcards: WildcardPolicy`, governs wildcard admission across all
resolvers and entity kinds. `Strict` is the default and requires a concrete molecule for an
aggregate `Determined` result. A retained noncovalent bond with open kind therefore yields
`Underdetermined` under `Strict`. `Admit` permits the wildcard entities described above.

A collection of generic compounds might accept unspecified substituents and require known bond
orders even when no hydrogen count depends on them. That is a requirement on records admitted to
the collection. Resolution uses the same rule for both entity kinds.

Each resolver that reads the setting receives it at construction. `IsotopeResolver` gets `new()`
and `with_config`. If a caller deliberately alters an intermediate result or makes phase settings
inconsistent after construction, the operation must still avoid panics. Its chemistry outcome in
that case is unspecified.

Check dependencies as each phase reads the relevant fields. Compare end-to-end resolution under
both policies against the current resolver on wildcard-free inputs with both valence sources.

### Names

Each of the eight entity forms gets `is_wildcard()`. `Molecule::has_wildcards()` is true when any stored
entity satisfies `is_wildcard()`. For a molecule containing just `C#h*` it is false because the
element is carbon. For one containing `[1*]` as read from SMILES it is true because the element
is undetermined, even though the isotope is known.

Python exposes `WildcardPolicy`, the `wildcards` keyword and attribute on `ResolveConfig`, and
the same predicates. An indexed molecule entity view delegates `is_wildcard()` to its stored
form. SMILES and MOL ingestion with a supplied config publish the aggregate resolver's
determined result, including an admitted partial molecule.

The [DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md) uses *wildcard* for `*`, and
the [fragment operation](../umol-graph-ir/src/ir/molecule/fragment.rs) uses *wildcard atom* for
an atom of undetermined element. *Principal field* describes where the wildcard occurs. Use
the same term in the compact DSL where the parser currently calls that field the *head*.
*Pattern* already names the role a molecule plays in matching.

## Earlier element-only design

This historical proposal predates the [agreed phase outcomes](#standalone-and-aggregate-outcomes).

### Determined result of `Resolver`

`resolve`, `resolve_with_report`, `resolve_into`, and `resolve_into_with_report` do not test
`Molecule::is_concrete`. After the last phase they return determined. As before, the isotope plan
is applied whatever its outcome, and a later phase that is underdetermined or contradictory ends
resolution with that outcome.

The determined molecule is the input with the edits of the phases. In the six situations of the
table above it holds what the phases leave: the isotope as given, the stereo entity with its open
coset, the noncovalent bond unchanged, the atom with the fields of its selected completion, the
bond with its order undetermined, the aromatic system with its attributes undetermined. A caller
that needs a concrete molecule tests `Molecule::is_concrete`.

The ingestion and parsing entry points return such a molecule where they returned the
underdetermined error. `ResolveUnderdetermined` remains the error for an outcome that a phase
reports. `IsotopeResolver::plan_resolve` called alone still reports a plan that leaves a
non-ground isotope as underdetermined.

### Policy

```rust
/// How resolution handles an atom whose element is Undetermined.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum UndeterminedElementPolicy {
    /// The atom leaves resolution underdetermined.
    #[default]
    Strict,
    /// The atom is admitted and left unchanged; the other atoms are resolved.
    Admit,
}

pub struct ResolveConfig {
    pub isotope: IsotopePolicy,
    pub undetermined_elements: UndeterminedElementPolicy,
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}
```

The type is defined in `ops/resolve.rs` beside `ResolveConfig`. It is one setting for the whole
resolution and is absent from `AromaticityResolveConfig` and `StereoResolveConfig`. It concerns
`ElementForm::Undetermined` only. An element given as a set or a variable is treated as before.
The policy belongs to resolution; `Resolver::project` and the `project` methods of the individual
resolvers do not read it.

A resolver that reads the policy receives it at construction and stores it, as every other
setting is stored. Those are `IsotopeResolver` and `ValenceResolver` with its two candidate
sources.

### Result under `Admit`

An atom of undetermined element in the determined molecule is equal to the input atom in every
field and assertion. A bond to it receives the defaults of any other bond.

`ingest_smiles_with` and the other callers that pass a `ResolveConfig` need no change: with
`Admit` they return the molecule that resolution returns.

### Isotope

`IsotopeResolver::plan_resolve` under `Admit` passes over an atom of undetermined element: it
plans no edit for it under `Natural`, and its isotope does not make the plan underdetermined.

### Valence

`AtomTypingValence::admit` and `CountsValence::admit` under `Admit` pass over an atom of
undetermined element in both of their loops: it does not make admission underdetermined and it
gets no completions. The check of `#v`, `#d`, and `#t` against the derived readings in
`ValenceResolver::admit` still covers every atom, since it reads bonds and not elements.

### Aromaticity

No change. An atom of undetermined element that reaches selection has no completions, so its
aromatic valence is read from its own evidence, as for any atom outside valence resolution. The
outcomes are those of the table above: with a literal `#a<n>` it takes part under Hückel when the
scope is any element; under an allow-list, and under HMO and Clar, the assertion fails and
`aromatic_valence_failure` applies; a non-literal count leaves selection underdetermined.

### Stereo

For a wildcard atom, stereo frame derivation only needs to check whether the element scope
is `Any`. For a literal element it can use the existing `ElementScope::contains`.

| Stereo element scope | Literal element | Wildcard (`*`) | Set, complement, variable |
| --- | --- | --- | --- |
| `Any` | In scope | In scope | In scope |
| `AllowList` | In scope if listed | Outside scope | Unresolved |

The outcomes of the stereo table become:

| Input | Scope | Outcome |
| --- | --- | --- |
| Stored stereo atom at `*` | any element | Determined; the entity is retained |
| The same | allow-list of carbon | `StereoAtomFailure`; `stereo_atom_failure` applies |
| `#T+` on `*` with literal counts | any element | Determined; the entity is created with an open coset |
| The same | allow-list of carbon | `TetrahedralStereoFailure`; `stereo_atom_failure` applies |

The same holds for stereo bonds with an end atom of undetermined element.

The reading is part of perception and does not depend on the policy. In `Resolver` an atom of
undetermined element reaches the stereo phase only under `Admit`. `StereoResolver` called alone
and `StereoConformanceValidator`, which also calls `derive`, get the outcomes of the table under
either policy. `StereoPerceiver::recognize` keeps its own test, which requires a literal element.

### Bonds and multicenter bonds

No change.

### Constructors

`new` takes what is required and uses default settings; `with_config` takes the same and the
settings.

| Resolver | Constructors |
| --- | --- |
| `IsotopeResolver` | `new()`, `with_config(policy: IsotopePolicy, undetermined_elements: UndeterminedElementPolicy)` |
| `ValenceResolver` | `new(model)`, `with_config(model: &'a ValenceModel, undetermined_elements: UndeterminedElementPolicy)` |
| `AtomTypingValence` | `new(registry)`, `with_config(registry: &'a AtomTypeRegistry, undetermined_elements: UndeterminedElementPolicy)` |
| `CountsValence` | `new(table)`, `with_config(table: &'a ValenceTable, undetermined_elements: UndeterminedElementPolicy)` |

`IsotopeResolver::new(policy)` is removed; its callers use `new()` or `with_config`. The other
`new` constructors keep their signatures and mean `Strict`. `Resolver::with_config` builds its
isotope and valence resolvers with `with_config` from `config.isotope` and
`config.undetermined_elements`. The constructors of `AromaticityResolver`, `StereoResolver`,
`BondsResolver`, and `MulticenterBondsResolver` are unchanged.

### Python

`umol-py` exposes `ResolveConfig`, so the policy is exposed with it:

- a class `UndeterminedElementPolicy` with `Strict` and `Admit`, bound like `IsotopePolicy`;
- `ResolveConfig(*, aromaticity, stereo, isotope=IsotopePolicy.Strict,
  undetermined_elements=UndeterminedElementPolicy.Strict)`;
- a getter `ResolveConfig.undetermined_elements`;
- `repr` of a non-default config names the field after `isotope`.

### Guides

- `docs/development/operations.md`: the paragraph on `Resolver::resolve` states the determined
  result without the concreteness test; the rows of `ValenceResolver::admit` and
  `IsotopeResolver::plan_resolve` state the behavior under the policy; the row of
  `StereoResolver::plan_resolve` states the scope reading.
- `docs/development/python-api.md`: the description of `ResolveConfig` names the field.

## Staged implementation plan

[Partial resolution](#partial-resolution) is the design for this plan. The [coverage table](#coverage-of-the-earlier-plan)
maps every subitem of the element-only plan to its implementation and tests below.

Work proceeds one subitem at a time. Run the touched crates' tests, report the changed files
and a suggested commit message, then stop for review. Leave Git operations to the user during
implementation. Show the old and new assertions before editing an existing expectation.
Keep cases beside the tests that use them, in source order, following the test-writing skill
and the [property-test guide](../docs/development/property-tests.md).

### S0 — Baseline, predicates, and shared vocabulary (done)

#### S0a — Record the resolver baseline (done)

`umol-graph/benches/resolve.rs`. Additive (green). [dep: none]

Before changing resolver behavior, save end-to-end Criterion timings for the SMILES counts and
registry atom-typing models. Use the existing `smiles_roundtrip/resolve` and
`valence_readings/typing_resolve` cases. Include a one-atom input, a branched molecule, a
64-carbon chain, a fused aromatic system and stereo on both backends. Add missing cases to the
existing harness. Record the expected determined molecules outside timing and keep parsing, raising,
resolver construction and cloning outside the timed operation. Record the revision, model,
build profile, benchmark IDs and batching so S5a can repeat the comparison.

Measured the resolver at `70ed9b1adbe3beaf87f8e8a3564f32e2006c56b0` with the S0a
benchmark harness, using rustc 1.99.1, Cargo's optimized bench profile with debug info,
and Criterion 0.8.1 on an Apple M1 Pro running macOS 15.7.3. Both backends used
`MostSaturated`, natural isotope completion and journal rollback. The counts backend used
the SMILES valence table. Atom typing used the default registry. Both used the remaining
chemistry and resolver defaults.

Each cell gives the Criterion slope estimate in µs per `Resolver::resolve_into` call,
with its 95% confidence interval. The benchmark ID is the column's group prefix followed
by the case name. The prefixes are `smiles_roundtrip/resolve/` and
`valence_readings/typing_resolve/`, respectively.

| Case | SMILES counts (µs) | Registry atom typing (µs) |
| --- | ---: | ---: |
| `charged` (ammonium) | 1.046 (1.044–1.049) | 1.354 (1.350–1.358) |
| `branched` | 7.085 (7.030–7.182) | 11.636 (11.593–11.695) |
| `chain_64` | 58.116 (58.060–58.171) | 115.054 (114.395–116.320) |
| `aromatic_fused` (naphthalene) | 25.570 (25.510–25.640) | 30.517 (30.330–30.860) |
| `tetra_four` | 4.421 (4.408–4.435) | 5.613 (5.587–5.643) |
| `alkene_four` | 5.069 (5.051–5.089) | 6.429 (6.411–6.448) |

Both models returned `Determined` and matched the six expected molecules in
[`resolve.rs`](../umol-graph/benches/resolve.rs), including the aromatic system and stereo
frames. Parsing, raising, resolver construction and these assertions run before timing.
`iter_batched_ref` with `BatchSize::SmallInput` clones each input before the timer starts
and drops the mutated molecules and returned reports afterward.

The raw samples and estimates are saved locally under the Criterion baseline name `s0a`.
From the workspace root, S5a can compare the same IDs with this command for each policy.
The original run used `--save-baseline s0a` in place of `--baseline s0a`.

```sh
CRITERION_HOME="$PWD/scratch/253-resolver-baseline/criterion" \
  cargo bench --locked --offline -p umol-graph --bench resolve -- \
  '^(smiles_roundtrip/resolve|valence_readings/typing_resolve)/(charged|branched|chain_64|aromatic_fused|tetra_four|alkene_four)$' \
  --baseline s0a --warm-up-time 3 --measurement-time 5 --sample-size 100 --noplot
```

#### S0b — Wildcard predicates on the eight entity forms (done)

`umol-graph-ir/src/ir/{atom,bond,dative,aromatic,multicenter,noncovalent,stereo}.rs`.
Additive (green). [dep: none]

Add `is_wildcard()` to each form in the [principal-fields table](#principal-fields). Use the
stored principal field. Both a kindless stereo configuration and a known kind with an open
coset qualify. Test each form with an open and a literal principal value, and all applicable
set, complement, variable and stereo-expression forms. Vary other inherent fields and
constraints independently. Groundness, normalization and a chemistry model are unnecessary.

All eight predicates inspect the stored principal field in constant time. The atom and bond
tests include an empty element complement and an empty conjunction over bond order. Both
return false for `is_wildcard()` even though normalization produces `Undetermined` in each case.

#### S0c — Molecule wildcard query (done)

`umol-graph-ir/src/ir/molecule.rs`. Additive (green). [dep: S0b]

Add `Molecule::has_wildcards()` as a short-circuit query over the eight entity stores. Test the
empty molecule, a molecule with each kind as its sole wildcard, `C#h*`, and a wildcard atom
with an explicit isotope as in `[1*]`. Keep `Molecule::is_concrete()` available to callers.

#### S0d — Shared admission policy (done)

`umol-graph/src/ops/resolve.rs`. Additive (green). [dep: none]

Define `WildcardPolicy::{Strict, Admit}` beside `ResolveConfig` with the value traits of
`IsotopePolicy`. Add `test_wildcard_policy_default` for `Strict`. S3b exposes the config field
together with working aggregate admission.

#### S0f — Record the nomenclature and use it in the compact DSL (done)

`docs/development/nomenclature.md`, `umol-graph-ir/spec/umol-dsl-spec.md`,
`umol-graph-ir/src/dsl/{electrons,stereo,aromatic,multicenter,edit}.rs` and
`umol-py/tests/test_stereo.py`. Terminology change (green). [dep: S0c]

Define *principal field* and *wildcard entity* alongside concreteness in the nomenclature guide
using the [eight entity forms](#principal-fields). Record the names `is_wildcard()` and `has_wildcards()`.
Use `C#h*` and `Th*` to illustrate which open fields make the entity a wildcard.

Replace *head* with *principal field* wherever it denotes that field in the compact DSL.
Update the documentation of `electrons::electron_counts`, `stereo::stereo_configuration`,
`stereo::fmt_stereo_configuration` and the two stereo update types, along with the stereo-string
specification. Rename the aromatic and multicenter parser case labels from `missing_head` to
`missing_principal` and the stereo configuration test parameters from `head` to `principal`
in the Rust edit and Python stereo tests.

Run the existing electron-count, aromatic, multicenter and stereo parsing, formatting and
edit tests and the Python stereo tests. Preserve their accepted strings, exact errors and
renderings, including the distinction between an omitted configuration and `*` in stereo updates.

### S1 — Isotope and valence resolution (done)

#### S1a — Isotope completion and projection (done)

`umol-graph/src/ops/{resolve/isotope.rs,resolve.rs}` and `umol-graph/tests/property/isotope.rs`.
Breaking (red→green within the subitem). [dep: S0b, S0d]

Store the wildcard policy. Replace `IsotopeResolver::new(policy)` with `new()` (both policies
strict) and `with_config(policy: IsotopePolicy, wildcards: WildcardPolicy)`. Migrate every
constructor call, including the aggregate, unit tests and properties, in the same change.
Under `Admit`, skip a wildcard atom during completion. In `plan_project` retain its natural
isotope under `Natural + Admit`. Preserve the ground-isotope precondition for projection.
Exercise the plan, consuming and in-place routes.

Test both constructors. In `test_isotope_resolver_plan_resolve` cover mixed `*` and carbon with
omitted isotopes, `*` alone, and supplied mass 1.
Under `Natural`, strict wildcard handling plans both defaults on the mixed input, and admission
plans only the carbon default. A lone wildcard under `Natural + Admit` gives a determined empty
plan. Cover omitted, set and variable isotopes on known elements.
Extend the exact project/resolve roundtrip and resolve idempotence properties across both
policies, with literal and wildcard atoms, natural composition and explicit masses.

#### S1b — Registry admission (done)

`umol-graph/src/ops/valence/atom_typing.rs`. Additive (green). [dep: S0b, S0d]

Add `AtomTypingValence::with_config(registry, wildcards)` and keep `new(registry)` strict.
Pass over an admitted wildcard element in the element check and candidate loop. At a known
atom, distinguish complete fields from missing work that needs an open localized or dative
order. Preserve established mismatches and the existing tie-break inputs. Check nonliteral
elements before extracting a literal.

Test construction and retain the strict `test_atom_typing_valence_admit_partial` cases. Extend
`test_atom_typing_valence_admit` under `Admit`. Here `*C(O)C` has
candidates for the three known atoms and none for `*`. A lone `*` gives determined empty
completions. Element sets and variables still block admission. Pair an open order needed to
complete a known atom with a case whose known fields and assertions need no choice of order.

#### S1c — Counts admission (done)

`umol-graph/src/ops/valence/counts.rs`. Additive (green). [dep: S0b, S0d]

Add `CountsValence::with_config(table, wildcards)` and keep `new(table)` strict. Apply the S1b
rules to counts enumeration. In particular, distinguish skipping an already complete atom
from skipping one whose missing hydrogen count needs the order. Test construction and keep
`test_counts_valence_admit_partial` for strict handling. In `test_counts_valence_admit` check
the exact candidates or outcome for `*C(O)C`, a lone wildcard, a symbolic element and the
two open-order cases.

#### S1d — Valence dispatch and incidence assertions (done)

`umol-graph/src/ops/resolve/valence.rs`. Additive (green). [dep: S1b, S1c]

Add `ValenceResolver::with_config(model: &ValenceModel, wildcards: WildcardPolicy)`, forwarding
to the selected source. `new(model)` selects `Strict`. Keep `#v`, `#d` and `#t` checks on every
atom, including a skipped wildcard. Report underdetermination for an asserted reading that
needs an open order and preserve an independently established contradiction. Test construction
through both dispatch branches. Extend `test_valence_resolver_admit` with `*C(O)C` under both
sources and entailed, blocked and contradictory assertions.

### S2 — Aromatic, stereo and bond phase outcomes

#### S2a — Aromatic selection and retained systems (done)

`umol-graph/src/ops/{aromaticity.rs,aromaticity/,resolve/aromaticity.rs}`.
Behavior change (green). [dep: S0b]

Apply the same evidence rules in standalone planning and joint selection, including the
reporting path. Under Hückel with any element allowed, a required ring member with `*#a1`
can supply one electron. A bare wildcard that would have to supply a missing contribution
leaves the requested system underdetermined. An exocyclic wildcard can be immaterial. Retain
the allow-list exclusion, the HMO and Clar model requirements, and the configured treatment of
an independently established failure.

For a newly established system, record its electron contributions and discharge the `#a1`
assertions it entails. For a stored system with independent member evidence, retain its
principal electron vector, including an open vector, and keep assertions that still supply a
contribution. Account for any required system charge or spin still open. Test new and stored
six-membered systems, retained member assertions, `*#a+`, a bare ring member versus an exocyclic
`*`, and permissive versus Daylight scope. Retain the HMO odd-electron and Clar nonbenzenoid
controls in the standalone table. Cover insufficient evidence under both `Error` and `Keep`
and exercise both atom-selection algorithms.

#### S2b — Stereo scope and frame derivation (done)

`umol-graph/src/ops/{stereo.rs,validate/stereo.rs}`. Behavior change (green).
[dep: none]

Apply the [stereo element scope table](#stereo) directly in `derive_stereo_atom`,
`derive_stereo_bond` and `identify_virtual_ligands`. Keep the literal-element requirements
of independent site recognition. Cover a stored stereo atom at `*` and a requested `#T+`
there under both scopes, a stereo bond with a wildcard endpoint, and a carbon site with a
wildcard ligand. Test literals inside and outside an allow-list and sets, complements and
variables under both scopes. Extend `test_stereo_perceiver_derive`, `test_stereo_perceiver_derive_stereo_atom`
and `test_stereo_perceiver_derive_stereo_bond`, and exercise the validator on a stored site.

#### S2c — Standalone stereo resolution and the configuration left open (done)

`umol-graph/src/ops/{resolve/stereo.rs,stereo.rs}`. Behavior change (green).
[dep: S0b, S2b]

Retain or create `Th*` or `Ct*` when a permitted request establishes its frame and kind.
Standalone resolution is determined in that case. Handle a stored kindless entity with an
established explicit frame separately from a request that requires choosing its kind.
Carry the configuration of retained and newly created entities into the aggregate decision
in S3a. A removed entity contributes no open field to the result.

In `test_stereo_resolver_plan_resolve` cover creation, retention and discharge, unique virtual-ligand identification,
incomplete frames, a known kind with a coset set, and a kindless row with and without a positive
kind assertion. Only an undetermined coset or kind qualifies for wildcard admission. Exercise
scope failures with the default policy and `Keep`, as well as existing `Remove` behavior.
Keep automatic discovery with `StereoAdder`.

#### S2d — Bond defaults and the remaining bond kinds

`umol-graph/src/ops/resolve/{bonds,multicenter}.rs`. Additive phase
accounting (green). [dep: S0b, S0d]

Account for the resulting localized-bond order, multicenter electron vector, charge and spin
while planning defaults. Preserve the scoped standalone return contracts, including the
`Edits` returned by `BondsResolver::plan_resolve`. Keep `#m` checks and distinguish a blocked
reading from a conflict. Include retained dative orders and noncovalent kinds in aggregate
bond accounting. Inspect those stores directly when populated. Distinguish an open dative
order, which `Admit` permits when unused, from a surviving symbolic order, which still
blocks publication. Valence accounts for uses of the dative order.

Test a bond with an open principal field and independently completed charge/spin, partially
specified spin, a residual non-wildcard field, and a `#m` assertion with complete or insufficient
evidence. Test the dative and noncovalent field readings here. S3b exercises both policies on
otherwise complete molecules with those bonds.

### S3 — Aggregate resolution and its callers

#### S3a — Derive the strict aggregate outcome from the phases

`umol-graph/src/ops/resolve.rs` and its phase planners. Internal rewire (green).
[dep: S1a, S1d, S2a, S2c, S2d]

In all four `resolve` entry points, account for the result left by each phase while that phase
has the entities in hand. Carry the isotope outcome alongside its applicable edits so a later
contradiction can still surface. At constitution, inspect the selected atom fields and account
for known atoms that received no candidate. At aromaticity, stereo and bond defaults, use the
fields the plan retains or produces, including open configurations created from assertions.
A field that only blocks publication must allow later phases to establish contradictions and
retain their report information. Compose these findings internally with the edits each phase
already plans. Preserve the scoped standalone results and the public plan and report interfaces.

Once every entity kind is accounted for, remove the final `Molecule::is_concrete()` calls
from `resolve`, `resolve_with_report`, `resolve_into` and `resolve_into_with_report` together.
Resolve never calls `has_wildcards()`. Preserve integrity checks, consuming ownership,
report contents, contradiction/error propagation and in-place rollback under `Journal` and
`Clone`. Update the module and four method rustdocs and the `ResolveUnderdetermined` rustdoc
to describe the policy-dependent result. Keep its error message and boundary mappings.

Use the six cases in [Final test of Resolver](#final-test-of-resolver) as strict regressions,
plus dative and multicenter bonds. Exercise consumed and borrowed molecules with and without a
report. Assert standalone stereo success and aggregate underdetermination for the same `Th*`.
Check later contradictions despite an unresolved isotope, and exclude an open stereo row removed
by a failure policy from the final decision.

#### S3b — Expose the policy through the aggregate config

`umol-graph/src/ops/resolve.rs`, Rust `ResolveConfig` literals across the workspace, and
`umol-py/src/resolve.rs`. Breaking Rust struct-literal change (red→green within the subitem).
[dep: S0d, S3a]

Add `wildcards: WildcardPolicy` after `isotope` on `ResolveConfig`. `Resolver::with_config`
passes the chosen policy to the isotope and valence constructors and uses it for the aggregate
phase outcomes. `Resolver::new` and `ResolveConfig::default` select `Strict`. Migrate exhaustive
literals in the same change, including the aggregate constructor tests, molecule and reaction
ingestion tests, and the Python constructor and its four Rust test expectations. The existing
Python signature supplies the strict default until S4b. Extend `test_resolve_config_default`
and `test_resolver_with_config` to verify the default and propagation into the isotope resolver
and both valence sources.

Run the following cases through the four entry points. Check exact molecules on success,
reports where requested, and restored receivers on rejection.

- `*C(O)C` under `Admit` keeps the wildcard and gives the central C, O and methyl C one, one
  and three hydrogens respectively. Its bond to `*` has zero charge and closed-shell spin.
  The wildcard contributes no report entry. `Strict` is underdetermined.
- `[1*]C` retains mass 1 and zero implicit H on `*`. Cover a lone `*`, two wildcard substituents
  on one carbon, a carbon stereo site with a wildcard ligand, and an element set or variable.
- `C~C` with missing hydrogen counts is underdetermined under both policies. Pair it with a
  fully specified example where no required reading depends on the open order. Exercise the
  corresponding dative, multicenter and noncovalent cases.
- Resolve the six-membered `*#a1` ring under `Admit` with permissive Hückel and Daylight
  models. Expect the six-atom system under the former and a contradiction under the latter's
  default failure policy. Exercise the S2c stereo cases with concrete surrounding atoms.
  A strict aggregate rejects a surviving open coset.
  `Admit` publishes it when the required work succeeds. Cover an input that gains its first
  wildcard during resolution.
- Check access to caller-altered intermediate results at the point of use. Exercise these
  consumers and inconsistent phase settings on the aggregate. Assert absence of panics. Normal
  configurations retain the full outcome, report and rollback assertions.

#### S3c — SMILES and MOL boundaries

`umol-graph/src/{ingest,parse,export}.rs` and `umol-graph/tests/property/project.rs`.
Integration coverage (green). [dep: S3b]

Follow `ResolveConfig` through molecule and reaction SMILES interpretation and configured MOL
parsing, including the byte entry points. Publish the determined partial molecule, retain
side-specific errors on reaction underdetermination, and keep unresolved ingestion available.
Exercise `test_ingest_smiles_with_resolution` on `*C(O)C` and `[1*]C` under `Natural + Admit`.
In `test_ingest_smiles_with_error` keep the custom carbon registry whose selected entry
leaves required fields open and add the default-model wildcard case. Keep the `*` case of
`test_ingest_smiles_error`.

In `test_export_smiles_with_mol` preserve the whole generated MOL domain, including both
valence sources and the existing render/reingest checks for definite geometries. For the
double bond marked “either”, test strict rejection with restoration and admitted resolution
with the stored open cis-trans coset. Check export against the format's own representability
rules. A successful resolve alone establishes no format roundtrip.

#### S3d — Cross-phase properties

`umol-graph/tests/property/{resolve,isotope,project}.rs`. Additive coverage (green).
[dep: S3b, S3c]

Retain the generators and assertions of `test_resolver_resolve_carbon` and
`test_resolver_resolve_into_heteroatoms`. These inputs have literal principal fields and no
requests that create wildcards. Require equal outcomes and determined molecules under both
policies. On the terminal-wildcard variants require the wildcard atom to be preserved and the
other atoms to be concrete on an admitted success. Those generators have no assertion requiring
discharge on the wildcard.

Check the refinement law by resolving a bounded family of literal assignments to admitted
elements or orders under the same model. For a determined partial result, every accepted
full result refines its retained fields and relations, allowing entailed assertions to be
discharged. Include inputs with accepted assignments and compare fields in corresponding
participant frames. Retain the strict implication that
every determined aggregate result is concrete. `is_concrete()` is a suitable test assertion.
Check parity of reports, outcomes and recovered input across the four aggregate routes.

### S4 — Python surfaces

#### S4a — Form, view and molecule queries

`umol-py/src/{atom,bond,dative,aromatic,multicenter,noncovalent,stereo,molecule}.rs`
and their Python tests. Additive (green). [dep: S0c]

Bind `is_wildcard()` on all eight owned entity forms and on the attached entity views returned
by indexing a molecule. Delegate to the underlying form, using the view's existing access and
invalidation checks. Bind `Molecule.has_wildcards()`. Test each kind from an owned form and
an indexed view, a changed principal field read through a live view, and ordinary invalidation.
Include `C#h*` and the isotope-labeled wildcard in molecule tests.

#### S4b — Policy, configuration and boundary outcomes

`umol-py/src/{resolve,lib,molecule,reaction}.rs`, `umol-py/python/umol/__init__.py`
and `umol-py/tests/`. Additive (green). [dep: S3b, S3c]

Bind and export `WildcardPolicy.Strict` and `WildcardPolicy.Admit` with `from_rust` and
`to_rust`. Add `wildcards=WildcardPolicy.Strict` to the keyword-only `ResolveConfig` constructor
and its read-only getter. Preserve `ResolveConfig.default()` for the default repr and place
`wildcards` after `isotope` in the non-default repr.

Exercise the Rust/Python conversions, default equality and repr in the binding and config tests.
Update `test_import.py` (including `__all__`) and `test_class_constructor_signatures.py`.
In `test_molecule.py` compare
`Molecule.from_smiles` under `Admit` with the exact `*C(O)C` result and check the default error.
Exercise partial success, underdetermination, reports, ownership and view invalidation through
the four resolution methods. Keep the lenient `C~C>>C.C` reaction case under both policies
and include a reaction with an admissible wildcard. The Python wrappers forward the config
and use the existing Rust outcomes.

### S5 — Verification and permanent documentation

#### S5a — Measure both policies on ordinary molecules

`umol-graph/benches/resolve.rs`. Additive benchmark coverage (green). [dep: S0a, S3b]

Run the S0a cases against the saved baseline for both `Strict` and `Admit` on both valence
sources. Hold chemistry, tie-breaks, rollback, batching and timed work constant. Check that
each wildcard-free fixture still produces the baseline molecule. Record end-to-end changes
per case and repeat noisy or regressing cases before attributing a cost to the phase checks.
Investigate and resolve a significant common-case regression before accepting this stage.
The earlier stand-alone scan timing cannot measure the implemented phases.

#### S5b — Publish the contracts

`docs/development/{operations,python-api}.md` and the affected rustdoc.
Documentation (green). [dep: S0f, S3d, S4a, S4b, S5a]

Update the phase, discharge, isotope projection and aggregate outcome descriptions, plus the
`ResolveConfig` and predicate entries in the Python guide. Cross-check examples against the
tested `*C`, `C~C` and `Th*` cases. Verify the guide and public rustdoc agree about caller
queries, the strict default, standalone outcomes, admitted results and rollback.

### Coverage of the earlier plan

| Earlier subitem | Updated subitems | Work carried forward |
| --- | --- | --- |
| S0a · final concreteness test | S2a–S2d, S3a–S3d, S4b | Four aggregate exits, rustdoc and boundary error, six strict regression cases, Rust and Python expectations below. Removal follows phase accounting. |
| S1a · policy and config | S0d, S3b, S4b | Enum, strict default, field, exhaustive Rust literals and Python constructor migration use the agreed `WildcardPolicy` name. |
| S1b · element scope | S2b | Stereo frame derivation tests scope membership, exclusion and unresolved cases directly. Hückel already admits a wildcard only under `Any`. |
| S2a · isotope resolver | S1a, S3a, S3b | Constructor migration and wildcard defaulting, plus projection and exact restoration. |
| S2b · valence sources | S1b, S1c | Both constructors, prior strict and new admitted candidate cases, and resolution of order-dependent work. |
| S2c · valence dispatch | S1d, S3b | Both sources receive the policy and incidence checks cover skipped atoms. |
| S2d · stereo scope | S2b, S2c, S3a, S3b | Perceiver, virtual ligands, both scopes, failure policies, validator, and standalone versus aggregate outcomes. |
| S3a · aggregate construction | S3b–S3d | Constructor propagation, exact integration fixtures, four routes, ingress and both original property domains. |
| S4a · Python | S4a, S4b | Bindings, exports, signature, repr, config and ingestion tests, plus predicates on forms, live views and molecules. |
| S5a · guides | S0f, S5b | S0f records the terms and aligns the compact DSL. S5b covers the operations and Python guides. |

The following expectations replace the unconditional partial successes in the earlier S0a.
Keep the original fixtures and controls when adding the new policy dimension.

| Existing tests | Expected coverage in the updated plan |
| --- | --- |
| `test_resolver_resolve_isotope`, `test_resolver_resolve_isotope_partial`, `test_resolver_resolve_with_report`, property `test_resolver_resolve_into_isotope` | On a known element, an omitted isotope under `IsotopePolicy::Strict` or a surviving set or variable yields `Underdetermined` under either wildcard policy. `Natural` still fills an omitted isotope. Retain the report, restoration and idempotence assertions. |
| `test_resolver_resolve_later_underdetermined`, `test_resolver_resolve_stereo_assertions` | Keep all strict cases. Under `Admit`, compare the three requested stereo rows to the exact created entity and discharged assertion, the stored open row to its retained entity, and the stored virtual row to the identified ligands. |
| `test_resolver_resolve_underdetermined` | Keep strict rollback for the open noncovalent kind. Under `Admit`, check the completed first atom and retained bond. |
| `test_ingest_smiles_with_error`, property `test_export_smiles_with_mol` | Retain the incomplete-registry and default wildcard errors and split the MOL “either” outcome by policy (S3c). |
| Python `test_molecule_resolve_isotope` and `test_molecule_resolve_with_report_empty` | A known carbon with omitted isotope under strict isotope handling remains underdetermined, including the empty report. |
| Python `test_molecule_resolve_into_aromaticity_partial` | Retain strict rollback with `Error` and `Keep`. Under `Admit` check the stored open vector and assertions that still supply contributions when the other required fields are complete. Also exercise missing required evidence. |
| Python `test_reaction_from_reaction_smiles_io_config` | Both wildcard policies report underdetermined reactants for `C~C>>C.C` when H completion needs the order. |

### Stage gates and critical path

Format with `cargo +nightly fmt --all`. After each stage run
`cargo test -p umol-graph --all-features` (includes the `resolution` and `property` targets)
and `cargo clippy -p umol-graph --all-features --all-targets -- -D warnings`. At the IR
stage also run `cargo test -p umol-graph-ir --features proptest`. For each stage touching
the bindings or their tests, run `cargo test -p umol-py --lib`, rebuild with
`maturin develop --locked --manifest-path umol-py/Cargo.toml` and run `pytest -q umol-py/tests`.
Activate `umol-py/.venv` in the same shell for each command that compiles or loads the
extension, including workspace commands, and confirm Python 3.13 first. Finish with the
workspace test and Clippy gates from `AGENTS.md`.

The valence path is S0b/S0d → S1b/S1c → S1d → S3a. The stereo path S2b → S2c also
joins S3a, together with S1a, S2a and S2d. From S3a the path runs through S3b and S3c to
S3d and S4b, then S5b. S0f and S4a can start after S0c, and the S0a baseline and S3b enable S5a.
All stages are required.
