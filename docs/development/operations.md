# Operation semantics

This guide defines shared contracts for molecular operations. Construction and validation tiers
are defined in [Data type contracts](data-types.md); representation prerequisites remain in
[Representation integrity](integrity.md).

Document intended contracts separately from observed implementation behavior. Missing documentation
does not by itself imply an open design choice. Check implementation against the intended contract
and record violations as implementation defects.

## Resolution

Valence, aromaticity, and stereo resolution use the same conflict and discharge contract.
Constraints restrict admissible results. Resolution may refine unresolved attributes and construct
compatible entities under the selected chemistry model. It must not overwrite explicit data or
delete a constraint to hide a disagreement. Model tie-breaks select among admissible candidates.

A determined resolution means that its plan can be applied under the configured policies. It does
not establish every validator's predicate or require every stored form to be ground. Stereo
configurations may retain undetermined cosets, and Keep may retain an unrealizable description.

### Refinement and conflicts

Refining an attribute narrows its admitted values while respecting its stored information and
constraints. A compatible constraint may still add information to the result and must then remain
stored. A definite conflict between independently meaningful descriptions is contradictory under
every failure policy. Neither deleting the constraint nor replacing the entity with a
constraint-derived candidate is permitted as mismatch recovery.

Perception and derivation classify findings without applying recovery policy. A constraint failure
means that the requested description is unrealizable under the selected topology and model; an
entity failure means that a stored description is unrealizable. A mismatch means that independently
meaningful descriptions disagree. Ambiguity is underdetermination, not a realizability failure.
An established conflict must not be suppressed by uncertainty elsewhere.

### Discharge

Discharge removes a constraint that adds no information to the planned result. Let C be its
normalized asserted value and D its normalized, independent structural reading from that result:

| Meet | Result |
| --- | --- |
| C ∧ D is bottom | Conflict; reject the resolution. |
| C ∧ D = D | Discharge C; the structure independently entails it. |
| C ∧ D is compatible but differs from D | Retain C; it still restricts the result. |

Use the open derived reading for D. In the keyed-reading API, None is bottom and
Some(undetermined(key)) is top. The evidenced readings include the assertion being tested and
therefore cannot establish its redundancy. A complete reading that closes an absent relation to a
negative value cannot justify discharge either. Candidate admission may use complete evidence
according to its own contract; discharge still requires independent open derivation.

The same rule applies to newly constructed and retained entities. Creation history does not
determine whether a constraint is redundant.

| Constraint | Independent structural reading | Result |
| --- | --- | --- |
| #v4 | Valence 4 | Discharge. |
| #v4 | Valence 3 | Conflict. |
| Atom #a1 | Aromatic-system contribution 1 | Discharge. |
| Atom #a1 | Aromatic-system contribution 2 | Conflict. |
| #T+ or #C+ | Entity of the asserted kind, with any coset | Discharge. |
| #T! or #C! | Entity of the excluded kind | Conflict. |
| #T+ or #C+ | No entity | Construct a compatible entity or retain the constraint; absence alone is not a conflict. |
| #T! or #C! | No entity | Retain; the open reading is undetermined. |
| #d0 or #t0 | No dative incidence | Retain; the open reading is undetermined. |

### Realizability failures and recovery

Recovery policies apply only to the classified failures they name. They do not choose a winner
between conflicting descriptions. Current resolver recovery actions are:

| Action | Meaning |
| --- | --- |
| Error | Report the classified failure as Contradictory, in the semantic Solution channel. |
| Keep | Retain the failed constraint or entity. |
| Remove | Remove the failed constraint or entity identified by the diagnostic, using ordinary reference compaction. |

AromaticityFailurePolicy provides Error and Keep for aromatic_valence_failure and
aromatic_system_failure in AromaticityResolveConfig. StereoFailurePolicy provides Error, Keep,
and Remove for stereo_atom_failure and stereo_bond_failure in StereoResolveConfig; each stereo
field covers both constraints and entities of its kind. These policies do not turn ambiguity into
a realizability failure. Operational errors remain separate from semantic contradictions.

### Resolver scope

| Domain | Owned constraints and planning |
| --- | --- |
| Valence | ValenceResolver::admit restricts atom completions. Resolver::plan_constitution applies the selected completions and discharges independently entailed #v, #d, and #t on all atoms, including already concrete atoms. |
| Aromaticity | AromaticityResolver discharges entailed atom aromatic-valence and bond-aromaticity assertions beside new and retained systems. Standalone and aggregate resolution use the same rule. |
| Stereo | StereoResolver discharges entailed atom #T and bond #C constraints beside new and retained entities. Kind agreement does not depend on the coset; stored configurations are preserved. |

Each resolver handles its owned keys. This is not a general constraint-validation or closing pass:
other entity assertions and molecule-level logical constraints remain stored for explicit
validation. A transformer removing an entity follows its own assertion-removal contract; that
removal is distinct from resolver discharge.

Planning, ownership, rejection, and publication guarantees are defined in
[Data type contracts](data-types.md#resolution-projection-and-transformations).

## Tier-3 conformance validators

A conformance validator checks molecular descriptions against a selected chemistry model.
The model determines chemical acceptance; operation configuration selects algorithms and execution
limits. Different models may legitimately give different conformance results for the same input.

Conformance validators operate on integrity-valid representations without modifying them. They may
assume tier-2 properties only when their public contract or producer establishes those properties.
Before a model-specific query, establish the model-independent prerequisites that query requires.
Report a failed prerequisite as such, rather than attributing it to the selected chemistry model.

### Outcomes

Validation uses Result<Solution<(), C>, E>, with the operation's contradiction and error types:

| Outcome | Meaning |
| --- | --- |
| Determined(()) | Every check within this validator's scope is established. This does not establish other validators' predicates or require every stored form to be ground. |
| Underdetermined(()) | Available information leaves at least one required check undecided, with no established contradiction. |
| Contradictory(C) | An applicable check establishes a semantic violation. |
| Err(E) | Setup or execution failed outside the semantic verdict, for example because required model parameters are unavailable. |

An established contradiction takes precedence over uncertainty elsewhere. An unresolved value
must not suppress an independently established rejection. Ambiguity is underdetermination,
not evidence of impossibility.

Validators may share read-only perception and derivation with resolvers. They report established
constraint failures, entity failures, and mismatches without applying resolver recovery policies.
Validation does not select an authoritative description, repair or remove an entity, materialize a
candidate, or discharge a constraint.

### Constraints and explicit entities

A constraint restricts admissible descriptions. Model-dependent validation can check whether a
positive assertion is realizable without storing the entity derived for that check. A successful
realizability check does not establish that the entity is already present in the supplied molecule.

Agreement with stored structure uses independent derived information. A stored entity has its own
conformance obligations, whether or not an additional scalar assertion accompanies it. For a
supported known kind, an entity's unrealizability must therefore be reported even when its site has
no corresponding scalar assertion. A definite disagreement between independently meaningful
descriptions is contradictory.

Model-independent constraint satisfaction remains a tier-2 check. Its existing open/complete
reading choice is separate from conformance and from the open reading required for resolver
discharge. For example, the incidence check of #T+ without a stereo atom is underdetermined with
Derived and contradictory with DerivedComplete. Model-dependent perception may separately establish
whether that requested description is realizable.

### Domain checks

| Domain | Conformance checks |
| --- | --- |
| Connectivity | Topology accepted by the selected connectivity model. |
| Valence | Atom states accepted by the selected atom-typing or counts model. |
| Aromaticity | Realizability of asserted aromatic contributions and stored aromatic systems, and their agreement under the selected model. |
| Stereo | Realizability of asserted stereo kinds and explicit stereo descriptions, model ring-size limits, and asserted stereogenicity, topicity, and ligand symmetry. |

Stereo scalar constraints #T/#C describe kind, independently of the entity's coset. An open coset
does not prevent deciding kind agreement, although configuration-dependent checks may remain
underdetermined. Neither a positive stereo constraint nor an explicit stereo entity asserts
stereogenicity.

Frame and configuration prerequisites are checked directly on explicit entities before queries
that require them. This includes agreement with the declared kind and admissible coset indices;
normalization alone does not validate coset ranges. These model-independent prerequisites remain
distinct from model acceptance. Virtual-ligand identification may be used for checking without
changing the stored frame; unresolved identification contributes underdetermination.

The stereo model supplies both ring-size limits. The general limit also applies to positive bond
stereo assertions. The E-specific check uses the configuration relative to the ring's continuation
ligands; a raw coset index without its frame cannot establish ring-relative E. Validation reports
violations without altering supplied descriptions.
