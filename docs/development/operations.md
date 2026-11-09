# Operation semantics

This guide defines shared contracts for molecular operations. Construction and validation tiers
are defined in [Data type contracts](data-types.md); representation prerequisites remain in
[Representation integrity](integrity.md).

Document intended contracts separately from observed implementation behavior. Missing documentation
does not by itself imply an open design choice. Check implementation against the intended contract
and record violations as implementation defects.

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
