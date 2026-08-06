# Development guides

These are living, normative guides for repository development. They are updated in place when a
rule changes and carry neither dates nor implementation status.

- [Data type contracts](data-types.md) defines construction, conversion, validation,
  transformation, provenance, fallibility, and public contract review.
- [Nomenclature](nomenclature.md) defines repository-specific terms and public naming conventions.
- [Property tests](property-tests.md) defines executable-specification, evidence, documentation,
  and property-suite organization policy.

Dated documents under `discussion/` retain design rationale, alternatives, implementation plans,
and historical snapshots. They may record the provenance of a settled rule but are not normative.
Source comments and public rustdoc must be self-contained and must not require a discussion document
to explain their contract.

When a discussion settles a general development rule, update the applicable living guide in the
same work item. Do not copy the rule into a second guide. `discussion/000-status.md` tracks dated
work; living guides do not appear in that status table.
