---
name: python-api-review
description: Adversarial review of umol's Python API against its owning Rust operations. Use for Python binding audits, semantic parity reviews, and reviews of ownership, lifecycle, snapshots, cloning, or Python-only API growth. Produces evidence and dispositions without implementing fixes. Separate from the Rust review-cycle skill; not required for ordinary binding edits.
---

# Python API review

Review the requested Python operation family against its Rust owners and consumers. Preserve sound
APIs; propose removal when an exposed operation is wrong, confusing, redundant, unnecessarily costly,
or forces an unjustified lifecycle. Do not expand the task into a general Rust hygiene review.

## Authority and scope

- Read AGENTS.md, `docs/development/python-api.md`, and the relevant data-type and nomenclature
  contracts. Use the status index to locate the owning discussion and distinguish current behavior
  from proposed changes. Apply data-type-contracts for invariant-bearing surfaces and test-writing
  when reviewing or constructing test cases.
- Every Python/Rust difference needs a specific justification: types, methods, names, construction,
  reference semantics, lifecycle, ownership, aliasing, errors, and copying. Existing code, tests,
  documentation, and previous acceptance are not exemptions. An omitted Rust feature is not
  automatically a defect or a request to expose it.
- Review a named family at a time. For a broad request, inventory the surface and divide it into
  bounded passes, beginning with reference-bearing operations, consuming operations, snapshots,
  and nested mutation. Report the boundary of each pass; do not imply whole-API coverage.
- Record the source revision and relevant working-tree changes, including governing policies.
  Give all reviewers the same fixed source and policy snapshot. Do not require a commit merely to
  review uncommitted policy or code. Keep instrumented reproductions in an isolated temporary copy.
  Review does not authorize production-code fixes, git mutations, or unrelated refactors.

## Inventory and tracing

Find the public surface through module registration, Python exports, generated bindings, and types
returned by operations even if they are not exported by name. For each type and operation, identify
the owning public Rust symbol and follow conversions into and out of it. Include properties,
operators, iterators, constructors, and lifecycle transitions.

Record the difference, its claimed justification, and where it is implemented. Methods generated
by one implementation may share evidence only when their contracts and relevant branches match;
list the covered members and check exceptions. A spot check does not clear unexamined members.

Trace enough of the producer and consumer to answer:

| Area | Questions |
| --- | --- |
| Operation ownership | Does Python delegate the Rust operation, or invent composition, normalization, mutation, or validation? Does familiar syntax change reference meaning? |
| Lifecycle and failure | Are consuming, borrowing, reusable, lazy, and transactional operations preserved? What can aliases observe after success, failure, abandonment, or a second call? |
| References and construction | Which namespace owns each id or handle? Are independently assembled objects combined correctly? Which boundary first requires contextual checks? |
| Nested access | Is the result owned, shared, or live? Do writes affect the expected owner? Are immutability, equality, and hashing consistent with aliasing? |
| Copying and execution | What is copied, when, and at whose request? Does a snapshot change visibility, eagerness, lifetime, or recovery semantics? Does the wrapper discard a Rust optimization? |
| Names and surface | Does a name preserve the Rust concept? Is the difference useful ergonomics, a second operation, or an unnecessary public type? |

## Evidence

- Reproduce semantic findings with a small input and observable consequence. Exercise the consumer:
  equal assembled lists do not establish equal application behavior. Distinguish a source-level
  trace from an executed reproduction and retain exact commands and results for the latter.
- Trace copies through conversions and subsequent mutations. Distinguish shared ownership, copying
  owned data, and later copy-on-write; counting clone calls is insufficient. Consider borrowing,
  ownership transfer with explicit consumption, or supported live access before accepting a copy.
  Do not replace justified copies with speculative wrapper machinery.
- Verify claims that PyO3 or the Python ABI requires a design against the relevant version and code
  path, using dependency sources or official documentation and a focused probe when needed. The
  implementer/reviewer supplies this evidence; the user need not disprove the claim.
- Read tests and guidance critically: they can encode the same mistake as the implementation.
  Check deferrals before flagging absent functionality; deferral does not justify a wrong exposed
  operation. Do not change Rust solely to make an unnecessary Python convenience appear supported.
- Measure only when results could decide between semantically viable alternatives. Use focused
  checks, not a workspace-wide gate for every finding. Apply python-build before PyO3 builds or
  Python tests and verify that the imported extension matches the reviewed sources.

## Adversarial review

Use a review agent for each bounded family and an independent challenger for its findings. Keep
roles focused on Python/Rust contracts rather than creating separate agents for formatting or file
organization. The coordinator checks coverage and resolves duplicated findings. If independent
review is unavailable, report that limitation; do not label a self-check independent.

The reviewer seeks concrete violations and records the strongest evidence-based defense of each.
The challenger checks both arguments independently: cited contracts, actual consumers, reference
frames, copying claims, asserted platform limitations, and the proposed correction. It should narrow
an overbroad finding rather than dismiss a valid core. Neither role must produce findings to meet a
quota, and neither may defend behavior merely because it exists.

Keep two decisions separate:

- **Finding verdict:** confirmed, narrowed, refuted, or unresolved. Refutation needs evidence;
  failure to reproduce a suspected defect is not proof of semantic parity.
- **Deviation status:** justified, unjustified, or unresolved. Uninspected surfaces stay unreviewed.
  An inconclusive challenge never approves an unexplained deviation. Missing justification alone
  does not establish a particular runtime bug; record precisely what remains unproved.

## Output and stopping point

Use discussion-doc-writing to record results in the owning review document. Lead with a short
decision summary, followed by the coverage inventory and evidence. Each finding contains one defect:

> Rust contract → Python difference → consequence → proposed correction.

Attach source locations, the governing rule, reproduction or cost evidence, the strongest defense,
and the challenger's disposition. Retain refuted findings with reasons so they are not repeatedly
rediscovered. Keep detailed evidence accessible without forcing it into the summary.

For each reviewed deviation, propose keeping a justified adaptation, delegating/simplifying, or
removing the exposed capability. Identify affected consumers and the migration; preserving a bad
pattern through a compatibility shim is not the default. Do not design a replacement hierarchy
without a demonstrated need.

Finish the pass when every in-scope surface has an explicit coverage status and every raised finding
has a disposition. List unresolved questions and unreviewed surfaces separately. The result is a
review and proposed corrections, not an implementation plan or authorization to apply fixes.
