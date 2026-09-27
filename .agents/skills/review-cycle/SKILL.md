---
name: review-cycle
description: Apply when asked to run a review cycle — a structured review-and-refutation pass over a named crate, module group, or type family that produces a dated review document and changes no code. Trigger on "review cycle", "adversarial review", or a request to audit construction and fallibility, nomenclature, module structure, tests and generators, or documentation against the repository guides. Covers the pinned-snapshot and worktree isolation rules, the review and refutation roles with their weighted stance, per-area agent obligations, and the review-document output.
---

# Review cycle

Run a structured, self-critical review of existing code without changing it. The normative
standard is `docs/development/code-reviews.md`; read it in full before acting, together with the
living guides and skills relevant to the reviewed area.

## Target and snapshot

- One cycle covers one named target: a crate, a module group, or a cross-module type family. The
  invoker names the target; this skill does not select one.
- Pin one commit for the whole cycle and record it in the review document. Create each agent's
  worktree explicitly with `git worktree add --detach <path> <pinned-commit>` and pass the path
  in the agent's instructions; do not rely on tooling that snapshots the current working state,
  which in a concurrently developed checkout is not the pinned commit. Remove the worktrees once
  the review document is written.
- Precondition: verify that every normative source the cycle will cite — the guides, the
  relevant skills, and this skill — exists at the pinned commit. If any is missing, stop and
  have it committed before launching agents.
- Never build, edit, or sample in the primary checkout; concurrent development there must remain
  untouched.

## Roles

- One review agent per review area of the guide (construction/integrity/fallibility,
  nomenclature, module structure and visibility, tests and generators, documentation); a small
  target may merge areas. Review agents read code, run tests, and sample generators in their own
  worktree only.
- Assign the production-use and method-boundary checks below to the structure/visibility area,
  and the name/behavior checks to nomenclature. Include their evidence in refutation; these
  checks do not require additional agents.
- Each review agent's instructions include: the must-read set (`code-reviews.md`, the guides and
  skills for its area, the governing discussion documents located through
  `discussion/000-status.md`, and `docs/umol-whitepaper.pdf` for design intent); the
  both-arguments obligation with the deferral check, whose defense names the documents
  consulted; the evidence requirements; and the one-defect-per-finding rule — a multi-part
  finding forces part-by-part adjudication and duplicate counting, so each defect is pushed
  separately.
- One refutation agent processes the pooled findings after the review agents finish, following
  the guide's steelman-first procedure and graded verdicts. It may consult git history as
  evidence, for example to check a discussion document's implementation record against what
  actually landed.
- Verification is premise-level: the refutation agent confirms that every cited normative source
  exists and supports the claim as cited, and checks the factual premises of both recorded
  arguments — claimed deferrals, claimed consumers, claimed history — rather than trusting them.
  A finding cannot survive on a citation that does not exist or does not support it. If a genuine
  normative basis exists, the verdict carries the corrected citation and records the correction;
  if none exists, the finding is refuted regardless of how plausible the defect reads.
- The orchestrating session synthesizes the surviving findings into the review document under
  the `discussion-doc-writing` skill and registers it as Proposed in
  `discussion/000-status.md`.

## Production use, names, and method boundaries

Apply these checks to private functions and methods as well as public APIs.

**Production callers.** Trace callers to production entry points across crates, bindings, trait
implementations, callbacks, macros, and relevant feature configurations. A wrapper counts only
when its own call chain has a production consumer. Separate this evidence from calls in tests,
benchmarks, fuzz targets, examples, and documentation; those calls do not establish production
use. Search hits and compiler dead-code warnings alone are insufficient.

For methods with no identified production caller, examine why they exist and whether removal
or inlining at the remaining callers would simplify the API. Tests alone do not justify a
production method or widened visibility. An intentionally supported external API or required
trait implementation can justify retention without an in-repository caller: cite its contract
and governing decision, and report unknown external use honestly. Do not invent future consumers.

**Names and behavior.** Compare each name with its implementation, inputs, output, mutation,
ownership, and failure behavior. Check the nomenclature guide and equivalent operations across
the repository, not just nearby spellings. Flag confusing names, misleading verbs or qualifiers,
and terms used outside their defined meaning even when there is no synonym or collision. For
example, normalization does not justify naming an operation canonical. Existing names, private
visibility, and historical usage are not defenses for a current rule violation. Propose a name
for the actual operation rather than a comment explaining away the mismatch.

**Methods split too finely.** Look for chains of small methods that production callers always
invoke together in the same order. Trace the intermediate values and obligations: does the split
expose a meaningless partial operation, force callers to coordinate invariants, or repeat
allocation, validation, or traversal? Consider one coherent operation where it removes that
burden. Co-use alone does not justify merging independent storage primitives or other methods
with a distinct semantic purpose.

**Methods combining unrelated work.** Examine large methods and shared helpers for independent
responsibilities bundled only to avoid duplication. Caller-specific flags, optional arguments,
discarded outputs, unrelated storage updates, and names listing several operations are evidence
to investigate, not automatic findings. Trace which responsibilities each caller actually needs.
Compare separating them and using direct caller composition, including modest duplication,
against keeping the bundle. Sharing or fewer lines is not a sufficient justification; a real
atomicity or invariant boundary can be. Preserve the semantics of one operation when proposing
a split, including its failure and rollback boundary.

For each finding, record the relevant production call chains, the claimed semantic boundary,
the concrete cost or caller burden, and the proposed removal, inline, merge, split, or rename.
The refutation pass verifies both caller classifications and the strongest reason to retain the
existing boundary. A long method or short helper is not a defect merely because of its size.

## Output

The review document contains scope, objective, per-area findings with verdicts and both recorded
arguments, refuted findings with their dismissing citations (so they are not re-flagged in later
cycles and each disposition can be audited), open items, and a proposed design. It contains no
staged implementation plan.

## Prohibitions

- No code edits anywhere, including the review worktrees; generator-sampling instrumentation is
  discarded with its worktree.
- No finding without its recorded defense argument and deferral check.
- No builds or test runs in the primary checkout.
