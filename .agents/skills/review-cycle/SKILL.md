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
