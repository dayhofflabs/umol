---
name: discussion-doc-writing
description: MANDATORY — use for material edits to dated discussion documents, lifecycle status, related-document links, closeout, implementation-plan insertion, or registration in discussion/000-status.md. Enforces naming, headers, relationships, lifecycle, indexing, and post-closeout rules.
---

# Discussion documents

`discussion/` contains dated design and implementation records, not living developer documentation.

## Before editing

- Follow `AGENTS.md`; do not reread it if already loaded this session.
- Read the target header, changed sections, and necessary adjacent or linked context. Read the whole
  document only for restructuring or closeout that requires a full audit.
- Read the status vocabulary and affected/directly related rows in `discussion/000-status.md`, not
  the whole index. Read related documents or permanent guides only as needed; never infer content
  from titles.

## Creation and index

Use `NNN-short-kebab-topic-YYYY-MM-DD.md`: allocate the next number from the status index, retain the
number and ISO date, use lowercase ASCII kebab case, and keep the complete basename at most 55
characters.

```markdown
# NNN — Concise title

Status: Proposed
Date: YYYY-MM-DD
Relates: [NNN](NNN-related-topic-YYYY-MM-DD.md),
[development guide](../docs/development/example.md)
```

- Keep header fields in that order. Use an exact status from the index vocabulary, without styling
  or commentary. The creation date never changes.
- Use relative links in `Relates`; wrap additional links without reflowing unrelated text. Omit the
  field only after confirming that no direct relationship exists.
- Register a new document immediately and keep header/index status synchronized.
- For open work, the index note states only the current action or blocker. Completed rows have an
  empty note; outcomes belong in the document.

## Lifecycle

Use one document per work unit unless a genuinely independent concern emerges.

1. **Scope:** record motivation, outcome, boundaries, exclusions, evidence, and open questions.
   Separate facts from proposals; do not add an implementation plan. Use `Proposed` for unresolved
   or future implementation and `Informational` only for analysis without tracked implementation.
2. **Design:** revise stale proposals into a coherent record of settled semantics/API, useful
   rejected alternatives, cross-crate/language/spec/test consequences, and remaining questions. Do
   not preserve a chat transcript or implementation diary.
3. **Plan:** only after semantics and names are settled, apply `staged-impl-plan` and append the plan
   to the same document. A plan does not change status; use `In Progress` only once implementation
   starts. Mark completed subitems. Use `Blocked` only for a concrete blocker and put it in the index
   note.
4. **Closeout:** verify all agreed non-deferred work; move deferred/new work to a linked `Proposed`
   document; set header and index to `Completed`; update `Last Checked` to the closeout date; clear
   the index note.

## Closed documents

- Preserve their historical semantics, names, and prose. Do not modernize or silently rewrite them.
- Current policy belongs in `docs/development/`; source comments and public rustdoc must not cite
  discussion documents.
- A narrow correction may be appended under a dated, explicit addendum that identifies the corrected
  conclusion and leaves the original intact. Broad redesign or reopened implementation needs a new
  document.
- A post-closeout `Relates` link is allowed; change only the link unless adding a permitted addendum.
- The index may mark a document `Superseded` or `Outdated` and link its replacement. Use
  `Superseded` only when the operative whole was replaced. If only part was superseded, retain the
  ordinary status and add a dated notice at the top identifying that part and its replacement.

## Editing

- Make narrow diffs; preserve untouched wrapping, heading style, and list style.
- Add reciprocal `Relates` links when either document should lead readers to the other, including
  permitted link-only changes to closed documents.
- Run `git diff --check`.
