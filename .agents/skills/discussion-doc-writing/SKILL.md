---
name: discussion-doc-writing
description: MANDATORY — load and apply before creating, restructuring, closing, or materially editing a dated document under discussion/, changing its lifecycle status, adding related-document links, or registering it in discussion/000-status.md. Also apply when turning an initial scope into a design record or preparing that settled design for a staged implementation plan. Enforces the repository header, filename, relationship, lifecycle, status-index, and post-closeout rules.
---

# Discussion document writing

Treat `discussion/` as a set of dated design and implementation records. These documents preserve
the reasoning and state of a work unit; they are not living developer documentation.

## Read before editing

1. Read `CLAUDE.md`.
2. Read the target document completely when revising it.
3. Read `discussion/000-status.md`; its vocabulary and index are authoritative for status.
4. Read directly related discussion documents and permanent guides only as needed to establish the
   current context. Do not infer their contents from titles.

## Create the file

Use this filename shape:

```text
NNN-short-kebab-topic-YYYY-MM-DD.md
```

- Allocate the next document number from `discussion/000-status.md`.
- Count the complete basename, including `.md`; it must be **shorter than 55 characters**.
- Keep the number and ISO date. Shorten the topic instead of omitting either.
- Use lowercase ASCII kebab case for the topic.

Start new documents with this header:

```markdown
# NNN — Concise title

Status: Proposed
Date: YYYY-MM-DD
Relates: [NNN](NNN-related-topic-YYYY-MM-DD.md),
[development guide](../docs/development/example.md)
```

Header rules:

- Keep the order `Status`, `Date`, `Relates`.
- Use one exact status from the vocabulary in `discussion/000-status.md`; do not decorate it with
  bold text or append explanatory prose.
- `Date` is the creation date and does not change during implementation or closeout.
- Include `Relates` when direct relationships exist. Use relative Markdown links, not bare document
  numbers or prose references. Wrap additional links onto following lines without reformatting
  unrelated text.
- Omit `Relates` only after checking that the document is genuinely independent.

Add the document to `discussion/000-status.md` when creating it. Keep the index row's status in sync
with the header. Use the index note for the concrete outstanding work or blocker, not for a running
history.

## Follow the document lifecycle

Use one document for the work unit unless the scope grows into a genuinely separate concern.

### 1. Scope

Begin with the motivation, required outcome, boundaries, exclusions, current evidence, and open
questions. Separate known facts from proposals. Do not present unsettled suggestions as decisions,
and do not add an implementation plan yet.

Use `Proposed` for future implementation or unresolved design work. Use `Informational` only when
the document records analysis or decisions without tracked implementation scope.

### 2. Discussion and technical design

Develop the same document as questions are resolved. Record:

- the selected semantics and public API;
- rejected alternatives only when their rationale remains useful;
- consequences across affected crates, languages, specifications, and tests; and
- genuinely open questions, clearly distinguished from settled decisions.

Write a coherent technical record, not a chat transcript or implementation diary. While the
document remains open, revise stale proposals so the operative design is unambiguous. Preserve
useful rationale without narrating every intermediate thought.

### 3. Implementation plan

Add a staged implementation plan only after the design, semantics, and names needed for the work
are settled. Load and apply the `staged-impl-plan` skill, then append the plan to the same document.
Do not use the plan to conceal unresolved design choices.

A completed plan does not make the document `In Progress`. Change the status to `In Progress` only
when implementation begins. Mark completed subitems as work proceeds so the next item is visible.
If work is genuinely blocked, use `Blocked` and state the concrete blocker in the status-index note.

### 4. Closeout

Before setting `Completed`:

1. Verify every agreed, non-deferred item is implemented and checked.
2. Move deferred or newly discovered work into a separate `Proposed` document and link it.
3. Mark the document header and status-index row `Completed`.
4. Set the index `Last Checked` date to the closeout date and make its note concise.

## Preserve closed documents

After closeout, treat the document as a historical artifact:

- Do not rewrite its semantics when later decisions differ.
- Do not migrate old crate, type, method, field, or terminology names after a rename.
- Do not perform cosmetic reflow or otherwise modernize the prose.
- Record current policy in a permanent guide or a new discussion document instead, and connect the
  records with links.

Adding a `Relates` link after closeout is allowed, including a forward link to superseding or
follow-up work. Make only the link change; do not update the surrounding historical account.
`discussion/000-status.md` may reclassify a document as `Superseded` or `Outdated` and identify the
replacement without rewriting the closed document.

Permanent guides under `docs/development/` are the living normative surface. Discussion documents
must not be cited from source comments or public rustdoc.

## Editing discipline

- Make narrow diffs. Do not wrap or reformat untouched lines.
- Preserve the document's established heading and list style unless restructuring is the task.
- Update related links reciprocally when the relationship is important for discovering either work
  unit, including by adding the permitted link to a closed document.
- Run `git diff --check` after edits.
