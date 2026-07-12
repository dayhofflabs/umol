---
description: Apply when asked to produce a staged implementation plan from a settled design — a discussion doc, a decisions list, a spec section, or the conclusion of a design conversation. Trigger on "staged impl plan", "plan the implementation", "break into stages/steps/subitems", "sequencing", "build order", or any request to turn settled design decisions into an ordered, buildable plan. Produces module-grouped subitems in dependence order, batched into stages that keep the build green at each stage boundary, using S0/S1 (stages) and S0a/S0b (subitems) notation. Consult before writing any multi-step implementation plan.
---

# Staged implementation plan

Turn a **settled** design into an ordered build plan. Do not use this to make design decisions — only to sequence decisions already made. Five steps.

## Method

**i. Collect + items by module (top-down).** Gather every settled decision. Convert to implementation *items* grouped by module, listed top-down (the consumer's needs make the items visible). The bottom-up ordering is step iii's job.

**ii. Split into subitems.** Break each item into *subitems* at the granularity of a **single struct / enum / trait**, or a **cohesive group of related functions or methods**. One subitem = one reviewable, independently-committable unit that carries its own tests.

**iii. Dependence order (bottom-up).** Order the subitems by dependence, foundation-first: a subitem precedes every subitem that uses its types or functions. Record each subitem's **explicit** dependencies. Foundation crates/modules precede their consumers precede the surface (e.g. graph-core → ast → dsl).

**iv. Group into stages.** Batch the dependence-ordered subitems into stages in rough dependence order, under one invariant:
- the tree **stays green after every subitem** wherever possible;
- a subitem **may go red** only when it is an unavoidably-breaking change — a public signature / return-type change, retiring a type, a rewire — red while mid-edit, green once its callers are migrated;
- the tree **MUST be green after every stage**. A breaking subitem and the caller-migration that restores green belong to the **same stage**.

**v. Write it down.** `S0`, `S1`, … for stages; `S0a`, `S0b`, … for subitems. Per subitem give: module, what it adds/changes, **additive (green)** vs **breaking (red→green)**, and explicit deps in `[dep: S0b, S1a]` form. State the critical path and any deferrable stages at the end.

## Rules

- **Additive first.** Batch all additive new-functionality (new type, new method, new test) before the breaking refactors — build the vocabulary first, localize risk to late stages.
- **Only breaking subitems may go red.** Changing a public signature, a return type, or retiring a type. Each must be followed *within its stage* by the caller-migration restoring green.
- **A stage is a green milestone.** If a stage can't end green, it's mis-cut — pull the caller-fixes in, or split differently.
- **Foundation before consumers.** graph-core before umol-ast before dsl; a shared kernel before every caller.
- **Mark deferrable stages explicitly**, placed last — unifications/optimizations not required for the core deliverable (the general path subsumes them). Note when the core works without them.
- **Every subitem carries its tests** (per the test-writing skill). "Green" means the suite passes, not just that it compiles.
- Prefer **one apply-surface change per stage** — batch a signature change and its migration; don't interleave two independent breaking edits to the same function across stages.
