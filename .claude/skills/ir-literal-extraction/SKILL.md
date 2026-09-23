---
name: ir-literal-extraction
description: MANDATORY — load and apply before creating or editing code in umol-graph or higher-level crates that extracts literal values from graph-IR forms, calls AsLit methods, adds a ground view, or introduces an operation-specific literal input type. Also apply when reviewing repeated as_lit calls, literal extraction followed by expect, non-literal failure handling, or a proposed macro/helper for graph-IR extraction. Classifies what non-literal values mean and selects the corresponding Rust control flow without hiding operation semantics.
---

# Graph-IR literal extraction

Read and follow [the canonical repository skill](../../../.agents/skills/ir-literal-extraction/SKILL.md)
in full before taking task actions. The `.agents` copy is authoritative.
