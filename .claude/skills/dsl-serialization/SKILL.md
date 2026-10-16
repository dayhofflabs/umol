---
name: dsl-serialization
description: Apply on ANY task that adds, extends, revises, or reviews DSL serialization/deserialization in a umol-workspace crate — `FromEdn`/`ToEdn` impls, `*Dsl` boundary types, compact-string DSL (`FromStr`/`Display`, winnow parsers), EDN readers/writers (tree or streaming), or the free functions that build/encode graph-IR types. Trigger whenever serde code is created or edited, or whenever the request mentions EDN, the DSL, `FromEdn`/`ToEdn`, `*Dsl` types, `parse_`/`fmt_`/`read_`/`render_` functions, `_to_edn`/`_from_edn`, streaming deserialization, or roundtrip. Covers which types must own their serde via traits vs which may use free functions, the mandatory function-naming scheme, the same-prefix rule for shared helpers, and the ban on single-use helpers. Consult before writing or editing any DSL/EDN ser/de, and re-check naming and structure on every such edit.
---

# umol DSL serialization conventions

Read and follow [the canonical repository skill](../../../.agents/skills/dsl-serialization/SKILL.md)
in full before taking task actions. The `.agents` copy is authoritative.
