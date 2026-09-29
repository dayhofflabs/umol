# 228 — SVG render configuration

Status: Completed
Date: 2026-09-15
Relates: [220](220-readable-depiction-2026-09-02.md),
[221](221-depiction-api-2026-09-03.md)

## Purpose

Doc 221 fixed `Depiction::render_svg() -> String` as the one rendering entry point, taking nothing.
Two rendering properties are hard-coded in the renderer and both bind a downstream consumer that
embeds umol SVG in an HTML page and also saves the same bytes as a standalone file.

- The atom-label mask is always `<mask id="umol-atom-label-mask">`, and every masked stroke group
  carries `mask="url(#umol-atom-label-mask)"`. Two depictions inlined in one document share the id,
  so the second depiction's strokes are clipped by the first depiction's mask (observed 2026-09-09).
- Mapping-index text — the correspondence-pair numbers beside reaction atoms — is drawn at 85% of
  the atom-label size (`0.3825` against `0.45`). A page stylesheet can override the presentation
  attribute on screen, but a saved SVG has no page stylesheet, so the file keeps the renderer's size.

Neither is a defect in the renderer's defaults. The mask id is a fixed name because a single
depiction needs no other, and the mapping-index ratio is a deliberate choice from doc 220. Both are
properties of the emitted text that a caller composing several depictions, or embedding at an
unusual scale, needs to control. Layout configuration (`DepictConfig`) is the wrong home: it
configures layout generation, and a `Depiction` is already laid out by the time it is rendered.

## Design

A rendering configuration type and a configured rendering variant, in Rust and Python.

**Rust.** `depict::SvgConfig` has public fields `mapping_index_text_size: f64` and
`mask_id: String`; `Default` yields the renderer's former constants. `Depiction::render_svg_with(&self,
config: &SvgConfig) -> String` renders with the configuration. `Depiction::render_svg()` is
`render_svg_with(&SvgConfig::default())` and emits byte-identical text to the previous
implementation. The renderer takes the configuration as a parameter; the mask id is written through
the existing attribute-value escaping.

**Python.** `SvgConfig(*, mapping_index_text_size=..., mask_id=...)` is a frozen, keyword-only
class with both defaults, `SvgConfig.default()`, getters, equality, and a repr that reads
`SvgConfig.default()` for the default value, following `DepictConfig`. The constructor raises
`ValueError` for a non-finite or non-positive size and for an empty mask id.
`Depiction.render_svg_with(config)` returns SVG text; `render_svg()` and `_repr_svg_()` are
unchanged.

Rejected alternatives for the mask id: a content-derived id (hash of the mask rectangles) is
deterministic and correct but changes every emitted SVG and gains nothing over an explicit id once a
configuration type exists; a per-process counter is non-deterministic. Rejected for the size: a
general text scale, which would also scale atom labels that the caller has no reason to change here.

## Scope

- `umol-io`: `SvgConfig`, `Depiction::render_svg_with`, renderer plumbing, exact tests that the
  default configuration reproduces `render_svg()` and that a custom configuration reaches the mask
  definition, every masked group, and only mapping-index text.
- `umol-py`: `SvgConfig` class, `Depiction.render_svg_with`, export, constructor rejection tests,
  and a rendered-output test.
- No default output changes. Release notes are written at release.
