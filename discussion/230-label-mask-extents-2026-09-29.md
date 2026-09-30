# 230 — Label mask extents

Status: Completed
Date: 2026-09-29
Relates: [220](220-readable-depiction-2026-09-02.md),
[225](225-depiction-problems-2026-09-09.md),
[228](228-svg-render-configuration-2026-09-15.md)

## Purpose

The SVG renderer knocks bonds out beneath each visible atom label with one rectangle in the label
mask. The rectangle was estimated from fixed constants introduced in `6a7290eb8` ("Improve atom
labels", 2026-09-03): every base character cost 0.36, every script character 0.252, and the base
half-height was 0.2475, all at label size 0.45 (script 0.315), plus 0.08 clearance on each side.

In a downstream canvas the knockout under `NH₂` visibly obscured the bond leading to it. Measured
against real sans-serif fonts, the rectangle was 1.4–2.8 times the glyph ink width and twice the
cap height, so bonds stopped well short of the letters. Vertical extent was the larger error: a
skeletal bond usually meets a label at about 30°, and the rectangle cuts it where it crosses the box
edge.

## Constraint kept from doc 220

Doc 220 rejected a glyph-shaped mask, because it leaves counters and inter-glyph gaps open to
underlying bonds, and required one continuous, conservative rectangle per structured label. The
estimate stays renderer-private, with no font measurement or label bounds in the public scene model
and no public configuration of mask size. This document keeps all of that. It changes only how the
rectangle is estimated, and it rounds the rectangle's corners.

## Measurement

Two reference fonts: DejaVu Sans 2.37 (the usual Linux `sans-serif`) and Liberation Sans 2.1.5
(metric-compatible with Arial). Both use 2048 units per em. Advances and outline extents were read
with fontTools 4.65.0. Chrome places `dominant-baseline="central"` midway between the hhea ascent
and descent lines: 709 units above the alphabetic baseline in DejaVu Sans and 710 in Liberation
Sans.

For each character of the label alphabet the renderer produces (ASCII letters, digits, `+`, U+2212,
U+2022), the new horizontal cost is the larger of the two fonts' advances. Advance is the right
measure because `text-anchor="middle"` centers by advance. Vertical extents are measured about the
central baseline and also take the larger of the two fonts:

| Extent (font units) | DejaVu Sans | Liberation Sans | Used |
| --- | ---: | ---: | ---: |
| Letter top above central | 847 | 774 | 847 |
| Letter bottom below central, no descender | 738 | 730 | 750 |
| Letter bottom below central, `g j p q y J Q` | 1135 | 1135 | 1160 |
| Script top above central | 811 | 720 | 811 |
| Script bottom below central | 738 | 730 | 750 |
| Script shift, half the line height | 1192 | 1144 | 1192 |

The 750 and 1160 depths come from Chrome 151 on macOS, whose `sans-serif` font is neither reference
font. Canvas `measureText` in that browser gave 0.366 em and 0.566 em, slightly deeper than both
reference fonts. Every character advance in that font was within the table, and its letter and script
ascents were below the table's.

The same browser, through `getBBox` on each tspan, placed a `baseline-shift` subscript or
superscript's central baseline 0.177–0.183 from the base's. That is half the script font's line
height, not the 0.1125 drop and 0.1575 rise the old estimate assumed. The old box therefore stopped
short of a subscript's ink.

## Design

- Width is the sum of per-character advances times the text size, plus clearance. A character
  outside the table costs the widest entry (`W`). This is a sizing rule, not an error path:
  rendering stays infallible, and the over-coverage matches the old behavior.
- Top is the base's letter extent, or the script shift plus script ascent when a superscript is
  present, whichever is larger. Bottom is the base's letter extent (deeper when the base has a
  descender letter), or the script shift plus script depth when a subscript is present, whichever is
  larger. Clearance is added to both. The box is no longer symmetric about the label anchor, because
  the ink is not.
- Clearance stays 0.08 on each side.
- The mask rectangle has `rx` and `ry` of 0.24, three times the clearance. A corner of the glyph box
  stays covered while the radius is at most clearance × √2 / (√2 − 1) ≈ 0.273, so rounding never
  uncovers ink. Visual review with sharp corners showed bond ends cut to the rectangle's corner. A
  radius equal to the clearance softened it too little; three times the clearance was chosen by eye.

## Evidence

The same fixtures were rendered with `depict_dsl` on `main` and on the proposing branch, then
reviewed side by side in Chrome 151 on macOS, on white and nonwhite backgrounds. The knockouts are
smaller and nothing else looked wrong, but bond ends still met the rectangle so that its corner was
visible. The rounded corners answer that. Label boxes before and after, width × height:

| Label | Before | After |
| --- | --- | --- |
| NH₂ | 1.132 × 0.693 | 1.035 × 0.645 |
| O | 0.520 × 0.655 | 0.514 × 0.511 |
| OH | 0.880 × 0.655 | 0.853 × 0.511 |
| NH | 0.880 × 0.655 | 0.835 × 0.511 |
| Cl | 0.880 × 0.655 | 0.610 × 0.511 |
| Br | 0.880 × 0.655 | 0.654 × 0.511 |
| Mg | 0.880 × 0.655 | 0.834 × 0.601 |
| O⁻ | 0.772 × 0.738 | 0.778 × 0.633 |
| NH₃⁺, NH₄⁺ | 1.384 × 0.776 | 1.299 × 0.767 |
| ¹³CH₃ | 1.636 × 0.776 | 1.425 × 0.767 |

`O⁻` and `W` become slightly wider: U+2212 and `W` are wider than the old flat per-character costs,
so the old box under-covered them. The SVG viewport follows the label box as before, so drawings
whose extreme item is a label have a slightly smaller view box.

The SVG benchmark was not rerun. The change replaces a multiplication with a table lookup per
character and adds two attributes per mask rectangle, with no new elements.

## Scope

- `umol-io`: `atom_label_box` estimates, the mask rectangle's corner radii, and exact tests for the
  boxes of `O`, `NH₂`, `Cl`, `F`, `Mg`, `¹³C`, `NH₄⁺`, `O⁻`, and an out-of-table character. A
  second test checks that the box is never narrower than either reference font's advance width plus
  clearance, and strictly narrower than before for labels with a character narrower than `N`. The
  mask-rectangle, molecule view box, and reaction mask-width pins move to the new values.
- No public type, configuration, or API changes.
- Out of scope: label orientation for bonds arriving from the right (doc 225 § 3), which also needs
  its bounds and masking updated when it is settled.
