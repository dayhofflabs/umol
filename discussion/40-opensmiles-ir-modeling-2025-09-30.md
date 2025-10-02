### OpenSMILES IR modeling and post-parse pipeline (2025-09-30)

This note captures the high-level IR design and post-parse workflow for SMILES in umol.

#### Goals
- Keep parsing simple; defer semantics to post-parse.
- Provide explicit, swappable models for valence/aromaticity, including a null model.
- Unify SMILES and MOL/SDF ingestion into a common Simple IR.
- Layer additional structure for algorithms and typed lints without entanglement.

#### IRs
- SimpleIR (SIR)
  - atoms: { element, charge, bracket_h_opt, aromatic_flag, chirality_opt, isotope_opt, class_opt, src_span }
  - bonds: { a_idx, b_idx, order, dir_opt, aromatic_flag, src_span }
  - Purpose: stable, format-agnostic exchange; easy to build from SMILES or MOL/SDF; no inference.

- GraphIR (GIR)
  - petgraph-backed graph built from SIR.
  - Carries SideChannel (optional): ring events/spans; bracket field order, etc.
  - Purpose: algorithms/topology; used by post-parse passes; convertible back to SIR.

- SMILES-typed IR (STIR)
  - Adds SMILES-scoped typed facts derived from SIR/GIR and a ModelProfile.
  - Per-atom: degree/order_sum, explicit_h (bracket), implicit_h_suggested_opt (profile), radical_opt (future), kekule_tag_opt (future), issue flags.
  - Per-bond: normalized dir for cis/trans checks, inferred_aromatic_opt (verify-only model).
  - Purpose: power lints/verification without full valence assignment.

#### Conversions
- sir_to_gir(sir) -> gir
- gir_to_sir(gir) -> sir
- annotate_smiles_typed({sir|gir}, profile) -> stir

#### Post-parse pipeline (check_smiles)
Inputs: GIR, optional SideChannel, ModelProfile, lint-config.

Order:
1) Normalize (bond dir storage, bond defaults; no semantics).
2) Topology checks: self-loop, parallel edges (use ring side-channel for spans when present).
3) Double-bond stereo validation: per end, require ≥1 marked substituent; per-end consistency; emit insufficient/conflict.
4) Aromaticity verification (model-driven):
   - Null: skip.
   - Verify-only: enforce aromatic atoms/bonds in rings; ring consistency; no kekulization initially.
5) Typed annotations: produce STIR fields per profile (no hard inference under Null).
6) Lint engine: run enabled rules (style/numeric) and merge with parser-mapped diagnostics.

Outputs: DiagnosticsReport; optionally STIR and/or updated GIR; artifacts (ring sequence, per-bond stereo view, etc.).

#### Profiles (ModelProfile)
- Null: no inference/validation; structural checks only; STIR fields only where explicitly specified.
- OrganicStrict: minimal inference (implicit_h suggestions), per-element sanity checks; still SMILES-scoped.
- Hybrid: per-element policy map selecting behaviors for valence/aromaticity.


#### Valence checking

Got it.

### Valence checking pass design

- Basic model
  - Valence model = map Element → ordered list of allowed valence states (u8), plus per-element policy.
  - Implicit H algorithm (strict): total_valence = sum(bond orders) − charge. For each atom:
    - Find the smallest allowed valence state ≥ total_valence; if none exists, emit VALENCE_EXCEEDS_MAX (error/warn per config).
    - implicit_H = allowed_valence − total_valence. If negative, emit VALENCE_NEGATIVE_IMPLICIT (error).
  - Use a compact table, similar to RDKit’s, but with per-element toggles: enabled, rounding policy, aromatic handling, and overflow policy. Reference: RDKit’s atomic valence states table [link](https://github.com/rdkit/rdkit/blob/23ffd85f60d5cbedc86c698933f0fbaeabc81437/Code/GraphMol/atomic_data.cpp).

- Bond order accounting
  - Default: Single=1, Double=2, Triple=3, Quadruple=4, Aromatic=1 (strict-lower-bound), configurable to 2 or 1.5 later.
  - Coordination/dative donation does not change order value in this pass (keep simple).
  - Ring/aromatic verification is a separate pass; valence does not attempt to resolve kekulization.

- Config surfaces
  - enable: on/off (skip pass entirely).
  - overflow_policy: Error | Warn | Off.
  - aromatic_as: 1 (default) | 2 | 1_5 (future).
  - elements: per-element switch; if an element has no valence states, skip implicit H inference and only sanity-check totals if requested (warn-only mode).
  - bracket_behavior:
    - check_bracket: validate bracketed atoms against model (default: on).
    - infer_bracket_implicit: off by default; if on, compute implicit_H for bracket atoms without H field; if bracket has Hn, verify consistency.
  - metals_policy: treat states as advisory only (Warn) unless explicitly flipped to Error.

- Diagnostics (Category=Valence)
  - VALENCE_EXCEEDS_MAX: total_valence > highest_state.
  - VALENCE_NEGATIVE_IMPLICIT: computed implicit_H < 0 (over-satisfied valence).
  - VALENCE_UNKNOWN_FOR_ELEMENT: no states and checks required by config.
  - VALENCE_BRACKET_H_MISMATCH: bracket Hn conflicts with computed implicit_H.
  - STYLE_PREFER_IMPLICIT_H_OR_EXPLICIT (optional): advisory if both implicit suggestion and explicit Hs appear inconsistent but legal.
  - STYLE_PREFER_NO_VALENCE_OVERFLOW (advisory when overflow_policy=Warn).

- Placement in pipeline
  - Run after topology, before double-bond stereo and aromaticity.
  - Outputs:
    - For STIR: per-atom implicit_h_suggested (u8), total_valence_observed (u8), chosen_valence_state (u8).
    - Diagnostics per above.

- Q1: Is “round up to next valence state” enough?
  - For organic subset: yes (C,N,O,P,S,halogens) under strict-aromatic-as-1 default it’s predictable; hypervalent S/P are covered by states 4/6 (and higher, model-dependent).
  - For metals/electronegatives and unusual chemistries, numeric states alone are sometimes insufficient; keep the algorithm numeric now and allow per-element overflow_policy=Warn (not Error) by default. Provide a per-element override to disable implicit H inference entirely.
  - Conclusion: model includes numeric states + small per-element policy. Keep strategy hooks for future specialized elements, but do not implement general strategies now.

- Q2: Preserving expressive valence states (full electron config)
  - Add a second model tier: TypedValenceModel backed by `umol-models-valence` `AtomSpec` profiles. Adapter: NumericValenceModel (today) and TypedValenceModel (future) both implement a `ValenceModel` trait:
    - fn pick_state(atom_view) -> ValenceDecision { allowed_state, implicit_h, diagnostics }
  - The typed model can match the full electron configuration (charge, lone pairs, unpaired e−, multiplicity, valence) and produce richer diagnostics. Keep this as an optional profile; default remains numeric.

- Q3: Bracketed atoms
  - Yes, run checks for bracketed atoms:
    - If Hn is present: verify Hn equals computed implicit_H; if mismatch → VALENCE_BRACKET_H_MISMATCH (Error).
    - If Hn absent:
      - Default: do not invent implicit Hs, only sanity-check total_valence vs model; emit VALENCE_EXCEEDS_MAX/NEGATIVE_IMPLICIT where applicable.
      - Optional (infer_bracket_implicit=true): set implicit_h_suggested as for bare atoms (round-up), but do not mutate IR, only report suggestion.
  - “[CH*]” extension:
    - Treat “*” as “fill to next state” only behind an extension flag; on success, record implicit_h_suggested and emit a STYLE note if requested. Absent flag → lexical/syntax error as per parser policy.

- Minimal implementation plan (non-code)
  - Add `pass_valence(mol, profile, ctx)`:
    - Compute per-atom total_valence.
    - Lookup states; apply rounding rule; compute implicit_h_suggested.
    - Validate bracket H field if present.
    - Emit diagnostics per policy.
    - Write STIR fields: implicit_h_suggested, total_valence_observed, chosen_valence_state.
  - Start with a small built-in numeric table; expose a builder to register/override states per element; ship a StrictOrganic profile and a Null profile.

If this looks good, I’ll add the `pass_valence` interface and wire it after topology in `checker`, then seed a minimal numeric table for C,N,O,P,S,F,Cl,Br,I,Si,B,H with conservative states.

#### APIs (initial)
- parse_smiles(bytes, flags) -> Result<GIR, ParseError>
- check_smiles(gir: &GIR, side: Option<&SideChannel>, profile: &ModelProfile, lint_cfg: &LintConfig) -> DiagnosticsReport
- smiles_to_ir(input, options) -> { SIR | GIR | STIR }
  - Options select target IR kind; SIR is default (enables MOL/SDF to target the same type).

---

### Interface design (sketch)

Note: Rust-like type signatures for design; exact modules/types may vary.

```rust
// Core IRs
pub struct SimpleIr {
    pub atoms: Vec<SirAtom>,
    pub bonds: Vec<SirBond>,
}

pub struct SirAtom {
    pub element: umol_data::Element,
    pub charge: Option<i32>,
    pub bracket_h: Option<u8>,
    pub aromatic: bool,
    pub chirality: Option<io::smiles::Chirality>,
    pub isotope: Option<u32>,
    pub class_num: Option<u32>,
    pub span: io::smiles::Span,
}

pub struct SirBond {
    pub a: u32,
    pub b: u32,
    pub order: io::smiles::BondOrder,
    pub dir: Option<io::smiles::BondDir>,
    pub aromatic: bool,
    pub span: io::smiles::Span,
}

// GraphIR is the parser’s molecule (petgraph-backed or equivalent)
pub type GraphIr = io::ir::Molecule; // existing parser IR

// Optional side-channel captured during parsing
pub struct SideChannel {
    pub ring_events: Vec<RingEvent>, // best-effort anchors for spans
    pub bracket_field_order: Option<Vec<(usize /*atom_idx*/, &'static str /*field*/)>>,
}

pub struct RingEvent {
    pub num: u32,
    pub open_pos: usize,
    pub close_pos: Option<usize>,
    pub a_idx: Option<u32>,
    pub b_idx: Option<u32>,
}

// SMILES-typed IR (annotations layered over SIR)
pub struct StIr {
    pub sir: SimpleIr,
    pub atoms: Vec<StAtom>,
    pub bonds: Vec<StBond>,
}

pub struct StAtom {
    pub degree: u8,
    pub order_sum: u8,
    pub explicit_h: Option<u8>,
    pub implicit_h_suggested: Option<u8>, // per profile
    pub radical: Option<u8>,              // future
    pub kekule_tag: Option<u8>,           // future
}

pub struct StBond {
    pub norm_dir: Option<io::smiles::BondDir>,
    pub inferred_aromatic: Option<bool>,
}

// Conversions
pub fn sir_to_gir(sir: &SimpleIr) -> GraphIr;
pub fn gir_to_sir(gir: &GraphIr) -> SimpleIr;
pub fn annotate_smiles_typed(source: &SimpleIr, gir: &GraphIr, profile: &ModelProfile) -> StIr;

// Profiles (strategy objects)
pub struct ModelProfile {
    pub valence: Box<dyn ValenceModel>,
    pub aromaticity: Box<dyn AromaticityModel>,
    pub element_policy: Option<std::collections::HashMap<umol_data::Element, ElementPolicy>>,
}

pub struct ElementPolicy {
    pub valence: Option<Box<dyn ValenceModel>>,      // overrides per element
    pub aromaticity: Option<Box<dyn AromaticityModel>>,
}

pub trait ValenceModel: Send + Sync {
    fn annotate(&self, gir: &GraphIr, stir: &mut StIr);
}

pub trait AromaticityModel: Send + Sync {
    fn verify(&self, gir: &GraphIr) -> Vec<io::smiles::Diagnostic>; // verify-only, no kekulization for now
}

// Lint configuration (placeholder)
pub struct LintConfig {
    pub enabled_codes: Vec<&'static str>,
    pub disabled_codes: Vec<&'static str>,
}

pub enum IrTarget { Sir, Gir, Stir }

pub enum IrOutput { Sir(SimpleIr), Gir(GraphIr), Stir(StIr) }

pub struct CheckOptions<'a> {
    pub profile: &'a ModelProfile,
    pub lint: &'a LintConfig,
    pub side: Option<&'a SideChannel>,
}

// Preferred top-level APIs
pub fn check_smiles(gir: &GraphIr, opts: &CheckOptions) -> io::smiles::DiagnosticsReport;

pub struct SmilesToIrOptions<'a> {
    pub flags: io::smiles::SmilesParseFlags,
    pub target: IrTarget,
    pub profile: Option<&'a ModelProfile>, // used only for STIR
}

pub fn smiles_to_ir(input: &str, opts: &SmilesToIrOptions) -> Result<IrOutput, io::smiles::ParseError>;
```

#### Execution outline
- smiles_to_ir("...", target=Sir) → parse → gir_to_sir
- smiles_to_ir("...", target=Gir) → parse → GraphIr
- smiles_to_ir("...", target=Stir) → parse → gir_to_sir → annotate_smiles_typed
- check_smiles(&gir, opts) → normalize → topology → double-bond stereo → aromaticity verify → lints

---

### API surface details

#### Modules and placement
- `umol_models_graph::io::smiles::api` (new): top-level helpers
  - `check_smiles`, `smiles_to_ir`, `sir_to_gir`, `gir_to_sir`, `annotate_smiles_typed`
- Keep existing parser entrypoints unchanged (`parse_smiles`, `parse_smiles_inner`).

#### Defaults & config
```rust
impl Default for LintConfig {
    fn default() -> Self { Self { enabled_codes: vec![], disabled_codes: vec![] } }
}

impl Default for ModelProfile {
    fn default() -> Self { ModelProfile { valence: Box::new(NullValenceModel), aromaticity: Box::new(NullAromModel), element_policy: None } }
}

pub struct ParseOptions {
    pub flags: io::smiles::SmilesParseFlags, // STRICT_OPENSMILES by default
}
```

Common call patterns:
```rust
// 1) Parse + check (strict, no typed annotations)
let gir = parse_smiles(bytes)?;
let report = check_smiles(&gir, &CheckOptions{ profile: &ModelProfile::default(), lint: &LintConfig::default(), side: None });

// 2) Produce SIR for snapshot/export
let IrOutput::Sir(sir) = smiles_to_ir(input, &SmilesToIrOptions{ flags: SmilesParseFlags::STRICT_OPENSMILES, target: IrTarget::Sir, profile: None })? else { unreachable!() };

// 3) Produce STIR with OrganicStrict profile
let profile = organic_strict_profile();
let IrOutput::Stir(stir) = smiles_to_ir(input, &SmilesToIrOptions{ flags: SmilesParseFlags::STRICT_OPENSMILES, target: IrTarget::Stir, profile: Some(&profile) })? else { unreachable!() };
```

#### Profiles
```rust
pub fn null_profile() -> ModelProfile { ModelProfile::default() }

pub fn organic_strict_profile() -> ModelProfile {
    ModelProfile {
        valence: Box::new(OrganicStrictValenceModel),
        aromaticity: Box::new(VerifyOnlyAromModel),
        element_policy: None,
    }
}

pub fn hybrid_profile(map: std::collections::HashMap<Element, ElementPolicy>) -> ModelProfile {
    ModelProfile { element_policy: Some(map), ..organic_strict_profile() }
}
```

#### Interop with MOL/SDF
- MOL/SDF parsers should output `SimpleIr` directly.
- Downstream checks reuse `sir_to_gir` + `check_smiles` for a unified pipeline.

#### Spans and side-channel
- `Span` fields are carried in SIR and preserved in GIR via node/edge metadata.
- `SideChannel` is optional and only produced when `SmilesParseFlags::LINT_SIDECHANNEL` is set; consumers must handle `None` gracefully.

#### Thread-safety and performance
- All models implement `Send + Sync`; API functions take shared references and avoid hidden global state.
- SIR/GIR conversions are O(V+E); STIR annotation adds O(V+E) overhead proportional to enabled models.

#### Error handling
- `smiles_to_ir` returns parser `ParseError` only; lints/semantic diagnostics are emitted via `check_smiles`.
- `check_smiles` never panics and returns a `DiagnosticsReport` (possibly empty).

#### Future extensions
- Kekulé assignment API (separate from `AromaticityModel::verify`).
- Optional geometry/seating inference hooks for enhanced stereo validation.
- Exporters: `gir_to_smiles` (canonicalizer) to be designed separately.

---

### Pass interfaces (normalization, topology, stereo)

We expose passes as functions with simple inputs/outputs; a common `PassCtx` carries config and emitters.

```rust
pub struct PassCtx<'a> {
    pub profile: &'a ModelProfile,
    pub lint: &'a LintConfig,
    pub emitter: &'a mut io::smiles::linter::Emitter<'a>,
}

// Optional artifacts each pass may return for later stages or debugging
pub struct NormArtifacts {
    pub normalized_bond_dirs: bool,
}

pub struct TopologyArtifacts {
    pub self_loops: usize,
    pub parallel_pairs: usize,
}

pub struct StereoArtifacts {
    pub checked_double_bonds: usize,
    pub insufficient_count: usize,
    pub conflict_count: usize,
}

// Stage 0: Normalize (O(V+E))
pub fn pass_normalize(gir: &mut GraphIr) -> NormArtifacts;

// Stage 1: Topology checks (O(V+E))
pub fn pass_topology(
    gir: &GraphIr,
    side: Option<&SideChannel>,
    ctx: &mut PassCtx,
) -> TopologyArtifacts;

// Stage 2: Double-bond stereo validation (local around double bonds)
pub fn pass_stereo_double(
    gir: &GraphIr,
    ctx: &mut PassCtx,
) -> StereoArtifacts;

// Stage 3: Aromaticity verification (model-driven); verification only
pub fn pass_aromaticity_verify(
    gir: &GraphIr,
    ctx: &mut PassCtx,
);

// Stage 4: Typed annotations (STIR) – optional, profile-driven
pub fn pass_typed_annotations(
    sir: &SimpleIr,
    gir: &GraphIr,
    profile: &ModelProfile,
) -> StIr;

// Pipeline runner used internally by check_smiles
pub fn run_pipeline(
    gir: &mut GraphIr,
    side: Option<&SideChannel>,
    ctx: &mut PassCtx,
) {
    let _norm = pass_normalize(gir);
    let _topo = pass_topology(gir, side, ctx);
    let _stereo = pass_stereo_double(gir, ctx);
    pass_aromaticity_verify(gir, ctx);
}
```

Diagnostics emission follows existing `Emitter` usage; each pass emits stable codes with spans.

Span policy:
- Prefer ring side-channel anchors for ring-index-originated diagnostics (open_pos/close_pos).
- Fallback to bond/atom endpoints when side-channel data is absent.

Configuration:
- Passes read `ctx.profile` and `ctx.lint` to decide behaviors (e.g., Null vs Strict profile, enabled codes) but should remain fast no-ops if corresponding diagnostics are disabled.

Threading:
- Passes are pure over inputs except for `Emitter` side effects; they can be made parallel over components later.

#### Notes on umol-models-valence
- Serves as inspiration for stricter typing; not a current target for automatic mapping.
- Future: optional export STIR→valence graph behind a dedicated profile/feature.

#### Next steps
- Define ModelProfile, ElementPolicy, and pass runner interfaces.
- Implement sir_to_gir/gir_to_sir and outline STIR annotation.
- Implement topology and stereo passes.
- Add parser ring side-channel (flag-gated) for spans.


