# MOL Parser Performance Regression Analysis

## Summary

The basic atom and bond parsers have regressed significantly, becoming slower than the extended parsers. This defeats the purpose of having optimized basic parsers.

**Key Findings:**
- Basic atom parser (`atom_input`): ~48% slower (135ns → 200ns for len69)
- Basic bond parser: Similar regression pattern
- Extended parsers actually got faster (~7% improvement)
- The "optimizations" in basic parsers are now hurting performance

**Critical Issues (in order of impact):**

1. **`cond()` overhead** - Multiple conditional parsers add cost even when false
2. **Dynamic runtime calculations** - Length calculations on every parse
3. **Complex tuple composition** - More combinators and nested structures
4. **Position parser abstraction** - `fixed_width_position()` slower than separate x,y,z
5. **Multiple flag checks** - Bitfield checks at function start

**Priority Fixes:**
1. Remove all `cond()` calls from basic parsers (use fixed structure)
2. Replace `fixed_width_position()` with separate x, y, z parsers
3. Remove dynamic length calculations (use compile-time constants)
4. Simplify tuple composition (fewer combinators)

## Root Causes

### 1. Excessive Conditional Logic (`cond()` calls)

**Old version:**
```rust
// Simple, direct parsing - no conditionals
let x = fixed_width_float::<f64>(10, 4);
let y = fixed_width_float::<f64>(10, 4);
let z = fixed_width_float::<f64>(10, 4);
// ... direct tuple composition
```

**New version:**
```rust
// Multiple cond() calls add overhead even when false
let hydrogen_count = cond(
    atom_map_hcount_fields,
    map_res(...),
);
let atom_map_num = cond(
    atom_map_hcount_fields,
    fixed_width_int_in_range_opt::<u32, _>(3, 1..=999),
);
```

**Impact:** `cond()` in nom has overhead even when the condition is false. For basic parsers with `BASIC` flags, these conditions are typically false, but we still pay the cost.

### 2. Dynamic Runtime Calculations

**New version:**
```rust
let count1 = if atom_map_hcount_fields { 1 } else { 2 };
let remaining_len = input
    .len()
    .min(69)
    .saturating_sub(if atom_map_hcount_fields { 63 } else { 60 });
let count3 = remaining_len / 3;
```

**Impact:** These calculations happen on every parse, adding overhead. The old version had fixed, compile-time known field positions.

### 3. More Complex Tuple Composition

**Old version:**
```rust
map(
    (x, y, z, symbol, mass_diff, charge_radical, stereo_parity, valence, atom_map_num),
    |(x, y, z, symbol, ...)| Atom { ... }
)
```

**New version:**
```rust
map(
    (
        position,  // Single parser instead of x, y, z
        symbol,
        mass_diff,
        charge_radical,
        stereo_parity,
        terminated(hydrogen_count, extended1),  // Conditional
        terminated(valence, unused2),
        terminated(atom_map_num, extended3),  // Conditional
    ),
    |(position, symbol, ...)| (Atom { ... }, position)  // Returns tuple
)
```

**Impact:** 
- More fields in tuple (9 → 8, but with conditionals)
- Returns `(Atom, Point3D)` instead of just `Atom`
- More `terminated()` combinators

### 4. Position Parser Complexity

**Old version:**
```rust
let x = fixed_width_float::<f64>(10, 4);
let y = fixed_width_float::<f64>(10, 4);
let z = fixed_width_float::<f64>(10, 4);
// ... in tuple: (x, y, z, ...)
```

**New version:**
```rust
let position = fixed_width_position(ignore_positions);
// ... in tuple: (position, ...)
```

**Impact:** `fixed_width_position()` likely has more overhead than three separate parsers, especially when `ignore_positions` is false (common case). The function needs to check the flag and potentially construct a zero Point3D.

### 5. Multiple Flag Checks at Function Start

**New version:**
```rust
let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
let atom_map_hcount_fields = flags.contains(CtabParseFlags::ATOM_MAP_HCOUNT_FIELDS);
let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
```

**Impact:** Multiple bitfield checks add up, especially when most are false for basic parsers.

### 6. Length-Dependent Dispatch Overhead

The length-based dispatch (`match len { ... }`) is still present, but the individual parser functions are now more complex, so the dispatch overhead is less significant relative to the parsing cost.

## Constraints

**Important context:**
1. **Correctness is paramount** - The added complexity (cond combinators, etc.) ensures basic and extended parsers behave consistently across all flag combinations. This is a strict requirement.
2. **TableIR unification** - The `(Atom, Point3D)` return type is required by the unified SMILES/MOL TableIR design. Cannot revert to old structure.
3. **Total conversion time matters** - Not just parsing speed, but the full `MOL string -> table_ir::Molecule` pipeline. Wrapper structs that get converted don't help overall performance.

## Recommendations

### Option 1: Non-Dispatching Basic Parser (Recommended)

Create a unified basic parser with the same structure as the extended parser, but simpler:

**Design:**
- Single parser function (no length-based dispatch)
- Same tuple structure as extended parser (for consistency)
- Simplified symbol parsing (basic atoms only)
- Fixed field handling (no conditionals for basic-only features)
- Use same `(Atom, Point3D)` return type

**Benefits:**
- Eliminates dispatch overhead (`match len { ... }`)
- Simpler code path (no branching on length)
- Maintains correctness (same structure as extended)
- Can share more code with extended parser

**Implementation approach:**
```rust
pub fn atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Atom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp> {
    // Single unified parser, no length dispatch
    atom_input_unified(flags)
}

fn atom_input_unified<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Atom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp> {
    // Similar structure to extended_atom_input_inner but:
    // - Uses atom_symbol() instead of extended_atom_symbol()
    // - Fixed field counts (no dynamic calculations)
    // - No cond() for extended-only features
    // - Same tuple composition pattern
}
```

### Option 2: Conditional Length Dispatch

Keep length-dispatching parsers but only use them when no non-default flags are set:

**Design:**
- Fast path: Default flags → use length-dispatching parsers (current implementation)
- General path: Non-default flags → use unified parser (Option 1)

**Benefits:**
- Optimizes common case (default flags)
- Maintains correctness for all flag combinations
- Best of both worlds

**Implementation approach:**
```rust
pub fn atom_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Atom, Point3D), Error = NomError<&'inp [u8]>> + use<'inp> {
    // Fast path: default flags only
    if flags == CtabParseFlags::BASIC {
        atom_input_dispatch(flags)  // Current length-dispatching version
    } else {
        atom_input_unified(flags)   // Unified parser
    }
}
```

### Option 3: Hybrid Approach (Best Performance)

Combine both strategies:

1. **Default flags**: Use optimized length-dispatching parsers (current, but simplified)
2. **Non-default flags**: Use unified parser (Option 1)
3. **Benchmark both**: Compare unified vs dispatching for default flags

This allows us to:
- Keep the fast path if it's actually faster
- Have a fallback unified parser for correctness
- Measure which is actually better

## Implementation Strategy

### Step 1: Create Unified Basic Parser

Create `atom_input_unified()` that:
- Uses same structure as `extended_atom_input_inner()`
- Uses `atom_symbol()` instead of `extended_atom_symbol()`
- Has fixed field handling (no dynamic calculations for basic case)
- Still uses `fixed_width_position()` and returns `(Atom, Point3D)` for consistency
- No length dispatch - single code path

**Key simplifications:**
- For `BASIC` flags, we know:
  - `atom_map_hcount_fields` is false → no hydrogen_count, no atom_map_num
  - `extended_range` is false → fixed hydrogen count range
  - `skip_unused_fields` is true → can skip validation
  - These can be compile-time constants, not runtime checks

### Step 2: Benchmark Unified vs Dispatching

Compare performance:
- Unified parser with `BASIC` flags
- Current length-dispatching parsers with `BASIC` flags
- Measure across all length variants (34, 36, 39, 42, 48, 60, 69)

### Step 3: Optimize Based on Results

**If unified is faster:**
- Replace length dispatch entirely with unified parser
- Simpler codebase, better performance

**If dispatching is faster:**
- Keep both:
  - Fast path: `BASIC` flags → length dispatch
  - General path: Other flags → unified parser
- Or optimize unified parser further (maybe it can match dispatching speed)

### Step 4: Apply Same Strategy to Bonds

Create `bond_input_unified()` following same pattern:
- Single parser, no length dispatch
- Fixed field handling for basic case
- Benchmark vs current dispatching version

## Potential Optimizations for Unified Parser

Even with the unified structure, we can optimize:

1. **Runtime specialization with separate functions:**
   ```rust
   pub fn atom_input<'inp>(
       flags: CtabParseFlags,
   ) -> impl Parser<...> {
       // Call specialized version for BASIC flags
       if flags == CtabParseFlags::BASIC {
           atom_input_unified_basic()  // Specialized: no cond(), uses constants
       } else {
           atom_input_unified_general(flags)  // General: handles all flag combinations
       }
   }
   
   fn atom_input_unified_basic<'inp>(
   ) -> impl Parser<...> {
       // For BASIC flags, we know:
       // - atom_map_hcount_fields = false → no hydrogen_count, no atom_map_num
       // - skip_unused_fields = true → can skip validation
       // - extended_range = false → fixed hydrogen count range
       // Use constants instead of runtime checks:
       let count1 = 2;  // Fixed, not: if atom_map_hcount_fields { 1 } else { 2 }
       // No cond() calls - use success(None) or skip fields entirely
   }
   ```
   
   **This means:** Yes, separate functions called conditionally. The specialized `atom_input_unified_basic()` is only called when `flags == CtabParseFlags::BASIC`.

2. **Reduce cond() overhead in specialized version:**
   - In `atom_input_unified_basic()`, `atom_map_hcount_fields` is always false
   - Instead of `cond(false, parser)` which has overhead, use `success(None)` or omit the field entirely
   - Use compile-time constants for field counts (e.g., `count1 = 2` instead of `if atom_map_hcount_fields { 1 } else { 2 }`)

3. **Optimize position parsing:**
   - For `BASIC` flags, `ignore_positions` is typically false
   - Can inline position parsing or use a specialized version
   - Or optimize `fixed_width_position()` to be faster when flag is false

4. **Symbol parsing:**
   - `atom_symbol()` is already simpler than `extended_atom_symbol()`
   - Ensure it's as fast as possible (no unnecessary allocations)
   - Consider inlining common paths

5. **Reduce tuple nesting:**
   - Flatten nested `terminated()` combinators where possible
   - Minimize intermediate tuple construction
   - Consider using a custom combinator if needed

## Additional Considerations

### Why Unified Parser Should Be Faster

1. **No dispatch overhead:** Single code path vs `match len { ... }`
2. **Better branch prediction:** Linear control flow
3. **Cache locality:** Single function, better instruction cache usage
4. **Simpler structure:** Can optimize the whole parser as one unit

### When Length Dispatch Might Still Win

Length dispatch could be faster if:
- Different length variants have very different parsing needs
- The dispatch overhead is negligible compared to parsing cost
- Specialized parsers can skip significant work

But given the current regression, the dispatch overhead likely outweighs benefits.

### Hybrid Approach Benefits

If we implement conditional dispatch:
- **Fast path (default flags):** Use whichever is faster (unified or dispatch)
- **General path (other flags):** Always use unified (ensures correctness)
- **Flexibility:** Can optimize each path independently
- **Future-proof:** Easy to add more optimizations later

## Expected Performance Improvement

**Unified parser approach:**
- Eliminates dispatch overhead (~5-10ns saved)
- Single code path is more cache-friendly
- Simpler control flow (better branch prediction)
- Expected: **~15-25% faster** than current dispatching version
- Target: Get basic parsers to **~150-160ns** for len69 (vs current 200ns)
- Goal: Make basic parsers **~15-20% faster** than extended parsers

**If keeping both (conditional dispatch):**
- Fast path (default flags): Current performance or better
- General path (other flags): Unified parser ensures correctness
- Best of both worlds

## Testing Strategy

1. **Create unified parser** - Implement `atom_input_unified()` and `bond_input_unified()`
2. **Benchmark comparison:**
   - Unified parser with `BASIC` flags
   - Current length-dispatching parsers with `BASIC` flags
   - Extended parser with `EXTENDED` flags (baseline)
   - Measure across all length variants
3. **Correctness verification:**
   - Ensure unified parser handles all flag combinations correctly
   - Verify basic parsers reject extended features appropriately
   - Test edge cases (various lengths, flag combinations)
4. **Performance validation:**
   - Verify basic parsers are faster than extended parsers
   - Check that performance gap is maintained across all variants
   - Measure total conversion time (parsing + TableIR construction)
5. **Decision point:**
   - If unified is faster: Replace dispatch entirely
   - If dispatch is faster: Implement conditional dispatch (Option 2)
   - If similar: Choose simpler unified approach
