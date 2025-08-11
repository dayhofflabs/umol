# Expected Results

This directory contains manually curated expected parsing results in YAML format.
These serve as reference outputs for complex test cases where the correct result
is not obvious or needs to be explicitly specified.

## File Format

Expected results are stored as YAML files with the same base name as the input file:
- `example.mol` → `example.mol.yaml`

## YAML Schema

```yaml
format_version: "1.0"
source_file: "example.mol" 
metadata:
  atom_count: 2
  bond_count: 1
  # ... other metadata

atoms:
  - index: 0
    element: "C"
    coordinates: [0.0, 0.0, 0.0]
    # ... atom properties

bonds:
  - atom1: 0
    atom2: 1
    bond_type: "Single"
    # ... bond properties

# ... additional sections
```

## When to Add Expected Results

Create expected result files for:
1. Complex molecules where structure verification is important
2. Test cases with subtle property combinations
3. Real-world files that represent important use cases
4. Cases where the "correct" result might be ambiguous

## Maintenance

Expected results should be:
- Reviewed when parser logic changes
- Updated when file format specifications evolve
- Validated against multiple implementations when possible
