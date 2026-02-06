# Plan: Data Structures for RG/LOG (Markush) and @:/@@:/THB/TLB/TEB (Extra Stereo)

Brief plan for storing the remaining CXSMILES tags. Recursive structures out of scope.

## 1. Markush (RG: / LOG:)

Store in `CxAnnotationData` using RGroup structs. Do **not** touch CtfileData. `RgroupsView` chains ctfile_data.rgroups and cx_data.rgroups (same pattern as SgroupsView). `rgroups_mut()` uses cx_data when no ctfile_data.

### LOG: — R-group logic

Maps to `cx_data.rgroups`. The `RGroup` struct has:
- `label`, `dependent_label`, `rgroup_or_h`, `occurrence`

### RG: — Member structures

CXSMILES embeds member structures as `{...}` blocks (SMILES strings). Store as raw strings for roundtripping.

```rust
// In CxAnnotationData
pub rgroups: BTreeMap<u32, RGroup>,
pub rgroup_members: BTreeMap<u32, Vec<String>>,
```

---

## 2. Extra stereo (@: / @@: / THB: / TLB: / TEB:)

Store in `CxAnnotationData` for roundtripping. Semantic mapping deferred.

### Local parity (@: / @@:)

One index for the chiral center and an ordered list of substituents that define @ or @@. Use `Chirality::Clockwise` and `Chirality::CounterClockwise` from atom.rs.

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct LocalParityEntry {
    pub center: u32,
    pub substituents: Vec<u32>,
    pub chirality: Chirality,  // Clockwise (@) or CounterClockwise (@@)
}

// In CxAnnotationData
pub local_parity: Option<Vec<LocalParityEntry>>,
```

**Global vs local parity** (ChemAxon):

- **Global parity**: CIP-like rules. Uses full molecular structure to assign priority. Nonzero when 4 different ligand types (or different stereo environments). Used for absolute configuration. ODD = clockwise, EVEN = counterclockwise when viewing with H behind the plane.
- **Local parity**: Uses atom block position only (no CIP). Numbering by atom index; H gets highest. Nonzero when implicit+explicit H < 3 and implicit H < 2. A center can have nonzero local but zero global parity (e.g. meso-like structure where CIP gives no chirality, but drawing order implies a direction).

### Bicyclo stereo (THB: / TLB: / TEB:)

In `atom.rs`, analogous to BondStereo. Variants: TowardsHigherBridge, TowardsLowerBridge, TowardsEitherBridge.

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BicycloStereoData {
    pub ligand_atom: u32,
    pub connection_atom: u32,
    pub lower_bridge_atoms: Vec<u32>,
    pub higher_bridge_atoms: Vec<u32>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BicycloStereo {
    TowardsHigherBridge(BicycloStereoData),
    TowardsLowerBridge(BicycloStereoData),
    TowardsEitherBridge(BicycloStereoData),
}

// In CxAnnotationData
pub bicyclo_stereo: Option<Vec<BicycloStereo>>,
```

---

## 3. Summary

| Tag   | Storage location        | Structure                                      |
|-------|-------------------------|------------------------------------------------|
| LOG:  | cx_data.rgroups         | RGroup (occurrence, rgroup_or_h)               |
| RG:   | cx_data.rgroup_members  | BTreeMap<u32, Vec<String>>                     |
| @:    | cx_data.local_parity    | LocalParityEntry { center, substituents, Chirality::Clockwise } |
| @@:   | cx_data.local_parity    | same, Chirality::CounterClockwise              |
| THB:  | cx_data.bicyclo_stereo  | BicycloStereo::TowardsHigherBridge             |
| TLB:  | cx_data.bicyclo_stereo  | BicycloStereo::TowardsLowerBridge              |
| TEB:  | cx_data.bicyclo_stereo  | BicycloStereo::TowardsEitherBridge             |
