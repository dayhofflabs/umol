## Proposed Design Sketch

### 1. Parsing Context Structure
```rust
pub struct PropertyParsingContext {
    // Key: (sgroup_index, field_name), Value: accumulated content parts
    sgroup_data_buffers: HashMap<(usize, String), Vec<String>>,
}

impl PropertyParsingContext {
    pub fn new() -> Self {
        Self {
            sgroup_data_buffers: HashMap::new(),
        }
    }
    
    pub fn finalize(self, molecule: &mut Molecule) -> Result<()> {
        // Validate that no incomplete SCD sequences remain
        if !self.sgroup_data_buffers.is_empty() {
            return Err(ValidationError::IncompleteDataSequences);
        }
        Ok(())
    }
}
```

### 2. Modified Apply Trait
```rust
pub trait Apply {
    fn apply(self, molecule: &mut Molecule, context: &mut PropertyParsingContext) -> Result<()>;
}
```

### 3. SGroupDataEntry Apply Implementation
```rust
impl Apply for SGroupDataEntry {
    fn apply(self, molecule: &mut Molecule, context: &mut PropertyParsingContext) -> Result<()> {
        ensure_sgroup(molecule, self.sgroup_index)?;
        
        // Get the field name from the existing SGroupData description
        let field_name = get_current_field_name(molecule, self.sgroup_index)?;
        let buffer_key = (self.sgroup_index, field_name.clone());
        
        if self.is_end {
            // SED entry - finalize the data
            let mut content_parts = context.sgroup_data_buffers
                .remove(&buffer_key)
                .unwrap_or_default();
            content_parts.push(self.data_content);
            
            // Concatenate all parts, trim to 200 chars, validate
            let full_content = content_parts.join("");
            let final_content = full_content.trim_end().chars().take(200).collect::<String>();
            
            // Add to the SGroup's data
            let sgroup = molecule.sgroups.get_mut(&self.sgroup_index).unwrap();
            if let Some(data) = sgroup.data.get_mut(&field_name) {
                data.data_content = Some(vec![final_content]);
            }
        } else {
            // SCD entry - buffer the content
            context.sgroup_data_buffers
                .entry(buffer_key)
                .or_default()
                .push(self.data_content);
        }
        
        Ok(())
    }
}
```

### 4. Parser Integration
```rust
// In the main parsing loop
let mut context = PropertyParsingContext::new();

for property_entry in property_entries {
    property_entry.apply(&mut molecule, &mut context)?;
}

// After all properties are processed
context.finalize(&mut molecule)?;
```

### 5. Other Apply Implementations
```rust
// Most other implementations just ignore the context
impl Apply for ChargeEntry {
    fn apply(self, molecule: &mut Molecule, _context: &mut PropertyParsingContext) -> Result<()> {
        // existing logic unchanged
    }
}

impl Apply for Vec<ChargeEntry> {
    fn apply(self, molecule: &mut Molecule, context: &mut PropertyParsingContext) -> Result<()> {
        for entry in self {
            entry.apply(molecule, context)?;
        }
        Ok(())
    }
}
```

### 6. Error Handling
```rust
pub enum ValidationError {
    // ... existing variants
    IncompleteDataSequences,
    MissingDataDescription(usize), // SGroup index
    DataSequenceWithoutDescription,
}
```

## Key Design Points

1. **State Management**: The context holds temporary buffers that are built up through SCD entries and finalized on SED entries.

2. **Field Name Resolution**: The field name comes from the preceding SDT (SGroupDataDescriptionEntry), which must have been processed first.

3. **Content Assembly**: Each SCD contributes up to 69 characters, SED provides the final chunk and triggers processing.

4. **Validation**: The `finalize()` method ensures no incomplete sequences remain, and individual applies validate proper sequencing.

5. **Memory Management**: Buffers are cleaned up as they're completed, and the entire context is discarded after parsing.

6. **Error Recovery**: Clear error messages for malformed data sequences, missing descriptions, etc.

This design keeps parsing state explicit and contained while handling the complex multi-line assembly requirements of the MOL format specification.