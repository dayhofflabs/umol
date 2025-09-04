
## WONTFIX

### **🔴 Group 1: $RGP Block Issues (44 files)**

Files containing R-group definitions that our parser doesn't support:

**CDK (1 file):**
- `cdk/rgfile.7.mol` - Contains $RGP blocks

**Indigo (22 files):**
- `indigo/Overlapping_Atoms.mol` - Contains $RGP blocks  
- `indigo/Query_for_CML.mol` - Contains $RGP blocks
- `indigo/Radical.mol` - Contains $RGP blocks
- `indigo/Rgroup_for_Dearomatize.mol` - Contains $RGP blocks
- `indigo/arom_rgroup_member.mol` - Contains $RGP blocks
- `indigo/c11100_3.mol` - Contains $RGP blocks
- `indigo/c_r1.mol` - Contains $RGP blocks
- `indigo/cistrans_r.mol` - Contains $RGP blocks
- `indigo/composition1.mol` - Contains $RGP blocks
- `indigo/composition2.mol` - Contains $RGP blocks
- `indigo/composition3.mol` - Contains $RGP blocks
- `indigo/q_11.mol` - Contains $RGP blocks
- `indigo/q_rg_recurs.mol` - Contains $RGP blocks
- `indigo/q_rg_recurs2.mol` - Contains $RGP blocks
- `indigo/r1-2ap-aal.mol` - Contains $RGP blocks
- `indigo/r1-2ap.mol` - Contains $RGP blocks
- `indigo/r1_2ap.mol` - Contains $RGP blocks
- `indigo/r1_2ap_aal.mol` - Contains $RGP blocks
- `indigo/r2.mol` - Contains $RGP blocks
- `indigo/r3.mol` - Contains $RGP blocks
- `indigo/r_occur.mol` - Contains $RGP blocks
- `indigo/r_resth.mol` - Contains $RGP blocks
- `indigo/recursive2.mol` - Contains $RGP blocks
- `indigo/rgroups.mol` - Contains $RGP blocks
- `indigo/sub_mar_q01.mol` - Contains $RGP blocks
- `indigo/sub_mar_q02.mol` - Contains $RGP blocks

**Ketcher (21 files):**
- `ketcher/R-Group-structure.mol` - Contains $RGP blocks
- `ketcher/R-fragment-structure.mol` - Contains $RGP blocks
- `ketcher/Rgroup.mol` - Contains $RGP blocks
- `ketcher/all-kind-of-r-group.mol` - Contains $RGP blocks
- `ketcher/clean-different-properties.mol` - Contains $RGP blocks
- `ketcher/clean-rgroups.mol` - Contains $RGP blocks
- `ketcher/complex-r-group-structure.mol` - Contains $RGP blocks
- `ketcher/markush-expected.mol` - Contains $RGP blocks
- `ketcher/markush.mol` - Contains $RGP blocks
- `ketcher/r-group-all-chain.mol` - Contains $RGP blocks
- `ketcher/r-group-expected.mol` - Contains $RGP blocks
- `ketcher/r-group-with-allkind-attachment-points-expectedV2000.mol` - Contains $RGP blocks
- `ketcher/r-group-with-allkind-attachment-points.mol` - Contains $RGP blocks
- `ketcher/r1-several-distorted.mol` - Contains $RGP blocks
- `ketcher/r1-several-structures-expected.mol` - Contains $RGP blocks
- `ketcher/r1-several-structures.mol` - Contains $RGP blocks
- `ketcher/structure-r-group-logic.mol` - Contains $RGP blocks

## FIXED

### **🔴 Group 2: Out-of-Bounds Values (25+ files)**

**Invalid Bond Orders:**
- `rdkit/bondorder0.mol` - **Bond order 0**
- `rdkit/bondorder9.mol` - **Bond order 9**  
- `rdkit/DativeBond2000.mol` - **Bond order 9**
- `indigo/hydrogen_test.mol` - **Bond order 10**
- `ketcher/saving-and-rendering-Dative-bond-(refactored).mol` - **Bond orders 9, 9**
- `ketcher/saving-and-rendering-Hydrogen-bond-(refactored).mol` - **Bond orders 10, 10** 
- `ketcher/four-bonds.mol` - **Bond orders 8, 9**
- `ketcher/all-kinds-of-bonds-test-file.mol` - **Bond orders 5, 6, 7, 9, 10**

**Invalid Element Symbols:**
- `cdk/mdlquery.mol` - Element `L` (query atom "list")
- `indigo/aniline_pol_psd.mol` - Elements `Psd`, `Pol` (pseudoatoms)
- `indigo/pseudo_target2.mol` - Element `Pol`
- `indigo/issue269test_PseudoatomWarning.mol` - Element `HAR` 
- `indigo/ind-295-biggy1v2000.mol` - Elements `T`, `D`, `Q`, `Asn`, `Asp` (isotopes + amino acids)
- `ketcher/Source.mol` - Elements `A`, `CYC` (query atoms)
- `ketcher/chain-with-group-generics-expected.mol` - Element `GH*`
- `ketcher/chain-with-group-generics.mol` - Element `GH*`
- `ketcher/generic-groups.mol` - Element `GH*` 
- `ketcher/s-group-features.mol` - Elements `A`, `L`, `T-4` (query atoms)

**Mass Data Files:** 
- `indigo/all2000.mol` - **280+ bond order 0**, Elements `SDF`, `R#`, `T`, `Rn`, `Rb`, `D`, `Ra`, `Re`, `Ru`, `Rf`, `Rh`, `PhN`
- `indigo/kconv.mol` - **280+ bond order 0**, same element issues
- `indigo/ketcher.mol` - **280+ bond order 0**, same element issues

### **🔴 Group 3: EOF Issues (13 files)**

**Missing V2000 Tag (7 files):**
- `cdk/Strychnine_nichtOK.mol` - Missing V2000 tag
- `cdk/bug1014344-1.mol` - Missing V2000 tag
- `indigo/01b4b097_1ab8_48c4_8b39_4b33574ff5e1.mol` - Missing V2000 tag
- `indigo/1e-0.mol` - Missing V2000 tag
- `indigo/Chirality.mol` - Missing V2000 tag
- `rdkit/issue148.mol` - Missing V2000 tag
- `rdkit/unsanitary2.mol` - Missing V2000 tag

**Truncated Files (6 files):**
- `cdk/rgroupsNumbered.mol` - Unexpected end of file
- `indigo/Row2.mol` - Unexpected end of file  
- `indigo/Row3.mol` - Unexpected end of file
- `indigo/recursive1.mol` - Unexpected end of file
- `indigo/test_molv2000_charge.mol` - Unexpected end of file
- `ketcher/three-structures.mol` - Unexpected end of file

---

## **📋 FILE CATALOG WITH PARSING ISSUES



### **🔴 Group 4: Digit Parsing Issues (13 files)**

**Spacing Problems in Counts Line (7 files):**
- `cdk/hisotopes.mol` - Non-digit in counts line (`"0999 V2000"` - missing space)
- `cdk/rgfile.1.mol` - Non-digit in counts line  
- `cdk/rgfile.2.mol` - Non-digit in counts line
- `cdk/rgfile.3.mol` - Non-digit in counts line
- `cdk/rgfile.4.mol` - Non-digit in counts line
- `cdk/rgfile.5.mol` - Non-digit in counts line  
- `cdk/rgfile.6.mol` - Non-digit in counts line

**Missing Space Before V2000 (4 files):**
- `indigo/1944-3D_Structure.mol` - Missing space before V2000
- `indigo/3D_Structure.mol` - Missing space before V2000  
- `indigo/Pseudoatom.mol` - Missing space before V2000
- `indigo/Stereochemistry.mol` - Missing space before V2000

**Other Digit Issues (2 files):**
- `cdk/mdlquery.mol` - Digit parsing error
- `openbabel/nsc2dmol.mol` - Digit parsing error

---

### **🔴 Group 5: Line Ending Issues (2 files)**
- `ketcher/Custom-expected.mol` - **0 bytes** (completely empty)
- `ketcher/empty-file.mol` - **2 bytes** (just LF characters)

---

### **🔴 Group 6: Other Issues (10 files)**

**RXN Files (3 files):**
- `cdk/ethylesterification.mol` - RXN file (not MOL)
- `indigo/empty_apid.mol` - RXN file (not MOL)  
- `indigo/ket-reaction-arrow.mol` - RXN file (not MOL)

**Unknown Mapping Errors (5 files):**
- `ketcher/alias-and-pseudoatoms-expected.mol` - Unknown mapping error
- `ketcher/alias-and-pseudoatoms.mol` - Unknown mapping error
- `ketcher/mol_1852_to_open-expected.mol` - Unknown mapping error
- `marvin/20.mol` - Unknown mapping error  
- `rdkit/Issue269.mol` - Unknown mapping error

**Other Data Errors (2 files):**
- `indigo/Row5.mol` - Unknown data error (ends with `$$$$%` instead of `$$$$`)
- `rdkit/rxn2.mol` - Unknown data error (RXN file that partially parses)

## **📊 SUMMARY STATISTICS:**
- **$RGP blocks:** 55 files
- **$RXN files:** 4 files  
- **Missing V2000:** 9 files
- **Empty/tiny files:** 2 files

The analysis shows clear patterns where extended range support could help with bond orders, isotopes, and query atoms.