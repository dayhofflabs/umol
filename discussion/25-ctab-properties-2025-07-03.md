# Property parsing functions for MOL files.

Implementation status for the Property Block MOL v2000 file
https://en.wikipedia.org/wiki/Chemical_table_file

|-------------------------------------------------------------------------------------------------------|
| Property   | Symbol | Implementation | Class   | RDKit* | ChemAxon** | CDK*** | Notes                 |
|------------|--------|----------------|---------|--------|------------|--------|-----------------------|
| Atom Alias | A      | x              | ISIS    | x      | x          | x      |                       |
| Atom Value | V      | x              | ISIS    | x      | x          | x      |                       |
! Group Abbr | G      | -              | ISIS    | -      | -          | -      | Outdated, use `M SUP` |
| Charge     | CHG    | x              | Generic | x      | x          | x      |                       |
| Radical    | RAD    | x              | Generic | x      | x          | x      |                       |
| Isotope    | ISO    | x              | Generic | x      | x          | x      |                       |
| Ring Bonds | RBC    | -              | Query   | x      | x          | x      |                       |
| Subs Count | SUB    | -              | Query   | x      | x          | x      |                       |
| Unsat Atom | UNS    | -              | Query   | x      | x          | x      |                       |
| Link Atom  | LIN    | -              | Query   | x      | x          | x      |                       |
| Atom List  | ALS    | -              | Query   | x      | x          | x      |                       |
| Att Point  | APO    | -              | RGroup  | x      | x          | x      |                       |
| Att Order  | AAL    | -              | RGroup  | -      | -          | x      |                       |
| Lab Loc    | RGP    | -              | RGroup  | x      | x          | x      |                       |
| Logic      | LOG    | -              | RGroup  | -      | x          | x      |                       |
| Sgrp Type  | STY    | x              | SGroup  | x      | x          | x      |                       |
| Sgrp Subt  | SST    | -              | SGroup  | x      | x          | x      |                       |
| Sgrp Label | SLB    | x              | SGroup  | x      | -          | x      |                       |
| Sgrp Conn  | SCN    | -              | SGroup  | x      | x          | x      |                       |
| Sgrp Expan | SDS    | -              | SGroup  | x      | x          | x      |                       |
| Sgrp Atoms | SAL    | x              | SGroup  | x      | x          | x      |                       |
| Sgrp Bonds | SBL    | x              | SGroup  | x      | x          | x      |                       |
| Sgrp Parnt | SPA    | -              | SGroup  | x      | x          | x      |                       |
| Sgrp Subs  | SMT    | -              | SGroup  | x      | x          | x      |                       |
| Sgrp Corr  | CRS    | -              | SGroup  | x      | -          | x      |                       |
| Sgrp Disp  | SDI    | -              | SGroup  | x      | x          | x      |                       |
| Sup Bd Vec | SBV    | -              | SGroup  | x      | -          | x      |                       |
| Data Flds  | SDT    | -              | SGroup  | x      | x          | x      |                       |
| Data Disp  | SDD    | -              | SGroup  | x      | x          | x      |                       |
| Data Sgrp  | SCD    | -              | SGroup  | x      | x          | x      | Continued data line   |
| Data Sgrp  | SED    | -              | SGroup  | x      | x          | x      | End of data line      |
| Sgrp Hier  | SPL    | -              | SGroup  | x      | x          | x      | Parent list           |
| Sgrp Comp# | SNC    | -              | SGroup  | x      | x          | x      |                       |
| 3D Feat    | $3D    | -              | 3D      | -      | -          | x      |                       |
| Phantom    | PXA    | -              | ISIS    | x      | -          | -      |                       |
| Sup Att Pt | SAP    | -              | ISIS    | x      | x          | -      |                       |
| Sup Class  | SCL    | -              | ISIS    | x      | -          | -      |                       |
| Regno      | REG    | -              | ISIS    | -      | -          | -      |                       |
| Sgrp Brkt  | SBT    | -              | ISIS    | x      | x          | x      |                       |
| 0-Order Bd | ZBO    | -              | Bd Ext  | x      | -          | -      | DOI:10.1021/ci200488k |
| Virt Hs    | ZCH    | -              | Bd Ext  | x      | -          | -      | DOI:10.1021/ci200488k |                     |
| Marvin SM  | MRV    | -              | Marvin  | x      | x          | -      |                       |
| Atom Label | ZZC    | -              | ADC     | -      | -          | x      |                       |
| Skip       | SKIP   | -              | Generic | ?      | -          | x      |                       |
| End        | END    | x              | Generic | x      | x          | x      |                       |
|-------------------------------------------------------------------------------------------------------|

* RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
** ChemAxon: https://docs.chemaxon.com/display/docs/formats_mdl-molfiles-rgfiles-sdfiles-rxnfiles-rdfiles-formats.md
*** CDK: https://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/io/MDLV2000Reader.html