# Property parsing functions for MOL files.

Implementation status for the Property Block MOL v2000 file
https://en.wikipedia.org/wiki/Chemical_table_file

|-----------------------------------------------------------------------------------------------------|
| Property   | Symbol | Impl | Class   | RDKit* | ChemAxon@ | CDK^   | Indigo+ |  Notes                |
|------------|--------|------|---------|--------|-----------|--------|--------------------------------|
| Atom Alias | A      | x S  | ISIS    | x      | x         | x      | x       |                       |
| Atom Value | V      | x S  | ISIS    | x      | x         | x      | -       |                       |
! Group Abbr | G      | -    | ISIS    | -      | -         | -      | -       | Outdated, use `M SUP` |
| Charge     | CHG    | x S  | Generic | x      | x         | x      | x       |                       |
| Radical    | RAD    | x S  | Generic | x      | x         | x      | x       |                       |
| Isotope    | ISO    | x S  | Generic | x      | x         | x      | x       |                       |
| Ring Bonds | RBC    | x    | Query   | x      | x         | x      | x       |                       |
| Subs Count | SUB    | x    | Query   | x      | x         | x      | x       |                       |
| Unsat Atom | UNS    | x    | Query   | x      | x         | x      | x       |                       |
| Link Atom  | LIN    | x    | Query   | x ->   | x         | -      | -       | >Only has 2-bond var. |
| Atom List  | ALS    | x    | Query   | x      | x         | x      | x       |                       |
| Att Point  | APO    | x    | RGroup  | x      | x         | x      | x       |                       |
| Att Order  | AAL    | x    | RGroup  | -      | -         | x      | x       |                       |
| RGrp Label | RGP    | x    | RGroup  | x      | x         | x      | x       |                       |
| Logic      | LOG    | x    | RGroup  | -      | x         | x      | x       |                       |
| Sgrp Type  | STY    | x S  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Subt  | SST    | x S  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Label | SLB    | x S  | SGroup  | x      | -         | x      | -       |                       |
| Sgrp Conn  | SCN    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Expan | SDS    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Atoms | SAL    | x S  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Bonds | SBL    | x S  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Parnt | SPA    | X    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Subs  | SMT    | x S  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Corr  | CRS    | x    | SGroup  | x      | -         | x      | -       |                       |
| Sgrp Disp  | SDI    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sup Bd Vec | SBV    | x    | SGroup  | x      | -         | x      | x       |                       |
| Data Flds  | SDT    | -    | SGroup  | x      | x         | x      | x       |                       |
| Data Disp  | SDD    | -    | SGroup  | x      | x         | x      | x       |                       |
| Data Sgrp  | SCD    | -    | SGroup  | x      | x         | x      | x       | Continued data line   |
| Data Sgrp  | SED    | -    | SGroup  | x      | x         | x      | x       | End of data line      |
| Sgrp Hier  | SPL    | -    | SGroup  | x      | x         | x      | x       | Parent list           |
| Sgrp Comp# | SNC    | -    | SGroup  | x      | x         | x      | -       |                       |
| 3D Feat    | $3D    | -    | 3D      | -      | -         | x      | x       |                       |
| Phantom    | PXA    | -    | ISIS    | x      | -         | -      | -       |                       |
| Sup Att Pt | SAP    | -    | ISIS    | x      | x         | -      | x       |                       |
| Sup Class  | SCL    | -    | ISIS    | x      | -         | -      | x       |                       |
| Regno      | REG    | -    | ISIS    | -      | -         | -      | -       |                       |
| Sgrp Brkt  | SBT    | -    | ISIS    | x      | x         | x      | x       |                       |
| 0-Order Bd | ZBO    | -    | Bd Ext  | x      | -         | -      | -       | DOI:10.1021/ci200488k |
| Virt Hs    | ZCH    | -    | Bd Ext  | x      | -         | -      | -       | DOI:10.1021/ci200488k |
| Marvin SM  | MRV    | -    | Marvin  | x      | x         | -      | x       |                       |
| Atom Label | ZZC    | -    | ADC     | -      | -         | x      | -       |                       |
| Skip       | SKIP   | -    | Generic | ?      | -         | x      | -       |                       |
| End        | END    | x    | Generic | x      | x         | x      | x       |                       |
|------------------------------------------------------------------------------------------------------|

S Parsed by standard parser
* RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
@ ChemAxon: https://docs.chemaxon.com/display/docs/formats_mdl-molfiles-rgfiles-sdfiles-rxnfiles-rdfiles-formats.md
^ CDK: https://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/io/MDLV2000Reader.html
+ Indigo: https://github.com/epam/Indigo/blob/master/core/indigo-core/molecule/src/molfile_loader.cpp