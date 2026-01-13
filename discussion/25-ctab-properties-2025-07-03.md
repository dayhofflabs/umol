# Property parsing functions for MOL files.

Implementation status for the Property Block MOL v2000 file
https://en.wikipedia.org/wiki/Chemical_table_file

|------------------------------------------------------------------------------------------------------|
| Property   | Symbol | Impl | Class   | RDKit* | ChemAxon@ | CDK^   | Indigo+ |  Notes                |
|------------|--------|------|---------|--------|-----------|--------|---------------------------------|
| Atom Alias | A      | x B  | ISIS    | x      | x         | x      | x       |                       |
| Atom Value | V      | x B  | ISIS    | x      | x         | x      | -       |                       |
! Group Abbr | G      | x    | ISIS    | -      | -         | -      | -       | Legacy, defer valid.  |
| Charge     | CHG    | x B  | Generic | x      | x         | x      | x       |                       |
| Radical    | RAD    | x B  | Generic | x      | x         | x      | x       |                       |
| Isotope    | ISO    | x B  | Generic | x      | x         | x      | x       |                       |
| Ring Bonds | RBC    | x    | Query   | x      | x         | x      | x       |                       |
| Subs Count | SUB    | x    | Query   | x      | x         | x      | x       |                       |
| Unsat Atom | UNS    | x    | Query   | x      | x         | x      | x       |                       |
| Link Atom  | LIN    | x    | Query   | x ->   | x         | -      | -       | >Only has 2-bond var. |
| Atom List  | ALS    | x    | Query   | x      | x         | x      | x       |                       |
| Att Point  | APO    | x    | RGroup  | x      | x         | x      | x       |                       |
| Att Order  | AAL    | x    | RGroup  | -      | -         | x      | x       |                       |
| RGrp Label | RGP    | x    | RGroup  | x      | x         | x      | x       |                       |
| Logic      | LOG    | x    | RGroup  | -      | x         | x      | x       |                       |
| Sgrp Type  | STY    | x B  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Subt  | SST    | x B  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Label | SLB    | x B  | SGroup  | x      | -         | x      | -       |                       |
| Sgrp Conn  | SCN    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Expan | SDS    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Atoms | SAL    | x B  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Bonds | SBL    | x B  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Parnt | SPA    | X    | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Subs  | SMT    | x B  | SGroup  | x      | x         | x      | x       |                       |
| Sgrp Corr  | CRS    | x    | SGroup  | x      | -         | x      | -       |                       |
| Sgrp Disp  | SDI    | x    | SGroup  | x      | x         | x      | x       |                       |
| Sup Bd Vec | SBV    | x    | SGroup  | x      | -         | x      | x       |                       |
| Data Flds  | SDT    | x    | SGroup  | x      | x         | x      | x       | See bel. for MRV ext. |
| Data Disp  | SDD    | x    | SGroup  | x      | x         | x      | x       |                       |
| Data Sgrp  | SCD    | x    | SGroup  | x      | x         | x      | x       | Continued data line   |
| Data Sgrp  | SED    | x    | SGroup  | x      | x         | x      | x       | End of data line      |
| Sgrp Hier  | SPL    | x    | SGroup  | x      | x         | x      | x       | Parent list           |
| Sgrp Comp# | SNC    | x    | SGroup  | x      | x         | x      | -       |                       |
| 3D Feat    | $3D    | -    | 3D      | -      | -         | x      | x       |                       |
| Phantom    | PXA    | -    | ISIS    | x      | -         | -      | -       |                       |
| Sup Att Pt | SAP    | -    | ISIS    | x      | x         | -      | x       |                       |
| Sup Class  | SCL    | -    | ISIS    | x      | -         | -      | x       |                       |
| Regno      | REG    | -    | ISIS    | -      | -         | -      | -       |                       |
| Sgrp Brkt  | SBT    | -    | ISIS    | x      | x         | x      | x       |                       |
| Bond order | ZBO    | x B  | Bd Ext  | x      | -         | -      | -       | DOI:10.1021/ci200488k |
| Atom charge| ZCH    | x B  | Bd Ext  | x      | -         | -      | -       | DOI:10.1021/ci200488k |
| Atom Hs    | HYD    | x B  | Bd Ext  | x      | -         | -      | -       | DOI:10.1021/ci200488k |
| Marvin SM  | MRV    | x    | Marvin  | x      | x         | -      | x       |                       |
| Atom Label | ZZC    | x B  | ADC     | -      | -         | x      | -       |                       |
| Skip       | SKIP   | -    | Generic | ?      | -         | x      | -       |                       |
| End        | END    | x B  | Generic | x      | x         | x      | x       |                       |
|------------------------------------------------------------------------------------------------------|

B Parsed by basic parser
* RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
@ ChemAxon: https://docs.chemaxon.com/display/docs/formats_mdl-molfiles-rgfiles-sdfiles-rxnfiles-rdfiles-formats.md
            https://docs.chemaxon.com/display/docs/formats_chemaxon-specific-information-in-mdl-mol-files.md
^ CDK: https://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/io/MDLV2000Reader.html
+ Indigo: https://github.com/epam/Indigo/blob/master/core/indigo-core/molecule/src/molfile_loader.cpp

# Chemaxon specific information in SDT property field

|-----------------|---------------------------------------|---------------------------------------|-----------------------|
| Property        | Symbol                                | Additional properties                 | Interpretation        |
|-----------------|---------------------------------------|---------------------------------------|-----------------------|
| Coordinate Bond | M  SDT   1 MRV_COORDINATE_BOND_TYPE   | M  STY  3   1 DAT   2 DAT   3 DAT     | Must be a data SGroup |
|                 |                                       | M  SED   1 50                         | Bond index*           |
|                 |                                       | M  SAL   1  2  16  22                 | Atom indices*@        |
| Implicit H      | M  SDT   1 MRV_IMPLICIT_H             | M  STY  2   1 DAT   2 DAT             | Must be a data SGroup |
|                 |                                       | M  SAL   2  1  12                     | Atom index*           |
|                 |                                       | M  SED   2 IMPL_H1                    | Number of hydrogens^  |
| Multicenter Bd  | M  SDT   1 MRV_MULTICENTER_ATOM_INDEX | M  STY  1   1 DAT                     | Must be a data SGroup |
|                 |                                       | M  SED   1 19                         | Atom index*+          |
|                 |                                       | M  SAL   1  6  12  13  14  15  16  17 | Atom indices*+        |
| Charge SGroup   | M  SDT   1 MRV_CHARGE_ON_GROUP        | ?                                     | No examples           |
|-----------------|---------------------------------------|---------------------------------------|-----------------------|

Marvin: https://docs.chemaxon.com/display/docs/formats_chemaxon-specific-information-in-mdl-mol-files.md

Multipage document annotation properties do not seem too useful.

* 1-based indexing
@ donor atom first, acceptor atom second?
^ IMPL_H<n>, n denotes the number of hydrogens
+ SED property denotes pseudoatom multicenter attachment site, SAL contains the atom indices contributing to the pseudoatom
