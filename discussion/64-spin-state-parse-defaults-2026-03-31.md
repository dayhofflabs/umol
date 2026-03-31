# Parsing modes for spin state data (unpaired electrons, multiplicity)

Parsing modes:

- unpaired: `enum UnpaireElectronsdMode { Zero, Required, Derived }`
- multiplicity: `enum MultiplicityMode { Required, Derived }`

| Mode                     | Input    | Output                       | Ground Valid? | Pattern Valid? | Note                                                                   |
| ------------------------ | -------- | ---------------------------- | ------------- | -------------- | ---------------------------------------------------------------------- |
| u: Zero, m: Derived      | ""       | #u = 0, #m = 0 + 1 = singlet | Yes           | Yes            | Implies singlet, default for ground terms                              |
| "                        | "#u0"    | #u = 0, #m = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #m = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#m1"    | #u = 0, #m = singlet         | Yes           | Yes            |                                                                        |
| "                        | "#m2"    | #u = 0, #m = doublet         | No            | No             |                                                                        |
| "                        | "#u1#m1" | #u = 1, #m = singlet         | No            | No             |                                                                        |
| "                        | "#u1#m2" | #u = 1, #m = doublet         | Yes           | Yes            |                                                                        |
| u: Zero, m: Required     | ""       | #u = 0, #m = any             | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #m = any             | No            | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #m = any             | No            | Yes            |                                                                        |
| "                        | "#m1"    | #u = 0, #m = singlet         | Yes           | Yes            |                                                                        |
| "                        | "#m2"    | #u = 0, #m = doublet         | No            | No             |                                                                        |
| "                        | "#u1m1"  | #u = 1, #m = singlet         | No            | No             |                                                                        |
| "                        | "#u1m2"  | #u = 1, #m = doublet         | Yes           | Yes            |                                                                        |
| u: Required, m: Derived  | ""       | #u = any, #m = any           | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #m = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #m = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#m1"    | #u = any, #m = singlet       | No            | Yes            |                                                                        |
| "                        | "#m2"    | #u = any, #m = doublet       | No            | Yes            |                                                                        |
| "                        | "#u1#m1" | #u = 1, #m = singlet         | No            | No             |                                                                        |
| "                        | "#u1#m2" | #u = 1, #m = doublet         | Yes           | Yes            |                                                                        |
| u: Required, m: Required | ""       | #u = any, #m = any           | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #m = any             | No            | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #m = any             | No            | Yes            |                                                                        |
| "                        | "#m1"    | #u = any, #m = singlet       | No            | Yes            |                                                                        |
| "                        | "#m2"    | #u = any, #m = doublet       | No            | Yes            |                                                                        |
| "                        | "#u1#m1" | #u = 1, #m = singlet         | No            | No             |                                                                        |
| "                        | "#u1#m2" | #u = 1, #m = doublet         | Yes           | Yes            |                                                                        |
| u: Derived, m: Derived   | ""       | #u = any, #m = any           | No            | Yes            | Requires one of the two parameters to be defined, default for patterns |
| "                        | "#u0"    | #u = 0, #m = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #m = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#m1"    | #u = 1 - 1 = 0, #m = singlet | Yes           | Yes            |                                                                        |
| "                        | "#m2"    | #u = 2 - 1 = 1, #m = doublet | Yes           | Yes            |                                                                        |
| "                        | "#u1#m1" | #u = 1, #m = singlet         | No            | No             |                                                                        |
| "                        | "#u1#m2" | #u = 1, #m = doublet         | Yes           | Yes            |                                                                        |
| u: Derived, m: Required  |          |                              |               |                | Reverse of u: Required, m: Derived                                     |
