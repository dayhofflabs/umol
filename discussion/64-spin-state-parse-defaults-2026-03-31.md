# Parsing modes for spin state data (unpaired electrons, multiplicity)

Parsing modes:

- unpaired: `enum UnpaireElectronsMode { Zero, Required, Derived }`
- multiplicity: `enum MultiplicityMode { Required, Derived }`

| Mode                     | Input    | Output                       | Ground Valid? | Pattern Valid? | Note                                                                   |
| ------------------------ | -------- | ---------------------------- | ------------- | -------------- | ---------------------------------------------------------------------- |
| u: Zero, m: Derived      | ""       | #u = 0, #s = 0 + 1 = singlet | Yes           | Yes            | Implies singlet, default for ground terms                              |
| "                        | "#u0"    | #u = 0, #s = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #s = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#s1"    | #u = 0, #s = singlet         | Yes           | Yes            |                                                                        |
| "                        | "#s2"    | #u = 0, #s = doublet         | No            | No             |                                                                        |
| "                        | "#u1#s1" | #u = 1, #s = singlet         | No            | No             |                                                                        |
| "                        | "#u1#s2" | #u = 1, #s = doublet         | Yes           | Yes            |                                                                        |
| u: Zero, m: Required     | ""       | #u = 0, #s = any             | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #s = any             | No            | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #s = any             | No            | Yes            |                                                                        |
| "                        | "#s1"    | #u = 0, #s = singlet         | Yes           | Yes            |                                                                        |
| "                        | "#s2"    | #u = 0, #s = doublet         | No            | No             |                                                                        |
| "                        | "#u1m1"  | #u = 1, #s = singlet         | No            | No             |                                                                        |
| "                        | "#u1m2"  | #u = 1, #s = doublet         | Yes           | Yes            |                                                                        |
| u: Required, m: Derived  | ""       | #u = any, #s = any           | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #s = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #s = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#s1"    | #u = any, #s = singlet       | No            | Yes            |                                                                        |
| "                        | "#s2"    | #u = any, #s = doublet       | No            | Yes            |                                                                        |
| "                        | "#u1#s1" | #u = 1, #s = singlet         | No            | No             |                                                                        |
| "                        | "#u1#s2" | #u = 1, #s = doublet         | Yes           | Yes            |                                                                        |
| u: Required, m: Required | ""       | #u = any, #s = any           | No            | Yes            |                                                                        |
| "                        | "#u0"    | #u = 0, #s = any             | No            | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #s = any             | No            | Yes            |                                                                        |
| "                        | "#s1"    | #u = any, #s = singlet       | No            | Yes            |                                                                        |
| "                        | "#s2"    | #u = any, #s = doublet       | No            | Yes            |                                                                        |
| "                        | "#u1#s1" | #u = 1, #s = singlet         | No            | No             |                                                                        |
| "                        | "#u1#s2" | #u = 1, #s = doublet         | Yes           | Yes            |                                                                        |
| u: Derived, m: Derived   | ""       | #u = any, #s = any           | No            | Yes            | Requires one of the two parameters to be defined, default for patterns |
| "                        | "#u0"    | #u = 0, #s = 0 + 1 = singlet | Yes           | Yes            |                                                                        |
| "                        | "#u1"    | #u = 1, #s = 1 + 1 = doublet | Yes           | Yes            |                                                                        |
| "                        | "#s1"    | #u = 1 - 1 = 0, #s = singlet | Yes           | Yes            |                                                                        |
| "                        | "#s2"    | #u = 2 - 1 = 1, #s = doublet | Yes           | Yes            |                                                                        |
| "                        | "#u1#s1" | #u = 1, #s = singlet         | No            | No             |                                                                        |
| "                        | "#u1#s2" | #u = 1, #s = doublet         | Yes           | Yes            |                                                                        |
| u: Derived, m: Required  |          |                              |               |                | Reverse of u: Required, m: Derived                                     |
