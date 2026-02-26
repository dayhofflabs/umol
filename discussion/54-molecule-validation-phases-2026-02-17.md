# Resolution phases

1. Topology
2. Valence
3. Aromaticity
4. Stereo

# Molecule information by phase

| Item                      | Input | Topology | Valence | Aromaticity | Stereo | Notes                                                                                     | TODO |
| ------------------------- | ----- | -------- | ------- | ----------- | ------ | ----------------------------------------------------------------------------------------- | ---- |
| Atom::element             | +     | +        | +       | +           | +      | CTab, SMILES: atom type                                                                   |      |
| Atom::isotope_mass        | +     | +        | +       | +           | +      | CTab: isotope, SMILES: bracket                                                            |      |
| Atom::charge              | +     | +        | +       | +           | +      | CTab: charge, SMILES: bracket                                                             |      |
| Atom::hydrogens           | (+)   | (+)      | +       | +           | +      | CTab: H code, SMILES: bracket (optional)                                                  |      |
| Atom::lone_pairs          | (+)   | (+)      | +       | +           | +      | CX: lone pairs (CTab LP only in extended)                                                 |      |
| Atom::donated_pairs       | -     | -        | +       | +           | +      |                                                                                           |      |
| Atom::accepted_pairs      | -     | -        | +       | +           | +      |                                                                                           |      |
| Atom::unpaired_electrons  | (+)   | (+)      | +       | +           | +      | CTab, CX: radical                                                                         |      |
| Atom::multiplicity        | (+)   | (+)      | +       | +           | +      | CTab, CX: radical                                                                         |      |
| Atom::valence             | (+)   | (+)      | [+]     | +           | +      | CTab: valence                                                                             |      |
| Atom::aromatic_valence    | ~     | ~        | [+]     | +           | +      | CTab: bond type, SMILES: atom, bond type                                                  |      |
| Atom::multicenter_valence | -     | -        | +       | +           | +      | CX: only atoms, not number of electrons                                                   |      |
| Bond::order               | (+)   | +        | +       | +           | +      | CTab, SMILES: bond type (incl. aromatic); only localized (sigma) bond order               |      |
| Bond::donation            | (+)   | (+)      | [+]     | +           | +      | CX: coordinate bonds                                                                      |      |
| Graph                     | -     | +        | +       | +           | +      | Edges = Atoms, Nodes = Ordinary bonds                                                     |      |
| Atom, Bond indices        | -     | +        | +       | +           | +      | AtomIndex, BondIndex                                                                      |      |
| Total charge              | (+)   | (+)      | (+)     | (+)         | +      | CTab, SMILES: no explicit total charge                                                    |      |
| Aromatic systems          | ~     | ~        | ~       | +           | +      | CTab: bond type, SMILES: atom, bond type                                                  |      |
| Multicenter bonds         | (+)   | +        | +       | +           | +      | CX: multicenter bond (no #e)                                                              | +    |
| Noncovalent bonds         | (+)   | +        | +       | +           | +      | CX: hydrogen bond                                                                         | +    |
| Chiral center             | ~     | ~        | ~       | ~           | +      | CTab, SMILES: parity, CX ext: relative stereo, stereo group, ligand order, (local parity) | +    |
| Stereo bond               | ~     | ~        | ~       | ~           | +      | CTab, SMILES: bond wedge, CX: cis, trans, unspecified, wiggly                             | +    |
| Allene stereo             | ~     | ~        | ~       | ~           | +      | SMILES: parity                                                                            | +    |
| Bicyclo stereo            | ~     | ~        | ~       | ~           | +      | CX: (bicyclo stereo)                                                                      | +    |
| Positions                 | +     | +        | +       | +           | +      |                                                                                           | +    |

Notes: +: available, -: not available, (+): optional, [+]: multiple variants possible, ~: hints provided
TODO: data structures to be implemented
