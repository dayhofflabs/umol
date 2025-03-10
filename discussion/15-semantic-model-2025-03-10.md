# Semantic Model

The semantic model consists of three sets of entities, each forming a graph structure.

* __Structure Graph__. Nodes: chemical structures, edges: transformations (incl. reactions)
  - Expresses different states of the chemical system and their interconversions.
  - Examples of structures: different chemical entities, conformers of the same molecule, 
    collections of molecules, sets of reactants, sets of products.
  - Examples of transformations: changes in spatial structure (conformers), changes in bonding 
    (reactions, transition states), alchimical changes (library generation).
  - Conversions between structures are usually not lossy?
* __Model Graph__. Nodes: models, edges: conversions between models
  - Represents different views of the chemical system.
  - Examples of models: sum formula, molecular graph with discrete bonds, three-dimensional
    structure, ensemble of molecular graphs (Kekule resonance structures), hypergraph representation
    of aromaticity (Green model), ensemble of conformers, aggregates of molecules (flasks).
  - Models form an algebra analogous to algebraic data types: ensemble models corresponds to
    a collection of models of alternative states of the same system (like sum data type),
    aggregate model is a combination of several distinct models (like product data type).
  - Inputs and outputs to external formats proceed from specific models, require first a
    conversion to the relevant model.
  - Models are characterized by their capabilities: Molecular properties representable or
    computable within the model. For examples, see below.
  - Conversions between models are fundamentally lossy.
  - Conversions can be accomplished using different methods or algorithms.
* __Property Graph__. Nodes: properties, edges: computations
  - Expresses different molecular properties and the methods for computing properties
    from others.
  - Examples of properties: Graph model has atoms, bonds, atom charges (discrete),
    bond orders (usually discrete), molecular graph as capabilities. Quantum chemical
    model (SCF) has set of nuclei, locations of nuclei, # of electrons, atom charges
    (fractional), bond orders (fractional), orbital energies, orbital shapes. 
  - Computations between properties may fail.
  - Computations can be accomplished using different methods or algorithms.