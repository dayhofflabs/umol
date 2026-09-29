import re

import pytest

from umol import (
    ChemistryModel,
    Correspondence,
    ModelConversionError,
    Molecule,
    Reaction,
)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("CCO", "CCO"),
        ("[CH3]", "[CH3]"),
        ("[13CH4]", "[13CH4]"),
        ("c1ccccc1", "c1ccccc1"),
        ("F[C@H](Cl)Br", "F[C@H](Cl)Br"),
        ("[H][C@](F)(Cl)Br", "[H][C@](F)(Cl)Br"),
        ("F/C=C/C=C/F", "F/C=C/C=C/F"),
    ],
)
def test_molecule_to_smiles(source, expected):
    assert Molecule.from_smiles(source).to_smiles() == expected


def test_molecule_to_smiles_chemistry_model():
    molecule = Molecule.from_smiles("CCO")
    assert molecule.to_smiles() == "CCO"
    assert (
        molecule.to_smiles(chemistry_model=ChemistryModel.default()) == "[CH3][CH2]O"
    )


@pytest.mark.parametrize(
    ("source", "message"),
    [
        (
            '{:atoms ["C#i=#c0#h256#n0#u0#s"]}',
            "Atom(AtomId(0)) implicit hydrogens cannot be represented in TableIR: "
            "Lit(256)",
        ),
        (
            '{:atoms ["C#i=#c0#h4#n0#u0#s" "C#i=#c0#h4#n0#u0#s"] '
            ':bonds [[0 1 "1#c+"]]}',
            "Bond(BondId(0)) requires concrete zero bond charge, found Lit(1)",
        ),
        (
            '{:atoms ["C#i=#c0#h10#n0#u0#s"]}',
            "unsupported atom 0 field implicit_hydrogens",
        ),
        ('{:atoms ["C#c0#h4#n0#u0#s"]}', "atom AtomId(0) has a non-ground isotope"),
    ],
)
def test_molecule_to_smiles_error(source, message):
    with pytest.raises(ModelConversionError, match=f"^{re.escape(message)}$"):
        Molecule.parse(source).to_smiles()


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (">>", ">>"),
        (">>C", ">>C"),
        ("[CH3:7][Cl:3]>>[CH3:7][OH:9]", "[CH3:1]Cl>>[CH3:1]O"),
        ("c1ccccc1>>c1ccccc1", "c1ccccc1>>c1ccccc1"),
        ("F[C@H](Cl)Br>>F[C@@H](Cl)Br", "F[C@H](Cl)Br>>F[C@@H](Cl)Br"),
        ("F/C=C/C=C/F>>F/C=C/C=C\\F", "F/C=C/C=C/F>>F/C=C/C=C\\F"),
    ],
)
def test_reaction_to_reaction_smiles(source, expected):
    assert Reaction.from_reaction_smiles(source).to_reaction_smiles() == expected


def test_reaction_to_reaction_smiles_from_sides():
    lhs = Molecule.from_smiles("CC(=O)O.OCC")
    rhs = Molecule.from_smiles("CC(=O)OCC")
    pairs = [(0, 0), (1, 1), (2, 2), (4, 3), (5, 4), (6, 5)]
    reaction = Reaction.from_sides(lhs, rhs, Correspondence(pairs, 7, 6))
    text = reaction.to_reaction_smiles()
    assert text.count(">>") == 1
    assert Reaction.from_reaction_smiles(text).to_reaction_smiles() == text


def test_reaction_to_reaction_smiles_error():
    atom = "C#i=#c0#h3#n0#u0#s"
    lhs = Molecule.parse(f'{{:atoms ["{atom}" "{atom}"] :bonds [[0 1 "1#c+"]]}}')
    reaction = Reaction.from_sides(lhs, lhs, Correspondence([(0, 0), (1, 1)], 2, 2))
    message = "reactants: Bond(BondId(0)) requires concrete zero bond charge, found Lit(1)"
    with pytest.raises(ModelConversionError, match=f"^{re.escape(message)}$"):
        reaction.to_reaction_smiles()
