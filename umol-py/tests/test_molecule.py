import re

import pytest

from umol import (
    AromaticSystemAst,
    AromaticValenceAst,
    AromaticityModel,
    AromaticityResolveConfig,
    AtomAst,
    AtomTypeRegistry,
    BondAst,
    ChemistryModel,
    ContradictionError,
    ElectronCountsAst,
    Element,
    ElementAst,
    ElementScope,
    InconsistencyPolicy,
    ModelConversionError,
    MoleculeAst,
    NoncovalentBondAst,
    NoncovalentBondKind,
    ParseError,
    ResolveConfig,
    RingLimits,
    SmilesIoConfig,
    SmilesSyntaxFlags,
    StereoCosetAst,
    StereoKind,
    StereoModel,
    StereoResolveConfig,
    UnderdeterminedError,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
    ValueAst,
)


def test_molecule_ast_new():
    assert len(MoleculeAst().atoms) == 0
    assert len(MoleculeAst().bonds) == 0


def test_molecule_ast_from_parts():
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("C"))],
        bonds=[(0, 1, BondAst(2))],
    )
    assert len(mol.atoms) == 2
    assert len(mol.bonds) == 1
    assert repr(mol) == "MoleculeAst(atoms=2, bonds=1)"


def test_molecule_ast_from_parts_default():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    assert len(mol.atoms) == 1
    assert len(mol.bonds) == 0


def test_molecule_ast_from_smiles():
    assert MoleculeAst.from_smiles("C") == MoleculeAst.from_parts(
        [AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
    )


@pytest.mark.parametrize(
    "io_config",
    [
        SmilesIoConfig.with_syntax_flags(
            syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS
        ),
        SmilesIoConfig.lenient(),
    ],
)
def test_molecule_ast_from_smiles_io_config(io_config):
    molecule = MoleculeAst.from_smiles("[se]1cccc1", io_config=io_config)

    assert [atom.element for atom in molecule.atoms] == [
        ElementAst.Lit(Element("Se")),
        *[ElementAst.Lit(Element("C")) for _ in range(4)],
    ]
    assert [
        (bond.atom_ids, bond.order) for bond in molecule.bonds
    ] == [
        ((0, 4), ValueAst.Lit(1)),
        ((0, 1), ValueAst.Lit(1)),
        ((1, 2), ValueAst.Lit(1)),
        ((2, 3), ValueAst.Lit(1)),
        ((3, 4), ValueAst.Lit(1)),
    ]
    assert [
        (system.atom_ids, system.electrons, system.charge)
        for system in molecule.aromatic_systems
    ] == [
        (
            (0, 1, 2, 3, 4),
            ElectronCountsAst.Lit([2, 1, 1, 1, 1]),
            ValueAst.Lit(0),
        )
    ]


def test_molecule_ast_from_smiles_io_config_error():
    with pytest.raises(ParseError, match="^Invalid token at position 2$"):
        MoleculeAst.from_smiles("C->N", io_config=SmilesIoConfig.opensmiles())
    with pytest.raises(
        ContradictionError,
        match="^no atom-typing match for AtomId\\(0\\) "
        "\\(element C, charge Some\\(0\\)\\)$",
    ):
        MoleculeAst.from_smiles("C->N", io_config=SmilesIoConfig.lenient())


@pytest.mark.parametrize(
    ("valence_model", "expected"),
    [
        (
            ValenceModel.AtomTyping(
                registry=AtomTypeRegistry.from_atoms(
                    [
                        AtomAst.parse(
                            "C#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"
                        )
                    ]
                )
            ),
            AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"),
        ),
        (
            ValenceModel.Counts(table=ValenceTable.default()),
            AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#a!"),
        ),
    ],
)
def test_molecule_ast_from_smiles_chemistry_model_valence(
    valence_model, expected
):
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=valence_model,
        aromaticity=default.aromaticity,
        stereo=default.stereo,
    )

    assert MoleculeAst.from_smiles(
        "C", chemistry_model=chemistry_model
    ) == MoleculeAst.from_parts([expected])


def test_molecule_ast_from_smiles_chemistry_model_aromaticity():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=default.valence,
        aromaticity=AromaticityModel.Hmo(
            scope=ElementScope.Any(), stabilization_threshold=0.375
        ),
        stereo=default.stereo,
    )

    molecule = MoleculeAst.from_smiles(
        "c1ccccc1", chemistry_model=chemistry_model
    )

    assert [atom.implicit_hydrogens for atom in molecule.atoms] == [
        ValueAst.Lit(1)
    ] * 6
    assert [
        atom.constraints.aromatic_valence for atom in molecule.atoms
    ] == [AromaticValenceAst.Aromatic(ValueAst.Lit(1))] * 6
    assert list(molecule.aromatic_systems) == []


def test_molecule_ast_from_smiles_chemistry_model_stereo():
    default = ChemistryModel.default()
    chemistry_model = ChemistryModel(
        valence=default.valence,
        aromaticity=default.aromaticity,
        stereo=StereoModel(
            kind_models={},
            para_stereo=False,
            max_iterations=16,
            inconsistency=InconsistencyPolicy.Strip,
        ),
    )

    molecule = MoleculeAst.from_smiles(
        "C[C@H](N)O", chemistry_model=chemistry_model
    )

    assert [atom.implicit_hydrogens for atom in molecule.atoms] == [
        ValueAst.Lit(3),
        ValueAst.Lit(1),
        ValueAst.Lit(2),
        ValueAst.Lit(1),
    ]
    assert [
        atom.constraints.tetrahedral_stereo for atom in molecule.atoms
    ] == [None] * 4
    assert list(molecule.stereo_atoms) == []


@pytest.mark.parametrize(
    ("source", "resolve_config", "expected"),
    [
        (
            "[cH+]1[cH][cH]1",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    delocalize_charge=False,
                    reset_aromatic_valence=False,
                ),
                stereo=StereoResolveConfig(),
            ),
            (
                [ValueAst.Lit(1), ValueAst.Lit(0), ValueAst.Lit(0)],
                [
                    AromaticValenceAst.Aromatic(ValueAst.Lit(0)),
                    AromaticValenceAst.Aromatic(ValueAst.Lit(1)),
                    AromaticValenceAst.Aromatic(ValueAst.Lit(1)),
                ],
                [None] * 3,
                [((0, 1, 2), ValueAst.Lit(0))],
                [],
            ),
        ),
        (
            "c1ccccc1",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(
                    delocalize_charge=True,
                    reset_aromatic_valence=True,
                ),
                stereo=StereoResolveConfig(),
            ),
            (
                [ValueAst.Lit(0)] * 6,
                [None] * 6,
                [None] * 6,
                [((0, 1, 2, 3, 4, 5), ValueAst.Lit(0))],
                [],
            ),
        ),
        (
            "C[C@H](N)O",
            ResolveConfig(
                aromaticity=AromaticityResolveConfig(),
                stereo=StereoResolveConfig(reset_stereo_constraints=True),
            ),
            (
                [ValueAst.Lit(0)] * 4,
                [AromaticValenceAst.NotAromatic()] * 4,
                [None] * 4,
                [],
                [(1, StereoKind.Tetrahedral, StereoCosetAst.Lit(0))],
            ),
        ),
    ],
)
def test_molecule_ast_from_smiles_resolve_config(
    source, resolve_config, expected
):
    molecule = MoleculeAst.from_smiles(source, resolve_config=resolve_config)

    assert (
        [atom.charge for atom in molecule.atoms],
        [atom.constraints.aromatic_valence for atom in molecule.atoms],
        [atom.constraints.tetrahedral_stereo for atom in molecule.atoms],
        [
            (system.atom_ids, system.charge)
            for system in molecule.aromatic_systems
        ],
        [
            (stereo.site_id, stereo.kind, stereo.coset)
            for stereo in molecule.stereo_atoms
        ],
    ) == expected


@pytest.mark.parametrize(
    ("source", "kwargs", "error_type", "message"),
    [
        (" C", {}, ParseError, "Leading whitespace"),
        (
            "C[S@]C",
            {},
            ModelConversionError,
            "tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands",
        ),
        (
            "[nH]1cccc1",
            {
                "chemistry_model": ChemistryModel(
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.Clar(
                        scope=ElementScope.Any(), ring_limits=RingLimits()
                    ),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "clar: non-benzenoid input: Clar model requires benzenoid input "
            "but non-carbon aromatic atoms are present",
        ),
        ("*", {}, UnderdeterminedError, "resolution underdetermined"),
        (
            "c1ccccc1",
            {
                "chemistry_model": ChemistryModel(
                    valence=ValenceModel.Counts(
                        table=ValenceTable(
                            entries={
                                Element("C"): ValenceEntry(
                                    target_covalences=[4],
                                    aromatic_valences=[0],
                                )
                            }
                        )
                    ),
                    aromaticity=AromaticityModel.Hmo(
                        scope=ElementScope.Any(),
                        stabilization_threshold=0.375,
                    ),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            RuntimeError,
            "hmo: missing parameters: no Van-Catledge parameters for C with "
            "0 pi-electrons",
        ),
    ],
)
def test_molecule_ast_from_smiles_error(
    source, kwargs, error_type, message
):
    with pytest.raises(error_type, match=f"^{re.escape(message)}$"):
        MoleculeAst.from_smiles(source, **kwargs)


def test_molecule_ast_from_smiles_keyword_error():
    with pytest.raises(TypeError):
        MoleculeAst.from_smiles("C", SmilesIoConfig.opensmiles())


def test_molecule_ast_from_smiles_ownership():
    io_config = SmilesIoConfig.opensmiles()
    chemistry_model = ChemistryModel.default()
    resolve_config = ResolveConfig.default()

    first = MoleculeAst.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    second = MoleculeAst.from_smiles(
        "C",
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )
    first.atoms[0].charge = 1

    assert second == MoleculeAst.from_parts(
        [AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")]
    )
    assert first != second
    assert io_config == SmilesIoConfig.opensmiles()
    assert chemistry_model == ChemistryModel.default()
    assert resolve_config == ResolveConfig.default()


def test_molecule_ast_bonds_error():
    with pytest.raises(IndexError):
        MoleculeAst.from_parts([AtomAst(Element("C"))]).bonds[0]


def test_molecule_ast_eq():
    assert MoleculeAst() == MoleculeAst()


@pytest.mark.parametrize(
    ("molecule", "expected"),
    [
        (MoleculeAst(), "MoleculeAst(atoms=0, bonds=0)"),
        (
            MoleculeAst.from_parts(
                [
                    AtomAst(Element("C")),
                    AtomAst(Element("C")),
                    AtomAst(Element("O")),
                ],
                bonds=[(0, 1, BondAst(1))],
                aromatic_systems=[
                    ([0, 1, 2], AromaticSystemAst([1, 1, 1]))
                ],
                noncovalent_bonds=[
                    (
                        [0, 2],
                        NoncovalentBondAst(
                            NoncovalentBondKind.HydrogenBond
                        ),
                    )
                ],
            ),
            "MoleculeAst(atoms=3, bonds=1, aromatic_systems=1, "
            "noncovalent_bonds=1)",
        ),
    ],
)
def test_molecule_ast_repr(molecule, expected):
    assert repr(molecule) == expected
