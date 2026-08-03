import pytest

from umol import (
    AtomAst,
    CisTransConfiguration,
    CisTransStereo,
    CisTransStereoAst,
    Element,
    MoleculeAst,
    Permutation,
    StereoAtomAst,
    StereoAtomConstraintAst,
    StereoAtomConstraintsAst,
    StereoAtomUpdate,
    StereoBondAst,
    StereoBondConstraintAst,
    StereoBondConstraintsAst,
    StereoBondUpdate,
    StereoConfigurationAst,
    StereoConfigurationUpdate,
    StereoCoset,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereoTerm,
    Stereogenicity,
    StereogenicityAst,
    TetrahedralConfiguration,
    TetrahedralStereo,
    TetrahedralStereoAst,
)


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (StereoConfigurationUpdate.Unchanged(), StereoConfigurationUpdate.Unchanged()),
        (StereoConfigurationUpdate.Undetermined(), StereoConfigurationUpdate.Undetermined()),
        (
            StereoConfigurationUpdate.Kinded(StereoKind.Tetrahedral, None),
            StereoConfigurationUpdate.Kinded(StereoKind.Tetrahedral, None),
        ),
        (
            StereoConfigurationUpdate.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
            StereoConfigurationUpdate.Kinded(
                StereoKind.Tetrahedral, StereoCoset.Undetermined()
            ),
        ),
        (
            StereoConfigurationUpdate.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(1)
            ),
            StereoConfigurationUpdate.Kinded(
                StereoKind.CisTrans, StereoCoset.Lit(1)
            ),
        ),
    ],
)
def test_stereo_configuration_update(update, expected):
    assert update == expected


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            StereoAtomUpdate(),
            (StereoConfigurationUpdate.Unchanged(), StereoAtomConstraintsAst([])),
        ),
        (
            StereoAtomUpdate(
                configuration=StereoConfigurationUpdate.Kinded(
                    StereoKind.Tetrahedral, None
                )
            ),
            (
                StereoConfigurationUpdate.Kinded(StereoKind.Tetrahedral, None),
                StereoAtomConstraintsAst([]),
            ),
        ),
        (
            StereoAtomUpdate(
                configuration=StereoConfigurationUpdate.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                )
            ),
            (
                StereoConfigurationUpdate.Kinded(
                    StereoKind.Tetrahedral, StereoCoset.Lit(1)
                ),
                StereoAtomConstraintsAst([]),
            ),
        ),
        (
            StereoAtomUpdate(configuration=StereoConfigurationUpdate.Undetermined()),
            (StereoConfigurationUpdate.Undetermined(), StereoAtomConstraintsAst([])),
        ),
        (
            StereoAtomUpdate(
                constraints=StereoAtomConstraintsAst(
                    [
                        StereoAtomConstraintAst.Stereogenicity(
                            StereogenicityAst.Undetermined()
                        )
                    ]
                )
            ),
            (
                StereoConfigurationUpdate.Unchanged(),
                StereoAtomConstraintsAst(
                    [
                        StereoAtomConstraintAst.Stereogenicity(
                            StereogenicityAst.Undetermined()
                        )
                    ]
                ),
            ),
        ),
    ],
)
def test_stereo_atom_update(update, expected):
    assert (update.configuration, update.constraints) == expected


@pytest.mark.parametrize(
    ("update", "expected"),
    [
        (
            StereoBondUpdate(),
            (StereoConfigurationUpdate.Unchanged(), StereoBondConstraintsAst([])),
        ),
        (
            StereoBondUpdate(
                configuration=StereoConfigurationUpdate.Kinded(
                    StereoKind.CisTrans, None
                )
            ),
            (
                StereoConfigurationUpdate.Kinded(StereoKind.CisTrans, None),
                StereoBondConstraintsAst([]),
            ),
        ),
        (
            StereoBondUpdate(
                configuration=StereoConfigurationUpdate.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                )
            ),
            (
                StereoConfigurationUpdate.Kinded(
                    StereoKind.CisTrans, StereoCoset.Lit(0)
                ),
                StereoBondConstraintsAst([]),
            ),
        ),
        (
            StereoBondUpdate(configuration=StereoConfigurationUpdate.Undetermined()),
            (StereoConfigurationUpdate.Undetermined(), StereoBondConstraintsAst([])),
        ),
        (
            StereoBondUpdate(
                constraints=StereoBondConstraintsAst(
                    [
                        StereoBondConstraintAst.Stereogenicity(
                            StereogenicityAst.Undetermined()
                        )
                    ]
                )
            ),
            (
                StereoConfigurationUpdate.Unchanged(),
                StereoBondConstraintsAst(
                    [
                        StereoBondConstraintAst.Stereogenicity(
                            StereogenicityAst.Undetermined()
                        )
                    ]
                ),
            ),
        ),
    ],
)
def test_stereo_bond_update(update, expected):
    assert (update.configuration, update.constraints) == expected


def test_permutation_image_degree():
    p = Permutation([1, 0, 2, 3])
    assert p.image() == [1, 0, 2, 3]
    assert p.degree == 4


@pytest.mark.parametrize(
    ("image", "message"),
    [
        ([0, 1, 2, 3, 4, 5, 6], "permutation image length 7 exceeds maximum 6"),
        ([0, 1, 3], "permutation image value 3 at position 2 is outside 0..3"),
        ([0, 1, 1], "permutation image value 1 occurs more than once"),
    ],
)
def test_permutation_error(image, message):
    with pytest.raises(ValueError) as error:
        Permutation(image)
    assert str(error.value) == message


def test_permutation_identity():
    assert Permutation.identity(4).image() == [0, 1, 2, 3]


def test_permutation_eq_hash():
    assert Permutation([1, 0, 2]) == Permutation([1, 0, 2])
    assert Permutation([1, 0, 2]) != Permutation([0, 1, 2])
    assert len({Permutation([1, 0, 2]), Permutation([1, 0, 2])}) == 1


def test_stereoterm_lit():
    assert StereoTerm.Lit(1) == StereoTerm.Lit(1)


def test_stereoterm_apply():
    term = StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))
    assert term == StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))


def test_stereocosetast_term():
    assert StereoCoset.Term(StereoTerm.Lit(2)) == StereoCoset.Term(StereoTerm.Lit(2))


def test_tetrahedralstereoast_stereo():
    assert TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)) == TetrahedralStereoAst.Stereo(
        StereoCoset.Lit(1)
    )


def test_tetrahedralstereoast_not_stereo():
    assert TetrahedralStereoAst.NotStereo() == TetrahedralStereoAst.NotStereo()


def test_tetrahedralstereoast_as_lit():
    assert TetrahedralStereoAst.NotStereo().as_lit() == TetrahedralStereo.NotStereo()
    assert TetrahedralStereoAst.Stereo(StereoCoset.Lit(1)).as_lit() == (
        TetrahedralStereo.Stereo(1)
    )
    assert TetrahedralStereoAst.Undetermined().as_lit() is None
    assert TetrahedralStereoAst.Stereo(StereoCoset.LitSet({0, 1})).as_lit() is None


def test_stereocosetast_eq_hash_repr():
    assert StereoCoset.Lit(1) == StereoCoset.Lit(1)
    assert StereoCoset.Lit(1) != StereoCoset.Lit(0)
    assert len({StereoCoset.Lit(1), StereoCoset.Lit(1)}) == 1
    assert repr(StereoCoset.Lit(1)) == "StereoCoset.Lit(1)"


def test_stereoterm_eq_repr():
    assert StereoTerm.Lit(0) == StereoTerm.Lit(0)
    assert StereoTerm.Lit(0) != StereoTerm.Swap(StereoTerm.Lit(0))
    assert repr(StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))) == (
        "StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))"
    )


def test_tetrahedralstereoast_stereo_repr():
    assert repr(TetrahedralStereoAst.Stereo(StereoCoset.Lit(1))) == (
        "TetrahedralStereoAst.Stereo(StereoCoset.Lit(1))"
    )


def test_cistransstereoast_stereo():
    assert CisTransStereoAst.Stereo(StereoCoset.Lit(1)) == CisTransStereoAst.Stereo(
        StereoCoset.Lit(1)
    )


def test_cistransstereoast_not_stereo():
    assert CisTransStereoAst.NotStereo() == CisTransStereoAst.NotStereo()


def test_cistransstereoast_as_lit():
    assert CisTransStereoAst.NotStereo().as_lit() == CisTransStereo.NotStereo()
    assert CisTransStereoAst.Stereo(StereoCoset.Lit(0)).as_lit() == CisTransStereo.Stereo(0)
    assert CisTransStereoAst.Undetermined().as_lit() is None
    assert CisTransStereoAst.Stereo(StereoCoset.LitSet({0, 1})).as_lit() is None


def test_cistransstereoast_stereo_repr():
    assert repr(CisTransStereoAst.Stereo(StereoCoset.Lit(1))) == (
        "CisTransStereoAst.Stereo(StereoCoset.Lit(1))"
    )


def test_cistransconfiguration_enum():
    assert CisTransConfiguration.Z == CisTransConfiguration.Z
    assert CisTransConfiguration.Z != CisTransConfiguration.E
    assert len(
        {CisTransConfiguration.E, CisTransConfiguration.E, CisTransConfiguration.Z}
    ) == 2


def stereo_atom_molecule():
    # a tetrahedral stereocenter on atom 0 with four atom ligands (atoms 1-4)
    return MoleculeAst.from_parts(
        [AtomAst(Element("C")) for _ in range(5)],
        stereo_atoms=[
            (
                0,
                [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
                StereoAtomAst(TetrahedralConfiguration.Ccw),
            )
        ],
    )


def test_stereoatomast_new():
    atom = StereoAtomAst(TetrahedralConfiguration.Ccw)
    assert str(atom) == "Th0"
    assert atom.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(0)
    )
    assert len(atom.constraints) == 0


@pytest.mark.parametrize("dsl", ["Th0", "Th*", "Sp2"])
def test_stereoatomast_parse_roundtrip(dsl):
    atom = StereoAtomAst.parse(dsl)
    assert str(atom) == dsl
    assert repr(atom) == f"StereoAtomAst.parse('{dsl}')"


def test_stereoatomast_configuration_setter():
    atom = StereoAtomAst(TetrahedralConfiguration.Ccw)
    atom.configuration = TetrahedralConfiguration.Cw
    assert str(atom) == "Th1"


def test_stereoatomast_constraints_kwarg():
    atom = StereoAtomAst(
        TetrahedralConfiguration.Ccw,
        constraints=StereoAtomConstraintsAst(
            [
                StereoAtomConstraintAst.Stereogenicity(
                    StereogenicityAst.Lit(Stereogenicity.Stereogenic)
                )
            ]
        ),
    )
    assert atom.constraints.stereogenicity() == StereogenicityAst.Lit(
        Stereogenicity.Stereogenic
    )


def test_stereoatomast_asdict():
    atom = StereoAtomAst(TetrahedralConfiguration.Ccw)
    d = atom.asdict()
    assert set(d.keys()) == {"configuration", "constraints"}
    assert d["configuration"] == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(0)
    )
    assert d["constraints"] == []


def test_stereobondast_parse_roundtrip():
    bond = StereoBondAst.parse("Ct0")
    assert str(bond) == "Ct0"
    assert bond.configuration == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCoset.Lit(0)
    )


def test_molecule_stereo_atoms_from_parts():
    views = stereo_atom_molecule().stereo_atoms
    assert len(views) == 1
    view = views[0]
    assert view.id == 0
    assert view.site_id == 0
    assert view.ligands == [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)]
    assert view.kind == StereoKind.Tetrahedral
    assert view.coset == StereoCoset.Lit(0)
    assert view.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCoset.Lit(0)
    )


def test_stereoatomview_configuration_write_through():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0].configuration = TetrahedralConfiguration.Cw
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.stereo_atoms[0].coset == StereoCoset.Lit(1)


def test_stereoatomview_constraints_write_through():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0].constraints.set(
        StereoAtomConstraintAst.Stereogenicity(
            StereogenicityAst.Lit(Stereogenicity.Stereogenic)
        )
    )
    assert mol.stereo_atoms[0].constraints.stereogenicity() == StereogenicityAst.Lit(
        Stereogenicity.Stereogenic
    )


def test_stereoatomviews_at():
    views = stereo_atom_molecule().stereo_atoms
    assert views.at(0).id == 0
    assert views.at(1) is None


def test_stereoatomviews_of():
    views = stereo_atom_molecule().stereo_atoms
    ligands = [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)]
    # order-independent full-ligand-set match
    assert views.of(0, list(reversed(ligands))).id == 0
    # a partial ligand set does not match
    assert views.of(0, ligands[:2]) is None


def test_stereoatomviews_setitem():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0] = StereoAtomAst.parse("Th1")
    view = mol.stereo_atoms[0]
    assert view.coset == StereoCoset.Lit(1)
    # site topology preserved
    assert view.site_id == 0


def test_molecule_repr_stereo():
    assert repr(stereo_atom_molecule()) == (
        "MoleculeAst(atoms=5, bonds=0, stereo_atoms=1)"
    )
