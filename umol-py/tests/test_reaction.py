import pytest
import umol

from umol import (
    AtomAst,
    AtomDelta,
    AtomFieldChange,
    CommonSubgraphEnumerationAlgorithm,
    ContradictionError,
    Correspondence,
    Delta,
    Deltas,
    Element,
    InvalidStructureError,
    MoleculeAst,
    MoleculeCorrespondence,
    ParseError,
    ReactionApplicationConfig,
    ReactionAst,
    ReactionDefaults,
    ReactionDerivation,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
    ValueAst,
)


def test_application_exports():
    assert {
        "Correspondence",
        "MoleculeCorrespondence",
        "ReactionApplicationConfig",
        "ReactionDerivation",
    } <= set(umol.__all__)
    assert umol.Correspondence is Correspondence
    assert umol.MoleculeCorrespondence is MoleculeCorrespondence
    assert umol.ReactionApplicationConfig is ReactionApplicationConfig
    assert umol.ReactionDerivation is ReactionDerivation


def test_reactionapplicationconfig_default():
    config = ReactionApplicationConfig()

    assert config == ReactionApplicationConfig.default()
    assert config.match_algorithm == SubstructureMatchAlgorithm.GraphAndOverlays()
    assert (
        config.subgraph_isomorphism_algorithm
        == SubgraphIsomorphismAlgorithm.Vf2Rdkit()
    )
    assert repr(config) == (
        "ReactionApplicationConfig("
        "match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(), "
        "subgraph_isomorphism_algorithm="
        "SubgraphIsomorphismAlgorithm.Vf2Rdkit())"
    )


@pytest.mark.parametrize(
    "match_algorithm,subiso_algorithm,expected_repr",
    [
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.Vf2(),
            "SubgraphIsomorphismAlgorithm.Vf2()",
            id="vf2",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.Ullmann(),
            "SubgraphIsomorphismAlgorithm.Ullmann()",
            id="ullmann",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.Ri(),
            "SubgraphIsomorphismAlgorithm.Ri()",
            id="ri",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6),
            "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6)",
            id="arc-match",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
            "SubgraphIsomorphismAlgorithm.Vf2Rdkit()",
            id="vf2-rdkit",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.RayKirsch(),
            "SubgraphIsomorphismAlgorithm.RayKirsch()",
            id="ray-kirsch",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.Incidence(),
            SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
            "SubgraphIsomorphismAlgorithm.Vf2Rdkit()",
            id="incidence",
        ),
    ],
)
def test_reactionapplicationconfig_value(
    match_algorithm, subiso_algorithm, expected_repr
):
    config = ReactionApplicationConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subiso_algorithm,
    )
    equal = ReactionApplicationConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subiso_algorithm,
    )

    assert config == equal
    assert config.match_algorithm == match_algorithm
    assert config.subgraph_isomorphism_algorithm == subiso_algorithm
    assert repr(config) == (
        "ReactionApplicationConfig("
        f"match_algorithm={match_algorithm!r}, "
        f"subgraph_isomorphism_algorithm={expected_repr})"
    )


@pytest.mark.parametrize(
    "value_type,message",
    [
        pytest.param(
            Correspondence,
            "cannot create 'builtins.Correspondence' instances",
            id="correspondence",
        ),
        pytest.param(
            MoleculeCorrespondence,
            "cannot create 'builtins.MoleculeCorrespondence' instances",
            id="molecule-correspondence",
        ),
        pytest.param(
            ReactionDerivation,
            "cannot create 'builtins.ReactionDerivation' instances",
            id="reaction-derivation",
        ),
    ],
)
def test_return_only_value_constructor_error(value_type, message):
    with pytest.raises(TypeError, match=f"^{message}$"):
        value_type()


def test_reactionast_constructor():
    empty = ReactionAst()
    populated = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]),
    )

    assert empty.lhs == MoleculeAst()
    assert empty.deltas == Deltas()
    assert populated.lhs == MoleculeAst.from_parts([AtomAst(Element("C"))])
    assert populated.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )


def test_reactionast_constructor_snapshot():
    lhs = MoleculeAst.from_parts([AtomAst(Element("C"))])
    deltas = Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))])
    reaction = ReactionAst(lhs, deltas)

    lhs.atoms[0].charge = 1
    deltas.append(Delta.Atom(AtomDelta.Add(id=2, ast=AtomAst(Element("N")))))

    assert reaction.lhs.atoms[0].charge == ValueAst.Undetermined()
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )
    assert reaction.lhs is not lhs
    assert reaction.deltas is not deltas


def test_reactionast_components():
    reaction = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas(),
    )

    lhs = reaction.lhs
    deltas = reaction.deltas
    lhs.atoms[0].charge = -1
    deltas.append(Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O")))))

    assert reaction.lhs is lhs
    assert reaction.deltas is deltas
    assert reaction.lhs.atoms[0].charge == ValueAst.Lit(-1)
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )


def test_reactionast_component_replacement():
    reaction = ReactionAst()
    lhs = MoleculeAst.from_parts([AtomAst(Element("C"))])
    deltas = Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))])

    reaction.lhs = lhs
    reaction.deltas = deltas
    lhs.atoms[0].charge = 1
    deltas.append(Delta.Atom(AtomDelta.Add(id=2, ast=AtomAst(Element("N")))))

    assert reaction.lhs.atoms[0].charge == ValueAst.Undetermined()
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )
    assert reaction.lhs is not lhs
    assert reaction.deltas is not deltas


def test_reactionast_component_replacement_self():
    reaction = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]),
    )
    expected = ReactionAst(reaction.lhs, reaction.deltas)

    reaction.lhs = reaction.lhs
    reaction.deltas = reaction.deltas

    assert reaction == expected


def test_reactionast_value():
    reaction = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]),
    )

    assert reaction == ReactionAst(reaction.lhs, reaction.deltas)
    assert reaction != ReactionAst()
    assert repr(reaction) == (
        "ReactionAst(lhs=MoleculeAst(atoms=1, bonds=0), "
        "deltas=Deltas([Delta.Atom(AtomDelta.Add("
        "id=1, ast=AtomAst.parse('O')))]))"
    )
    with pytest.raises(TypeError):
        hash(reaction)


@pytest.mark.parametrize(
    ("value", "expected", "expected_repr"),
    [
        (ReactionDefaults(), ReactionDefaults(), "ReactionDefaults()"),
        (
            ReactionDefaults.ground(),
            ReactionDefaults.ground(),
            "ReactionDefaults.ground()",
        ),
    ],
)
def test_reactiondefaults_value(value, expected, expected_repr):
    assert value == expected
    assert repr(value) == expected_repr


@pytest.mark.parametrize(
    "text",
    [
        pytest.param(
            '{:lhs {:atoms ["C" "O"]} :deltas '
            '[{:atom {:add "N"}} {:atom {:remove 1}}]}',
            id="atom-add-remove",
        ),
        pytest.param(
            '{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}',
            id="atom-modify",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
            ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
            ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} '
            ":deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}",
            id="stereo-mirror",
        ),
        pytest.param(
            '{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}',
            id="molecule-constraint",
        ),
    ],
)
def test_reactionast_parse(text):
    first = ReactionAst.parse(text)

    canonical = str(first)
    second = ReactionAst.parse(canonical)

    assert second == first
    assert str(second) == canonical
    assert second.lhs is not first.lhs
    assert second.deltas is not first.deltas


def test_reactionast_parse_defaults():
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"]} '
        ':deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}]}',
        defaults=ReactionDefaults.ground(),
    )

    assert reaction == ReactionAst.parse(
        '{:lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]} '
        ":deltas [{:atom {:add "
        '"O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"}}]}'
    )


def test_reactionast_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        ReactionAst.parse("not edn")


def test_reactionast_parse_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^ReactionAst.parse\\(\\) takes 1 positional arguments but 2 were given$"
        ),
    ):
        ReactionAst.parse(
            '{:lhs {:atoms ["C"]} :deltas []}',
            ReactionDefaults.ground(),
        )


def test_reactionast_from_sides():
    lhs = MoleculeAst.from_parts([AtomAst(Element("C")), AtomAst(Element("O"))])
    rhs = MoleculeAst.from_parts([AtomAst(Element("C")), AtomAst(Element("N"))])
    lhs_snapshot = MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("O"))]
    )
    rhs_snapshot = MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("N"))]
    )

    reaction = ReactionAst.from_sides(
        lhs,
        rhs,
        ((left, right) for left, right in [(0, 0)]),
    )

    assert reaction == ReactionAst(
        lhs_snapshot,
        Deltas(
            [
                Delta.Atom(AtomDelta.Remove(id=1, ast=AtomAst(Element("O")))),
                Delta.Atom(AtomDelta.Add(id=2, ast=AtomAst(Element("N")))),
            ]
        ),
    )
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot
    assert reaction.lhs is not lhs


def test_reactionast_from_sides_snapshot():
    lhs = MoleculeAst.from_parts([AtomAst(Element("C")), AtomAst(Element("O"))])
    rhs = MoleculeAst.from_parts([AtomAst(Element("C")), AtomAst(Element("N"))])
    reaction = ReactionAst.from_sides(lhs, rhs, [(0, 0)])
    expected = ReactionAst(reaction.lhs, reaction.deltas)

    lhs.atoms[0].charge = 1
    rhs.atoms[0].charge = -1

    assert reaction == expected
    assert reaction.lhs is not lhs

    reaction.lhs.atoms[0].charge = 2
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=3, ast=AtomAst(Element("F")))))

    assert reaction.lhs.atoms[0].charge == ValueAst.Lit(2)
    assert reaction.deltas[-1] == Delta.Atom(
        AtomDelta.Add(id=3, ast=AtomAst(Element("F")))
    )
    assert lhs.atoms[0].charge == ValueAst.Lit(1)
    assert rhs.atoms[0].charge == ValueAst.Lit(-1)


@pytest.mark.parametrize(
    "lhs_count,rhs_count,atom_pairs,message",
    [
        pytest.param(
            2,
            2,
            [(0, 0), (0, 1)],
            "duplicate left atom id 0",
            id="duplicate-left",
        ),
        pytest.param(
            2,
            2,
            [(0, 1), (1, 1)],
            "duplicate right atom id 1",
            id="duplicate-right",
        ),
        pytest.param(
            2,
            1,
            [(2, 0)],
            "left atom id 2 out of range for 2 atoms",
            id="left-out-of-range",
        ),
        pytest.param(
            1,
            1,
            [(0, 1)],
            "right atom id 1 out of range for 1 atoms",
            id="right-out-of-range",
        ),
    ],
)
def test_reactionast_from_sides_error(lhs_count, rhs_count, atom_pairs, message):
    lhs = MoleculeAst.from_parts([AtomAst(Element("C")) for _ in range(lhs_count)])
    rhs = MoleculeAst.from_parts([AtomAst(Element("C")) for _ in range(rhs_count)])

    with pytest.raises(ValueError, match=f"^{message}$"):
        ReactionAst.from_sides(
            lhs,
            rhs,
            (pair for pair in atom_pairs),
        )


def test_reactionast_str_components():
    reaction = ReactionAst.parse('{:lhs {:atoms ["C"]} :deltas []}')

    reaction.lhs.atoms[0].charge = 1
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O")))))

    assert str(reaction) == (
        '{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C#c+"] :bonds []}}'
    )


def test_reactionast_parse_repr():
    reaction = ReactionAst.parse('{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}')

    assert repr(reaction) == (
        "ReactionAst(lhs=MoleculeAst(atoms=1, bonds=0), "
        "deltas=Deltas([Delta.Atom(AtomDelta.Add("
        "id=1, ast=AtomAst.parse('O')))]))"
    )


def test_reactionast_canonicalize():
    source = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"), charge=0)]),
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.ModifyField(
                        id=0,
                        change=AtomFieldChange.Charge(
                            old=ValueAst.Lit(0), new=ValueAst.Lit(1)
                        ),
                    )
                ),
                Delta.Atom(
                    AtomDelta.ModifyField(
                        id=0,
                        change=AtomFieldChange.Charge(
                            old=ValueAst.Lit(1), new=ValueAst.Lit(2)
                        ),
                    )
                ),
            ]
        ),
    )
    snapshot = ReactionAst(source.lhs, source.deltas)

    canonical = source.canonicalize()

    assert canonical.deltas == Deltas(
        [
            Delta.Atom(
                AtomDelta.ModifyField(
                    id=0,
                    change=AtomFieldChange.Charge(
                        old=ValueAst.Lit(0), new=ValueAst.Lit(2)
                    ),
                )
            )
        ]
    )
    assert canonical.canonicalize() == canonical
    assert source == snapshot
    assert canonical.lhs is not source.lhs
    assert canonical.deltas is not source.deltas

    canonical.lhs.atoms[0].charge = 3
    canonical.deltas.append(Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O")))))
    assert canonical.lhs.atoms[0].charge == ValueAst.Lit(3)
    assert len(canonical.deltas) == 2


def test_reactionast_canonicalize_error():
    source = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"), charge=0)]),
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.ModifyField(
                        id=0,
                        change=AtomFieldChange.Charge(
                            old=ValueAst.Lit(0), new=ValueAst.Lit(1)
                        ),
                    )
                ),
                Delta.Atom(
                    AtomDelta.ModifyField(
                        id=0,
                        change=AtomFieldChange.Charge(
                            old=ValueAst.Lit(2), new=ValueAst.Lit(3)
                        ),
                    )
                ),
            ]
        ),
    )
    snapshot = ReactionAst(source.lhs, source.deltas)

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        source.canonicalize()

    assert source == snapshot


def test_reactionast_reverse():
    source = ReactionAst.parse(
        '{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}'
    )
    snapshot = ReactionAst(source.lhs, source.deltas)

    reversed_reaction = source.reverse()
    roundtrip = reversed_reaction.reverse()

    assert reversed_reaction.lhs == MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("N"))]
    )
    assert roundtrip.canonicalize() == source.canonicalize()
    assert source == snapshot
    assert reversed_reaction.lhs is not source.lhs
    assert reversed_reaction.deltas is not source.deltas

    reversed_reaction.lhs.atoms[0].charge = 1
    reversed_reaction.deltas.append(
        Delta.Atom(AtomDelta.Add(id=2, ast=AtomAst(Element("F"))))
    )
    assert reversed_reaction.lhs.atoms[0].charge == ValueAst.Lit(1)
    assert len(reversed_reaction.deltas) == 3


@pytest.mark.parametrize(
    "algorithm",
    [
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
            id="modular-product-backtracking",
        ),
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
            id="direct-backtracking",
        ),
    ],
)
@pytest.mark.parametrize(
    "first,second,expected",
    [
        pytest.param(
            '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}',
            '{:lhs {:atoms ["N#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}',
            [
                '{:lhs {:atoms ["C#c0" "N#c0"]} :deltas '
                '[{:atom {:modify [0 "#c+"]}} '
                '{:atom {:modify [1 "#c+"]}}]}'
            ],
            id="no-match",
        ),
        pytest.param(
            '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}',
            '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}',
            [
                '{:lhs {:atoms ["C#c0" "C#c+"]} :deltas '
                '[{:atom {:modify [0 "#c+"]}} '
                '{:atom {:modify [1 "#c+2"]}}]}',
                '{:lhs {:atoms ["C#c0"]} '
                ':deltas [{:atom {:modify [0 "#c+2"]}}]}',
            ],
            id="admissible",
        ),
        pytest.param(
            '{:lhs {:atoms ["C"]} :deltas [{:atom {:remove 0}}]}',
            '{:lhs {:atoms ["N#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}',
            [
                '{:lhs {:atoms ["N#c0" "C"]} :deltas '
                '[{:atom {:remove 1}} '
                '{:atom {:modify [0 "#c+"]}}]}'
            ],
            id="deletion-only",
        ),
    ],
)
def test_reactionast_compose(algorithm, first, second, expected):
    first = ReactionAst.parse(first)
    second = ReactionAst.parse(second)

    composites = first.compose(
        second,
        algorithm=algorithm,
    )

    assert composites == [ReactionAst.parse(reaction) for reaction in expected]


def test_reactionast_compose_default():
    first = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second = ReactionAst.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    fused = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    disjoint = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0" "C#c+"]} :deltas '
        '[{:atom {:modify [0 "#c+"]}} '
        '{:atom {:modify [1 "#c+2"]}}]}'
    )

    omitted = first.compose(second)
    explicit = first.compose(
        second,
        algorithm=CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
    )

    assert omitted == [disjoint, fused]
    assert explicit == [disjoint, fused]
    with pytest.raises(TypeError):
        first.compose(
            second,
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
        )


def test_reactionast_compose_snapshot():
    first = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second = ReactionAst.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    first_snapshot = ReactionAst(first.lhs, first.deltas)
    second_snapshot = ReactionAst(second.lhs, second.deltas)

    first.compose(first)
    composites = first.compose(second)

    assert first == first_snapshot
    assert second == second_snapshot
    assert len(composites) == 2
    assert composites[0].lhs is not first.lhs
    assert composites[0].lhs is not second.lhs
    assert composites[0].deltas is not first.deltas
    assert composites[0].deltas is not second.deltas
    assert composites[0].lhs is not composites[1].lhs
    assert composites[0].deltas is not composites[1].deltas

    for composite in composites:
        composite.lhs.atoms[0].charge = 7
        composite.deltas.append(
            Delta.Atom(AtomDelta.Add(id=8, ast=AtomAst(Element("Cl"))))
        )

        assert composite.lhs.atoms[0].charge == ValueAst.Lit(7)
        assert composite.deltas[-1] == Delta.Atom(
            AtomDelta.Add(id=8, ast=AtomAst(Element("Cl")))
        )

    assert first == first_snapshot
    assert second == second_snapshot


def test_reactionast_apply():
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    host = ReactionAst.parse('{:lhs {:atoms ["C#c0" "C#c0"]} :deltas []}').lhs
    reaction_snapshot = ReactionAst(reaction.lhs, reaction.deltas)
    host_snapshot = ReactionAst(host).lhs
    first_product = ReactionAst.parse('{:lhs {:atoms ["C#c+" "C#c0"]} :deltas []}').lhs
    second_product = ReactionAst.parse('{:lhs {:atoms ["C#c0" "C#c+"]} :deltas []}').lhs

    application = reaction.apply(host)

    assert iter(application) is application
    assert reaction == reaction_snapshot
    assert host == host_snapshot

    reaction.lhs.atoms[0].charge = 7
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("N")))))
    host.atoms[0].charge = 8

    first = next(application)
    remaining = list(application)

    assert len(remaining) == 1
    second = remaining[0]
    assert first.lhs == host_snapshot
    assert first.rhs == first_product
    assert second.lhs == host_snapshot
    assert second.rhs == second_product
    with pytest.raises(StopIteration):
        next(application)
    with pytest.raises(StopIteration):
        next(application)

    comap = first.comap
    atom_map = first.atom_map

    assert first.comap == comap
    assert first.comap is not comap
    assert first.atom_map == atom_map
    assert first.atom_map is not atom_map
    assert atom_map == comap.atoms
    assert comap.atoms.mates == [(0, 0), (1, 1)]
    assert comap.atoms.left_count == 2
    assert comap.atoms.right_count == 2
    assert comap.atoms.left_exposed == []
    assert comap.atoms.right_exposed == []
    for entity_map in (
        comap.bonds,
        comap.dative_bonds,
        comap.aromatic_systems,
        comap.multicenter_bonds,
        comap.noncovalent_bonds,
        comap.stereo_atoms,
        comap.stereo_bonds,
    ):
        assert entity_map.mates == []
        assert entity_map.left_count == 0
        assert entity_map.right_count == 0
        assert entity_map.left_exposed == []
        assert entity_map.right_exposed == []

    detached_lhs = first.lhs
    detached_rhs = first.rhs
    detached_lhs.atoms[0].charge = 5
    detached_rhs.atoms[0].charge = 6

    assert first.lhs == host_snapshot
    assert first.rhs == first_product
    assert second.lhs == host_snapshot
    assert second.rhs == second_product

    reversed_first = first.reverse()

    assert reversed_first.lhs == first_product
    assert reversed_first.rhs == host_snapshot
    assert reversed_first.reverse() == first

    second_step = ReactionAst.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    config = ReactionApplicationConfig(
        match_algorithm=SubstructureMatchAlgorithm.Incidence(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Vf2(),
    )
    following = next(second_step.apply(first.rhs, config=config))
    chained = first.chain(following)
    chained_product = ReactionAst.parse(
        '{:lhs {:atoms ["C#c+2" "C#c0"]} :deltas []}'
    ).lhs

    assert chained.lhs == host_snapshot
    assert chained.rhs == chained_product
    assert first.rhs == first_product
    assert following.lhs == first_product

    expected_reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0" "C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    recovered = chained.to_reaction()
    independent = chained.to_reaction()

    assert recovered == expected_reaction
    assert independent == expected_reaction
    assert recovered.lhs is not independent.lhs
    assert recovered.deltas is not independent.deltas

    recovered.lhs.atoms[1].charge = 9
    recovered.deltas.append(Delta.Atom(AtomDelta.Add(id=2, ast=AtomAst(Element("O")))))

    assert independent == expected_reaction
    assert chained.lhs == host_snapshot
    assert chained.rhs == chained_product

    zero = ReactionAst.parse('{:lhs {:atoms ["N"]} :deltas []}').apply(host_snapshot)

    assert iter(zero) is zero
    assert list(zero) == []
    with pytest.raises(StopIteration):
        next(zero)


@pytest.mark.parametrize(
    "config",
    [
        pytest.param(None, id="default"),
        pytest.param(
            ReactionApplicationConfig(
                match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(),
                subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Ullmann(),
            ),
            id="graph-and-overlays-ullmann",
        ),
        pytest.param(
            ReactionApplicationConfig(
                match_algorithm=SubstructureMatchAlgorithm.Incidence(),
                subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Vf2(),
            ),
            id="incidence-vf2",
        ),
    ],
)
def test_reaction_ast_apply_config(config):
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    host = ReactionAst.parse('{:lhs {:atoms ["C#c0" "C#c0"]} :deltas []}').lhs
    expected = [
        ReactionAst.parse('{:lhs {:atoms ["C#c+" "C#c0"]} :deltas []}').lhs,
        ReactionAst.parse('{:lhs {:atoms ["C#c0" "C#c+"]} :deltas []}').lhs,
    ]

    application = (
        reaction.apply(host) if config is None else reaction.apply(host, config=config)
    )

    assert [derivation.rhs for derivation in application] == expected


def test_reaction_ast_apply_config_error():
    with pytest.raises(TypeError):
        ReactionAst().apply(MoleculeAst(), ReactionApplicationConfig())


def test_reaction_ast_apply_rejection():
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C"]} :deltas [{:atom {:remove 0}}]}'
    )
    host = ReactionAst.parse(
        '{:lhs {:atoms ["C" "C" "C" "O"] :bonds [[1 3 "1"]]} :deltas []}'
    ).lhs
    expected = [
        ReactionAst.parse(
            '{:lhs {:atoms ["C" "C" "O"] :bonds [[0 2 "1"]]} :deltas []}'
        ).lhs,
        ReactionAst.parse(
            '{:lhs {:atoms ["C" "C" "O"] :bonds [[1 2 "1"]]} :deltas []}'
        ).lhs,
    ]

    assert [derivation.rhs for derivation in reaction.apply(host)] == expected


def test_reaction_ast_apply_precondition_error():
    reaction = ReactionAst()
    host = ReactionAst.parse(
        '{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"] [0 1 "2"]]} :deltas []}'
    ).lhs

    with pytest.raises(
        InvalidStructureError,
        match=(
            r"^invalid host: bond: parallel bonds on atoms "
            r"\[AtomId\(0\), AtomId\(1\)\]$"
        ),
    ):
        reaction.apply(host)


def test_reaction_ast_apply_iteration_error():
    reaction = ReactionAst.parse(
        "{:lhs {:atoms [\"C\"] "
        ":constraints [{:charge-sum {:atoms [0] :sum 0}}]} "
        ":deltas [{:constraint {:remove {:charge-sum {:atoms [0] :sum 0}}}}]}"
    )
    host = ReactionAst.parse('{:lhs {:atoms ["C"]} :deltas []}').lhs

    application = reaction.apply(host)

    with pytest.raises(
        RuntimeError,
        match=r"^apply transaction failed: missing constraint entry on remove$",
    ):
        next(application)
    with pytest.raises(StopIteration):
        next(application)
    with pytest.raises(StopIteration):
        next(application)


def test_reactionast_workflow():
    lhs = ReactionAst.parse(
        '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]} '
        ":deltas []}"
    ).lhs
    rhs = ReactionAst.parse(
        '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
        ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th0"}] '
        ":constraints [{:connected {}}]} :deltas []}"
    ).lhs
    lhs_snapshot = ReactionAst(lhs).lhs
    rhs_snapshot = ReactionAst(rhs).lhs
    expected_forward = ReactionAst.parse(
        "{:deltas ["
        "{:stereo-atom {:add {:ligands [1 2 3 4] :site 0 :type :ccw}}} "
        "{:constraint {:add {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]]}}"
    )

    reaction = ReactionAst.from_sides(
        lhs,
        rhs,
        ((atom_id, atom_id) for atom_id in range(5)),
    )

    assert reaction == expected_forward
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot
    assert reaction.lhs is not lhs

    normalized = reaction.canonicalize()

    assert normalized == expected_forward
    assert reaction == expected_forward
    assert normalized.lhs is not reaction.lhs
    assert normalized.deltas is not reaction.deltas

    reaction.lhs.atoms[0].charge = 1
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=5, ast=AtomAst(Element("Xe")))))

    assert normalized == expected_forward
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot

    rendered = str(normalized)

    assert rendered == (
        "{:deltas [{:stereo-atom {:add {:ligands [1 2 3 4] "
        ":site 0 :type :ccw}}} {:constraint {:add {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]]}}"
    )
    assert normalized == expected_forward

    parsed = ReactionAst.parse(rendered)

    assert parsed == expected_forward
    assert normalized == expected_forward
    assert parsed.lhs is not normalized.lhs
    assert parsed.deltas is not normalized.deltas

    normalized.lhs.atoms[0].charge = 2
    normalized.deltas.append(
        Delta.Atom(AtomDelta.Add(id=5, ast=AtomAst(Element("Ne"))))
    )

    assert parsed == expected_forward

    expected_reverse = ReactionAst.parse(
        "{:deltas ["
        "{:stereo-atom {:remove 0}} "
        "{:constraint {:remove {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]] :constraints [{:connected {}}] "
        ":stereo-atoms [{:ligands [1 2 3 4] :site 0 :type :ccw}]}}"
    )
    reversed_reaction = parsed.reverse()

    assert reversed_reaction == expected_reverse
    assert parsed == expected_forward
    assert reversed_reaction.lhs is not parsed.lhs
    assert reversed_reaction.deltas is not parsed.deltas

    parsed.lhs.atoms[0].charge = 3
    parsed.deltas.append(Delta.Atom(AtomDelta.Add(id=5, ast=AtomAst(Element("Ar")))))

    assert reversed_reaction == expected_reverse

    second = ReactionAst.parse(
        '{:lhs {:atoms ["Xe#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second_snapshot = ReactionAst(second.lhs, second.deltas)
    expected_composite = ReactionAst.parse(
        "{:deltas ["
        '{:atom {:modify [5 "#c+"]}} '
        "{:stereo-atom {:remove 0}} "
        "{:constraint {:remove {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I" "Xe#c0"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]] :constraints [{:connected {}}] "
        ":stereo-atoms [{:ligands [1 2 3 4] :site 0 :type :ccw}]}}"
    )

    composites = reversed_reaction.compose(second)

    assert composites == [expected_composite]
    assert reversed_reaction == expected_reverse
    assert second == second_snapshot
    assert composites[0].lhs is not reversed_reaction.lhs
    assert composites[0].lhs is not second.lhs
    assert composites[0].deltas is not reversed_reaction.deltas
    assert composites[0].deltas is not second.deltas

    reversed_reaction.lhs.atoms[0].charge = 4
    reversed_reaction.deltas.append(
        Delta.Atom(AtomDelta.Add(id=5, ast=AtomAst(Element("Kr"))))
    )

    assert composites == [expected_composite]
    assert second == second_snapshot

    composite_lhs = composites[0].lhs
    composite_deltas = composites[0].deltas
    composite_lhs.atoms[5].charge = 2
    composite_deltas.append(Delta.Atom(AtomDelta.Add(id=6, ast=AtomAst(Element("O")))))

    assert composites[0].lhs is composite_lhs
    assert composites[0].deltas is composite_deltas
    assert composites[0].lhs.atoms[5].charge == ValueAst.Lit(2)
    assert composites[0].deltas[-1] == Delta.Atom(
        AtomDelta.Add(id=6, ast=AtomAst(Element("O")))
    )
    assert second == second_snapshot
