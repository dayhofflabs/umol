import pytest

from umol import (
    AtomAst,
    AtomDelta,
    AtomFieldChange,
    ContradictionError,
    Delta,
    Deltas,
    Element,
    MoleculeAst,
    ParseError,
    ReactionAst,
    ValueAst,
)


def test_reactionast_constructor():
    empty = ReactionAst()
    populated = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.Add(id=1, ast=AtomAst(Element("O")))
                )
            ]
        ),
    )

    assert empty.lhs == MoleculeAst()
    assert empty.deltas == Deltas()
    assert populated.lhs == MoleculeAst.from_parts([AtomAst(Element("C"))])
    assert populated.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )


def test_reactionast_constructor_snapshot():
    lhs = MoleculeAst.from_parts([AtomAst(Element("C"))])
    deltas = Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )
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
    deltas = Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))]
    )

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
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.Add(id=1, ast=AtomAst(Element("O")))
                )
            ]
        ),
    )
    expected = ReactionAst(reaction.lhs, reaction.deltas)

    reaction.lhs = reaction.lhs
    reaction.deltas = reaction.deltas

    assert reaction == expected


def test_reactionast_value():
    reaction = ReactionAst(
        MoleculeAst.from_parts([AtomAst(Element("C"))]),
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.Add(id=1, ast=AtomAst(Element("O")))
                )
            ]
        ),
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
    "text",
    [
        pytest.param(
            '{:lhs {:atoms ["C" "O"]} :deltas '
            '[{:atom {:add "N"}} {:atom {:remove 1}}]}',
            id="atom-add-remove",
        ),
        pytest.param(
            '{:lhs {:atoms ["Br#c0"]} :deltas '
            '[{:atom {:modify [0 "#c-1"]}}]}',
            id="atom-modify",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
            ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
            ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} '
            ':deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}',
            id="stereo-mirror",
        ),
        pytest.param(
            '{:lhs {:atoms ["C"]} :deltas '
            '[{:constraint {:add {:connected {}}}}]}',
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


def test_reactionast_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        ReactionAst.parse("not edn")


def test_reactionast_str_components():
    reaction = ReactionAst.parse('{:lhs {:atoms ["C"]} :deltas []}')

    reaction.lhs.atoms[0].charge = 1
    reaction.deltas.append(
        Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))
    )

    assert str(reaction) == (
        '{:deltas [{:atom {:add "O"}}] '
        ':lhs {:atoms ["C#c+"] :bonds []}}'
    )


def test_reactionast_parse_repr():
    reaction = ReactionAst.parse(
        '{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}'
    )

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
    canonical.deltas.append(
        Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst(Element("O"))))
    )
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
        '{:lhs {:atoms ["C" "O"]} :deltas '
        '[{:atom {:add "N"}} {:atom {:remove 1}}]}'
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
