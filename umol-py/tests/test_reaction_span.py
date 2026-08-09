import pytest
from umol import (
    AtomForm,
    BondForm,
    BondDelta,
    BondFieldChange,
    ContradictionError,
    Delta,
    Deltas,
    Element,
    Entity,
    MetadataError,
    MoleculeAst,
    MoleculeDefaults,
    MoleculeMetadata,
    ParseError,
    ReactionAst,
    ReactionSpanAst,
    NumForm,
)


def test_reaction_span_ast_parse():
    span = ReactionSpanAst.parse(r'{:atoms ["C" {:add "O"}]}')

    assert span == ReactionSpanAst.from_entries(
        [
            (AtomForm(Element("C")), AtomForm(Element("C"))),
            (None, AtomForm(Element("O"))),
        ]
    )


def test_reaction_span_ast_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        ReactionSpanAst.parse("not edn")


def test_reaction_span_ast_parse_defaults():
    span = ReactionSpanAst.parse(
        r'{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.ground(),
    )

    atom = AtomForm.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")
    assert span == ReactionSpanAst.from_entries([(atom, atom)])


def test_reaction_span_ast_parse_with_metadata_roundtrip():
    span, metadata = ReactionSpanAst.parse_with_metadata(
        r'{:atoms [[:carbon :x] {:add "O"}] '
        r':bonds [{:add [0 1 :single]}] :atom-aliases [:x "C"]}'
    )

    rendered = span.render_with_metadata(metadata)
    reparsed, reparsed_metadata = ReactionSpanAst.parse_with_metadata(rendered)

    assert rendered == (
        r'{:atom-aliases [:x "C"] :atoms [[:carbon :x] {:add "O"}] '
        r':bonds [{:add [:carbon 1 :single]}]}'
    )
    assert (reparsed, reparsed_metadata) == (span, metadata)


def test_reaction_span_ast_render():
    span = ReactionSpanAst.from_entries(
        [
            (AtomForm(Element("C")), AtomForm(Element("C"))),
            (None, AtomForm(Element("O"))),
        ]
    )

    assert span.render() == r'{:atoms ["C" {:add "O"}]}'
    assert str(span) == span.render()


def test_reaction_span_ast_render_with_metadata_error():
    span = ReactionSpanAst.from_entries(
        [(AtomForm(Element("C")), AtomForm(Element("C")))]
    )
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(1), "outside")

    with pytest.raises(
        MetadataError,
        match="^metadata entity is out of range: atom 1$",
    ):
        span.render_with_metadata(metadata)


def test_reaction_span_ast_from_entries():
    lhs = AtomForm(Element("C"), charge=NumForm.Lit(1))
    rhs = AtomForm(Element("C"), charge=NumForm.LitSet({1}))

    span = ReactionSpanAst.from_entries([(lhs, rhs)])

    assert span == ReactionSpanAst.from_entries([(lhs, lhs)])


def test_reaction_span_ast_from_entries_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entry is absent from both sides$",
    ):
        ReactionSpanAst.from_entries([(None, None)])


def test_reaction_span_ast_from_entries_reference_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entries reference unavailable atom 1$",
    ):
        ReactionSpanAst.from_entries(
            [(AtomForm(Element("C")), AtomForm(Element("C")))],
            bonds=[(0, 1, (BondForm(1), BondForm(1)))],
        )


def test_reaction_span_ast_lhs():
    span = ReactionSpanAst.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    assert span.lhs() == MoleculeAst.parse(
        r'{:atoms ["C" "O"] :bonds [[0 1 :single]]}'
    )


def test_reaction_span_ast_rhs():
    span = ReactionSpanAst.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    assert span.rhs() == MoleculeAst.parse(
        r'{:atoms ["C" "N"] :bonds [[0 1 :double]]}'
    )


def test_reaction_span_ast_correspondence():
    span = ReactionSpanAst.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    correspondence = span.correspondence()

    assert correspondence.atoms.matched_pairs == [(0, 0)]
    assert correspondence.atoms.left_count == 2
    assert correspondence.atoms.right_count == 2
    assert correspondence.atoms.left_unmatched == [1]
    assert correspondence.atoms.right_unmatched == [1]
    assert correspondence.bonds.matched_pairs == []
    assert correspondence.bonds.left_count == 1
    assert correspondence.bonds.right_count == 1
    assert correspondence.bonds.left_unmatched == [0]
    assert correspondence.bonds.right_unmatched == [0]


def test_reaction_span_ast_to_reaction_roundtrip():
    span = ReactionSpanAst.parse(
        r'{:atoms ["C" {:add "O"}] :bonds [{:add [0 1 :single]}]}'
    )

    reaction = span.to_reaction()

    assert reaction.lhs == span.lhs()
    assert reaction.to_reaction_span() == span


def test_reaction_ast_to_reaction_span_error():
    reaction = ReactionAst(
        MoleculeAst.from_entries(
            [AtomForm(Element("C")), AtomForm(Element("C"))],
            bonds=[(0, 1, BondForm(1))],
        ),
        Deltas(
            [
                Delta.Bond(
                    BondDelta.ModifyField(
                        id=0,
                        change=BondFieldChange.Order(
                            old=NumForm.Lit(2),
                            new=NumForm.Lit(3),
                        ),
                    )
                )
            ]
        ),
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        reaction.to_reaction_span()
