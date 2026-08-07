import pytest
from umol import (
    AtomAst,
    BondAst,
    Element,
    Entity,
    MetadataError,
    MoleculeDefaults,
    MoleculeMetadata,
    ParseError,
    ReactionSpanAst,
    ValueAst,
)


def test_reaction_span_ast_parse():
    span = ReactionSpanAst.parse(r'{:atoms ["C" {:add "O"}]}')

    assert span == ReactionSpanAst.from_entries(
        [
            (AtomAst(Element("C")), AtomAst(Element("C"))),
            (None, AtomAst(Element("O"))),
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

    atom = AtomAst.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")
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
            (AtomAst(Element("C")), AtomAst(Element("C"))),
            (None, AtomAst(Element("O"))),
        ]
    )

    assert span.render() == r'{:atoms ["C" {:add "O"}]}'
    assert str(span) == span.render()


def test_reaction_span_ast_render_with_metadata_error():
    span = ReactionSpanAst.from_entries(
        [(AtomAst(Element("C")), AtomAst(Element("C")))]
    )
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(1), "outside")

    with pytest.raises(
        MetadataError,
        match="^metadata entity is out of range: atom 1$",
    ):
        span.render_with_metadata(metadata)


def test_reaction_span_ast_from_entries():
    lhs = AtomAst(Element("C"), charge=ValueAst.Lit(1))
    rhs = AtomAst(Element("C"), charge=ValueAst.LitSet({1}))

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
            [(AtomAst(Element("C")), AtomAst(Element("C")))],
            bonds=[(0, 1, (BondAst(1), BondAst(1)))],
        )
