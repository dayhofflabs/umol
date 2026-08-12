import re

import pytest

from umol import (
    AromaticityFailurePolicy,
    AromaticityModel,
    AromaticityResolveConfig,
    AtomForm,
    AtomDelta,
    AtomFieldChange,
    ChemistryModel,
    CommonSubgraphEnumerationAlgorithm,
    Correspondence,
    ContradictionError,
    Delta,
    Deltas,
    Element,
    ElementScope,
    Entity,
    MetadataError,
    ModelConversionError,
    Molecule,
    MoleculeCorrespondence,
    ParseError,
    ReactionApplicationConfig,
    Reaction,
    ReactionCompositionConfig,
    ReactionDefaults,
    ReactionDerivation,
    ReactionMetadata,
    RelevantCycleEnumerationAlgorithm,
    ResolveConfig,
    RingLimits,
    SmilesIoConfig,
    StereoResolveConfig,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
    TransactionError,
    UnderdeterminedError,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
    NumForm,
)


def test_reaction_composition_config_default():
    config = ReactionCompositionConfig()

    assert config == ReactionCompositionConfig.default()
    assert (
        config.common_subgraph_enumeration_algorithm
        == CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
    )
    assert repr(config) == (
        "ReactionCompositionConfig("
        "common_subgraph_enumeration_algorithm="
        "CommonSubgraphEnumerationAlgorithm.DirectBacktracking())"
    )


@pytest.mark.parametrize(
    "algorithm",
    [
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
            id="direct-backtracking",
        ),
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
            id="modular-product-backtracking",
        ),
    ],
)
def test_reaction_composition_config_value(algorithm):
    config = ReactionCompositionConfig(
        common_subgraph_enumeration_algorithm=algorithm,
    )
    equal = ReactionCompositionConfig(
        common_subgraph_enumeration_algorithm=algorithm,
    )

    assert config == equal
    assert config.common_subgraph_enumeration_algorithm == algorithm
    assert repr(config) == (
        "ReactionCompositionConfig("
        f"common_subgraph_enumeration_algorithm={algorithm!r})"
    )
    assert config != ReactionCompositionConfig(
        common_subgraph_enumeration_algorithm=(
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking()
            if algorithm
            == CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
            else CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
        )
    )
    with pytest.raises(AttributeError):
        config.common_subgraph_enumeration_algorithm = algorithm


def test_reaction_application_config_default():
    config = ReactionApplicationConfig()

    assert config == ReactionApplicationConfig.default()
    assert config.match_algorithm == SubstructureMatchAlgorithm.GraphAndOverlays()
    assert (
        config.subgraph_isomorphism_algorithm
        == SubgraphIsomorphismAlgorithm.Vf2Rdkit()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert repr(config) == (
        "ReactionApplicationConfig("
        "match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(), "
        "subgraph_isomorphism_algorithm="
        "SubgraphIsomorphismAlgorithm.Vf2Rdkit(), "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara())"
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
def test_reaction_application_config_value(
    match_algorithm, subiso_algorithm, expected_repr
):
    config = ReactionApplicationConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subiso_algorithm,
        relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara(),
    )
    equal = ReactionApplicationConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subiso_algorithm,
        relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara(),
    )

    assert config == equal
    assert config.match_algorithm == match_algorithm
    assert config.subgraph_isomorphism_algorithm == subiso_algorithm
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert repr(config) == (
        "ReactionApplicationConfig("
        f"match_algorithm={match_algorithm!r}, "
        f"subgraph_isomorphism_algorithm={expected_repr}, "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara())"
    )


@pytest.mark.parametrize(
    "value_type,message",
    [
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


def test_reaction_constructor():
    empty = Reaction()
    populated = Reaction(
        Molecule.from_entries([AtomForm(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]),
    )

    assert empty.lhs == Molecule()
    assert empty.deltas == Deltas()
    assert populated.lhs == Molecule.from_entries([AtomForm(Element("C"))])
    assert populated.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]
    )


def test_reaction_constructor_snapshot():
    lhs = Molecule.from_entries([AtomForm(Element("C"))])
    deltas = Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))])
    reaction = Reaction(lhs, deltas)

    lhs.atoms[0].charge = 1
    deltas.append(Delta.Atom(AtomDelta.Add(id=2, attributes=AtomForm(Element("N")))))

    assert reaction.lhs.atoms[0].charge == NumForm.Undetermined()
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]
    )
    assert reaction.lhs is not lhs
    assert reaction.deltas is not deltas


def test_reaction_components():
    reaction = Reaction(
        Molecule.from_entries([AtomForm(Element("C"))]),
        Deltas(),
    )

    lhs = reaction.lhs
    deltas = reaction.deltas
    lhs.atoms[0].charge = -1
    deltas.append(Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O")))))

    assert reaction.lhs is lhs
    assert reaction.deltas is deltas
    assert reaction.lhs.atoms[0].charge == NumForm.Lit(-1)
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]
    )


def test_reaction_component_replacement():
    reaction = Reaction()
    lhs = Molecule.from_entries([AtomForm(Element("C"))])
    deltas = Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))])

    reaction.lhs = lhs
    reaction.deltas = deltas
    lhs.atoms[0].charge = 1
    deltas.append(Delta.Atom(AtomDelta.Add(id=2, attributes=AtomForm(Element("N")))))

    assert reaction.lhs.atoms[0].charge == NumForm.Undetermined()
    assert reaction.deltas == Deltas(
        [Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]
    )
    assert reaction.lhs is not lhs
    assert reaction.deltas is not deltas


def test_reaction_component_replacement_self():
    reaction = Reaction(
        Molecule.from_entries([AtomForm(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]),
    )
    expected = Reaction(reaction.lhs, reaction.deltas)

    reaction.lhs = reaction.lhs
    reaction.deltas = reaction.deltas

    assert reaction == expected


def test_reaction_value():
    reaction = Reaction(
        Molecule.from_entries([AtomForm(Element("C"))]),
        Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O"))))]),
    )

    assert reaction == Reaction(reaction.lhs, reaction.deltas)
    assert reaction != Reaction()
    assert repr(reaction) == (
        "Reaction(lhs=Molecule(atoms=1, bonds=0), "
        "deltas=Deltas([Delta.Atom(AtomDelta.Add("
        "id=1, attributes=AtomForm.parse('O')))]))"
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
            ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]} '
            ':deltas [{:stereo-atom {:modify [0 "Th0"]}}]}',
            id="stereo-modify",
        ),
        pytest.param(
            '{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}',
            id="molecule-constraint",
        ),
    ],
)
def test_reaction_parse(text):
    first = Reaction.parse(text)

    canonical = str(first)
    second = Reaction.parse(canonical)

    assert second == first
    assert str(second) == canonical
    assert second.lhs is not first.lhs
    assert second.deltas is not first.deltas


def test_reaction_parse_defaults():
    reaction = Reaction.parse(
        '{:lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"]} '
        ':deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}]}',
        defaults=ReactionDefaults.ground(),
    )

    assert reaction == Reaction.parse(
        '{:lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]} '
        ":deltas [{:atom {:add "
        '"O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"}}]}'
    )


def test_reaction_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        Reaction.parse("not edn")


def test_reaction_parse_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Reaction.parse\\(\\) takes 1 positional arguments but 2 were given$"
        ),
    ):
        Reaction.parse(
            '{:lhs {:atoms ["C"]} :deltas []}',
            ReactionDefaults.ground(),
        )


@pytest.mark.parametrize(
    ("source", "lhs_entity", "delta_entity"),
    [
        pytest.param(
            '{:lhs {:atoms [[:lhs "C"]]} '
            ':deltas [{:atom {:add [:delta "O"]}}]}',
            Entity.Atom(0),
            Entity.Atom(1),
            id="atom",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "C"] '
            ':bonds [{:id :lhs :atoms [0 1] :attrs "1"}]} '
            ':deltas [{:bond {:add '
            '{:id :delta :atoms [0 1] :attrs "2"}}}]}',
            Entity.Bond(0),
            Entity.Bond(1),
            id="bond",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "N"] '
            ':dative-bonds [{:id :lhs :donors [0] '
            ':acceptor 1 :attrs "1#R"}]} '
            ':deltas [{:dative-bond {:add '
            '{:id :delta :donors [0] :acceptor 1 :attrs "1#R"}}}]}',
            Entity.DativeBond(0),
            Entity.DativeBond(1),
            id="dative-bond",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "C"] '
            ':aromatic-systems [{:id :lhs :atoms [0 1] :attrs "*#e2"}]} '
            ':deltas [{:aromatic-system {:add '
            '{:id :delta :atoms [0 1] :attrs "*#e2"}}}]}',
            Entity.AromaticSystem(0),
            Entity.AromaticSystem(1),
            id="aromatic-system",
        ),
        pytest.param(
            '{:lhs {:atoms ["B" "H" "B"] '
            ':multicenter-bonds [{:id :lhs :atoms [0 1 2] '
            ':attrs "[1,0,1]#e2"}]} '
            ':deltas [{:multicenter-bond {:add '
            '{:id :delta :atoms [0 1 2] :attrs "[1,0,1]#e2"}}}]}',
            Entity.MulticenterBond(0),
            Entity.MulticenterBond(1),
            id="multicenter-bond",
        ),
        pytest.param(
            '{:lhs {:atoms ["N" "H"] '
            ':noncovalent-bonds [{:id :lhs :atoms [0 1] :attrs "Hbd"}]} '
            ':deltas [{:noncovalent-bond {:add '
            '{:id :delta :atoms [0 1] :attrs "Hbd"}}}]}',
            Entity.NoncovalentBond(0),
            Entity.NoncovalentBond(1),
            id="noncovalent-bond",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
            ':stereo-atoms [{:id :lhs :site 0 :ligands [1 2 3 4] '
            ':attrs "Th1"}]} '
            ':deltas [{:stereo-atom {:add '
            '{:id :delta :site 0 :ligands [1 2 3 4] :attrs "Th2"}}}]}',
            Entity.StereoAtom(0),
            Entity.StereoAtom(1),
            id="stereo-atom",
        ),
        pytest.param(
            '{:lhs {:atoms ["C" "C" "C" "C"] '
            ':bonds [{:id :first :atoms [0 1] :attrs "2"} '
            '{:id :second :atoms [2 3] :attrs "2"}] '
            ':stereo-bonds [{:id :lhs :site :first '
            ':ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct1"}]} '
            ':deltas [{:stereo-bond {:add '
            '{:id :delta :site :second :ligands [0 [:h 2] 1 [:h 3]] :attrs "Ct2"}}}]}',
            Entity.StereoBond(0),
            Entity.StereoBond(1),
            id="stereo-bond",
        ),
    ],
)
def test_reaction_parse_with_metadata(source, lhs_entity, delta_entity):
    reaction, metadata = Reaction.parse_with_metadata(source)

    rendered = reaction.render_with_metadata(metadata)
    reparsed, reparsed_metadata = Reaction.parse_with_metadata(rendered)

    assert metadata.lhs.keyword(lhs_entity) == "lhs"
    assert metadata.lhs.entity("lhs") == lhs_entity
    assert metadata.delta_keyword(delta_entity) == "delta"
    assert metadata.delta_entity("delta") == delta_entity
    assert reparsed == reaction
    assert reparsed_metadata == metadata
    assert reparsed.render_with_metadata(reparsed_metadata) == rendered


def test_reaction_parse_with_metadata_aliases():
    source = (
        '{:lhs {:atoms [:lhs-c] :atom-aliases [:lhs-c "C"]} '
        ':atom-aliases [:delta-o "O"] '
        ':deltas [{:atom {:add [:added :delta-o]}}]}'
    )

    reaction, metadata = Reaction.parse_with_metadata(source)
    rendered = reaction.render_with_metadata(metadata)
    reparsed, reparsed_metadata = Reaction.parse_with_metadata(rendered)

    assert repr(metadata.lhs) == (
        "MoleculeMetadata(keywords=[], atom_alias_count=1)"
    )
    assert repr(metadata) == (
        "ReactionMetadata(lhs=MoleculeMetadata(keywords=[], "
        "atom_alias_count=1), "
        'delta_keywords=[(Entity.Atom(1), "added")], '
        "reaction_atom_alias_count=1)"
    )
    assert reparsed == reaction
    assert reparsed_metadata == metadata
    assert reparsed.render_with_metadata(reparsed_metadata) == rendered


def test_reaction_parse_with_metadata_defaults():
    reaction, metadata = Reaction.parse_with_metadata(
        '{:lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"]} '
        ':deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}]}',
        defaults=ReactionDefaults.ground(),
    )

    assert reaction == Reaction(
        Molecule.from_entries(
            [
                AtomForm.parse(
                    "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"
                )
            ]
        ),
        Deltas(
            [
                Delta.Atom(
                    AtomDelta.Add(
                        id=1,
                        attributes=AtomForm.parse(
                            "O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"
                        ),
                    )
                )
            ]
        ),
    )
    assert metadata == ReactionMetadata()


def test_reaction_parse_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Reaction.parse_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        Reaction.parse_with_metadata(
            '{:lhs {:atoms ["C"]} :deltas []}',
            ReactionDefaults.ground(),
        )


@pytest.mark.parametrize(
    ("reaction", "defaults", "expected"),
    [
        pytest.param(
            Reaction(
                Molecule.from_entries([AtomForm(Element("C"))]),
                Deltas(
                    [
                        Delta.Atom(
                            AtomDelta.Add(
                                id=1,
                                attributes=AtomForm(Element("O")),
                            )
                        )
                    ]
                ),
            ),
            ReactionDefaults(),
            '{:deltas [{:atom {:add "O"}}] '
            ':lhs {:atoms ["C"] :bonds []}}',
            id="required",
        ),
        pytest.param(
            Reaction(
                Molecule.from_entries(
                    [
                        AtomForm.parse(
                            "C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"
                        )
                    ]
                ),
                Deltas(
                    [
                        Delta.Atom(
                            AtomDelta.Add(
                                id=1,
                                attributes=AtomForm.parse(
                                    "O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"
                                ),
                            )
                        )
                    ]
                ),
            ),
            ReactionDefaults.ground(),
            '{:deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}] '
            ':lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}}',
            id="ground",
        ),
    ],
)
def test_reaction_render(reaction, defaults, expected):
    assert reaction.render(defaults=defaults) == expected


def test_reaction_render_keyword_error():
    with pytest.raises(
        TypeError,
        match="^Reaction.render\\(\\) takes 0 positional arguments but 1 was given$",
    ):
        Reaction().render(ReactionDefaults())


def test_reaction_render_with_metadata():
    reaction, metadata = Reaction.parse_with_metadata(
        '{:lhs {:atoms [[:carbon "C"] [:oxygen "O"]]} '
        ':deltas [{:atom {:add [:nitrogen "N"]}}]}'
    )

    assert reaction.render_with_metadata(metadata) == (
        '{:deltas [{:atom {:add [:nitrogen "N"]}}] '
        ':lhs {:atoms [[:carbon "C"] [:oxygen "O"]] :bonds []}}'
    )
    assert reaction.render() == (
        '{:deltas [{:atom {:add "N"}}] '
        ':lhs {:atoms ["C" "O"] :bonds []}}'
    )


def test_reaction_render_with_metadata_error():
    metadata = ReactionMetadata()
    metadata.set_delta_keyword(Entity.Atom(1), "absent")

    with pytest.raises(
        MetadataError,
        match=(
            "^metadata entity is not introduced by an add delta: atom 1$"
        ),
    ):
        Reaction(
            Molecule.from_entries([AtomForm(Element("C"))]),
        ).render_with_metadata(metadata)


def test_reaction_render_with_metadata_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Reaction.render_with_metadata\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        Reaction().render_with_metadata(
            ReactionMetadata(),
            ReactionDefaults(),
        )


def test_reaction_from_sides():
    lhs = Molecule.from_entries([AtomForm(Element("C")), AtomForm(Element("O"))])
    rhs = Molecule.from_entries([AtomForm(Element("C")), AtomForm(Element("N"))])
    lhs_snapshot = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("O"))]
    )
    rhs_snapshot = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("N"))]
    )

    atom_correspondence = Correspondence([(0, 0)], 2, 2)
    reaction = Reaction.from_sides(lhs, rhs, atom_correspondence)

    assert reaction == Reaction(
        lhs_snapshot,
        Deltas(
            [
                Delta.Atom(AtomDelta.Remove(id=1, attributes=AtomForm(Element("O")))),
                Delta.Atom(AtomDelta.Add(id=2, attributes=AtomForm(Element("N")))),
            ]
        ),
    )
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot
    assert reaction.lhs is not lhs
    assert Reaction.from_sides(lhs, rhs, atom_correspondence) == reaction


def test_reaction_from_sides_snapshot():
    lhs = Molecule.from_entries([AtomForm(Element("C")), AtomForm(Element("O"))])
    rhs = Molecule.from_entries([AtomForm(Element("C")), AtomForm(Element("N"))])
    reaction = Reaction.from_sides(lhs, rhs, Correspondence([(0, 0)], 2, 2))
    expected = Reaction(reaction.lhs, reaction.deltas)

    lhs.atoms[0].charge = 1
    rhs.atoms[0].charge = -1

    assert reaction == expected
    assert reaction.lhs is not lhs

    reaction.lhs.atoms[0].charge = 2
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm(Element("F")))))

    assert reaction.lhs.atoms[0].charge == NumForm.Lit(2)
    assert reaction.deltas[-1] == Delta.Atom(
        AtomDelta.Add(id=3, attributes=AtomForm(Element("F")))
    )
    assert lhs.atoms[0].charge == NumForm.Lit(1)
    assert rhs.atoms[0].charge == NumForm.Lit(-1)


def test_reaction_from_sides_error():
    lhs = Molecule.from_entries([AtomForm(Element("C"))])
    rhs = Molecule.from_entries([AtomForm(Element("C"))])
    atom_correspondence = Correspondence([(0, 0)], 2, 1)

    with pytest.raises(
        ValueError,
        match="^atom correspondence is incompatible with the reaction sides$",
    ):
        Reaction.from_sides(lhs, rhs, atom_correspondence)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        pytest.param(
            "[CH4:1]>>[CH4:1]",
            '{:deltas [] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] '
            ":bonds []}}",
            id="mapped",
        ),
        pytest.param(
            "[CH4:1]>>[CH4:1].[OH2:2]",
            '{:deltas [{:atom {:add '
            '"O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"}}] '
            ':lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}',
            id="one-sided",
        ),
        pytest.param(
            "C>>O",
            '{:deltas [{:atom {:remove 0}} {:atom {:add '
            '"O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"}}] '
            ':lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}',
            id="unmapped",
        ),
    ],
)
def test_reaction_from_reaction_smiles(source, expected):
    assert Reaction.from_reaction_smiles(source) == Reaction.parse(expected)


def test_reaction_from_reaction_smiles_io_config():
    with pytest.raises(
        UnderdeterminedError,
        match="^reactants: resolution underdetermined$",
    ):
        Reaction.from_reaction_smiles(
            "C~C>>C.C",
            io_config=SmilesIoConfig.lenient(),
        )


def test_reaction_from_reaction_smiles_resolve_config():
    reaction = Reaction.from_reaction_smiles(
        "[cH+:1]1[cH:2][cH:3]1>>[cH+:1]1[cH:2][cH:3]1",
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                reset_aromatic_valence=False,
            ),
            stereo=StereoResolveConfig(),
        ),
    )

    assert reaction == Reaction.parse(
        '{:deltas [] :lhs {:aromatic-systems '
        '[{:atoms [0 1 2] :attrs "[0,1,1]#c0#u0#s"}] '
        ':atoms ["C#i=#c+#h#n0#u0#s#v2#d0#t0#a0#m!" '
        '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
        '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
        ':bonds [[0 2 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
        '[1 2 "1#c0#u0#s#a"]]}}'
    )


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        pytest.param(
            "[o:1]1[cH:2][cH:3][cH:4][cH:5]1>>"
            "[o:1]1[cH:2][cH:3][cH:4][cH:5]1",
            '{:deltas [] :lhs {:atoms '
            '["O#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}}',
            id="furan",
        ),
        pytest.param(
            "[s:1]1[cH:2][cH:3][cH:4][cH:5]1>>"
            "[s:1]1[cH:2][cH:3][cH:4][cH:5]1",
            '{:deltas [] :lhs {:atoms '
            '["S#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}}',
            id="thiophene",
        ),
        pytest.param(
            "[nH:1]1[cH:2][cH:3][cH:4][cH:5]1>>"
            "[nH:1]1[cH:2][cH:3][cH:4][cH:5]1",
            '{:deltas [] :lhs {:atoms '
            '["N#i=#c0#h#n0#u0#s#v2#d0#t0#a2#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" '
            '"C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] '
            ':bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] '
            '[1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"] '
            '[3 4 "1#c0#u0#s#a"]]}}',
            id="pyrrole",
        ),
    ],
)
def test_reaction_from_reaction_smiles_aromaticity_policy(source, expected):
    default = ChemistryModel.default()

    assert Reaction.from_reaction_smiles(
        source,
        chemistry_model=ChemistryModel(
            connectivity=ChemistryModel.default().connectivity,
            valence=default.valence,
            aromaticity=AromaticityModel.mdl(),
            stereo=default.stereo,
        ),
        resolve_config=ResolveConfig(
            aromaticity=AromaticityResolveConfig(
                aromatic_valence_failure=AromaticityFailurePolicy.Keep
            ),
            stereo=StereoResolveConfig(),
        ),
    ) == Reaction.parse(expected)


@pytest.mark.parametrize(
    ("source", "kwargs", "error_type", "message"),
    [
        pytest.param(" C>>C", {}, ParseError, "Leading whitespace", id="syntax"),
        pytest.param(
            "C[S@]C>>",
            {},
            ModelConversionError,
            "reactants: tetrahedral stereo at atom 1 with 2 ligands, "
            "expected 3 or 4 ligands",
            id="model-conversion",
        ),
        pytest.param(
            "[nH]1cccc1>>",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.Clar(
                        scope=ElementScope.Any(),
                        ring_limits=RingLimits(),
                    ),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "reactants: clar: non-benzenoid input: Clar model requires "
            "benzenoid input but non-carbon aromatic atoms are present",
            id="contradiction",
        ),
        pytest.param(
            "o1cccc1>>C",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "reactants: aromaticity inconsistency: aromatic valence at atom "
            "AtomId(0) cannot produce a valid aromatic system",
            id="mdl-furan",
        ),
        pytest.param(
            "s1cccc1>>C",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "reactants: aromaticity inconsistency: aromatic valence at atom "
            "AtomId(0) cannot produce a valid aromatic system",
            id="mdl-thiophene",
        ),
        pytest.param(
            "[nH]1cccc1>>C",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
                    valence=ChemistryModel.default().valence,
                    aromaticity=AromaticityModel.mdl(),
                    stereo=ChemistryModel.default().stereo,
                )
            },
            ContradictionError,
            "reactants: aromaticity inconsistency: aromatic valence at atom "
            "AtomId(0) cannot produce a valid aromatic system",
            id="mdl-pyrrole",
        ),
        pytest.param(
            "*>>",
            {},
            UnderdeterminedError,
            "reactants: resolution underdetermined",
            id="underdetermined",
        ),
        pytest.param(
            "c1ccccc1>>",
            {
                "chemistry_model": ChemistryModel(
                    connectivity=ChemistryModel.default().connectivity,
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
            "reactants: hmo: missing parameters: no Van-Catledge parameters "
            "for C with 0 pi-electrons",
            id="execution",
        ),
        pytest.param(
            "[C:1].[O:1]>>[C:1]",
            {},
            ModelConversionError,
            "atom-map class 1 cannot be projected into one correspondence "
            "(reactant atoms: 2, product atoms: 1)",
            id="ambiguous-map",
        ),
        pytest.param(
            "C>O>C",
            {},
            ModelConversionError,
            "reaction agents cannot be represented in Reaction",
            id="agents",
        ),
    ],
)
def test_reaction_from_reaction_smiles_error(
    source,
    kwargs,
    error_type,
    message,
):
    with pytest.raises(error_type, match=f"^{re.escape(message)}$"):
        Reaction.from_reaction_smiles(source, **kwargs)


def test_reaction_from_reaction_smiles_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^Reaction.from_reaction_smiles\\(\\) takes 1 positional "
            "arguments but 2 were given$"
        ),
    ):
        Reaction.from_reaction_smiles(
            "C>>C",
            SmilesIoConfig.opensmiles(),
        )


def test_reaction_from_reaction_smiles_ownership():
    source = "[CH4:1]>>[CH4:1]"
    io_config = SmilesIoConfig.opensmiles()
    chemistry_model = ChemistryModel.default()
    resolve_config = ResolveConfig(
        aromaticity=AromaticityResolveConfig(),
        stereo=StereoResolveConfig(),
    )
    reaction = Reaction.from_reaction_smiles(
        source,
        io_config=io_config,
        chemistry_model=chemistry_model,
        resolve_config=resolve_config,
    )

    del source, io_config, chemistry_model, resolve_config

    assert reaction == Reaction.parse(
        '{:deltas [] :lhs {:atoms '
        '["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}'
    )


def test_reaction_str_components():
    reaction = Reaction.parse('{:lhs {:atoms ["C"]} :deltas []}')

    reaction.lhs.atoms[0].charge = 1
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("O")))))

    assert str(reaction) == reaction.render()
    assert reaction.render() == (
        '{:deltas [{:atom {:add "O"}}] '
        ':lhs {:atoms ["C#c+"] :bonds []}}'
    )


def test_reaction_parse_repr():
    reaction = Reaction.parse('{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}')

    assert repr(reaction) == (
        "Reaction(lhs=Molecule(atoms=1, bonds=0), "
        "deltas=Deltas([Delta.Atom(AtomDelta.Add("
        "id=1, attributes=AtomForm.parse('O')))]))"
    )


def test_reaction_reverse():
    source = Reaction.parse(
        '{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}'
    )
    snapshot = Reaction(source.lhs, source.deltas)

    reversed_reaction = source.reverse()
    roundtrip = reversed_reaction.reverse()

    assert reversed_reaction.lhs == Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("N"))]
    )
    assert roundtrip.lhs == source.lhs
    assert roundtrip.deltas.normalize() == source.deltas.normalize()
    assert source == snapshot
    assert reversed_reaction.lhs is not source.lhs
    assert reversed_reaction.deltas is not source.deltas

    reversed_reaction.lhs.atoms[0].charge = 1
    reversed_reaction.deltas.append(
        Delta.Atom(AtomDelta.Add(id=2, attributes=AtomForm(Element("F"))))
    )
    assert reversed_reaction.lhs.atoms[0].charge == NumForm.Lit(1)
    assert len(reversed_reaction.deltas) == 3


@pytest.mark.parametrize(
    "config",
    [
        pytest.param(
            ReactionCompositionConfig(
                common_subgraph_enumeration_algorithm=(
                    CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking()
                )
            ),
            id="modular-product-backtracking",
        ),
        pytest.param(
            ReactionCompositionConfig(
                common_subgraph_enumeration_algorithm=(
                    CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
                )
            ),
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
def test_reaction_compose(config, first, second, expected):
    first = Reaction.parse(first)
    second = Reaction.parse(second)

    composites = first.compose(
        second,
        config=config,
    )

    expected = [Reaction.parse(reaction) for reaction in expected]
    assert len(composites) == len(expected)
    for reaction in expected:
        assert reaction in composites


def test_reaction_compose_default():
    first = Reaction.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second = Reaction.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    fused = Reaction.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    disjoint = Reaction.parse(
        '{:lhs {:atoms ["C#c0" "C#c+"]} :deltas '
        '[{:atom {:modify [0 "#c+"]}} '
        '{:atom {:modify [1 "#c+2"]}}]}'
    )

    omitted = first.compose(second)
    explicit = first.compose(
        second,
        config=ReactionCompositionConfig(
            common_subgraph_enumeration_algorithm=(
                CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
            )
        ),
    )

    expected = [disjoint, fused]
    assert omitted == explicit
    assert len(omitted) == len(expected)
    for reaction in expected:
        assert reaction in omitted


def test_reaction_compose_error():
    first = Reaction()
    second = Reaction()

    with pytest.raises(
        TypeError,
        match="^Reaction.compose\\(\\) got an unexpected keyword argument 'algorithm'$",
    ):
        first.compose(
            second,
            algorithm=CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
        )


def test_reaction_compose_snapshot():
    first = Reaction.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second = Reaction.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    first_snapshot = Reaction(first.lhs, first.deltas)
    second_snapshot = Reaction(second.lhs, second.deltas)

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
            Delta.Atom(AtomDelta.Add(id=8, attributes=AtomForm(Element("Cl"))))
        )

        assert composite.lhs.atoms[0].charge == NumForm.Lit(7)
        assert composite.deltas[-1] == Delta.Atom(
            AtomDelta.Add(id=8, attributes=AtomForm(Element("Cl")))
        )

    assert first == first_snapshot
    assert second == second_snapshot


def test_reaction_apply():
    reaction = Reaction.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    host = Reaction.parse('{:lhs {:atoms ["C#c0" "C#c0"]} :deltas []}').lhs
    reaction_snapshot = Reaction(reaction.lhs, reaction.deltas)
    host_snapshot = Reaction(host).lhs
    first_product = Reaction.parse('{:lhs {:atoms ["C#c+" "C#c0"]} :deltas []}').lhs
    second_product = Reaction.parse('{:lhs {:atoms ["C#c0" "C#c+"]} :deltas []}').lhs

    application = reaction.apply(host)

    assert iter(application) is application
    assert reaction == reaction_snapshot
    assert host == host_snapshot

    reaction.lhs.atoms[0].charge = 7
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm(Element("N")))))
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
    atom_correspondence = first.atom_correspondence

    assert first.comap == comap
    assert first.comap is not comap
    assert first.atom_correspondence == atom_correspondence
    assert first.atom_correspondence is not atom_correspondence
    assert atom_correspondence == comap.atoms
    assert (
        Reaction.from_sides(first.lhs, first.rhs, atom_correspondence)
        == first.to_reaction()
    )
    assert comap.atoms.matched_pairs == [(0, 0), (1, 1)]
    assert comap.atoms.left_count == 2
    assert comap.atoms.right_count == 2
    assert comap.atoms.left_unmatched == []
    assert comap.atoms.right_unmatched == []
    for entity_map in (
        comap.bonds,
        comap.dative_bonds,
        comap.aromatic_systems,
        comap.multicenter_bonds,
        comap.noncovalent_bonds,
        comap.stereo_atoms,
        comap.stereo_bonds,
    ):
        assert entity_map.matched_pairs == []
        assert entity_map.left_count == 0
        assert entity_map.right_count == 0
        assert entity_map.left_unmatched == []
        assert entity_map.right_unmatched == []

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

    second_step = Reaction.parse(
        '{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    config = ReactionApplicationConfig(
        match_algorithm=SubstructureMatchAlgorithm.Incidence(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Vf2(),
    )
    following = next(second_step.apply(first.rhs, config=config))
    chained = first.chain(following)
    chained_product = Reaction.parse(
        '{:lhs {:atoms ["C#c+2" "C#c0"]} :deltas []}'
    ).lhs

    assert chained.lhs == host_snapshot
    assert chained.rhs == chained_product
    assert first.rhs == first_product
    assert following.lhs == first_product

    expected_reaction = Reaction.parse(
        '{:lhs {:atoms ["C#c0" "C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}'
    )
    recovered = chained.to_reaction()
    independent = chained.to_reaction()

    assert recovered == expected_reaction
    assert independent == expected_reaction
    assert recovered.lhs is not independent.lhs
    assert recovered.deltas is not independent.deltas

    recovered.lhs.atoms[1].charge = 9
    recovered.deltas.append(Delta.Atom(AtomDelta.Add(id=2, attributes=AtomForm(Element("O")))))

    assert independent == expected_reaction
    assert chained.lhs == host_snapshot
    assert chained.rhs == chained_product

    zero = Reaction.parse('{:lhs {:atoms ["N"]} :deltas []}').apply(host_snapshot)

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
def test_reaction_apply_config(config):
    reaction = Reaction.parse(
        '{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    host = Reaction.parse('{:lhs {:atoms ["C#c0" "C#c0"]} :deltas []}').lhs
    expected = [
        Reaction.parse('{:lhs {:atoms ["C#c+" "C#c0"]} :deltas []}').lhs,
        Reaction.parse('{:lhs {:atoms ["C#c0" "C#c+"]} :deltas []}').lhs,
    ]

    application = (
        reaction.apply(host) if config is None else reaction.apply(host, config=config)
    )

    assert [derivation.rhs for derivation in application] == expected


def test_reaction_apply_config_error():
    with pytest.raises(TypeError):
        Reaction().apply(Molecule(), ReactionApplicationConfig())


def test_reaction_apply_rejection():
    reaction = Reaction.parse(
        '{:lhs {:atoms ["C"]} :deltas [{:atom {:remove 0}}]}'
    )
    host = Reaction.parse(
        '{:lhs {:atoms ["C" "C" "C" "O"] :bonds [[1 3 "1"]]} :deltas []}'
    ).lhs
    expected = [
        Reaction.parse(
            '{:lhs {:atoms ["C" "C" "O"] :bonds [[0 2 "1"]]} :deltas []}'
        ).lhs,
        Reaction.parse(
            '{:lhs {:atoms ["C" "C" "O"] :bonds [[1 2 "1"]]} :deltas []}'
        ).lhs,
    ]

    assert [derivation.rhs for derivation in reaction.apply(host)] == expected


def test_reaction_apply_iteration_error():
    reaction = Reaction.parse(
        "{:lhs {:atoms [\"C\"] "
        ":constraints [{:charge-sum {:atoms [0] :sum 0}}]} "
        ":deltas [{:constraint {:remove {:charge-sum {:atoms [0] :sum 0}}}}]}"
    )
    host = Reaction.parse('{:lhs {:atoms ["C"]} :deltas []}').lhs

    application = reaction.apply(host)

    with pytest.raises(
        TransactionError,
        match=r"^missing constraint entry on remove$",
    ):
        next(application)
    with pytest.raises(StopIteration):
        next(application)
    with pytest.raises(StopIteration):
        next(application)


def test_reaction_workflow():
    lhs = Reaction.parse(
        '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]} '
        ":deltas []}"
    ).lhs
    rhs = Reaction.parse(
        '{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] '
        ':stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}] '
        ":constraints [{:connected {}}]} :deltas []}"
    ).lhs
    lhs_snapshot = Reaction(lhs).lhs
    rhs_snapshot = Reaction(rhs).lhs
    expected_forward = Reaction.parse(
        "{:deltas ["
        "{:stereo-atom {:add {:ligands [1 2 3 4] :site 0 :attrs :ccw}}} "
        "{:constraint {:add {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]]}}"
    )

    reaction = Reaction.from_sides(
        lhs,
        rhs,
        Correspondence([(atom_id, atom_id) for atom_id in range(5)], 5, 5),
    )

    assert reaction == expected_forward
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot
    assert reaction.lhs is not lhs

    normalized = Reaction(reaction.lhs, reaction.deltas)

    assert normalized == expected_forward
    assert reaction == expected_forward
    assert normalized.lhs is not reaction.lhs
    assert normalized.deltas is not reaction.deltas

    reaction.lhs.atoms[0].charge = 1
    reaction.deltas.append(Delta.Atom(AtomDelta.Add(id=5, attributes=AtomForm(Element("Xe")))))

    assert normalized == expected_forward
    assert lhs == lhs_snapshot
    assert rhs == rhs_snapshot

    rendered = str(normalized)

    assert rendered == (
        "{:deltas [{:stereo-atom {:add {:attrs :ccw :ligands [1 2 3 4] "
        ":site 0}}} {:constraint {:add {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]]}}"
    )
    assert normalized == expected_forward

    parsed = Reaction.parse(rendered)

    assert parsed == expected_forward
    assert normalized == expected_forward
    assert parsed.lhs is not normalized.lhs
    assert parsed.deltas is not normalized.deltas

    normalized.lhs.atoms[0].charge = 2
    normalized.deltas.append(
        Delta.Atom(AtomDelta.Add(id=5, attributes=AtomForm(Element("Ne"))))
    )

    assert parsed == expected_forward

    expected_reverse = Reaction.parse(
        "{:deltas ["
        "{:stereo-atom {:remove 0}} "
        "{:constraint {:remove {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]] :constraints [{:connected {}}] "
        ":stereo-atoms [{:ligands [1 2 3 4] :site 0 :attrs :ccw}]}}"
    )
    reversed_reaction = parsed.reverse()

    assert reversed_reaction == expected_reverse
    assert parsed == expected_forward
    assert reversed_reaction.lhs is not parsed.lhs
    assert reversed_reaction.deltas is not parsed.deltas

    parsed.lhs.atoms[0].charge = 3
    parsed.deltas.append(Delta.Atom(AtomDelta.Add(id=5, attributes=AtomForm(Element("Ar")))))

    assert reversed_reaction == expected_reverse

    second = Reaction.parse(
        '{:lhs {:atoms ["Xe#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}'
    )
    second_snapshot = Reaction(second.lhs, second.deltas)
    expected_composite = Reaction.parse(
        "{:deltas ["
        '{:atom {:modify [5 "#c+"]}} '
        "{:stereo-atom {:remove 0}} "
        "{:constraint {:remove {:connected {}}}}] "
        ':lhs {:atoms ["C" "F" "Cl" "Br" "I" "Xe#c0"] '
        ":bonds [[0 1 :single] [0 2 :single] [0 3 :single] "
        "[0 4 :single]] :constraints [{:connected {}}] "
        ":stereo-atoms [{:ligands [1 2 3 4] :site 0 :attrs :ccw}]}}"
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
        Delta.Atom(AtomDelta.Add(id=5, attributes=AtomForm(Element("Kr"))))
    )

    assert composites == [expected_composite]
    assert second == second_snapshot

    composite_lhs = composites[0].lhs
    composite_deltas = composites[0].deltas
    composite_lhs.atoms[5].charge = 2
    composite_deltas.append(Delta.Atom(AtomDelta.Add(id=6, attributes=AtomForm(Element("O")))))

    assert composites[0].lhs is composite_lhs
    assert composites[0].deltas is composite_deltas
    assert composites[0].lhs.atoms[5].charge == NumForm.Lit(2)
    assert composites[0].deltas[-1] == Delta.Atom(
        AtomDelta.Add(id=6, attributes=AtomForm(Element("O")))
    )
    assert second == second_snapshot
