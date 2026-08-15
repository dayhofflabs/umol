import re

import pytest

from umol import (
    AromaticityConfig,
    AromaticityModel,
    AromaticityRule,
    AromaticityTieBreak,
    AtomForm,
    AtomTypeRegistry,
    ChemistryModel,
    ConnectedComponentsAlgorithm,
    ConnectivityModel,
    Element,
    ElementScope,
    MaximumIndependentSetAlgorithm,
    RingConfig,
    RingLimits,
    StereoKind,
    StereoKindModel,
    StereoModel,
    ValenceCandidateSource,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
    ValenceTieBreak,
)


def test_atom_type_registry_default():
    registry = AtomTypeRegistry.default()

    assert registry == AtomTypeRegistry.default()
    assert repr(registry) == "AtomTypeRegistry.default()"
    assert registry.content_hash_hex == f"{registry.content_hash:016x}"


def test_atom_type_registry_from_atoms():
    atom = AtomForm.parse("C#c0#v4")
    registry = AtomTypeRegistry.from_atoms(
        [atom, AtomForm.parse("C#c+#v3"), AtomForm.parse("O#c0#v2")]
    )

    assert registry == AtomTypeRegistry.from_atoms(
        [
            AtomForm.parse("C#c0#v4"),
            AtomForm.parse("C#c+#v3"),
            AtomForm.parse("O#c0#v2"),
        ]
    )
    assert registry != AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])

    atom.charge = 1
    assert registry.patterns_for_element_and_charge(Element("C"), 0) == [
        AtomForm.parse("C#c0#v4")
    ]


@pytest.mark.parametrize(
    ("atom", "message"),
    [
        (
            AtomForm.parse("*#c0"),
            "atom type registry entry 0 must have a literal element",
        ),
        (
            AtomForm.parse("C"),
            "atom type registry entry 0 must have a literal charge",
        ),
        (
            AtomForm.parse("C#c128"),
            "atom type registry entry 0 charge 128 is outside -128..=127",
        ),
    ],
)
def test_atom_type_registry_from_atoms_error(atom, message):
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        AtomTypeRegistry.from_atoms([atom])


def test_atom_type_registry_from_toml():
    registry = AtomTypeRegistry.from_toml(
        '[C]\n0 = ["C#c0#v4"]\n1 = ["C#c+#v3"]\n[O]\n0 = ["O#c0#v2"]'
    )

    assert registry.patterns_for_element(Element("C")) == [
        AtomForm.parse("C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!"),
        AtomForm.parse("C#i=#c+#h0#n0#u0#s#v3#d0#t0#a!#m!"),
    ]
    assert registry.patterns_for_element_and_charge(Element("O"), 0) == [
        AtomForm.parse("O#i=#c0#h0#n0#u0#s#v2#d0#t0#a!#m!")
    ]


def test_atom_type_registry_from_toml_error():
    with pytest.raises(
        ValueError,
        match="^invalid atom type registry: unknown element: X$",
    ):
        AtomTypeRegistry.from_toml('[X]\n0 = ["X#c0"]')


def test_atom_type_registry_patterns():
    registry = AtomTypeRegistry.from_atoms(
        [AtomForm.parse("C#c0#v4"), AtomForm.parse("C#c+#v3")]
    )
    patterns = registry.patterns_for_element(Element("C"))

    patterns[0].charge = -1
    patterns.pop()

    assert registry.patterns_for_element(Element("C")) == [
        AtomForm.parse("C#c0#v4"),
        AtomForm.parse("C#c+#v3"),
    ]


def test_atom_type_registry_repr():
    registry = AtomTypeRegistry.from_atoms(
        [AtomForm.parse("C#c0#v4"), AtomForm.parse("O#c0#v2")]
    )

    assert repr(registry) == (
        'AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4"), '
        'AtomForm.parse("O#c0#v2")])'
    )


def test_atom_type_registry_mutation():
    registry = AtomTypeRegistry.default()

    with pytest.raises(AttributeError):
        registry.content_hash = 0
    assert not hasattr(registry, "add")


@pytest.mark.parametrize(
    ("target_covalences", "aromatic_valences", "expected"),
    [
        ([], [], ValenceEntry()),
        (
            [6, 2, 4, 2],
            [4, 2],
            ValenceEntry(
                target_covalences=[2, 2, 4, 6], aromatic_valences=[2, 4]
            ),
        ),
    ],
)
def test_valence_entry_new(target_covalences, aromatic_valences, expected):
    assert (
        ValenceEntry(
            target_covalences=target_covalences,
            aromatic_valences=aromatic_valences,
        )
        == expected
    )


def test_valence_entry_properties():
    entry = ValenceEntry(target_covalences=[6, 2, 4], aromatic_valences=[4, 2])

    target_covalences = entry.target_covalences
    aromatic_valences = entry.aromatic_valences
    target_covalences.append(8)
    aromatic_valences.clear()

    assert entry.target_covalences == [2, 4, 6]
    assert entry.aromatic_valences == [2, 4]


@pytest.mark.parametrize(
    ("entry", "expected"),
    [
        (
            ValenceEntry(),
            "ValenceEntry(target_covalences=[], aromatic_valences=[])",
        ),
        (
            ValenceEntry(target_covalences=[6, 2, 4], aromatic_valences=[4, 2]),
            "ValenceEntry(target_covalences=[2, 4, 6], aromatic_valences=[2, 4])",
        ),
    ],
)
def test_valence_entry_repr(entry, expected):
    assert repr(entry) == expected


def test_valence_entry_mutation():
    entry = ValenceEntry(target_covalences=[2, 4], aromatic_valences=[2])

    with pytest.raises(AttributeError):
        entry.target_covalences = [6]


def test_valence_table_new():
    entries = {
        Element("C"): ValenceEntry(
            target_covalences=[6, 2, 4], aromatic_valences=[3, 2]
        ),
        Element("O"): ValenceEntry(target_covalences=[2]),
    }
    table = ValenceTable(entries=entries)

    entries.clear()

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4, 6], aromatic_valences=[2, 3]
    )
    assert table.entry(Element("O")) == ValenceEntry(target_covalences=[2])
    assert table.entry(Element("N")) is None


def test_valence_table_default():
    table = ValenceTable.default()

    assert table == ValenceTable.default()
    assert repr(table) == "ValenceTable.default()"
    assert table.content_hash_hex == f"{table.content_hash:016x}"


def test_valence_table_from_toml():
    table = ValenceTable.from_toml(
        "[C]\ntarget_covalences = [6, 2, 4]\naromatic_valences = [3, 2]"
    )

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4, 6], aromatic_valences=[2, 3]
    )
    assert table.entry(Element("O")) is None


def test_valence_table_from_toml_error():
    with pytest.raises(ValueError, match="^invalid valence table: unknown element: X$"):
        ValenceTable.from_toml("[X]\ntarget_covalences = [1]")


def test_valence_table_entry():
    table = ValenceTable(
        entries={
            Element("C"): ValenceEntry(
                target_covalences=[2, 4], aromatic_valences=[2]
            )
        }
    )
    entry = table.entry(Element("C"))
    target_covalences = entry.target_covalences
    aromatic_valences = entry.aromatic_valences

    target_covalences.append(6)
    aromatic_valences.clear()

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4], aromatic_valences=[2]
    )


@pytest.mark.parametrize(
    ("table", "expected"),
    [
        (ValenceTable(entries={}), "ValenceTable(entries={})"),
        (
            ValenceTable(
                entries={
                    Element("C"): ValenceEntry(target_covalences=[4, 2]),
                    Element("O"): ValenceEntry(target_covalences=[2]),
                }
            ),
            "ValenceTable(entries={Element('C'): "
            "ValenceEntry(target_covalences=[2, 4], aromatic_valences=[]), "
            "Element('O'): ValenceEntry(target_covalences=[2], aromatic_valences=[])})",
        ),
    ],
)
def test_valence_table_repr(table, expected):
    assert repr(table) == expected


def test_valence_table_mutation():
    table = ValenceTable(entries={})

    with pytest.raises(AttributeError):
        table.content_hash = 0
    assert not hasattr(table, "insert")


def test_valence_candidate_source_atom_typing():
    registry = AtomTypeRegistry.from_atoms(
        [AtomForm.parse("C#c0#v4"), AtomForm.parse("O#c0#v2")]
    )
    source = ValenceCandidateSource.AtomTyping(registry=registry)

    assert isinstance(source, ValenceCandidateSource)
    assert isinstance(source, ValenceCandidateSource.AtomTyping)
    assert source.registry == registry
    assert source == ValenceCandidateSource.AtomTyping(registry=registry)
    assert source != ValenceCandidateSource.Counts(table=ValenceTable(entries={}))


def test_valence_candidate_source_counts():
    table = ValenceTable(
        entries={Element("C"): ValenceEntry(target_covalences=[4, 2])}
    )
    source = ValenceCandidateSource.Counts(table=table)

    assert isinstance(source, ValenceCandidateSource)
    assert isinstance(source, ValenceCandidateSource.Counts)
    assert source.table == table
    assert source == ValenceCandidateSource.Counts(table=table)
    assert source != ValenceCandidateSource.AtomTyping(
        registry=AtomTypeRegistry.from_atoms([])
    )


@pytest.mark.parametrize(
    ("variant", "payload"),
    [
        (ValenceCandidateSource.AtomTyping, AtomTypeRegistry.from_atoms([])),
        (ValenceCandidateSource.Counts, ValenceTable(entries={})),
    ],
)
def test_valence_candidate_source_new_error(variant, payload):
    with pytest.raises(TypeError):
        variant(payload)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (
            ValenceCandidateSource.AtomTyping(
                registry=AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])
            ),
            "ValenceCandidateSource.AtomTyping("
            "registry=AtomTypeRegistry.from_atoms(["
            'AtomForm.parse("C#c0#v4")]))',
        ),
        (
            ValenceCandidateSource.Counts(
                table=ValenceTable(
                    entries={
                        Element("C"): ValenceEntry(target_covalences=[4, 2])
                    }
                )
            ),
            "ValenceCandidateSource.Counts(table=ValenceTable(entries="
            "{Element('C'): "
            "ValenceEntry(target_covalences=[2, 4], aromatic_valences=[])}))",
        ),
    ],
)
def test_valence_candidate_source_repr(source, expected):
    assert repr(source) == expected


@pytest.mark.parametrize(
    ("source", "field", "value"),
    [
        (
            ValenceCandidateSource.AtomTyping(
                registry=AtomTypeRegistry.from_atoms([])
            ),
            "registry",
            AtomTypeRegistry.default(),
        ),
        (
            ValenceCandidateSource.Counts(table=ValenceTable(entries={})),
            "table",
            ValenceTable.default(),
        ),
    ],
)
def test_valence_candidate_source_mutation(source, field, value):
    with pytest.raises(AttributeError):
        setattr(source, field, value)


def test_valence_tie_break_equality():
    assert ValenceTieBreak.Strict == ValenceTieBreak.Strict
    assert ValenceTieBreak.Strict != ValenceTieBreak.MostSaturated


@pytest.mark.parametrize(
    ("tie_break", "expected"),
    [
        (ValenceTieBreak.Strict, "ValenceTieBreak.Strict"),
        (ValenceTieBreak.MostSaturated, "ValenceTieBreak.MostSaturated"),
    ],
)
def test_valence_tie_break_repr(tie_break, expected):
    assert repr(tie_break) == expected


def test_valence_model_new():
    registry = AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])
    model = ValenceModel(
        candidates=ValenceCandidateSource.AtomTyping(registry=registry),
        tie_break=ValenceTieBreak.MostSaturated,
    )

    assert model.candidates == ValenceCandidateSource.AtomTyping(
        registry=registry
    )
    assert model.tie_break == ValenceTieBreak.MostSaturated
    assert model == ValenceModel(
        candidates=ValenceCandidateSource.AtomTyping(registry=registry),
        tie_break=ValenceTieBreak.MostSaturated,
    )
    assert model != ValenceModel(
        candidates=ValenceCandidateSource.AtomTyping(registry=registry)
    )


def test_valence_model_new_default_tie_break():
    model = ValenceModel(
        candidates=ValenceCandidateSource.Counts(table=ValenceTable(entries={}))
    )

    assert model.tie_break == ValenceTieBreak.Strict


def test_valence_model_new_error():
    with pytest.raises(TypeError):
        ValenceModel(ValenceCandidateSource.Counts(table=ValenceTable(entries={})))


def test_valence_model_atom_typing():
    registry = AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])

    assert ValenceModel.atom_typing(registry) == ValenceModel(
        candidates=ValenceCandidateSource.AtomTyping(registry=registry),
        tie_break=ValenceTieBreak.Strict,
    )


def test_valence_model_counts():
    table = ValenceTable(
        entries={Element("C"): ValenceEntry(target_covalences=[4, 2])}
    )

    assert ValenceModel.counts(table) == ValenceModel(
        candidates=ValenceCandidateSource.Counts(table=table),
        tie_break=ValenceTieBreak.Strict,
    )


def test_valence_model_smiles():
    model = ValenceModel.smiles()
    table = model.candidates.table

    assert isinstance(model.candidates, ValenceCandidateSource.Counts)
    assert model.tie_break == ValenceTieBreak.MostSaturated
    assert table != ValenceTable.default()
    assert table.entry(Element("Cl")) == ValenceEntry(
        target_covalences=[1], aromatic_valences=[]
    )
    assert table.entry(Element("N")) == ValenceEntry(
        target_covalences=[3, 5], aromatic_valences=[1, 2]
    )


def test_valence_model_mdl():
    model = ValenceModel.mdl()
    table = model.candidates.table

    assert isinstance(model.candidates, ValenceCandidateSource.Counts)
    assert model.tie_break == ValenceTieBreak.MostSaturated
    assert table != ValenceTable.default()
    assert table.entry(Element("N")) == ValenceEntry(
        target_covalences=[3], aromatic_valences=[1, 2]
    )
    assert table.entry(Element("Cl")) == ValenceEntry(
        target_covalences=[1, 3, 5, 7], aromatic_valences=[]
    )


def test_valence_model_repr():
    model = ValenceModel(
        candidates=ValenceCandidateSource.AtomTyping(
            registry=AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])
        ),
    )

    assert repr(model) == (
        "ValenceModel(candidates=ValenceCandidateSource.AtomTyping("
        'registry=AtomTypeRegistry.from_atoms([AtomForm.parse("C#c0#v4")])), '
        "tie_break=ValenceTieBreak.Strict)"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("candidates", ValenceCandidateSource.Counts(table=ValenceTable.default())),
        ("tie_break", ValenceTieBreak.MostSaturated),
    ],
)
def test_valence_model_mutation(field, value):
    model = ValenceModel.smiles()

    with pytest.raises(AttributeError):
        setattr(model, field, value)


def test_element_scope_any():
    scope = ElementScope.Any()

    assert isinstance(scope, ElementScope)
    assert isinstance(scope, ElementScope.Any)
    assert scope == ElementScope.Any()
    assert scope != ElementScope.AllowList([])
    assert repr(scope) == "ElementScope.Any()"


@pytest.mark.parametrize(
    ("elements", "expected"),
    [
        ([], []),
        (
            [Element("C"), Element("N"), Element("C")],
            [Element("C"), Element("N"), Element("C")],
        ),
    ],
)
def test_element_scope_allow_list(elements, expected):
    scope = ElementScope.AllowList(elements)

    elements.clear()

    assert isinstance(scope, ElementScope)
    assert isinstance(scope, ElementScope.AllowList)
    assert scope.elements == expected
    assert scope == ElementScope.AllowList(expected)


def test_element_scope_elements():
    scope = ElementScope.AllowList([Element("C"), Element("N")])
    elements = scope.elements

    elements.append(Element("O"))

    assert scope.elements == [Element("C"), Element("N")]


@pytest.mark.parametrize(
    ("scope", "expected"),
    [
        (ElementScope.AllowList([]), "ElementScope.AllowList([])"),
        (
            ElementScope.AllowList([Element("C"), Element("N")]),
            "ElementScope.AllowList([Element('C'), Element('N')])",
        ),
    ],
)
def test_element_scope_repr(scope, expected):
    assert repr(scope) == expected


def test_element_scope_mutation():
    scope = ElementScope.AllowList([Element("C")])

    with pytest.raises(AttributeError):
        scope.elements = [Element("N")]


def test_ring_limits_default():
    limits = RingLimits()

    assert limits.min_ring_size == 3
    assert limits.max_ring_size == 22
    assert limits.include_fused is True
    assert limits.max_fused_combination == 6
    assert limits.max_fused_search == 10_000


def test_ring_limits_new():
    limits = RingLimits(
        min_ring_size=5,
        max_ring_size=18,
        include_fused=False,
        max_fused_combination=4,
        max_fused_search=2_500,
    )

    assert limits.min_ring_size == 5
    assert limits.max_ring_size == 18
    assert limits.include_fused is False
    assert limits.max_fused_combination == 4
    assert limits.max_fused_search == 2_500
    assert limits == RingLimits(
        min_ring_size=5,
        max_ring_size=18,
        include_fused=False,
        max_fused_combination=4,
        max_fused_search=2_500,
    )
    assert limits != RingLimits()


@pytest.mark.parametrize(
    "field",
    [
        "min_ring_size",
        "max_ring_size",
        "max_fused_combination",
        "max_fused_search",
    ],
)
def test_ring_limits_zero(field):
    limits = RingLimits(**{field: 0})

    assert getattr(limits, field) == 0


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("min_ring_size", -1),
        ("min_ring_size", 1 << 64),
        ("max_ring_size", -1),
        ("max_ring_size", 1 << 64),
        ("max_fused_combination", -1),
        ("max_fused_combination", 1 << 64),
        ("max_fused_search", -1),
        ("max_fused_search", 1 << 64),
    ],
)
def test_ring_limits_new_integer_error(field, value):
    with pytest.raises(OverflowError):
        RingLimits(**{field: value})


def test_ring_limits_new_positional_error():
    with pytest.raises(TypeError):
        RingLimits(3)


@pytest.mark.parametrize(
    ("limits", "expected"),
    [
        (
            RingLimits(),
            "RingLimits(min_ring_size=3, max_ring_size=22, include_fused=True, "
            "max_fused_combination=6, max_fused_search=10000)",
        ),
        (
            RingLimits(
                min_ring_size=5,
                max_ring_size=18,
                include_fused=False,
                max_fused_combination=4,
                max_fused_search=2_500,
            ),
            "RingLimits(min_ring_size=5, max_ring_size=18, include_fused=False, "
            "max_fused_combination=4, max_fused_search=2500)",
        ),
    ],
)
def test_ring_limits_repr(limits, expected):
    assert repr(limits) == expected


def test_ring_limits_mutation():
    limits = RingLimits()

    with pytest.raises(AttributeError):
        limits.min_ring_size = 5


def test_aromaticity_config_default():
    config = AromaticityConfig()

    assert config.ring_config == RingConfig()
    assert config.ring_config is not config.ring_config
    assert (
        config.connected_components_algorithm
        == ConnectedComponentsAlgorithm.Bfs()
    )
    assert (
        config.maximum_independent_set_algorithm
        == MaximumIndependentSetAlgorithm.BranchAndBound()
    )
    assert config == AromaticityConfig()


def test_aromaticity_config_new():
    ring_config = RingConfig()
    config = AromaticityConfig(
        ring_config=ring_config,
        connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(),
        maximum_independent_set_algorithm=(
            MaximumIndependentSetAlgorithm.BranchAndBound()
        ),
    )

    assert config.ring_config == ring_config
    assert config.ring_config is not ring_config
    assert (
        config.connected_components_algorithm
        == ConnectedComponentsAlgorithm.Bfs()
    )
    assert (
        config.maximum_independent_set_algorithm
        == MaximumIndependentSetAlgorithm.BranchAndBound()
    )


def test_aromaticity_config_new_error():
    with pytest.raises(TypeError):
        AromaticityConfig(
            RingConfig(),
            ConnectedComponentsAlgorithm.Bfs(),
            MaximumIndependentSetAlgorithm.BranchAndBound(),
        )


def test_aromaticity_config_repr():
    assert repr(AromaticityConfig()) == (
        "AromaticityConfig("
        "ring_config=RingConfig("
        "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara()), "
        "connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), "
        "maximum_independent_set_algorithm="
        "MaximumIndependentSetAlgorithm.BranchAndBound())"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("ring_config", RingConfig()),
        (
            "connected_components_algorithm",
            ConnectedComponentsAlgorithm.Bfs(),
        ),
        (
            "maximum_independent_set_algorithm",
            MaximumIndependentSetAlgorithm.BranchAndBound(),
        ),
    ],
)
def test_aromaticity_config_mutation(field, value):
    config = AromaticityConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


@pytest.mark.parametrize(
    ("variant", "ring_limits"),
    [
        (
            AromaticityRule.Hueckel,
            RingLimits(
                min_ring_size=4,
                max_ring_size=18,
                include_fused=False,
                max_fused_combination=3,
                max_fused_search=2_000,
            ),
        ),
        (
            AromaticityRule.Clar,
            RingLimits(
                min_ring_size=6,
                max_ring_size=14,
                include_fused=True,
                max_fused_combination=4,
                max_fused_search=1_500,
            ),
        ),
    ],
)
def test_aromaticity_rule_ring_variant(variant, ring_limits):
    rule = variant(ring_limits=ring_limits)

    assert isinstance(rule, AromaticityRule)
    assert isinstance(rule, variant)
    assert rule.ring_limits == ring_limits
    assert rule.ring_limits is not ring_limits
    assert rule == variant(ring_limits=ring_limits)
    assert rule != variant(ring_limits=RingLimits())


def test_aromaticity_rule_hmo():
    rule = AromaticityRule.Hmo(stabilization_threshold=0.375)

    assert isinstance(rule, AromaticityRule)
    assert isinstance(rule, AromaticityRule.Hmo)
    assert rule.stabilization_threshold == 0.375
    assert rule == AromaticityRule.Hmo(stabilization_threshold=0.375)
    assert rule != AromaticityRule.Hmo(stabilization_threshold=0.5)


@pytest.mark.parametrize(
    ("variant", "args"),
    [
        (AromaticityRule.Hueckel, (RingLimits(),)),
        (AromaticityRule.Hmo, (0.375,)),
        (AromaticityRule.Clar, (RingLimits(),)),
    ],
)
def test_aromaticity_rule_new_error(variant, args):
    with pytest.raises(TypeError):
        variant(*args)


@pytest.mark.parametrize(
    ("rule", "expected"),
    [
        (
            AromaticityRule.Hueckel(
                ring_limits=RingLimits(
                    min_ring_size=4,
                    max_ring_size=18,
                    include_fused=False,
                    max_fused_combination=3,
                    max_fused_search=2_000,
                )
            ),
            "AromaticityRule.Hueckel(ring_limits=RingLimits(min_ring_size=4, "
            "max_ring_size=18, include_fused=False, max_fused_combination=3, "
            "max_fused_search=2000))",
        ),
        (
            AromaticityRule.Hmo(stabilization_threshold=0.375),
            "AromaticityRule.Hmo(stabilization_threshold=0.375)",
        ),
        (
            AromaticityRule.Clar(
                ring_limits=RingLimits(
                    min_ring_size=6,
                    max_ring_size=14,
                    max_fused_combination=4,
                    max_fused_search=1_500,
                )
            ),
            "AromaticityRule.Clar(ring_limits=RingLimits(min_ring_size=6, "
            "max_ring_size=14, include_fused=True, max_fused_combination=4, "
            "max_fused_search=1500))",
        ),
    ],
)
def test_aromaticity_rule_repr(rule, expected):
    assert repr(rule) == expected


@pytest.mark.parametrize(
    ("rule", "field", "value"),
    [
        (
            AromaticityRule.Hueckel(ring_limits=RingLimits()),
            "ring_limits",
            RingLimits(min_ring_size=4),
        ),
        (
            AromaticityRule.Hmo(stabilization_threshold=0.375),
            "stabilization_threshold",
            0.5,
        ),
        (
            AromaticityRule.Clar(ring_limits=RingLimits()),
            "ring_limits",
            RingLimits(min_ring_size=6),
        ),
    ],
)
def test_aromaticity_rule_mutation(rule, field, value):
    with pytest.raises(AttributeError):
        setattr(rule, field, value)


def test_aromaticity_tie_break_equality():
    assert AromaticityTieBreak.Strict == AromaticityTieBreak.Strict
    assert AromaticityTieBreak.Strict != AromaticityTieBreak.MinElectronCount


@pytest.mark.parametrize(
    ("tie_break", "expected"),
    [
        (AromaticityTieBreak.Strict, "AromaticityTieBreak.Strict"),
        (AromaticityTieBreak.MinElectronCount, "AromaticityTieBreak.MinElectronCount"),
    ],
)
def test_aromaticity_tie_break_repr(tie_break, expected):
    assert repr(tie_break) == expected


def test_aromaticity_model_new():
    scope = ElementScope.AllowList([Element("C"), Element("N")])
    rule = AromaticityRule.Hmo(stabilization_threshold=0.375)
    model = AromaticityModel(
        scope=scope, rule=rule, tie_break=AromaticityTieBreak.MinElectronCount
    )

    assert model.scope == scope
    assert model.scope is not scope
    assert model.rule == rule
    assert model.rule is not rule
    assert model.tie_break == AromaticityTieBreak.MinElectronCount
    assert model == AromaticityModel(
        scope=scope, rule=rule, tie_break=AromaticityTieBreak.MinElectronCount
    )
    assert model != AromaticityModel(scope=scope, rule=rule)
    assert model != AromaticityModel(
        scope=ElementScope.Any(),
        rule=rule,
        tie_break=AromaticityTieBreak.MinElectronCount,
    )
    assert model != AromaticityModel(
        scope=scope,
        rule=AromaticityRule.Hmo(stabilization_threshold=0.5),
        tie_break=AromaticityTieBreak.MinElectronCount,
    )
    assert model != AromaticityModel(
        scope=scope,
        rule=AromaticityRule.Hueckel(ring_limits=RingLimits()),
        tie_break=AromaticityTieBreak.MinElectronCount,
    )


def test_aromaticity_model_new_default_tie_break():
    model = AromaticityModel(
        scope=ElementScope.Any(),
        rule=AromaticityRule.Hueckel(ring_limits=RingLimits()),
    )

    assert model.tie_break == AromaticityTieBreak.Strict


def test_aromaticity_model_new_error():
    with pytest.raises(TypeError):
        AromaticityModel(
            ElementScope.Any(), AromaticityRule.Hueckel(ring_limits=RingLimits())
        )


@pytest.mark.parametrize(
    ("preset", "expected"),
    [
        (
            AromaticityModel.daylight,
            AromaticityModel(
                scope=ElementScope.AllowList(
                    [
                        Element("C"),
                        Element("N"),
                        Element("O"),
                        Element("S"),
                        Element("Se"),
                        Element("As"),
                    ]
                ),
                rule=AromaticityRule.Hueckel(ring_limits=RingLimits()),
                tie_break=AromaticityTieBreak.MinElectronCount,
            ),
        ),
        (
            AromaticityModel.mdl,
            AromaticityModel(
                scope=ElementScope.AllowList([Element("C"), Element("N")]),
                rule=AromaticityRule.Hueckel(
                    ring_limits=RingLimits(min_ring_size=6)
                ),
                tie_break=AromaticityTieBreak.MinElectronCount,
            ),
        ),
        (
            AromaticityModel.permissive,
            AromaticityModel(
                scope=ElementScope.Any(),
                rule=AromaticityRule.Hueckel(ring_limits=RingLimits()),
                tie_break=AromaticityTieBreak.MinElectronCount,
            ),
        ),
    ],
)
def test_aromaticity_model_preset(preset, expected):
    assert preset() == expected


def test_aromaticity_model_repr():
    model = AromaticityModel(
        scope=ElementScope.AllowList([Element("C")]),
        rule=AromaticityRule.Hmo(stabilization_threshold=0.375),
    )

    assert repr(model) == (
        "AromaticityModel(scope=ElementScope.AllowList([Element('C')]), "
        "rule=AromaticityRule.Hmo(stabilization_threshold=0.375), "
        "tie_break=AromaticityTieBreak.Strict)"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("scope", ElementScope.AllowList([Element("C")])),
        ("rule", AromaticityRule.Hmo(stabilization_threshold=0.5)),
        ("tie_break", AromaticityTieBreak.MinElectronCount),
    ],
)
def test_aromaticity_model_mutation(field, value):
    model = AromaticityModel.permissive()

    with pytest.raises(AttributeError):
        setattr(model, field, value)


@pytest.mark.parametrize(
    ("scope", "fluxionality"),
    [
        (ElementScope.Any(), False),
        (ElementScope.AllowList([Element("C"), Element("N")]), True),
    ],
)
def test_stereo_kind_model_new(scope, fluxionality):
    model = StereoKindModel(scope=scope, fluxionality=fluxionality)

    assert model.scope == scope
    assert model.scope is not scope
    assert model.fluxionality is fluxionality
    assert model == StereoKindModel(scope=scope, fluxionality=fluxionality)


def test_stereo_kind_model_new_error():
    with pytest.raises(TypeError):
        StereoKindModel(ElementScope.Any(), False)


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        (
            StereoKindModel(scope=ElementScope.Any(), fluxionality=False),
            "StereoKindModel(scope=ElementScope.Any(), fluxionality=False)",
        ),
        (
            StereoKindModel(
                scope=ElementScope.AllowList([Element("C"), Element("N")]),
                fluxionality=True,
            ),
            "StereoKindModel(scope=ElementScope.AllowList([Element('C'), "
            "Element('N')]), fluxionality=True)",
        ),
    ],
)
def test_stereo_kind_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("scope", ElementScope.AllowList([Element("C")])),
        ("fluxionality", True),
    ],
)
def test_stereo_kind_model_mutation(field, value):
    model = StereoKindModel(scope=ElementScope.Any(), fluxionality=False)

    with pytest.raises(AttributeError):
        setattr(model, field, value)


def test_stereo_model_default():
    model = StereoModel.default()

    assert model.kind_models == {
        StereoKind.Tetrahedral: StereoKindModel(
            scope=ElementScope.Any(), fluxionality=False
        ),
        StereoKind.CisTrans: StereoKindModel(
            scope=ElementScope.Any(), fluxionality=False
        ),
    }
    assert model.para_stereo is False
    assert model == StereoModel.default()


def test_stereo_model_new():
    kind_models = {
        StereoKind.Tetrahedral: StereoKindModel(
            scope=ElementScope.Any(), fluxionality=False
        ),
        StereoKind.CisTrans: StereoKindModel(
            scope=ElementScope.AllowList([Element("C")]), fluxionality=True
        ),
        StereoKind.Axial: StereoKindModel(
            scope=ElementScope.AllowList([Element("N")]), fluxionality=False
        ),
        StereoKind.SquarePlanar: StereoKindModel(
            scope=ElementScope.AllowList([Element("O")]), fluxionality=True
        ),
        StereoKind.TrigonalBipyramidal: StereoKindModel(
            scope=ElementScope.AllowList([Element("S")]), fluxionality=False
        ),
        StereoKind.Octahedral: StereoKindModel(
            scope=ElementScope.AllowList([Element("Fe")]), fluxionality=True
        ),
    }
    expected = kind_models.copy()
    model = StereoModel(
        kind_models=kind_models,
        para_stereo=True,
    )
    kind_models.clear()

    assert model.kind_models == expected
    assert model.para_stereo is True
    for kind, kind_model in expected.items():
        assert model.kind_models[kind] is not kind_model


def test_stereo_model_kind_models():
    model = StereoModel.default()
    kind_models = model.kind_models
    kind_models.clear()

    assert model.kind_models == {
        StereoKind.Tetrahedral: StereoKindModel(
            scope=ElementScope.Any(), fluxionality=False
        ),
        StereoKind.CisTrans: StereoKindModel(
            scope=ElementScope.Any(), fluxionality=False
        ),
    }


def test_stereo_model_new_error():
    with pytest.raises(TypeError):
        StereoModel({}, False)


@pytest.mark.parametrize(
    ("field", "value"),
    [("max_iterations", 8), ("inconsistency", "strip")],
)
def test_stereo_model_new_field_error(field, value):
    with pytest.raises(TypeError):
        StereoModel(kind_models={}, para_stereo=False, **{field: value})


@pytest.mark.parametrize(
    "other",
    [
        StereoModel(
            kind_models={},
            para_stereo=False,
        ),
        StereoModel(
            kind_models={
                StereoKind.Tetrahedral: StereoKindModel(
                    scope=ElementScope.Any(), fluxionality=False
                ),
                StereoKind.CisTrans: StereoKindModel(
                    scope=ElementScope.Any(), fluxionality=False
                ),
            },
            para_stereo=True,
        ),
    ],
)
def test_stereo_model_equality(other):
    assert StereoModel.default() != other


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        (StereoModel.default(), "StereoModel.default()"),
        (
            StereoModel(
                kind_models={
                    StereoKind.Octahedral: StereoKindModel(
                        scope=ElementScope.AllowList([Element("Fe")]),
                        fluxionality=True,
                    ),
                    StereoKind.Tetrahedral: StereoKindModel(
                        scope=ElementScope.Any(), fluxionality=False
                    ),
                },
                para_stereo=True,
            ),
            "StereoModel(kind_models={StereoKind.Tetrahedral: "
            "StereoKindModel(scope=ElementScope.Any(), fluxionality=False), "
            "StereoKind.Octahedral: StereoKindModel(scope="
            "ElementScope.AllowList([Element('Fe')]), fluxionality=True)}, "
            "para_stereo=True)",
        ),
    ],
)
def test_stereo_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("kind_models", {}),
        ("para_stereo", True),
    ],
)
def test_stereo_model_mutation(field, value):
    model = StereoModel.default()

    with pytest.raises(AttributeError):
        setattr(model, field, value)


def test_connectivity_model_default():
    model = ConnectivityModel.default()

    assert model.allow_disconnected is True
    assert model.allow_disconnected_dative is True
    assert model.allow_disconnected_aromatic is False
    assert model.allow_disconnected_multicenter is True
    assert model.allow_disconnected_noncovalent is True
    assert model.allow_disconnected_stereo_atom is False
    assert model.allow_disconnected_stereo_bond is False
    assert model.allow_disconnected_constraints is True
    assert repr(model) == "ConnectivityModel.default()"


def test_connectivity_model_new():
    model = ConnectivityModel(
        allow_disconnected=False,
        allow_disconnected_dative=False,
        allow_disconnected_aromatic=False,
        allow_disconnected_multicenter=False,
        allow_disconnected_noncovalent=False,
        allow_disconnected_stereo_atom=False,
        allow_disconnected_stereo_bond=False,
        allow_disconnected_constraints=False,
    )

    assert model.allow_disconnected is False
    assert model != ConnectivityModel.default()


def test_chemistry_model_default():
    model = ChemistryModel.default()

    assert model.connectivity == ConnectivityModel.default()
    assert model.valence == ValenceModel.atom_typing(AtomTypeRegistry.default())
    assert model.aromaticity == AromaticityModel.daylight()
    assert model.stereo == StereoModel.default()
    assert model == ChemistryModel.default()


def test_chemistry_model_new():
    valence = ValenceModel.counts(
        ValenceTable(
            entries={
                Element("C"): ValenceEntry(target_covalences=[4]),
                Element("O"): ValenceEntry(target_covalences=[2]),
            }
        )
    )
    aromaticity = AromaticityModel(
        scope=ElementScope.AllowList([Element("C")]),
        rule=AromaticityRule.Clar(ring_limits=RingLimits(min_ring_size=6)),
    )
    stereo = StereoModel(
        kind_models=StereoModel.default().kind_models,
        para_stereo=True,
    )
    connectivity = ConnectivityModel.default()
    model = ChemistryModel(
        connectivity=connectivity,
        valence=valence,
        aromaticity=aromaticity,
        stereo=stereo,
    )

    assert model.connectivity == connectivity
    assert model.connectivity is not connectivity
    assert model.valence == valence
    assert model.valence is not valence
    assert model.aromaticity == aromaticity
    assert model.aromaticity is not aromaticity
    assert model.stereo == stereo
    assert model.stereo is not stereo


def test_chemistry_model_new_error():
    with pytest.raises(TypeError):
        ChemistryModel(
            ValenceModel.atom_typing(AtomTypeRegistry.default()),
            AromaticityModel.daylight(),
            StereoModel.default(),
        )


@pytest.mark.parametrize(
    "other",
    [
        ChemistryModel(
            connectivity=ConnectivityModel(
                allow_disconnected=False,
                allow_disconnected_dative=True,
                allow_disconnected_aromatic=False,
                allow_disconnected_multicenter=True,
                allow_disconnected_noncovalent=True,
                allow_disconnected_stereo_atom=False,
                allow_disconnected_stereo_bond=False,
                allow_disconnected_constraints=True,
            ),
            valence=ChemistryModel.default().valence,
            aromaticity=AromaticityModel.daylight(),
            stereo=StereoModel.default(),
        ),
        ChemistryModel(
            connectivity=ChemistryModel.default().connectivity,
            valence=ValenceModel.counts(ValenceTable(entries={})),
            aromaticity=AromaticityModel.daylight(),
            stereo=StereoModel.default(),
        ),
        ChemistryModel(
            connectivity=ChemistryModel.default().connectivity,
            valence=ValenceModel.atom_typing(AtomTypeRegistry.default()),
            aromaticity=AromaticityModel.permissive(),
            stereo=StereoModel.default(),
        ),
        ChemistryModel(
            connectivity=ChemistryModel.default().connectivity,
            valence=ValenceModel.atom_typing(AtomTypeRegistry.default()),
            aromaticity=AromaticityModel.daylight(),
            stereo=StereoModel(
                kind_models={},
                para_stereo=False,
            ),
        ),
    ],
)
def test_chemistry_model_equality(other):
    assert ChemistryModel.default() != other


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        (ChemistryModel.default(), "ChemistryModel.default()"),
        (
            ChemistryModel(
                connectivity=ChemistryModel.default().connectivity,
                valence=ValenceModel.counts(
                    ValenceTable(
                        entries={
                            Element("C"): ValenceEntry(target_covalences=[4])
                        }
                    )
                ),
                aromaticity=AromaticityModel(
                    scope=ElementScope.Any(),
                    rule=AromaticityRule.Hmo(stabilization_threshold=0.375),
                ),
                stereo=StereoModel(
                    kind_models=StereoModel.default().kind_models,
                    para_stereo=True,
                ),
            ),
            "ChemistryModel(connectivity=ConnectivityModel.default(), "
            "valence=ValenceModel(candidates=ValenceCandidateSource.Counts("
            "table=ValenceTable(entries={Element('C'): ValenceEntry("
            "target_covalences=[4], aromatic_valences=[])})), "
            "tie_break=ValenceTieBreak.Strict), aromaticity="
            "AromaticityModel(scope=ElementScope.Any(), "
            "rule=AromaticityRule.Hmo(stabilization_threshold=0.375), "
            "tie_break=AromaticityTieBreak.Strict), "
            "stereo=StereoModel(kind_models={"
            "StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any(), "
            "fluxionality=False), StereoKind.CisTrans: StereoKindModel(scope="
            "ElementScope.Any(), fluxionality=False)}, para_stereo=True))",
        ),
    ],
)
def test_chemistry_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        (
            "connectivity",
            ConnectivityModel(
                allow_disconnected=False,
                allow_disconnected_dative=True,
                allow_disconnected_aromatic=False,
                allow_disconnected_multicenter=True,
                allow_disconnected_noncovalent=True,
                allow_disconnected_stereo_atom=False,
                allow_disconnected_stereo_bond=False,
                allow_disconnected_constraints=True,
            ),
        ),
        (
            "valence",
            ValenceModel.counts(ValenceTable(entries={})),
        ),
        ("aromaticity", AromaticityModel.permissive()),
        (
            "stereo",
            StereoModel(
                kind_models={},
                para_stereo=False,
            ),
        ),
    ],
)
def test_chemistry_model_mutation(field, value):
    model = ChemistryModel.default()

    with pytest.raises(AttributeError):
        setattr(model, field, value)
