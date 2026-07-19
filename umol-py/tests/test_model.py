import re

import pytest

from umol import (
    AromaticityModel,
    AtomAst,
    AtomTypeRegistry,
    ChemistryModel,
    Element,
    ElementScope,
    InconsistencyPolicy,
    RingLimits,
    StereoKind,
    StereoKindModel,
    StereoModel,
    ValenceEntry,
    ValenceModel,
    ValenceTable,
)


def test_atom_type_registry_default():
    registry = AtomTypeRegistry.default()

    assert registry == AtomTypeRegistry.default()
    assert repr(registry) == "AtomTypeRegistry.default()"
    assert registry.content_hash_hex == f"{registry.content_hash:016x}"


def test_atom_type_registry_from_atoms():
    atom = AtomAst.parse("C#c0#v4")
    registry = AtomTypeRegistry.from_atoms(
        [atom, AtomAst.parse("C#c+#v3"), AtomAst.parse("O#c0#v2")]
    )

    assert registry == AtomTypeRegistry.from_atoms(
        [
            AtomAst.parse("C#c0#v4"),
            AtomAst.parse("C#c+#v3"),
            AtomAst.parse("O#c0#v2"),
        ]
    )
    assert registry != AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4")])

    atom.charge = 1
    assert registry.patterns_for_element_and_charge(Element("C"), 0) == [
        AtomAst.parse("C#c0#v4")
    ]


@pytest.mark.parametrize(
    ("atom", "message"),
    [
        (
            AtomAst.parse("*#c0"),
            "atom type registry entry 0 must have a literal element",
        ),
        (
            AtomAst.parse("C"),
            "atom type registry entry 0 must have a literal charge",
        ),
        (
            AtomAst.parse("C#c128"),
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
        AtomAst.parse("C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!"),
        AtomAst.parse("C#i=#c+#h0#n0#u0#s#v3#d0#t0#a!#m!"),
    ]
    assert registry.patterns_for_element_and_charge(Element("O"), 0) == [
        AtomAst.parse("O#i=#c0#h0#n0#u0#s#v2#d0#t0#a!#m!")
    ]


def test_atom_type_registry_from_toml_error():
    with pytest.raises(
        ValueError,
        match="^invalid atom type registry: unknown element: X$",
    ):
        AtomTypeRegistry.from_toml('[X]\n0 = ["X#c0"]')


def test_atom_type_registry_patterns():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("C#c+#v3")]
    )
    patterns = registry.patterns_for_element(Element("C"))

    patterns[0].charge = -1
    patterns.pop()

    assert registry.patterns_for_element(Element("C")) == [
        AtomAst.parse("C#c0#v4"),
        AtomAst.parse("C#c+#v3"),
    ]


def test_atom_type_registry_repr():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("O#c0#v2")]
    )

    assert repr(registry) == (
        'AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4"), '
        'AtomAst.parse("O#c0#v2")])'
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


def test_valence_model_atom_typing():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("O#c0#v2")]
    )
    model = ValenceModel.AtomTyping(registry=registry)

    assert isinstance(model, ValenceModel)
    assert isinstance(model, ValenceModel.AtomTyping)
    assert model.registry == registry
    assert model == ValenceModel.AtomTyping(registry=registry)
    assert model != ValenceModel.Counts(table=ValenceTable(entries={}))


def test_valence_model_counts():
    table = ValenceTable(
        entries={Element("C"): ValenceEntry(target_covalences=[4, 2])}
    )
    model = ValenceModel.Counts(table=table)

    assert isinstance(model, ValenceModel)
    assert isinstance(model, ValenceModel.Counts)
    assert model.table == table
    assert model == ValenceModel.Counts(table=table)
    assert model != ValenceModel.AtomTyping(
        registry=AtomTypeRegistry.from_atoms([])
    )


@pytest.mark.parametrize(
    ("variant", "payload"),
    [
        (ValenceModel.AtomTyping, AtomTypeRegistry.from_atoms([])),
        (ValenceModel.Counts, ValenceTable(entries={})),
    ],
)
def test_valence_model_new_error(variant, payload):
    with pytest.raises(TypeError):
        variant(payload)


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        (
            ValenceModel.AtomTyping(
                registry=AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4")])
            ),
            'ValenceModel.AtomTyping(registry=AtomTypeRegistry.from_atoms(['
            'AtomAst.parse("C#c0#v4")]))',
        ),
        (
            ValenceModel.Counts(
                table=ValenceTable(
                    entries={
                        Element("C"): ValenceEntry(target_covalences=[4, 2])
                    }
                )
            ),
            "ValenceModel.Counts(table=ValenceTable(entries={Element('C'): "
            "ValenceEntry(target_covalences=[2, 4], aromatic_valences=[])}))",
        ),
    ],
)
def test_valence_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("model", "field", "value"),
    [
        (
            ValenceModel.AtomTyping(registry=AtomTypeRegistry.from_atoms([])),
            "registry",
            AtomTypeRegistry.default(),
        ),
        (
            ValenceModel.Counts(table=ValenceTable(entries={})),
            "table",
            ValenceTable.default(),
        ),
    ],
)
def test_valence_model_mutation(model, field, value):
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


@pytest.mark.parametrize(
    ("variant", "scope", "ring_limits"),
    [
        (
            AromaticityModel.HueckelRule,
            ElementScope.Any(),
            RingLimits(
                min_ring_size=4,
                max_ring_size=18,
                include_fused=False,
                max_fused_combination=3,
                max_fused_search=2_000,
            ),
        ),
        (
            AromaticityModel.Clar,
            ElementScope.AllowList([Element("C")]),
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
def test_aromaticity_model_ring_variant(variant, scope, ring_limits):
    model = variant(scope=scope, ring_limits=ring_limits)

    assert isinstance(model, AromaticityModel)
    assert isinstance(model, variant)
    assert model.scope == scope
    assert model.scope is not scope
    assert model.ring_limits == ring_limits
    assert model.ring_limits is not ring_limits
    assert model == variant(scope=scope, ring_limits=ring_limits)


def test_aromaticity_model_hmo():
    scope = ElementScope.AllowList([Element("C"), Element("N")])
    model = AromaticityModel.Hmo(
        scope=scope,
        stabilization_threshold=0.375,
    )

    assert isinstance(model, AromaticityModel)
    assert isinstance(model, AromaticityModel.Hmo)
    assert model.scope == scope
    assert model.scope is not scope
    assert model.stabilization_threshold == 0.375
    assert model == AromaticityModel.Hmo(
        scope=scope,
        stabilization_threshold=0.375,
    )


@pytest.mark.parametrize(
    ("variant", "args"),
    [
        (AromaticityModel.HueckelRule, (ElementScope.Any(), RingLimits())),
        (AromaticityModel.Hmo, (ElementScope.Any(), 0.375)),
        (AromaticityModel.Clar, (ElementScope.Any(), RingLimits())),
    ],
)
def test_aromaticity_model_new_error(variant, args):
    with pytest.raises(TypeError):
        variant(*args)


@pytest.mark.parametrize(
    ("preset", "expected"),
    [
        (
            AromaticityModel.daylight,
            AromaticityModel.HueckelRule(
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
                ring_limits=RingLimits(),
            ),
        ),
        (
            AromaticityModel.mdl,
            AromaticityModel.HueckelRule(
                scope=ElementScope.AllowList([Element("C"), Element("N")]),
                ring_limits=RingLimits(min_ring_size=6),
            ),
        ),
        (
            AromaticityModel.permissive,
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(),
                ring_limits=RingLimits(),
            ),
        ),
    ],
)
def test_aromaticity_model_preset(preset, expected):
    assert preset() == expected


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            AromaticityModel.Clar(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
        ),
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            AromaticityModel.HueckelRule(
                scope=ElementScope.AllowList([Element("C")]),
                ring_limits=RingLimits(),
            ),
        ),
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(),
                ring_limits=RingLimits(min_ring_size=4),
            ),
        ),
        (
            AromaticityModel.Hmo(
                scope=ElementScope.Any(), stabilization_threshold=0.375
            ),
            AromaticityModel.Hmo(
                scope=ElementScope.AllowList([Element("C")]),
                stabilization_threshold=0.375,
            ),
        ),
        (
            AromaticityModel.Hmo(
                scope=ElementScope.Any(), stabilization_threshold=0.375
            ),
            AromaticityModel.Hmo(
                scope=ElementScope.Any(), stabilization_threshold=0.5
            ),
        ),
        (
            AromaticityModel.Clar(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            AromaticityModel.Clar(
                scope=ElementScope.AllowList([Element("C")]),
                ring_limits=RingLimits(),
            ),
        ),
        (
            AromaticityModel.Clar(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            AromaticityModel.Clar(
                scope=ElementScope.Any(),
                ring_limits=RingLimits(min_ring_size=6),
            ),
        ),
    ],
)
def test_aromaticity_model_equality(left, right):
    assert left != right


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(),
                ring_limits=RingLimits(
                    min_ring_size=4,
                    max_ring_size=18,
                    include_fused=False,
                    max_fused_combination=3,
                    max_fused_search=2_000,
                ),
            ),
            "AromaticityModel.HueckelRule(scope=ElementScope.Any(), "
            "ring_limits=RingLimits(min_ring_size=4, max_ring_size=18, "
            "include_fused=False, max_fused_combination=3, "
            "max_fused_search=2000))",
        ),
        (
            AromaticityModel.Hmo(
                scope=ElementScope.AllowList([Element("C"), Element("N")]),
                stabilization_threshold=0.375,
            ),
            "AromaticityModel.Hmo(scope=ElementScope.AllowList([Element('C'), "
            "Element('N')]), stabilization_threshold=0.375)",
        ),
        (
            AromaticityModel.Clar(
                scope=ElementScope.AllowList([Element("C")]),
                ring_limits=RingLimits(
                    min_ring_size=6,
                    max_ring_size=14,
                    max_fused_combination=4,
                    max_fused_search=1_500,
                ),
            ),
            "AromaticityModel.Clar(scope=ElementScope.AllowList([Element('C')]), "
            "ring_limits=RingLimits(min_ring_size=6, max_ring_size=14, "
            "include_fused=True, max_fused_combination=4, "
            "max_fused_search=1500))",
        ),
    ],
)
def test_aromaticity_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("model", "field", "value"),
    [
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            "scope",
            ElementScope.AllowList([Element("C")]),
        ),
        (
            AromaticityModel.HueckelRule(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            "ring_limits",
            RingLimits(min_ring_size=4),
        ),
        (
            AromaticityModel.Hmo(
                scope=ElementScope.Any(), stabilization_threshold=0.375
            ),
            "stabilization_threshold",
            0.5,
        ),
        (
            AromaticityModel.Clar(
                scope=ElementScope.Any(), ring_limits=RingLimits()
            ),
            "ring_limits",
            RingLimits(min_ring_size=6),
        ),
    ],
)
def test_aromaticity_model_mutation(model, field, value):
    with pytest.raises(AttributeError):
        setattr(model, field, value)


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (InconsistencyPolicy.Keep, InconsistencyPolicy.Keep, True),
        (InconsistencyPolicy.Strip, InconsistencyPolicy.Strip, True),
        (InconsistencyPolicy.Error, InconsistencyPolicy.Error, True),
        (InconsistencyPolicy.Keep, InconsistencyPolicy.Strip, False),
        (InconsistencyPolicy.Strip, InconsistencyPolicy.Error, False),
        (InconsistencyPolicy.Error, InconsistencyPolicy.Keep, False),
    ],
)
def test_inconsistency_policy_equality(left, right, expected):
    assert (left == right) is expected


def test_inconsistency_policy_hash():
    policies = {
        InconsistencyPolicy.Keep: "keep",
        InconsistencyPolicy.Strip: "strip",
        InconsistencyPolicy.Error: "error",
    }

    assert policies[InconsistencyPolicy.Keep] == "keep"
    assert policies[InconsistencyPolicy.Strip] == "strip"
    assert policies[InconsistencyPolicy.Error] == "error"


@pytest.mark.parametrize(
    ("policy", "expected"),
    [
        (InconsistencyPolicy.Keep, "InconsistencyPolicy.Keep"),
        (InconsistencyPolicy.Strip, "InconsistencyPolicy.Strip"),
        (InconsistencyPolicy.Error, "InconsistencyPolicy.Error"),
    ],
)
def test_inconsistency_policy_repr(policy, expected):
    assert repr(policy) == expected


@pytest.mark.parametrize(
    "policy",
    [
        InconsistencyPolicy.Keep,
        InconsistencyPolicy.Strip,
        InconsistencyPolicy.Error,
    ],
)
def test_inconsistency_policy_mutation(policy):
    with pytest.raises(AttributeError):
        policy.value = "changed"


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
    assert model.max_iterations == 16
    assert model.inconsistency == InconsistencyPolicy.Error
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
        max_iterations=8,
        inconsistency=InconsistencyPolicy.Strip,
    )
    kind_models.clear()

    assert model.kind_models == expected
    assert model.para_stereo is True
    assert model.max_iterations == 8
    assert model.inconsistency == InconsistencyPolicy.Strip
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
        StereoModel({}, False, 16, InconsistencyPolicy.Error)


@pytest.mark.parametrize("max_iterations", [-1, 1 << 64])
def test_stereo_model_new_integer_error(max_iterations):
    with pytest.raises(OverflowError):
        StereoModel(
            kind_models={},
            para_stereo=False,
            max_iterations=max_iterations,
            inconsistency=InconsistencyPolicy.Error,
        )


def test_stereo_model_max_iterations():
    model = StereoModel(
        kind_models={},
        para_stereo=False,
        max_iterations=0,
        inconsistency=InconsistencyPolicy.Error,
    )

    assert model.max_iterations == 0


@pytest.mark.parametrize(
    "other",
    [
        StereoModel(
            kind_models={},
            para_stereo=False,
            max_iterations=16,
            inconsistency=InconsistencyPolicy.Error,
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
            max_iterations=16,
            inconsistency=InconsistencyPolicy.Error,
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
            para_stereo=False,
            max_iterations=8,
            inconsistency=InconsistencyPolicy.Error,
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
            para_stereo=False,
            max_iterations=16,
            inconsistency=InconsistencyPolicy.Keep,
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
                max_iterations=8,
                inconsistency=InconsistencyPolicy.Keep,
            ),
            "StereoModel(kind_models={StereoKind.Tetrahedral: "
            "StereoKindModel(scope=ElementScope.Any(), fluxionality=False), "
            "StereoKind.Octahedral: StereoKindModel(scope="
            "ElementScope.AllowList([Element('Fe')]), fluxionality=True)}, "
            "para_stereo=True, max_iterations=8, "
            "inconsistency=InconsistencyPolicy.Keep)",
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
        ("max_iterations", 8),
        ("inconsistency", InconsistencyPolicy.Keep),
    ],
)
def test_stereo_model_mutation(field, value):
    model = StereoModel.default()

    with pytest.raises(AttributeError):
        setattr(model, field, value)


def test_chemistry_model_default():
    model = ChemistryModel.default()

    assert model.valence == ValenceModel.AtomTyping(
        registry=AtomTypeRegistry.default()
    )
    assert model.aromaticity == AromaticityModel.daylight()
    assert model.stereo == StereoModel.default()
    assert model == ChemistryModel.default()


def test_chemistry_model_new():
    valence = ValenceModel.Counts(
        table=ValenceTable(
            entries={
                Element("C"): ValenceEntry(target_covalences=[4]),
                Element("O"): ValenceEntry(target_covalences=[2]),
            }
        )
    )
    aromaticity = AromaticityModel.Clar(
        scope=ElementScope.AllowList([Element("C")]),
        ring_limits=RingLimits(min_ring_size=6),
    )
    stereo = StereoModel(
        kind_models=StereoModel.default().kind_models,
        para_stereo=True,
        max_iterations=8,
        inconsistency=InconsistencyPolicy.Strip,
    )
    model = ChemistryModel(
        valence=valence,
        aromaticity=aromaticity,
        stereo=stereo,
    )

    assert model.valence == valence
    assert model.valence is not valence
    assert model.aromaticity == aromaticity
    assert model.aromaticity is not aromaticity
    assert model.stereo == stereo
    assert model.stereo is not stereo


def test_chemistry_model_new_error():
    with pytest.raises(TypeError):
        ChemistryModel(
            ValenceModel.AtomTyping(registry=AtomTypeRegistry.default()),
            AromaticityModel.daylight(),
            StereoModel.default(),
        )


@pytest.mark.parametrize(
    "other",
    [
        ChemistryModel(
            valence=ValenceModel.Counts(table=ValenceTable(entries={})),
            aromaticity=AromaticityModel.daylight(),
            stereo=StereoModel.default(),
        ),
        ChemistryModel(
            valence=ValenceModel.AtomTyping(
                registry=AtomTypeRegistry.default()
            ),
            aromaticity=AromaticityModel.permissive(),
            stereo=StereoModel.default(),
        ),
        ChemistryModel(
            valence=ValenceModel.AtomTyping(
                registry=AtomTypeRegistry.default()
            ),
            aromaticity=AromaticityModel.daylight(),
            stereo=StereoModel(
                kind_models={},
                para_stereo=False,
                max_iterations=16,
                inconsistency=InconsistencyPolicy.Error,
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
                valence=ValenceModel.Counts(
                    table=ValenceTable(
                        entries={
                            Element("C"): ValenceEntry(target_covalences=[4])
                        }
                    )
                ),
                aromaticity=AromaticityModel.Hmo(
                    scope=ElementScope.Any(),
                    stabilization_threshold=0.375,
                ),
                stereo=StereoModel(
                    kind_models=StereoModel.default().kind_models,
                    para_stereo=True,
                    max_iterations=8,
                    inconsistency=InconsistencyPolicy.Keep,
                ),
            ),
            "ChemistryModel(valence=ValenceModel.Counts(table="
            "ValenceTable(entries={Element('C'): ValenceEntry("
            "target_covalences=[4], aromatic_valences=[])})), aromaticity="
            "AromaticityModel.Hmo(scope=ElementScope.Any(), "
            "stabilization_threshold=0.375), stereo=StereoModel(kind_models={"
            "StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any(), "
            "fluxionality=False), StereoKind.CisTrans: StereoKindModel(scope="
            "ElementScope.Any(), fluxionality=False)}, para_stereo=True, "
            "max_iterations=8, inconsistency=InconsistencyPolicy.Keep))",
        ),
    ],
)
def test_chemistry_model_repr(model, expected):
    assert repr(model) == expected


@pytest.mark.parametrize(
    ("field", "value"),
    [
        (
            "valence",
            ValenceModel.Counts(table=ValenceTable(entries={})),
        ),
        ("aromaticity", AromaticityModel.permissive()),
        (
            "stereo",
            StereoModel(
                kind_models={},
                para_stereo=False,
                max_iterations=16,
                inconsistency=InconsistencyPolicy.Error,
            ),
        ),
    ],
)
def test_chemistry_model_mutation(field, value):
    model = ChemistryModel.default()

    with pytest.raises(AttributeError):
        setattr(model, field, value)
