import inspect

import pytest
import umol


@pytest.mark.parametrize(
    ("constructor", "expected"),
    [
        (umol.RefinementRounds.Fixed, "(*, rounds)"),
        (umol.RefinementRounds.ToFixpoint, "()"),
        (
            umol.HashedFingerprintConfig.Morgan,
            "(*, radius=2, ring_config=Ellipsis)",
        ),
        (
            umol.HashedFingerprintConfig.Ecfp,
            "(*, radius, hashing_scheme=Ellipsis, ring_config=Ellipsis)",
        ),
        (
            umol.HashedFingerprintConfig.Wl,
            "(*, rounds, hashing_scheme=Ellipsis)",
        ),
        (
            umol.ReactionCombinedFingerprintConfig.Difference,
            "(*, molecule)",
        ),
        (
            umol.ReactionCombinedFingerprintConfig.DisjointUnion,
            "(*, molecule)",
        ),
        (
            umol.PatternFingerprintConfig,
            "(*, width=2048, match_algorithm=Ellipsis, "
            "subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
        ),
        (
            umol.StructuralFingerprintConfig,
            "(*, max_bonds, subgraph_enumeration_algorithm=Ellipsis, "
            "automorphism_algorithm=Ellipsis)",
        ),
        (
            umol.RingConfig,
            "(*, simple_cycle_algorithm=None, relevant_cycle_algorithm=None)",
        ),
        (
            umol.AromaticityConfig,
            "(*, ring_config=Ellipsis, connected_components_algorithm=Ellipsis, "
            "maximum_independent_set_algorithm=Ellipsis)",
        ),
        (
            umol.AromaticityResolveConfig,
            "(*, perception=Ellipsis, aromatic_valence_failure=Ellipsis, "
            "aromatic_system_failure=Ellipsis, "
            "aromatic_valence_mismatch=Ellipsis, "
            "aromatic_bond_constraint_mismatch=Ellipsis, "
            "reset_aromatic_valence=False)",
        ),
        (
            umol.StereoResolveConfig,
            "(*, tetrahedral_stereo_failure=Ellipsis, "
            "stereo_atom_failure=Ellipsis, "
            "tetrahedral_stereo_mismatch=Ellipsis, "
            "cis_trans_stereo_failure=Ellipsis, "
            "stereo_bond_failure=Ellipsis, "
            "cis_trans_stereo_mismatch=Ellipsis, "
            "reset_stereo_constraints=False)",
        ),
        (umol.ResolveConfig, "(*, aromaticity, stereo)"),
        (
            umol.ConnectivityModel,
            "(*, allow_disconnected, allow_disconnected_dative, "
            "allow_disconnected_aromatic, allow_disconnected_multicenter, "
            "allow_disconnected_noncovalent, allow_disconnected_stereo_atom, "
            "allow_disconnected_stereo_bond, allow_disconnected_constraints)",
        ),
        (
            umol.ChemistryModel,
            "(*, connectivity, valence, aromaticity, stereo)",
        ),
        (
            umol.ReactionApplicationConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
        ),
        (
            umol.ReactionCompositionConfig,
            "(*, common_subgraph_enumeration_algorithm=Ellipsis)",
        ),
        (
            umol.CanonicalizeConfig,
            "(*, automorphism_algorithm=Ellipsis)",
        ),
        (
            umol.SubstructureSearchConfig,
            "(*, match_algorithm=Ellipsis, subgraph_isomorphism_algorithm=Ellipsis, "
            "relevant_cycle_algorithm=Ellipsis)",
        ),
    ],
)
def test_public_class_constructor_signature(constructor, expected):
    assert str(inspect.signature(constructor)) == expected
