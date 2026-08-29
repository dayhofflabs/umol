//! Complete participant-frame actions for overlay aggregates.

use std::collections::BTreeMap;

use umol_perm::{ClassKey, DynPermutation, Permutation};

use super::id::{
    AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
    StereoBondId,
};

trait FramePermutation: Clone {
    fn degree(&self) -> usize;
    fn identity(&self) -> Self;
    fn inverse(&self) -> Self;
    fn compose(&self, other: &Self) -> Option<Self>;
}

impl FramePermutation for DynPermutation {
    fn degree(&self) -> usize {
        self.degree()
    }

    fn identity(&self) -> Self {
        Self::identity(self.degree())
    }

    fn inverse(&self) -> Self {
        self.inverse()
    }

    fn compose(&self, other: &Self) -> Option<Self> {
        self.compose(other)
    }
}

impl FramePermutation for Permutation {
    fn degree(&self) -> usize {
        (*self).degree()
    }

    fn identity(&self) -> Self {
        Self::identity((*self).degree())
    }

    fn inverse(&self) -> Self {
        (*self).inverse()
    }

    fn compose(&self, other: &Self) -> Option<Self> {
        ((*self).degree() == (*other).degree()).then(|| (*self).compose(*other))
    }
}

macro_rules! overlay_frame_action {
    (
        $(#[$meta:meta])*
        $name:ident, $id:ty, $action:ty, $allows:expr
    ) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name(BTreeMap<$id, $action>);

        impl $name {
            pub(crate) fn from_vec(actions: Vec<$action>) -> Option<Self> {
                Self::from_action_map(
                    actions
                        .into_iter()
                        .enumerate()
                        .map(|(index, action)| (<$id>::from(index), action))
                        .collect(),
                )
            }

            pub(crate) fn from_action_map(actions: BTreeMap<$id, $action>) -> Option<Self> {
                actions.values().all($allows).then_some(Self(actions))
            }

            /// Number of entity ids in this exact action domain.
            pub fn count(&self) -> usize {
                self.0.len()
            }

            /// Whether the action domain covers `id`.
            pub fn contains(&self, id: $id) -> bool {
                self.0.contains_key(&id)
            }

            /// Typed entity ids in the action domain.
            pub fn ids(&self) -> impl ExactSizeIterator<Item = $id> + '_ {
                self.0.keys().copied()
            }

            /// Local participant-frame action for `id`, if the domain covers it.
            pub fn action(&self, id: $id) -> Option<&$action> {
                self.0.get(&id)
            }

            /// Identity action on exactly this typed id-and-degree domain.
            pub fn identity(&self) -> Self {
                Self(
                    self.0
                        .iter()
                        .map(|(&id, action)| (id, action.identity()))
                        .collect(),
                )
            }

            /// Componentwise inverse, preserving the exact typed id-and-degree domain.
            pub fn inverse(&self) -> Self {
                Self(
                    self.0
                        .iter()
                        .map(|(&id, action)| (id, action.inverse()))
                        .collect(),
                )
            }

            /// Componentwise composition on equal typed id-and-degree domains.
            pub fn compose(&self, other: &Self) -> Option<Self> {
                if self.0.keys().ne(other.0.keys()) {
                    return None;
                }
                self.0
                    .iter()
                    .zip(&other.0)
                    .map(|((&id, left), (&other_id, right))| {
                        debug_assert!(id == other_id);
                        (left.degree() == right.degree())
                            .then(|| left.compose(right))
                            .flatten()
                            .map(|action| (id, action))
                    })
                    .collect::<Option<BTreeMap<_, _>>>()
                    .map(Self)
            }
        }
    };
}

overlay_frame_action! {
    /// Complete participant-frame action for the dative-bond aggregate.
    DativeBondsFrameAction, DativeBondId, DynPermutation, |_| true
}

overlay_frame_action! {
    /// Complete participant-frame action for the aromatic-system aggregate.
    AromaticSystemsFrameAction, AromaticSystemId, DynPermutation, |_| true
}

overlay_frame_action! {
    /// Complete participant-frame action for the multicenter-bond aggregate.
    MulticenterBondsFrameAction, MulticenterBondId, DynPermutation, |_| true
}

overlay_frame_action! {
    /// Complete participant-frame action for the noncovalent-bond aggregate.
    NoncovalentBondsFrameAction, NoncovalentBondId, DynPermutation,
    |action: &DynPermutation| action.degree() == 2
}

overlay_frame_action! {
    /// Complete participant-frame action for the stereo-atom aggregate.
    StereoAtomsFrameAction, StereoAtomId, Permutation, |_| true
}

overlay_frame_action! {
    /// Complete participant-frame action for the stereo-bond aggregate.
    StereoBondsFrameAction, StereoBondId, Permutation,
    |action: &Permutation| ClassKey::CisTrans.space().allows(*action)
}

/// Complete participant-frame action over all six overlays.
///
/// Fields are private because this is an operation-issued witness. Producers establish the typed
/// domains and local action groups; consumers may use a covering witness and ignore entries they do
/// not reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OverlaysFrameAction {
    dative_bonds: DativeBondsFrameAction,
    aromatic_systems: AromaticSystemsFrameAction,
    multicenter_bonds: MulticenterBondsFrameAction,
    noncovalent_bonds: NoncovalentBondsFrameAction,
    stereo_atoms: StereoAtomsFrameAction,
    stereo_bonds: StereoBondsFrameAction,
}

impl OverlaysFrameAction {
    pub(crate) fn new(
        dative_bonds: DativeBondsFrameAction,
        aromatic_systems: AromaticSystemsFrameAction,
        multicenter_bonds: MulticenterBondsFrameAction,
        noncovalent_bonds: NoncovalentBondsFrameAction,
        stereo_atoms: StereoAtomsFrameAction,
        stereo_bonds: StereoBondsFrameAction,
    ) -> Self {
        Self {
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        }
    }

    pub fn dative_bonds(&self) -> &DativeBondsFrameAction {
        &self.dative_bonds
    }

    pub fn aromatic_systems(&self) -> &AromaticSystemsFrameAction {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &MulticenterBondsFrameAction {
        &self.multicenter_bonds
    }

    pub fn noncovalent_bonds(&self) -> &NoncovalentBondsFrameAction {
        &self.noncovalent_bonds
    }

    pub fn stereo_atoms(&self) -> &StereoAtomsFrameAction {
        &self.stereo_atoms
    }

    pub fn stereo_bonds(&self) -> &StereoBondsFrameAction {
        &self.stereo_bonds
    }

    /// Identity action on exactly this six-overlay domain.
    pub fn identity(&self) -> Self {
        Self {
            dative_bonds: self.dative_bonds.identity(),
            aromatic_systems: self.aromatic_systems.identity(),
            multicenter_bonds: self.multicenter_bonds.identity(),
            noncovalent_bonds: self.noncovalent_bonds.identity(),
            stereo_atoms: self.stereo_atoms.identity(),
            stereo_bonds: self.stereo_bonds.identity(),
        }
    }

    /// Componentwise inverse, preserving every entity-kind domain and local degree.
    pub fn inverse(&self) -> Self {
        Self {
            dative_bonds: self.dative_bonds.inverse(),
            aromatic_systems: self.aromatic_systems.inverse(),
            multicenter_bonds: self.multicenter_bonds.inverse(),
            noncovalent_bonds: self.noncovalent_bonds.inverse(),
            stereo_atoms: self.stereo_atoms.inverse(),
            stereo_bonds: self.stereo_bonds.inverse(),
        }
    }

    /// Componentwise composition when every component has the same typed id-and-degree domain.
    pub fn compose(&self, other: &Self) -> Option<Self> {
        Some(Self {
            dative_bonds: self.dative_bonds.compose(&other.dative_bonds)?,
            aromatic_systems: self.aromatic_systems.compose(&other.aromatic_systems)?,
            multicenter_bonds: self.multicenter_bonds.compose(&other.multicenter_bonds)?,
            noncovalent_bonds: self.noncovalent_bonds.compose(&other.noncovalent_bonds)?,
            stereo_atoms: self.stereo_atoms.compose(&other.stereo_atoms)?,
            stereo_bonds: self.stereo_bonds.compose(&other.stereo_bonds)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    fn dynamic(image: Vec<usize>) -> DynPermutation {
        DynPermutation::try_from(image).expect("case image is a permutation")
    }

    fn overlays_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction {
            dative_bonds: DativeBondsFrameAction::from_vec(vec![dynamic(vec![1, 0])])
                .expect("dative action is admissible"),
            aromatic_systems: AromaticSystemsFrameAction::from_vec(vec![dynamic(vec![2, 0, 1])])
                .expect("aromatic action is admissible"),
            multicenter_bonds: MulticenterBondsFrameAction::from_vec(vec![dynamic(vec![0, 1])])
                .expect("multicenter action is admissible"),
            noncovalent_bonds: NoncovalentBondsFrameAction::from_vec(vec![dynamic(vec![1, 0])])
                .expect("noncovalent action is admissible"),
            stereo_atoms: StereoAtomsFrameAction::from_vec(vec![Permutation::from_image(&[
                1, 0, 2, 3,
            ])])
            .expect("stereo-atom action is admissible"),
            stereo_bonds: StereoBondsFrameAction::from_vec(vec![Permutation::from_image(&[
                2, 3, 0, 1,
            ])])
            .expect("stereo-bond action is admissible"),
        }
    }

    #[rstest]
    fn test_aromatic_systems_frame_action_typed_lookup() {
        let action =
            AromaticSystemsFrameAction::from_vec(vec![dynamic(vec![1, 0]), dynamic(vec![2, 0, 1])])
                .expect("actions are admissible");

        assert_eq!(action.count(), 2);
        assert!(action.contains(AromaticSystemId(1)));
        assert_eq!(
            action.ids().collect::<Vec<_>>(),
            [AromaticSystemId(0), AromaticSystemId(1)]
        );
        assert_eq!(
            action
                .action(AromaticSystemId(1))
                .map(DynPermutation::image),
            Some([2, 0, 1].as_slice()),
        );
        assert_eq!(action.action(AromaticSystemId(2)), None);
    }

    #[rstest]
    fn test_aromatic_systems_frame_action_from_action_map() {
        let action = AromaticSystemsFrameAction::from_action_map(BTreeMap::from([
            (AromaticSystemId(2), dynamic(vec![1, 0])),
            (AromaticSystemId(7), dynamic(vec![2, 0, 1])),
        ]))
        .expect("actions are admissible");

        assert_eq!(
            action.ids().collect::<Vec<_>>(),
            [AromaticSystemId(2), AromaticSystemId(7)],
        );
        assert_eq!(
            action
                .action(AromaticSystemId(7))
                .map(DynPermutation::image),
            Some([2, 0, 1].as_slice()),
        );
        assert_eq!(action.action(AromaticSystemId(0)), None);
    }

    #[rstest]
    fn test_overlays_frame_action_identity_preserves_domain() {
        let action = overlays_frame_action();
        let identity = action.identity();

        assert_eq!(
            identity.dative_bonds().ids().collect::<Vec<_>>(),
            [DativeBondId(0)]
        );
        assert_eq!(
            identity.aromatic_systems().action(AromaticSystemId(0)),
            Some(&DynPermutation::identity(3)),
        );
        assert_eq!(
            identity.stereo_bonds().action(StereoBondId(0)),
            Some(&Permutation::identity(4)),
        );
    }

    #[rstest]
    fn test_overlays_frame_action_inverse_is_two_sided() {
        let action = overlays_frame_action();
        let inverse = action.inverse();

        assert_eq!(action.compose(&inverse), Some(action.identity()));
        assert_eq!(inverse.compose(&action), Some(action.identity()));
    }

    #[rstest]
    #[case::compatible(
        vec![dynamic(vec![1, 2, 0])],
        vec![dynamic(vec![2, 0, 1])],
        Some(vec![dynamic(vec![0, 1, 2])]),
    )]
    #[case::id_domain(
        vec![dynamic(vec![1, 0])],
        vec![dynamic(vec![1, 0]), dynamic(vec![0, 1])],
        None,
    )]
    #[case::degree(
        vec![dynamic(vec![1, 0])],
        vec![dynamic(vec![1, 2, 0])],
        None,
    )]
    fn test_dative_bonds_frame_action_compose(
        #[case] left: Vec<DynPermutation>,
        #[case] right: Vec<DynPermutation>,
        #[case] expected: Option<Vec<DynPermutation>>,
    ) {
        let left = DativeBondsFrameAction::from_vec(left).expect("left actions are admissible");
        let right = DativeBondsFrameAction::from_vec(right).expect("right actions are admissible");

        assert_eq!(
            left.compose(&right).map(|action| {
                action
                    .ids()
                    .map(|id| action.action(id).expect("id came from domain").clone())
                    .collect::<Vec<_>>()
            }),
            expected,
        );
    }

    #[rstest]
    fn test_noncovalent_bonds_frame_action_reject_wrong_degree() {
        assert_eq!(
            NoncovalentBondsFrameAction::from_vec(vec![dynamic(vec![0, 1, 2])]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_bonds_frame_action_reject_cross_endpoint_action() {
        assert_eq!(
            StereoBondsFrameAction::from_vec(vec![Permutation::from_image(&[1, 2, 0, 3])]),
            None,
        );
    }
}
