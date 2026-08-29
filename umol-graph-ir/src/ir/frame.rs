//! Complete participant-frame actions for overlay aggregates.

use std::collections::BTreeMap;

use umol_perm::{ClassKey, DynPermutation, Permutation};

use super::id::{
    AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
    StereoBondId,
};

trait LocalFrameAction: Clone {
    fn degree(&self) -> usize;
    fn identity(&self) -> Self;
    fn inverse(&self) -> Self;
    fn compose(&self, other: &Self) -> Option<Self>;
}

impl LocalFrameAction for DynPermutation {
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

impl LocalFrameAction for Permutation {
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

#[derive(Clone, Debug, PartialEq, Eq)]
struct FamilyFrameActions<I, A> {
    actions: BTreeMap<I, A>,
}

impl<I, A> FamilyFrameActions<I, A>
where
    I: Copy + Ord + From<usize>,
    A: LocalFrameAction,
{
    fn from_dense(actions: Vec<A>) -> Self {
        Self {
            actions: actions
                .into_iter()
                .enumerate()
                .map(|(index, action)| (I::from(index), action))
                .collect(),
        }
    }

    fn count(&self) -> usize {
        self.actions.len()
    }

    fn contains(&self, id: I) -> bool {
        self.actions.contains_key(&id)
    }

    fn ids(&self) -> impl ExactSizeIterator<Item = I> + '_ {
        self.actions.keys().copied()
    }

    fn action(&self, id: I) -> Option<&A> {
        self.actions.get(&id)
    }

    fn identity(&self) -> Self {
        Self {
            actions: self
                .actions
                .iter()
                .map(|(&id, action)| (id, action.identity()))
                .collect(),
        }
    }

    fn inverse(&self) -> Self {
        Self {
            actions: self
                .actions
                .iter()
                .map(|(&id, action)| (id, action.inverse()))
                .collect(),
        }
    }

    fn compose(&self, other: &Self) -> Option<Self> {
        if self.actions.keys().ne(other.actions.keys()) {
            return None;
        }
        let actions = self
            .actions
            .iter()
            .zip(&other.actions)
            .map(|((&id, left), (&other_id, right))| {
                debug_assert!(id == other_id);
                (left.degree() == right.degree())
                    .then(|| left.compose(right))
                    .flatten()
                    .map(|action| (id, action))
            })
            .collect::<Option<BTreeMap<_, _>>>()?;
        Some(Self { actions })
    }
}

macro_rules! family_frame_actions {
    (
        $(#[$meta:meta])*
        $name:ident, $id:ty, $action:ty, $allows:expr
    ) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name(FamilyFrameActions<$id, $action>);

        impl $name {
            pub(crate) fn from_dense(actions: Vec<$action>) -> Option<Self> {
                actions
                    .iter()
                    .all($allows)
                    .then(|| Self(FamilyFrameActions::from_dense(actions)))
            }

            /// Number of entity ids in this exact action domain.
            pub fn count(&self) -> usize {
                self.0.count()
            }

            /// Whether the action domain covers `id`.
            pub fn contains(&self, id: $id) -> bool {
                self.0.contains(id)
            }

            /// Typed entity ids in the action domain.
            pub fn ids(&self) -> impl ExactSizeIterator<Item = $id> + '_ {
                self.0.ids()
            }

            /// Local participant-frame action for `id`, if the domain covers it.
            pub fn action(&self, id: $id) -> Option<&$action> {
                self.0.action(id)
            }

            /// Identity action on exactly this typed id-and-degree domain.
            pub fn identity(&self) -> Self {
                Self(self.0.identity())
            }

            /// Componentwise inverse, preserving the exact typed id-and-degree domain.
            pub fn inverse(&self) -> Self {
                Self(self.0.inverse())
            }

            /// Componentwise composition on equal typed id-and-degree domains.
            pub fn compose(&self, other: &Self) -> Option<Self> {
                self.0.compose(&other.0).map(Self)
            }
        }
    };
}

family_frame_actions! {
    /// Complete participant-frame action for one dative-bond family aggregate.
    DativeBondFrameActions, DativeBondId, DynPermutation, |_| true
}

family_frame_actions! {
    /// Complete participant-frame action for one aromatic-system family aggregate.
    AromaticSystemFrameActions, AromaticSystemId, DynPermutation, |_| true
}

family_frame_actions! {
    /// Complete participant-frame action for one multicenter-bond family aggregate.
    MulticenterBondFrameActions, MulticenterBondId, DynPermutation, |_| true
}

family_frame_actions! {
    /// Complete participant-frame action for one noncovalent-bond family aggregate.
    NoncovalentBondFrameActions, NoncovalentBondId, DynPermutation,
    |action: &DynPermutation| action.degree() == 2
}

family_frame_actions! {
    /// Complete participant-frame action for one stereo-atom family aggregate.
    StereoAtomFrameActions, StereoAtomId, Permutation, |_| true
}

family_frame_actions! {
    /// Complete participant-frame action for one stereo-bond family aggregate.
    StereoBondFrameActions, StereoBondId, Permutation,
    |action: &Permutation| ClassKey::CisTrans.space().allows(*action)
}

/// Complete participant-frame action over all six overlay entity families.
///
/// Fields are private because this is an operation-issued witness. Producers establish the typed
/// domains and local action groups; consumers may use a covering witness and ignore entries they do
/// not reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OverlayFrameActions {
    dative_bonds: DativeBondFrameActions,
    aromatic_systems: AromaticSystemFrameActions,
    multicenter_bonds: MulticenterBondFrameActions,
    noncovalent_bonds: NoncovalentBondFrameActions,
    stereo_atoms: StereoAtomFrameActions,
    stereo_bonds: StereoBondFrameActions,
}

impl OverlayFrameActions {
    pub fn dative_bonds(&self) -> &DativeBondFrameActions {
        &self.dative_bonds
    }

    pub fn aromatic_systems(&self) -> &AromaticSystemFrameActions {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &MulticenterBondFrameActions {
        &self.multicenter_bonds
    }

    pub fn noncovalent_bonds(&self) -> &NoncovalentBondFrameActions {
        &self.noncovalent_bonds
    }

    pub fn stereo_atoms(&self) -> &StereoAtomFrameActions {
        &self.stereo_atoms
    }

    pub fn stereo_bonds(&self) -> &StereoBondFrameActions {
        &self.stereo_bonds
    }

    /// Identity action on exactly this six-family domain.
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

    /// Componentwise inverse, preserving every family domain and local degree.
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

    /// Componentwise composition when every family has the same typed id-and-degree domain.
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

    fn overlay_actions() -> OverlayFrameActions {
        OverlayFrameActions {
            dative_bonds: DativeBondFrameActions::from_dense(vec![dynamic(vec![1, 0])])
                .expect("dative action is admissible"),
            aromatic_systems: AromaticSystemFrameActions::from_dense(vec![dynamic(vec![2, 0, 1])])
                .expect("aromatic action is admissible"),
            multicenter_bonds: MulticenterBondFrameActions::from_dense(vec![dynamic(vec![0, 1])])
                .expect("multicenter action is admissible"),
            noncovalent_bonds: NoncovalentBondFrameActions::from_dense(vec![dynamic(vec![1, 0])])
                .expect("noncovalent action is admissible"),
            stereo_atoms: StereoAtomFrameActions::from_dense(vec![Permutation::from_image(&[
                1, 0, 2, 3,
            ])])
            .expect("stereo-atom action is admissible"),
            stereo_bonds: StereoBondFrameActions::from_dense(vec![Permutation::from_image(&[
                2, 3, 0, 1,
            ])])
            .expect("stereo-bond action is admissible"),
        }
    }

    #[rstest]
    fn test_aromatic_system_frame_actions_typed_lookup() {
        let actions = AromaticSystemFrameActions::from_dense(vec![
            dynamic(vec![1, 0]),
            dynamic(vec![2, 0, 1]),
        ])
        .expect("actions are admissible");

        assert_eq!(actions.count(), 2);
        assert!(actions.contains(AromaticSystemId(1)));
        assert_eq!(
            actions.ids().collect::<Vec<_>>(),
            [AromaticSystemId(0), AromaticSystemId(1)]
        );
        assert_eq!(
            actions
                .action(AromaticSystemId(1))
                .map(DynPermutation::image),
            Some([2, 0, 1].as_slice()),
        );
        assert_eq!(actions.action(AromaticSystemId(2)), None);
    }

    #[rstest]
    fn test_overlay_frame_actions_identity_preserves_domain() {
        let actions = overlay_actions();
        let identity = actions.identity();

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
    fn test_overlay_frame_actions_inverse_is_two_sided() {
        let actions = overlay_actions();
        let inverse = actions.inverse();

        assert_eq!(actions.compose(&inverse), Some(actions.identity()));
        assert_eq!(inverse.compose(&actions), Some(actions.identity()));
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
    fn test_dative_bond_frame_actions_compose(
        #[case] left: Vec<DynPermutation>,
        #[case] right: Vec<DynPermutation>,
        #[case] expected: Option<Vec<DynPermutation>>,
    ) {
        let left = DativeBondFrameActions::from_dense(left).expect("left actions are admissible");
        let right =
            DativeBondFrameActions::from_dense(right).expect("right actions are admissible");

        assert_eq!(
            left.compose(&right).map(|actions| {
                actions
                    .ids()
                    .map(|id| actions.action(id).expect("id came from domain").clone())
                    .collect::<Vec<_>>()
            }),
            expected,
        );
    }

    #[rstest]
    fn test_noncovalent_bond_frame_actions_reject_wrong_degree() {
        assert_eq!(
            NoncovalentBondFrameActions::from_dense(vec![dynamic(vec![0, 1, 2])]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_bond_frame_actions_reject_cross_endpoint_action() {
        assert_eq!(
            StereoBondFrameActions::from_dense(vec![Permutation::from_image(&[1, 2, 0, 3])]),
            None,
        );
    }
}
