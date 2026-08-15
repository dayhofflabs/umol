//! Per-atom completion carrier for the resolve pipeline: each entry maps an
//! atom to the disjunction of completions surviving a resolution phase.

use std::collections::BTreeMap;

use smallvec::SmallVec;
use umol_graph_ir::ir::{AtomForm, AtomId};

/// Per-atom disjunctions of surviving completions, in ascending atom order.
///
/// An entry is non-empty — an empty candidate set is a contradiction, not an
/// entry. This is the one representation invariant.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomCompletions {
    entries: BTreeMap<AtomId, SmallVec<[AtomForm; 1]>>,
}

impl AtomCompletions {
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the surviving completions of `atom`, replacing any previous entry.
    ///
    /// # Panics
    ///
    /// Panics if `completions` is empty.
    pub fn insert(&mut self, atom: AtomId, completions: SmallVec<[AtomForm; 1]>) {
        assert!(
            !completions.is_empty(),
            "an empty candidate set is a contradiction, not an entry"
        );
        self.entries.insert(atom, completions);
    }

    pub fn get(&self, atom: AtomId) -> Option<&[AtomForm]> {
        self.entries.get(&atom).map(|entry| entry.as_slice())
    }

    pub fn remove(&mut self, atom: AtomId) -> Option<SmallVec<[AtomForm; 1]>> {
        self.entries.remove(&atom)
    }

    /// Entries in ascending atom order.
    pub fn iter(&self) -> impl Iterator<Item = (AtomId, &[AtomForm])> {
        self.entries
            .iter()
            .map(|(&atom, entry)| (atom, entry.as_slice()))
    }

    /// Number of atoms carrying an entry.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl FromIterator<(AtomId, SmallVec<[AtomForm; 1]>)> for AtomCompletions {
    /// Collects entries through `insert`: later duplicates replace earlier
    /// ones.
    ///
    /// # Panics
    ///
    /// Panics if any entry is empty.
    fn from_iter<I: IntoIterator<Item = (AtomId, SmallVec<[AtomForm; 1]>)>>(iter: I) -> Self {
        let mut completions = Self::new();
        for (atom, entry) in iter {
            completions.insert(atom, entry);
        }
        completions
    }
}

/// Resolution report: a descriptive record with no invariant.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ResolveReport {
    /// Per-atom plural survivors; empty under a `Determined` solution.
    pub unresolved: AtomCompletions,
    /// Atoms whose completion was selected by the valence-preference
    /// tie-break; sorted and deduplicated.
    pub tie_breaks: Vec<AtomId>,
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use rstest::rstest;
    use smallvec::{smallvec, SmallVec};
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomForm, AtomId};

    use super::AtomCompletions;

    #[rstest]
    #[case::single_into_empty(
        AtomCompletions::new(),
        AtomId(0),
        smallvec![AtomForm::from_element(Element::C)],
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        }
    )]
    #[case::plural_alongside_existing(
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        },
        AtomId(3),
        smallvec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::N).with_implicit_hydrogens(1),
        ],
        AtomCompletions {
            entries: BTreeMap::from([
                (
                    AtomId(0),
                    SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
                ),
                (
                    AtomId(3),
                    SmallVec::from_vec(vec![
                        AtomForm::from_element(Element::N),
                        AtomForm::from_element(Element::N).with_implicit_hydrogens(1),
                    ]),
                ),
            ]),
        }
    )]
    #[case::replaces_previous_entry(
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        },
        AtomId(0),
        smallvec![AtomForm::from_element(Element::O)],
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::O)]),
            )]),
        }
    )]
    fn test_atom_completions_insert(
        #[case] mut completions: AtomCompletions,
        #[case] atom: AtomId,
        #[case] entry: SmallVec<[AtomForm; 1]>,
        #[case] expected: AtomCompletions,
    ) {
        completions.insert(atom, entry);
        assert_eq!(completions, expected);
    }

    #[rstest]
    #[should_panic(expected = "an empty candidate set is a contradiction, not an entry")]
    fn test_atom_completions_insert_error() {
        AtomCompletions::new().insert(AtomId(0), SmallVec::new());
    }

    #[rstest]
    #[case::empty(vec![], AtomCompletions::new())]
    #[case::entries(
        vec![
            (AtomId(4), smallvec![AtomForm::from_element(Element::N)]),
            (
                AtomId(0),
                smallvec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
            ),
        ],
        AtomCompletions {
            entries: BTreeMap::from([
                (
                    AtomId(0),
                    SmallVec::from_vec(vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::O),
                    ]),
                ),
                (
                    AtomId(4),
                    SmallVec::from_vec(vec![AtomForm::from_element(Element::N)]),
                ),
            ]),
        }
    )]
    #[case::duplicate_replaces(
        vec![
            (AtomId(0), smallvec![AtomForm::from_element(Element::C)]),
            (AtomId(0), smallvec![AtomForm::from_element(Element::O)]),
        ],
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::O)]),
            )]),
        }
    )]
    fn test_atom_completions_from_iter(
        #[case] entries: Vec<(AtomId, SmallVec<[AtomForm; 1]>)>,
        #[case] expected: AtomCompletions,
    ) {
        assert_eq!(AtomCompletions::from_iter(entries), expected);
    }

    #[rstest]
    #[should_panic(expected = "an empty candidate set is a contradiction, not an entry")]
    fn test_atom_completions_from_iter_error() {
        AtomCompletions::from_iter([(AtomId(0), SmallVec::new())]);
    }

    #[rstest]
    #[case::present(
        AtomId(1),
        Some(vec![AtomForm::from_element(Element::C)])
    )]
    #[case::absent(AtomId(2), None)]
    fn test_atom_completions_get(#[case] atom: AtomId, #[case] expected: Option<Vec<AtomForm>>) {
        let completions = AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(1),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        };
        assert_eq!(completions.get(atom), expected.as_deref());
    }

    #[rstest]
    #[case::present(
        AtomId(1),
        Some(smallvec![AtomForm::from_element(Element::C)]),
        AtomCompletions::new()
    )]
    #[case::absent(
        AtomId(2),
        None,
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(1),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        }
    )]
    fn test_atom_completions_remove(
        #[case] atom: AtomId,
        #[case] expected: Option<SmallVec<[AtomForm; 1]>>,
        #[case] expected_remaining: AtomCompletions,
    ) {
        let mut completions = AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(1),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        };
        assert_eq!(completions.remove(atom), expected);
        assert_eq!(completions, expected_remaining);
    }

    #[rstest]
    fn test_atom_completions_iter() {
        let completions = AtomCompletions {
            entries: BTreeMap::from([
                (
                    AtomId(4),
                    SmallVec::from_vec(vec![AtomForm::from_element(Element::N)]),
                ),
                (
                    AtomId(0),
                    SmallVec::from_vec(vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::O),
                    ]),
                ),
            ]),
        };
        assert_eq!(
            completions.iter().collect::<Vec<_>>(),
            vec![
                (
                    AtomId(0),
                    [
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::O),
                    ]
                    .as_slice(),
                ),
                (AtomId(4), [AtomForm::from_element(Element::N)].as_slice(),),
            ]
        );
    }

    #[rstest]
    #[case::empty(AtomCompletions::new(), 0)]
    #[case::two_atoms(
        AtomCompletions {
            entries: BTreeMap::from([
                (
                    AtomId(0),
                    SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
                ),
                (
                    AtomId(1),
                    SmallVec::from_vec(vec![
                        AtomForm::from_element(Element::N),
                        AtomForm::from_element(Element::O),
                    ]),
                ),
            ]),
        },
        2
    )]
    fn test_atom_completions_len(#[case] completions: AtomCompletions, #[case] expected: usize) {
        assert_eq!(completions.len(), expected);
    }

    #[rstest]
    #[case::empty(AtomCompletions::new(), true)]
    #[case::populated(
        AtomCompletions {
            entries: BTreeMap::from([(
                AtomId(0),
                SmallVec::from_vec(vec![AtomForm::from_element(Element::C)]),
            )]),
        },
        false
    )]
    fn test_atom_completions_is_empty(
        #[case] completions: AtomCompletions,
        #[case] expected: bool,
    ) {
        assert_eq!(completions.is_empty(), expected);
    }
}
