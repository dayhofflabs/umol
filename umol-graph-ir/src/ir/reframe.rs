//! Acting on a frame-relative value by a relabeling of its frame positions.
//!
//! `FrameAction` is the action itself; `find_reframed` searches the actions carrying one frame to
//! another. The family-level frame quotient is `Reframe` in [`super::traits`].

use std::ops::ControlFlow;

use umol_perm::Permutation;

/// A value read against an ordered participant frame, restatable under a relabeling of that
/// frame's positions. `None` when the permutation is not an admissible action for the value.
pub trait FrameAction: Sized {
    fn reframe_by(self, permutation: Permutation) -> Option<Self>;
}

/// The first admissible restatement of `value` carrying frame `from` to frame `to` for which
/// `select` yields a value, together with nothing else.
///
/// Equal repeated participants leave several restatements. They differ by a stabilizer element of
/// `to`, so the values they produce denote one arrangement and the first `select` accepts stands.
/// `select` receives the action as well, for a caller that must carry a second value under the
/// same one. An action the value declines is skipped rather than ending the walk. `None` when the
/// two frames do not hold the same participants, or when `select` accepts none of the
/// restatements.
pub fn find_reframed<L, T, B, F>(value: &T, from: &[L], to: &[L], mut select: F) -> Option<B>
where
    L: Eq,
    T: FrameAction + Clone,
    F: FnMut(Permutation, T) -> Option<B>,
{
    match Permutation::visit_between(from, to, |action| {
        match value
            .clone()
            .reframe_by(action)
            .and_then(|restated| select(action, restated))
        {
            Some(found) => ControlFlow::Break(found),
            None => ControlFlow::Continue(()),
        }
    }) {
        ControlFlow::Break(found) => Some(found),
        ControlFlow::Continue(()) => None,
    }
}
