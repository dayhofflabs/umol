//! Binary coset actions from permutation parity and the two Ct endpoint blocks.
//! No coset-space, orbit, or graph-symmetry algorithm is used by the reference.

use std::array;

use rstest::rstest;
use umol_graph_ir::ir::StereoKind;
use umol_perm::PermutationGroup;

pub(super) fn binary_coset_action(kind: StereoKind, coset: u32, image: &[usize; 4]) -> Option<u32> {
    let flip = match kind {
        StereoKind::Tetrahedral => {
            let inversions = (0..4)
                .flat_map(|a| (a + 1..4).map(move |b| (a, b)))
                .filter(|&(a, b)| image[a] > image[b])
                .count();
            (inversions % 2) as u32
        }
        StereoKind::CisTrans => {
            if image[0] / 2 != image[1] / 2 || image[2] / 2 != image[3] / 2 {
                return None;
            }
            ((image[0] % 2) ^ (image[2] % 2)) as u32
        }
        _ => return None,
    };
    Some(coset ^ flip)
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
fn test_stereo_kind_act(#[case] kind: StereoKind) {
    for &permutation in PermutationGroup::symmetric(4).elements() {
        let image = array::from_fn(|position| permutation.apply(position));
        for coset in 0..2 {
            assert_eq!(
                kind.act(coset, permutation),
                binary_coset_action(kind, coset, &image),
                "kind {kind:?}, coset {coset}, image {image:?}"
            );
        }
    }
}
