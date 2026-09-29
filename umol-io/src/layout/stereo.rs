//! Selection of definite cis/trans sites and their displayed ligands.

use umol_graph_ir::ir::{
    AtomId, BondId, CisTransConfiguration, Molecule, StereoBondView, StereoCoset, StereoKind,
    StereoLigand, StereoLigandKind,
};

/// A stereo bond with a literal cis/trans coset, read over one actual ligand per site atom.
///
/// `ligands[i]` is bonded to `site[i]`. `configuration` is the stored coset in the frame
/// `[ligands[0], other_0, ligands[1], other_1]`: `Z` places the two selected ligands on the same
/// side of the site axis, `E` on opposite sides.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CisTransSite {
    pub(crate) bond: BondId,
    pub(crate) site: [AtomId; 2],
    pub(crate) ligands: [AtomId; 2],
    pub(crate) configuration: CisTransConfiguration,
}

pub(crate) fn cis_trans_site(
    molecule: &Molecule,
    stereo: StereoBondView<'_>,
) -> Option<CisTransSite> {
    if stereo.attributes.configuration.kind() != Some(StereoKind::CisTrans) {
        return None;
    }
    let ligands = stereo
        .ligands()
        .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
        .collect::<Vec<_>>();
    let [first_0, first_1, second_0, second_1] = ligands.as_slice() else {
        return None;
    };
    let [site_0, site_1] = stereo.site().atom_ids();
    let (first_pair, second_pair) = if [first_0, first_1]
        .into_iter()
        .all(|ligand| ligand_matches_endpoint(molecule, *ligand, site_0, site_1))
    {
        ([*first_0, *first_1], [*second_0, *second_1])
    } else {
        ([*second_0, *second_1], [*first_0, *first_1])
    };
    let (first_ligand, first_other) = select_actual_ligand(first_pair)?;
    let (second_ligand, second_other) = select_actual_ligand(second_pair)?;
    let requested = [first_ligand, first_other, second_ligand, second_other];
    let StereoCoset::Lit(coset) = stereo.coset_for(requested)? else {
        return None;
    };
    let configuration = match coset {
        0 => CisTransConfiguration::Z,
        1 => CisTransConfiguration::E,
        _ => return None,
    };

    Some(CisTransSite {
        bond: stereo.site_id(),
        site: [site_0, site_1],
        ligands: [first_ligand.atom_id, second_ligand.atom_id],
        configuration,
    })
}

fn ligand_matches_endpoint(
    molecule: &Molecule,
    ligand: StereoLigand,
    endpoint: AtomId,
    other_endpoint: AtomId,
) -> bool {
    match ligand.kind {
        StereoLigandKind::Atom => {
            ligand.atom_id != other_endpoint
                && molecule
                    .neighbors(endpoint)
                    .any(|neighbor| neighbor.atom_id() == ligand.atom_id)
        }
        StereoLigandKind::ImplicitHydrogen | StereoLigandKind::LonePair => {
            ligand.atom_id == endpoint
        }
    }
}

fn select_actual_ligand(pair: [StereoLigand; 2]) -> Option<(StereoLigand, StereoLigand)> {
    match pair {
        [first, second] if first.kind == StereoLigandKind::Atom => Some((first, second)),
        [first, second] if second.kind == StereoLigandKind::Atom => Some((second, first)),
        _ => None,
    }
}
