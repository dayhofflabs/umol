//! Transformations between explicit and implicit hydrogen representations.

use std::iter;

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AddBond, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomHandle, AtomId,
    BondConstraintForm, BondConstraintKey, BondForm, BondHandle, BondId, CisTransStereoForm,
    Constraint, Edit, Edits, ElementForm, Entity, Lattice, Molecule, MoleculeConstraint, NumForm,
    RelationalConstraint, StereoAtomHandle, StereoBondHandle, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, TetrahedralStereoForm,
};
use umol_perm::Permutation;

use crate::ops::transform::Transformer;

/// A hydrogen collapse that cannot preserve its affected count or stereo frame.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum FoldHydrogensError {
    #[error("atom {atom} has negative implicit-hydrogen count {count}")]
    NegativeImplicitHydrogenCount { atom: AtomId, count: i64 },
    #[error("folding hydrogen {hydrogen} would repeat an implicit-hydrogen ligand in {entity}")]
    StereoHydrogenConflict { entity: Entity, hydrogen: AtomId },
}

/// Folds ordinary explicit hydrogens into their partners' implicit-hydrogen counts.
///
/// An eligible hydrogen has concrete default attributes, one ordinary single bond,
/// a known non-hydrogen partner, and no constraints or non-stereo relation membership.
/// Isotope labels, non-default atom or bond state, stereo sites, and explicit
/// molecule-constraint references prevent collapse. Non-literal partner counts
/// leave their explicit hydrogens unchanged.
///
/// Localized valence and degree assertions decrease with the removed bonds, using
/// the existing numeric-form arithmetic. Total and ring assertions are unchanged.
/// Actual stereo ligands become virtual hydrogens in their stored positions;
/// configurations and stereo constraints are preserved. Tetrahedral and cis/trans
/// assertions are transported to the changed implicit neighbor frame.
/// Molecule-level predicates receive reference compaction only.
///
/// # Errors
///
/// Reports a negative literal count on an otherwise eligible partner, or a fold
/// that would require repeated implicit-hydrogen ligands in a stereo frame.
/// All rejections precede mutation. Consuming failure drops the input; in-place
/// failure leaves it unchanged. Iteration lazily yields one independent result
/// on success and no result on rejection.
///
/// # Semantic properties
///
/// Every successful result preserves molecule representation integrity and is
/// unchanged by another application under exact representation equality.
/// With literal elements and hydrogen counts, the number of hydrogen atoms plus
/// the sum of implicit-hydrogen counts is conserved. Retained atoms preserve their
/// total valence, degree, and hydrogen readings under the surviving atom correspondence.
///
/// Atom and bond identities are removed, so no general inverse or lossless
/// roundtrip with hydrogen unfolding is promised.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HydrogenFolder;

impl HydrogenFolder {
    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, FoldHydrogensError> {
        let mut protected_atoms = Vec::new();
        let mut protected_bonds = Vec::new();
        for constraint in molecule.constraints().iter() {
            collect_topology_references(constraint, &mut protected_atoms, &mut protected_bonds);
        }
        protected_atoms.sort_unstable();
        protected_atoms.dedup();
        protected_bonds.sort_unstable();
        protected_bonds.dedup();

        let ordinary_hydrogen = AtomForm::from_element(Element::H).into_concrete();
        let ordinary_bond = BondForm::from_order(1).into_concrete();
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut counts = Vec::new();
        let mut valences = Vec::new();
        let mut degrees = Vec::new();
        let mut tetrahedral = Vec::new();
        for atom in molecule.atoms().iter() {
            match atom.element() {
                ElementForm::Undetermined | ElementForm::Lit(Element::H) => continue,
                ElementForm::Lit(_) => {}
                element if !ElementForm::not(Element::H).matches(element) => continue,
                _ => {}
            }
            let NumForm::Lit(count) = atom.implicit_hydrogens() else {
                continue;
            };
            let start = atoms.len();
            for neighbor in atom.neighbors() {
                let hydrogen = neighbor.atom_id();
                let bond = neighbor.bond_id();
                if neighbor.atom().attributes() != &ordinary_hydrogen
                    || neighbor.bond().attributes() != &ordinary_bond
                    || neighbor.atom().neighbors().len() != 1
                    || protected_atoms.binary_search(&hydrogen).is_ok()
                    || protected_bonds.binary_search(&bond).is_ok()
                    || molecule.dative_bonds().has_incident(hydrogen)
                    || molecule.aromatic_systems().has_incident(hydrogen)
                    || molecule.multicenter_bonds().has_incident(hydrogen)
                    || molecule.noncovalent_bonds().has_incident(hydrogen)
                    || molecule.stereo_atoms().is_at(hydrogen)
                    || molecule.stereo_bonds().is_at(bond)
                {
                    continue;
                }
                atoms.push(hydrogen);
                bonds.push(BondHandle::Id(bond));
            }
            let added = (atoms.len() - start) as i64;
            if added == 0 {
                continue;
            }
            if *count < 0 {
                return Err(FoldHydrogensError::NegativeImplicitHydrogenCount {
                    atom: atom.id(),
                    count: *count,
                });
            }
            if count + added > 1 {
                let entity = molecule
                    .stereo_atoms()
                    .at_id(atom.id())
                    .map(Entity::StereoAtom)
                    .or_else(|| {
                        molecule
                            .stereo_bonds()
                            .incident_as_site_ids(atom.id())
                            .next()
                            .map(Entity::StereoBond)
                    });
                if let Some(entity) = entity {
                    return Err(FoldHydrogensError::StereoHydrogenConflict {
                        entity,
                        hydrogen: atoms[start],
                    });
                }
            }
            counts.push((
                AtomHandle::Id(atom.id()),
                NumForm::Lit(*count),
                NumForm::Lit(count + added),
            ));
            let constraints = &atom.attributes().constraints;
            if let Some(old @ AtomConstraintForm::Valence(value)) =
                constraints.get(AtomConstraintKey::Valence)
            {
                valences.push((
                    AtomHandle::Id(atom.id()),
                    Some(old.clone()),
                    Some(AtomConstraintForm::valence(value - added)),
                ));
            }
            if let Some(old @ AtomConstraintForm::Degree(value)) =
                constraints.get(AtomConstraintKey::Degree)
            {
                degrees.push((
                    AtomHandle::Id(atom.id()),
                    Some(old.clone()),
                    Some(AtomConstraintForm::degree(value - added)),
                ));
            }
            if let Some(
                old @ AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(coset)),
            ) = constraints.get(AtomConstraintKey::TetrahedralStereo)
            {
                let mut image = [0, 1, 2, 3];
                fold_neighbor_frame(molecule, atom.id(), None, &atoms[start..], &mut image);
                if image != [0, 1, 2, 3] {
                    let new = AtomConstraintForm::tetrahedral_stereo(
                        TetrahedralStereoForm::Stereo(transport_coset(
                            coset,
                            StereoKind::Tetrahedral,
                            Permutation::from_image(&image),
                        )),
                    );
                    if &new != old {
                        tetrahedral.push((AtomHandle::Id(atom.id()), Some(old.clone()), Some(new)));
                    }
                }
            }
        }
        let mut edits = Edits::new();
        if atoms.is_empty() {
            return Ok(edits);
        }
        atoms.sort_unstable();
        edits.modify_atom_implicit_hydrogens(counts);
        for (key, changes) in [
            (AtomConstraintKey::Valence, valences),
            (AtomConstraintKey::Degree, degrees),
            (AtomConstraintKey::TetrahedralStereo, tetrahedral),
        ] {
            if !changes.is_empty() {
                edits.modify_atom_constraints(key, changes);
            }
        }
        let mut cis_trans = Vec::new();
        for bond in molecule.bonds().iter() {
            let Some(old @ BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) =
                bond.attributes()
                    .constraints
                    .get(BondConstraintKey::CisTransStereo)
            else {
                continue;
            };
            let [a, b] = bond.atom_ids();
            let mut image = [0, 1, 2, 3];
            fold_neighbor_frame(molecule, a, Some(b), &atoms, &mut image[..2]);
            fold_neighbor_frame(molecule, b, Some(a), &atoms, &mut image[2..]);
            if image != [0, 1, 2, 3] {
                let new = BondConstraintForm::cis_trans_stereo(CisTransStereoForm::Stereo(
                    transport_coset(coset, StereoKind::CisTrans, Permutation::from_image(&image)),
                ));
                if &new != old {
                    cis_trans.push((BondHandle::Id(bond.id()), Some(old.clone()), Some(new)));
                }
            }
        }
        if !cis_trans.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::CisTransStereo, cis_trans);
        }
        for stereo in molecule.stereo_atoms().iter() {
            if let Some(new) = fold_ligands(
                molecule,
                Entity::StereoAtom(stereo.id()),
                stereo.ligand_ids(),
                &atoms,
            )? {
                edits.push(Edit::ReplaceStereoAtomLigands {
                    id: StereoAtomHandle::Id(stereo.id()),
                    old: stereo
                        .ligand_ids()
                        .iter()
                        .map(|l| (AtomHandle::Id(l.atom_id), l.kind))
                        .collect(),
                    new,
                });
            }
        }
        for stereo in molecule.stereo_bonds().iter() {
            if let Some(new) = fold_ligands(
                molecule,
                Entity::StereoBond(stereo.id()),
                stereo.ligand_ids(),
                &atoms,
            )? {
                edits.push(Edit::ReplaceStereoBondLigands {
                    id: StereoBondHandle::Id(stereo.id()),
                    old: stereo
                        .ligand_ids()
                        .iter()
                        .map(|l| (AtomHandle::Id(l.atom_id), l.kind))
                        .collect(),
                    new,
                });
            }
        }
        edits.remove_topology(atoms.into_iter().map(AtomHandle::Id).collect(), bonds);
        Ok(edits)
    }
}

impl Transformer for HydrogenFolder {
    type Error = FoldHydrogensError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error> {
        let edits = self.plan_transform(&molecule)?;
        if edits.is_empty() {
            return Ok(molecule);
        }
        Ok(molecule
            .apply(edits)
            .expect("hydrogen folding plan applies and preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

fn fold_ligands(
    molecule: &Molecule,
    entity: Entity,
    ligands: &[StereoLigand],
    removed: &[AtomId],
) -> Result<Option<Vec<(AtomHandle, StereoLigandKind)>>, FoldHydrogensError> {
    if !ligands
        .iter()
        .any(|l| l.kind == StereoLigandKind::Atom && removed.binary_search(&l.atom_id).is_ok())
    {
        return Ok(None);
    }
    let mut new = Vec::with_capacity(ligands.len());
    for ligand in ligands {
        if ligand.kind != StereoLigandKind::Atom || removed.binary_search(&ligand.atom_id).is_err()
        {
            new.push((AtomHandle::Id(ligand.atom_id), ligand.kind));
            continue;
        }
        let bearer = molecule
            .atom(ligand.atom_id)
            .neighbors()
            .next()
            .expect("eligible hydrogen has one neighbor")
            .atom_id();
        let replacement = (AtomHandle::Id(bearer), StereoLigandKind::ImplicitHydrogen);
        if new.contains(&replacement)
            || ligands.contains(&StereoLigand::new(
                bearer,
                StereoLigandKind::ImplicitHydrogen,
            ))
        {
            return Err(FoldHydrogensError::StereoHydrogenConflict {
                entity,
                hydrogen: ligand.atom_id,
            });
        }
        new.push(replacement);
    }
    Ok(Some(new))
}

fn fold_neighbor_frame(
    molecule: &Molecule,
    atom: AtomId,
    opposite: Option<AtomId>,
    removed: &[AtomId],
    image: &mut [usize],
) {
    let mut folded = [false; 4];
    let mut count = 0;
    for neighbor in molecule
        .atom(atom)
        .neighbors()
        .filter(|n| Some(n.atom_id()) != opposite)
    {
        if count == image.len() {
            return;
        }
        folded[count] = removed.binary_search(&neighbor.atom_id()).is_ok();
        count += 1;
    }
    let offset = image[0];
    image[..count].sort_unstable_by_key(|&i| (folded[i - offset], i));
}

fn transport_coset(coset: &StereoCoset, kind: StereoKind, permutation: Permutation) -> StereoCoset {
    let configuration = StereoConfigurationForm::kinded(kind, coset.clone());
    match configuration.apply(permutation) {
        Some(StereoConfigurationForm::Kinded(_, transported)) => transported,
        // An already contradictory coset stays contradictory in the new frame.
        _ => coset.clone(),
    }
}

fn collect_topology_references(
    constraint: &Constraint,
    atoms: &mut Vec<AtomId>,
    bonds: &mut Vec<BondId>,
) {
    match constraint {
        Constraint::Atom(atom, _) => atoms.push(*atom),
        Constraint::Bond(bond, _) => bonds.push(*bond),
        Constraint::And(items) | Constraint::Or(items) => {
            for item in items {
                collect_topology_references(item, atoms, bonds);
            }
        }
        Constraint::Not(item) => collect_topology_references(item, atoms, bonds),
        Constraint::Molecule(item) => match item {
            MoleculeConstraint::ChargeSum {
                atoms: Some(ids), ..
            }
            | MoleculeConstraint::UnpairedElectronCoupling {
                atoms: Some(ids), ..
            }
            | MoleculeConstraint::Connected { atoms: Some(ids) } => atoms.extend(ids),
            MoleculeConstraint::BondOrderSum {
                bonds: Some(ids), ..
            } => bonds.extend(ids),
            MoleculeConstraint::ChargeSum { atoms: None, .. }
            | MoleculeConstraint::UnpairedElectronCoupling { atoms: None, .. }
            | MoleculeConstraint::Connected { atoms: None }
            | MoleculeConstraint::BondOrderSum { bonds: None, .. } => {}
        },
        Constraint::Relational(item) => match item {
            RelationalConstraint::DativeBondDonors { atoms: ids, .. }
            | RelationalConstraint::DativeBondContainsAllDonors { atoms: ids, .. }
            | RelationalConstraint::AromaticSystemAtoms { atoms: ids, .. }
            | RelationalConstraint::AromaticSystemContainsAll { atoms: ids, .. }
            | RelationalConstraint::MulticenterBondAtoms { atoms: ids, .. }
            | RelationalConstraint::MulticenterBondContainsAll { atoms: ids, .. }
            | RelationalConstraint::StereoAtomLigands { atoms: ids, .. }
            | RelationalConstraint::StereoBondLigands { atoms: ids, .. } => atoms.extend(ids),
            RelationalConstraint::NoncovalentBondEnds { atoms: ids, .. } => atoms.extend(ids),
            RelationalConstraint::DativeBondDonor { atom, .. }
            | RelationalConstraint::DativeBondAcceptor { atom, .. }
            | RelationalConstraint::AromaticSystemContains { atom, .. }
            | RelationalConstraint::MulticenterBondContains { atom, .. }
            | RelationalConstraint::NoncovalentBondContains { atom, .. }
            | RelationalConstraint::StereoAtomSite { atom, .. }
            | RelationalConstraint::StereoAtomContains { atom, .. }
            | RelationalConstraint::StereoBondContains { atom, .. } => atoms.push(*atom),
            RelationalConstraint::DativeBondParallels { parallel: bond, .. }
            | RelationalConstraint::StereoBondSite { bond, .. } => bonds.push(*bond),
            RelationalConstraint::DativeBondAllDonors { .. }
            | RelationalConstraint::DativeBondAnyDonor { .. }
            | RelationalConstraint::DativeBondAcceptorSatisfies { .. }
            | RelationalConstraint::AromaticSystemAllAtoms { .. }
            | RelationalConstraint::AromaticSystemAnyAtom { .. }
            | RelationalConstraint::MulticenterBondAllAtoms { .. }
            | RelationalConstraint::MulticenterBondAnyAtom { .. }
            | RelationalConstraint::NoncovalentBondEndsSatisfy { .. }
            | RelationalConstraint::StereoAtomAllLigands { .. }
            | RelationalConstraint::StereoAtomAnyLigand { .. }
            | RelationalConstraint::StereoBondAllLigands { .. }
            | RelationalConstraint::StereoBondAnyLigand { .. } => {}
        },
        Constraint::DativeBond(..)
        | Constraint::AromaticSystem(..)
        | Constraint::MulticenterBond(..)
        | Constraint::NoncovalentBond(..)
        | Constraint::StereoAtom(..)
        | Constraint::StereoBond(..) => {}
    }
}

/// A hydrogen expansion that cannot preserve its affected count or stereo frame.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum UnfoldHydrogensError {
    #[error("atom {atom} has negative implicit-hydrogen count {count}")]
    NegativeImplicitHydrogenCount { atom: AtomId, count: i64 },
    #[error("implicit-hydrogen ligand on atom {atom} in {entity} has a non-concrete count")]
    NonConcreteStereoHydrogenCount {
        entity: Entity,
        atom: AtomId,
        value: NumForm,
    },
    #[error(
        "atom {atom} has implicit-hydrogen count {implicit_hydrogens} but {virtual_hydrogens} virtual hydrogens in {entity}"
    )]
    StereoHydrogenCountMismatch {
        entity: Entity,
        atom: AtomId,
        implicit_hydrogens: i64,
        virtual_hydrogens: usize,
    },
}

/// Expands literal implicit-hydrogen counts into ordinary hydrogen atoms and single bonds.
///
/// Each positive count produces that many concrete default hydrogens and bonds,
/// appended in bearer-atom order, and becomes literal zero. Non-literal counts
/// remain unchanged unless they back a virtual hydrogen in a stereo frame.
/// All existing entities retain their identities; unrelated fields are unchanged.
///
/// Localized valence and degree assertions increase with the added bonds, using
/// the existing numeric-form arithmetic. Total and ring assertions are unchanged.
/// Virtual hydrogen ligands become actual ligands in the same stored positions;
/// configurations and frame-relative constraints are preserved. Appending the
/// hydrogens preserves the implicit neighbor ordering of tetrahedral and cis/trans
/// assertions. Molecule-level predicates are unchanged.
///
/// # Errors
///
/// Reports any negative literal implicit-hydrogen count, a non-literal count
/// backing a virtual hydrogen, or a disagreement between a stereo bearer's
/// literal count and its frame's virtual hydrogen census. Each stereo-atom frame
/// and stereo-bond endpoint supplies the complete census for its bearer.
/// All rejections precede mutation. Consuming failure drops the input; in-place
/// failure leaves it unchanged. Iteration lazily yields one independent result
/// on success and no result on rejection.
///
/// # Semantic properties
///
/// Every successful result preserves molecule representation integrity and is
/// unchanged by another application under exact representation equality.
/// With literal elements and hydrogen counts, the number of hydrogen atoms plus
/// the sum of implicit-hydrogen counts is conserved. Existing atoms preserve their
/// total valence, degree, and hydrogen readings at their original ids.
///
/// Unfolding followed by folding restores the exact stored representation when
/// no explicit hydrogen is initially eligible for folding, all hydrogen counts
/// are nonnegative literals, positive counts belong to known non-hydrogen atoms,
/// affected localized valence/degree assertions are literal or absent, and stereo
/// frames and assertions are compatible with the hydrogen counts.
/// No general inverse or lossless roundtrip with hydrogen folding is promised.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HydrogenUnfolder;

impl HydrogenUnfolder {
    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, UnfoldHydrogensError> {
        let ordinary_hydrogen = AtomForm::from_element(Element::H).into_concrete();
        let ordinary_bond = BondForm::from_order(1).into_concrete();
        let has_stereo = molecule.has_stereo_atoms() || molecule.has_stereo_bonds();
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut hydrogens = Vec::new();
        let mut counts = Vec::new();
        let mut valences = Vec::new();
        let mut degrees = Vec::new();
        for atom in molecule.atoms().iter() {
            let NumForm::Lit(count) = atom.implicit_hydrogens() else {
                continue;
            };
            if *count < 0 {
                return Err(UnfoldHydrogensError::NegativeImplicitHydrogenCount {
                    atom: atom.id(),
                    count: *count,
                });
            }
            if *count == 0 {
                continue;
            }
            if has_stereo && *count == 1 {
                hydrogens.push((atom.id(), atoms.len()));
            }
            for _ in 0..*count {
                bonds.push(AddBond {
                    endpoints: [AtomHandle::Id(atom.id()), AtomHandle::New(atoms.len())],
                    attributes: ordinary_bond.clone(),
                });
                atoms.push(ordinary_hydrogen.clone());
            }
            counts.push((
                AtomHandle::Id(atom.id()),
                NumForm::Lit(*count),
                NumForm::Lit(0),
            ));
            let constraints = &atom.attributes().constraints;
            if let Some(old @ AtomConstraintForm::Valence(value)) =
                constraints.get(AtomConstraintKey::Valence)
            {
                valences.push((
                    AtomHandle::Id(atom.id()),
                    Some(old.clone()),
                    Some(AtomConstraintForm::valence(value + *count)),
                ));
            }
            if let Some(old @ AtomConstraintForm::Degree(value)) =
                constraints.get(AtomConstraintKey::Degree)
            {
                degrees.push((
                    AtomHandle::Id(atom.id()),
                    Some(old.clone()),
                    Some(AtomConstraintForm::degree(value + *count)),
                ));
            }
        }
        let mut edits = Edits::new();
        if !atoms.is_empty() {
            edits.push(Edit::AddAtoms { atoms });
            edits.push(Edit::AddBonds { bonds });
            edits.modify_atom_implicit_hydrogens(counts);
            for (key, changes) in [
                (AtomConstraintKey::Valence, valences),
                (AtomConstraintKey::Degree, degrees),
            ] {
                if !changes.is_empty() {
                    edits.modify_atom_constraints(key, changes);
                }
            }
        }
        for stereo in molecule.stereo_atoms().iter() {
            if let Some(new) = unfold_ligands(
                molecule,
                Entity::StereoAtom(stereo.id()),
                &[stereo.site_id()],
                stereo.ligand_ids(),
                &hydrogens,
            )? {
                edits.push(Edit::ReplaceStereoAtomLigands {
                    id: StereoAtomHandle::Id(stereo.id()),
                    old: stereo
                        .ligand_ids()
                        .iter()
                        .map(|l| (AtomHandle::Id(l.atom_id), l.kind))
                        .collect(),
                    new,
                });
            }
        }
        for stereo in molecule.stereo_bonds().iter() {
            if let Some(new) = unfold_ligands(
                molecule,
                Entity::StereoBond(stereo.id()),
                &stereo.site().atom_ids(),
                stereo.ligand_ids(),
                &hydrogens,
            )? {
                edits.push(Edit::ReplaceStereoBondLigands {
                    id: StereoBondHandle::Id(stereo.id()),
                    old: stereo
                        .ligand_ids()
                        .iter()
                        .map(|l| (AtomHandle::Id(l.atom_id), l.kind))
                        .collect(),
                    new,
                });
            }
        }
        Ok(edits)
    }
}

impl Transformer for HydrogenUnfolder {
    type Error = UnfoldHydrogensError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error> {
        let edits = self.plan_transform(&molecule)?;
        if edits.is_empty() {
            return Ok(molecule);
        }
        Ok(molecule
            .apply(edits)
            .expect("hydrogen unfolding plan applies and preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

fn unfold_ligands(
    molecule: &Molecule,
    entity: Entity,
    sites: &[AtomId],
    ligands: &[StereoLigand],
    hydrogens: &[(AtomId, usize)],
) -> Result<Option<Vec<(AtomHandle, StereoLigandKind)>>, UnfoldHydrogensError> {
    let mut changed = false;
    for &atom in sites {
        let virtual_hydrogens = ligands
            .iter()
            .filter(|l| l.atom_id == atom && l.kind == StereoLigandKind::ImplicitHydrogen)
            .count();
        let value = molecule.atom(atom).implicit_hydrogens();
        match value {
            NumForm::Lit(count) if *count != virtual_hydrogens as i64 => {
                return Err(UnfoldHydrogensError::StereoHydrogenCountMismatch {
                    entity,
                    atom,
                    implicit_hydrogens: *count,
                    virtual_hydrogens,
                });
            }
            NumForm::Lit(_) => {}
            _ if virtual_hydrogens != 0 => {
                return Err(UnfoldHydrogensError::NonConcreteStereoHydrogenCount {
                    entity,
                    atom,
                    value: value.clone(),
                });
            }
            _ => {}
        }
        changed |= virtual_hydrogens != 0;
    }
    if !changed {
        return Ok(None);
    }
    Ok(Some(
        ligands
            .iter()
            .map(|ligand| {
                if ligand.kind == StereoLigandKind::ImplicitHydrogen {
                    let index = hydrogens
                        .binary_search_by_key(&ligand.atom_id, |&(atom, _)| atom)
                        .expect("the checked stereo census has one planned hydrogen per bearer");
                    (AtomHandle::New(hydrogens[index].1), StereoLigandKind::Atom)
                } else {
                    (AtomHandle::Id(ligand.atom_id), ligand.kind)
                }
            })
            .collect(),
    ))
}

#[cfg(test)]
mod tests;
