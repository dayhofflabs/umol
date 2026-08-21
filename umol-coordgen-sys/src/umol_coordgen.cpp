#include "umol_coordgen.h"

#include <cstddef>
#include <mutex>
#include <new>
#include <vector>

#include "sketcherMinimizer.h"

namespace {

std::mutex coordgen_mutex;

void delete_unowned_molecule(sketcherMinimizerMolecule* molecule)
{
    for (auto* bond : molecule->_bonds) {
        delete bond;
    }
    for (auto* atom : molecule->_atoms) {
        delete atom;
    }
    delete molecule;
}

} // namespace

extern "C" umol_coordgen_error umol_coordgen_generate(
    size_t atom_count,
    const uint16_t* atomic_numbers,
    size_t bond_count,
    const umol_coordgen_bond* bonds,
    umol_coordgen_point* points)
{
    if ((atom_count != 0 && (atomic_numbers == nullptr || points == nullptr)) ||
        (bond_count != 0 && bonds == nullptr)) {
        return UMOL_COORDGEN_NULL_POINTER;
    }
    for (size_t index = 0; index < bond_count; ++index) {
        if (bonds[index].atom_0 >= atom_count ||
            bonds[index].atom_1 >= atom_count) {
            return UMOL_COORDGEN_ATOM_OUT_OF_BOUNDS;
        }
    }
    if (atom_count == 0) {
        return UMOL_COORDGEN_OK;
    }

    std::lock_guard<std::mutex> lock(coordgen_mutex);
    auto* molecule = new (std::nothrow) sketcherMinimizerMolecule();
    if (molecule == nullptr) {
        return UMOL_COORDGEN_ALLOCATION_FAILED;
    }

    bool ownership_transferred = false;
    try {
        molecule->_atoms.reserve(atom_count);
        molecule->_bonds.reserve(bond_count);

        std::vector<sketcherMinimizerAtom*> input_atoms;
        input_atoms.reserve(atom_count);
        for (size_t index = 0; index < atom_count; ++index) {
            auto* atom = molecule->addNewAtom();
            atom->setAtomicNumber(static_cast<int>(atomic_numbers[index]));
            input_atoms.push_back(atom);
        }
        for (size_t index = 0; index < bond_count; ++index) {
            auto* bond = molecule->addNewBond(input_atoms[bonds[index].atom_0],
                                               input_atoms[bonds[index].atom_1]);
            bond->setBondOrder(static_cast<int>(bonds[index].order));
        }

        sketcherMinimizer minimizer;
        ownership_transferred = true;
        minimizer.initialize(molecule);
        minimizer.runGenerateCoordinates();

        for (size_t index = 0; index < atom_count; ++index) {
            const auto coordinates = input_atoms[index]->getCoordinates();
            points[index].x = static_cast<double>(coordinates.x());
            points[index].y = static_cast<double>(coordinates.y());
        }
        return UMOL_COORDGEN_OK;
    } catch (const std::bad_alloc&) {
        if (!ownership_transferred) {
            delete_unowned_molecule(molecule);
        }
        return UMOL_COORDGEN_ALLOCATION_FAILED;
    } catch (...) {
        if (!ownership_transferred) {
            delete_unowned_molecule(molecule);
        }
        return UMOL_COORDGEN_BACKEND_EXCEPTION;
    }
}
