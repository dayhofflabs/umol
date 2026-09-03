#ifndef UMOL_COORDGEN_H
#define UMOL_COORDGEN_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t atom_0;
    size_t atom_1;
    uint8_t order;
} umol_coordgen_bond;

enum {
    UMOL_COORDGEN_SAME_SIDE = 0,
    UMOL_COORDGEN_OPPOSITE_SIDE = 1
};

typedef struct {
    size_t bond;
    size_t first_ligand;
    size_t second_ligand;
    uint8_t relation;
} umol_coordgen_cis_trans_bond;

typedef struct {
    double x;
    double y;
} umol_coordgen_point;

enum {
    UMOL_COORDGEN_OK = 0,
    UMOL_COORDGEN_NULL_POINTER = 1,
    UMOL_COORDGEN_ATOM_OUT_OF_BOUNDS = 2,
    UMOL_COORDGEN_ALLOCATION_FAILED = 3,
    UMOL_COORDGEN_BACKEND_EXCEPTION = 4,
    UMOL_COORDGEN_CIS_TRANS_SITE_OUT_OF_BOUNDS = 5,
    UMOL_COORDGEN_CIS_TRANS_LIGAND_OUT_OF_BOUNDS = 6,
    UMOL_COORDGEN_INVALID_SIDE_RELATION = 7
};

typedef int32_t umol_coordgen_error;

/* Generate one point per input atom while preserving the input atom order.
 * Atomic number zero denotes a generic atom. Empty input is valid. */
umol_coordgen_error umol_coordgen_generate(
    size_t atom_count,
    const uint16_t *atomic_numbers,
    size_t bond_count,
    const umol_coordgen_bond *bonds,
    size_t cis_trans_bond_count,
    const umol_coordgen_cis_trans_bond *cis_trans_bonds,
    umol_coordgen_point *points
);

#ifdef __cplusplus
}
#endif

#endif
