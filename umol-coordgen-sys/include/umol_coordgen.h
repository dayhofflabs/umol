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

typedef struct {
    double x;
    double y;
} umol_coordgen_point;

enum {
    UMOL_COORDGEN_OK = 0,
    UMOL_COORDGEN_NULL_POINTER = 1,
    UMOL_COORDGEN_ATOM_OUT_OF_BOUNDS = 2,
    UMOL_COORDGEN_ALLOCATION_FAILED = 3,
    UMOL_COORDGEN_BACKEND_EXCEPTION = 4
};

typedef int32_t umol_coordgen_error;

/* Generate one point per input atom while preserving the input atom order.
 * Atomic number zero denotes a generic atom. Empty input is valid. */
umol_coordgen_error umol_coordgen_generate(
    size_t atom_count,
    const uint16_t *atomic_numbers,
    size_t bond_count,
    const umol_coordgen_bond *bonds,
    umol_coordgen_point *points
);

#ifdef __cplusplus
}
#endif

#endif
