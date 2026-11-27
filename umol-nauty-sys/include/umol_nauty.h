#ifndef UMOL_NAUTY_H
#define UMOL_NAUTY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    UMOL_NAUTY_OK = 0,
    UMOL_NAUTY_NULL_POINTER = 1,
    UMOL_NAUTY_INVALID_VERTEX_COUNT = 2,
    UMOL_NAUTY_INVALID_CSR = 3,
    UMOL_NAUTY_INTEGER_OVERFLOW = 4,
    UMOL_NAUTY_ALLOCATION_FAILED = 5,
    UMOL_NAUTY_REENTRANT_CALL = 6,
    UMOL_NAUTY_N_TOO_BIG = 7,
    UMOL_NAUTY_M_TOO_BIG = 8,
    UMOL_NAUTY_CANON_GRAPH_MISSING = 9,
    UMOL_NAUTY_ABORTED = 10,
    UMOL_NAUTY_KILLED = 11,
    UMOL_NAUTY_UNKNOWN_ERROR = 12,
    UMOL_NAUTY_INVALID_PARTITION = 13,
    UMOL_NAUTY_INVALID_WEIGHT = 14
} umol_nauty_error;

/* Called synchronously while umol_nauty_run or umol_nauty_run_weighted is
 * active. permutation[v] is the image of vertex v and is valid only for the
 * duration of the callback. */
typedef void (*umol_nauty_generator_fn)(
    void *context,
    const uint32_t *permutation,
    uint32_t vertex_count
);

/* Canonicalize a nonempty, vertex-colored sparse graph and report generators.
 *
 * offsets has vertex_count + 1 entries, starts at zero, and ends at the
 * directed edge count. neighbors contains that many entries. The caller owns
 * every input and output buffer. canonical_labels and orbits each have
 * vertex_count entries. partition is a permutation of the vertices ordered by
 * nondecreasing color. Colors need not be contiguous.
 *
 * This operation is thread-safe but not reentrant on the same thread. */
umol_nauty_error umol_nauty_run(
    uint32_t vertex_count,
    const size_t *offsets,
    const uint32_t *neighbors,
    const uint32_t *colors,
    const uint32_t *partition,
    uint32_t *canonical_labels,
    uint32_t *orbits,
    umol_nauty_generator_fn report_generator,
    void *callback_context,
    double *group_mantissa,
    int32_t *group_exponent
);

/* Run weighted sparse nauty on a nonempty, vertex- and edge-colored simple
 * undirected sparse graph.
 *
 * weights has directed_edge_count entries aligned with neighbors; each weight
 * must fit in [0, INT32_MAX - 1] and opposite directed edges must carry equal
 * weights. When get_canon is true, canonical_labels must be non-NULL and
 * receives vertex_count entries. When get_canon is false, canonical_labels is
 * ignored and may be NULL, and nauty computes only orbits, group order, and
 * generators. All other contracts match umol_nauty_run. */
umol_nauty_error umol_nauty_run_weighted(
    uint32_t vertex_count,
    const size_t *offsets,
    const uint32_t *neighbors,
    const uint32_t *weights,
    const uint32_t *colors,
    const uint32_t *partition,
    bool get_canon,
    uint32_t *canonical_labels,
    uint32_t *orbits,
    umol_nauty_generator_fn report_generator,
    void *callback_context,
    double *group_mantissa,
    int32_t *group_exponent
);

#ifdef __cplusplus
}
#endif

#endif
