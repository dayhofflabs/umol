#include "umol_nauty.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "nausparse.h"

_Static_assert(sizeof(int) == sizeof(uint32_t), "nauty requires 32-bit int");
_Static_assert(sizeof(int) == sizeof(int32_t), "nauty exponent must fit int32_t");

static TLS_ATTR umol_nauty_generator_fn active_report_generator;
static TLS_ATTR void *active_callback_context;
static TLS_ATTR int active_call;

static void report_automorphism(
    int count,
    int *permutation,
    int *orbits,
    int orbit_count,
    int stabilizer_vertex,
    int vertex_count)
{
    (void)count;
    (void)orbits;
    (void)orbit_count;
    (void)stabilizer_vertex;

    if (active_report_generator != NULL) {
        active_report_generator(
            active_callback_context,
            (const uint32_t *)permutation,
            (uint32_t)vertex_count
        );
    }
}

static umol_nauty_error map_nauty_error(int status)
{
    switch (status) {
    case 0: return UMOL_NAUTY_OK;
    case NTOOBIG: return UMOL_NAUTY_N_TOO_BIG;
    case MTOOBIG: return UMOL_NAUTY_M_TOO_BIG;
    case CANONGNIL: return UMOL_NAUTY_CANON_GRAPH_MISSING;
    case NAUABORTED: return UMOL_NAUTY_ABORTED;
    case NAUKILLED: return UMOL_NAUTY_KILLED;
    default: return UMOL_NAUTY_UNKNOWN_ERROR;
    }
}

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
    int32_t *group_exponent)
{
    umol_nauty_error result = UMOL_NAUTY_OK;
    sparsegraph graph;
    sparsegraph canonical_graph;
    statsblk stats = {0};
    int *lab = NULL;
    int *ptn = NULL;
    int *nauty_orbits = NULL;
    size_t directed_edge_count;
    size_t position;
    int n;
    int m;
    DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(graph);
    SG_INIT(canonical_graph);

    if (active_call) return UMOL_NAUTY_REENTRANT_CALL;
    if (vertex_count == 0 || vertex_count > (uint32_t)INT_MAX) {
        return UMOL_NAUTY_INVALID_VERTEX_COUNT;
    }
    if (offsets == NULL || colors == NULL || partition == NULL || canonical_labels == NULL ||
        orbits == NULL || group_mantissa == NULL || group_exponent == NULL) {
        return UMOL_NAUTY_NULL_POINTER;
    }
    if (offsets[0] != 0) return UMOL_NAUTY_INVALID_CSR;

    directed_edge_count = offsets[vertex_count];
    if (directed_edge_count > 0 && neighbors == NULL) {
        return UMOL_NAUTY_NULL_POINTER;
    }
    if (directed_edge_count > SIZE_MAX / sizeof(int)) {
        return UMOL_NAUTY_INTEGER_OVERFLOW;
    }

    n = (int)vertex_count;
    if ((size_t)n > SIZE_MAX / sizeof(*lab) ||
        (size_t)n > SIZE_MAX / sizeof(*graph.v) ||
        (size_t)n > SIZE_MAX / sizeof(*graph.d)) {
        return UMOL_NAUTY_INTEGER_OVERFLOW;
    }
    lab = malloc((size_t)n * sizeof(*lab));
    ptn = malloc((size_t)n * sizeof(*ptn));
    nauty_orbits = calloc((size_t)n, sizeof(*nauty_orbits));
    graph.v = malloc((size_t)n * sizeof(*graph.v));
    graph.d = malloc((size_t)n * sizeof(*graph.d));
    if (directed_edge_count > 0) {
        graph.e = malloc(directed_edge_count * sizeof(*graph.e));
    }
    if (lab == NULL || ptn == NULL || nauty_orbits == NULL ||
        graph.v == NULL || graph.d == NULL ||
        (directed_edge_count > 0 && graph.e == NULL)) {
        result = UMOL_NAUTY_ALLOCATION_FAILED;
        goto cleanup;
    }

    graph.nv = n;
    graph.nde = directed_edge_count;
    graph.vlen = (size_t)n;
    graph.dlen = (size_t)n;
    graph.elen = directed_edge_count;

    for (position = 0; position < (size_t)n; ++position) {
        size_t begin = offsets[position];
        size_t end = offsets[position + 1];
        size_t edge_position;

        if (end < begin || end > directed_edge_count || end - begin > (size_t)INT_MAX) {
            result = UMOL_NAUTY_INVALID_CSR;
            goto cleanup;
        }
        graph.v[position] = begin;
        graph.d[position] = (int)(end - begin);
        for (edge_position = begin; edge_position < end; ++edge_position) {
            uint32_t neighbor = neighbors[edge_position];
            if (neighbor >= vertex_count) {
                result = UMOL_NAUTY_INVALID_CSR;
                goto cleanup;
            }
            graph.e[edge_position] = (int)neighbor;
        }
    }

    for (position = 0; position < (size_t)n; ++position) {
        uint32_t vertex = partition[position];
        if (vertex >= vertex_count) {
            result = UMOL_NAUTY_INVALID_PARTITION;
            goto cleanup;
        }
        lab[position] = (int)vertex;
        if (nauty_orbits[vertex] != 0) {
            result = UMOL_NAUTY_INVALID_PARTITION;
            goto cleanup;
        }
        nauty_orbits[vertex] = 1;
        if (position > 0 && colors[partition[position - 1]] > colors[vertex]) {
            result = UMOL_NAUTY_INVALID_PARTITION;
            goto cleanup;
        }
    }
    for (position = 0; position < (size_t)n; ++position) {
        ptn[position] = position + 1 < (size_t)n &&
            colors[partition[position]] == colors[partition[position + 1]];
    }

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.userautomproc = report_generator == NULL ? NULL : report_automorphism;
    m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

    active_call = 1;
    active_report_generator = report_generator;
    active_callback_context = callback_context;
    sparsenauty(
        &graph,
        lab,
        ptn,
        nauty_orbits,
        &options,
        &stats,
        &canonical_graph
    );
    active_callback_context = NULL;
    active_report_generator = NULL;
    active_call = 0;

    result = map_nauty_error(stats.errstatus);
    if (result != UMOL_NAUTY_OK) goto cleanup;

    for (position = 0; position < (size_t)n; ++position) {
        canonical_labels[position] = (uint32_t)lab[position];
        orbits[position] = (uint32_t)nauty_orbits[position];
    }
    *group_mantissa = stats.grpsize1;
    *group_exponent = (int32_t)stats.grpsize2;

cleanup:
    active_callback_context = NULL;
    active_report_generator = NULL;
    active_call = 0;
    free(lab);
    free(ptn);
    free(nauty_orbits);
    SG_FREE(graph);
    SG_FREE(canonical_graph);
    return result;
}
