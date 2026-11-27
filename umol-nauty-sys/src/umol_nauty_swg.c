/* Weighted sparse-graph dispatch callbacks for nauty over sparsegraph.w.
 *
 * Adapted from nauty 2.9.3 `nausparse.c` (Copyright (1984-) Brendan McKay,
 * Apache License, Version 2.0; see ../nauty/COPYRIGHT and ../LICENSE-APACHE),
 * extending the sparse dispatch callbacks to operate over non-negative edge
 * weights in `sparsegraph.w` without modifying the vendored upstream sources. */

#include "umol_nauty_swg.h"

#include <limits.h>
#include <stdlib.h>

#define MASH(l, i) ((long)((((unsigned long)(l) ^ 065435UL) + (unsigned long)(i)) & 077777UL))
#define CLEANUP(l) ((int)((l) % 077777L))

typedef struct {
    int weight;
    int target;
} swg_hit_edge;

DYNALLSTAT(short, swg_mark1, swg_mark1_sz);
DYNALLSTAT(short, swg_mark2, swg_mark2_sz);
DYNALLSTAT(int, swg_wmark1, swg_wmark1_sz);
DYNALLSTAT(int, swg_work1, swg_work1_sz);
DYNALLSTAT(int, swg_work2, swg_work2_sz);
DYNALLSTAT(int, swg_work3, swg_work3_sz);
DYNALLSTAT(int, swg_work4, swg_work4_sz);
DYNALLSTAT(int, swg_work5, swg_work5_sz);
DYNALLSTAT(swg_hit_edge, swg_edges, swg_edges_sz);
DYNALLSTAT(set, swg_snwork, swg_snwork_sz);

static TLS_ATTR short swg_mark1_val = 32000;
static TLS_ATTR short swg_mark2_val = 32000;

#define MARK1(i) (swg_mark1[i] = swg_mark1_val)
#define UNMARK1(i) (swg_mark1[i] = 0)
#define ISMARKED1(i) (swg_mark1[i] == swg_mark1_val)
#define ISNOTMARKED1(i) (swg_mark1[i] != swg_mark1_val)

#define MARK2(i) (swg_mark2[i] = swg_mark2_val)
#define ISMARKED2(i) (swg_mark2[i] == swg_mark2_val)

#define RESETMARKS1 do { \
    if (swg_mark1_val++ >= 32000) { \
        size_t ij_; \
        for (ij_ = 0; ij_ < swg_mark1_sz; ++ij_) swg_mark1[ij_] = 0; \
        swg_mark1_val = 1; \
    } \
} while (0)

#define RESETMARKS2 do { \
    if (swg_mark2_val++ >= 32000) { \
        size_t ij_; \
        for (ij_ = 0; ij_ < swg_mark2_sz; ++ij_) swg_mark2[ij_] = 0; \
        swg_mark2_val = 1; \
    } \
} while (0)

static void prepare_marks1(size_t nn)
{
    size_t oldsize = swg_mark1_sz;
    short *oldpos = swg_mark1;
    DYNALLOC1(short, swg_mark1, swg_mark1_sz, nn, "swg prepare_marks1");
    if (swg_mark1_sz != oldsize || swg_mark1 != oldpos) swg_mark1_val = 32000;
}

static void prepare_marks2(size_t nn)
{
    size_t oldsize = swg_mark2_sz;
    short *oldpos = swg_mark2;
    DYNALLOC1(short, swg_mark2, swg_mark2_sz, nn, "swg prepare_marks2");
    if (swg_mark2_sz != oldsize || swg_mark2 != oldpos) swg_mark2_val = 32000;
}

static void sort_ints_asc(int *x, int n)
{
    int i, j, h, xi;
    if (n <= 1) return;
    if (n <= 16) {
        for (i = 1; i < n; ++i) {
            xi = x[i];
            j = i;
            while (j > 0 && x[j - 1] > xi) {
                x[j] = x[j - 1];
                --j;
            }
            x[j] = xi;
        }
        return;
    }
    h = 1;
    while (h < n / 3) h = 3 * h + 1;
    while (h > 0) {
        for (i = h; i < n; ++i) {
            xi = x[i];
            j = i;
            while (j >= h && x[j - h] > xi) {
                x[j] = x[j - h];
                j -= h;
            }
            x[j] = xi;
        }
        h /= 3;
    }
}

static void sort_indirect_asc(int *lab, const int *key, int n)
{
    int i, j, h, li, ki;
    if (n <= 1) return;
    if (n <= 16) {
        for (i = 1; i < n; ++i) {
            li = lab[i];
            ki = key[li];
            j = i;
            while (j > 0 && key[lab[j - 1]] > ki) {
                lab[j] = lab[j - 1];
                --j;
            }
            lab[j] = li;
        }
        return;
    }
    h = 1;
    while (h < n / 3) h = 3 * h + 1;
    while (h > 0) {
        for (i = h; i < n; ++i) {
            li = lab[i];
            ki = key[li];
            j = i;
            while (j >= h && key[lab[j - h]] > ki) {
                lab[j] = lab[j - h];
                j -= h;
            }
            lab[j] = li;
        }
        h /= 3;
    }
}

static void sort_hit_edges_by_weight(swg_hit_edge *edges, size_t n)
{
    size_t i, j, h;
    swg_hit_edge cur;
    if (n <= 1) return;
    if (n <= 24) {
        for (i = 1; i < n; ++i) {
            cur = edges[i];
            j = i;
            while (j > 0 && edges[j - 1].weight > cur.weight) {
                edges[j] = edges[j - 1];
                --j;
            }
            edges[j] = cur;
        }
        return;
    }
    h = 1;
    while (h < n / 3) h = 3 * h + 1;
    while (h > 0) {
        for (i = h; i < n; ++i) {
            cur = edges[i];
            j = i;
            while (j >= h && edges[j - h].weight > cur.weight) {
                edges[j] = edges[j - h];
                j -= h;
            }
            edges[j] = cur;
        }
        h /= 3;
    }
}

static boolean isautom_swg(graph *g, int *p, boolean digraph, int m, int n)
{
    size_t *v, vi, vpi, j;
    int *d, *e;
    sg_weight *wt;
    int i, pi, di, u;

    (void)m;
    SWG_VDE(g, v, d, e, wt);
    prepare_marks1((size_t)n);
    DYNALLOC1(int, swg_wmark1, swg_wmark1_sz, (size_t)n, "isautom_swg");

    for (i = 0; i < n; ++i) {
        if (p[i] != i || digraph) {
            pi = p[i];
            di = d[i];
            if (d[pi] != di) return FALSE;

            vi = v[i];
            vpi = v[pi];
            RESETMARKS1;
            for (j = 0; j < (size_t)di; ++j) {
                u = p[e[vi + j]];
                MARK1(u);
                swg_wmark1[u] = wt[vi + j];
            }
            for (j = 0; j < (size_t)di; ++j) {
                u = e[vpi + j];
                if (ISNOTMARKED1(u) || swg_wmark1[u] != wt[vpi + j]) {
                    return FALSE;
                }
            }
        }
    }
    return TRUE;
}

/* Compare g^lab against canong by row degree, then lexicographically over
 * canonical columns 0..n-1 with edge absence sorting before any present weight. */
static int testcanlab_swg(graph *g, graph *canong, int *lab, int *samerows, int m, int n)
{
    size_t *v, *cv, vi, vli, j;
    int *d, *e, *cd, *ce;
    sg_weight *wt, *cwt;
    int i, di, dli, k, min_g_col, min_g_wt, min_c_col, min_c_wt;

    (void)m;
    SWG_VDE(g, v, d, e, wt);
    SWG_VDE(canong, cv, cd, ce, cwt);

    DYNALLOC1(int, swg_work1, swg_work1_sz, (size_t)n, "testcanlab_swg");
    DYNALLOC1(int, swg_wmark1, swg_wmark1_sz, (size_t)n, "testcanlab_swg");
    prepare_marks1((size_t)n);

    for (i = 0; i < n; ++i) swg_work1[lab[i]] = i;

    for (i = 0; i < n; ++i) {
        vi = cv[i];
        di = cd[i];
        vli = v[lab[i]];
        dli = d[lab[i]];

        if (dli != di) {
            *samerows = i;
            return dli < di ? -1 : 1;
        }

        RESETMARKS1;
        for (j = 0; j < (size_t)di; ++j) {
            k = ce[vi + j];
            MARK1(k);
            swg_wmark1[k] = cwt[vi + j];
        }

        min_g_col = n;
        min_g_wt = 0;
        for (j = 0; j < (size_t)di; ++j) {
            int gw = wt[vli + j];
            k = swg_work1[e[vli + j]];
            if (ISMARKED1(k) && swg_wmark1[k] == gw) {
                UNMARK1(k);
            } else if (k < min_g_col) {
                min_g_col = k;
                min_g_wt = gw;
            }
        }

        if (min_g_col != n) {
            *samerows = i;
            min_c_col = n;
            min_c_wt = 0;
            for (j = 0; j < (size_t)di; ++j) {
                k = ce[vi + j];
                if (ISMARKED1(k) && k < min_c_col) {
                    min_c_col = k;
                    min_c_wt = swg_wmark1[k];
                }
            }
            if (min_c_col < min_g_col) return -1;
            if (min_g_col < min_c_col) return 1;
            return min_g_wt < min_c_wt ? -1 : 1;
        }
    }

    *samerows = n;
    return 0;
}

static void updatecan_swg(graph *g, graph *canong, int *lab, int samerows, int m, int n)
{
    size_t *v, *cv, vli, j, k;
    int *d, *e, *cd, *ce;
    sg_weight *wt, *cwt;
    int i, dli;

    (void)m;
    SWG_VDE(g, v, d, e, wt);
    SWG_VDE(canong, cv, cd, ce, cwt);

    DYNALLOC1(int, swg_work1, swg_work1_sz, (size_t)n, "updatecan_swg");

    ((sparsegraph *)canong)->nv = n;
    ((sparsegraph *)canong)->nde = ((sparsegraph *)g)->nde;

    for (i = 0; i < n; ++i) swg_work1[lab[i]] = i;

    if (samerows == 0) k = 0;
    else k = cv[samerows - 1] + (size_t)cd[samerows - 1];

    for (i = samerows; i < n; ++i) {
        cv[i] = k;
        cd[i] = dli = d[lab[i]];
        vli = v[lab[i]];
        for (j = 0; j < (size_t)dli; ++j) {
            ce[k] = swg_work1[e[vli + j]];
            cwt[k] = wt[vli + j];
            ++k;
        }
    }
}

static void init_swg(
    graph *gin,
    graph **gout,
    graph *hin,
    graph **hout,
    int *lab,
    int *ptn,
    set *active,
    struct optionstruct *options,
    int *status,
    int m,
    int n
)
{
    (void)gout;
    (void)hout;
    (void)lab;
    (void)ptn;
    (void)active;
    (void)m;
    (void)n;

    if (options->getcanon && hin != NULL) {
        sparsegraph *sg = (sparsegraph *)gin;
        sparsegraph *sh = (sparsegraph *)hin;
        SWG_ALLOC(*sh, (size_t)sg->nv, sg->nde, "init_swg");
        sh->nv = sg->nv;
        sh->nde = sg->nde;
    }
    *status = 0;
}

/* Conservative sufficient condition for edge-colored undirected simple graphs:
 * every non-trivial cell has size 2, with at most one cell of size 3. */
static boolean cheapautom_swg(int *ptn, int level, boolean digraph, int n)
{
    int i, k, nnt;

    if (digraph) return FALSE;

    k = n;
    nnt = 0;
    for (i = 0; i < n; ++i) {
        --k;
        if (ptn[i] > level) {
            ++nnt;
            while (ptn[++i] > level) {}
        }
    }

    return k <= nnt + 1;
}

static void split_cell_by_key(
    int *lab,
    int *ptn,
    int level,
    int *numcells,
    set *active,
    int *active_queue,
    int *nactive,
    int *cellstart,
    const int *key,
    long *longcode,
    int v1,
    int n
)
{
    int v2, v3, j, k, l, lj;
    int w1, w2, w3;
    int size, bigsize, bigpos;

    w1 = key[lab[v1]];
    *longcode = MASH(*longcode, v1);

    v2 = v1 + 1;
    while (ptn[v2 - 1] > level && key[lab[v2]] == w1) ++v2;
    if (ptn[v2 - 1] <= level) return;

    w2 = INT_MAX;
    v3 = j = v2;
    do {
        lj = lab[j];
        w3 = key[lj];
        if (w3 == w1) {
            lab[j] = lab[v3];
            lab[v3] = lab[v2];
            lab[v2] = lj;
            ++v2;
            ++v3;
        } else if (w3 == w2) {
            lab[j] = lab[v3];
            lab[v3] = lj;
            ++v3;
        } else if (w3 < w1) {
            lab[j] = lab[v2];
            lab[v2] = lab[v1];
            lab[v1] = lj;
            v3 = v2 + 1;
            v2 = v1 + 1;
            w2 = w1;
            w1 = w3;
        } else if (w3 < w2) {
            lab[j] = lab[v2];
            lab[v2] = lj;
            v3 = v2 + 1;
            w2 = w3;
        }
    } while (ptn[j++] > level);

    *longcode = MASH(*longcode, w1);
    *longcode = MASH(*longcode, v2);
    if (j == v2) return;

    if (v2 == v1 + 1) {
        cellstart[lab[v1]] = n;
    }

    if (v3 == v2 + 1) {
        cellstart[lab[v2]] = n;
    } else {
        for (k = v2; k < v3; ++k) cellstart[lab[k]] = v2;
    }
    ++(*numcells);
    ptn[v2 - 1] = level;

    if (j == v3) {
        *longcode = MASH(*longcode, w2);
        if (v2 - v1 <= v3 - v2 && !ISELEMENT(active, v1)) {
            ADDELEMENT(active, v1);
            active_queue[(*nactive)++] = v1;
        } else {
            ADDELEMENT(active, v2);
            active_queue[(*nactive)++] = v2;
        }
    } else {
        *longcode = MASH(*longcode, w2);
        *longcode = MASH(*longcode, v3);
        sort_indirect_asc(lab + v3, key, j - v3);
        active_queue[(*nactive)++] = v2;
        ADDELEMENT(active, v2);
        if (v2 - v1 >= v3 - v2) {
            bigpos = -1;
            bigsize = v2 - v1;
        } else {
            bigpos = *nactive - 1;
            bigsize = v3 - v2;
            *longcode = MASH(*longcode, bigsize);
        }
        for (k = v3 - 1; k < j - 1;) {
            ptn[k] = level;
            ++(*numcells);
            l = k + 1;
            ADDELEMENT(active, l);
            active_queue[(*nactive)++] = l;
            w3 = key[lab[l]];
            *longcode = MASH(*longcode, w3);
            for (k = l; k < j - 1 && key[lab[k + 1]] == w3; ++k) {
                cellstart[lab[k + 1]] = l;
            }
            size = k - l + 1;
            if (size == 1) {
                cellstart[lab[l]] = n;
            } else {
                cellstart[lab[l]] = l;
                if (size > bigsize) {
                    bigsize = size;
                    bigpos = *nactive - 1;
                }
            }
        }
        if (bigpos >= 0 && !ISELEMENT(active, v1)) {
            DELELEMENT(active, active_queue[bigpos]);
            ADDELEMENT(active, v1);
            active_queue[bigpos] = v1;
        }
    }
}

static void refine_swg(
    graph *g,
    int *lab,
    int *ptn,
    int level,
    int *numcells,
    int *count,
    set *active,
    int *code,
    int m,
    int n
)
{
    size_t *v, vi, ii;
    int *d, *e;
    sg_weight *wt;
    int i, j, k, isplit, splitv, di;
    int nactive, hitcells;
    boolean trivsplit;
    long longcode;

    (void)count;
    SWG_VDE(g, v, d, e, wt);

    DYNALLOC1(int, swg_work1, swg_work1_sz, (size_t)n, "refine_swg");
    DYNALLOC1(int, swg_work2, swg_work2_sz, (size_t)n, "refine_swg");
    DYNALLOC1(int, swg_work3, swg_work3_sz, (size_t)n, "refine_swg");
    DYNALLOC1(int, swg_work4, swg_work4_sz, (size_t)n, "refine_swg");
    prepare_marks1((size_t)n);
    prepare_marks2((size_t)n);

    int *cellstart = swg_work1;
    int *active_queue = swg_work2;
    int *hits = swg_work3;
    int *hitcell = swg_work4;

    longcode = *numcells;

    nactive = 0;
    for (i = -1; (i = nextelement(active, m, i)) >= 0;) {
        active_queue[nactive++] = i;
    }
    if (nactive == 0) {
        *code = CLEANUP(longcode);
        return;
    }

    for (i = 0; i < n;) {
        if (ptn[i] <= level) {
            cellstart[lab[i]] = n;
            ++i;
        } else {
            j = i;
            do {
                cellstart[lab[i]] = j;
            } while (ptn[i++] > level);
        }
    }

    while (nactive > 0 && *numcells < n) {
        for (i = 0; i < nactive && i < 10; ++i) {
            if (ptn[active_queue[i]] <= level) break;
        }

        if (i < nactive && i < 10) {
            trivsplit = TRUE;
            isplit = active_queue[i];
            active_queue[i] = active_queue[--nactive];
        } else {
            isplit = active_queue[--nactive];
            trivsplit = ptn[isplit] <= level;
        }

        DELELEMENT(active, isplit);
        longcode = MASH(longcode, isplit);

        if (trivsplit) {
            RESETMARKS1;
            RESETMARKS2;
            hitcells = 0;
            splitv = lab[isplit];
            vi = v[splitv];
            di = d[splitv];
            for (ii = 0; ii < (size_t)di; ++ii) {
                j = e[vi + ii];
                k = cellstart[j];
                if (k != n) {
                    if (ISNOTMARKED1(k)) {
                        MARK1(k);
                        hitcell[hitcells++] = k;
                    }
                    MARK2(j);
                    hits[j] = wt[vi + ii] + 1;
                }
            }

            if (hitcells > 1) sort_ints_asc(hitcell, hitcells);
            longcode = MASH(longcode, hitcells);

            for (i = 0; i < hitcells; ++i) {
                int v1 = hitcell[i];
                int p = v1;
                do {
                    int lp = lab[p];
                    if (!ISMARKED2(lp)) hits[lp] = 0;
                } while (ptn[p++] > level);

                split_cell_by_key(
                    lab,
                    ptn,
                    level,
                    numcells,
                    active,
                    active_queue,
                    &nactive,
                    cellstart,
                    hits,
                    &longcode,
                    v1,
                    n
                );
            }
        } else {
            size_t total_hit_edges = 0;
            size_t capacity_needed = 0;
            size_t group_start = 0;
            int first_wt = -1;
            boolean multi_weight = FALSE;
            int scan = isplit;

            do {
                capacity_needed += (size_t)d[lab[scan]];
            } while (ptn[scan++] > level);

            if (capacity_needed == 0) continue;
            if (swg_edges_sz < capacity_needed) {
                DYNALLOC1(swg_hit_edge, swg_edges, swg_edges_sz, capacity_needed, "refine_swg edges");
            }

            scan = isplit;
            do {
                int sv = lab[scan];
                vi = v[sv];
                di = d[sv];
                for (ii = 0; ii < (size_t)di; ++ii) {
                    j = e[vi + ii];
                    if (cellstart[j] != n) {
                        int ew = wt[vi + ii];
                        if (first_wt < 0) first_wt = ew;
                        else if (ew != first_wt) multi_weight = TRUE;
                        swg_edges[total_hit_edges].weight = ew;
                        swg_edges[total_hit_edges].target = j;
                        ++total_hit_edges;
                    }
                }
            } while (ptn[scan++] > level);

            if (total_hit_edges == 0) continue;

            if (multi_weight) {
                sort_hit_edges_by_weight(swg_edges, total_hit_edges);
            }

            while (group_start < total_hit_edges && *numcells < n) {
                int cur_wt = swg_edges[group_start].weight;
                size_t group_end = total_hit_edges;
                if (multi_weight) {
                    group_end = group_start + 1;
                    while (group_end < total_hit_edges && swg_edges[group_end].weight == cur_wt) {
                        ++group_end;
                    }
                }

                longcode = MASH(longcode, cur_wt);
                RESETMARKS1;
                hitcells = 0;
                for (ii = group_start; ii < group_end; ++ii) {
                    j = swg_edges[ii].target;
                    k = cellstart[j];
                    if (k != n) {
                        if (ISNOTMARKED1(k)) {
                            int p = k;
                            MARK1(k);
                            hitcell[hitcells++] = k;
                            do {
                                hits[lab[p]] = 0;
                            } while (ptn[p++] > level);
                        }
                        ++hits[j];
                    }
                }

                if (hitcells > 1) sort_ints_asc(hitcell, hitcells);
                longcode = MASH(longcode, hitcells);
                for (i = 0; i < hitcells; ++i) {
                    split_cell_by_key(
                        lab,
                        ptn,
                        level,
                        numcells,
                        active,
                        active_queue,
                        &nactive,
                        cellstart,
                        hits,
                        &longcode,
                        hitcell[i],
                        n
                    );
                }

                group_start = group_end;
            }
        }
    }

    longcode = MASH(longcode, *numcells);
    *code = CLEANUP(longcode);
}

static int bestcell_swg(graph *g, int *lab, int *ptn, int level, int n)
{
    size_t *v, vi, j;
    int *d, *e;
    sg_weight *wt;
    int i, k, di, nnt, best_idx, maxcnt;

    SWG_VDE(g, v, d, e, wt);

    DYNALLOC1(int, swg_work1, swg_work1_sz, (size_t)n, "bestcell_swg");
    DYNALLOC1(int, swg_work2, swg_work2_sz, (size_t)n, "bestcell_swg");
    DYNALLOC1(int, swg_work3, swg_work3_sz, (size_t)n, "bestcell_swg");
    DYNALLOC1(int, swg_work4, swg_work4_sz, (size_t)n, "bestcell_swg");
    DYNALLOC1(int, swg_work5, swg_work5_sz, (size_t)n, "bestcell_swg");
    DYNALLOC1(int, swg_wmark1, swg_wmark1_sz, (size_t)n, "bestcell_swg");

    int *start = swg_work1;
    int *size = swg_work1 + (n / 2);
    int *nntcell = swg_work2;
    int *hits = swg_work3;
    int *cnt = swg_work4;
    int *first_wt = swg_work5;
    int *mixed_wt = swg_wmark1;

    i = 0;
    nnt = 0;
    while (i < n) {
        if (ptn[i] > level) {
            int p = i;
            start[nnt] = i;
            do {
                nntcell[lab[p]] = nnt;
            } while (ptn[p++] > level);
            size[nnt] = p - i;
            ++nnt;
            i = p;
        } else {
            nntcell[lab[i]] = n;
            ++i;
        }
    }

    if (nnt == 0) return n;

    for (i = 0; i < nnt; ++i) {
        hits[i] = 0;
        cnt[i] = 0;
        first_wt[i] = -1;
        mixed_wt[i] = 0;
    }

    for (i = 0; i < nnt; ++i) {
        vi = v[lab[start[i]]];
        di = d[lab[start[i]]];
        for (j = 0; j < (size_t)di; ++j) {
            k = nntcell[e[vi + j]];
            if (k != n) {
                int ew = wt[vi + j];
                if (hits[k] == 0) first_wt[k] = ew;
                else if (first_wt[k] != ew) mixed_wt[k] = 1;
                ++hits[k];
            }
        }
        for (j = 0; j < (size_t)di; ++j) {
            k = nntcell[e[vi + j]];
            if (k != n) {
                if (hits[k] > 0 && (hits[k] < size[k] || mixed_wt[k])) {
                    ++cnt[i];
                }
                hits[k] = 0;
                first_wt[k] = -1;
                mixed_wt[k] = 0;
            }
        }
    }

    best_idx = 0;
    maxcnt = cnt[0];
    for (i = 1; i < nnt; ++i) {
        if (cnt[i] > maxcnt) {
            best_idx = i;
            maxcnt = cnt[i];
        }
    }
    return start[best_idx];
}

static int targetcell_swg(
    graph *g,
    int *lab,
    int *ptn,
    int level,
    int tc_level,
    boolean digraph,
    int hint,
    int m,
    int n
)
{
    int i;
    (void)digraph;
    (void)m;

    if (hint >= 0 && ptn[hint] > level && (hint == 0 || ptn[hint - 1] <= level)) {
        return hint;
    }
    if (level <= tc_level) {
        return bestcell_swg(g, lab, ptn, level, n);
    }
    for (i = 0; i < n && ptn[i] <= level; ++i) {}
    return i == n ? 0 : i;
}

static void swg_freedyn(void)
{
    DYNFREE(swg_mark1, swg_mark1_sz);
    DYNFREE(swg_mark2, swg_mark2_sz);
    DYNFREE(swg_wmark1, swg_wmark1_sz);
    DYNFREE(swg_work1, swg_work1_sz);
    DYNFREE(swg_work2, swg_work2_sz);
    DYNFREE(swg_work3, swg_work3_sz);
    DYNFREE(swg_work4, swg_work4_sz);
    DYNFREE(swg_work5, swg_work5_sz);
    DYNFREE(swg_edges, swg_edges_sz);
    DYNFREE(swg_snwork, swg_snwork_sz);
}

dispatchvec umol_dispatch_swg = {
    isautom_swg,
    testcanlab_swg,
    updatecan_swg,
    refine_swg,
    refine_swg,
    cheapautom_swg,
    targetcell_swg,
    swg_freedyn,
    nausparse_check,
    init_swg,
    NULL
};

void umol_sparsenauty_weighted(
    sparsegraph *sg,
    int *lab,
    int *ptn,
    int *orbits,
    optionblk *options,
    statsblk *stats,
    sparsegraph *h
)
{
    int n = sg->nv;
    int m = SETWORDSNEEDED(n);
    boolean saved_digraph = options->digraph;

    if (!saved_digraph && sg->nde > 0) {
        size_t *v;
        int *d, *e;
        int i;
        SG_VDE(sg, v, d, e);
        for (i = 0; i < n; ++i) {
            size_t vi = v[i];
            int di = d[i];
            size_t j;
            for (j = 0; j < (size_t)di; ++j) {
                if (e[vi + j] == i) {
                    options->digraph = TRUE;
                    break;
                }
            }
            if (options->digraph) break;
        }
    }

    DYNALLOC1(set, swg_snwork, swg_snwork_sz, (size_t)(2 * 500 * m), "umol_sparsenauty_weighted");
    nauty(
        (graph *)sg,
        lab,
        ptn,
        NULL,
        orbits,
        options,
        stats,
        swg_snwork,
        2 * 500 * m,
        m,
        n,
        (graph *)h
    );
    options->digraph = saved_digraph;
}
