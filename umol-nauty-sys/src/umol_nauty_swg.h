/* Weighted sparse-graph dispatch interface for nauty over sparsegraph.w.
 *
 * Adapted from nauty 2.9.3 `nausparse.h` (Copyright (1984-) Brendan McKay,
 * Apache License, Version 2.0; see ../nauty/COPYRIGHT and ../LICENSE-APACHE). */

#ifndef UMOL_NAUTY_SWG_H
#define UMOL_NAUTY_SWG_H

#include "nausparse.h"

extern dispatchvec umol_dispatch_swg;

/* Run nauty with the weighted sparse dispatch vector on an undirected simple
 * sparsegraph whose edge weights are non-negative integers in [0, INT_MAX - 1]. */
void umol_sparsenauty_weighted(
    sparsegraph *sg,
    int *lab,
    int *ptn,
    int *orbits,
    optionblk *options,
    statsblk *stats,
    sparsegraph *h
);

#endif
