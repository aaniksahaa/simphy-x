/* Species tree simulation: port of SimPhy's NewBDSTree (BDSA with a fixed
 * number of leaves and optional fixed height; SSA with a fixed height) plus
 * outgroup addition.  Random draws happen in exactly SimPhy's order. */
#ifndef SX_STREE_H
#define SX_STREE_H
#include "sx_tree.h"
#include "sx_rng.h"

long sx_new_bds_tree(sx_stree **out_tree, int leaves, double time, double b_rate, double d_rate, double gen_time, int Ne, double mu,
                     int ind_per_sp, double outgroup, int complete, int mrca_time, gsl_rng *seed, int verbosity);

/* Rateheter_lineagespec: one gamma draw per node in array order. */
void sx_stree_rateheter(sx_stree *t, double alpha, gsl_rng *seed);

#endif
