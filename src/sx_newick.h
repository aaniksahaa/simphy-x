/* Fixed species trees given as (extended) Newick strings: port of SimPhy's
 * CheckNewickSTree / ParseNewickSTree / RefineSNodes.  Supported annotations
 * after a node: :length (generations), #Ne, /individuals, *mu_mult, ~gtime_mult. */
#ifndef SX_NEWICK_H
#define SX_NEWICK_H
#include "sx_tree.h"
#include "sx_write.h"

sx_names *sx_names_new(int n_names);
void sx_names_free(sx_names *n);

/* Returns a spread (uncollapsed) species tree with post-order index labels, or aborts on error. */
sx_stree *sx_parse_newick_stree(const char *newick, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp);

#endif
