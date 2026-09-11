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
/* Nexus flavour: annotations as [&pop_size=,n_ind=,u_mult=,g_mult=] comments. */
sx_stree *sx_parse_nexus_stree(const char *nexus, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp);

/* Fixed locus trees (SimPhy's ParseNewickLTree / ParseNexusLTree).  Extended Newick
 * codes: :length *mu ~gtime #Ne /n_ind %kind _paralog ; Nexus comments:
 * [&pop_size= kind_n= paralog= n_ind= u_mult= g_mult=].  The returned tree is
 * post-order collapsed with the lineage-count probabilities computed
 * (compat arithmetic).  n_dup/n_loss/n_trans/n_gc are incremented as SimPhy does
 * (losses also count as gene conversions, a fall-through in RefineLNodes). */
sx_ltree *sx_parse_newick_ltree(const char *newick, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc);
sx_ltree *sx_parse_nexus_ltree(const char *nexus, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc);

/* Nexus tree files: count the trees (NNexusTrees), position after "BEGIN TREES;"
 * (InitNexusParser), read the next tree string from "(" to end of line (NextNexusTree). */
long sx_nexus_count(FILE *f, int *n_trees);
long sx_nexus_init(FILE *f);
long sx_nexus_next(FILE *f, char **tree_str);

#endif
