/* Output writers reproducing SimPhy's file formats byte for byte. */
#ifndef SX_WRITE_H
#define SX_WRITE_H
#include "sx_tree.h"

typedef struct { char **names; int n_names; } sx_names;   /* NULL when trees are simulated */

void sx_write_stree(sx_buf *b, const sx_stree *t, const sx_names *names, int time_units, int int_labels);
void sx_write_ltree(sx_buf *b, const sx_ltree *t, const sx_names *names, int time_units, int int_labels);
void sx_write_gtree(sx_buf *b, const sx_gtree *t, const sx_names *names, int int_labels);
/* bounded_locus_subtrees.out line for one locus (WriteDaughtersFile incl. the
 * trailing-comma → newline fix-up). */
void sx_write_daughters(sx_buf *b, const sx_ltree *t, const sx_names *names);
long sx_write_mapping_sl(sx_buf *b, sx_stree *st, const sx_ltree *lt, const sx_names *names);
long sx_write_mapping_lg(sx_buf *b, const sx_gtree *gt, const sx_names *names);

#endif
