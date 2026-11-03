/* SPDX-License-Identifier: GPL-3.0-only */
#ifndef SX_DB_H
#define SX_DB_H
#include <sqlite3.h>
#include "sx_common.h"
long sx_db_init(sqlite3 **db, const char *filename);
long sx_db_close(sqlite3 **db);
long sx_db_write_stree(sqlite3 **db, int n_leaves, double sb_rate, double sd_rate, int bds_leaves, double bds_length, double height, double length, double outgroup, int ind_per_sp, int n_loci, double b_rate, double d_rate, double t_rate, double gc_rate, double alpha_s, double alpha_l, double salpha_g, int Ne, double mu, double gen_time);
long sx_db_write_ltree(sqlite3 **db, int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, int n_trials, double gamma, double lalpha_g);
long sx_db_write_gtree(sqlite3 **db, int n_gtree, unsigned long long LID, int n_ltree, int SID, double alpha_g, int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss);
/* Text form of the INSERT statements (for deferred/parallel writing). */
void sx_db_ltree_sql(char *query, size_t n, int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, int n_trials, double gamma, double lalpha_g);
void sx_db_gtree_sql(char *query, size_t n, int n_gtree, unsigned long long LID, int n_ltree, int SID, double alpha_g, int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss);
long sx_db_exec(sqlite3 **db, const char *query);
#endif
