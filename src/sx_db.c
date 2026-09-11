#include "sx_db.h"

long sx_db_exec(sqlite3 **db, const char *query)
{
    char *zErrMsg = NULL;
    if (sqlite3_exec(*db, query, NULL, NULL, &zErrMsg) != SQLITE_OK) {
        fprintf(stderr, "\n\nCan't write into the db: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        sqlite3_close(*db);
        return SX_DB_ERROR;
    }
    return SX_NO_ERROR;
}

long sx_db_init(sqlite3 **database, const char *db_filename)
{
    char *zErrMsg = NULL;
    if (sqlite3_open(db_filename, database) != SQLITE_OK) {
        fprintf(stderr, "\n\nCan't open database: %s\n", sqlite3_errmsg(*database));
        sqlite3_close(*database);
        return SX_DB_ERROR;
    }
    sqlite3_exec(*database, "DROP TABLE IF EXISTS Species_Trees; DROP TABLE IF EXISTS Locus_Trees; DROP TABLE IF EXISTS Gene_Trees; PRAGMA foreign_keys=yes", NULL, NULL, &zErrMsg);
    const char *q[] = {
        "CREATE TABLE IF NOT EXISTS Species_Trees(SID INTEGER PRIMARY KEY,Leaves INTEGER,SB_rate REAL,SD_rate REAL,F_length_gen REAL,Height_cu REAL,Length_cu REAL, Outgroup_intlength_gen REAL, Ind_per_sp INTEGER,N_loci INTEGER,GB_rate REAL,GD_rate REAL,GT_rate REAL,GGC_rate REAL,Alpha_s REAL,Alpha_l REAL,Salpha_g REAL,Ne INTEGER,Mu REAL,Gen_time REAL)",
        "CREATE TABLE IF NOT EXISTS Locus_Trees(LID INTEGER PRIMARY KEY,n_ltree INTEGER NOT NULL,SID INTEGER NOT NULL,b_rate REAL,d_rate REAL,t_rate REAL,gc_rate REAL,n_leaves INTEGER,n_dup INTEGER,n_loss INTEGER, n_transf INTEGER, n_gc INTEGER, n_ltrials INTEGER, gamma REAL, Lalpha_g REAL, FOREIGN KEY(SID) REFERENCES Species_Trees(SID))",
        "CREATE TABLE IF NOT EXISTS Gene_Trees(GID INTEGER PRIMARY KEY,n_gtree INTEGER NOT NULL,LID INTEGER NOT NULL,n_ltree INTEGER NOT NULL, SID INTEGER NOT NULL,Alpha_g REAL,n_leaves INTEGER,Extra_l INTEGER,Tree_h_cu REAL,Tree_l_bl REAL,FOREIGN KEY(LID) REFERENCES Locus_Trees(LID),FOREIGN KEY(SID) REFERENCES Species_Trees(SID))",
        "BEGIN IMMEDIATE TRANSACTION"
    };
    for (int i = 0; i < 4; ++i) {
        if (sqlite3_exec(*database, q[i], NULL, NULL, &zErrMsg) != SQLITE_OK) {
            fprintf(stderr, "\n\nCan't create the tables: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
            sqlite3_close(*database);
            return SX_DB_ERROR;
        }
    }
    return SX_NO_ERROR;
}

long sx_db_close(sqlite3 **database)
{
    char *zErrMsg = NULL;
    if (sqlite3_exec(*database, "COMMIT TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        fprintf(stderr, "\n\nCan't commit the results to the db: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        sqlite3_close(*database);
        return SX_DB_ERROR;
    }
    sqlite3_close(*database);
    return SX_NO_ERROR;
}

long sx_db_write_stree(sqlite3 **db, int n_leaves, double sb_rate, double sd_rate, int bds_leaves, double bds_length, double height, double length, double outgroup, int ind_per_sp, int n_loci, double b_rate, double d_rate, double t_rate, double gc_rate, double alpha_s, double alpha_l, double salpha_g, int Ne, double mu, double gen_time)
{
    char query[1000];
    (void)bds_leaves;
    snprintf(query, sizeof query, "INSERT INTO Species_Trees VALUES(NULL,'%d','%e','%e','%e','%e','%e','%e','%d','%d','%e','%e','%e','%e','%e','%e','%e','%d','%e','%e')", n_leaves, sb_rate, sd_rate, bds_length, height, length, outgroup, ind_per_sp, n_loci, b_rate, d_rate, t_rate, gc_rate, alpha_s, alpha_l, salpha_g, Ne, mu, gen_time);
    return sx_db_exec(db, query);
}

void sx_db_ltree_sql(char *query, size_t n, int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, int n_trials, double gamma, double lalpha_g)
{
    snprintf(query, n, "INSERT INTO Locus_Trees VALUES(NULL,'%d','%d','%e','%e','%e','%e','%d','%d','%d','%d','%d','%d','%e','%e')", n_ltree, SID, b_rate, d_rate, t_rate, gc_rate, n_leaves, n_dup, n_loss, n_transf, n_gc, n_trials, gamma, lalpha_g);
}
long sx_db_write_ltree(sqlite3 **db, int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, int n_trials, double gamma, double lalpha_g)
{
    char query[1000];
    sx_db_ltree_sql(query, sizeof query, n_ltree, SID, b_rate, d_rate, t_rate, gc_rate, n_leaves, n_dup, n_loss, n_transf, n_gc, n_trials, gamma, lalpha_g);
    return sx_db_exec(db, query);
}
void sx_db_gtree_sql(char *query, size_t n, int n_gtree, unsigned long long LID, int n_ltree, int SID, double alpha_g, int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss)
{
    snprintf(query, n, "INSERT INTO Gene_Trees VALUES(NULL,'%d','%lld','%d','%d','%e','%d','%d','%e','%e')", n_gtree, LID, n_ltree, SID, alpha_g, n_leaves, extra_lineages, tree_h_cu, tree_l_ss);
}
long sx_db_write_gtree(sqlite3 **db, int n_gtree, unsigned long long LID, int n_ltree, int SID, double alpha_g, int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss)
{
    char query[1000];
    sx_db_gtree_sql(query, sizeof query, n_gtree, LID, n_ltree, SID, alpha_g, n_leaves, extra_lineages, tree_h_cu, tree_l_ss);
    return sx_db_exec(db, query);
}
