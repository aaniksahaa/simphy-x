#!/usr/bin/env python3
"""Statistical equivalence of simphy-x fast mode vs the reference (compat).

Runs the same scenario with many replicates in both modes (different random
streams by construction), collects per-locus and per-gene-tree summary
statistics from the SQLite databases and the Newick files, and compares the
distributions with two-sample Kolmogorov-Smirnov tests (continuous) and
chi-square tests (counts).  Prints a table and exits non-zero if any test
fails after a Bonferroni correction.

Usage: tests/verify_fast.py [--taxa N] [--loci L] [--reps R] [--pop P] [--dup D] [--loss-ratio r] [--work DIR]
"""
import argparse, os, subprocess, sqlite3, sys, math, re, shutil
import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.environ.get("SIMPHY_REF", os.path.join(ROOT, "upstream/SimPhy/bin/simphy_ref"))
NEW = os.environ.get("SIMPHY_X", os.path.join(ROOT, "bin/simphy-x"))

def run(binary, args, out, extra=(), env=''):
    if os.path.isdir(out): shutil.rmtree(out)
    cmd = [binary] + args + ["-o", out] + list(extra)
    e = dict(os.environ)
    for kv in env.split(','):
        if '=' in kv: k, v = kv.split('=', 1); e[k] = v
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=e)
    if r.returncode != 0:
        print(r.stdout[-3000:]); raise SystemExit(f"{binary} failed")

def newick_stats(path):
    """per gene tree: n_leaves, height (max root-to-tip bl), total bl, n_cherries"""
    out = []
    for line in open(path):
        line = line.strip()
        if not line: continue
        # tokenise
        toks = re.findall(r"\(|\)|,|[^(),:;]+:[-+0-9.eE]+|;", line)
        depth = 0; stack = []  # (children_count, cherry flag)
        heights = []; total = 0.0
        # iterative parse with explicit stack of node records: [n_children, all_children_leaves, height_children_max]
        nodes = [[0, True, 0.0]]
        leaves = 0; cherries = 0
        i = 0
        cur_len = None
        for t in toks:
            if t == '(':
                nodes.append([0, True, 0.0])
            elif t == ',':
                pass
            elif t == ')':
                pass
            elif t == ';':
                pass
            else:
                name, ln = t.rsplit(':', 1)
                ln = float(ln); total += ln
                if name.startswith(')'):
                    node = nodes.pop()
                    if node[0] == 2 and node[1]: cherries += 1
                    h = node[2] + ln
                    par = nodes[-1]; par[0] += 1; par[1] = False; par[2] = max(par[2], h)
                else:
                    leaves += 1
                    par = nodes[-1]; par[0] += 1; par[2] = max(par[2], ln)
        # root closing ");" has no length: handle remaining node
        # tokens like ")" followed by ":len" were merged as name starting with ')'
        root = nodes.pop() if len(nodes) > 1 else nodes[0]
        if root[0] == 2 and root[1]: cherries += 1
        out.append((leaves, root[2], total, cherries))
    return out

def collect(outdir, reps):
    db = [f for f in os.listdir(outdir) if f.endswith('.db')][0]
    c = sqlite3.connect(os.path.join(outdir, db))
    lt = np.array(c.execute("select n_leaves, n_dup, n_loss from Locus_Trees").fetchall(), dtype=float)
    gt = np.array(c.execute("select Extra_l, Tree_h_cu, Tree_l_bl, n_leaves from Gene_Trees").fetchall(), dtype=float)
    nw = []
    for d in sorted(os.listdir(outdir)):
        p = os.path.join(outdir, d)
        if not os.path.isdir(p): continue
        for f in sorted(os.listdir(p)):
            if f.startswith('g_trees'): nw.extend(newick_stats(os.path.join(p, f)))
    nw = np.array(nw, dtype=float)
    return lt, gt, nw

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--taxa', type=int, default=20); ap.add_argument('--loci', type=int, default=40); ap.add_argument('--reps', type=int, default=25)
    ap.add_argument('--pop', type=float, default=290000000); ap.add_argument('--dup', type=float, default=2.046e-10); ap.add_argument('--loss-ratio', type=float, default=0.1)
    ap.add_argument('--work', default='/tmp/simphy-x-verify-fast'); ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--alpha', type=float, default=0.01)
    ap.add_argument('--ref-mode', default='compat'); ap.add_argument('--fast-mode', default='fast')
    ap.add_argument('--ref-env', default=''); ap.add_argument('--fast-env', default='')
    ap.add_argument('--label', default='')
    ap.add_argument('--stree', default='', help='Newick file: use this fixed species tree (-S) instead of simulating one')
    ap.add_argument('--simple-rates', action='store_true', help='no lineage/gene-family rate heterogeneity and fixed mu, so substitution-unit statistics are iid per locus')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    args = ["-sl", f"f:{a.taxa}", "-rs", str(a.reps), "-rl", f"f:{a.loci}", "-rg", "1", "-sb", "f:0.000001", "-sd", "f:0", "-st", "ln:21.25,0.2", "-so", "f:1", "-si", "f:1",
            "-sp", f"f:{int(a.pop)}", "-su", "ln:-21.9,0.1", "-hh", "f:1", "-hs", "ln:1.5,1", "-hl", "ln:1.551533,0.6931472", "-hg", "ln:1.5,1",
            "-lb", f"f:{a.dup}", "-ld", f"f:{a.dup * a.loss_ratio}", "-lt", "f:0", "-cs", "9644", "-v", "0", "-ot", "0", "-od", "1", "-oc", "0"]
    if a.simple_rates:
        na = []; skip = False
        for tok in args:
            if skip: skip = False; continue
            if tok in {"-hs", "-hl", "-su"}: skip = True; continue
            na.append(tok)
        args = na + ["-su", "f:1e-9"]
    if a.stree:
        nw = open(a.stree).read().strip().split('\n')[0].strip()
        # drop species-tree simulation options, add the fixed tree
        drop = {"-sl", "-sb", "-sd", "-st", "-so"}
        na = []; skip = False
        for tok in args:
            if skip: skip = False; continue
            if tok in drop: skip = True; continue
            na.append(tok)
        args = ["-s", nw] + na
    ref_out = os.path.join(a.work, "ref"); fast_out = os.path.join(a.work, "fast")
    print(f"running A: mode {a.ref_mode} env [{a.ref_env}] ..."); run(NEW, args, ref_out, ["-XM", a.ref_mode, "-XT", str(a.threads)], a.ref_env)
    print(f"running B: mode {a.fast_mode} env [{a.fast_env}] ..."); run(NEW, args, fast_out, ["-XM", a.fast_mode, "-XT", str(a.threads)], a.fast_env)
    lt1, gt1, nw1 = collect(ref_out, a.reps); lt2, gt2, nw2 = collect(fast_out, a.reps)
    tests = []
    def ks(name, x, y): tests.append((name, 'KS', stats.ks_2samp(x, y).pvalue, len(x), len(y), np.mean(x), np.mean(y)))
    def chi(name, x, y):
        vals = np.union1d(x, y)
        # pool sparse categories
        cx = np.array([np.sum(x == v) for v in vals]); cy = np.array([np.sum(y == v) for v in vals])
        keep = (cx + cy) >= 10
        if keep.sum() < 2: return ks(name, x, y)
        cx2 = np.append(cx[keep], cx[~keep].sum()); cy2 = np.append(cy[keep], cy[~keep].sum())
        tab = np.array([cx2, cy2]); tab = tab[:, tab.sum(0) > 0]
        p = stats.chi2_contingency(tab)[1] if tab.shape[1] > 1 else 1.0
        tests.append((name, 'chi2', p, len(x), len(y), np.mean(x), np.mean(y)))
    chi('locus n_leaves', lt1[:, 0], lt2[:, 0]); chi('locus n_dup', lt1[:, 1], lt2[:, 1]); chi('locus n_loss', lt1[:, 2], lt2[:, 2])
    chi('gene extra lineages', gt1[:, 0], gt2[:, 0])
    ks('gene height (cu)', gt1[:, 1], gt2[:, 1]); ks('gene total length (bl)', gt1[:, 2], gt2[:, 2])
    ks('newick height (bl)', nw1[:, 1], nw2[:, 1]); chi('newick cherries', nw1[:, 3], nw2[:, 3])
    # conditional on duplications present: heights of MLC trees
    m1 = lt1[:, 1] > 0; m2 = lt2[:, 1] > 0
    if m1.sum() > 20 and m2.sum() > 20:
        ks('gene height (cu) | dups>0', gt1[m1, 1], gt2[m2, 1]); ks('gene total length | dups>0', gt1[m1, 2], gt2[m2, 2]); chi('extra lineages | dups>0', gt1[m1, 0], gt2[m2, 0])
    bonf = a.alpha / len(tests)
    fail = 0
    print(f"\n{'statistic':32s} {'test':5s} {'p-value':>9s} {'n_ref':>6s} {'n_fast':>6s} {'mean_ref':>12s} {'mean_fast':>12s}")
    for name, kind, p, n1, n2, m1_, m2_ in tests:
        flag = '' if p >= bonf else '  <-- FAIL'
        if p < bonf: fail += 1
        print(f"{name:32s} {kind:5s} {p:9.4f} {n1:6d} {n2:6d} {m1_:12.6g} {m2_:12.6g}{flag}")
    print(f"\n{a.label} {len(tests)} tests, Bonferroni threshold {bonf:.4g}, failures: {fail}")
    sys.exit(1 if fail else 0)

if __name__ == '__main__':
    main()
