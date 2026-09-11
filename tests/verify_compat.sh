#!/usr/bin/env bash
# Byte-exact verification of simphy-x (compat mode) against the reference
# SimPhy 1.0.2 build.  Usage: tests/verify_compat.sh [workdir] [threads...]
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${SIMPHY_REF:-$ROOT/upstream/SimPhy/bin/simphy_ref}"
NEW="${SIMPHY_X:-$ROOT/bin/simphy-x}"
WORK="${1:-/tmp/simphy-x-verify}"
shift || true
THREADS=("$@"); [[ ${#THREADS[@]} -eq 0 ]] && THREADS=(1)
mkdir -p "$WORK"
pass=0; fail=0; failed_cases=()

# name|args
CASES=(
"msc_small|-sb f:0.000001 -ld f:0 -lb f:0 -lt f:0 -rs 2 -rl f:10 -rg 1 -sp f:10000 -su f:0.00001 -sg f:1 -sl f:20 -st f:1000000 -om 1 -v 2 -od 1 -op 1 -oc 1 -on 1 -cs 22 -ll 20 -ls 20"
"msc_500x50|-sb f:0.000001 -ld f:0 -lb f:0 -lt f:0 -rs 1 -rl f:50 -rg 1 -sp f:10000 -su f:0.00001 -sg f:1 -sl f:500 -st f:1000000 -om 1 -v 1 -od 1 -op 1 -oc 1 -on 1 -cs 22 -ll 500 -ls 500"
"msc_distr|-sb f:0.000001 -ld f:0 -lb f:0 -lt f:0 -rs 3 -rl u:5,12 -rg 2 -sp u:500000,1500000 -su ln:-17.27461,0.6931472 -sg f:1 -sl f:30 -st ln:16.2,1 -om 1 -v 1 -od 1 -op 1 -oc 1 -on 1 -cs 42"
"msc_outgroup_ind|-sb f:0.000001 -rs 2 -rl f:6 -rg 1 -sp f:20000 -su f:0.00001 -sl f:12 -st f:500000 -so f:0.5 -si f:3 -ol 1 -ot 1 -os 1 -om 1 -od 1 -op 1 -v 1 -cs 7"
"msc_ssa|-sb f:0.00001 -sd f:0.000002 -rs 2 -rl f:4 -rg 1 -sp f:20000 -su f:0.00001 -st f:600000 -so f:1 -v 1 -om 1 -os 1 -od 1 -cs 11"
"msc_ssa_noog|-sb f:0.00001 -sd f:0.000002 -rs 2 -rl f:4 -rg 1 -sp f:20000 -su f:0.00001 -st f:600000 -v 1 -om 1 -od 1 -cs 13"
"msc_origin|-sb f:0.000001 -sd f:0.0000002 -rs 2 -rl f:5 -rg 1 -sp f:20000 -su f:0.00001 -sl f:15 -v 1 -om 1 -od 1 -cs 5"
"msc_hetero|-sb f:0.000001 -rs 2 -rl f:6 -rg 2 -sp f:20000 -su f:0.00001 -sl f:12 -st f:500000 -hs f:2 -hl f:1.5 -hg f:1 -gp f:2 -hh ln:0.5,0.3 -v 1 -od 1 -os 1 -cs 3"
"gdl_small|-sl f:10 -rs 2 -rl f:5 -rg 1 -sb f:0.000001 -sd f:0 -st ln:21.25,0.2 -so f:1 -si f:1 -sp f:50000000 -su ln:-21.9,0.1 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -lb f:2.046e-10 -ld f:2.046e-11 -lt f:0 -cs 9644 -v 3 -ot 0 -op 1"
"gdl_maps|-sl f:12 -rs 2 -rl f:6 -rg 2 -sb f:0.000001 -st ln:21.25,0.2 -so f:1 -sp f:50000000 -su ln:-21.9,0.1 -lb f:4e-10 -ld f:4e-10 -om 1 -ol 1 -on 1 -os 1 -od 1 -op 1 -cs 17 -v 1"
"gdl_ind2|-sl f:8 -rs 2 -rl f:6 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:200000 -su f:0.00001 -si f:2 -lb f:0.0000005 -ld f:0.0000003 -om 1 -on 1 -cs 23 -v 1"
"gdl_hgt_k1|-sl f:10 -rs 2 -rl f:8 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -lb f:0.0000005 -ld f:0.0000003 -lt f:0.0000005 -lk 1 -om 1 -on 1 -ol 1 -od 1 -os 1 -cs 29 -v 1"
"gdl_hgt_k0|-sl f:10 -rs 2 -rl f:8 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -lb f:0.0000005 -ld f:0.0000003 -lt f:0.0000005 -lk 0 -om 1 -on 1 -cs 31 -v 1"
"gdl_gc|-sl f:10 -rs 2 -rl f:8 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -lb f:0.0000005 -ld f:0.0000003 -lg f:0.0000005 -lt f:0.0000002 -om 1 -on 1 -od 1 -cs 37 -v 1"
"gdl_genomewide|-sl f:10 -rs 2 -rl f:6 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -gb ln:-14.5,0.5 -gd f:0.0000002 -lb sl:0,0.2,gb -ld f:gd -om 1 -on 1 -od 1 -cs 41 -v 1"
"gdl_minleaves|-sl f:10 -rs 2 -rl f:6 -rg 1 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -lb f:0.0000005 -ld f:0.0000009 -ll 5 -ls 3 -om 1 -on 1 -od 1 -cs 43 -v 1"
"fixed_msc|-s (((A:1000000,B:1000000):500000,(C:800000,D:800000):700000):300000,E:1800000); -rs 2 -rl f:8 -rg 1 -sp f:20000 -su f:0.00001 -om 1 -od 1 -op 1 -on 1 -v 1 -cs 61"
"fixed_msc_annot|-s (((A:1000000#30000/2,B:1000000):500000#40000,(C:800000,D:800000/3):700000):300000,E:1800000); -rs 1 -rl f:6 -rg 1 -sp f:20000 -su f:0.00001 -si f:2 -om 1 -od 1 -op 1 -ol 1 -v 1 -cs 67"
"fixed_gdl|-s (((A:1000000,B:1000000):500000,(C:800000,D:800000):700000):300000,E:1800000); -rs 2 -rl f:10 -rg 1 -sp f:100000 -su f:0.00001 -lb f:0.0000008 -ld f:0.0000004 -lt f:0.0000003 -om 1 -on 1 -od 1 -os 1 -op 1 -v 1 -cs 71"
"gdl_v4|-sl f:8 -rs 2 -rl f:4 -rg 2 -sb f:0.000001 -st f:1000000 -sp f:100000 -su f:0.00001 -lb f:0.0000006 -ld f:0.0000003 -hs f:2 -hl f:1.5 -hg f:1 -om 1 -on 1 -cs 53 -v 4"
"msc_v0|-sb f:0.000001 -rs 2 -rl f:5 -rg 1 -sp f:20000 -su f:0.00001 -sl f:12 -st f:500000 -v 0 -cs 59"
"gdl_ot1|-sl f:8 -rs 1 -rl f:5 -rg 1 -sb f:0.000002 -st f:800000 -sg f:2.5 -sp f:100000 -su f:0.00001 -lb f:0.0000005 -ld f:0.0000003 -ot 1 -ol 1 -cs 47 -v 1"
)

compare_dir() {
    local a="$1" b="$2" na="$3" nb="$4"
    # tree/mapping/etc files: exact
    if ! diff -r -q -x "*.command" -x "*.params" -x "*.db" "$a" "$b" >/dev/null; then
        diff -r -q -x "*.command" -x "*.params" -x "*.db" "$a" "$b" | head -5
        return 1
    fi
    # params: identical except the prefix line
    if [[ -f "$a/$na.params" ]]; then
        if ! diff <(grep -v 'Output files prefix' "$a/$na.params") <(grep -v 'Output files prefix' "$b/$nb.params") >/dev/null; then
            echo "params differ"; return 1
        fi
    fi
    if [[ -f "$a/$na.db" ]]; then
        python3 - "$a/$na.db" "$b/$nb.db" <<'PY' || return 1
import sqlite3,sys
def dump(p):
    c=sqlite3.connect(p)
    return [(t,[tuple(r) for r in c.execute('select * from %s order by 1'%t)]) for t in ['Species_Trees','Locus_Trees','Gene_Trees']]
a,b=dump(sys.argv[1]),dump(sys.argv[2])
if a!=b:
    print("db differs"); sys.exit(1)
PY
    fi
    return 0
}

for entry in "${CASES[@]}"; do
    name="${entry%%|*}"; args="${entry#*|}"
    rm -rf "$WORK/${name}_ref" "$WORK/${name}_new"*
    ( cd "$WORK" && $REF $args -o "${name}_ref" >"${name}_ref.log" 2>&1 ); rc_ref=$?
    ok=1
    for t in "${THREADS[@]}"; do
        ( cd "$WORK" && $NEW $args -o "${name}_new" -XT "$t" >"${name}_new_t$t.log" 2>&1 ); rc_new=$?
        if [[ $rc_ref -ne $rc_new ]]; then echo "[$name t=$t] exit codes differ: ref=$rc_ref new=$rc_new"; ok=0; fi
        if [[ $rc_ref -eq 0 ]] && ! compare_dir "$WORK/${name}_ref" "$WORK/${name}_new" "${name}_ref" "${name}_new"; then ok=0; echo "[$name t=$t] MISMATCH"; fi
        if [[ -n "${SX_CHECK_STDOUT:-}" && $rc_ref -eq 0 ]]; then
            # program text after the banner separator; the output-directory name differs by construction
            strip() { awk 'f{print} /^-{20,}$/{f=1}' "$1" | sed "s/${name}_ref/OUT/g; s/${name}_new/OUT/g"; }
            if ! diff <(strip "$WORK/${name}_ref.log") <(strip "$WORK/${name}_new_t$t.log") > "$WORK/${name}_stdout_t$t.diff"; then ok=0; echo "[$name t=$t] STDOUT DIFFERS ($(wc -l < "$WORK/${name}_stdout_t$t.diff") diff lines)"; fi
        fi
        rm -rf "$WORK/${name}_new"
    done
    if [[ $ok -eq 1 ]]; then echo "[$name] OK (ref rc=$rc_ref)"; pass=$((pass+1)); else fail=$((fail+1)); failed_cases+=("$name"); fi
done
echo "passed=$pass failed=$fail ${failed_cases[*]:-}"
[[ $fail -eq 0 ]]
