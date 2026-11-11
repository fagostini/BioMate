#!/usr/bin/env bash
# Benchmark `biomate nspector` on real FASTQ files.
#
# Usage: scripts/bench_nspector.sh OUTDIR FASTQ [FASTQ ...]
# Env:   ENGINES="rust python" (default: rust)   THREADS=0   EXTRA="--no-multiqc"
#
# For each engine prints wall time, reads scanned, reads/s, peak RSS and
# report size. Results go under OUTDIR/<engine>/.
set -euo pipefail
export LC_NUMERIC=C

if [[ $# -lt 2 ]]; then
    sed -n '2,8p' "$0"
    exit 1
fi
outdir=$1
shift
engines=${ENGINES:-rust}
threads=${THREADS:-0}
biomate=${BIOMATE:-biomate}

printf '%-7s %10s %14s %12s %12s %12s\n' engine wall_s reads Mreads_s peak_rss_MB report_MB
for engine in $engines; do
    out="$outdir/$engine"
    mkdir -p "$out"
    log="$out/time.log"
    # shellcheck disable=SC2086
    /usr/bin/time -v -o "$log" "$biomate" nspector --input "$@" --output "$out" \
        --engine "$engine" --threads "$threads" ${EXTRA:-} > "$out/nspector.log" 2>&1
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $2}' "$log" |
        awk -F: '{ if (NF == 3) print $1*3600 + $2*60 + $3; else print $1*60 + $2 }')
    rss_kb=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$log")
    reads=$(python3 -c 'import polars, sys; print(polars.read_parquet(sys.argv[1])["total"].sum())' \
        "$out/nspector_sample_lane.parquet")
    report=0
    [[ -f "$out/nspector_report.html" ]] && report=$(stat -c %s "$out/nspector_report.html")
    # a missing biomate_kernels makes --engine rust fall back silently to Python
    if grep -q "falling back to the Python scan" "$out/nspector.log"; then
        engine="$engine(->python!)"
    fi
    awk -v e="$engine" -v w="$wall" -v r="$reads" -v m="$rss_kb" -v s="$report" 'BEGIN {
        printf "%-7s %10.1f %14d %12.2f %12.0f %12.1f\n", e, w, r, r / w / 1e6, m / 1024, s / 1048576 }'
done
