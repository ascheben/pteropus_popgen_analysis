#!/usr/bin/env bash
# Generate fastsimcoal2 input files (.tpl, .est, _jointMAFpop1_0.obs) for every
# population pair x demographic model, for a given mutation rate.
#
# Usage (from the 5_Fastsimcoal directory):
#   scripts/make_fsc_inputs.sh <mutation_rate> <out_dir>
#
# Examples (inputs used in the manuscript):
#   scripts/make_fsc_inputs.sh 2.5e-8 runs/human_rate    # human rate (main results)
#   scripts/make_fsc_inputs.sh 8.8e-9 runs/mammal_rate   # mammal rate (2.2e-9/yr x 4 yr generation)
#
# Sources:
#   pairs.tsv               pair name and haploid sample sizes (n0, n1) as used in the .tpl
#   templates/<model>.tpl   model skeletons; @N0@, @N1@ and @MU@ are substituted
#   est/<model>.est         parameter search ranges (m3 and m4 share m3_m4_migration.est)
#   sfs/<pair>_jointMAFpop1_0.obs   observed joint SFS from easySFS
#
# fastsimcoal2 expects the .obs file to share the prefix of the .tpl file, so the
# SFS is copied once per model. Output is byte-identical to the files that were run.
set -euo pipefail

if [ $# -ne 2 ]; then
    sed -n '5,12p' "$0"
    exit 1
fi

MU=$1
OUT=$2
FSC_DIR=$(cd "$(dirname "$0")/.." && pwd)
MODELS="m1_strict_isolation m2_constant_migration m3_ancient_migration m4_recent_migration"

mkdir -p "$OUT"
tail -n +2 "$FSC_DIR/pairs.tsv" | while IFS=$'\t' read -r PAIR N0 N1; do
    for MODEL in $MODELS; do
        PREFIX="$OUT/${PAIR}_${MODEL}"
        sed -e "s/^@N0@$/${N0}/" -e "s/^@N1@$/${N1}/" -e "s/@MU@/${MU}/" \
            "$FSC_DIR/templates/${MODEL}.tpl" > "${PREFIX}.tpl"
        case $MODEL in
            m3_*|m4_*) EST=m3_m4_migration ;;
            *)         EST=$MODEL ;;
        esac
        cp "$FSC_DIR/est/${EST}.est" "${PREFIX}.est"
        cp "$FSC_DIR/sfs/${PAIR}_jointMAFpop1_0.obs" "${PREFIX}_jointMAFpop1_0.obs"
    done
done
