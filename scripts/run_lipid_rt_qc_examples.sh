#!/usr/bin/env bash
# Reproduce results/lipid_rt_qc: QC reports + error-injection benchmark on public
# MetaboLights reverse-phase lipidomics datasets (MS-DIAL annotations).
set -euo pipefail
cd "$(dirname "$0")/.."
RAW=data/raw/rt
BASE=https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public
mkdir -p "$RAW" results/lipid_rt_qc
while read -r st maf; do
  f="$RAW/$maf"
  [ -s "$f" ] || curl -sfL -o "$f" "$BASE/$st/$maf"
  lipid-rt-qc run "$f" --out "results/lipid_rt_qc/$st"
  lipid-rt-qc benchmark "$f" --out "results/lipid_rt_qc/$st/benchmark"
done <<'LIST'
MTBLS5861 m_MTBLS5861_LC-MS_cid_positive_reverse-phase_metabolite_profiling_v2_maf.tsv
MTBLS14841 m_MTBLS14841_LC-MS_positive_reverse-phase_v2_maf.tsv
MTBLS14692 m_MTBLS14692_LC-MS_alternating_reverse-phase_v2_maf.tsv
LIST
