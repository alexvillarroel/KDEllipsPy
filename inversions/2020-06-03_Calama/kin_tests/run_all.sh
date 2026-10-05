#!/usr/bin/env bash
# Lanza run_kin_multiseed.py en todos los casos (o los que se pasen), uno por
# proceso: cada corrida es n_jobs=1, así que el paralelismo va por caso.
set -euo pipefail
cd "$(dirname "$0")"
PY=/home/alex/.conda/envs/kdellipspy/bin/python
MAXJOBS=${MAXJOBS:-6}
mkdir -p logs
CASES=("$@")
(( ${#CASES[@]} )) || mapfile -t CASES < <(ls -d np[12]_*/ | sed 's#/##')
for c in "${CASES[@]}"; do
  echo "=== $c"
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 KIN_CASE="$c" \
    nohup "$PY" -u run_kin_multiseed.py > "logs/$c.log" 2>&1 &
  sleep 3
  while (( $(jobs -rp | wc -l) >= MAXJOBS )); do wait -n; done
done
wait
