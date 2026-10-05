#!/usr/bin/env bash
# NA dinámico A (M0 libre) y B (M0 de agencia) para los dos planos nodales,
# con varias semillas. Cada corrida lleva su propio DYN_WORK (fd3d escribe en
# disco: dos corridas compartiendo directorio se pisan) y su propio log.
#
# El preflight de na_dynamic.py corre antes del NA en cada proceso y aborta si
# ningún modelo de referencia ajusta (DYN_PREFLIGHT_MAX): mejor perder minutos
# que una noche.
#
# Uso:  ./run_AB.sh [semillas...]      (por defecto 0 1 2)
set -euo pipefail
cd "$(dirname "$0")"
PY=/home/alex/.conda/envs/kdellipspy/bin/python
M0_USGS=2.29e19          # USGS Mww 6.8 (us6000a4yi, tensor W-phase)
# Casos dinámicos: derivados del cinemático que ganó el paso 1 (hipocentro
# ISOLA, modelo de Potin, banda 0.04-0.15 Hz, sin A09F ni la transversal de
# A08F), con 20x20 subfallas. Ver make_dyn_case.py.
SUF=_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F
SEEDS=("${@:-0 1 2}")
read -ra SEEDS <<< "${SEEDS[*]}"
MAXJOBS=${MAXJOBS:-6}    # 8 núcleos, cada corrida es n_jobs=1; dejar aire para el resto
rm -f na_AB.pids

for plane in "np1${SUF}" "np2${SUF}"; do
  for seed in "${SEEDS[@]}"; do
    for mode in A B; do
      short="${plane%%_*}"              # np1 / np2, para nombres de log y work
      tag="${short}_${mode}_s${seed}"
      work="tsn_work_${tag}"
      [[ $mode == B ]] && m0env="DYN_M0=$M0_USGS" || m0env="DYN_M0="
      echo "=== lanzando $tag (work $work)"
      env DYN_CASE="dyn_cases/${plane}" DYN_WORK="$work" NA_SEED="$seed" $m0env \
        nohup "$PY" -u na_dynamic.py > "na_${tag}.log" 2>&1 &
      echo "$!" >> na_AB.pids
      sleep 5   # el arranque de axitra compite por disco si salen todos a la vez
      while (( $(jobs -rp | wc -l) >= MAXJOBS )); do wait -n; done
    done
  done
done
wait
