#!/usr/bin/env bash
# Espera a A y B; luego reportes y mapas (A: M0 libre, B: M0 CSN).
cd "$(dirname "$0")"
PY=/home/alex/.conda/envs/kdellipspy/bin/python
C=isc_tomo/b0.05-0.30_noM14L
while ps -p $(paste -sd, na_AB.pids) >/dev/null 2>&1; do sleep 60; done
date
for f in na_A.log na_B.log; do echo "== $f"; grep -E "best misfit|Traceback|Error|M0 impuesto" $f | tail -3; done
OA=na_output_isc_tomo__b0.05-0.30_noM14L; OB=${OA}_M0
DYN_CASE=$C DYN_WORK=tsn_work_rep_A $PY report_dynamic.py $OA 2>&1 | grep -vE "Moment|^\s*$|convmPy|openMp|Reading|^ +[0-9]" &
DYN_CASE=$C DYN_WORK=tsn_work_rep_B DYN_M0=1.0e19 $PY report_dynamic.py $OB 2>&1 | grep -vE "Moment|^\s*$|convmPy|openMp|Reading|^ +[0-9]" &
wait
DYN_CASE=$C DYN_WORK=tsn_work_rep_A $PY map_dynamic.py $OA 2>&1 | grep -E "^ok|Error|Trace"
DYN_CASE=$C DYN_WORK=tsn_work_rep_B DYN_M0=1.0e19 $PY map_dynamic.py $OB 2>&1 | grep -E "^ok|Error|Trace"
