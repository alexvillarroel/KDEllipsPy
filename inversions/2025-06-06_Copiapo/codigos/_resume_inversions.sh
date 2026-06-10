#!/usr/bin/env bash
# Supervisor desacoplado para terminar DC + MT aunque se cierre la terminal.
# Se lanza con:  setsid nohup bash codigos/_resume_inversions.sh > output/supervisor.log 2>&1 < /dev/null &
set -u
cd /home/alex/KDEllipsPy/inversions/2025-06-06_Copiapo
PY=/home/alex/.conda/envs/kdellipspy/bin/python
log(){ echo "[$(date +%F\ %H:%M:%S)] $*"; }

log "supervisor iniciado (PID $$)"

# 1) Esperar a que termine cualquier corrida DC en vuelo (la actual va al 79%).
while pgrep -f "run_inversion.py --input-ctl input_dc.ctl" >/dev/null 2>&1; do
    sleep 20
done

# 2) Asegurar que DC quedo completa; si la terminal la mato a medias, re-correr.
if ! grep -q "Inversion completada" output_dc/run.log 2>/dev/null; then
    log "DC incompleta -> re-corriendo desde cero"
    $PY codigos/run_inversion.py --input-ctl input_dc.ctl -o output_dc > output_dc/run.log 2>&1
fi
log "DC lista (misfit=$(grep -oE 'misfit minimo = [0-9.eE+-]+' output_dc/run.log | tail -1))"

# 3) MT (solo si no esta ya completa).
if ! grep -q "Inversion completada" output_mt/run.log 2>/dev/null; then
    log "corriendo MT"
    $PY codigos/run_inversion.py --input-ctl input.ctl -o output_mt > output_mt/run.log 2>&1
fi
log "MT lista (misfit=$(grep -oE 'misfit minimo = [0-9.eE+-]+' output_mt/run.log | tail -1))"

log "=== AMBAS COMPLETADAS ==="
