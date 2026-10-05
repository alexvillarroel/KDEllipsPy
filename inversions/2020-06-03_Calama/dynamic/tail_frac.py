"""Fraccion del momento liberada despues de T, para varios modelos del prior.

Una sola corrida de fd3d por modelo (NT=2759, 40 s), leida DESPUES de que el
proceso termina, y borrando los .res antes de cada corrida. Truncar a T
equivale a descartar esa cola, asi que esto responde la pregunta sin comparar
dos corridas (que fue donde se me enredo el test anterior).
"""
import os
import subprocess
import sys

import numpy as np

os.environ.setdefault("DYN_CASE", "dyn_cases/np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F")
sys.path.insert(0, "/home/alex/KDEllipsPy/inversions/2020-06-03_Calama/dynamic")
import kdellipspy as kde  # noqa: E402
import forward_dynamic as F  # noqa: E402
from kdellipspy.core.geometry import build_tsn_dynamic_fields, tsn_hypocentre_coarse  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import (  # noqa: E402
    read_tsn_fault_field, write_fd3d_tsn_forwardmodel,
)

cfg = kde.ConfigParser(str(F.CASE_DIR / "input.ctl"))
fp = cfg.fault_plane
hypo = tsn_hypocentre_coarse(fp)
WORK = F.HERE / "tsn_work_tail"
F.WORK = WORK
F.NT = 2759
F.setup_work_dir(cfg)

MODELS = {
    "lento  (a=b=8, Te=2,  cte1=1.8,  Dc=4.0)": [8.0, 8.0, hypo[0], hypo[1], 0.0, 2.0, 1.8, 1.1, 1.5, 4.0],
    "rapido (a=b=8, Te=20, cte1=1.05, Dc=0.4)": [8.0, 8.0, hypo[0], hypo[1], 0.0, 20.0, 1.05, 1.1, 1.5, 0.4],
    "medio  (a=b=5, Te=8,  cte1=1.3,  Dc=1.5)": [5.0, 5.0, hypo[0], hypo[1], 0.0, 8.0, 1.3, 1.1, 1.5, 1.5],
    "chico  (a=b=3, Te=15, cte1=1.15, Dc=1.0)": [3.0, 3.0, hypo[0], hypo[1], 0.0, 15.0, 1.15, 1.1, 1.5, 1.0],
}

print(f"{'modelo':<42}{'t95':>7}{'t99':>7}{'cola>15s':>10}{'cola>20s':>10}{'cola>25s':>10}", flush=True)
for label, mv in MODELS.items():
    m = np.array(mv, dtype=np.float32)
    t0, ts, dc = build_tsn_dynamic_fields(m, fp.nx, fp.ny, F.tsn_grid(cfg), hypo=hypo)
    write_fd3d_tsn_forwardmodel(WORK / "forwardmodel.dat", t0, ts, dc)
    for f in ("sliprateZ.res", "sliprateX.res"):
        (WORK / "result" / f).unlink(missing_ok=True)
    subprocess.run([str(WORK / F.BINARY)], cwd=WORK, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    srz = read_tsn_fault_field(WORK / "result" / "sliprateZ.res", F.NXTT, F.NZTT)
    mr = np.abs(srz).sum(axis=(1, 2))
    t = np.arange(mr.size) * F.DT
    c = np.cumsum(mr) / mr.sum()
    t95, t99 = t[np.searchsorted(c, 0.95)], t[np.searchsorted(c, 0.99)]
    tails = [1.0 - c[np.searchsorted(t, T)] for T in (15.0, 20.0, 25.0)]
    print(f"{label:<42}{t95:7.2f}{t99:7.2f}" + "".join(f"{100 * x:9.3f}%" for x in tails), flush=True)
