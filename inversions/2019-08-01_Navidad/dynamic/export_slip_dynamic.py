"""Archivo de slip por subfalla de la solución dinámica (mismo formato que
loc_mid_mt/output/slip_loc_mid_mt.txt) + grillas .grd en coordenadas de falla y geográficas.

Uso:  python export_slip_dynamic.py [A|B]      (B por defecto: M0 CSN, mejor ajuste)

Tiempo de ruptura: primer instante en que el slip-rate medio de la subfalla supera
1 mm/s (criterio SCEC). Rise time: tiempo entre el 10 % y el 90 % del slip final.
B se escala por autosemejanza a M0 = 1e19 N·m, como en todos los análisis.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
os.environ.setdefault("DYN_CASE", "isc_tomo/b0.05-0.30_noM14L")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path[:0] = [str(NAV / "okada"), str(HERE)]

import joblib  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

import okada_models as om  # noqa: E402
from kdellipspy.core.geometry import tsn_hypocentre_coarse  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import run_tsn_forward  # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else "B"
RT_THRESHOLD = 1e-3  # m/s


def main():
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    _, geom = om.subfaults(cfg)
    fp, sp = cfg.fault_plane, cfg.source_position
    r = joblib.load(om.DYN_OUT[WHICH] / "na_result.joblib")
    vals = np.asarray(r.best_model.model)
    full = om.nd.to_full_model_named(r.param_names, vals)
    rate = run_tsn_forward(full, fp.nx, fp.ny, om.fd.tsn_grid(cfg), om.fd.tsn_run_cfg(), hypo=tsn_hypocentre_coarse(fp))
    fx, fz, dt = om.fd.NXTT // fp.nx, om.fd.NZTT // fp.ny, om.fd.DT
    mu = np.array([sf.mu_pa for sf in geom.subfaults])

    def to_sub(rc):  # (t, nxtT, nztT), manteo desde abajo -> (t, nsub) en el orden de la malla (idip desde arriba, istk rápido)
        nt = rc.shape[0]
        s = rc.reshape(nt, fp.nx, fx, fp.ny, fz).mean(axis=(2, 4))[:, :, ::-1]
        return s.transpose(0, 2, 1).reshape(nt, -1)

    vx, vz = to_sub(rate["sliprateX"]), to_sub(rate["sliprateZ"])
    k = 1.0
    if WHICH == "B":
        m0_raw = mu.mean() * np.hypot(rate["sliprateX"].sum(0), rate["sliprateZ"].sum(0)).sum() * dt * om.fd.DH**2
        k = om.M0_B / m0_raw
    vx, vz = vx * k, vz * k
    ss, ds = vx.sum(0) * dt, vz.sum(0) * dt
    slip = np.hypot(ss, ds)
    v = np.hypot(vx, vz)
    t = np.arange(v.shape[0]) * dt
    on = v > RT_THRESHOLD
    rt = np.where(on.any(0), t[on.argmax(0)], np.nan)
    cum = np.cumsum(v, axis=0) * dt
    frac = cum / np.where(cum[-1] > 0, cum[-1], 1.0)
    rise = np.where(slip > 0.01, t[(frac >= 0.9).argmax(0)] - t[(frac >= 0.1).argmax(0)], np.nan)
    rake = np.where(slip > 0, np.degrees(np.arctan2(ds, ss)), np.nan)
    area = (fp.lx / fp.nx) * (fp.ly / fp.ny)
    m0 = float((mu * area * slip).sum())
    mw = (np.log10(m0) - 9.1) / 1.5

    dstk, ddip = fp.lx / fp.nx / 1e3, fp.ly / fp.ny / 1e3
    rows = []
    for i, sf in enumerate(geom.subfaults):
        istk, idip = i % fp.nx + 1, i // fp.nx + 1
        x = (istk - 0.5) * dstk - fp.hx / 1e3
        y = (idip - 0.5) * ddip - fp.hy / 1e3
        rows.append(f"{i + 1:4d} {istk:3d} {idip:3d} {x:8.3f} {y:8.3f} {sf.lat:10.5f} {sf.lon:10.5f} {abs(sf.z_m) / 1e3:8.3f} "
                    f"{slip[i]:8.4f} {ss[i]:8.4f} {ds[i]:8.4f} {np.nan_to_num(rake[i]):7.1f} "
                    f"{np.nan_to_num(rt[i]):7.2f} {np.nan_to_num(rise[i]):6.2f} {v[:, i].max():7.4f} {mu[i]:.4e}")
    dt0 = dict(zip(r.param_names, vals)).get("dt0 (s)", float("nan"))
    pv = {n.split(" ")[0]: v for n, v in zip(r.param_names, vals)}
    te = pv.get("Te", om.nd.TE_REF)  # B no invierte Te: usa TE_REF y se escala por k
    phys = {"Te_MPa": te * k, "peak_MPa": pv["cte1"] * te * k, "Dc_m": pv["Dc"] * k}
    tag = f"dyn{WHICH}"
    head = [
        f"# Distribucion de slip - Navidad 2019-08-01, solucion DINAMICA {WHICH} (KDEllipsPy + fd3d_TSN, NA)",
        f"# Caso: {os.environ['DYN_CASE']} | hipocentro: lat {sp.latitude} lon {sp.longitude} prof {sp.depth} km | "
        f"origen {cfg.source_position.origin_time if hasattr(cfg.source_position, 'origin_time') else ''} + dt0 = {dt0:+.2f} s",
        f"# Plano: strike {sp.strike} dip {sp.dip} | malla {fp.nx}x{fp.ny} subfallas de {dstk:g}x{ddip:g} km (fd3d: celdas de {om.fd.DH:g} m promediadas)",
        "# Modelo: " + "  ".join(f"{n.split(' ')[0]}={v:.4f}" for n, v in zip(r.param_names, vals)),
        f"# misfit = {r.best_model.misfit:.6f} | escalamiento por autosemejanza k = {k:.3f}" + (" (M0 fijado al CSN)" if WHICH == "B" else ""),
        "# Parametros fisicos (escalados): " + "  ".join(f"{n}={v:.3f}" for n, v in phys.items()),
        f"# M0 = {m0:.4e} N.m  Mw = {mw:.3f} | slip max = {slip.max():.3f} m | subfallas con slip > 1 cm: {(slip > 0.01).sum()} de {slip.size}",
        "# x_strike/y_dip: centro de la subfalla respecto al hipocentro sobre el plano (x positivo = en el sentido del rumbo, y positivo = pendiente abajo)",
        f"# rupture_time_s: slip-rate > {RT_THRESHOLD * 1e3:g} mm/s desde el inicio de la ruptura; rise_time_s: 10-90 % del slip; rake desde (ss, ds)",
        "# idx istk idip x_strike_km y_dip_km lat lon depth_km slip_m ss_m ds_m rake_deg rupture_time_s rise_time_s peak_sliprate_ms mu_pa",
    ]
    out = om.DYN_OUT[WHICH]
    txt = out / f"slip_{tag}.txt"
    txt.write_text("\n".join(head + rows) + "\n")

    xs = (np.arange(fp.nx) + 0.5) * dstk - fp.hx / 1e3
    ys = (np.arange(fp.ny) + 0.5) * ddip - fp.hy / 1e3
    xr.DataArray(slip.reshape(fp.ny, fp.nx), coords={"y": ys, "x": xs}, dims=("y", "x"), name="slip").to_netcdf(out / f"slip_{tag}_fault.grd")
    la = np.array([sf.lat for sf in geom.subfaults]); lo = np.array([sf.lon for sf in geom.subfaults])
    glon = np.arange(lo.min(), lo.max() + 1e-9, 0.005); glat = np.arange(la.min(), la.max() + 1e-9, 0.005)
    gz = griddata(np.c_[lo, la], slip, tuple(np.meshgrid(glon, glat)), method="linear")
    xr.DataArray(gz, coords={"lat": glat, "lon": glon}, dims=("lat", "lon"), name="slip").to_netcdf(out / f"slip_{tag}_geo.grd")
    print(f"{txt}\nMw {mw:.3f}  slip max {slip.max():.2f} m  k {k:.3f}  rupture time max {np.nanmax(np.where(slip > 0.01, rt, np.nan)):.1f} s")


if __name__ == "__main__":
    main()
