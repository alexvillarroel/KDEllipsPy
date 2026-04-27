from dataclasses import dataclass, field
from math import cos, log10, pi as math_pi, radians, sin, sqrt
from typing import List, Optional, Tuple

import numpy as np
from pyproj import CRS, Transformer

# Handle both relative imports (when used as package) and direct imports (from notebooks)
try:
    from .config_parser import ConfigParser
except ImportError:
    from config_parser import ConfigParser


class UTMProjection:
    """
    UTM (Universal Transverse Mercator) projection using pyproj.
    
    Automatically detects the UTM zone based on the center coordinates.
    Uses WGS84 ellipsoid for accurate geodetic transformations.
    This replaces the legacy Lambert conformal projection with a
    faster (C-native) and more standard approach.
    """

    def __init__(self, lat0_deg: float, lon0_deg: float):
        """Initialize UTM projection with automatic zone detection."""
        # Validate coordinates
        if abs(lon0_deg) > 180 or abs(lat0_deg) > 90:
            raise ValueError(f"Coordinates out of range: lat={lat0_deg}, lon={lon0_deg}")
        
        self.lat0_deg = float(lat0_deg)
        self.lon0_deg = float(lon0_deg)
        
        # Auto-detect UTM zone from central longitude
        zone = self._calculate_utm_zone(lon0_deg)
        hemisphere = "north" if lat0_deg >= 0 else "south"
        
        self.crs_wgs84 = CRS.from_epsg(4326)  # WGS84 (lat/lon)
        self.crs_utm = CRS.from_proj4(f"+proj=utm +zone={zone} +{hemisphere}")
        
        # Create bidirectional transformers
        self.to_utm = Transformer.from_crs(self.crs_wgs84, self.crs_utm, always_xy=True)
        self.to_wgs84 = Transformer.from_crs(self.crs_utm, self.crs_wgs84, always_xy=True)

    @staticmethod
    def _calculate_utm_zone(lon_deg: float) -> int:
        """Calculate UTM zone from longitude (1-60)."""
        zone = int((lon_deg + 180) / 6) + 1
        return max(1, min(60, zone))

    def latlon_to_xy(self, lat_deg: float, lon_deg: float) -> Tuple[float, float]:
        """Convert lat/lon (degrees) to UTM x/y (meters)."""
        if abs(lon_deg) > 180 or abs(lat_deg) > 90:
            raise ValueError(f"Coordinates out of range: lat={lat_deg}, lon={lon_deg}")
        
        x_m, y_m = self.to_utm.transform(lon_deg, lat_deg)
        return float(x_m), float(y_m)

    def xy_to_latlon(self, x_m: float, y_m: float) -> Tuple[float, float]:
        """Convert UTM x/y (meters) to lat/lon (degrees)."""
        lon_deg, lat_deg = self.to_wgs84.transform(x_m, y_m)
        return float(lat_deg), float(lon_deg)


@dataclass
class Station:
    index: int
    name: str
    x_m: float
    y_m: float
    z_m: float
    lat: float
    lon: float


@dataclass
class StationGeometry:
    ref_lat: float
    ref_lon: float
    stations: List[Station]
    _projection: UTMProjection = field(default=None, init=False, repr=False)

    def __post_init__(self):
        self._projection = UTMProjection(self.ref_lat, self.ref_lon)

    @property
    def nstations(self) -> int:
        return len(self.stations)

    def _latlon_to_local(self, lat: float, lon: float) -> Tuple[float, float]:
        x_m, y_m = self._projection.latlon_to_xy(lat, lon)
        return x_m, y_m

    def _local_to_latlon(self, x_m: float, y_m: float) -> Tuple[float, float]:
        return self._projection.xy_to_latlon(x_m, y_m)

    def to_axitra_stations(self, latlon: bool = False) -> np.ndarray:
        rows = []
        for st in self.stations:
            if latlon:
                rows.append([st.index, st.lat, st.lon, st.z_m])
            else:
                rows.append([st.index, st.x_m, st.y_m, st.z_m])
        return np.array(rows, dtype="float64")


@dataclass
class Subfault:
    index: int
    x_m: float
    y_m: float
    z_m: float
    rupture_time_s: float
    mu_pa: float = 0.0
    area_m2: float = 0.0


@dataclass
class SourcePoint:
    index: int
    subfault_index: int
    x_m: float
    y_m: float
    z_m: float
    rupture_time_s: float
    moment: float = 0.0
    displacement: float = 0.0
    strike_deg: float = 0.0
    dip_deg: float = 0.0
    rake_deg: float = 0.0
    width: float = 0.0
    length: float = 0.0
    mu_pa: float = 0.0
    basis_slot: int = 0


@dataclass
class FaultGeometry:
    length_strike_m: float
    length_dip_m: float
    hypo_strike_m: float
    hypo_dip_m: float
    nx: int
    ny: int
    strike_deg: float
    dip_deg: float
    rake_deg: float
    source_depth_m: float
    source_lat: float
    source_lon: float
    rupture_velocity_km_s: float
    mt_enabled: bool
    subfaults: List[Subfault]
    source_points: List[SourcePoint]
    _projection: UTMProjection = field(default=None, init=False, repr=False)

    def __post_init__(self):
        self._projection = UTMProjection(self.source_lat, self.source_lon)

    @property
    def nsubfaults(self) -> int:
        return len(self.subfaults)

    @property
    def nsources(self) -> int:
        return len(self.source_points)

    def _local_to_latlon(self, x_m: float, y_m: float) -> Tuple[float, float]:
        return self._projection.xy_to_latlon(x_m, y_m)

    def _latlon_to_local(self, lat: float, lon: float) -> Tuple[float, float]:
        return self._projection.latlon_to_xy(lat, lon)

    def to_axitra_sources(self, latlon: bool = True) -> np.ndarray:
        rows = []
        for sp in self.source_points:
            if latlon:
                lat, lon = self._local_to_latlon(sp.x_m, sp.y_m)
                rows.append([sp.index, lat, lon, sp.z_m])
            else:
                rows.append([sp.index, sp.x_m, sp.y_m, sp.z_m])
        return np.array(rows, dtype="float64")

    def to_axitra_hist(self) -> np.ndarray:
        rows = []
        for sp in self.source_points:
            if self.mt_enabled:
                rows.append(
                    [
                        sp.index,
                        sp.moment,
                        sp.strike_deg,
                        sp.dip_deg,
                        sp.rake_deg,
                        sp.width,
                        sp.length,
                        sp.rupture_time_s,
                    ]
                )
            else:
                rows.append(
                    [
                        sp.index,
                        sp.displacement,
                        sp.strike_deg,
                        sp.dip_deg,
                        sp.rake_deg,
                        sp.width,
                        sp.length,
                        sp.rupture_time_s,
                    ]
                )
        return np.array(rows, dtype="float64")

    def total_moment_nm(self) -> float:
        """Return total scalar seismic moment in N.m for current source set."""
        if self.mt_enabled:
            return float(sum(abs(float(sp.moment)) for sp in self.source_points))

        sf_by_index = {int(sf.index): sf for sf in self.subfaults}
        total = 0.0
        for sp in self.source_points:
            sf = sf_by_index.get(int(sp.subfault_index))
            if sf is None:
                continue
            total += abs(float(sf.mu_pa) * float(sf.area_m2) * float(sp.displacement))
        return float(total)

    def moment_magnitude_mw(self) -> float:
        """Return moment magnitude Mw from total moment (Hanks & Kanamori, M0 in N.m)."""
        m0 = self.total_moment_nm()
        if m0 <= 0.0:
            return float("-inf")
        return float((2.0 / 3.0) * (log10(m0) - 9.1))


class GeometryBuilder:
    """Build fault geometry from input.ctl through ConfigParser."""

    def __init__(self, config: ConfigParser):
        self.config = config

    @classmethod
    def from_input_ctl(cls, input_ctl_path: str) -> "GeometryBuilder":
        return cls(ConfigParser(input_ctl_path))

    def _mt_basis_and_amplitudes(self, nsubfaults: int) -> Tuple[List[float], List[float], List[float], List[float], List[float], List[float]]:
        mt = self.config.moment_tensor

        scale = 10.0 ** float(mt.exponent)
        mrr = float(mt.mrr) * scale
        mtt = float(mt.mtt) * scale
        mpp = float(mt.mpp) * scale
        mrt = float(mt.mrt) * scale
        mrp = float(mt.mrp) * scale
        mtp = float(mt.mtp) * scale

        c1 = mrt
        c2 = mrp
        c3 = -mtp
        c5 = mpp - mrr
        c6 = mtt + mpp - mrr
        c4 = mtt + mpp - 2.0 * mrr

        vals = [
            abs(c1),
            abs(c2),
            abs(c3),
            abs(c4),
            abs(c5),
            abs(c6) * sqrt(1.5),
        ]
        vals = [v / float(nsubfaults) for v in vals]

        b_str = [0.0, 270.0, 0.0, 90.0, 0.0, 0.0]
        b_dip = [90.0, 90.0, 90.0, 45.0, 45.0, 0.0]
        b_rak = [0.0, -90.0, 90.0, 90.0, 90.0, 0.0]
        b_wd = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
        b_len = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        return vals, b_str, b_dip, b_rak, b_wd, b_len

    def build(
        self,
        slip_geom: float = 1.0,
    ) -> FaultGeometry:
        fp = self.config.fault_plane
        src = self.config.source_position
        mt = self.config.moment_tensor

        phi = radians(src.strike)
        delta = radians(src.dip)

        dstk = fp.lx / float(fp.nx)
        ddip = fp.ly / float(fp.ny)
        area_subfault_m2 = float(dstk * ddip)

        dxs = dstk * cos(phi)
        dxd = -ddip * cos(delta) * sin(phi)
        dys = dstk * sin(phi)
        dyd = ddip * cos(delta) * cos(phi)
        dzd = ddip * sin(delta)

        depth_m = float(src.depth) * 1000.0

        layers = self.config.velocity_model.layers

        # Support both conventions documented by AXITRA:
        # - layer thickness per row
        # - upper-interface depth per row (first row must be 0)
        is_interface_depth = False
        if len(layers) > 1:
            vals = [float(l.thickness) for l in layers]
            if abs(vals[0]) < 1e-9 and all(vals[i + 1] >= vals[i] for i in range(len(vals) - 1)):
                is_interface_depth = True

        def get_rigidity_at_depth(z_m: float) -> float:
            z_local = max(float(z_m), 0.0)

            if is_interface_depth:
                current_layer = layers[-1]
                for i, layer in enumerate(layers):
                    z_top = float(layer.thickness)
                    z_bot = float(layers[i + 1].thickness) if i + 1 < len(layers) else np.inf
                    if z_local >= z_top and z_local < z_bot:
                        current_layer = layer
                        break
            else:
                cum_depth = 0.0
                current_layer = layers[-1]
                for layer in layers:
                    cum_depth += float(layer.thickness)
                    current_layer = layer
                    if z_local <= cum_depth:
                        break

            vs = float(current_layer.vs)
            rho = float(current_layer.rho)
            # If vs is in km/s (small values), convert to m/s.
            if abs(vs) < 100.0:
                vs = vs * 1000.0
            # If rho is in kg/km^3 (very large values), convert to kg/m^3.
            if abs(rho) > 1e8:
                rho = rho / 1e9
            return float(rho * (vs ** 2))

        x0 = (dxs + dxd) / 2.0 - fp.hx * cos(phi) + fp.hy * cos(delta) * sin(phi)
        y0 = (dys + dyd) / 2.0 - fp.hx * sin(phi) - fp.hy * cos(delta) * cos(phi)
        z0 = dzd / 2.0 - fp.hy * sin(delta) + depth_m

        subfaults: List[Subfault] = []
        for idip in range(1, fp.ny + 1):
            for istk in range(1, fp.nx + 1):
                idx = (idip - 1) * fp.nx + istk
                x = x0 + (istk - 1) * dxs + (idip - 1) * dxd
                y = y0 + (istk - 1) * dys + (idip - 1) * dyd
                z = z0 + (idip - 1) * dzd
                mu_pa = get_rigidity_at_depth(z)
                subfaults.append(
                    Subfault(
                        index=idx,
                        x_m=x,
                        y_m=y,
                        z_m=z,
                        rupture_time_s=0.0,
                        mu_pa=mu_pa,
                        area_m2=area_subfault_m2,
                    )
                )

        source_points: List[SourcePoint] = []
        mt_enabled = int(mt.flag) == 1

        if mt_enabled:
            _, b_str, b_dip, b_rak, b_wd, b_len = self._mt_basis_and_amplitudes(len(subfaults))
            sid = 1
            for sf in subfaults:
                for k in range(6):
                    source_points.append(
                        SourcePoint(
                            index=sid,
                            subfault_index=sf.index,
                            x_m=sf.x_m,
                            y_m=sf.y_m,
                            z_m=sf.z_m,
                            rupture_time_s=0.0,
                            moment=0.0,
                            displacement=0.0,
                            strike_deg=b_str[k],
                            dip_deg=b_dip[k],
                            rake_deg=b_rak[k],
                            width=b_wd[k],
                            length=b_len[k],
                            mu_pa=sf.mu_pa,
                            basis_slot=k,
                        )
                    )
                    sid += 1
        else:
            sid = 1
            for sf in subfaults:
                source_points.append(
                    SourcePoint(
                        index=sid,
                        subfault_index=sf.index,
                        x_m=sf.x_m,
                        y_m=sf.y_m,
                        z_m=sf.z_m,
                        rupture_time_s=0.0,
                        moment=0.0,
                        displacement=0.0,
                        strike_deg=float(src.strike),
                        dip_deg=float(src.dip),
                        rake_deg=float(src.rake),
                        width=float(ddip),
                        length=float(dstk),
                        mu_pa=sf.mu_pa,
                        basis_slot=0,
                    )
                )
                sid += 1

        return FaultGeometry(
            length_strike_m=fp.lx,
            length_dip_m=fp.ly,
            hypo_strike_m=fp.hx,
            hypo_dip_m=fp.hy,
            nx=fp.nx,
            ny=fp.ny,
            strike_deg=src.strike,
            dip_deg=src.dip,
            rake_deg=src.rake,
            source_depth_m=depth_m,
            source_lat=src.latitude,
            source_lon=src.longitude,
            rupture_velocity_km_s=0.0,
            mt_enabled=mt_enabled,
            subfaults=subfaults,
            source_points=source_points,
        )


class EllipticalSlipMapper:
    """Apply ellipse-based slip factors to fault geometry source points."""

    def __init__(self, config: ConfigParser):
        self.config = config

    def _prepare(self, model: np.ndarray) -> dict:
        if len(model) < 7:
            raise ValueError("Model must include 7 parameters: a1,a2,theta,np,tp,dmax,vr")

        nx = int(self.config.fault_plane.nx)
        ny = int(self.config.fault_plane.ny)
        dstk = float(self.config.fault_plane.lx) / float(nx)
        ddip = float(self.config.fault_plane.ly) / float(ny)

        a1_m = float(model[0]) * 1000.0
        a2_m = float(model[1]) * 1000.0
        alpha = float(model[2]) * math_pi
        np_frac_val = float(model[3])
        tp_angle = float(model[4]) * 2.0 * math_pi
        dmax_val = float(model[5])
        vr_km_s = float(model[6])
        vr_m_s = max(vr_km_s * 1000.0, 1.0)

        estk = float(self.config.fault_plane.hx)
        edip = float(self.config.fault_plane.hy)
        slip_shape = int(self.config.ellipse.slip_shape)

        x01 = a1_m * np_frac_val * cos(tp_angle)
        y01 = a2_m * np_frac_val * sin(tp_angle)
        xe = x01 * cos(alpha) + y01 * sin(alpha) + estk
        ye = -x01 * sin(alpha) + y01 * cos(alpha) + edip

        return {
            "nx": nx,
            "dstk": dstk,
            "ddip": ddip,
            "a1_m": a1_m,
            "a2_m": a2_m,
            "alpha": alpha,
            "dmax": dmax_val,
            "vr_km_s": vr_km_s,
            "vr_m_s": vr_m_s,
            "xe": xe,
            "ye": ye,
            "slip_shape": slip_shape,
        }

    @staticmethod
    def _subfault_fault_plane_xy(sf_idx: int, nx: int, dstk: float, ddip: float) -> Tuple[float, float]:
        istk = ((sf_idx - 1) % nx) + 1
        idip = ((sf_idx - 1) // nx) + 1
        xpos = (float(istk) - 0.5) * dstk
        ypos = (float(idip) - 0.5) * ddip
        return xpos, ypos

    def _slip_factor_for_subfault(self, sf_idx: int, prepared: dict) -> float:
        xpos, ypos = self._subfault_fault_plane_xy(
            sf_idx=sf_idx,
            nx=int(prepared["nx"]),
            dstk=float(prepared["dstk"]),
            ddip=float(prepared["ddip"]),
        )

        alpha = float(prepared["alpha"])
        xe = float(prepared["xe"])
        ye = float(prepared["ye"])
        a1_m = float(prepared["a1_m"])
        a2_m = float(prepared["a2_m"])

        xx = (xpos - xe) * cos(alpha) - (ypos - ye) * sin(alpha)
        yy = (xpos - xe) * sin(alpha) + (ypos - ye) * cos(alpha)

        if a1_m > 0 and a2_m > 0:
            d = (xx / a1_m) ** 2 + (yy / a2_m) ** 2
        else:
            d = np.inf

        if d > 1.0:
            return 0.0

        slip_shape = int(prepared["slip_shape"])
        if slip_shape == 0:
            return 1.0
        if slip_shape == 1:
            return float(np.exp(-d))
        return float(np.sqrt(max(0.0, 1.0 - d)))

    def _mt_component_weights(self) -> List[float]:
        mt = self.config.moment_tensor

        scale = 10.0 ** float(mt.exponent)
        mrr = float(mt.mrr) * scale
        mtt = float(mt.mtt) * scale
        mpp = float(mt.mpp) * scale
        mrt = float(mt.mrt) * scale
        mrp = float(mt.mrp) * scale
        mtp = float(mt.mtp) * scale

        c1 = mrt
        c2 = mrp
        c3 = -mtp
        c5 = mpp - mrr
        c6 = mtt + mpp - mrr
        c4 = mtt + mpp - 2.0 * mrr

        vals = [
            abs(c1),
            abs(c2),
            abs(c3),
            abs(c4),
            abs(c5),
            abs(c6) * sqrt(1.5),
        ]
        total = float(sum(vals))
        if total <= 0.0:
            return [1.0 / 6.0] * 6
        return [v / total for v in vals]

    def apply_to_geometry(
        self,
        geom: FaultGeometry,
        model: np.ndarray,
    ) -> FaultGeometry:
        prepared = self._prepare(model)
        dmax = float(prepared["dmax"])
        vr_m_s = float(prepared["vr_m_s"])
        source_depth_m = float(self.config.source_position.depth) * 1000.0

        sf_by_index = {int(sf.index): sf for sf in geom.subfaults}

        for sf in sf_by_index.values():
            tr = sqrt(
                sf.x_m * sf.x_m
                + sf.y_m * sf.y_m
                + (sf.z_m - source_depth_m) * (sf.z_m - source_depth_m)
            ) / vr_m_s
            sf.rupture_time_s = float(tr)

        slip_factor_by_subfault = {}
        for sf_idx in sf_by_index.keys():
            slip_factor_by_subfault[sf_idx] = self._slip_factor_for_subfault(sf_idx, prepared)

        mt_weights = self._mt_component_weights() if geom.mt_enabled else None
        for sp in geom.source_points:
            sf_idx = int(sp.subfault_index)
            sf = sf_by_index.get(sf_idx)
            if sf is None:
                continue

            sp.rupture_time_s = float(sf.rupture_time_s)
            slip_factor = slip_factor_by_subfault[sf_idx]
            slip_real = float(dmax * slip_factor)

            if geom.mt_enabled and mt_weights is not None:
                moment_total = float(sf.mu_pa * sf.area_m2 * slip_real)
                k = int(sp.basis_slot) % 6
                sp.moment = float(moment_total * mt_weights[k])
                sp.displacement = 0.0
            else:
                sp.moment = 0.0
                sp.displacement = slip_real

        # Filter source points with near-zero source terms that do not contribute to synthetics.
        if geom.mt_enabled:
            threshold_moment = 1e14
            geom.source_points = [sp for sp in geom.source_points if abs(sp.moment) > threshold_moment]
        else:
            threshold_disp = 1e-14
            geom.source_points = [sp for sp in geom.source_points if abs(sp.displacement) > threshold_disp]

        # AXITRA expects sequential source indices 1..N.
        for i, sp in enumerate(geom.source_points, start=1):
            sp.index = i

        geom.rupture_velocity_km_s = float(prepared["vr_km_s"])

        return geom


def build_geometry_from_input_ctl(
    input_ctl_path: str,
    slip_geom: float = 1.0,
) -> FaultGeometry:
    builder = GeometryBuilder.from_input_ctl(input_ctl_path)
    return builder.build(
        slip_geom=slip_geom,
    )


def build_station_geometry(
    ref_lat: float,
    ref_lon: float,
    station_data: List[Tuple[str, float, float, float]],
) -> StationGeometry:
    stations: List[Station] = []
    ref_geometry = StationGeometry(ref_lat=ref_lat, ref_lon=ref_lon, stations=[])

    # UTM absolute coordinates: x=east, y=north
    epi_east_utm, epi_north_utm = ref_geometry._latlon_to_local(ref_lat, ref_lon)

    for idx, (name, lat, lon, elev_m) in enumerate(station_data, start=1):
        # Get station in UTM coordinates (absolute)
        sta_east_utm, sta_north_utm = ref_geometry._latlon_to_local(lat, lon)

        # Convert to local convention used by forward/geometry:
        # x -> north, y -> east (relative to event location).
        x_local = sta_north_utm - epi_north_utm
        y_local = sta_east_utm - epi_east_utm

        stations.append(
            Station(
                index=idx,
                name=name,
                x_m=x_local,
                y_m=y_local,
                z_m=elev_m,
                lat=lat,
                lon=lon,
            )
        )

    return StationGeometry(ref_lat=ref_lat, ref_lon=ref_lon, stations=stations)


@dataclass
class EllipseDiagnosticsResult:
    x_m: np.ndarray
    y_m: np.ndarray
    d: np.ndarray
    slip_factor: np.ndarray
    inside_mask: np.ndarray
    xe_m: float
    ye_m: float
    a1_m: float
    a2_m: float
    alpha_rad: float
    active_subfault_ratio: float


class EllipseDiagnostics:
    """Evaluate and plot the ellipse/slip mapping that feeds the forward model."""

    def __init__(self, config: ConfigParser):
        self.config = config
        self.builder = GeometryBuilder(config)

    @classmethod
    def from_input_ctl(cls, input_ctl_path: str) -> "EllipseDiagnostics":
        return cls(ConfigParser(input_ctl_path))

    def evaluate(self, model: np.ndarray) -> EllipseDiagnosticsResult:
        if len(model) < 7:
            raise ValueError("Model must include 7 parameters: a1,a2,theta,np,tp,dmax,vr")

        geom = self.builder.build()

        a1_m = float(model[0]) * 1000.0
        a2_m = float(model[1]) * 1000.0
        alpha = float(model[2]) * math_pi
        np_frac = float(model[3])
        tp_angle = float(model[4]) * 2.0 * math_pi

        estk = float(self.config.fault_plane.hx)
        edip = float(self.config.fault_plane.hy)
        slip_shape = int(self.config.ellipse.slip_shape)

        x01 = a1_m * np_frac * cos(tp_angle)
        y01 = a2_m * np_frac * sin(tp_angle)
        xe = x01 * cos(alpha) + y01 * sin(alpha) + estk
        ye = -x01 * sin(alpha) + y01 * cos(alpha) + edip

        x_vals = []
        y_vals = []
        d_vals = []
        s_vals = []
        inside_vals = []

        nx = int(self.config.fault_plane.nx)
        dstk = float(self.config.fault_plane.lx) / float(nx)
        ddip = float(self.config.fault_plane.ly) / float(self.config.fault_plane.ny)

        for sf in geom.subfaults:
            sf_idx = int(sf.index)
            istk = ((sf_idx - 1) % nx) + 1
            idip = ((sf_idx - 1) // nx) + 1
            xpos = (float(istk) - 0.5) * dstk
            ypos = (float(idip) - 0.5) * ddip
            xx = (xpos - xe) * cos(alpha) - (ypos - ye) * sin(alpha)
            yy = (xpos - xe) * sin(alpha) + (ypos - ye) * cos(alpha)

            if a1_m > 0 and a2_m > 0:
                d = (xx / a1_m) ** 2 + (yy / a2_m) ** 2
            else:
                d = np.inf

            if d <= 1.0:
                if slip_shape == 0:
                    slip_factor = 1.0
                elif slip_shape == 1:
                    slip_factor = float(np.exp(-d))
                else:
                    slip_factor = float(np.sqrt(max(0.0, 1.0 - d)))
                inside = True
            else:
                slip_factor = 0.0
                inside = False

            x_vals.append(xpos)
            y_vals.append(ypos)
            d_vals.append(float(d))
            s_vals.append(float(slip_factor))
            inside_vals.append(bool(inside))

        x_arr = np.asarray(x_vals, dtype=float)
        y_arr = np.asarray(y_vals, dtype=float)
        d_arr = np.asarray(d_vals, dtype=float)
        s_arr = np.asarray(s_vals, dtype=float)
        in_arr = np.asarray(inside_vals, dtype=bool)

        active_ratio = float(np.count_nonzero(s_arr > 1e-14)) / float(len(s_arr)) if len(s_arr) else 0.0

        return EllipseDiagnosticsResult(
            x_m=x_arr,
            y_m=y_arr,
            d=d_arr,
            slip_factor=s_arr,
            inside_mask=in_arr,
            xe_m=float(xe),
            ye_m=float(ye),
            a1_m=float(a1_m),
            a2_m=float(a2_m),
            alpha_rad=float(alpha),
            active_subfault_ratio=active_ratio,
        )

    def plot(
        self,
        result: EllipseDiagnosticsResult,
        title: str = "Ellipse diagnostics (fault plane)",
        save_path: Optional[str] = None,
        show: bool = True,
    ):
        try:
            import matplotlib.pyplot as plt
        except Exception as exc:
            raise ImportError("matplotlib is required for EllipseDiagnostics.plot") from exc

        fig, ax = plt.subplots(figsize=(8, 7))

        sc = ax.scatter(
            result.x_m / 1000.0,
            result.y_m / 1000.0,
            c=result.slip_factor,
            s=70,
            cmap="inferno_r",
            edgecolors="k",
            linewidths=0.25,
        )

        t = np.linspace(0.0, 2.0 * math_pi, 361)
        ex = result.a1_m * np.cos(t)
        ey = result.a2_m * np.sin(t)
        bx = ex * np.cos(result.alpha_rad) + ey * np.sin(result.alpha_rad) + result.xe_m
        by = -ex * np.sin(result.alpha_rad) + ey * np.cos(result.alpha_rad) + result.ye_m
        ax.plot(bx / 1000.0, by / 1000.0, "c-", lw=2.0, label="Ellipse boundary")

        ax.plot(result.xe_m / 1000.0, result.ye_m / 1000.0, "c*", ms=12, label="Ellipse center")
        hx = float(self.config.fault_plane.hx) / 1000.0
        hy = float(self.config.fault_plane.hy) / 1000.0
        ax.plot(hx, hy, marker="D", color="royalblue", ms=10, label="Hypocenter")

        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("Slip factor")
        ax.set_xlabel("Along strike (km)")
        ax.set_ylabel("Along dip (km)")
        ax.set_title(f"{title}\nactive_subfault_ratio={result.active_subfault_ratio:.3f}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=180)

        if show:
            plt.show()
        return fig, ax
