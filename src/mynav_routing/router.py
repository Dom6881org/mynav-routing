"""Façade de routage météo pour MyNavApp.

route() calcule la trajectoire la plus rapide par méthode des isochrones :
à chaque pas, le front des positions atteignables est propagé sur tous les caps
candidats (vent local -> polaire -> vitesse), élagué, puis on trace la route
optimale depuis le point le plus proche de l'arrivée.

Contrat MyNavApp : départ = position NMEA courante, vent = GRIB téléchargé
à terre. Exécution 100 % offline.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .engine.isochrones import haversine_nm, move_nm
from .polar import Polar


@dataclass
class RouteResult:
    found: bool
    coords: list          # [(lat, lon), ...] du départ à l'arrivée
    times_s: list         # temps cumulé (s) à chaque point
    sog_kts: list         # vitesse surface (kt) sur chaque tronçon
    headings_deg: list    # cap suivi sur chaque tronçon
    eta_s: float = 0.0
    min_dist_nm: float = float("inf")


def route(
    start,
    finish,
    polar: Polar,
    wind,
    dt=1800.0,
    horizon_h=72.0,
    heading_step=5.0,
    bin_deg=0.25,
    arrival_nm=5.0,
    tws_max=None,
    land=None,
):
    """Calcule la route la plus rapide de start vers finish.

    start, finish : (lat, lon) en degrés décimaux.
    polar : Polar (.pol) ; wind : WindField (voir grib.py).
    dt : pas des isochrones (s) ; horizon_h : durée max (h) ;
    heading_step : pas des caps candidats (degrés) ;
    bin_deg : maille d'élagage (degrés) ; arrival_nm : rayon d'arrivée (mn) ;
    tws_max : contrainte de confort (kt), None = illimité ;
    land : LandMask optionnel ; les tronçons qui touchent la terre sont écartés.
    """
    fronts = [
        {
            "lat": np.array([float(start[0])]),
            "lon": np.array([float(start[1])]),
            "hdg": np.array([0.0]),
            "bs": np.array([0.0]),
            "parent": np.array([-1]),
        }
    ]
    headings = np.arange(0.0, 360.0, heading_step)
    n_steps = int(horizon_h * 3600.0 / dt)
    for k in range(n_steps):
        prev = fronts[-1]
        n = len(prev["lat"])
        if n == 0:
            break
        lat_r = np.repeat(prev["lat"], len(headings))
        lon_r = np.repeat(prev["lon"], len(headings))
        hdg = np.tile(headings, n)
        parent = np.repeat(np.arange(n), len(headings))
        u, v = wind.sample(lat_r, lon_r, k * dt)
        u = np.asarray(u)
        v = np.asarray(v)
        wind_from = np.degrees(np.arctan2(-u, -v)) % 360.0
        tws = np.hypot(u, v)
        twa = (hdg - wind_from) % 360.0
        bs = polar.speed_grid(tws, twa)
        keep = bs > 0.1
        if tws_max is not None:
            keep &= tws <= tws_max
        if not keep.any():
            break
        lat_r, lon_r, hdg, bs, parent = (
            lat_r[keep],
            lon_r[keep],
            hdg[keep],
            bs[keep],
            parent[keep],
        )
        lat2, lon2 = move_nm(lat_r, lon_r, hdg, bs * dt / 3600.0)
        if land is not None:
            ok = ~land.crosses_land(lat_r, lon_r, lat2, lon2)
            lat2, lon2, hdg, bs, parent = lat2[ok], lon2[ok], hdg[ok], bs[ok], parent[ok]
            if len(lat2) == 0:
                break
        # Élagage : un seul point par cellule de bin_deg degrés
        key = np.round(lat2 / bin_deg) * 4096.0 + np.round(lon2 / bin_deg)
        _, idx = np.unique(key, return_index=True)
        lat2, lon2, hdg, bs, parent = (
            lat2[idx],
            lon2[idx],
            hdg[idx],
            bs[idx],
            parent[idx],
        )
        d = haversine_nm(lat2, lon2, float(finish[0]), float(finish[1]))
        reached = d <= np.maximum(arrival_nm, bs * dt / 3600.0)
        if land is not None and reached.any():
            reached &= ~land.crosses_land(lat2, lon2, float(finish[0]), float(finish[1]))
        fronts.append(
            {"lat": lat2, "lon": lon2, "hdg": hdg, "bs": bs, "parent": parent}
        )
        if reached.any():
            i = int(np.argmin(np.where(reached, d, np.inf)))
            coords, bs_list, hdg_list = _traceback(fronts, i)
            if len(bs_list) > 1 and bs_list[0] == 0.0:
                bs_list[0] = bs_list[1]
            times = [x * dt for x in range(len(coords))]
            bs_final = max(bs_list[-1], 0.1)
            eta = times[-1] + float(d[i]) / bs_final
            coords.append((float(finish[0]), float(finish[1])))
            bs_list.append(bs_final)
            hdg_list.append(hdg_list[-1])
            times.append(eta)
            return RouteResult(
                found=True,
                coords=coords,
                times_s=times,
                sog_kts=bs_list,
                headings_deg=hdg_list,
                eta_s=eta,
                min_dist_nm=float(d[i]),
            )
    # Horizon dépassé : retourner la meilleure progression atteinte
    best_d, best_step, best_i = None, 0, 0
    for step, f in enumerate(fronts):
        if len(f["lat"]):
            dd = haversine_nm(f["lat"], f["lon"], float(finish[0]), float(finish[1]))
            m = float(np.min(dd))
            if best_d is None or m < best_d:
                best_d, best_step, best_i = m, step, int(np.argmin(dd))
    coords, bs_list, hdg_list = _traceback(fronts[: best_step + 1], best_i)
    if len(bs_list) > 1 and bs_list[0] == 0.0:
        bs_list[0] = bs_list[1]
    times = [x * dt for x in range(len(coords))]
    return RouteResult(
        found=False,
        coords=coords,
        times_s=times,
        sog_kts=bs_list,
        headings_deg=hdg_list,
        eta_s=0.0,
        min_dist_nm=best_d if best_d is not None else float("inf"),
    )


def _traceback(fronts, i):
    coords, bs_list, hdg_list = [], [], []
    step, j = len(fronts) - 1, int(i)
    while step >= 0:
        f = fronts[step]
        coords.append((float(f["lat"][j]), float(f["lon"][j])))
        bs_list.append(float(f["bs"][j]))
        hdg_list.append(float(f["hdg"][j]))
        j = int(f["parent"][j])
        step -= 1
    coords.reverse()
    bs_list.reverse()
    hdg_list.reverse()
    return coords, bs_list, hdg_list


def to_geojson(res: RouteResult):
    """Export GeoJSON (LineString [lon, lat]) pour la carte Leaflet de MyNavApp."""
    if not res.found:
        raise ValueError("Aucune route trouvée (horizon dépassé)")
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[lon, lat] for lat, lon in res.coords],
                },
                "properties": {
                    "eta_s": res.eta_s,
                    "sog_kts": [round(s, 2) for s in res.sog_kts],
                },
            }
        ],
    }
