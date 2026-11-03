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
from time import perf_counter

import numpy as np

from .bathy import DepthMask
from .engine.isochrones import haversine_nm, move_nm
from .land import LandMask
from .polar import Polar

_AUTO_LAND = object()

# B1 voile/moteur : vitesse moteur constante (kt) prise en relais dès que la
# vitesse voile issue de la polaire passe sous le seuil engine_below_kts.
ENGINE_SPEED_KTS = 4.0


@dataclass
class RouteResult:
    found: bool
    coords: list          # [(lat, lon), ...] du départ à l'arrivée
    times_s: list         # temps cumulé (s) à chaque point
    sog_kts: list         # vitesse surface (kt) d'arrivée à chaque point
    headings_deg: list    # cap suivi pour arriver à chaque point
    eta_s: float = 0.0
    min_dist_nm: float = float("inf")
    corridor_margin: float | None = None  # marge corridor effective du résultat
    # B1 : une entrée bool par point (même convention que sog_kts) — tronçon
    # d'arrivée au point passé en moteur ; None = pur voile (back-compat).
    motor_flags: list | None = None
    # Départ hors zone du champ GRIB (issue MyNavApp #141) : False = le champ
    # nearest-clippe le départ (vent du bord de grille, sans erreur) ; True =
    # couvert ; None = vent sans emprise connue (ex. ConstantWind).
    start_in_grib: bool | None = None


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
    land=_AUTO_LAND,
    bathy: DepthMask | None = None,
    min_depth_m=None,
    phase_timings=None,
    corridor_margin=0.2,
    engine_below_kts=2.0,
):
    """Calcule la route la plus rapide de start vers finish.

    start, finish : (lat, lon) en degrés décimaux.
    polar : Polar (.pol) ; wind : WindField (voir grib.py).
    dt : pas des isochrones (s) ; horizon_h : durée max (h) ;
    heading_step : pas des caps candidats (degrés) ;
    bin_deg : maille d'élagage (degrés) ; arrival_nm : rayon d'arrivée (mn) ;
    tws_max : contrainte de confort (kt), None = illimité ;
    land : LandMask optionnel ; chargé depuis Natural Earth par défaut. None désactive le masque.
    bathy : DepthMask optionnel (GEBCO) pour éviter terre et hauts-fonds ;
        None = pas de contrainte de profondeur.
    min_depth_m : profondeur minimale acceptée (m) ; défaut = bathy.default_min_depth_m.
    phase_timings : dict optionnel, cumule les temps sample_wind_s, land_mask_s et
        pruning_s, ainsi que le nombre de candidats explorés.
    corridor_margin : marge relative du corridor start-finish (défaut 0.2 = +20 %) ;
        les candidats hors de l'ellipse de foyers start/finish, dont le grand axe
        vaut (1 + marge) x distance directe, sont élagués avant tout autre test.
        None désactive l'élagage par corridor. Fallback : si le calcul échoue
        (found=False) avec un corridor actif, route() retente automatiquement
        sans corridor et rend le meilleur des deux résultats — une route
        trouvable (ex. contournement de cap dont le détour excède la marge)
        ne doit jamais être refusée à cause de cette optimisation. Le résultat
        porte la marge effectivement utilisée (corridor_margin) ; en cas de
        re-run, phase_timings reflète le second calcul, celui du résultat.
    engine_below_kts : seuil de vitesse voile (kt, défaut 2.0) sous lequel le
        moteur prend le relais (B1) : la vitesse polaire du tronçon est
        remplacée par ENGINE_SPEED_KTS (4.0 kt), cap inchangé, et le tronçon
        est marqué moteur dans RouteResult.motor_flags (une entrée par point,
        même convention que sog_kts) ainsi que dans properties.motor du
        GeoJSON. None = pur voile, aucun relais (back-compat, motor_flags
        vaut alors None).
    """
    # Départ hors zone GRIB (#141) : détection par duck-typing — un WindField
    # qui connaît son emprise expose covers() (GRIBWindField) ; ConstantWind
    # n'en a pas et start_in_grib reste None.
    covers = getattr(wind, "covers", None)
    start_in_grib = None if covers is None else bool(covers(start[0], start[1]))
    result = _route_impl(
        start,
        finish,
        polar=polar,
        wind=wind,
        dt=dt,
        horizon_h=horizon_h,
        heading_step=heading_step,
        bin_deg=bin_deg,
        arrival_nm=arrival_nm,
        tws_max=tws_max,
        land=land,
        bathy=bathy,
        min_depth_m=min_depth_m,
        phase_timings=phase_timings,
        corridor_margin=corridor_margin,
        engine_below_kts=engine_below_kts,
    )
    result.corridor_margin = corridor_margin
    result.start_in_grib = start_in_grib
    if result.found or corridor_margin is None:
        return result
    # Fallback corridor : 2e essai sans corridor, meilleur des deux résultats.
    retry = _route_impl(
        start,
        finish,
        polar=polar,
        wind=wind,
        dt=dt,
        horizon_h=horizon_h,
        heading_step=heading_step,
        bin_deg=bin_deg,
        arrival_nm=arrival_nm,
        tws_max=tws_max,
        land=land,
        bathy=bathy,
        min_depth_m=min_depth_m,
        phase_timings=phase_timings,
        corridor_margin=None,
        engine_below_kts=engine_below_kts,
    )
    retry.corridor_margin = None
    retry.start_in_grib = start_in_grib
    if retry.found or retry.min_dist_nm <= result.min_dist_nm:
        return retry
    return result


def _route_impl(
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
    land=_AUTO_LAND,
    bathy: DepthMask | None = None,
    min_depth_m=None,
    phase_timings=None,
    corridor_margin=0.2,
    engine_below_kts=2.0,
):
    """Corps du calcul par isochrones (fallback corridor : route(), moteur B1)."""
    if phase_timings is not None:
        phase_timings.update(
            sample_wind_s=0.0,
            land_mask_s=0.0,
            pruning_s=0.0,
            explored_nodes=0,
        )
    if land is _AUTO_LAND:
        land = LandMask.from_open_data()
    if bathy is not None and min_depth_m is None:
        min_depth_m = bathy.default_min_depth_m
    flat, flon = float(finish[0]), float(finish[1])
    direct_nm = float(haversine_nm(float(start[0]), float(start[1]), flat, flon))
    fronts = [
        {
            "lat": np.array([float(start[0])]),
            "lon": np.array([float(start[1])]),
            "hdg": np.array([0.0]),
            "bs": np.array([0.0]),
            "mot": np.array([False]),
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
        phase_start = perf_counter() if phase_timings is not None else None
        u, v = wind.sample(lat_r, lon_r, k * dt)
        if phase_timings is not None:
            phase_timings["sample_wind_s"] += perf_counter() - phase_start
            phase_timings["explored_nodes"] += len(lat_r)
        u = np.asarray(u)
        v = np.asarray(v)
        wind_from = np.degrees(np.arctan2(-u, -v)) % 360.0
        tws = np.hypot(u, v)
        twa = (hdg - wind_from) % 360.0
        bs = polar.speed_grid(tws, twa)
        # B1 : relais moteur — sous le seuil de vitesse voile, le tronçon
        # passe en moteur à vitesse constante, cap inchangé.
        mot = np.zeros(len(bs), dtype=bool)
        if engine_below_kts is not None:
            mot = bs < engine_below_kts
            bs = np.where(mot, ENGINE_SPEED_KTS, bs)
        keep = bs > 0.1
        if tws_max is not None:
            keep &= tws <= tws_max
        if not keep.any():
            break
        lat_r, lon_r, hdg, bs, mot, parent = (
            lat_r[keep],
            lon_r[keep],
            hdg[keep],
            bs[keep],
            mot[keep],
            parent[keep],
        )
        lat2, lon2 = move_nm(lat_r, lon_r, hdg, bs * dt / 3600.0)
        if corridor_margin is not None:
            in_corridor = corridor_mask(lat2, lon2, start, finish, corridor_margin, direct_nm)
            if not in_corridor.all():
                lat2, lon2, hdg, bs, mot, parent = (
                    lat2[in_corridor],
                    lon2[in_corridor],
                    hdg[in_corridor],
                    bs[in_corridor],
                    mot[in_corridor],
                    parent[in_corridor],
                )
                lat_r, lon_r = lat_r[in_corridor], lon_r[in_corridor]
                if len(lat2) == 0:
                    break
        ok = None
        if land is not None:
            phase_start = perf_counter() if phase_timings is not None else None
            crosses_land = land.crosses_land(lat_r, lon_r, lat2, lon2)
            if phase_timings is not None:
                phase_timings["land_mask_s"] += perf_counter() - phase_start
            ok = ~crosses_land
        if bathy is not None:
            shallow = bathy.crosses_shallow(lat_r, lon_r, lat2, lon2, min_depth_m)
            ok = ~shallow if ok is None else (ok & ~shallow)
        if ok is not None:
            lat2, lon2, hdg, bs, mot, parent = (
                lat2[ok],
                lon2[ok],
                hdg[ok],
                bs[ok],
                mot[ok],
                parent[ok],
            )
            if len(lat2) == 0:
                break
        # Élagage : un seul point par cellule de bin_deg degrés, celui de temps
        # d'arrivée minimal (à temps égal, le premier rencontré : tri stable)
        phase_start = perf_counter() if phase_timings is not None else None
        d = haversine_nm(lat2, lon2, flat, flon)
        t = np.full(len(lat2), (k + 1) * dt)
        key = cell_keys(lat2, lon2, bin_deg)
        order = np.lexsort((t, key))
        ks = key[order]
        first = np.empty(len(ks), dtype=bool)
        first[0] = True
        first[1:] = ks[1:] != ks[:-1]
        idx = order[first]
        lat2, lon2, hdg, bs, mot, parent, d = (
            lat2[idx],
            lon2[idx],
            hdg[idx],
            bs[idx],
            mot[idx],
            parent[idx],
            d[idx],
        )
        if phase_timings is not None:
            phase_timings["pruning_s"] += perf_counter() - phase_start
        reached = d <= np.maximum(arrival_nm, bs * dt / 3600.0)
        if reached.any() and land is not None:
            phase_start = perf_counter() if phase_timings is not None else None
            reached &= ~land.crosses_land(lat2, lon2, float(finish[0]), float(finish[1]))
            if phase_timings is not None:
                phase_timings["land_mask_s"] += perf_counter() - phase_start
        if reached.any() and bathy is not None:
            reached &= ~bathy.crosses_shallow(
                lat2, lon2, float(finish[0]), float(finish[1]), min_depth_m
            )
        fronts.append(
            {"lat": lat2, "lon": lon2, "hdg": hdg, "bs": bs, "mot": mot, "parent": parent}
        )
        if reached.any():
            i = int(np.argmin(np.where(reached, d, np.inf)))
            coords, bs_list, hdg_list, mot_list = _traceback(fronts, i)
            if len(bs_list) > 1 and bs_list[0] == 0.0:
                bs_list[0] = bs_list[1]
            if len(mot_list) > 1:
                mot_list[0] = mot_list[1]
            times = [x * dt for x in range(len(coords))]
            bs_final = max(bs_list[-1], 0.1)
            eta = times[-1] + float(d[i]) / bs_final
            coords.append((float(finish[0]), float(finish[1])))
            bs_list.append(bs_final)
            hdg_list.append(hdg_list[-1])
            mot_list.append(mot_list[-1])
            times.append(eta)
            return RouteResult(
                found=True,
                coords=coords,
                times_s=times,
                sog_kts=bs_list,
                headings_deg=hdg_list,
                eta_s=eta,
                min_dist_nm=float(d[i]),
                motor_flags=_flags(mot_list, engine_below_kts),
            )
    # Horizon dépassé : retourner la meilleure progression atteinte
    best_d, best_step, best_i = None, 0, 0
    for step, f in enumerate(fronts):
        if len(f["lat"]):
            dd = haversine_nm(f["lat"], f["lon"], float(finish[0]), float(finish[1]))
            m = float(np.min(dd))
            if best_d is None or m < best_d:
                best_d, best_step, best_i = m, step, int(np.argmin(dd))
    coords, bs_list, hdg_list, mot_list = _traceback(fronts[: best_step + 1], best_i)
    if len(bs_list) > 1 and bs_list[0] == 0.0:
        bs_list[0] = bs_list[1]
    if len(mot_list) > 1:
        mot_list[0] = mot_list[1]
    times = [x * dt for x in range(len(coords))]
    return RouteResult(
        found=False,
        coords=coords,
        times_s=times,
        sog_kts=bs_list,
        headings_deg=hdg_list,
        eta_s=0.0,
        min_dist_nm=best_d if best_d is not None else float("inf"),
        motor_flags=_flags(mot_list, engine_below_kts),
    )


def corridor_mask(lat, lon, start, finish, margin, direct_nm=None):
    """Vrai pour les points dans l'ellipse de foyers start/finish dont le grand axe
    vaut (1 + margin) x distance directe (d(p, start) + d(p, finish) <= grand axe)."""
    if direct_nm is None:
        direct_nm = float(haversine_nm(start[0], start[1], finish[0], finish[1]))
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    total = haversine_nm(lat, lon, float(start[0]), float(start[1])) + haversine_nm(
        lat, lon, float(finish[0]), float(finish[1])
    )
    return total <= direct_nm * (1.0 + margin) + 1e-9


def cell_keys(lat, lon, bin_deg):
    """Clé entière de cellule de bin_deg degrés (int32 si elle y tient, sinon int64)."""
    n_lon = int(np.ceil(360.0 / bin_deg)) + 2
    n_lat = int(np.ceil(180.0 / bin_deg)) + 2
    dtype = np.int32 if n_lat * n_lon < np.iinfo(np.int32).max else np.int64
    ilat = np.floor((np.asarray(lat) + 90.0) / bin_deg).astype(dtype)
    ilon = np.floor((np.asarray(lon) + 180.0) / bin_deg).astype(dtype)
    return ilat * dtype(n_lon) + ilon


def _traceback(fronts, i):
    coords, bs_list, hdg_list, mot_list = [], [], [], []
    step, j = len(fronts) - 1, int(i)
    while step >= 0:
        f = fronts[step]
        coords.append((float(f["lat"][j]), float(f["lon"][j])))
        bs_list.append(float(f["bs"][j]))
        hdg_list.append(float(f["hdg"][j]))
        mot_list.append(bool(f["mot"][j]))
        j = int(f["parent"][j])
        step -= 1
    coords.reverse()
    bs_list.reverse()
    hdg_list.reverse()
    mot_list.reverse()
    return coords, bs_list, hdg_list, mot_list


def _flags(mot_list, engine_below_kts):
    """motor_flags exposé : None en pur voile (back-compat), sinon une entrée par point."""
    return list(mot_list) if engine_below_kts is not None else None


def to_geojson(res: RouteResult):
    """Export GeoJSON (LineString [lon, lat]) pour la carte Leaflet de MyNavApp."""
    if not res.found:
        raise ValueError("Aucune route trouvée (horizon dépassé)")
    props = {
        "eta_s": res.eta_s,
        "sog_kts": [round(s, 2) for s in res.sog_kts],
    }
    if res.motor_flags is not None:
        props["motor"] = [bool(m) for m in res.motor_flags]
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[lon, lat] for lat, lon in res.coords],
                },
                "properties": props,
            }
        ],
    }
