"""Campagne d'acceptation du routage sur routes réelles (issue #30).

Le routage est le point le plus sensible de MyNavApp (constat terrain du 04/11 :
le corridor historique refusait des routes contournant un cap, mynav-routing #28).
Ces tests verrouillent les comportements attendus sur des géométries côtières
réelles (masque Natural Earth + polaire Sugar) :

- suite rapide (chaque push) : les routes terrain sensibles au vent de tramontane
  doivent toujours être trouvées, sans jamais couper la terre ;
- matrice complète (vents variés + traversée longue) : job opt-in
  ``routing-acceptance`` (workflow_dispatch), et script local de campagne
  (issue #30) pour les relevés de performance sur machine de référence.
"""

import time
from pathlib import Path

import numpy as np
import pytest

from mynav_routing import ConstantWind, LandMask, Polar, route
from mynav_routing.engine.isochrones import haversine_nm, move_nm

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"

# Routes terrain : (départ, arrivée, étiquette).
ROUTES_COURTES = [
    # Côte ouverte ~2,5 mn : contrôle, ne doit jamais régresser.
    ((42.481, 3.139), (42.477, 3.083), "cote-ouverte-banyuls-collioure"),
    # Contournement du Cap Béar (~3 mn) : le cas court signalé sur le terrain.
    ((42.500, 3.100), (42.535, 3.150), "contournement-cap-bear"),
    # Contournement du Cap de Creus (22,5 mn) : le cas mynav-routing #28.
    ((42.61932, 3.04720), (42.25284, 3.15549), "contournement-cap-de-creus"),
]

# Traversée longue du golfe du Lion (~190 mn), dt/bin_deg comme l'app (>150 mn).
TRAVERSEE = ((42.619, 3.047), (42.568, 8.759), "banyuls-calvi")

# Vents de campagne : tramontane modérée, levant au près, tramontane forte.
VENTS = [
    (270, 12, "tramontane-12"),
    (90, 12, "levant-12"),
    (270, 25, "tramontane-25"),
]


@pytest.fixture(scope="module")
def land():
    return LandMask.from_open_data()


@pytest.fixture(scope="module")
def polar():
    return Polar.from_file(str(EX))


def _nearest_water(land, lat, lon, max_nm=2.0, step_nm=0.1):
    """Ramène une extrémité « sur la terre » au point d'eau le plus proche.

    Même logique que le snappage de l'app : anneaux croissants, caps de 5 deg,
    2 mn maximum.
    """
    if not land.is_land([lat], [lon])[0]:
        return lat, lon
    bearings = np.arange(0.0, 360.0, 5.0)
    for radius_nm in np.arange(step_nm, max_nm + step_nm / 2.0, step_nm):
        plats, plons = move_nm(
            np.full(bearings.shape, float(lat)),
            np.full(bearings.shape, float(lon)),
            bearings,
            np.full(bearings.shape, float(radius_nm)),
        )
        water = ~land.is_land(plats, plons)
        if water.any():
            idx = int(np.flatnonzero(water)[0])
            return float(plats[idx]), float(plons[idx])
    pytest.fail(f"aucune eau à moins de {max_nm} mn de ({lat}, {lon}) selon le masque")


def _route_acceptation(land, polar, start, finish, direction, speed, dt, bin_deg, budget_s):
    """Route réelle avec extrémités snappées ; assertions communes à la campagne."""
    snapped_start = _nearest_water(land, *start)
    snapped_finish = _nearest_water(land, *finish)
    started = time.perf_counter()
    res = route(
        start=snapped_start,
        finish=snapped_finish,
        polar=polar,
        wind=ConstantWind(direction_from_deg=direction, speed_kts=speed),
        horizon_h=72.0,
        dt=dt,
        bin_deg=bin_deg,
        arrival_nm=2.0,
        land=land,
    )
    duration = time.perf_counter() - started
    print(
        f"\nACCEPTANCE {snapped_start} -> {snapped_finish} vent {direction}/{speed} kt "
        f"dt={dt} bin={bin_deg} : found={res.found} "
        f"min_dist={res.min_dist_nm:.2f} mn corridor={res.corridor_margin} "
        f"eta_h={res.eta_s / 3600.0:.1f} duree_s={duration:.1f}"
    )
    assert res.found, (
        f"route non trouvée : min_dist={res.min_dist_nm:.2f} mn, "
        f"corridor effectif={res.corridor_margin}"
    )
    assert 0.0 < res.eta_s <= 72.0 * 3600.0
    lat = [c[0] for c in res.coords]
    lon = [c[1] for c in res.coords]
    assert not land.crosses_land(lat[:-1], lon[:-1], lat[1:], lon[1:]).any(), (
        "la trace traverse la terre"
    )
    assert duration < budget_s, f"calcul trop long : {duration:.1f} s"
    return res, duration


@pytest.mark.parametrize(
    "start,finish,nom", ROUTES_COURTES, ids=[r[2] for r in ROUTES_COURTES]
)
def test_acceptation_route_courte_tramontane(land, polar, start, finish, nom):
    """Suite rapide : chaque route terrain doit être trouvée (tramontane 12 kt)."""
    _route_acceptation(land, polar, start, finish, 270, 12, 1800.0, 0.02, 300.0)


@pytest.mark.slow
@pytest.mark.parametrize("direction,speed,vent", VENTS, ids=[v[2] for v in VENTS])
@pytest.mark.parametrize(
    "start,finish,nom", ROUTES_COURTES, ids=[r[2] for r in ROUTES_COURTES]
)
def test_acceptance_matrice_routes_courtes(
    land, polar, start, finish, nom, direction, speed, vent
):
    """Matrice opt-in : routes courtes par trois régimes de vent."""
    _route_acceptation(land, polar, start, finish, direction, speed, 1800.0, 0.02, 600.0)


@pytest.mark.slow
@pytest.mark.parametrize("direction,speed,vent", VENTS, ids=[v[2] for v in VENTS])
def test_acceptance_traversee_longue(land, polar, direction, speed, vent):
    """Matrice opt-in : traversée Banyuls -> Calvi (~190 mn) dans l'horizon."""
    start, finish, _nom = TRAVERSEE
    res, _duration = _route_acceptation(
        land, polar, start, finish, direction, speed, 3600.0, 0.04, 1200.0
    )
    assert 24.0 * 3600.0 < res.eta_s, f"ETA invraisemblable : {res.eta_s / 3600.0:.1f} h"
    distance_nm = haversine_nm(*start, *finish)
    assert 170.0 <= distance_nm <= 210.0, f"distance inattendue : {distance_nm:.1f} mn"
