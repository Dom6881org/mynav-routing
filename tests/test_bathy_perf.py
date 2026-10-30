"""Tests de performance v0.3 : bathymétrie GEBCO (DepthMask) dans le routeur.

Budgets volontairement larges (runners CI partagés) : l'important est de détecter
une régression d'ordre de grandeur — par ex. la perte de l'index spatial, qui
ferait échantillonner tous les tronçons à chaque pas d'isochrone.
"""

import time
from pathlib import Path

import numpy as np
import pytest

from mynav_routing import ConstantWind, DepthMask, Polar, route
from mynav_routing.engine.isochrones import move_nm

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"

# grille régionale 15° x 25° à 0.02° (~940k mailles) : côte ondulée au nord,
# île et haut-fond entre le départ et l'arrivée
LATS = np.arange(35.0, 50.001, 0.02)
LONS = np.arange(-10.0, 15.001, 0.02)


def _elev():
    gy, gx = np.meshgrid(LATS, LONS, indexing="ij")
    elev = np.full(gy.shape, -4000.0, dtype=np.float32)
    elev[(gy >= 42.55) & (gy <= 42.70) & (gx >= 3.20) & (gx <= 3.30)] = 5.0   # île
    elev[(gy >= 42.58) & (gy <= 42.66) & (gx >= 3.35) & (gx <= 3.45)] = -3.0  # haut-fond
    coast = 43.2 + 0.5 * np.sin(gx * 2.0)                                       # côte au nord
    elev[gy > coast] = 10.0
    return elev


@pytest.fixture(scope="module")
def sea():
    t0 = time.perf_counter()
    mask = DepthMask(LATS, LONS, _elev(), cell_deg=1.0)
    build_s = time.perf_counter() - t0
    print(f"\nconstruction DepthMask {mask.elev.shape} : {build_s:.3f}s")
    assert build_s < 5.0  # index spatial calculé une fois au chargement
    return mask


def test_perf_elevation_at_vectorise(sea):
    # interpolation bilinéaire : O(1) par point, ~10^7 points/s attendu
    rng = np.random.default_rng(0)
    lat = rng.uniform(35.0, 50.0, 1_000_000)
    lon = rng.uniform(-10.0, 15.0, 1_000_000)
    t0 = time.perf_counter()
    e = sea.elevation_at(lat, lon)
    dt = time.perf_counter() - t0
    assert np.isfinite(e).all()
    print(f"\nelevation_at 1M points : {dt:.3f}s")
    assert dt < 2.0, f"elevation_at 1M points trop lent : {dt:.3f}s"


def test_perf_crosses_shallow_au_large_pre_filtre(sea):
    # 200k tronçons de ~2 mn au large, loin de toute cellule dangereuse :
    # le pré-filtre spatial doit éviter tout échantillonnage
    rng = np.random.default_rng(1)
    n = 200_000
    lat1 = rng.uniform(36.0, 40.0, n)
    lon1 = rng.uniform(-8.0, 0.0, n)
    lat2, lon2 = move_nm(lat1, lon1, 90.0, 2.0)
    t0 = time.perf_counter()
    hit = sea.crosses_shallow(lat1, lon1, lat2, lon2, 10.0)
    dt = time.perf_counter() - t0
    assert not hit.any()
    print(f"\ncrosses_shallow 200k tronçons au large : {dt:.3f}s")
    assert dt < 5.0, f"pré-filtre inefficace : {dt:.3f}s pour 200k tronçons au large"


def test_perf_crosses_shallow_pres_du_haut_fond(sea):
    # 50k tronçons qui coupent le haut-fond : échantillonnés, tous détectés
    n = 50_000
    lat1 = np.full(n, 42.60)
    lon1 = np.full(n, 3.30)
    lat2, lon2 = move_nm(lat1, lon1, 90.0, 8.0)
    t0 = time.perf_counter()
    hit = sea.crosses_shallow(lat1, lon1, lat2, lon2, 10.0)
    dt = time.perf_counter() - t0
    assert hit.all()
    print(f"\ncrosses_shallow 50k tronçons sur le haut-fond : {dt:.3f}s")
    assert dt < 10.0, f"échantillonnage côtier trop lent : {dt:.3f}s"


def test_perf_route_avec_bathy(sea):
    # routage complet avec contrainte de profondeur : l'île et le haut-fond
    # sont entre le départ et l'arrivée, la route doit les contourner
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 3.5),
        "polar": polar,
        "wind": wind,
        "horizon_h": 48.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "arrival_nm": 5.0,
        "land": None,
    }
    t0 = time.perf_counter()
    res_free = route(bathy=None, **kw)
    t_free = time.perf_counter() - t0
    t0 = time.perf_counter()
    res = route(bathy=sea, min_depth_m=10.0, **kw)
    t_bathy = time.perf_counter() - t0
    print(f"\nroute libre : {t_free:.2f}s (eta {res_free.eta_s / 3600:.2f}h) — "
          f"route avec bathy : {t_bathy:.2f}s (eta {res.eta_s / 3600:.2f}h)")
    assert res.found
    lat = np.array([c[0] for c in res.coords])
    lon = np.array([c[1] for c in res.coords])
    assert not sea.crosses_shallow(lat[:-1], lon[:-1], lat[1:], lon[1:], 10.0).any()
    assert (sea.elevation_at(lat, lon) <= -10.0).all()
    assert res.eta_s >= res_free.eta_s - 1e-6
    # la contrainte bathymétrique ne doit pas coûter un ordre de grandeur
    assert t_bathy < 60.0, f"routage avec bathy trop lent : {t_bathy:.1f}s"
    assert t_bathy < 10.0 * max(t_free, 0.5), (
        f"surcoût bathy excessif : {t_bathy:.2f}s vs {t_free:.2f}s sans contrainte"
    )
