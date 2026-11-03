import io
import json
from pathlib import Path

import numpy as np
import pytest

from mynav_routing import ENGINE_SPEED_KTS, ConstantWind, Polar, route, to_geojson
from mynav_routing.engine.isochrones import haversine_nm, move_nm

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"


def test_geodesie():
    assert haversine_nm(0.0, 0.0, 0.0, 1.0) == pytest.approx(60.0, abs=0.5)
    lat, lon = move_nm(0.0, 0.0, 0.0, 60.0)
    assert lat == pytest.approx(1.0, abs=0.01)
    assert lon == pytest.approx(0.0, abs=1e-6)


def test_route_vent_arriere():
    # vent d'ouest 12 kt, cap vers l'est (~21 mn) : vent arrière sur la polaire Sugar
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=polar,
        wind=wind,
        horizon_h=24.0,
        dt=1800.0,
        bin_deg=0.05,
        arrival_nm=5.0,
        land=None,
    )
    assert res.found
    # ~21 mn à ~4.7 kt -> 3 à 6 h
    assert 3 * 3600 < res.eta_s < 6 * 3600
    assert res.coords[0] == (42.62, 3.02)
    geo = to_geojson(res)
    assert geo["type"] == "FeatureCollection"
    assert geo["features"][0]["geometry"]["type"] == "LineString"


def test_route_horizon_depasse():
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=polar,
        wind=wind,
        horizon_h=2.0,
        dt=1800.0,
        land=None,
    )
    assert not res.found
    assert res.min_dist_nm > 0.0
    with pytest.raises(ValueError):
        to_geojson(res)


def test_route_evite_la_terre():
    from mynav_routing import LandMask

    # île rectangulaire entre départ et arrivée
    island = [(42.55, 3.2), (42.55, 3.3), (42.70, 3.3), (42.70, 3.2)]
    land = LandMask([island], sample_nm=0.25)
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
    }
    res = route(land=land, **kw)
    assert res.found
    lat = [c[0] for c in res.coords]
    lon = [c[1] for c in res.coords]
    assert not land.crosses_land(lat[:-1], lon[:-1], lat[1:], lon[1:]).any()
    assert res.eta_s > route(**kw).eta_s


def test_landmask_geojson():
    from mynav_routing import LandMask

    ring = [[3.2, 42.5], [3.3, 42.5], [3.3, 42.7], [3.2, 42.7], [3.2, 42.5]]
    gj = {"type": "Polygon", "coordinates": [ring]}
    land = LandMask.from_geojson(gj)
    assert land.is_land(42.6, 3.25)
    assert not land.is_land(42.6, 3.1)


def test_landmask_charge_les_donnees_ouvertes(monkeypatch, tmp_path):
    import mynav_routing.land as land_module
    from mynav_routing.land import CACHE_DIR_ENV
    monkeypatch.setenv(CACHE_DIR_ENV, str(tmp_path / "cache"))

    geojson = {
        "type": "Polygon",
        "coordinates": [[[3.2, 42.5], [3.3, 42.5], [3.3, 42.7], [3.2, 42.5]]],
    }
    calls = []

    def open_url(url, timeout):
        calls.append((url, timeout))
        return io.BytesIO(json.dumps(geojson).encode())

    land_module._load_open_land_geojson.cache_clear()
    monkeypatch.setattr(land_module, "urlopen", open_url)
    land = land_module.LandMask.from_open_data()
    second_land = land_module.LandMask.from_open_data()
    assert land.is_land(42.6, 3.25)
    assert second_land.is_land(42.6, 3.25)
    assert calls == [(land_module.OPEN_LAND_URL, 30)]
    land_module._load_open_land_geojson.cache_clear()


def test_route_charge_automatiquement_le_masque(monkeypatch):
    from mynav_routing.land import LandMask

    calls = []

    def load(cls):
        calls.append(cls)
        return cls([])

    monkeypatch.setattr(LandMask, "from_open_data", classmethod(load))
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=Polar.from_file(str(EX)),
        wind=ConstantWind(direction_from_deg=270, speed_kts=12),
        horizon_h=24.0,
        dt=1800.0,
        bin_deg=0.05,
        arrival_nm=5.0,
    )
    assert res.found
    assert calls == [LandMask]


def test_corridor_mask_margin():
    from mynav_routing.router import corridor_mask

    start, finish = (42.0, 3.0), (42.0, 6.0)
    lat = [42.0, 42.0, 42.0, 44.0, 42.0]
    lon = [4.5, 3.0, 6.0, 4.5, 1.0]
    inside = corridor_mask(lat, lon, start, finish, 0.2)
    assert inside.tolist() == [True, True, True, False, False]
    assert corridor_mask(lat, lon, start, finish, 5.0).all()


def test_cell_keys_int32_and_distinct():
    from mynav_routing.router import cell_keys

    k = cell_keys([42.001, 42.002, 42.5], [3.001, 3.002, 3.0], 0.05)
    assert k.dtype == np.int32
    assert k[0] == k[1] != k[2]


def test_route_corridor_keeps_route():
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 3.5),
        "polar": polar,
        "wind": wind,
        "horizon_h": 24.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "land": None,
    }
    a = route(**kw)
    b = route(corridor_margin=None, **kw)
    assert a.found and b.found
    assert a.eta_s == pytest.approx(b.eta_s, rel=0.1)


def test_route_corridor_fallback_contournement():
    from mynav_routing import LandMask

    # Mur de terre coupant la ligne directe : le seul passage (au sud du mur)
    # sort de l'ellipse corridor (+20 %), comme un contournement de cap.
    # Le fallback doit retenter sans corridor et trouver la route.
    wall = [(42.3, 3.2), (42.3, 3.3), (43.0, 3.3), (43.0, 3.2)]
    land = LandMask([wall], sample_nm=0.25)
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=polar,
        wind=wind,
        horizon_h=24.0,
        dt=1800.0,
        bin_deg=0.05,
        arrival_nm=5.0,
        land=land,
    )
    assert res.found
    # La route vient du second essai : corridor effectif désactivé.
    assert res.corridor_margin is None
    # Elle passe au sud du mur.
    assert min(c[0] for c in res.coords) < 42.3


def test_route_corridor_effectif_quand_trouvee_au_premier_essai():
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 3.5),
        "polar": polar,
        "wind": wind,
        "horizon_h": 24.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "arrival_nm": 5.0,
        "land": None,
    }
    res = route(**kw)
    assert res.found
    assert res.corridor_margin == 0.2
    res_sans = route(corridor_margin=None, **kw)
    assert res_sans.found
    assert res_sans.corridor_margin is None


def test_route_moteur_seuil_parametrable():
    # B1 : vent arrière 4 kt -> polaire 2.1 kn au vent arrière. Seuil défaut
    # (2.0 kt) : la route reste en pur voile ; à 2.5 kt le moteur prend le
    # relais à 4.0 kn sur tous les tronçons ; None = pur voile (back-compat).
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=4)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 3.5),
        "polar": polar,
        "wind": wind,
        "horizon_h": 24.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "arrival_nm": 5.0,
        "land": None,
    }
    res_voile = route(engine_below_kts=None, **kw)
    assert res_voile.found
    assert res_voile.motor_flags is None  # pur voile : pas de flags (back-compat)
    assert "motor" not in to_geojson(res_voile)["features"][0]["properties"]

    res_def = route(**kw)
    assert res_def.found
    assert res_def.motor_flags is not None and not any(res_def.motor_flags)

    res_mot = route(engine_below_kts=2.5, **kw)
    assert res_mot.found
    assert res_mot.motor_flags is not None and all(res_mot.motor_flags)
    # alignement : une entrée par point, comme sog_kts et coords
    assert len(res_mot.motor_flags) == len(res_mot.sog_kts) == len(res_mot.coords)
    assert all(s == pytest.approx(ENGINE_SPEED_KTS) for s in res_mot.sog_kts)
    geo = to_geojson(res_mot)
    assert geo["features"][0]["properties"]["motor"] == res_mot.motor_flags


def test_route_moteur_zone_morte():
    # B1 défaut : au vent, la zone morte de la polaire (TWA < ~40°, vitesse 0)
    # passe en moteur à 4.0 kn : la route au vent devient directe, plus rapide
    # que le près serré en pur voile.
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=4)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 2.6),
        "polar": polar,
        "wind": wind,
        "horizon_h": 24.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "arrival_nm": 5.0,
        "land": None,
    }
    res = route(**kw)
    assert res.found
    assert res.motor_flags is not None and all(res.motor_flags)
    assert all(s == pytest.approx(ENGINE_SPEED_KTS) for s in res.sog_kts)
    res_voile = route(engine_below_kts=None, **kw)
    assert res_voile.found
    assert res_voile.motor_flags is None
    assert res.eta_s < res_voile.eta_s


def test_route_moteur_vent_navigable():
    # B1 défaut : vent navigable (12 kt, vent arrière ~4.7 kn) -> aucun
    # tronçon moteur, flags présents et tous faux.
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=polar,
        wind=wind,
        horizon_h=24.0,
        dt=1800.0,
        bin_deg=0.05,
        arrival_nm=5.0,
        land=None,
    )
    assert res.found
    assert res.motor_flags is not None and not any(res.motor_flags)
