from pathlib import Path

import pytest

from mynav_routing import ConstantWind, Polar, route, to_geojson
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
    )
    assert not res.found
    assert res.min_dist_nm > 0.0
    with pytest.raises(ValueError):
        to_geojson(res)
