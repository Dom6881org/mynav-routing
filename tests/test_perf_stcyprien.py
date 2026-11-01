"""Benchmark reproductible du routage Saint-Cyprien → Ajaccio."""

import time
from pathlib import Path

import pytest

from mynav_routing import ConstantWind, LandMask, Polar, route
from mynav_routing.engine.isochrones import haversine_nm

pytestmark = pytest.mark.slow

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"
START = (42.62, 3.03)
FINISH = (41.919, 8.738)


def test_perf_stcyprien_ajaccio():
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    # Contours fixes et grossiers pour garder le benchmark indépendant du réseau
    # et des versions de Natural Earth ; la mer du corridor reste navigable.
    land = LandMask(
        [
            [(42.9, 0.0), (44.0, 0.0), (44.0, 10.0), (42.9, 10.0)],
            [(41.95, 9.0), (42.0, 9.5), (42.7, 9.4), (42.8, 9.0)],
        ],
        sample_nm=0.5,
        cell_deg=1.0,
        raster_deg=0.02,
    )
    distance_nm = haversine_nm(*START, *FINISH)
    assert 250.0 <= distance_nm <= 265.0

    phase_timings = {}
    started = time.perf_counter()
    result = route(
        start=START,
        finish=FINISH,
        polar=polar,
        wind=wind,
        dt=1800.0,
        horizon_h=72.0,
        heading_step=10.0,
        bin_deg=0.1,
        arrival_nm=5.0,
        land=land,
        phase_timings=phase_timings,
    )
    total_s = time.perf_counter() - started

    print(
        "\nBENCH_STCYPRIEN_AJACCIO "
        f"total_s={total_s:.3f} "
        f"distance_nm={distance_nm:.1f} "
        f"explored_nodes={phase_timings['explored_nodes']} "
        f"sample_wind_s={phase_timings['sample_wind_s']:.3f} "
        f"land_mask_s={phase_timings['land_mask_s']:.3f} "
        f"pruning_s={phase_timings['pruning_s']:.3f}"
    )
    assert result.found
    assert phase_timings["explored_nodes"] > 0
    assert phase_timings["sample_wind_s"] > 0
    assert phase_timings["land_mask_s"] > 0
    assert phase_timings["pruning_s"] > 0
    assert not land.crosses_land(
        [point[0] for point in result.coords[:-1]],
        [point[1] for point in result.coords[:-1]],
        [point[0] for point in result.coords[1:]],
        [point[1] for point in result.coords[1:]],
    ).any()
    assert total_s < 60.0, f"routage Saint-Cyprien → Ajaccio trop lent : {total_s:.1f}s"
