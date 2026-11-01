"""Performance checks for vectorized land-mask segment filtering."""

import time

import numpy as np
import pytest

from mynav_routing import LandMask
from mynav_routing.engine.isochrones import NM_PER_RADIAN


def _legacy_crosses_land(land, lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = np.broadcast_arrays(lat1, lon1, lat2, lon2)
    shape = lat1.shape
    la1, lo1, la2, lo2 = (a.ravel() for a in (lat1, lon1, lat2, lon2))
    dlon = (lo2 - lo1 + 180.0) % 360.0 - 180.0
    seg_la0 = np.minimum(la1, la2)
    seg_la1 = np.maximum(la1, la2)
    seg_lo0 = np.minimum(lo1, lo1 + dlon)
    seg_lo1 = np.maximum(lo1, lo1 + dlon)
    i_min = int(np.floor(seg_la0.min() / land.cell_deg))
    i_max = int(np.floor(seg_la1.max() / land.cell_deg))
    j_min = int(np.floor(seg_lo0.min() / land.cell_deg))
    j_max = int(np.floor(seg_lo1.max() / land.cell_deg))
    candidates = {
        idx
        for i in range(i_min, i_max + 1)
        for j in range(j_min, j_max + 1)
        for idx in land._grid.get((i, j), ())
    }
    near = np.zeros(la1.shape, dtype=bool)
    for idx in candidates:
        _ring, (ra0, ra1, ro0, ro1) = land._rings[idx]
        near |= (
            (seg_la1 >= ra0)
            & (seg_la0 <= ra1)
            & (seg_lo1 >= ro0)
            & (seg_lo0 <= ro1)
        )
    hit = np.zeros(la1.shape, dtype=bool)
    selected = np.flatnonzero(near)
    if selected.size:
        dla = la2[selected] - la1[selected]
        dlo = dlon[selected]
        dist = np.hypot(dla, dlo * np.cos(np.radians((la1[selected] + la2[selected]) / 2.0)))
        n = max(1, int(np.ceil(float(dist.max()) * NM_PER_RADIAN * np.pi / 180.0 / land.sample_nm)))
        fs = (np.arange(1, n + 1, dtype=float) / n)[:, None]
        points_lat = (la1[selected] + dla * fs).ravel()
        points_lon = (lo1[selected] + dlo * fs).ravel()
        hit[selected] = land.is_land(points_lat, points_lon).reshape(n, selected.size).any(axis=0)
    return hit.reshape(shape)


@pytest.mark.slow
def test_crosses_land_rich_mask_benchmark():
    rng = np.random.default_rng(26)
    rings = []
    for row in range(40):
        for col in range(20):
            lat = 35.0 + row * 0.24
            lon = -5.0 + col * 0.5
            rings.append([(lat, lon), (lat, lon + 0.06), (lat + 0.06, lon + 0.06), (lat + 0.06, lon)])
    land = LandMask(rings, sample_nm=0.25, cell_deg=0.5, raster_deg=0.1)

    lat1 = rng.uniform(39.8, 40.2, 100_000)
    lon1 = rng.uniform(-0.2, 0.2, 100_000)
    lat2 = lat1 + 0.0001
    lon2 = lon1 + 0.0001
    lat1 = np.r_[lat1, [34.8, 34.8, 44.8, 44.8]]
    lon1 = np.r_[lon1, [-5.2, 5.2, -5.2, 5.2]]
    lat2 = np.r_[lat2, lat1[-4:] + 0.0001]
    lon2 = np.r_[lon2, lon1[-4:] + 0.0001]
    args = (lat1, lon1, lat2, lon2)

    reference = _legacy_crosses_land(land, *args)
    actual = land.crosses_land(*args)
    np.testing.assert_array_equal(actual, reference)

    legacy_times = []
    optimized_times = []
    for _ in range(3):
        started = time.perf_counter()
        _legacy_crosses_land(land, *args)
        legacy_times.append(time.perf_counter() - started)

        started = time.perf_counter()
        land.crosses_land(*args)
        optimized_times.append(time.perf_counter() - started)

    legacy_seconds = float(np.median(legacy_times))
    optimized_seconds = float(np.median(optimized_times))
    speedup = legacy_seconds / optimized_seconds
    print(
        f"LAND_CROSSES_BENCH rings={len(rings)} segments={lat1.size} "
        f"legacy_s={legacy_seconds:.3f} optimized_s={optimized_seconds:.3f} "
        f"speedup={speedup:.1f}x"
    )
    assert speedup >= 3.0
