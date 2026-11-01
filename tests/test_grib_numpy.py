"""Tests for numpy-backed GRIB wind sampling."""

import time

import numpy as np
import pytest

from mynav_routing import GRIBWindField

MS_TO_KNOTS = 1.94384


def make_dataset(times, latitudes, longitudes, u10, v10):
    xr = pytest.importorskip("xarray")
    return xr.Dataset(
        {
            "u10": (("time", "latitude", "longitude"), u10),
            "v10": (("time", "latitude", "longitude"), v10),
        },
        coords={
            "time": times,
            "latitude": latitudes,
            "longitude": longitudes,
        },
    )


def test_sample_uses_numpy_nearest_after_initialization(monkeypatch):
    times = np.array(["2026-01-01T00", "2026-01-01T06"], dtype="datetime64[h]")
    latitudes = np.array([50.0, 45.0, 40.0])
    longitudes = np.array([-10.0, 0.0, 10.0])
    u10 = np.arange(18, dtype=float).reshape(2, 3, 3)
    v10 = -u10
    dataset = make_dataset(times, latitudes, longitudes, u10, v10)
    field = GRIBWindField(dataset)

    def fail_xarray_access(*args, **kwargs):
        raise AssertionError("sample() must not access xarray")

    monkeypatch.setattr(type(dataset), "isel", fail_xarray_access)
    monkeypatch.setattr(type(dataset), "sel", fail_xarray_access)

    u, v = field.sample([49.8, 40.2], [-9.8, 9.7], 5 * 3600)

    np.testing.assert_allclose(u, [u10[1, 0, 0] * MS_TO_KNOTS, u10[1, 2, 2] * MS_TO_KNOTS])
    np.testing.assert_allclose(v, [v10[1, 0, 0] * MS_TO_KNOTS, v10[1, 2, 2] * MS_TO_KNOTS])

    scalar_u, scalar_v = field.sample(49.8, -9.8, 0)
    np.testing.assert_allclose(scalar_u, [u10[0, 0, 0] * MS_TO_KNOTS])
    np.testing.assert_allclose(scalar_v, [v10[0, 0, 0] * MS_TO_KNOTS])


@pytest.mark.slow
def test_sample_numpy_is_at_least_ten_times_faster_than_xarray():
    xr = pytest.importorskip("xarray")
    rng = np.random.default_rng(42)
    times = np.arange(8).astype("timedelta64[h]") + np.datetime64("2026-01-01")
    latitudes = np.linspace(35, 55, 80)
    longitudes = np.linspace(-15, 20, 100)
    u10 = rng.random((8, 80, 100))
    v10 = rng.random((8, 80, 100))
    dataset = make_dataset(times, latitudes, longitudes, u10, v10)
    field = GRIBWindField(dataset)
    lat = rng.uniform(latitudes.min(), latitudes.max(), 512)
    lon = rng.uniform(longitudes.min(), longitudes.max(), 512)

    def xarray_sample():
        target = times[0] + np.timedelta64(3600, "s")
        time_idx = int(np.argmin(np.abs(times - target)))
        points = dataset.isel(time=time_idx).sel(
            latitude=xr.DataArray(lat, dims="pts"),
            longitude=xr.DataArray(lon, dims="pts"),
            method="nearest",
        )
        return points["u10"].values * MS_TO_KNOTS, points["v10"].values * MS_TO_KNOTS

    for _ in range(5):
        field.sample(lat, lon, 3600)
        xarray_sample()

    numpy_samples = []
    xarray_samples = []
    for _ in range(5):
        started = time.perf_counter()
        for _ in range(20):
            field.sample(lat, lon, 3600)
        numpy_samples.append((time.perf_counter() - started) / 20)

        started = time.perf_counter()
        for _ in range(20):
            xarray_sample()
        xarray_samples.append((time.perf_counter() - started) / 20)

    numpy_seconds = float(np.median(numpy_samples))
    xarray_seconds = float(np.median(xarray_samples))
    speedup = xarray_seconds / numpy_seconds
    print(
        f"GRIB_SAMPLE_BENCH points=512 numpy_ms={numpy_seconds * 1000:.3f} "
        f"xarray_ms={xarray_seconds * 1000:.3f} speedup={speedup:.1f}x"
    )
    assert speedup >= 10
