"""Tests du champ de vent GRIBWindField (unités m/s -> nœuds)."""

import numpy as np
import pytest

from mynav_routing import GRIBWindField

MS_TO_KNOTS = 1.94384


def make_field(u10_ms_by_step, v10_ms_by_step):
    """Champ GRIB in-memory : 1 point, un pas de temps par valeur fournie."""
    xr = pytest.importorskip("xarray")
    times = np.array(
        ["2026-01-01T00:00:00", "2026-01-01T06:00:00"], dtype="datetime64[ns]"
    )
    u10 = np.array([[[value]] for value in u10_ms_by_step])
    v10 = np.array([[[value]] for value in v10_ms_by_step])
    return xr.Dataset(
        {
            "u10": (("time", "latitude", "longitude"), u10),
            "v10": (("time", "latitude", "longitude"), v10),
        },
        coords={
            "time": times,
            "latitude": np.array([46.5]),
            "longitude": np.array([-4.5]),
        },
    )


def test_sample_convertit_ms_en_noeuds():
    # 5 m/s de vent (composante est) = 9.72 nœuds : même unité que ConstantWind.
    field = GRIBWindField(make_field([5.0, 5.0], [0.0, 0.0]))
    u, v = field.sample([46.5], [-4.5], 0)
    assert u[0] == pytest.approx(5.0 * MS_TO_KNOTS)
    assert v[0] == pytest.approx(0.0, abs=1e-9)


def test_sample_choisit_le_pas_de_temps_le_plus_proche():
    # t=0 -> premier pas (2 m/s), t=6 h -> second pas (10 m/s).
    field = GRIBWindField(make_field([2.0, 10.0], [0.0, 0.0]))
    u_first, _ = field.sample([46.5], [-4.5], 0)
    u_second, _ = field.sample([46.5], [-4.5], 6 * 3600)
    assert u_first[0] == pytest.approx(2.0 * MS_TO_KNOTS)
    assert u_second[0] == pytest.approx(10.0 * MS_TO_KNOTS)
