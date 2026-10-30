"""Champs de vent pour le moteur de routage.

Le moteur dépend uniquement de l'interface sample(lat, lon, t) -> (u, v)
en nœuds (composantes est/nord du vecteur de déplacement du vent).
t est en secondes depuis le départ du routage.
"""
from __future__ import annotations

import numpy as np


class WindField:
    def sample(self, lat, lon, t):
        raise NotImplementedError


class ConstantWind(WindField):
    """Vent uniforme : tests, scénarios, dégradation offline."""

    def __init__(self, direction_from_deg, speed_kts):
        d = np.radians(direction_from_deg)
        self.u = -speed_kts * np.sin(d)
        self.v = -speed_kts * np.cos(d)

    def sample(self, lat, lon, t):
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        return (np.full(lat.shape, self.u), np.full(lat.shape, self.v))


class GRIBWindField(WindField):
    """Champ GRIB (u10/v10, m/s) lu via xarray + cfgrib (extra 'grib').

    Expérimental v1 : pas de temps le plus proche, voisin spatial le plus proche.
    """

    def __init__(self, ds):
        self.ds = ds

    @classmethod
    def from_grib(cls, path):
        import xarray as xr

        return cls(xr.open_dataset(path, engine="cfgrib"))

    def sample(self, lat, lon, t):
        import xarray as xr

        ds = self.ds
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        times = ds["time"].values.astype("datetime64[ns]")
        target = times[0] + np.timedelta64(int(t * 1e9), "ns")
        idx = int(np.argmin(np.abs(times - target)))
        pt = ds.isel(time=idx).sel(
            latitude=xr.DataArray(lat, dims="pts"),
            longitude=xr.DataArray(lon, dims="pts"),
            method="nearest",
        )
        u = np.atleast_1d(pt["u10"].values) / 1.94384
        v = np.atleast_1d(pt["v10"].values) / 1.94384
        return u, v
