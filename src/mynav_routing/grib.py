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
    """Champ GRIB (u10/v10 en m/s) lu via xarray + cfgrib (extra 'grib').

    Expérimental v1 : pas de temps le plus proche, voisin spatial le plus proche.
    Les données sont extraites en tableaux numpy à l'initialisation pour éviter
    les accès xarray dans sample(). Les composantes m/s sont converties en nœuds
    (x 1.94384), comme ConstantWind.
    """

    def __init__(self, ds):
        self.times = np.array(ds["time"].values, dtype="datetime64[ns]", copy=True)
        self.latitude = np.array(ds["latitude"].values, copy=True)
        self.longitude = np.array(ds["longitude"].values, copy=True)
        self.u = np.array(
            ds["u10"].transpose("time", "latitude", "longitude").values, copy=True
        )
        self.v = np.array(
            ds["v10"].transpose("time", "latitude", "longitude").values, copy=True
        )

    @classmethod
    def from_grib(cls, path):
        import xarray as xr

        return cls(xr.open_dataset(path, engine="cfgrib"))

    def covers(self, lat, lon):
        """True si la position est dans l'emprise (lat/lon) du champ.

        sample() nearest-clippe silencieusement toute position hors emprise
        (vent du bord de grille, sans erreur) : covers() permet aux appelants
        de détecter ce cas avant d'interpréter un résultat (MyNavApp #141).
        Vectorisée comme sample() ; scalaires acceptés.
        """
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        return (
            (lat >= self.latitude.min())
            & (lat <= self.latitude.max())
            & (lon >= self.longitude.min())
            & (lon <= self.longitude.max())
        )

    def sample(self, lat, lon, t):
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        target = self.times[0] + np.timedelta64(int(t * 1e9), "ns")
        time_idx = int(np.argmin(np.abs(self.times - target)))
        lat_idx = self._nearest_indices(self.latitude, lat)
        lon_idx = self._nearest_indices(self.longitude, lon)
        return self.u[time_idx, lat_idx, lon_idx] * 1.94384, self.v[
            time_idx, lat_idx, lon_idx
        ] * 1.94384

    @staticmethod
    def _nearest_indices(coordinates, values):
        ascending = coordinates[0] <= coordinates[-1]
        ordered = coordinates if ascending else coordinates[::-1]
        positions = np.searchsorted(ordered, values)
        left = np.clip(positions - 1, 0, len(ordered) - 1)
        right = np.clip(positions, 0, len(ordered) - 1)
        choose_left = np.abs(values - ordered[left]) < np.abs(values - ordered[right])
        indices = np.where(choose_left, left, right)
        return indices if ascending else len(ordered) - 1 - indices
