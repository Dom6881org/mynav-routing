"""Masque de terre : empêche le routeur de traverser la côte.

Données attendues : polygones de terre (GeoJSON Polygon / MultiPolygon, coordonnées
[lon, lat]), par exemple OpenStreetMap land polygons ou GSHHG.

Gestion hors ligne (exécution 100 % offline à bord) : le GeoJSON téléchargé est
mis en cache sur disque ; hors connexion, le cache est utilisé en repli. Un fichier
local peut être imposé via MYNAV_ROUTING_LAND_FILE, le répertoire de cache via
MYNAV_ROUTING_CACHE_DIR (défaut : ~/.cache/mynav-routing). Si aucune donnée n'est
disponible, LandDataError est levée.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from functools import lru_cache
from pathlib import Path
from urllib.request import urlopen

import numpy as np

from .engine.isochrones import NM_PER_RADIAN

OPEN_LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_10m_land.geojson"
)

LAND_FILE_ENV = "MYNAV_ROUTING_LAND_FILE"
CACHE_DIR_ENV = "MYNAV_ROUTING_CACHE_DIR"


class LandDataError(ValueError):
    """Masque de terre indisponible (hors ligne et aucun cache local).

    Sous-classe de ValueError : le backend MyNavApp traduit déjà les ValueError
    levées par route() en message utilisateur (RoutingError).
    """


def _cache_dir() -> Path:
    override = os.environ.get(CACHE_DIR_ENV)
    if override:
        return Path(override)
    return Path.home() / ".cache" / "mynav-routing"


def _cache_path(url: str) -> Path:
    digest = hashlib.sha1(url.encode("utf-8")).hexdigest()[:16]
    return _cache_dir() / f"{digest}.geojson"


def _write_cache(path: Path, data) -> None:
    """Écrit le cache de façon atomique ; best effort (répertoire en lecture seule = ignoré)."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name + ".", suffix=".tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(data, handle)
            os.replace(tmp, path)
        except BaseException:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
    except OSError:
        pass


@lru_cache(maxsize=4)
def _load_open_land_geojson(url):
    # 1. Fichier local explicite (prioritaire : ex. GeoJSON embarqué à bord).
    local = os.environ.get(LAND_FILE_ENV)
    if local:
        try:
            return json.loads(Path(local).read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise LandDataError(f"{LAND_FILE_ENV}={local} illisible : {error}") from error
    cache = _cache_path(url)
    # 2. Téléchargement, puis mise à jour du cache disque.
    try:
        with urlopen(url, timeout=30) as response:
            data = json.loads(response.read())
    except (OSError, ValueError):
        data = None
    if data is not None:
        _write_cache(cache, data)
        return data
    # 3. Hors ligne : repli sur le cache disque.
    if cache.is_file():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    raise LandDataError(
        "Masque de terre indisponible : hors connexion et aucun cache local "
        f"({cache}). Lancez un premier routage en ligne pour remplir le cache, "
        f"ou définissez {LAND_FILE_ENV} vers un fichier GeoJSON de terres."
    )


class LandMask:
    """Ensemble de polygones de terre ; is_land / crosses_land sont vectorisés."""

    def __init__(self, polygons, sample_nm=0.25):
        """polygons : liste d'anneaux extérieurs, chacun [(lat, lon), ...].
        sample_nm : pas d'échantillonnage le long d'un tronçon (mn)."""
        self.sample_nm = float(sample_nm)
        self._rings = []
        for ring in polygons:
            arr = np.asarray(ring, dtype=float)
            if len(arr) >= 3:
                bbox = (arr[:, 0].min(), arr[:, 0].max(), arr[:, 1].min(), arr[:, 1].max())
                self._rings.append((arr, bbox))

    @classmethod
    def from_geojson(cls, data, sample_nm=0.25):
        """Construit un masque depuis un dict GeoJSON (Feature, FeatureCollection ou géométrie)."""
        rings = []

        def walk(obj):
            t = obj.get("type")
            if t == "FeatureCollection":
                for f in obj["features"]:
                    walk(f)
            elif t == "Feature":
                walk(obj["geometry"])
            elif t == "GeometryCollection":
                for g in obj["geometries"]:
                    walk(g)
            elif t == "Polygon":
                rings.append([(lat, lon) for lon, lat, *_ in obj["coordinates"][0]])
            elif t == "MultiPolygon":
                for poly in obj["coordinates"]:
                    rings.append([(lat, lon) for lon, lat, *_ in poly[0]])

        walk(data)
        return cls(rings, sample_nm=sample_nm)

    @classmethod
    def from_open_data(cls, sample_nm=0.25, url=OPEN_LAND_URL):
        """Charge les terres Natural Earth : fichier local MYNAV_ROUTING_LAND_FILE,
        sinon téléchargement mis en cache, sinon repli sur le cache (hors ligne)."""
        return cls.from_geojson(_load_open_land_geojson(url), sample_nm=sample_nm)

    def is_land(self, lat, lon):
        """Vrai pour les points situés sur la terre (tableaux numpy acceptés)."""
        lat = np.asarray(lat, dtype=float)
        lon = np.asarray(lon, dtype=float)
        out = np.zeros(lat.shape, dtype=bool)
        for ring, (la0, la1, lo0, lo1) in self._rings:
            cand = ~out & (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)
            if not cand.any():
                continue
            y, x = lat[cand], lon[cand]
            inside = np.zeros(y.shape, dtype=bool)
            ys, xs = ring[:, 0], ring[:, 1]
            for i in range(len(ring)):
                j = i - 1
                crosses = (ys[i] > y) != (ys[j] > y)
                with np.errstate(divide="ignore", invalid="ignore"):
                    xint = (xs[j] - xs[i]) * (y - ys[i]) / (ys[j] - ys[i]) + xs[i]
                inside ^= crosses & (x < xint)
            out[cand] = inside
        return out

    def crosses_land(self, lat1, lon1, lat2, lon2):
        """Vrai si le tronçon (lat1, lon1) -> (lat2, lon2) touche la terre.

        Le point de départ n'est pas testé ; l'arrivée et les points intermédiaires le sont."""
        lat1, lon1, lat2, lon2 = np.broadcast_arrays(
            *(np.asarray(a, dtype=float) for a in (lat1, lon1, lat2, lon2))
        )
        if not self._rings or lat1.size == 0:
            return np.zeros(lat1.shape, dtype=bool)
        dlon = (lon2 - lon1 + 180.0) % 360.0 - 180.0
        dist = np.hypot(lat2 - lat1, dlon * np.cos(np.radians((lat1 + lat2) / 2.0)))
        dist_nm = float(dist.max()) * NM_PER_RADIAN * np.pi / 180.0
        n = max(1, int(np.ceil(dist_nm / self.sample_nm)))
        hit = np.zeros(lat1.shape, dtype=bool)
        for s in range(1, n + 1):
            f = s / n
            hit |= self.is_land(lat1 + (lat2 - lat1) * f, lon1 + dlon * f)
        return hit
