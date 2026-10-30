"""Masque de terre : empêche le routeur de traverser la côte.

Données attendues : polygones de terre (GeoJSON Polygon / MultiPolygon, coordonnées
[lon, lat]), par exemple OpenStreetMap land polygons ou GSHHG.

Gestion hors ligne (exécution 100 % offline à bord) : le GeoJSON téléchargé est
mis en cache sur disque ; hors connexion, le cache est utilisé en repli. Un fichier
local peut être imposé via MYNAV_ROUTING_LAND_FILE, le répertoire de cache via
MYNAV_ROUTING_CACHE_DIR (défaut : ~/.cache/mynav-routing). Si aucune donnée n'est
disponible, LandDataError est levée.

Performance : un index spatial (grille de cell_deg degrés) restreint chaque requête
aux anneaux dont la bounding box chevauche les cellules interrogées — sinon is_land
balaie les ~7000 anneaux de Natural Earth à chaque appel d'isochrone.
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
    """Ensemble de polygones de terre ; is_land / crosses_land sont vectorisés.

    L'index spatial (grille de cell_deg degrés) accélère is_land : seules les
    cellules interrogées sont consultées, pas les milliers d'anneaux du masque.
    """

    def __init__(self, polygons, sample_nm=0.25, cell_deg=2.0, raster_deg=0.01):
        """polygons : liste d'anneaux extérieurs, chacun [(lat, lon), ...].
        sample_nm : pas d'échantillonnage le long d'un tronçon (mn).
        cell_deg : taille des cellules de l'index spatial (degrés).
        raster_deg : résolution de la rasterisation paresseuse par cellule,
            utilisée par is_land (lookup O(1) par point) ; None = calcul exact.
            Tolérance induite : ~raster_deg/2 près des côtes."""
        self.sample_nm = float(sample_nm)
        self.cell_deg = float(cell_deg)
        self.raster_deg = None if raster_deg is None else float(raster_deg)
        self._raster = {}
        self._rings = []
        for ring in polygons:
            arr = np.asarray(ring, dtype=float)
            if len(arr) >= 3:
                bbox = (arr[:, 0].min(), arr[:, 0].max(), arr[:, 1].min(), arr[:, 1].max())
                self._rings.append((arr, bbox))
        # index spatial : cellule (i, j) -> indices des anneaux dont la bbox la chevauche
        self._grid = {}
        for idx, (_, (la0, la1, lo0, lo1)) in enumerate(self._rings):
            i0 = int(np.floor(la0 / self.cell_deg))
            i1 = int(np.floor(la1 / self.cell_deg))
            j0 = int(np.floor(lo0 / self.cell_deg))
            j1 = int(np.floor(lo1 / self.cell_deg))
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    self._grid.setdefault((i, j), []).append(idx)
        # index des arêtes par cellule de latitude (1°) : le ray casting d'un point
        # ne teste que les arêtes de sa cellule de latitude — un anneau géant
        # (ex. Eurasie, ~80k sommets dans Natural Earth) est ainsi évalué en
        # ~V/180 arêtes par point au lieu de V, sans changer le résultat.
        self._edge_index = []
        for arr, _bbox in self._rings:
            ys, xs = arr[:, 0], arr[:, 1]
            yp, xp = np.roll(ys, 1), np.roll(xs, 1)  # sommet précédent (j = i - 1)
            lo = np.minimum(ys, yp)
            hi = np.maximum(ys, yp)
            by_cell = {}
            for i in range(len(ys)):
                for c in range(int(np.floor(lo[i])), int(np.floor(hi[i])) + 1):
                    by_cell.setdefault(c, []).append(i)
            self._edge_index.append({
                c: (ys[np.asarray(lst)], xs[np.asarray(lst)],
                    yp[np.asarray(lst)], xp[np.asarray(lst)])
                for c, lst in by_cell.items()
            })

    @classmethod
    def from_geojson(cls, data, sample_nm=0.25, cell_deg=2.0, raster_deg=0.01):
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
        return cls(rings, sample_nm=sample_nm, cell_deg=cell_deg, raster_deg=raster_deg)

    @classmethod
    def from_open_data(cls, sample_nm=0.25, url=OPEN_LAND_URL, cell_deg=2.0, raster_deg=0.01):
        """Charge les terres Natural Earth : fichier local MYNAV_ROUTING_LAND_FILE,
        sinon téléchargement mis en cache, sinon repli sur le cache (hors ligne)."""
        return cls.from_geojson(
            _load_open_land_geojson(url), sample_nm=sample_nm, cell_deg=cell_deg, raster_deg=raster_deg
        )

    def _candidate_rings(self, lat, lon):
        """Indices des anneaux dont une cellule spatiale peut contenir les points."""
        ci = np.floor(lat / self.cell_deg).astype(np.int64)
        cj = np.floor(lon / self.cell_deg).astype(np.int64)
        cells = set(zip(ci.ravel().tolist(), cj.ravel().tolist()))
        seen = {}
        for cell in cells:
            for idx in self._grid.get(cell, ()):
                seen[idx] = True
        return list(seen)

    def is_land(self, lat, lon):
        """Vrai pour les points situés sur la terre (tableaux numpy acceptés).

        Par défaut : rasterisation paresseuse par cellule (self._raster), lookup
        O(1) par point. raster_deg=None : ray casting exact sur les anneaux candidats."""
        lat = np.asarray(lat, dtype=float)
        lon = np.asarray(lon, dtype=float)
        out = np.zeros(lat.shape, dtype=bool)
        if lat.size == 0:
            return out
        if self.raster_deg is not None:
            flat_out = out.ravel()
            flat_lat = lat.ravel()
            flat_lon = lon.ravel()
            ci = np.floor(flat_lat / self.cell_deg).astype(np.int64)
            cj = np.floor(flat_lon / self.cell_deg).astype(np.int64)
            for cell in set(zip(ci.tolist(), cj.tolist())):
                raster = self._raster.get(cell)
                if raster is None:
                    raster = self._rasterize(cell)
                sel = (ci == cell[0]) & (cj == cell[1])
                n_lat, n_lon = raster.shape
                li = np.clip(
                    np.floor((flat_lat[sel] - cell[0] * self.cell_deg) / self.raster_deg).astype(np.int64),
                    0, n_lat - 1,
                )
                lj = np.clip(
                    np.floor((flat_lon[sel] - cell[1] * self.cell_deg) / self.raster_deg).astype(np.int64),
                    0, n_lon - 1,
                )
                flat_out[sel] = raster[li, lj]
            return out
        flat_out = out.ravel()
        flat_lat = lat.ravel()
        flat_lon = lon.ravel()
        for idx in self._candidate_rings(lat, lon):
            _ring, (la0, la1, lo0, lo1) = self._rings[idx]
            cand = ~flat_out & (flat_lat >= la0) & (flat_lat <= la1) & (flat_lon >= lo0) & (flat_lon <= lo1)
            if not cand.any():
                continue
            positions = np.flatnonzero(cand)
            y, x = flat_lat[positions], flat_lon[positions]
            inside = _ring_contains(self._edge_index[idx], y, x)
            flat_out[positions[inside]] = True
        return out

    def _rasterize(self, cell):
        """Rasterise une cellule de cell_deg degrés à raster_deg : ray casting exact
        des points de la grille contre les anneaux de la cellule, mis en cache."""
        i, j = cell
        n = max(1, round(self.cell_deg / self.raster_deg))
        lats = i * self.cell_deg + (np.arange(n) + 0.5) * self.raster_deg
        lons = j * self.cell_deg + (np.arange(n) + 0.5) * self.raster_deg
        gy, gx = np.meshgrid(lats, lons, indexing="ij")
        inside = np.zeros(gy.shape, dtype=bool)
        for idx in self._grid.get(cell, ()):
            inside |= _ring_contains(self._edge_index[idx], gy.ravel(), gx.ravel()).reshape(gy.shape)
        self._raster[cell] = inside
        return inside

    def crosses_land(self, lat1, lon1, lat2, lon2):
        """Vrai si le tronçon (lat1, lon1) -> (lat2, lon2) touche la terre.

        Le point de départ n'est pas testé ; l'arrivée et les points intermédiaires le sont.
        Pré-filtre : seuls les tronçons dont la bbox touche un anneau candidat sont
        échantillonnés (au large, aucun tronçon ne l'est -> coût nul)."""
        lat1, lon1, lat2, lon2 = np.broadcast_arrays(
            *(np.asarray(a, dtype=float) for a in (lat1, lon1, lat2, lon2))
        )
        if not self._rings or lat1.size == 0:
            return np.zeros(lat1.shape, dtype=bool)
        shape = lat1.shape
        la1, lo1, la2, lo2 = (a.ravel() for a in (lat1, lon1, lat2, lon2))
        dlon = (lo2 - lo1 + 180.0) % 360.0 - 180.0
        seg_la0 = np.minimum(la1, la2)
        seg_la1 = np.maximum(la1, la2)
        seg_lo0 = np.minimum(lo1, lo1 + dlon)
        seg_lo1 = np.maximum(lo1, lo1 + dlon)
        near = np.zeros(la1.shape, dtype=bool)
        # anneaux candidats : cellules de la grille couvrant le rectangle englobant
        # de tous les tronçons (un anneau qui coupe un tronçon est forcément
        # enregistré dans une de ces cellules)
        i_min = int(np.floor(seg_la0.min() / self.cell_deg))
        i_max = int(np.floor(seg_la1.max() / self.cell_deg))
        j_min = int(np.floor(seg_lo0.min() / self.cell_deg))
        j_max = int(np.floor(seg_lo1.max() / self.cell_deg))
        cand = {}
        for i in range(i_min, i_max + 1):
            for j in range(j_min, j_max + 1):
                for idx in self._grid.get((i, j), ()):
                    cand[idx] = True
        for idx in cand:
            _ring, (ra0, ra1, ro0, ro1) = self._rings[idx]
            near |= (seg_la1 >= ra0) & (seg_la0 <= ra1) & (seg_lo1 >= ro0) & (seg_lo0 <= ro1)
        hit = np.zeros(la1.shape, dtype=bool)
        sel = np.flatnonzero(near)
        if sel.size:
            dla = la2[sel] - la1[sel]
            dlo = dlon[sel]
            dist = np.hypot(dla, dlo * np.cos(np.radians((la1[sel] + la2[sel]) / 2.0)))
            dist_nm = float(dist.max()) * NM_PER_RADIAN * np.pi / 180.0
            n = max(1, int(np.ceil(dist_nm / self.sample_nm)))
            fs = (np.arange(1, n + 1, dtype=float) / n)[:, None]  # (n, 1)
            plat = (la1[sel] + dla * fs).ravel()  # (n, S) aplati
            plon = (lo1[sel] + dlo * fs).ravel()
            hit[sel] = self.is_land(plat, plon).reshape(n, sel.size).any(axis=0)
        return hit.reshape(shape)


def _ring_contains(edge_index, y, x, _edges=512):
    """Ray casting exact : un point de latitude y ne peut être traversé que par
    les arêtes de sa cellule de latitude (1°) — mêmes résultats que la boucle
    sur toutes les arêtes, mais vectorisée par blocs d'arêtes."""
    inside = np.zeros(y.shape, dtype=bool)
    cells = np.floor(y).astype(np.int64)
    for cell in np.unique(cells):
        entry = edge_index.get(int(cell))
        if entry is None:
            continue
        sel = cells == cell
        yb = y[sel][:, None]
        xb = x[sel][:, None]
        yi, xi, yj, xj = entry
        acc = np.zeros(yb.shape[0], dtype=bool)
        for es in range(0, yi.size, _edges):
            e_yi = yi[es:es + _edges][None, :]
            e_xi = xi[es:es + _edges][None, :]
            e_yj = yj[es:es + _edges][None, :]
            e_xj = xj[es:es + _edges][None, :]
            crosses = (e_yi > yb) != (e_yj > yb)
            with np.errstate(divide="ignore", invalid="ignore"):
                xint = (e_xj - e_xi) * (yb - e_yi) / (e_yj - e_yi) + e_xi
            acc ^= np.bitwise_xor.reduce(crosses & (xb < xint), axis=1)
        inside[sel] = acc
    return inside
