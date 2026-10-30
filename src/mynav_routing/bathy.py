"""Bathymétrie GEBCO : éviter les hauts-fonds (et la terre) en routage.

Grille régulière (lat x lon) d'élévation en mètres, convention GEBCO :
valeur négative = profondeur sous la surface, positive = terre.
Sources acceptées :
  - NetCDF GEBCO (grille globale 15 arc-secondes ou extrait régional :
    variable 'elevation', dimensions 'lat' / 'lon') via xarray + netCDF4
    (extra 'bathy') ;
  - fichier .npz léger (clés 'elev' (2D), 'lats', 'lons'), pratique pour un
    extrait embarqué à bord.

Gestion hors ligne (exécution 100 % offline à bord) : MYNAV_ROUTING_BATHY_FILE
impose un fichier local (prioritaire) ; sinon MYNAV_ROUTING_BATHY_URL (ou le
paramètre url=) est téléchargé une fois puis mis en cache sur disque (même
répertoire de cache que le masque de terre, MYNAV_ROUTING_CACHE_DIR) ; hors
connexion, le cache est utilisé en repli. Sans aucune donnée : BathyDataError.
La grille globale GEBCO fait ~1 Go : à bord, préférer un extrait régional.

Performance : l'interpolation de l'élévation est bilinéaire et vectorisée
(O(1) par point, ~10^7 points/s) ; un index spatial (max d'élévation par cellule
de cell_deg degrés, calculé une fois au chargement) permet à crosses_shallow de
ne jamais échantillonner les tronçons dont la bbox ne touche aucune cellule
"dangereuse" — au large, le coût est nul. Une cellule est sûre si même son
point le plus haut reste sous -min_depth_m.

Hors couverture de la grille, aucune contrainte n'est appliquée (océan profond
supposé) ; la terre est gérée par LandMask, mais is_shallow / crosses_shallow
traitent aussi l'élévation >= 0 comme bloquante.
"""

from __future__ import annotations

import hashlib
import io
import os
import tempfile
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np

from .engine.isochrones import NM_PER_RADIAN
from .land import _cache_dir

BATHY_FILE_ENV = "MYNAV_ROUTING_BATHY_FILE"
BATHY_URL_ENV = "MYNAV_ROUTING_BATHY_URL"

# Portail de téléchargement GEBCO (grille globale NetCDF + appli de sous-ensemble
# régional en netCDF / GeoTIFF / ESRI ASCII).
GEBCO_PORTAL_URL = "https://download.gebco.net/"

DEFAULT_MIN_DEPTH_M = 10.0


class BathyDataError(ValueError):
    """Bathymétrie indisponible (hors ligne et aucun cache/fichier local).

    Sous-classe de ValueError : le backend MyNavApp traduit déjà les ValueError
    levées par route() en message utilisateur (RoutingError).
    """


def _bathy_cache_path(url: str) -> Path:
    suffix = Path(urlparse(url).path).suffix.lower()
    if suffix not in (".nc", ".nc4", ".npz"):
        suffix = ".nc"
    digest = hashlib.sha1(url.encode("utf-8")).hexdigest()[:16]
    return _cache_dir() / f"{digest}{suffix}"


def _write_cache_bytes(path: Path, data: bytes) -> None:
    """Écrit le cache de façon atomique ; best effort (répertoire en lecture seule = ignoré)."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name + ".", suffix=".tmp")
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
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
def _fetch_bathy(url: str) -> bytes:
    """Télécharge l'URL et met le brut en cache ; hors ligne, repli sur le cache."""
    cache = _bathy_cache_path(url)
    try:
        with urlopen(url, timeout=60) as response:
            data = response.read()
    except (OSError, ValueError):
        data = None
    if data is not None:
        _write_cache_bytes(cache, data)
        return data
    if cache.is_file():
        try:
            return cache.read_bytes()
        except OSError:
            pass
    raise BathyDataError(
        "Bathymétrie indisponible : hors connexion et aucun cache local "
        f"({cache}). Définissez {BATHY_FILE_ENV} vers un fichier GEBCO local "
        f"(NetCDF ou .npz), ou {BATHY_URL_ENV} vers une URL ; portail GEBCO : "
        f"{GEBCO_PORTAL_URL} (grille globale ~1 Go : préférer un extrait régional)."
    )


def _extract_grid(lats, lons, elev) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Normalise une grille (lat montante, lon croissante, elev 2D float32)."""
    lats = np.asarray(lats, dtype=float).ravel()
    lons = np.asarray(lons, dtype=float).ravel()
    elev = np.asarray(elev)
    if elev.ndim > 2:
        elev = np.squeeze(elev)
    if elev.ndim != 2 or elev.shape != (lats.size, lons.size):
        raise BathyDataError(
            f"Grille invalide : elev {elev.shape} vs lat {lats.size} x lon {lons.size}"
        )
    if lats.size < 2 or lons.size < 2:
        raise BathyDataError("Grille trop petite (moins de 2 points par axe)")
    if lats[0] > lats[-1]:
        lats = lats[::-1]
        elev = elev[::-1, :]
    if lons[0] > lons[-1]:
        lons = lons[::-1]
        elev = elev[:, ::-1]
    dlat = np.diff(lats)
    dlon = np.diff(lons)
    if not (np.allclose(dlat, dlat[0], rtol=1e-3, atol=1e-9)
            and np.allclose(dlon, dlon[0], rtol=1e-3, atol=1e-9)):
        raise BathyDataError("Grille non régulière (lat/lon doivent être uniformes)")
    return lats, lons, np.ascontiguousarray(elev, dtype=np.float32)


def _grid_from_npz(fileobj) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(fileobj) as z:
        elev = z["elev"] if "elev" in z else z["elevation"]
        return _extract_grid(z["lats"], z["lons"], elev)


def _open_netcdf(source):
    """Ouvre un NetCDF GEBCO (chemin ou objet binaire) : netcdf4 d'abord, sinon auto."""
    import xarray as xr

    try:
        return xr.open_dataset(source, engine="netcdf4")
    except (ValueError, ImportError, OSError):
        return xr.open_dataset(source)


def _grid_from_netcdf(source) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with _open_netcdf(source) as ds:
        name = next((n for n in ("elevation", "elev", "z", "Band1") if n in ds), None)
        if name is None:
            raise BathyDataError(f"Variable d'élévation introuvable dans le NetCDF : {list(ds)}")
        lat_name = next((n for n in ("lat", "latitude") if n in ds.coords), None)
        lon_name = next((n for n in ("lon", "longitude") if n in ds.coords), None)
        if lat_name is None or lon_name is None:
            raise BathyDataError("Coordonnées lat/lon introuvables dans le NetCDF")
        return _extract_grid(ds[lat_name].values, ds[lon_name].values, ds[name].values)


class DepthMask:
    """Grille régulière d'élévation GEBCO ; requêtes vectorisées.

    elevation_at : interpolation bilinéaire (NaN hors couverture).
    is_shallow / crosses_shallow : contrainte de profondeur minimale —
    bloquent l'eau peu profonde (hauts-fonds) et la terre (élevation >= 0).
    """

    def __init__(self, lats, lons, elev, sample_nm=0.25, cell_deg=1.0,
                 default_min_depth_m=DEFAULT_MIN_DEPTH_M):
        """lats / lons : coordonnées 1D des centres de maille (degrés, toute convention
        de longitude) ; elev : tableau 2D (nlat, nlon) en mètres (GEBCO : négatif = fond).
        sample_nm : pas d'échantillonnage le long d'un tronçon (mn).
        cell_deg : taille des cellules de l'index spatial (degrés).
        default_min_depth_m : profondeur minimale utilisée si non précisée aux requêtes."""
        self.lats, self.lons, self.elev = _extract_grid(lats, lons, elev)
        self.sample_nm = float(sample_nm)
        self.default_min_depth_m = float(default_min_depth_m)
        self._nlat, self._nlon = self.elev.shape
        self._lat0 = float(self.lats[0])
        self._lon0 = float(self.lons[0])
        self._dlat = float(self.lats[1] - self.lats[0])
        self._dlon = float(self.lons[1] - self.lons[0])
        # grille globale en longitude (recouvre les 360°) -> bouclage au méridien 0/360
        self._lon_global = (self._nlon - 1) * self._dlon + self._dlon >= 360.0 - 0.5 * self._dlon
        self._flat = self.elev.ravel()
        # index spatial : max d'élévation par cellule (une cellule n'est "sûre" que
        # si son point le plus haut reste sous le seuil), élargi d'une maille grille
        # de chaque côté (conservateur : une cellule marquée sûre l'est vraiment)
        k_lat = max(1, round(float(cell_deg) / self._dlat))
        k_lon = max(1, round(float(cell_deg) / self._dlon))
        self._cell_lat = k_lat * self._dlat
        self._cell_lon = k_lon * self._dlon
        rs = np.arange(0, self._nlat, k_lat)
        cs = np.arange(0, self._nlon, k_lon)
        cell_max = np.maximum.reduceat(np.maximum.reduceat(self.elev, rs, axis=0), cs, axis=1)
        cell_max = cell_max.astype(np.float32)
        p = np.pad(cell_max, 1, mode="constant", constant_values=-np.inf)
        self._cell_max = np.maximum.reduce([
            p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:],
            p[1:-1, :-2], p[1:-1, 1:-1], p[1:-1, 2:],
            p[2:, :-2], p[2:, 1:-1], p[2:, 2:],
        ])
        self._danger_cache = {}  # seuil (m) -> booléen des cellules dangereuses

    # --- chargement ---

    @classmethod
    def from_arrays(cls, lats, lons, elev, **kw):
        """Depuis des tableaux numpy (lats, lons 1D ; elev 2D)."""
        return cls(lats, lons, elev, **kw)

    @classmethod
    def from_npz(cls, path, **kw):
        """Fichier .npz : clés 'elev' (2D), 'lats', 'lons' (extrait embarqué léger)."""
        try:
            with open(path, "rb") as handle:
                lats, lons, elev = _grid_from_npz(handle)
        except (OSError, ValueError, KeyError) as error:
            raise BathyDataError(f"{path} illisible comme .npz : {error}") from error
        return cls(lats, lons, elev, **kw)

    @classmethod
    def from_gebco_nc(cls, path, **kw):
        """Fichier NetCDF GEBCO (variable 'elevation', dims lat/lon) — extra 'bathy'."""
        try:
            lats, lons, elev = _grid_from_netcdf(str(path))
        except (OSError, ValueError) as error:
            raise BathyDataError(f"{path} illisible comme NetCDF GEBCO : {error}") from error
        return cls(lats, lons, elev, **kw)

    @classmethod
    def _from_file(cls, path: Path, **kw):
        if path.suffix.lower() == ".npz":
            return cls.from_npz(path, **kw)
        return cls.from_gebco_nc(path, **kw)

    @classmethod
    def _from_bytes(cls, data: bytes, suffix: str, **kw):
        if suffix == ".npz":
            try:
                lats, lons, elev = _grid_from_npz(io.BytesIO(data))
            except (ValueError, KeyError) as error:
                raise BathyDataError(f"contenu .npz illisible : {error}") from error
        else:
            try:
                lats, lons, elev = _grid_from_netcdf(io.BytesIO(data))
            except (OSError, ValueError) as error:
                raise BathyDataError(f"contenu NetCDF illisible : {error}") from error
        return cls(lats, lons, elev, **kw)

    @classmethod
    def from_open_data(cls, url=None, **kw):
        """Charge la bathymétrie open data : fichier local MYNAV_ROUTING_BATHY_FILE
        (prioritaire), sinon téléchargement de url (ou MYNAV_ROUTING_BATHY_URL)
        mis en cache, sinon repli sur le cache hors ligne, sinon BathyDataError."""
        local = os.environ.get(BATHY_FILE_ENV)
        if local:
            try:
                return cls._from_file(Path(local), **kw)
            except (OSError, BathyDataError) as error:
                raise BathyDataError(f"{BATHY_FILE_ENV}={local} illisible : {error}") from error
        url = url or os.environ.get(BATHY_URL_ENV)
        if url is None:
            raise BathyDataError(
                f"Aucune bathymétrie configurée : définissez {BATHY_FILE_ENV} (fichier "
                f"GEBCO local, ex. extrait régional embarqué) ou {BATHY_URL_ENV} / url= "
                f"(téléchargement mis en cache). Portail GEBCO : {GEBCO_PORTAL_URL} — la "
                "grille globale fait ~1 Go, préférer un extrait régional."
            )
        return cls._from_bytes(_fetch_bathy(url), Path(urlparse(url).path).suffix.lower(), **kw)

    # --- requêtes ---

    def _norm_lon(self, lon):
        """Ramène les longitudes dans la convention de la grille [lon0, lon0 + 360)."""
        return (np.asarray(lon, dtype=float) - self._lon0) % 360.0 + self._lon0

    def elevation_at(self, lat, lon):
        """Élévation (m, convention GEBCO) par interpolation bilinéaire ; NaN hors grille."""
        lat = np.asarray(lat, dtype=float)
        lon = self._norm_lon(lon)
        fi = (lat - self._lat0) / self._dlat
        fj = (lon - self._lon0) / self._dlon
        valid = (fi >= 0.0) & (fi <= self._nlat - 1)
        if self._lon_global:
            fj = fj % self._nlon  # toute longitude est couverte (bouclage au méridien)
        else:
            valid &= (fj >= 0.0) & (fj <= self._nlon - 1)
        i0 = np.clip(np.floor(fi), 0, self._nlat - 2).astype(np.int64)
        wi = np.clip(fi - i0, 0.0, 1.0)
        j0 = np.floor(fj).astype(np.int64)
        if self._lon_global:
            j0 = j0 % self._nlon
            j1 = (j0 + 1) % self._nlon
        else:
            j0 = np.clip(j0, 0, self._nlon - 2)
            j1 = j0 + 1
        wj = np.clip(fj - j0, 0.0, 1.0)
        flat = self._flat
        n = self._nlon
        e00 = flat[i0 * n + j0]
        e01 = flat[i0 * n + j1]
        e10 = flat[(i0 + 1) * n + j0]
        e11 = flat[(i0 + 1) * n + j1]
        out = (e00 * (1.0 - wi) * (1.0 - wj) + e01 * (1.0 - wi) * wj
               + e10 * wi * (1.0 - wj) + e11 * wi * wj)
        return np.where(valid, out, np.nan)

    def depth_at(self, lat, lon):
        """Profondeur (m, positive sous la surface) ; NaN hors grille, négative sur terre."""
        return -self.elevation_at(lat, lon)

    def is_shallow(self, lat, lon, min_depth_m=None):
        """Vrai si le point est couvert par la grille et moins profond que min_depth_m
        (ou sur la terre). Hors couverture : False (aucune contrainte)."""
        m = self.default_min_depth_m if min_depth_m is None else float(min_depth_m)
        e = self.elevation_at(lat, lon)
        return np.isfinite(e) & (e > -m)

    def _danger_mask(self, min_depth_m):
        """Booléen (nci, ncj) des cellules dont le max d'élévation peut être au-dessus
        de -min_depth_m (potentiellement peu profondes) — mis en cache par seuil."""
        m = float(min_depth_m)
        mask = self._danger_cache.get(m)
        if mask is None:
            mask = self._cell_max > -m
            self._danger_cache[m] = mask
        return mask

    def crosses_shallow(self, lat1, lon1, lat2, lon2, min_depth_m=None):
        """Vrai si le tronçon (lat1, lon1) -> (lat2, lon2) touche une eau de profondeur
        < min_depth_m ou la terre (élevation > -min_depth_m, NaN hors grille exclu).

        Le point de départ n'est pas testé ; l'arrivée et les points intermédiaires le sont.
        Pré-filtre : seuls les tronçons dont la bbox touche une cellule dangereuse sont
        échantillonnés (au large, aucune cellule dangereuse -> coût nul)."""
        m = self.default_min_depth_m if min_depth_m is None else float(min_depth_m)
        lat1, lon1, lat2, lon2 = np.broadcast_arrays(
            *(np.asarray(a, dtype=float) for a in (lat1, lon1, lat2, lon2))
        )
        if lat1.size == 0:
            return np.zeros(lat1.shape, dtype=bool)
        shape = lat1.shape
        la1, lo1, la2, lo2 = (a.ravel() for a in (lat1, lon1, lat2, lon2))
        danger = self._danger_mask(m)
        hit = np.zeros(la1.shape, dtype=bool)
        if not danger.any():
            return hit.reshape(shape)
        nci, ncj = danger.shape
        dlon = (lo2 - lo1 + 180.0) % 360.0 - 180.0
        seg_la0 = np.minimum(la1, la2)
        seg_la1 = np.maximum(la1, la2)
        lo_a = self._norm_lon(lo1)
        seg_lo0 = np.minimum(lo_a, lo_a + dlon)
        seg_lo1 = np.maximum(lo_a, lo_a + dlon)
        # cellules de la bbox de chaque tronçon : un tronçon court (routage) ne
        # couvre que quelques cellules -> test O(1) par tronçon via la grille
        # de danger ; les tronçons longs (peu nombreux) sont échantillonnés d'office
        ci0 = np.clip(np.floor((seg_la0 - self._lat0) / self._cell_lat), 0, nci - 1).astype(np.int64)
        ci1 = np.clip(np.floor((seg_la1 - self._lat0) / self._cell_lat), 0, nci - 1).astype(np.int64)
        span_i = ci1 - ci0
        span_j = np.floor((seg_lo1 - self._lon0) / self._cell_lon).astype(np.int64) \
            - np.floor((seg_lo0 - self._lon0) / self._cell_lon).astype(np.int64)
        k_max = 4
        near = (span_i > k_max) | (span_j > k_max)
        for shift in (0.0, -360.0, 360.0):
            a = np.floor((seg_lo0 + shift - self._lon0) / self._cell_lon).astype(np.int64)
            if self._lon_global:
                a %= ncj
            else:
                a = np.clip(a, 0, ncj - 1)
            for di in range(k_max + 1):
                ii = np.clip(ci0 + di, 0, nci - 1)
                mi = di <= span_i
                for dj in range(k_max + 1):
                    jj = (a + dj) % ncj if self._lon_global else np.clip(a + dj, 0, ncj - 1)
                    near |= mi & (dj <= span_j) & danger[ii, jj]
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
            elev = self.elevation_at(plat, plon)
            hit[sel] = (elev > -m).reshape(n, sel.size).any(axis=0)  # NaN > -m = False
        return hit.reshape(shape)
