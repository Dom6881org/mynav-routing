import io
import json
from urllib.error import URLError

import numpy as np
import pytest

import mynav_routing.land as land_module
from mynav_routing import LandMask
from mynav_routing.land import CACHE_DIR_ENV, LAND_FILE_ENV, LandDataError

URL = "https://example.test/ne_10m_land.geojson"

# carré de terre : lat 0..10, lon 0..10 (GeoJSON en [lon, lat])
GEOJSON = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]]
                ],
            },
        }
    ],
}


@pytest.fixture(autouse=True)
def _cache_isole(tmp_path, monkeypatch):
    monkeypatch.setenv(CACHE_DIR_ENV, str(tmp_path / "cache"))
    monkeypatch.delenv(LAND_FILE_ENV, raising=False)
    land_module._load_open_land_geojson.cache_clear()
    yield
    land_module._load_open_land_geojson.cache_clear()


def _offline(url, timeout):
    raise URLError("offline")


def _online(url, timeout):
    return io.BytesIO(json.dumps(GEOJSON).encode("utf-8"))


def test_masque_geojson():
    land = LandMask.from_geojson(GEOJSON)
    assert land.is_land(5.0, 5.0)
    assert not land.is_land(20.0, 20.0)
    # tronçon lat 5 de lon -5 à lon 15 : coupe le carré de terre
    assert land.crosses_land(5.0, -5.0, 5.0, 15.0)
    assert not land.crosses_land(20.0, -5.0, 20.0, 15.0)


def test_telechargement_puis_cache_sur_disque(monkeypatch):
    monkeypatch.setattr(land_module, "urlopen", _online)
    data = land_module._load_open_land_geojson(URL)
    assert data == GEOJSON
    assert land_module._cache_path(URL).is_file()


def test_hors_ligne_repli_sur_le_cache(monkeypatch):
    monkeypatch.setattr(land_module, "urlopen", _online)
    land_module._load_open_land_geojson(URL)  # remplit le cache disque
    land_module._load_open_land_geojson.cache_clear()
    monkeypatch.setattr(land_module, "urlopen", _offline)
    data = land_module._load_open_land_geojson(URL)
    assert data == GEOJSON
    land = LandMask.from_open_data(url=URL)
    assert land.is_land(5.0, 5.0)
    assert not land.is_land(-5.0, -5.0)


def test_hors_ligne_sans_cache_leve_land_data_error(monkeypatch):
    monkeypatch.setattr(land_module, "urlopen", _offline)
    with pytest.raises(LandDataError):
        land_module._load_open_land_geojson(URL)
    # ValueError : le backend MyNavApp traduit déjà cela en RoutingError
    assert issubclass(LandDataError, ValueError)


def test_fichier_local_impose_par_env(monkeypatch, tmp_path):
    local = tmp_path / "terres.geojson"
    local.write_text(json.dumps(GEOJSON), encoding="utf-8")
    monkeypatch.setenv(LAND_FILE_ENV, str(local))
    monkeypatch.setattr(land_module, "urlopen", _offline)  # ne doit pas être appelé
    land = LandMask.from_open_data(url=URL)
    assert land.is_land(5.0, 5.0)
    assert not land.is_land(-5.0, -5.0)


def test_fichier_local_illisible(monkeypatch, tmp_path):
    monkeypatch.setenv(LAND_FILE_ENV, str(tmp_path / "absent.geojson"))
    with pytest.raises(LandDataError):
        land_module._load_open_land_geojson(URL)


# --- index spatial (perf v0.3.2) ---

def _point_in_ring(y, x, ring):
    """Référence naïve (ray casting scalaire) pour valider l'index spatial."""
    inside = False
    j = len(ring) - 1
    for i in range(len(ring)):
        yi, xi = ring[i]
        yj, xj = ring[j]
        if (yi > y) != (yj > y):
            xint = (xj - xi) * (y - yi) / (yj - yi) + xi
            if x < xint:
                inside = not inside
        j = i
    return inside


def _polygones_test():
    """60 petits carrés éparpillés + 1 grand anneau (style continent)."""
    rng = np.random.default_rng(42)
    rings = []
    for _ in range(60):
        la = -80.0 + rng.uniform(0, 160)
        lo = -180.0 + rng.uniform(0, 359)
        s = 0.5 + rng.uniform(0, 2.0)
        rings.append([(la, lo), (la, lo + s), (la + s, lo + s), (la + s, lo)])
    rings.append([(-90.0, -180.0), (-90.0, 180.0), (-60.0, 180.0), (-60.0, -180.0)])
    return rings


def test_index_spatial_exact_vs_reference():
    rings = _polygones_test()
    land = LandMask(rings, cell_deg=2.0, raster_deg=None)
    rng = np.random.default_rng(7)
    lat = rng.uniform(-90, 90, 500)
    lon = rng.uniform(-180, 180, 500)
    got = land.is_land(lat, lon)
    expected = np.array([
        any(_point_in_ring(y, x, r) for r in rings) for y, x in zip(lat, lon)
    ])
    assert (got == expected).all()


def test_index_spatial_invariant_selon_taille_cellule():
    rings = _polygones_test()
    rng = np.random.default_rng(3)
    lat = rng.uniform(-90, 90, 300)
    lon = rng.uniform(-180, 180, 300)
    ref = LandMask(rings, cell_deg=0.5, raster_deg=None).is_land(lat, lon)
    for cell in (1.0, 2.0, 8.0, 30.0):
        got = LandMask(rings, cell_deg=cell, raster_deg=None).is_land(lat, lon)
        assert (got == ref).all(), f"cell_deg={cell}"


def test_index_spatial_reduit_les_candidats():
    # 50 carrés tous dans le Pacifique Sud, requête locale au large de la France
    rings = []
    for k in range(50):
        la = -50.0 + 0.2 * k
        lo = 100.0 + 0.5 * k
        rings.append([(la, lo), (la, lo + 0.4), (la + 0.4, lo + 0.4), (la + 0.4, lo)])
    land = LandMask(rings, cell_deg=2.0)
    lat = np.full(100, 42.6)
    lon = np.linspace(3.0, 3.6, 100)
    assert land._candidate_rings(lat, lon) == []  # aucune cellule ne contient d'anneau
    assert not land.is_land(lat, lon).any()
    # requête au cœur du Pacifique Sud : des candidats, et des terres détectées
    lat2 = np.full(100, -45.0)
    lon2 = np.linspace(105.0, 125.0, 100)
    assert 0 < len(land._candidate_rings(lat2, lon2)) < len(land._rings)
    assert land.is_land(lat2, lon2).any()


def test_crosses_land_forme_et_echantillonnage():
    land = LandMask.from_geojson(GEOJSON)  # carré lat 0..10, lon 0..10
    # scalaire -> scalaire
    assert land.crosses_land(5.0, -5.0, 5.0, 15.0)
    assert not land.crosses_land(20.0, -5.0, 20.0, 15.0)
    # tableau 2D : forme préservée (ligne 0 lat=5 coupe le carré, ligne 1 lat=20 non)
    a = np.array([[5.0, 5.0], [20.0, 20.0]])
    b = np.array([[-5.0, 15.0], [-5.0, 15.0]])
    r = land.crosses_land(a, b, a, b + 10.0)
    assert r.shape == (2, 2)
    assert r[0, 0] and not r[0, 1] and not r[1, 0] and not r[1, 1]
    # l'arrivée est testée même sur un tronçon court
    assert land.crosses_land(5.0, 9.9, 5.0, 10.1)


def test_crosses_land_cellules_et_echantillonnage_par_segment():
    rng = np.random.default_rng(26)
    rings = []
    for k in range(300):
        lat = 30.0 + (k % 20) * 0.5
        lon = -5.0 + (k // 20) * 0.5
        rings.append([(lat, lon), (lat, lon + 0.12), (lat + 0.12, lon + 0.12), (lat + 0.12, lon)])
    land = LandMask(rings, sample_nm=0.25, cell_deg=0.5, raster_deg=None)

    lat1 = 30.0 + rng.uniform(0, 10, 500)
    lon1 = -5.0 + rng.uniform(0, 7.5, 500)
    lat2 = lat1 + rng.uniform(-0.2, 0.2, 500)
    lon2 = lon1 + rng.uniform(-0.8, 0.8, 500)
    lat1[:300] = 30.06 + (np.arange(300) % 20) * 0.5
    lon1[:300] = -5.0 + (np.arange(300) // 20) * 0.5 - 0.2
    lat2[:300] = lat1[:300]
    lon2[:300] = lon1[:300] + rng.uniform(0.3, 0.8, 300)

    got = land.crosses_land(lat1, lon1, lat2, lon2)
    dlon = (lon2 - lon1 + 180.0) % 360.0 - 180.0
    dist = np.hypot(
        lat2 - lat1,
        dlon * np.cos(np.radians((lat1 + lat2) / 2.0)),
    )
    sample_counts = np.maximum(
        1, np.ceil(dist * land_module.NM_PER_RADIAN * np.pi / 180.0 / land.sample_nm).astype(int)
    )
    expected = np.zeros(lat1.shape, dtype=bool)
    for n in np.unique(sample_counts):
        group = np.flatnonzero(sample_counts == n)
        fs = (np.arange(1, n + 1, dtype=float) / n)[:, None]
        points_lat = (lat1[group] + (lat2[group] - lat1[group]) * fs).ravel()
        points_lon = (lon1[group] + dlon[group] * fs).ravel()
        expected[group] = land.is_land(points_lat, points_lon).reshape(n, group.size).any(axis=0)

    np.testing.assert_array_equal(got, expected)


def test_raster_fidèle_gros_anneau():
    # cercle de 3600 sommets, rayon 3° centré en (42, 3)
    t = np.linspace(0, 2 * np.pi, 3601)
    ring = list(zip(42.0 + 3.0 * np.sin(t), 3.0 + 3.0 * np.cos(t)))
    land = LandMask([ring])  # raster_deg=0.01 par défaut
    assert land.is_land(42.0, 3.0)          # centre
    assert not land.is_land(42.0, 10.0)      # loin à l'est
    assert not land.is_land(50.0, 3.0)       # loin au nord
    assert land.is_land(42.0, 3.0 + 2.5)     # à 2.5° du centre, dans le cercle
    assert not land.is_land(42.0, 3.0 + 3.5) # à 3.5°, hors du cercle
    # cohérence raster vs exact loin des bords (tolérance raster)
    exact = LandMask([ring], raster_deg=None)
    rng = np.random.default_rng(11)
    lat = rng.uniform(38, 46, 400)
    lon = rng.uniform(-1, 7, 400)
    got = land.is_land(lat, lon)
    ref = exact.is_land(lat, lon)
    d = np.hypot(lat - 42.0, (lon - 3.0) * np.cos(np.radians(42.0)))
    loin_du_bord = np.abs(d - 3.0) > 0.05  # > tolérance raster
    assert (got[loin_du_bord] == ref[loin_du_bord]).all()
