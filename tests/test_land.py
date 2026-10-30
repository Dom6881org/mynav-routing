import io
import json
from urllib.error import URLError

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
