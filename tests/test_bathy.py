"""Tests de la bathymétrie GEBCO (DepthMask) : chargement, interpolation, contrainte
de profondeur minimale, intégration dans route()."""

import io
from pathlib import Path
from urllib.error import URLError

import numpy as np
import pytest

import mynav_routing.bathy as bathy_module
from mynav_routing import BathyDataError, ConstantWind, DepthMask, Polar, route
from mynav_routing.bathy import BATHY_FILE_ENV, BATHY_URL_ENV
from mynav_routing.land import CACHE_DIR_ENV

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"
URL = "https://example.test/gebco_extrait.npz"

# grille régionale test : lat 41.5..43.5, lon 2..5 à 0.01°
LATS = np.arange(41.5, 43.5001, 0.01)
LONS = np.arange(2.0, 5.0001, 0.01)


def _elev():
    gy, gx = np.meshgrid(LATS, LONS, indexing="ij")
    elev = np.full(gy.shape, -4000.0, dtype=np.float32)  # grand large
    # île (terre) : lat 42.55..42.70, lon 3.20..3.30
    elev[(gy >= 42.55) & (gy <= 42.70) & (gx >= 3.20) & (gx <= 3.30)] = 5.0
    # haut-fond (3 m) : lat 42.58..42.66, lon 3.35..3.45
    elev[(gy >= 42.58) & (gy <= 42.66) & (gx >= 3.35) & (gx <= 3.45)] = -3.0
    return elev


def _npz_bytes():
    buf = io.BytesIO()
    np.savez(buf, elev=_elev(), lats=LATS, lons=LONS)
    return buf.getvalue()


@pytest.fixture
def sea():
    return DepthMask(LATS, LONS, _elev())


@pytest.fixture(autouse=True)
def _cache_isole(tmp_path, monkeypatch):
    monkeypatch.setenv(CACHE_DIR_ENV, str(tmp_path / "cache"))
    monkeypatch.delenv(BATHY_FILE_ENV, raising=False)
    monkeypatch.delenv(BATHY_URL_ENV, raising=False)
    bathy_module._fetch_bathy.cache_clear()
    yield
    bathy_module._fetch_bathy.cache_clear()


# --- chargement ---


def test_npz_roundtrip(tmp_path):
    p = tmp_path / "extrait.npz"
    np.savez(p, elev=_elev(), lats=LATS, lons=LONS)
    m = DepthMask.from_npz(p)
    assert m.elev.shape == (LATS.size, LONS.size)
    assert m.is_shallow(42.62, 3.40, 10.0)
    assert m.elevation_at(42.0, 2.5) == pytest.approx(-4000.0, abs=1.0)


def test_gebco_netcdf(tmp_path):
    nc4 = pytest.importorskip("netCDF4")
    # GEBCO : lat décroissant, élévation int16, variable 'elevation', dims lat/lon
    lats = LATS[::-1].copy()
    lons = LONS.copy()
    elev = _elev()[::-1, :].astype(np.int16)
    p = tmp_path / "gebco_2025_extrait.nc"
    with nc4.Dataset(p, "w", format="NETCDF4") as ds:
        ds.createDimension("lat", lats.size)
        ds.createDimension("lon", lons.size)
        vlat = ds.createVariable("lat", "f8", ("lat",))
        vlat[:] = lats
        vlon = ds.createVariable("lon", "f8", ("lon",))
        vlon[:] = lons
        ve = ds.createVariable("elevation", "i2", ("lat", "lon"))
        ve[:] = elev
    m = DepthMask.from_gebco_nc(p)
    assert m.lats[0] < m.lats[-1]  # remis en croissant
    assert m.elevation_at(42.6, 3.4) == pytest.approx(-3.0, abs=0.01)  # haut-fond
    assert m.is_shallow(42.6, 3.25, 10.0)  # île
    assert not m.is_shallow(42.0, 2.5, 10.0)


def test_fichier_local_impose_par_env(monkeypatch, tmp_path):
    p = tmp_path / "gebco.npz"
    np.savez(p, elev=_elev(), lats=LATS, lons=LONS)
    monkeypatch.setenv(BATHY_FILE_ENV, str(p))

    def _boom(url, timeout):
        raise AssertionError("urlopen ne doit pas être appelé")

    monkeypatch.setattr(bathy_module, "urlopen", _boom)
    m = DepthMask.from_open_data(url=URL)
    assert m.is_shallow(42.62, 3.40, 10.0)


def test_fichier_local_illisible(monkeypatch, tmp_path):
    monkeypatch.setenv(BATHY_FILE_ENV, str(tmp_path / "absent.npz"))
    with pytest.raises(BathyDataError):
        DepthMask.from_open_data()


def test_sans_config_leve_bathy_data_error():
    with pytest.raises(BathyDataError):
        DepthMask.from_open_data()
    # ValueError : le backend MyNavApp traduit déjà cela en RoutingError
    assert issubclass(BathyDataError, ValueError)


def test_telechargement_puis_cache_sur_disque(monkeypatch):
    monkeypatch.setattr(
        bathy_module, "urlopen", lambda url, timeout: io.BytesIO(_npz_bytes())
    )
    m = DepthMask.from_open_data(url=URL)
    assert m.is_shallow(42.62, 3.40, 10.0)
    assert bathy_module._bathy_cache_path(URL).is_file()


def test_hors_ligne_repli_sur_le_cache(monkeypatch):
    monkeypatch.setattr(
        bathy_module, "urlopen", lambda url, timeout: io.BytesIO(_npz_bytes())
    )
    DepthMask.from_open_data(url=URL)  # remplit le cache disque
    bathy_module._fetch_bathy.cache_clear()

    def _offline(url, timeout):
        raise URLError("offline")

    monkeypatch.setattr(bathy_module, "urlopen", _offline)
    m = DepthMask.from_open_data(url=URL)
    assert m.is_shallow(42.62, 3.40, 10.0)


def test_hors_ligne_sans_cache_leve_bathy_data_error(monkeypatch):
    def _offline(url, timeout):
        raise URLError("offline")

    monkeypatch.setattr(bathy_module, "urlopen", _offline)
    with pytest.raises(BathyDataError):
        DepthMask.from_open_data(url=URL)


def test_grille_invalide():
    with pytest.raises(BathyDataError):
        DepthMask([0.0, 1.0, 3.0], [0.0, 1.0], np.zeros((3, 2)))  # lat non régulière
    with pytest.raises(BathyDataError):
        DepthMask([0.0, 1.0], [0.0, 1.0], np.zeros((3, 2)))  # forme incohérente


# --- interpolation ---


def test_interpolation_bilineaire():
    m = DepthMask([0.0, 1.0], [0.0, 1.0], np.array([[-10.0, -20.0], [-30.0, -40.0]]))
    assert m.elevation_at(0.0, 0.0) == pytest.approx(-10.0)
    assert m.elevation_at(0.5, 0.5) == pytest.approx(-25.0)
    assert m.elevation_at(0.25, 0.75) == pytest.approx(-22.5)
    assert m.depth_at(0.5, 0.5) == pytest.approx(25.0)
    assert np.isnan(m.elevation_at(5.0, 5.0))  # hors grille
    assert np.isnan(m.elevation_at(-1.0, 0.0))


def test_elevation_at_hors_grille(sea):
    assert np.isnan(sea.elevation_at(50.0, 10.0))
    assert sea.elevation_at(42.6, 3.4) == pytest.approx(-3.0, abs=0.01)
    assert sea.elevation_at(42.0, 2.5) == pytest.approx(-4000.0, abs=1.0)


def test_traversee_du_meridien_360():
    # grille globale en longitude 0..360 : haut-fond de part et d'autre du méridien
    lats = np.arange(-10.0, 10.001, 0.5)
    lons = np.arange(0.0, 360.0, 1.0)
    gy, gx = np.meshgrid(lats, lons, indexing="ij")
    elev = np.full(gy.shape, -4000.0, np.float32)
    elev[(gx >= 359.0) | (gx <= 1.0)] = -2.0
    m = DepthMask(lats, lons, elev)
    assert m.elevation_at(0.0, 359.5) == pytest.approx(-2.0, abs=0.6)  # à cheval sur 0/360
    assert m.elevation_at(0.0, -0.5) == pytest.approx(-2.0, abs=0.6)
    assert m.elevation_at(0.0, 180.0) == pytest.approx(-4000.0, abs=1.0)
    assert m.crosses_shallow(0.0, -0.5, 0.0, 0.5, 10.0)  # tronçon qui franchit le méridien
    assert not m.crosses_shallow(0.0, 10.0, 0.0, 11.0, 10.0)


# --- contrainte de profondeur ---


def test_is_shallow(sea):
    assert sea.is_shallow(42.62, 3.40, 10.0)       # haut-fond 3 m
    assert not sea.is_shallow(42.62, 3.40, 2.0)    # 3 m >= 2 m : autorisé
    assert sea.is_shallow(42.60, 3.25, 2.0)        # île = terre, bloquée à tout seuil
    assert not sea.is_shallow(42.0, 2.5, 10.0)     # -4000 m
    assert not sea.is_shallow(50.0, 10.0, 10.0)    # hors grille : aucune contrainte
    assert sea.is_shallow(42.62, 3.40)             # défaut = default_min_depth_m (10 m)


def test_crosses_shallow(sea):
    # tronçon lat 42.62 de lon 3.30 à 3.50 : coupe le haut-fond (3 m)
    assert sea.crosses_shallow(42.62, 3.30, 42.62, 3.50, 10.0)
    assert not sea.crosses_shallow(42.62, 3.30, 42.62, 3.50, 2.0)
    # tronçon qui coupe l'île : bloqué même avec un seuil faible
    assert sea.crosses_shallow(42.62, 3.10, 42.62, 3.25, 2.0)
    # au large et hors grille : rien
    assert not sea.crosses_shallow(42.0, 2.2, 42.0, 2.5, 10.0)
    assert not sea.crosses_shallow(50.0, 10.0, 50.0, 11.0, 10.0)
    # l'arrivée est testée même sur un tronçon court (point final dans le haut-fond)
    assert sea.crosses_shallow(42.62, 3.46, 42.62, 3.44, 10.0)
    # entrée vide
    assert sea.crosses_shallow([], [], [], []).shape == (0,)


def test_crosses_shallow_forme_tableau(sea):
    a = np.array([[42.62, 42.62], [42.0, 42.0]])
    b = np.array([[3.30, 4.0], [2.5, 4.0]])
    r = sea.crosses_shallow(a, b, a, b + 0.2)
    assert r.shape == (2, 2)
    assert r[0, 0] and not r[0, 1] and not r[1, 0] and not r[1, 1]


def test_index_spatial_cellules_dangereuses(sea):
    # seules les cellules qui couvrent l'île / le haut-fond (± 1 cellule) sont dangereuses
    danger = sea._danger_mask(10.0)
    assert danger.any()
    assert danger.sum() < danger.size
    # loin à l'est (dernière colonne de cellules, lon ~5°), aucune cellule dangereuse
    assert not danger[:, -1].any()


# --- intégration routeur ---


def test_route_evite_terres_et_hauts_fonds(sea):
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    kw = {
        "start": (42.62, 3.02),
        "finish": (42.62, 3.5),
        "polar": polar,
        "wind": wind,
        "horizon_h": 48.0,
        "dt": 1800.0,
        "bin_deg": 0.05,
        "arrival_nm": 5.0,
        "land": None,
    }
    res_free = route(bathy=None, **kw)
    res = route(bathy=sea, min_depth_m=10.0, **kw)
    assert res.found
    lat = np.array([c[0] for c in res.coords])
    lon = np.array([c[1] for c in res.coords])
    # aucun tronçon de la route ne touche terre / haut-fond (élevation <= -10 m partout)
    assert not sea.crosses_shallow(lat[:-1], lon[:-1], lat[1:], lon[1:], 10.0).any()
    assert (sea.elevation_at(lat, lon) <= -10.0).all()
    # la route contrainte fait un détour : ETA >= route libre (qui traverse l'île)
    assert res.eta_s >= res_free.eta_s - 1e-6


def test_route_min_depth_par_defaut(sea):
    # sans min_depth_m explicite, le défaut du masque s'applique (10 m)
    polar = Polar.from_file(str(EX))
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(
        start=(42.62, 3.02),
        finish=(42.62, 3.5),
        polar=polar,
        wind=wind,
        horizon_h=48.0,
        dt=1800.0,
        bin_deg=0.05,
        arrival_nm=5.0,
        land=None,
        bathy=sea,
    )
    assert res.found
    lat = np.array([c[0] for c in res.coords])
    lon = np.array([c[1] for c in res.coords])
    assert not sea.crosses_shallow(lat[:-1], lon[:-1], lat[1:], lon[1:], sea.default_min_depth_m).any()
