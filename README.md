# mynav-routing

Moteur de routage météo pour MyNavApp — méthode des isochrones, polaires .pol,
champs de vent GRIB. 100 % open data (cible : Raspberry Pi à bord).

Package Python autonome : le backend MyNavApp l'importe comme bibliothèque.
Sens de la dépendance : MyNavApp -> mynav-routing, jamais l'inverse.

## v0.1

- Parseur de polaires .pol (format OpenCPN / qtVlm), tolérant (lignes Up/Down = VMG ignorées)
- Interpolation bilinéaire TWS x TWA, symétrique bâbord/tribord
- Moteur d'isochrones vectorisé numpy (quelques secondes sur un Pi 4)
- Champs de vent : ConstantWind (tests, dégradation offline) et GRIBWindField (extra "grib", expérimental)
- Contraintes : TWS max (confort), rayon d'arrivée
- Masque terrestre Natural Earth téléchargé automatiquement, puis calcul local
- Export GeoJSON pour la carte Leaflet de MyNavApp

## v0.3 : masque terre + bathymétrie GEBCO

- Masque de terre Natural Earth (voir v0.1) accéléré : index spatial par cellule + raster,
  routage avec masque ~x10 plus rapide (v0.3.2)
- Bathymétrie GEBCO (grille 15 arc-secondes, NetCDF) : `DepthMask` évite la terre ET les
  hauts-fonds — tout tronçon dont la profondeur passe sous `min_depth_m` (défaut 10 m) est
  rejeté, via `route(..., bathy=DepthMask.from_open_data(), min_depth_m=10)`
- Chargement : fichier NetCDF GEBCO (extra `bathy`) ou `.npz` léger (`elev`, `lats`, `lons`)
  pour un extrait régional embarqué ; `MYNAV_ROUTING_BATHY_FILE` (fichier local prioritaire),
  `MYNAV_ROUTING_BATHY_URL` (téléchargement mis en cache, repli hors ligne dans le même cache
  que les terres). Sans donnée : `BathyDataError`. La grille globale GEBCO fait ~1 Go : à
  bord, préférer un extrait régional. Hors couverture de la grille : aucune contrainte
- Performance : interpolation bilinéaire vectorisée (O(1) par point) + index spatial
  (min d'élévation par cellule) : au large, les tronçons ne sont pas échantillonnés ;
  test de performance : `pytest -q tests/test_bathy_perf.py`

## Installation

    pip install -e .
    pip install -e '.[grib]'   # lecture GRIB (xarray + cfgrib + eccodes)
    pip install -e '.[bathy]'  # lecture NetCDF GEBCO (xarray + netcdf4)
    pip install -e '.[dev]'    # pytest + ruff

## Exemple

    from mynav_routing import ConstantWind, Polar, route, to_geojson

    polar = Polar.from_file("examples/polaire-sugar.pol")
    wind = ConstantWind(direction_from_deg=270, speed_kts=12)
    res = route(start=(42.62, 3.02), finish=(42.62, 3.5), polar=polar, wind=wind)
    if res.found:
        print("ETA (h):", res.eta_s / 3600.0)
        geojson = to_geojson(res)   # FeatureCollection LineString (lon, lat)

## Contrat d'intégration MyNavApp

- Départ : position courante NMEA (RMC/VTG) fournie par le backend
- Polaire : examples/polaire-sugar.pol (à recalibrer avec les logs enregistrés par MyNavApp)
- Vent : GRIB téléchargé à terre (GFS / ICON / AROME, open data) ou ConstantWind
- Terres : les polygones Natural Earth sont téléchargés automatiquement au premier routage,
  puis mis en cache sur disque (~/.cache/mynav-routing, ou MYNAV_ROUTING_CACHE_DIR).
  Hors ligne, le cache est utilisé en repli ; MYNAV_ROUTING_LAND_FILE permet d'imposer
  un fichier GeoJSON local (ex. embarqué à bord). Sans aucune donnée : LandDataError.
- Bathymétrie : `bathy=` est optionnel (None par défaut = pas de contrainte de profondeur) ;
  passer un `DepthMask` (GEBCO) pour éviter terre et hauts-fonds, avec `min_depth_m` (m).
- Sortie : GeoJSON affiché sur la carte, recalcul au fil de la traversée

Le calcul du masque est local. Pour désactiver explicitement la contrainte terrestre,
passer `land=None` à `route()` ; la contrainte de profondeur se désactive avec `bathy=None`.

## Feuille de route

- v0.2 : GRIBWindField stabilisé + tests avec GRIB réel (GFS 0.25°)
- v0.3 : masque terre / bathymétrie GEBCO (éviter terre et hauts-fonds) — fait
- v0.6 : **qualité de route** — anti-zigzag et côtes (voir "Diagnostic et améliorations" ci-dessous)
- v0.4 : courants de surface (Copernicus CMEMS / HYCOM) composés au SOG
- v0.5 : endpoint FastAPI montable dans le backend MyNavApp (POST /api/routing)
- v1.0 : calibration automatique de la polaire depuis les logs NMEA (feature record de MyNavApp)

## Sources open data

- Vent : GFS (NOAA), ICON (DWD), ARPEGE/AROME (Météo-France), ECMWF Open Data
- Courants : Copernicus Marine (CMEMS), HYCOM / RTOFS
- Bathymétrie : GEBCO 2026, EMODnet, MNT HOMONIM (Shom)

## Tests

    pytest -q
