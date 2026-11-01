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
- v0.4 : courants de surface (Copernicus CMEMS / HYCOM) composés au SOG
- v0.5 : endpoint FastAPI montable dans le backend MyNavApp (POST /api/routing)
- v0.6 : **qualité de route** — anti-zigzag et côtes (voir "Diagnostic et améliorations" ci-dessous)
- v1.0 : calibration automatique de la polaire depuis les logs NMEA (feature record de MyNavApp)

## Diagnostic v0.3 et améliorations prévues (v0.6 : qualité de route)

Symptômes observés dans MyNavApp : la trace passe sur la terre et présente des
zigzags serrés. Causes identifiées dans le code actuel :

### Traces sur terre

1. **MyNavApp épingle un commit antérieur au masque terre** : `requirements-routing.txt`
   pointe sur `9e85fe7` (pré-v0.3). La version consommée ne connaît ni `LandMask` ni
   `DepthMask`. Action : réépingler sur le dernier main (ou tag) et vérifier que le
   backend passe `land=` par défaut (pas `land=None`).
2. **Natural Earth 10m est trop grossier près des côtes** : presqu'îles fines, îlots,
   ports et digues sont absents du 1:10 000 000. Améliorations : GSHHG full/high
   resolution ou polygones OSM (fichier local imposé via `MYNAV_ROUTING_LAND_FILE`,
   extrait régional embarqué à bord), et marge de sécurité par dilatation du masque
   (~0,5 MN) plutôt qu'un simple test d'intersection.
3. **Tolérance du raster** (`raster_deg=0.01`, ~0,3 MN près des côtes) et pas
   d'échantillonnage `sample_nm=0.25` : un tronçon qui longe la côte peut être
   accepté à cheval sur le trait de côte. Prévoir `raster_deg=None` (ray casting
   exact) dans un rayon côtier, et densifier `sample_nm` près des terres.

### Zigzag

1. **Élagage arbitraire par cellule lat/lon** : `np.unique(key, return_index=True)`
   garde le premier point rencontré dans chaque bin de 0,25°, pas le meilleur temps.
   Les implémentations de référence (qtVlm, Expedition, MaxSea) élaguent par **secteur
   de relèvement** depuis le départ (ou le long de l'orthodromie) en conservant le
   point de temps minimal par secteur — et stockent le temps d'arrivée par point.
2. **Caps quantifiés à 5°** : l'optimum oscille entre deux caps adjacents selon le
   pas de temps du GRIB → dents de scie. Améliorations : raffinement local (±5°
   autour du cap de l'étape précédente), ou caps exacts vers les bins voisins.
3. **Pas de post-traitement de la trace** : la sortie brute de l'isochrone est une
   chaîne de segments de `bs*dt` (1–3 MN à dt=30 min). Ajouter : simplification
   Douglas-Peucker **sous contraintes** (terre + bathymétrie re-vérifiées après
   simplification), puis resampling en segments de cap constant (loomings).
4. **Fin de parcours brutale** : le dernier tronçon vers l'arrivée est ajouté sans
   lissage ; prévoir un raccord orthodromique vérifié.

### Alternatives algorithmiques (écartées ou complémentaires)

- **Isochrones** (actuel) : standard du routage voile (qtVlm, Expedition, MaxSea) —
  à conserver, avec les raffinements ci-dessus.
- **A*/Dijkstra sur graphe en temps** : robuste, mais peine avec les caps
  dépendant du vent et coûteux en mémoire sur zone large.
- **Fast Marching anisotrope** (ordered upwind) : élégant sur grille continue,
  mais la contrainte de cap par rapport au vent (pas de remontée du vent au près
  sous le TWA cible) est difficile à imposer proprement.
- Recommandation : rester sur les isochrones, investir dans l'élagage par secteur,
  le meilleur-temps par cellule et le post-traitement de trace.

## Sources open data

- Vent : GFS (NOAA), ICON (DWD), ARPEGE/AROME (Météo-France), ECMWF Open Data
- Courants : Copernicus Marine (CMEMS), HYCOM / RTOFS
- Bathymétrie : GEBCO 2026, EMODnet, MNT HOMONIM (Shom)

## Tests

    pytest -q

## Benchmark Saint-Cyprien → Ajaccio

Le benchmark déterministe utilise les coordonnées fixes de Saint-Cyprien
`(42.62, 3.03)` et Ajaccio `(41.919, 8.738)` (environ 257 MN orthodromiques),
la polaire d'exemple, un vent constant et un masque terrestre grossier intégré au
test. Il ne télécharge ni GRIB, ni données de terre ou bathymétrie ; cette dernière
reste optionnelle. Il est marqué `slow` et exclu de `pytest -q` pour préserver le
CI rapide.

Exécution locale :

    pytest -o addopts= -m slow tests/test_perf_stcyprien.py -s

Le benchmark affiche une ligne `BENCH_STCYPRIEN_AJACCIO` avec le temps total,
les nœuds explorés et les temps de sampling du vent, de masque terrestre et
d'élagage. Le budget bloquant est de 60 s sur le runner CI x86. Pour le suivi sur
Raspberry Pi 4, utiliser provisoirement un ratio indicatif de 3× (180 s), non
bloquant et à confirmer par mesure sur le matériel. Dans Actions, le job de
benchmark est lancé manuellement via **Run workflow** avec
`run_routing_benchmark` activé ; il ne rallonge pas les runs CI ordinaires.

Note : l'estimation d'environ 190 MN dans l'issue ne correspond pas à ces ports :
leur distance orthodromique calculée avec le moteur est d'environ 257 MN.

Après chaque PR de performance, reporter dans le commentaire de référence de
[l'issue #14](https://github.com/Dom6881org/mynav-routing/issues/14) les lignes
de sortie avant/après, avec le commit testé et l'environnement. Les optimisations
doivent comparer les mêmes paramètres et runner.
