# mynav-routing v0.6.0 — Départ hors zone GRIB signalé

## Nouveautés depuis v0.5.0
- **`GRIBWindField.covers(lat, lon)`** : True si la position est dans l'emprise du
  champ (vectorisée, scalaires acceptés). `sample()` nearest-clippe silencieusement
  toute position hors emprise (vent du bord de grille, sans erreur) — `covers()`
  permet de le détecter : fin du clip silencieux (MyNavApp #141).
- **`RouteResult.start_in_grib`** : `True` (départ couvert) / `False` (départ hors
  zone : le vent du départ est celui du bord de grille) / `None` (vent sans emprise
  connue, ex. ConstantWind — back-compat complet). Détection par duck-typing dans
  `route()` : tout `WindField` exposant `covers()` est interrogé, aucun autre
  changement de comportement (corridor, fallback, moteur B1 inchangés).
- L'app MyNavApp surface l'avertissement côté serveur/panneau de routage dans une
  PR séparée (politique #56 : moteur d'abord, puis intégration).

## Technique
- Tests : `covers` (dedans/hors emprise/vectorisé), départ hors zone →
  `start_in_grib` False avec route trouvée au vent du bord, départ couvert → True,
  vent constant → None.

# mynav-routing v0.5.0 — Voile/moteur (B1)

## Nouveautés depuis v0.4.0
- **Relais moteur sous seuil de vitesse voile** (cadrage MyNavApp #133) :
  `route(engine_below_kts=...)` — dès que la vitesse polaire d'un tronçon passe sous le
  seuil (défaut 2.0 kt, paramétrable, `None` = pur voile), le tronçon passe en moteur à
  `ENGINE_SPEED_KTS` (4.0 kt), cap inchangé. La zone morte de la polaire (TWA < ~40°,
  vitesse 0) ne bloque plus le routage : la route au vent devient directe.
- `RouteResult.motor_flags` : une entrée bool par point (même convention que `sog_kts`,
  index 0 artificiel), `None` en pur voile (back-compat).
- `to_geojson()` expose `properties.motor` par point quand le relais est actif.
- `ENGINE_SPEED_KTS` exporté depuis le paquet.

## Technique
- Tests : seuil paramétrable (voile 2.1 kn vs moteur 4.0 kn), zone morte au vent
  (route introuvable en pur voile, trouvée en relais moteur), vent navigable
  (aucun moteur), alignement flags/sog_kts/coords, tronçons moteur à vitesse
  constante, GeoJSON avec et sans `motor`.

## Cadrage
Dom6881org/MyNavApp#133 (livrable 1/2 — l'app MyNavApp suit dans une PR séparée :
serveur + champ « Motor below (kt) » + trace moteur + re-pin v0.5.0).

---

# mynav-routing v0.4.0 — Performances & robustesse du routage côtier

## Nouveautés depuis v0.3.3
- **Sampling GRIB numpy** (PR #24, issue #16) : pré-extraction du champ + `sample()` vectorisé — 17,1× mesuré sur le micro-bench dédié `GRIBWindField.sample()` ; interface publique inchangée.
- **Corridor d'élagage** (PR #25, issue #17) : corridor ellipse départ-arrivée (marge 20 %, configurable via `corridor_margin`, `None` = off) + clés de cellules int32 — bench ~2× (5,56 M → 2,05 M nœuds).
- **Perf crosses_land** (PR #27, issue #26) : pré-filtre spatial (bucketing par cellule + échantillonnage par segment) — ≥3×, égalité stricte avec l'implémentation legacy.
- **Fallback sans corridor** (PR #29, issue #28) : relance automatique sans corridor quand `found=False` — les routes contournant un cap (Cap Béar, Cap de Creus) ne sont plus refusées à tort ; `corridor_margin` effectif exposé dans `RouteResult`.
- **Campagne d'acceptance** (PR #31, issue #30) : à chaque push CI, 3 routes réelles (côte ouverte Banyuls→Collioure, contournements Cap Béar ~3 MN et Cap de Creus 22,5 MN) avec assertions `found` + **aucun segment de trace sur terre** + ETA ≤ horizon. Matrice opt-in `routing-acceptance` (`workflow_dispatch`) : 3 routes × 3 vents + traversée Banyuls→Calvi (252,3 MN).

## Technique
- Bench de référence opt-in `BENCH_STCYPRIEN_AJACCIO` : 1,9 s bout en bout (sample 0,8 % du total).
- CI : ruff (line-length 100, I + RUF) + pytest suite rapide.
- Le paquet porte désormais la vraie version `0.4.0` (jusqu'ici figé à `0.3.3` à tous les SHA).

## Issues fermées
#16, #17, #19, #26, #28, #30
