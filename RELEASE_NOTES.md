# mynav-routing v0.4.0 — Performances & robustesse du routage côtier

## Nouveautés depuis v0.3.3
- **Sampling GRIB numpy** (PR #24, issue #16) : pré-extraction du champ + `sample()` vectorisé — 17,1× mesuré sur le micro-bench dédié `GRIBWindField.sample()` ; interface publique inchangée.
- **Corridor d'élagage** (PR #25, issue #17) : corridor ellipse départ-arrivée (marge 20 %, configurable via `corridor_margin`, `None` = off) + clés de cellules int32 — bench ~2× (5,56 M → 2,05 M nœuds).
- **Perf crosses_land** (PR #27, issue #26) : pré-filtre spatial (bucketing par cellule + échantillonnage par segment) — ≥3×, égalité stricte avec l'implémentation legacy.
- **Fallback sans corridor** (PR #29, issue #28) : relance automatique sans corridor quand `found=False` — les routes contournant un cap (Cap Béar, Cap de Creus) ne sont plus refusées à tort ; `corridor_margin` effectif exposé dans `RouteResult`.
- **Campagne d'acceptation** (PR #31, issue #30) : à chaque push CI, 3 routes réelles (côte ouverte Banyuls→Collioure, contournements Cap Béar ~3 MN et Cap de Creus 22,5 MN) avec assertions `found` + **aucun segment de trace sur terre** + ETA ≤ horizon. Matrice opt-in `routing-acceptance` (`workflow_dispatch`) : 3 routes × 3 vents + traversée Banyuls→Calvi (252,3 MN).

## Technique
- Bench de référence opt-in `BENCH_STCYPRIEN_AJACCIO` : 1,9 s bout en bout (sample 0,8 % du total).
- CI : ruff (line-length 100, I + RUF) + pytest suite rapide.
- Le paquet porte désormais la vraie version `0.4.0` (jusqu'ici figé à `0.3.3` à tous les SHA).

## Issues fermées
#16, #17, #19, #26, #28, #30
