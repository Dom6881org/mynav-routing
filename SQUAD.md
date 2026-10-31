# Squad IA MyNavApp

Ce dépôt (**mynav-routing**, moteur de routing) fait partie de l'écosystème MyNavApp piloté par la [squad IA](https://github.com/Dom6881org/MyNavApp/blob/main/SQUAD.md) :

- 🎨 **Agent UI** — écrans et design (dépôt MyNavApp uniquement)
- ⚙️ **Agent Dev** — moteur de routing : GRIB, polaires, calculs (TWA, TWA/ETA), bathymétrie, unités
- ✅ **Agent Validation** — revue de chaque PR avant soumission à Dom
- 👑 **Dom Latge** — unique approbateur des PR

## Règles spécifiques à ce dépôt

1. Toute évolution du moteur passe par une PR sur ce dépôt (branch `feature/<slug>` ou `fix/<slug>` depuis `main`).
2. Après merge ici, l'intégration côté MyNavApp se fait par une PR séparée de mise à jour de la dépendance.
3. Aucun merge sans revue de l'Agent Validation puis **approbation de Dom**.
4. Une PR = un périmètre minimal, tests et CI au vert avant revue.
