# Squad IA MyNavApp

Ce dépôt (**mynav-routing**, moteur de routing) fait partie de l'écosystème MyNavApp piloté par la [squad IA](https://github.com/Dom6881org/MyNavApp/blob/main/SQUAD.md) :

- 🧭 **Agent Dev Routing** — agent **principal de ce dépôt** : moteur de routing (isochrones, polaires, calculs TWA/ETA/SOG, bathymétrie, unités) et performance du moteur (issues #14, #16-#20)
- 🎨 **Agent Dev Front-end** — écrans, composants, hooks, design (dépôt MyNavApp uniquement)
- ⚙️ **Agent Dev Back-end** — serveur Python de MyNavApp : API, NMEA, AIS, GRIB côté app (dépôt MyNavApp uniquement)
- 🔗 **Agent Intégrateur** — intégration inter-dépôts : après chaque release ici, PR de mise à jour de la dépendance côté MyNavApp ; tests E2E
- 🏛️ **Agent Architecte** — décisions d'architecture et cohérence de l'écosystème (issues, ADR, revues transversales) ; ne code pas
- ✅ **Agent Validation** — revue de chaque PR avant soumission à Dom
- 👑 **Dom Latge** — unique approbateur des PR

Le **responsable d'un ticket** peut faire appel aux autres spécialistes si le périmètre l'exige (contribution à sa PR, ou découpage en sous-issues par l'Architecte).

## Règles spécifiques à ce dépôt

1. Toute évolution du moteur passe par une PR sur ce dépôt (branch `feature/<slug>` ou `fix/<slug>` depuis `main`) ; l'agent principal est le **Dev Routing**.
2. Après merge ici, l'**Agent Intégrateur** ouvre la PR côté MyNavApp de mise à jour de la dépendance.
3. Aucun merge sans revue de l'Agent Validation puis **approbation de Dom**.
4. Une PR = un périmètre minimal, tests et CI au vert avant revue.
