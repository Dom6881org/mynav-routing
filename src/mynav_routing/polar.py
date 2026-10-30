"""Polaires de vitesse (format .pol OpenCPN/qtVlm) + interpolation bilinéaire.

Format attendu : 1re ligne = TWS (nœuds), lignes suivantes = TWA puis vitesses.
Les lignes annotées non numériques (en-tête, Up/Down = VMG) sont ignorées.
"""
from __future__ import annotations

import numpy as np


def _is_num(token):
    try:
        float(token)
        return True
    except ValueError:
        return False


class Polar:
    """Polaire de vitesse, interpolée bilinéairement.

    tws en nœuds, twa = angle relatif au vent en degrés (0 = au près,
    180 = vent arrière), symétrique bâbord/tribord.
    """

    def __init__(self, tws, twa, speeds):
        self.tws = np.asarray(tws, dtype=float)
        self.twa = np.asarray(twa, dtype=float)
        self.speeds = np.asarray(speeds, dtype=float)
        if self.speeds.shape != (len(self.twa), len(self.tws)):
            raise ValueError("Dimensions incohérentes : attendu (len(twa), len(tws))")

    @classmethod
    def from_file(cls, path):
        tws, rows = None, []
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.replace(",", " ").replace(";", " ").split()
                if not parts:
                    continue
                if not _is_num(parts[0]):
                    if tws is None:
                        # ligne d'en-tête : les nombres restants = TWS
                        vals = [float(x) for x in parts[1:] if _is_num(x)]
                        if vals:
                            tws = np.array(vals)
                    # sinon : ligne annotée (Up/Down = VMG) -> ignorée
                    continue
                vals = [float(x) for x in parts if _is_num(x)]
                if tws is None:
                    tws = np.array(vals)
                    continue
                if len(vals) >= 2:
                    rows.append(vals)
        if tws is None or not rows:
            raise ValueError("Format .pol non reconnu : " + str(path))
        n = len(tws)
        twa = np.array([r[0] for r in rows])
        speeds = np.zeros((len(rows), n))
        for i, r in enumerate(rows):
            data = r[1 : 1 + n]
            speeds[i, : len(data)] = data
        order = np.argsort(twa)
        return cls(tws, twa[order], speeds[order])

    def _clip(self, tws, twa):
        twa = np.abs((np.asarray(twa, dtype=float) + 180.0) % 360.0 - 180.0)
        tws = np.clip(np.asarray(tws, dtype=float), self.tws[0], self.tws[-1])
        twa = np.clip(twa, self.twa[0], self.twa[-1])
        return tws, twa

    def speed_grid(self, tws, twa):
        """Version vectorisée : retourne un tableau de vitesses (nœuds)."""
        tws, twa = self._clip(tws, twa)
        i = np.clip(np.searchsorted(self.tws, tws) - 1, 0, len(self.tws) - 2)
        j = np.clip(np.searchsorted(self.twa, twa) - 1, 0, len(self.twa) - 2)
        ft = (tws - self.tws[i]) / np.maximum(self.tws[i + 1] - self.tws[i], 1e-9)
        fa = (twa - self.twa[j]) / np.maximum(self.twa[j + 1] - self.twa[j], 1e-9)
        s00, s10 = self.speeds[j, i], self.speeds[j + 1, i]
        s01, s11 = self.speeds[j, i + 1], self.speeds[j + 1, i + 1]
        return (
            s00 * (1 - fa) * (1 - ft)
            + s10 * fa * (1 - ft)
            + s01 * (1 - fa) * ft
            + s11 * fa * ft
        )

    def speed(self, tws, twa):
        """Scalaire : vitesse (kt) pour un vent donné."""
        return float(np.asarray(self.speed_grid([tws], [twa]))[0])
