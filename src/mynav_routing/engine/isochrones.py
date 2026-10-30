"""Géodésie et propagation : utilitaires partagés par le moteur d'isochrones."""
from __future__ import annotations

import numpy as np

NM_PER_RADIAN = 3437.75  # rayon terrestre moyen en milles nautiques


def haversine_nm(lat1, lon1, lat2, lon2):
    """Distance orthodromique en milles nautiques (entrées en degrés)."""
    p1 = np.radians(np.asarray(lat1, dtype=float))
    p2 = np.radians(np.asarray(lat2, dtype=float))
    dl = np.radians(np.asarray(lon2, dtype=float) - np.asarray(lon1, dtype=float))
    a = np.sin((p2 - p1) / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    return 2.0 * NM_PER_RADIAN * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def move_nm(lat, lon, heading_deg, dist_nm):
    """Avance le long d'un grand cercle. Cap 0 = nord. Retourne (lat, lon) en degrés."""
    d = np.radians(np.asarray(dist_nm, dtype=float) / NM_PER_RADIAN)
    tc = np.radians(np.asarray(heading_deg, dtype=float))
    p1 = np.radians(np.asarray(lat, dtype=float))
    l1 = np.radians(np.asarray(lon, dtype=float))
    p2 = np.arcsin(
        np.clip(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(tc), -1.0, 1.0)
    )
    l2 = l1 + np.arctan2(
        np.sin(tc) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2)
    )
    return np.degrees(p2), (np.degrees(l2) + 540.0) % 360.0 - 180.0
