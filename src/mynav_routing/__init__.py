"""mynav-routing : moteur de routage météo pour MyNavApp (isochrones, polaires, GRIB)."""

from .bathy import BathyDataError, DepthMask
from .grib import ConstantWind, GRIBWindField, WindField
from .land import LandDataError, LandMask
from .polar import Polar
from .router import ENGINE_SPEED_KTS, RouteResult, route, to_geojson

__all__ = [
    "ENGINE_SPEED_KTS",
    "BathyDataError",
    "ConstantWind",
    "DepthMask",
    "GRIBWindField",
    "LandDataError",
    "LandMask",
    "Polar",
    "RouteResult",
    "WindField",
    "route",
    "to_geojson",
]
__version__ = "0.6.0"
