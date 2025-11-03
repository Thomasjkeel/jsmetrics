"""Top-level package for jsmetrics."""

__author__ = """Tom Keel"""
__email__ = "thomasjames.keel@gmail.com"
__version__ = "0.3.1"

from jsmetrics import details_for_all_metrics
from jsmetrics.metrics import (
    jet_core_algorithms,
    jet_statistics,
    waviness_metrics,
)


__all__ = [
    "details_for_all_metrics",
    "jet_core_algorithms",
    "jet_statistics",
    "waviness_metrics",
]
