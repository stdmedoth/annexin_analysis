"""
Annexin A11 Conformational Analysis Package

A comprehensive toolkit for analyzing the conformational space of
Annexin A11 protein and its mutations using molecular dynamics trajectories.

This package provides object-oriented tools for:
- Trajectory loading and preprocessing
- RMSF (Root Mean Square Fluctuation) analysis
- PCA (Principal Component Analysis) of conformational space
- Contact map computation and comparison
- Convergence analysis
- Publication-ready visualizations
"""

from .config import AnnexinConfig, ProteinRegions
from .trajectory import TrajectoryLoader
from .analysis import ConformationalAnalyzer
from .contact_map import ContactMapAnalyzer
from .visualization import ConformationalVisualizer
from .utils import AminoAcidConverter, NPZConverter

__version__ = "1.0.0"
__author__ = "Calistu"

__all__ = [
    "AnnexinConfig",
    "ProteinRegions",
    "TrajectoryLoader",
    "ConformationalAnalyzer",
    "ContactMapAnalyzer",
    "ConformationalVisualizer",
    "AminoAcidConverter",
    "NPZConverter",
]
