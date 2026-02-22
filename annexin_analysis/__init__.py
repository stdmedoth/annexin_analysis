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

from .config import AnnexinConfig, ProteinRegions, VariantConfig
from .trajectory import TrajectoryLoader, ProcessedTrajectory
from .analysis import ConformationalAnalyzer, RMSFResult, PCAResult, ConvergenceResult
from .contact_map import ContactMapAnalyzer, ContactMapResult, ContactComparisonResult
from .visualization import ConformationalVisualizer
from .utils import AminoAcidConverter, NPZConverter, TrajectoryFilter

__version__ = "1.0.0"
__author__ = "Calistu"

__all__ = [
    "AnnexinConfig",
    "ProteinRegions",
    "VariantConfig",
    "TrajectoryLoader",
    "ProcessedTrajectory",
    "ConformationalAnalyzer",
    "RMSFResult",
    "PCAResult",
    "ConvergenceResult",
    "ContactMapAnalyzer",
    "ContactMapResult",
    "ContactComparisonResult",
    "ConformationalVisualizer",
    "AminoAcidConverter",
    "NPZConverter",
    "TrajectoryFilter",
]
