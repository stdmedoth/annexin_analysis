#!/usr/bin/env python3
"""
RG Analysis for Annexin A11 Ensemble


Output:
- Radius of Gyration calculation result
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer
)


def compute_rg(variant_name: str):
    """
    Analyze convergence for a specific variant.

    Args:
        variant_name: 'wt', 'p36r', ...
    """
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)

    # Get variant
    if variant_name.lower() == "wt":
        variant = config.wt_variant
    elif variant_name.lower() == "wt2":
        variant = config.wt2_variant
    elif variant_name.lower() == "p36r":
        variant = config.p36r_variant
    else:
        raise ValueError(f"Unknown variant: {variant_name}")

    print(f"Radius of Gyration Analysis: {variant.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    processed = loader.process_variant(variant)

    print(f"  Total frames: {processed.n_frames}")

    # Run convergence analysis
    print("\nAnalyzing...")
    rg = analyzer.compute_rg(processed)

    return rg


def main():
    """Run convergence analysis for both WT, P36R..."""

    print("=" * 50)
    print("RADIUS OF GYRATION ANALYSIS - ANNEXIN A11")
    print("=" * 50)

    # Analyze WT
    print("\n")
    wt_rg = compute_rg("wt")

    # Analyze WT 2 (luis felipe/ geraldo data)
    print("\n")
    wt2_rg = compute_rg("wt2")

    # Analyze P36R
    print("\n")
    p36r_rg = compute_rg("p36r")

    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"WT RADIUS OF GYRATION: {wt_rg}")
    print(f"WT2 RADIUS OF GYRATION: {wt2_rg}")
    print(f"P36R RADIUS OF GYRATION: {p36r_rg}")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
