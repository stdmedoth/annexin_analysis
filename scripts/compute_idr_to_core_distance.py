#!/usr/bin/env python3
"""
Compute the distances between the IDR residues and the Core carbons for Annexin A11 Ensemble


Output:
- Plot for IDR residues
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ConformationalVisualizer,
)


def compute_idr_to_core_mean_distance(variant_name: str):
    """
    Return the DSSP data for a specific variant.

    Args:
        variant_name: 'wt', 'p36r', ...
    """
    # Initialize components
    config = AnnexinConfig()
    loader = TrajectoryLoader(config)
    analyzer = ConformationalAnalyzer(config)
    visualizer = ConformationalVisualizer(config)


    # Get variant
    if variant_name.lower() == "wt":
        variant = config.wt_variant
    elif variant_name.lower() == "wt2":
        variant = config.wt2_variant
    elif variant_name.lower() == "p36r":
        variant = config.p36r_variant
    else:
        raise ValueError(f"Unknown variant: {variant_name}")

    print(f"SASA Analysis: {variant.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    processed = loader.process_variant(variant, extract_ca_atoms=False)

    print(f"  Total frames: {processed.n_frames}")

    # Run convergence analysis
    print("\nFinding exposed residues...")
    exposed_residues = analyzer.compute_exposed_residues(processed)
    
    print("\nFinding distances...")
    distances = analyzer.compute_idr_to_core_mean_distance(processed, exposed_residues)

    print("\nPloting...")
    visualizer.plot_core_idr_mean_distances(
        distances,
        title=f"IDR residues distance to CORE - {variant.label}",
        filename=f"idr_core_distances_{variant_name.lower()}.png",
        show=True
    )

    return distances





def main():
    """Find IDR to CORE Distance for both WT, P36R..."""

    print("=" * 50)
    print("IDR to CORE Distance Analysis - ANNEXIN A11")
    print("=" * 50)

    # Analyze WT
    print("\n")
    wt_distances = compute_idr_to_core_mean_distance("wt")

    # Analyze P36R
    print("\n")
    p36r_distances = compute_idr_to_core_mean_distance("p36r")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
