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


def compute_idr_to_core_distance_comparison():
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


    variant1 = config.wt_variant
    variant2 = config.p36r_variant

    print("IDR to CORE Distance Comparison: WT vs P36R")
    print("=" * 40)
    print(f"Variant 1: {variant1.xtc_path}")
    print(f"Variant 2:    {variant2.xtc_path}")

    print(f"Analysis: {variant1.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    variant1_traj = loader.process_variant(variant1, extract_ca_atoms=False, convert_units=False)

    print(f"  Total frames: {variant1_traj.n_frames}")

    # Run convergence analysis
    print("\nFinding exposed residues for ...")
    exposed_residues = analyzer.compute_exposed_residues(variant1_traj)
    
    print("\nFinding distances...")
    distances1 = analyzer.compute_idr_to_core_mean_distance(variant1_traj, exposed_residues)

    
    print(f"Analysis: {variant2.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    variant2_traj = loader.process_variant(variant2, extract_ca_atoms=False)

    print(f"  Total frames: {variant2_traj.n_frames}")

    # Run convergence analysis
    print("\nFinding exposed residues for ...")
    exposed_residues = analyzer.compute_exposed_residues(variant2_traj)
    
    print("\nFinding distances...")
    distances2 = analyzer.compute_idr_to_core_mean_distance(variant2_traj, exposed_residues)


    print("\nPloting...")
    visualizer.plot_core_idr_mean_distances_comparisson(
        label1=variant1.label,
        distances1=distances1,
        label2=variant2.label,
        distances2=distances2,
        title=f"IDR residues distance to CORE - {variant1.label} X {variant2.label}",
        filename=f"idr_core_distances_{variant1.name.lower()}_{variant2.name.lower()}.png",
        show=True
    )

    return [distances1, distances2]





def main():
    """Find IDR to CORE Distance for both WT, P36R..."""

    print("=" * 50)
    print("IDR to CORE Distance Analysis - ANNEXIN A11")
    print("=" * 50)

    # Analyze WT
    print("\n")
    distances = compute_idr_to_core_distance_comparison()

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
