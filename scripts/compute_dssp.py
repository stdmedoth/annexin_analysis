#!/usr/bin/env python3
"""
DSSP algorithm for Annexin A11 Ensemble


Output:
- DSSP algorithm for an specific trajectory
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


def compute_dssp(variant_name: str):
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

    print(f"DSSP Analysis: {variant.label}")
    print("=" * 40)

    # Process trajectory
    print("Loading trajectory...")
    processed = loader.process_variant(variant, extract_ca_atoms=False)

    print(f"  Total frames: {processed.n_frames}")

    # Run convergence analysis
    print("\nAnalyzing...")
    dssp_data = analyzer.compute_dssp(processed)

    visualizer.plot_dssp(
        dssp_data,
        title=f"DSSP Analysis - {variant.label}",
        filename=f"dssp_{variant_name.lower()}.png"
    )

    return dssp_data





def main():
    """Run DSSP for both WT, P36R..."""

    print("=" * 50)
    print("DSSP Analysis - ANNEXIN A11")
    print("=" * 50)

    # Analyze WT
    print("\n")
    wt_dssp = compute_dssp("wt")

    # Analyze P36R
    print("\n")
    p36r_dssp = compute_dssp("p36r")


    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
