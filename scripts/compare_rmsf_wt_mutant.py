#!/usr/bin/env python3
"""
RMSF Comparison: Wild-Type vs Mutant P36R

This script compares the RMSF profiles of wild-type Annexin A11
and the P36R mutant to identify regions of altered flexibility.

The comparison helps prove that the P36R mutation affects the
conformational dynamics of Annexin A11.

Output:
- Overlaid RMSF profiles for WT and P36R
- Statistical comparison (mean RMSF for each variant)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from annexin_analysis import (
    AnnexinConfig,
    VariantConfig,
    TrajectoryLoader,
    ConformationalAnalyzer,
    ConformationalVisualizer
)
from annexin_analysis.analysis import ComparativeAnalyzer


def main():
    """Compare RMSF profiles between WT and P36R mutant."""

    # Initialize components
    config = AnnexinConfig()
    comparative = ComparativeAnalyzer(config)
    visualizer = ConformationalVisualizer(config)

    # Define variants to compare
    variants = [config.wt_variant, config.p36r_variant]

    print("RMSF Comparison: WT vs P36R")
    print("=" * 40)

    # Compute RMSF for all variants
    rmsf_results = comparative.compare_rmsf(variants)

    # Create comparison dictionary with nice labels
    labeled_results = {}
    for variant in variants:
        if variant.name in rmsf_results:
            labeled_results[variant.label] = rmsf_results[variant.name]

    # Visualize comparison
    print("\nGenerating comparison plot...")
    visualizer.plot_rmsf_comparison(
        labeled_results,
        title="RMSF Comparison: Wild-Type vs P36R - Annexin A11",
        filename="rmsf_comparison_wt_p36r.png"
    )

    # Print summary statistics
    print("\n" + "=" * 40)
    print("Summary:")
    print("=" * 40)

    if len(rmsf_results) >= 2:
        wt_mean = rmsf_results["wt"].mean_rmsf
        p36r_mean = rmsf_results["p36r"].mean_rmsf
        diff = p36r_mean - wt_mean
        diff_percent = (diff / wt_mean) * 100

        print(f"WT Mean RMSF:   {wt_mean:.4f} Å")
        print(f"P36R Mean RMSF: {p36r_mean:.4f} Å")
        print(f"Difference:     {diff:+.4f} Å ({diff_percent:+.1f}%)")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
